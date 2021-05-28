#include "precomp.h"
#include "mcviewer.h"

#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <math.h>

Game* game = new MCViewer();

/*
Notes on the data structure used in this application:

- The world is 8192 x 256 x 8192 voxels in size, 1 byte per voxel,
  so a total of 16GB of raw data.
- This data is subdivided in 8x8x8 bricks. There are thus
  1024 x 32 x 1024 bricks in the world.
- A brick may be used multiple times, to save disk space. Therefore,
  a top-level grid contains 1024 x 32 x 1024 brick indices. For the
  test world, these are all unique, for now. The indices are stored
  in 'index.bin', which is a 128MB file.
- The actual brick data is stored in 256 'block000.bin' files. Each
  of these contains 131072 bricks, for a total of 64MB. Brick x
  can be found in block file x / 131072.

Special note:

- Brick idx 0 is a special one: this denotes an empty 8x8x8 area in
  the map. An obvious optimization of the data set would be to
  replace any empty 8x8x8 area with a zero index, to save space,
  and processing time.

Your task:

1. For the first 5 points:
   Improve the speed of the application using caching. You may cache
   parts of the index table, as well as parts of block files.
   IMPORTANT: The total amount of cached data may not exceed 8MB.
   This does not include any administratative data, just the actual
   cached data. Note that this means you cannot cache the full index
   table. Test your cache using the OptimizeWorld function.

2. Extra challenge 1 (1.5pts):
   Only applicable if you completed the first task.
   Prepare your cache for more generic access patterns. A test
   function will be provided later.

3. Extra challenge 2 (1.5pts):
   Only applicable if you completed the first and second task.
   Prepare your cache for reading and writing. Modify the 'WriteVoxel' 
   for this. A test function will be provided later.

4. Extra challenge 3a (1 pt):
   Only applicable if you completed the first and second task.
   Implement a cache hierarchy. Test the hierarchy using a 2MB L1$,
   a 4MB L2$ and a 16MB L3$.

   Extra challenge 3b (1 pt):
   Only applicable if you completed the first and second task.
   Implement a cache hierarchy that supports reading and writing. 
   Test the hierarchy using a 2MB L1$, a 4MB L2$ and a 16MB L3$.

5. Hardcore challenge (5 pts):
   Make your cache thread-safe.
   WARNING: This may very well be unreasonably hard.
   
6. For an extra point, regardless of accepted challenges:
   - Analyze the caching scheme: hits, misses, behavior over time.
   - Experiment with cache sizes, cache lines sizes, and
     and associativity, and report on this.
   Depending on the quality of your analysis, the full point or
   part of it may be awarded.

   DEADLINE: MAY 28, 17:00.
                                     "May the Light be with you!"
*/

struct IndexCacheLine
{
	int idx;
	uint res;
	bool valid, dirty;
};

struct BrickCacheLine
{
	int idx;
	uchar bricks[512];
	bool valid, dirty;
};

static int brickCount = 1024 * 1024 * 32; // all bricks are unique when we start
static int deletedBricks = 0;
static int blockFileCount = 256;
static int startCell = 1;

static long int cacheChecks = 0;
static long int cacheHits = 0;

static long int brickChecks = 0;
static long int brickHits = 0;

#define IndexSetSize 256
#define BrickSetSize 512

#define Index4Bytes (262144 / 4)
#define Brick4Bytes ((8388608 / 4) - Index4Bytes)

#define PrefetchIndex 8
#define PrefetchBricks 8

IndexCacheLine IndexCache[IndexSetSize][Index4Bytes / IndexSetSize];
BrickCacheLine BrickCache[BrickSetSize][Brick4Bytes / BrickSetSize];

FILE* indexFile;

int EvictFromIndexCache()
{
	return rand() % (Index4Bytes / IndexSetSize);
}

int EvictFromBrickCache()
{
	return rand() % (Brick4Bytes / BrickSetSize);
}

uint GetBrickIndexFromCache(int index)
{
	//return UINT32_MAX;
	cacheChecks++;
	int set = index % (IndexSetSize);
	auto cacheSet = IndexCache[set];
	for (int i = 0; i < 32 * 128; i++)
	{
		if (cacheSet[i].valid && cacheSet[i].idx == index)
		{
			cacheHits++;
			return cacheSet[i].res;
		}
	}
	return UINT32_MAX;
}


void WriteBrickIndexFromCacheToDisk(int set, int cacheIndex)
{
	int tlidx = IndexCache[set][cacheIndex].idx;
	uint brickIdx = IndexCache[set][cacheIndex].res;
	// write the brick index for the specified location to the index file
	fseek(indexFile, tlidx * 4, SEEK_SET);
	fwrite(&brickIdx, 1, 4, indexFile);
}

// SetBrickIndexInCache(place on disk in brickindex, place on disk in brickfiles)
// call only when placing clean data in cache
void SetBrickIndexInCache(int index, uint result)
{
	int set = index % (IndexSetSize);

	for (int i = 0; i < Index4Bytes / IndexSetSize; i++)
	{
		if (IndexCache[set][i].valid && IndexCache[set][i].idx == index)
		{
			return;
		}
	}

	int rando = EvictFromIndexCache();
	if (IndexCache[set][rando].dirty)
		WriteBrickIndexFromCacheToDisk(set, rando);
	IndexCache[set][rando].idx = index;
	IndexCache[set][rando].res = result;
	IndexCache[set][rando].valid = true;
	IndexCache[set][rando].dirty = false;
}

void WriteBrickIndexInCache(int index, uint result)
{
	cacheChecks++;
	int set = index % (IndexSetSize);
	
	for (int i = 0; i < IndexSetSize; i++)
	{
		if (IndexCache[set][i].valid && IndexCache[set][i].idx == index)
		{
			IndexCache[set][i].res = result;
			IndexCache[set][i].dirty = true;
			cacheHits++;
			return;
		}
	}

	int rando = EvictFromIndexCache();

	if (IndexCache[set][rando].dirty)
		WriteBrickIndexFromCacheToDisk(set, rando);
	IndexCache[set][rando].idx = index;
	IndexCache[set][rando].res = result;
	IndexCache[set][rando].dirty = true;
	IndexCache[set][rando].valid = true;
}

void WriteBrickFromCacheToDisk(int set, int cacheIndex)
{
	uint index = BrickCache[set][cacheIndex].idx;
	uchar brickData[512];
	for (int i = 0; i < 512; i++)
		brickData[i] = BrickCache[set][cacheIndex].bricks[i];

	int dstRegionIdx = index / 131072;
	char dstBinFile[128];
	sprintf(dstBinFile, "assets/block%03i.bin", dstRegionIdx);
	FILE* d = fopen(dstBinFile, "r+b" /* open for updating */);
	fseek(d, (index& 131071) * 512, SEEK_SET);	
	fwrite(brickData, 1, 512, d);
	fclose(d);
}

int GetBrickFromCache(uchar buffer[512], uint index)
{
	brickChecks++;
	int set = index % (BrickSetSize);

	auto cacheSet = BrickCache[set];

	for (int i = 0; i < Brick4Bytes / BrickSetSize; i++)
	{
		if (cacheSet[i].valid && cacheSet[i].idx == index)
		{
			for (int j = 0; j < 512; j++)
			{
				buffer[j] = cacheSet[i].bricks[j];
			}
			brickHits++;
			return 1;
		}
	}

	return -1;
}

// call only when writing clean data to cache
void SetBrickInCache(uchar brick[512], uint index)
{
	int set = index % (BrickSetSize);

	for (int i = 0; i < Brick4Bytes / BrickSetSize; i++)
	{
		if (BrickCache[set][i].valid && BrickCache[set][i].idx == index)
			return;
	}

	int rando = EvictFromBrickCache();

	if (BrickCache[set][rando].dirty)
		WriteBrickFromCacheToDisk(set, rando);

	BrickCache[set][rando].idx = index;
	for (int i = 0; i < 512; i++)
	{
		BrickCache[set][rando].bricks[i] = brick[i];
	}
	BrickCache[set][rando].valid = true;
	BrickCache[set][rando].dirty = false;
}

void WriteBrickToCache(uchar brick[512], uint index)
{
	brickChecks++;
	int set = index % (BrickSetSize);

	for (int i = 0; i < Brick4Bytes / BrickSetSize; i++)
	{
		if ( BrickCache[set][i].idx == index)
		{
			for (int j = 0; j < 512; j++)
			{
				BrickCache[set][i].bricks[j] = brick[j];
			}
			BrickCache[set][i].dirty = true;
			BrickCache[set][i].valid = true;
			brickHits++;
			return;
		}
	}

	int rando = EvictFromBrickCache();

	if (BrickCache[set][rando].dirty)
		WriteBrickFromCacheToDisk(set, rando);

	
	for (int i = 0; i < 512; i++)
	{
		BrickCache[set][rando].bricks[i] = brick[i];
	}
	BrickCache[set][rando].idx = index;
	BrickCache[set][rando].dirty = true;
	BrickCache[set][rando].valid = true;
	
}

void FlushCaches()
{
	for(int i = 0; i < IndexSetSize; i++)
		for (int j = 0; j < Index4Bytes / IndexSetSize; j++)
		{
			if (IndexCache[i][j].valid && IndexCache[i][j].dirty)
				WriteBrickIndexFromCacheToDisk(i, j);
		}
	for (int i = 0; i < BrickSetSize; i++)
		for (int j = 0; j < Brick4Bytes / BrickSetSize; j++)
		{
			if (BrickCache[i][j].valid && BrickCache[i][j].dirty)
				WriteBrickFromCacheToDisk(i, j);
		}
}

// =================================================================
// VOXEL DATA INTERACTION
// The five functions below interact with the raw voxel data.
// Profiling indicates that disk i/o completely dominates runtime.
// Fix this by preventing disk i/o: use a cache to meet the vast
// majority of requests from a small amount of RAM.
// Note that the code that requests data (OptimizeWorld) never
// touches the same data twice. It does however exhibit a rather
// predictable access pattern. Exploit the fact that reading
// 4 bytes from a file takes just as long as e.g. 64 bytes. This is
// conceptually very similar to how a CPU cache operates.
// =================================================================

uint GetBrickIndex( int x, int y = 0, int z = 0 )
{
	// calculate position of voxel in top level grid
	int tlidx = x + z * 1024 + y * 1024 * 1024;
	uint brickIdx[PrefetchIndex];
	brickIdx[0] = GetBrickIndexFromCache(tlidx);
	if (brickIdx[0] == UINT32_MAX)
	{
		// read the brick index for the specified voxel from the index file
		fseek(indexFile, tlidx * 4, SEEK_SET);
		int maxRead = min(4 * PrefetchIndex, 256 * 131072 - 1);
		fread(&brickIdx, 1, maxRead, indexFile);
		for(int i = 0; i < PrefetchIndex; i++)
			SetBrickIndexInCache(tlidx + i, brickIdx[i]);
	}
	return brickIdx[0];
}

void SetBrickIndex( int x, int y, int z, uint brickIdx )
{
	// calculate position of voxel in top level grid
	int tlidx = x + z * 1024 + y * 1024 * 1024;

	WriteBrickIndexInCache(tlidx, brickIdx);	
}

void CopyBrick( uint dst, uint src )
{
	int srcRegionIdx = src / 131072;
	int dstRegionIdx = dst / 131072;

	uchar dstBrick[512];
	uchar srcBrick[512];
	
	if (GetBrickFromCache(srcBrick, src) == -1)
	{
		uchar prefetchedBricks[512 * PrefetchBricks];
		char srcBinFile[128];
		sprintf(srcBinFile, "assets/block%03i.bin", srcRegionIdx);
		FILE* s = fopen(srcBinFile, "rb"); // TODO: will this work if they are the same? Probably yes.

		//int seekStart = max( (src % 131072) - PrefetchBricks, (uint)0 ) * 512;

		int seekStart = max(0, (int) (src % 131072) - (PrefetchBricks - 1) ) * 512;

		fseek(s, seekStart, SEEK_SET);

		//int maxRead = min(PrefetchBricks, (int)(131072 - (src % 131072))) * 512;
		int maxRead = PrefetchBricks * 512;
		fread(prefetchedBricks, 1, maxRead, s);
		fclose(s);

		uchar brickybrick[512];
		int bricksRead = min(PrefetchBricks, (int)(131072 - (src % 131072) )) - 1;
		for (int i = 0; i < PrefetchBricks; i++)
		{
			if ((src & 131071) - bricksRead + i > 0 && (src % 131072) - bricksRead + i < 131072)
			{
				for (int j = 0; j < 512; j++)
					brickybrick[j] = prefetchedBricks[i * 512 + j];
				
				SetBrickInCache(brickybrick, (src - bricksRead) + i);
			}
		}
		for (int i = 0; i < 512; i++)
		{
			srcBrick[i] = prefetchedBricks[(bricksRead * 512) + i];
		}

	}
	WriteBrickToCache(srcBrick, dst);
}

uchar ReadVoxel( int x, int y, int z )
{
	static uchar brick[512]; // thread-safety down the drain
	// read the brick index from the index file
	uint brickIdx = GetBrickIndex( x / 8, y / 8, z / 8 );
	// handle the special 'empty' brick
	if (brickIdx == 0) return 0;
	// find the block file that contains the specified brick
	int blockIdx = brickIdx / 131072;
	if (GetBrickFromCache(brick, brickIdx) == -1)
	{
		uchar prefetchedBricks[512 * PrefetchBricks];
		char blockFileName[128];
		sprintf(blockFileName, "assets/block%03i.bin", blockIdx);
		// find the brick in the block file
		FILE* r = fopen(blockFileName, "rb");
		if (!r) return 0;
		fseek(r, (brickIdx & 131071) * 512, SEEK_SET);
		int maxRead = min(PrefetchBricks, (int)(  131072 - (brickIdx % 131072))) * 512;
		fread(prefetchedBricks, 1, maxRead, r);
		fclose(r);
		// find the specified voxel in the brick
		uchar brickybrick[512];
		for (int i = 0; i < PrefetchBricks; i++)
		{
			
			if ((brickIdx & 131071) + i < 131072)
			{
				for (int j = 0; j < 512; j++)
					brickybrick[j] = prefetchedBricks[i * 512 + j];
				SetBrickInCache(brickybrick, brickIdx + i);
			}
		}
		for (int i = 0; i < 512; i++)
		{
			brick[i] = prefetchedBricks[i];
		}
	}
	int localx = x & 7, localy = y & 7, localz = z & 7;
	return brick[localx + localy * 8 + localz * 8 * 8];
}

// todo
void WriteVoxel( int x, int y, int z, uchar v )
{
	// read the brick index from the index file
	uint brickIdx = GetBrickIndex( x / 8, y / 8, z / 8 );
	// handle the special 'empty' brick
	if (brickIdx == 0)
	{
		// brick is empty; claim a fresh one at the end
		brickIdx = brickCount++;
		// update the new brick index in the index file
		SetBrickIndex( x / 8, y / 8, z / 8, brickIdx );
		// we can now use this new brick index.
	}
	// find the block file that contains the specified brick
	uchar brick[512];
	int localx = x & 7, localy = y & 7, localz = z & 7;
	int posInBrick = localx + localy * 8 + localz * 8 * 8;
	if (GetBrickFromCache(brick, brickIdx) == -1)
	{
		int blockIdx = brickIdx / 131072;
		char blockFileName[128];
		sprintf(blockFileName, "assets/block%03i.bin", blockIdx);
		// write the voxel in the brick in the block file
		
		FILE* r = fopen(blockFileName, "r+b");
		if (!r)
		{
			// we have too few block files available; create a new one
			r = fopen(blockFileName, "wb");
			uchar emptyBrick[512];
			memset(emptyBrick, 0, 512);
			for (int i = 0; i < 131072; i++) fwrite(emptyBrick, 1, 512, r);
			fclose(r);
			blockFileCount++;
			// open the newly created file
			r = fopen(blockFileName, "r+b");
		}

		uchar prefetchedBricks[512 * PrefetchBricks];

		fseek(r, (brickIdx & 131071) * 512, SEEK_SET);
		int maxRead = min(PrefetchBricks, (int)(131072 - (brickIdx % 131072))) * 512;
		fread(prefetchedBricks, 1, maxRead, r);
		fclose(r);
		uchar brickybrick[512];
		for (int i = 0; i < PrefetchBricks; i++)
		{

			if ((brickIdx & 131071) + i < 131072)
			{
				for (int j = 0; j < 512; j++)
					brickybrick[j] = prefetchedBricks[i * 512 + j];
				SetBrickInCache(brickybrick, brickIdx + i);
			}
		}
		for (int i = 0; i < 512; i++)
		{
			brick[i] = prefetchedBricks[i];
		}
	}
	brick[posInBrick] = v;
	WriteBrickToCache(brick, brickIdx);
	
}

// =================================================================
// WORLD OPTIMIZATION
// Note: consider this a black box; this code is not an optimization
// target. FYI: optimization exploits the fact that the raw world
// only contains unique bricks. Once an empty brick is detected, it
// is replaced by '0'. The 8x8x8 brick data thus becomes unused, but
// this is useless, since it is in the middle of the set. Therefore,
// the now unoccupied data is replaced by the voxels of the last
// brick. The top-level grid cell that points to this brick is then
// redirected to the new copy. Once we have done this 131072 times,
// a full 'block255.bin' file can be deleted.
// =================================================================

// GetBrickCoordinatesForIdx:
// Returns the location in the world of a brick with the specified
// index. Only valid for bricks that are still in their original
// block file. After optimization this no longer works.
// -----------------------------------------------------------------
void GetBrickCoordinatesForIdx( const int idx, int& bx, int& by, int& bz )
{
	// get region coordinates (512x256x512)
	int region = idx / 131072;
	int rx = region & 15, rz = region / 16;
	// get coordinates inside region
	int lx = idx & 63;
	int lz = (idx >> 6) & 63;
	int ly = (idx >> 12) & 31;
	// finalize coordinates
	bx = rx * 64 + lx, by = ly, bz = rz * 64 + lz;
}

// BrickIsEmpty:
// Loads the brick with the specified index from a block file and 
// tests whether it is empty.
// -----------------------------------------------------------------
bool BrickIsEmpty( int brickIdx )
{
	uchar brick[512];

	if (GetBrickFromCache(brick, brickIdx) == -1)
	{
		uchar prefetchedBricks[512 * PrefetchBricks];
		int srcRegionIdx = brickIdx / 131072;

		if (srcRegionIdx > 0)
			int w = 5;

		char srcBinFile[128];
		sprintf(srcBinFile, "assets/block%03i.bin", srcRegionIdx);
		FILE* s = fopen(srcBinFile, "rb");
		fseek(s, (brickIdx & 131071) * 512, SEEK_SET);
		int maxRead = min(PrefetchBricks, (int)(131072 - (brickIdx % 131072))) * 512;
		fread(prefetchedBricks, 1, maxRead, s);
		fclose(s);
		
		uchar brickybrick[512];

		for (int i = 0; i < PrefetchBricks; i++)
		{
			
			if ((brickIdx & 131071) + i < 131072)
			{
				for (int j = 0; j < 512; j++)
					brickybrick[j] = prefetchedBricks[i * 512 + j];
				SetBrickInCache(brickybrick, brickIdx + i);
			}
		}

		for (int i = 0; i < 512; i++)
		{
			brick[i] = prefetchedBricks[i];
		}
	}
	for (int i = 0; i < 512; i++) if (brick[i]) return false;
	return true;
}

// OptimizeWorld:
// World optimization algorithm. Has two aims:
// 1. In the top-level grid (index.bin), every reference to an empty 
//    8x8x8 brick is replaced by 0 to optimize subsequent queries.
// 2. The set of block files (blockXXX.bin) is reduced by
//    'defragmentation'. Each time 131072 bricks have been cleared, 
//    one block file is deleted.
// You can interrupt the process at any time by pressing ESCAPE.
// Restarting the application continues the optimization process
// where it left of.
// -----------------------------------------------------------------
void OptimizeWorld()
{
	// optimize the full world, will take ages (well, 22 days on my system)
	int optimized = 0, countdown = 256;
	// scan the cells of the top-level grid for empty bricks
	Timer timer;
	timer.reset();
	for (int b = startCell; b < 1024 * 32 * 1024; b++)
	{
		// iterate through the 64x32x64 regions
		int region = b >> 17;
		int rx = region & 15, rz = region / 16;
		int posInRegion = b & 131071;
		int lx = posInRegion & 63, lz = (posInRegion >> 6) & 63, ly = posInRegion >> 12;
		// get the index of a brick from the cell
		uint brickIdx = GetBrickIndex( rx * 64 + lx, ly, rz * 64 + lz );
		if (brickIdx == 0) continue; // skip the cell if it was already cleared
		// scan the brick for 8x8x8 zeroes
		if (BrickIsEmpty( brickIdx ))
		{
			// this brick is empty; replace cell value by 0
			SetBrickIndex( rx * 64 + lx, ly, rz * 64 + lz, 0 );
			// halt when we are optimizing the last region
			if ((b >> 17) == (blockFileCount - 1)) break;
			// move the last brick from the last region file
			CopyBrick( brickIdx, --brickCount );
			optimized++;
			// find out where that brick came from
			int ox, oy, oz;
			GetBrickCoordinatesForIdx( brickCount, ox, oy, oz );
			SetBrickIndex( ox, oy, oz, brickIdx );
			deletedBricks++;
			// once we have removed 131072 bricks, we can kill one block file
			while (deletedBricks >= 131072)
			{
				printf( "killed a 64MB block file.\n" );
				deletedBricks -= 131072;
				// do the actual region file killing
				char regionFileName[128];
				sprintf( regionFileName, "assets/block%03i.bin", --blockFileCount );
				remove( regionFileName );
			}
			// decrease b; we may have pulled in an empty brick from the end
			b--;
		}
		// report
		if (--countdown == 0)
		{
			printf( "removed %i of 256 bricks at b=%i; deleted bricks: %i\n", optimized, b, deletedBricks );
			//cout << "removed " << optimized << " of 256 bricks at b=" << b << "; deleted bricks: " << deletedBricks << "\n";
			optimized = 0, countdown = 256;
		}
		// graceful exit when ESC is pressed
		if (GetAsyncKeyState( VK_ESCAPE ))
		{
			// save optimization state; we'll load this when we restart the program.
			FILE* f = fopen( "optimstate.bin", "wb" );
			fwrite( &blockFileCount, 1, 4, f );
			fwrite( &deletedBricks, 1, 4, f );
			fwrite( &brickCount, 1, 4, f );
			fwrite( &b, 1, 4, f );
			fclose( f );
			
			FlushCaches();
			fclose(indexFile);
			exit( 0 );
		}
		if ((b % 1024) == 0)
		{
			//cout << timer.elapsed() << "\n";
			//cout << "cache check " << cacheChecks << " cache hits " << cacheHits << " percent " << (float)cacheHits / float(cacheChecks) << "\n";
			//cout << "brick check " << brickChecks << " brick hits " << brickHits << " percent " << (float)brickHits / float(brickChecks) << "\n";
			cacheHits = 0;
			cacheChecks = 0;
			brickChecks = 0;
			brickHits = 0;
			timer.reset();
		}
	}
}

// MikadoWorld:
// Draws a set of random lines in the disk-based copy of the world.
// -----------------------------------------------------------------
void MikadoWorld()
{
	for( int i = 0; i < 100; i++ )
	{
		// place the line in a 512 x 256 x 512 box
	#if 1
		int bx = 7200;
		int bz = 7200;
	#else
		int bx = RandomUInt() % (8192 - 512);
		int bz = RandomUInt() % (8192 - 512);
	#endif
		// actual line start and end
		int x1 = bx + (RandomUInt() & 511), y1 = RandomUInt() & 255, z1 = bz + (RandomUInt() & 511);
		int x2 = bx + (RandomUInt() & 511), y2 = RandomUInt() & 255, z2 = bz + (RandomUInt() & 511);
		// draw line
		printf( "Line: (%i,%i,%i) - (%i,%i,%i)... ", x1, y1, z1, x2, y2, z2 );
		Timer t;
		float l = (float)max( max( abs( x2 - x1 ), abs( y2 - y1 ) ), abs( z2 - z1 ) );
		printf( "voxels: %i ", (int)l );
		float dx = (x2 - x1) / l, dy = (y2 - y1) / l, dz = (z2 - z1) / l;
		float x = (float)x1, y = (float)y1, z = (float)z1;
		for( int i = 0; i < l; i++, x += dx, y += dy, z += dz )
			WriteVoxel( (int)x, (int)y, (int)z, 7 << 5 );
		printf( "%1.2fms.\n", t.elapsed() * 1000 );
	}
}

// =================================================================
// TEMPLATE CODE
// The usual Init and Tick.
// =================================================================

// Init:
// Prepare the world, initializes the 3D camera, restore the state 
// of the optimization algorithm and starts the actual optimizaion.
// -----------------------------------------------------------------
void MCViewer::Init()
{
	//std::ofstream out("out.txt");
	//std::streambuf* coutbuf = std::cout.rdbuf(); //save old buf
	//std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
	indexFile = fopen("assets/index.bin", "r+b" /* open for updating */);
	// initialize the world
	ClearWorld();
	// position the camera
	GetWorld()->SetCameraMatrix( mat4::LookAt( make_float3( 600, 250, 800 ), make_float3( 452, 10, 212 ) ) );
	// restore optimization state if this is not our first rodeo.
	FILE* f = fopen( "optimstate.bin", "rb" );
	if (f) // otherwise we'll use the defaults.
	{
		fread( &blockFileCount, 1, 4, f );
		fread( &deletedBricks, 1, 4, f );
		fread( &brickCount, 1, 4, f );
		fread( &startCell, 1, 4, f );
		fclose( f );
	}
	// ENABLE ME FOR TESTING THE OPTIMIZATION FUNCTION:
	OptimizeWorld();
	
	// ENABLE ME FOR TESTING LINE DRAWING
	//MikadoWorld();
	FlushCaches();
}

// Tick:
// NOTE: We only ever get here if Init completes. Since init does
// the optimization, this may take a while (days). Disable
// OptimizeWorld to view a slice of the landscape. The slice that is
// displayed is intentionally at the end of the set; any error in
// the optimization process should therefore show up quickly.
// -----------------------------------------------------------------
void MCViewer::Tick( float deltaTime )
{
	// This function gets called once per frame by the template code.
	const int firstx = 128, lastx = 896; // big slice
	// visualize a slice of the landscape
	static int x = firstx, y = 0, z = 380 /* same as x */, zstep = 1;
	Timer timer;
	timer.reset();
	for (int i = 0; i < 1024; i++)
	{
		uchar v = ReadVoxel( x + 7292, y, z + 7292 );
		Plot( x, y, z, v );
		if (++y == 256)
		{
			y = 0; if (++x == lastx) x = firstx, z = (z + zstep) & 1023;
		}
	}
	cout << timer.elapsed() << "\n";
}