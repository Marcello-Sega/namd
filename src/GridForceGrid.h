/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GRIDFORCEGRID_H
#define GRIDFORCEGRID_H

#include <set>

#include "Vector.h"
#include "Tensor.h"
#include "SimParameters.h"
#include "NamdTypes.h"
#include "MStream.h"
#include "charm++.h"

#include "MGridforceParams.h"

#define GRIDBLOCKX 32
#define GRIDBLOCKY 32
#define GRIDBLOCKZ 32

class GridforceGrid {
private:
    
public:
    GridforceGrid(int gridnum);
    //GridforceGrid(GridforceGrid g);
    ~GridforceGrid();
    
    struct _box
    {
	Vector loc;	// Position within grid
	float b[64];	// b[0 .. 7] = V (eight corners)
			// b[8 ..15] = dV/dx
			// b[16..23] = dV/dy
			// b[24..31] = dV/dz
			// b[32..39] = d2V/dxdy
			// b[40..47] = d2V/dxdz
			// b[48..55] = d2V/dydz
			// b[56..63] = d3V/dxdydz
	Tensor scale;
    };
    typedef struct _box Box;
    
//    void initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams);
    int request_box(Vector pos);
    int get_box(Box *box, Vector pos) const;
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    void init1(char *potfilename, 
               SimParameters *simParams,
               MGridforceParams *mgridParams);
    void init2();
    void init3(float *grid, int *start, int *count);
    void init4(float *grid, int start, int count);
    inline Position get_center(void) const { return center; }
    inline Position get_origin(void) const { return origin; }
    inline Tensor get_e (void) const { return e; }
    inline Tensor get_inv(void) const { return inv; }
    inline Vector get_scale(void) const { return scale; }

    inline void req_grid(const int myrank, const int i) {
      int bi, boffset;
      
      getGridBlockIdxOffset(i,bi,boffset);
      if (!isGridValAvail(bi,boffset)) {
        addGridProcRequest(myrank,bi);
      }
      return; 
    };

    inline float get_grid(const int i) const {
      int bi, boffset;
      
      getGridBlockIdxOffset(i,bi,boffset);
      
//      if (grid_cache[bi] == 0) {
//          iout << iINFO << "ELEM[" << CkMyPe() << "][" << i
//             << "] not in grid cache "
//             << " gridnum=" << mygridnum
//             << " block=" << bi
//             << " boff" << boffset
//             << "\n" << endi;
//      }
      return grid_cache[bi][boffset];
    }
    
    void allocateGridCache();
    void allocateGridEntry(int idx);
    
    inline int getGridBlockIdx(const int idx) const {
      const int iz = idx % k[2];
      const int iy = (idx / k[2]) % k[1];
      const int ix = idx / (k[1] * k[2]);
      
      const int bx = block_index[0][ix];
      const int by = block_index[1][iy];
      const int bz = block_index[2][iz];
      
      const int bi = (bx * bysz + by) * bzsz + bz;
      
//      iout << iINFO << "Idx " << idx 
//           << " ix=" << ix << "," << iy << "," << iz
//           << " bx=" << bx << "," << by << "," << bz
//           << " bi=" << bi
//           << "\n" << endi;
     
      return bi;
    }

    inline int getGridBlockOffset(const int idx) const {
      int bl, off;
      getGridBlockIdxOffset(idx,bl,off);
      return off;
    }

    inline void getGridBlockIdxOffset(const int idx, 
                 int &block, int &offset) const 
    {
      const int iz = idx % k[2];
      const int iy = (idx / k[2]) % k[1];
      const int ix = idx / (k[1] * k[2]);
      
      const int bx = block_index[0][ix];
      const int by = block_index[1][iy];
      const int bz = block_index[2][iz];
      
      const int ox = block_offset[0][ix];
      const int oy = block_offset[1][iy];
      const int oz = block_offset[2][iz];

      block = (bx * bysz + by) * bzsz + bz;
      offset = (ox * GRIDBLOCKY +  oy) * GRIDBLOCKZ + oz;
      
//      iout << iINFO << "Idx " << idx 
//           << " ix=" << ix << "," << iy << "," << iz
//           << " bx=" << bx << "," << by << "," << bz
//           << " ox=" << ox << "," << oy << "," << oz
//           << " ox=" << ox << "," << oy << "," << oz
//           << " bxsz=" << bysz << "," << bzsz
//           << " block=" << block
//           << " offset=" << offset
//           << "\n" << endi;
      
      return;
    }
    

    inline int getBlockSize(int blk) const {
//      CkPrintf("Blk %d sz %d %d %d %d %d %d\n",
//                blk,k[0],k[1],k[2],bxsz,bysz,bzsz);
      const int bz = blk % bzsz;
      const int by = (blk / bzsz) % bysz;
      const int bx = blk / (bysz * bzsz);
      int sx, sy, sz;
      
      if (bx < bxsz - 1)
        sx = GRIDBLOCKX;
      else {
        sx = k[0] % GRIDBLOCKX;
        if (sx == 0) sx = GRIDBLOCKX;
      }

      if (by < bysz - 1)
        sy = GRIDBLOCKY;
      else {
        sy = k[1] % GRIDBLOCKY;
        if (sy == 0) sy = GRIDBLOCKY;
      }

      if (bz < bzsz - 1)
        sz = GRIDBLOCKZ;
      else {
        sz = k[2] % GRIDBLOCKZ;
        if (sz == 0) sz = GRIDBLOCKZ;
      }
      
//      if (sx * sy * sz != GRIDBLOCKX * GRIDBLOCKY * GRIDBLOCKZ) {
//        CkPrintf("Returning partial block %d (%d,%d,%d) = %d,%d,%d\n",
//                  blk, bx,by,bz,sx,sy,sz);
//      }
      
      return sx * sy * sz;
    }

    inline int getIndexForBlock(int blk,int bidx) const {
      const int bz = blk % bzsz;
      const int by = (blk / bzsz) % bysz;
      const int bx = blk / (bysz * bzsz);
      int sx, sy, sz;
      
      if (bx < bxsz - 1)
        sx = GRIDBLOCKX;
      else {
        sx = k[0] % GRIDBLOCKX;
        if (sx == 0) sx = GRIDBLOCKX;
      }

      if (by < bysz - 1)
        sy = GRIDBLOCKY;
      else {
        sy = k[1] % GRIDBLOCKY;
        if (sy == 0) sy = GRIDBLOCKY;
      }

      if (bz < bzsz - 1)
        sz = GRIDBLOCKZ;
      else {
        sz = k[2] % GRIDBLOCKZ;
        if (sz == 0) sz = GRIDBLOCKZ;
      }
         
      const int biz = bidx % sz;
      const int biy = (bidx / sz) % sy;
      const int bix = bidx / (sy * sz);
      
      const int idx = ((bx * GRIDBLOCKX + bix) * k[1] 
                       + (by * GRIDBLOCKY + biy)) * k[2] 
                       + (bz * GRIDBLOCKZ + biz);
                       
//      CkPrintf("Returning index for block %d (%d,%d,%d) = %d (%d,%d,%d) %d\n",
//                  blk, bx,by,bz,bidx,bix,biy,biz,idx);
      return idx;
    }
    
    inline int getGridValHomeNode(int blk) const {
//      const int blksz = 17;
//      const int bk0 = (k[0]+blksz-1) / blksz;
//      const int bk1 = (k[1]+blksz-1) / blksz;
      
//      const int i0 = ind / dk[0];
//      const int i1 = (ind - i0*dk[0]) / dk[1];
//      const int i2 = (ind - i0*dk[0] - i1*dk[1]) / dk[2];
      
//      const int new_ind = ((i0 / blksz * bk0) + i1/blksz) * bk1 + i2/blksz;
      
      return blk % CkNumNodes();
    }

    
    inline Bool isGridValCached(int ind) const {
      return (getGridValHomeNode(getGridBlockIdx(ind)) == CkMyNode());
    }

    inline Bool isGridValAvail(int bi,int boffset) const {
//      CkPrintf("%d %d\n",bi,boffset);
      return (grid_cache[bi] != 0);
    }

    void addGridValue(int index, float val) {
      int bi, boffset;
      getGridBlockIdxOffset(index,bi,boffset);
      
      if (grid_cache[bi] == 0)  {
//        CkPrintf("[%d] adding index %d %d %d %f\n",
//                 CkMyPe(),index,bi,boffset,val);
        allocateGridEntry(bi);
      }
       
      grid_cache[bi][boffset] = val;
    }
    
    void allocateGridRequestTables();
    void consolidateGridRequests();
    static void getGridIndices(int *gridStartIndex, 
                               int *gridIndexList,
                               int *gridProcList,
                               int *gridProcCount);
    
    inline void addGridProcRequest(int rank, int bi) {
      grid_proc_request[rank][bi] = true;
    };
    
    inline int getNumGridRequests() const { 
      return num_requests;
    };
    
    inline int get_k0(void) const { return k[0]; }
    inline int get_k1(void) const { return k[1]; }
    inline int get_k2(void) const { return k[2]; }
    
private:
    struct GridIndices {
      int inds2;
      int dk_hi;
      int dk_lo;
      Bool zero_derivs;      
    };
   
    void precomputeIndexing();
    
    FILE *poten_fp;
    //    float *grid;	// Actual grid
    
    int k[3];		// Grid dimensions
    int k_nopad[3];		// Grid dimensions
    int size;
    int size_nopad;
    int dk[3];
    int dk_nopad[3];
    float factor;
    
    Position origin;	// Grid origin
    Position center;	// Center of grid (for wrapping)
    Tensor e;		// Grid unit vectors
    Tensor inv;		// Inverse of unit vectors
    
    float p_sum[3];     // Accumulators for sums
    float n_sum[3];
    float pad_p[3];	// Pad values (p = positive side, n = negative side) for each dimension
    float pad_n[3];
    Bool cont[3];	// Whether grid is continuous in each dimension
    float offset[3];	// Potential offset in each dimension
    float gap[3];	// Gap between images of grid in grid units for each dimension
    float gapinv[3];	// 1.0/gap

    Vector scale;
    
    static const int border = 1; // Size of border (for generality)
    float **grid_cache;
    Bool *grid_request;
    Bool **grid_proc_request;
    int num_requests;
    
    Bool *gridvals_needed;
    int grid_elems_read;
    int mygridnum;
    GridIndices *grid_index_table[8][3];
    int *block_index[3];
    int *block_offset[3];
    int bxsz, bysz, bzsz;  // Dimensions of grid in grid blocks
};


#endif
