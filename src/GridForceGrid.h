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
//#include "ComputeGridForce.h"

#include "MGridforceParams.h"

class GridDepositMsg;
class GridforceSubGrid;

#define GRIDBLOCKX 32
#define GRIDBLOCKY 32
#define GRIDBLOCKZ 32

class GridforceGrid {
    friend class GridforceMainGrid;
    friend class GridforceSubGrid;
    
public:
    GridforceGrid(void);
    virtual ~GridforceGrid();
    
    int request_box(Vector pos);
//    int get_box(Box *box, Vector pos) const;
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    virtual void init1(char *potfilename, 
		       SimParameters *simParams,
		       MGridforceParams *mgridParams) = 0;
    virtual void init2();
    void init3(float *grid, int *start, int *count);
    virtual void init4(float *grid, int start, int count) = 0;
    
    inline Position get_center(void) const { return center; }
    inline Position get_origin(void) const { return origin; }
    inline Tensor get_e (void) const { return e; }
    inline Tensor get_inv(void) const { return inv; }
    inline Vector get_scale(void) const { return scale; }

    int compute_VdV(Position pos, float &V, Vector &dV) const;
    
    inline int get_k0(void) const { return k[0]; }
    inline int get_k1(void) const { return k[1]; }
    inline int get_k2(void) const { return k[2]; }
    
    inline void req_grid(const int myrank, const int i) {
      int bi, boffset;
      
      //CkPrintf("req_grid: i = %d myrank = %d\n", i, myrank);
      
      getGridBlockIdxOffset(i,bi,boffset);
      if (!isGridValAvail(bi,boffset)) {
	//CkPrintf("calling addGridProcRequest(%d, %d)\n", myrank, bi);
        addGridProcRequest(myrank,bi);
      }
      return; 
    };

    inline float get_grid(const int i) const {
      int bi, boffset;
      
      getGridBlockIdxOffset(i,bi,boffset);
      //CkPrintf("get_grid: i = %d bi = %d boffset = %d pointer = %d\n", i, bi, boffset, (int64)this);
      
//       if (grid_cache[bi] == 0) {
// 	  iout << iINFO << "ELEM[" << CkMyPe() << "][" << i
// 	       << "] not in grid cache "
// 	       << " gridnum=" << mygridnum
// 	       << " block=" << bi
// 	       << " boff=" << boffset
// 	       << " pointer=" << (int64)this
// 	       << "\n" << endi;
//       }
      //CkPrintf("grid_cache[bi][boffset] = %f\n", grid_cache[bi][boffset]);
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
      
//       iout << iINFO << "Idx " << idx 
// 	   << " ix=" << ix << "," << iy << "," << iz
// 	   << " bx=" << bx << "," << by << "," << bz
// 	   << " ox=" << ox << "," << oy << "," << oz
// 	   << " ox=" << ox << "," << oy << "," << oz
// 	   << " bxsz=" << bysz << "," << bzsz
// 	   << " block=" << block
// 	   << " offset=" << offset
// 	   << "\n" << endi;
      
      return;
    }
    

    inline int getBlockSize(int blk) const {
      //CkPrintf("Blk %d sz %d %d %d %d %d %d\n",blk,k[0],k[1],k[2],bxsz,bysz,bzsz);
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
      
//       if (sx * sy * sz != GRIDBLOCKX * GRIDBLOCKY * GRIDBLOCKZ) {
// 	  CkPrintf("Returning partial block %d (%d,%d,%d) = %d,%d,%d\n",
// 		   blk, bx,by,bz,sx,sy,sz);
//       }
      
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
      return blk % CkNumNodes();
    }

    
    inline Bool isGridValCached(int ind) const {
      return (getGridValHomeNode(getGridBlockIdx(ind)) == CkMyNode());
    }

    inline Bool isGridValAvail(int bi,int boffset) const {
	//CkPrintf("isGridValAvail(%d, %d) called\n",bi,boffset);
      return (grid_cache[bi] != 0);
    }

    void addGridValue(int index, float val) {
      int bi, boffset;
      getGridBlockIdxOffset(index,bi,boffset);
      
      if (grid_cache[bi] == 0)  {
// 	CkPrintf("[%d] adding index %d %d %d %f\n",
// 		   CkMyPe(),index,bi,boffset,val);
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
      //CkPrintf("addGridProcRequest(%d, %d)\n", rank, bi);
      grid_proc_request[rank][bi] = true;
    };
    
    inline int getNumGridRequests() const { 
      //CkPrintf("num_requests = %d\n", num_requests);
      return num_requests;
    };
    
protected:
    struct GridIndices {
      int inds2;
      int dk_hi;
      int dk_lo;
      Bool zero_derivs;
    };
   
    void precomputeIndexing();
    
    // Utility functions
    void readHeader(SimParameters *simParams, MGridforceParams *mgridParams);
    
    //virtual int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const = 0;
    int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const;
    void compute_a(float *a, float *b) const;
    virtual void compute_b(float *b, int *inds, Vector gapscale)  const = 0;
    float compute_V(float *a, float *x, float *y, float *z) const;
    Vector compute_dV(float *a, float *x, float *y, float *z) const;
    Vector compute_d2V(float *a, float *x, float *y, float *z) const;
    float compute_d3V(float *a, float *x, float *y, float *z) const;
    
    void addDepositMsg(GridDepositMsg **outmsgs, int &arrayIdx, int startIdx, int grandTotalGrids);
    void readSubgridHierarchy(FILE *poten, int &totalGrids);
    
    FILE *poten_fp;
    //    float *grid;	// Actual grid
    
    GridforceSubGrid **subgrids;
    int numSubgrids;
    int generation;	// Subgrid level (0 = main grid)
    
    int k[3];		// Grid dimensions
    int k_nopad[3];	// Grid dimensions
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
    
    float **grid_cache;
    Bool *grid_request;
    Bool **grid_proc_request;
    int num_requests;
    
    Bool *gridvals_needed;
    int grid_elems_read;
    GridIndices *grid_index_table[8][3];
    int *block_index[3];
    int *block_offset[3];
    int bxsz, bysz, bzsz;  // Dimensions of grid in grid blocks
};


class GridforceMainGrid : public GridforceGrid {
    friend class GridforceGrid;
    friend class GridforceSubGrid;

public:
    GridforceMainGrid(int gridnum);
    
    void init1(char *potfilename,
	       SimParameters *simParams,
	       MGridforceParams *mgridParams);
    void init2();
    void init4(float *grid, int start, int count);
    
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    inline int getTotalGrids(void) const { return totalGrids; }
    int buildDepositMsgs(GridDepositMsg **outmsgs, int startIdx, int grandTotalGrids);
    
protected:
    //int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const;
    void compute_b(float *b, int *inds, Vector gapscale)  const;
    void buildSubgridsFlat(void);
    
    int totalGrids;
    GridforceSubGrid **subgrids_flat;
    int mygridnum;
    
    static const int border = 1;
};


class GridforceSubGrid : public GridforceGrid {
    friend class GridforceGrid;
    friend class GridforceMainGrid;

public:
    GridforceSubGrid(GridforceGrid *parent_in);
    
    void init1(char *potfilename,
	       SimParameters *simParams,
	       MGridforceParams *mgridParams);
    void init2();
    void init4(float *grid, int start, int count);
    
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    inline Tensor tensorMult (const Tensor &t1, const Tensor &t2) {
	Tensor tmp;
	tmp.xx = t1.xx * t2.xx + t1.xy * t2.yx + t1.xz * t2.zx;
	tmp.xy = t1.xx * t2.xy + t1.xy * t2.yy + t1.xz * t2.zy;
	tmp.xz = t1.xx * t2.xz + t1.xy * t2.yz + t1.xz * t2.zz;
	tmp.yx = t1.yx * t2.xx + t1.yy * t2.yx + t1.yz * t2.zx;
	tmp.yy = t1.yx * t2.xy + t1.yy * t2.yy + t1.yz * t2.zy;
	tmp.yz = t1.yx * t2.xz + t1.yy * t2.yz + t1.yz * t2.zz;
	tmp.zx = t1.zx * t2.xx + t1.zy * t2.yx + t1.zz * t2.zx;
	tmp.zy = t1.zx * t2.xy + t1.zy * t2.yy + t1.zz * t2.zy;
	tmp.zz = t1.zx * t2.xz + t1.zy * t2.yz + t1.zz * t2.zz;
	return tmp;
    }
    
protected:
    //int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const;
    void compute_b(float *b, int *inds, Vector gapscale) const;
    void addToSubgridsFlat(void);
    
    // Utility numbers
    Tensor scale_dV;
    Tensor scale_d2V;
    float scale_d3V;
    
    GridforceGrid *parent;
    int pmin[3], pmax[3];
    GridforceMainGrid *maingrid;
    int subgridIdx;
};

#endif
