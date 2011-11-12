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

class GridforceSubGrid;

/*
#define GRIDBLOCKX 32
#define GRIDBLOCKY 32
#define GRIDBLOCKZ 32
*/

class GridforceGrid {
    friend class GridforceMainGrid;
    friend class GridforceSubGrid;
    
public:
    GridforceGrid(void);
    virtual ~GridforceGrid();
    
//    int request_box(Vector pos);
//    int get_box(Box *box, Vector pos) const;
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
  
    inline Position get_center(void) const { return center; }
    inline Position get_origin(void) const { return origin; }
    inline Tensor get_e (void) const { return e; }
    inline Tensor get_inv(void) const { return inv; }
    inline Vector get_scale(void) const { return scale; }
    
    float get_grid(int i0, int i1, int i2) const;
    void set_grid(int i0, int i1, int i2, float V);
    
    int compute_VdV(Position pos, float &V, Vector &dV) const;
    
    inline int get_k0(void) const { return k[0]; }
    inline int get_k1(void) const { return k[1]; }
    inline int get_k2(void) const { return k[2]; }
    
protected:
    struct GridIndices {
      int inds2;
      int dk_hi;
      int dk_lo;
      Bool zero_derivs;
    };
   
    // Utility functions
    void readHeader(SimParameters *simParams, MGridforceParams *mgridParams);
    
    inline int grid_index(int i0, int i1, int i2) const
    {
	return i0*dk[0] + i1*dk[1] + i2*dk[2];
    }
    
    //virtual int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const = 0;
    int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const;
    void compute_a(float *a, float *b) const;
    virtual void compute_b(float *b, int *inds, Vector gapscale)  const = 0;
    float compute_V(float *a, float *x, float *y, float *z) const;
    Vector compute_dV(float *a, float *x, float *y, float *z) const;
    Vector compute_d2V(float *a, float *x, float *y, float *z) const;
    float compute_d3V(float *a, float *x, float *y, float *z) const;
    
    void readSubgridHierarchy(FILE *poten, int &totalGrids);
    
    FILE *poten_fp;
    float *grid;	// Actual grid
    
    GridforceSubGrid **subgrids;
    int numSubgrids;
    int generation;	// Subgrid level (0 = main grid)
    
    // should move 'nopad' versions to maingrid only ... or not have them as ivars at all, why are they here?
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
};


class GridforceMainGrid : public GridforceGrid {
    friend class GridforceGrid;
    friend class GridforceSubGrid;

public:
    GridforceMainGrid(int gridnum);
    virtual ~GridforceMainGrid();
    
    void initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams);
    void reinitialize(SimParameters *simParams, MGridforceParams *mgridParams);
    
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    int get_all_gridvals(float** all_gridvals) const;
    void set_all_gridvals(float* all_gridvals, int sz);
    
    inline int getTotalGrids(void) const { return totalGrids; }
    
protected:
    //int get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const;
    void compute_b(float *b, int *inds, Vector gapscale)  const;
    void buildSubgridsFlat(void);
    
    char filename[129];
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
    
    void initialize(SimParameters *simParams, MGridforceParams *mgridParams);
    
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
