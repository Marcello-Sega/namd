/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GRIDFORCEGRID_H
#define GRIDFORCEGRID_H

#include "Vector.h"
#include "Tensor.h"
#include "SimParameters.h"
#include "NamdTypes.h"
#include "MStream.h"


class GridforceGrid {
public:
    GridforceGrid(void);
    //GridforceGrid(GridforceGrid g);
    ~GridforceGrid();
    
    struct __box
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
    typedef struct __box Box;
    
    void initialize(char *potfilename, SimParameters *simParams);
    int get_box(Box *box, Vector pos) const;
    void pack(MOStream *msg) const;
    void unpack(MIStream *msg);
    
    inline Position get_center(void) const { return center; }
    inline Tensor get_inv(void) const { return inv; }
    
private:
    float *grid;	// Actual grid
    
    int k[3];		// Grid dimensions
    int size;
    int dk[3];
    
    Position origin;	// Grid origin
    Position center;	// Center of grid (for wrapping)
    Tensor e;		// Grid unit vectors
    Tensor inv;		// Inverse of unit vectors
    
    float pad_p[3];	// Pad values (p = positive side, n = negative side) for each dimension
    float pad_n[3];
    Bool cont[3];	// Whether grid is continuous in each dimension
    float offset[3];	// Potential offset in each dimension
    float gap[3];	// Gap between images of grid in grid units for each dimension
    float gapinv[3];	// 1.0/gap
    
    enum {border = 3}; // Size of border (for generality)
};


#endif
