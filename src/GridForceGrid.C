/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <iostream>

#include "GridForceGrid.h"
#include "Vector.h"
#include "SimParameters.h"
#include "InfoStream.h"
#include "common.h"

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"


GridforceGrid::GridforceGrid(void)
{
    grid = NULL;
}

// GridforceGrid::GridforceGrid(GridforceGrid g) :
//     size(g.size), origin(g.origin), center(g.center), e(g.e), inv(g.inv)
// {
//     for (int i = 0; i < 3; i++) {
// 	k[i] = g.k[i];
// 	dk[i] = g.dk[i];
// 	pad_p[i] = g.pad_p[i];
// 	pad_n[i] = g.pad_n[i];
// 	cont[i] = g.cont[i];
// 	offset[i] = g.offset[i];
// 	gap[i] = g.gap[i];
//     }
    
//     grid = new float[size];
//     memcpy(grid, g.grid, size*sizeof(float));
// }

GridforceGrid::~GridforceGrid() {
    delete [] grid;
}

void GridforceGrid::pack(MOStream *msg) const
{
    DebugM(2, "Packing message\n" << endi);
    
    msg->put(3*sizeof(int), (char*)k);
    msg->put(size);
    msg->put(3*sizeof(int), (char*)dk);
    
    msg->put(sizeof(Vector), (char*)&origin);
    msg->put(sizeof(Vector), (char*)&center);
    msg->put(sizeof(Tensor), (char*)&e);
    msg->put(sizeof(Tensor), (char*)&inv);
	     
    msg->put(3*sizeof(float), (char*)pad_p);
    msg->put(3*sizeof(float), (char*)pad_n);
    msg->put(3*sizeof(Bool), (char*)cont);
    msg->put(3*sizeof(float), (char*)offset);
    msg->put(3*sizeof(float), (char*)gap);
    
    msg->put(size*sizeof(float), (char*)grid);
}

void GridforceGrid::unpack(MIStream *msg)
{
    DebugM(2, "Unpacking message\n" << endi);
    
    msg->get(3*sizeof(int), (char*)k);
    msg->get(size);
    msg->get(3*sizeof(int), (char*)dk);
    
    msg->get(sizeof(Vector), (char*)&origin);
    msg->get(sizeof(Vector), (char*)&center);
    msg->get(sizeof(Tensor), (char*)&e);
    msg->get(sizeof(Tensor), (char*)&inv);
	     
    msg->get(3*sizeof(float), (char*)pad_p);
    msg->get(3*sizeof(float), (char*)pad_n);
    msg->get(3*sizeof(Bool), (char*)cont);
    msg->get(3*sizeof(float), (char*)offset);
    msg->get(3*sizeof(float), (char*)gap);
    
    if (size) {
	delete [] grid;
	grid = new float[size];
	msg->get(size*sizeof(float), (char*)grid);
    }
    
//     DebugM(2, "grid[0] = " << grid[0] << "\n");
//     DebugM(2, "k = " << k[0] << " " << k[1] << " " << k[2] << "\n");
//     DebugM(2, "size = " << size << "\n");
//     DebugM(2, "dk = " << dk[0] << " " << dk[1] << " " << dk[2] << "\n");
//     DebugM(2, "origin = " << origin << "\n");
//     DebugM(2, "center = " << center << "\n");
//     DebugM(2, "e = " << e << "\n");
//     DebugM(2, "inv = " << inv << "\n");
//     DebugM(2, "pad_p = " << pad_p[0] << " " << pad_p[1] << " " << pad_p[2] << "\n");
//     DebugM(2, "pad_n = " << pad_n[0] << " " << pad_n[1] << " " << pad_n[2] << "\n");
//     DebugM(2, "cont = " << cont[0] << " " << cont[1] << " " << cont[2] << "\n");
//     DebugM(2, "offset = " << offset[0] << " " << offset[1] << " " << offset[2] << "\n");
//     DebugM(2, "gap = " << gap[0] << " " << gap[1] << " " << gap[2] << endi << "\n");
}

void GridforceGrid::initialize(char *potfilename, SimParameters *simParams)
{
    FILE *poten = Fopen(potfilename, "r");
    if (!poten) {
	NAMD_die("Problem reading grid force potential file");
    }
    
    char line[256];
    do {
	fgets(line, 256, poten);	// Read comment lines
    } while (line[0] == '#');
    
    int k_nopad[3];
    sscanf(line, "object 1 class gridpositions counts %d %d %d\n", &k_nopad[0], &k_nopad[1], &k_nopad[2]);
    
    // Read origin
    fscanf(poten, "origin %lf %lf %lf\n", &origin.x, &origin.y, &origin.z);
        
    // Read delta (unit vectors)
    // These are column vectors, so must fill gridfrcE tensor to reflect this
    fscanf(poten, "delta %lf %lf %lf\n", &e.xx, &e.yx, &e.zx);
    fscanf(poten, "delta %lf %lf %lf\n", &e.xy, &e.yy, &e.zy);
    fscanf(poten, "delta %lf %lf %lf\n", &e.xz, &e.yz, &e.zz);
    
    center = origin + e * 0.5 * Position(k_nopad[0], k_nopad[1], k_nopad[2]);

    BigReal tmp[3];   // Temporary storage
    fscanf(poten, "object 2 class gridconnections counts %lf %lf %lf\n", tmp, tmp+1, tmp+2);
    fscanf(poten, "object 3 class array type double rank 0 items %lf data follows\n", tmp);
    
    // Calculate inverse tensor
    BigReal det;
    det = e.xx*(e.yy*e.zz - e.yz*e.zy) - e.xy*(e.yx*e.zz - e.yz*e.zx) + e.xz*(e.yx*e.zy - e.yy*e.zx);
    inv.xx =  (e.yy*e.zz - e.yz*e.zy)/det;
    inv.xy = -(e.xy*e.zz - e.xz*e.zy)/det;
    inv.xz =  (e.xy*e.yz - e.xz*e.yy)/det;
    inv.yx = -(e.yx*e.zz - e.yz*e.zx)/det;
    inv.yy =  (e.xx*e.zz - e.xz*e.zx)/det;
    inv.yz = -(e.xx*e.yz - e.xz*e.yx)/det;
    inv.zx =  (e.yx*e.zy - e.yy*e.zx)/det;
    inv.zy = -(e.xx*e.zy - e.xy*e.zx)/det;
    inv.zz =  (e.xx*e.yy - e.xy*e.yx)/det;
    
    // Allocate storage for potential
    int size_nopad = k_nopad[0] * k_nopad[1] * k_nopad[2];
    float *grid_nopad = new float[size_nopad];
    
    float factor;
    if (simParams->gridforceVolts) {
	factor = 1.0/0.0434;  // convert V -> kcal/mol*e
    } else {
	factor = 1.0;
    }
    
    float tmp2;
    for (int count = 0; count < size_nopad; count++) {
	int err = fscanf(poten, "%f", &tmp2);
	if (err == EOF || err == 0) {
	    NAMD_die("Grid force potential file incorrectly formatted");
	}
	grid_nopad[count] = tmp2 * factor;
    }
    
    DebugM(4, "e = " << e << "\n" << endi);
    DebugM(4, "inv = " << inv << "\n" << endi);
    
    int dk_nopad[3];
    
    // Shortcuts for accessing 1-D array with four indices
    dk_nopad[0] = k_nopad[1] * k_nopad[2];
    dk_nopad[1] = k_nopad[2];
    dk_nopad[2] = 1;
    
    // PARAMETERS which probably ought to be somewhere else
    BigReal modThresh = 1.0;
    
    Vector Kvec[3];
    Kvec[0] = e * Position(k_nopad[0]-1, 0, 0);
    Kvec[1] = e * Position(0, k_nopad[1]-1, 0);
    Kvec[2] = e * Position(0, 0, k_nopad[2]-1);
    Vector Avec[3];
    Avec[0] = simParams->lattice.a();
    Avec[1] = simParams->lattice.b();
    Avec[2] = simParams->lattice.c();
    
    // Decide whether we're wrapping
    for (int i0 = 0; i0 < 3; i0++) {
	if (simParams->gridforceCont[i0]) {
	    for (int i1 = 0; i1 < 3; i1++) {
		if (cross(Avec[i0], Kvec[i1]) == 0) {
		    cont[i1] = TRUE;
		    offset[i1] = simParams->gridforceVOffset[i0] * factor;
		    gap[i1] = (inv * (Avec[i0] - Kvec[i1])).length();	// want in grid-point units (normal = 1)
		    
		    if (gap[i1] < 0) {
			NAMD_die("Gridforce Grid overlap!");
		    }
		    
		    DebugM(3, "cont[" << i1 << "] = " << cont[i1] << "\n" << endi);
		    
		}
	    }
	    
	    if (!(cont[0] || cont[1] || cont[2])) {
		NAMD_die("No Gridforce unit vector found parallel to requested continuous grid direction!");
	    }
	}
	else {
	    // TO DO: must check for grid overlap in non-wrapping
	    // dimensions taking non-orthogonal unit vectors into
	    // consideration
	}
    }
    
    // Figure out size of true grid (padded on non-periodic sides)
    Vector delta = 0;
    for (int i = 0; i < 3; i++) {
	if (cont[i]) {
	    k[i] = k_nopad[i];
	} else {
	    k[i] = k_nopad[i] + 2*border;
	    delta[i] -= border;
	}
    }
    DebugM(3, "delta = " << e * delta << " (" << delta << ")\n");
    origin += e * delta;
    
    size = k[0] * k[1] * k[2];
    dk[0] = k[1] * k[2];
    dk[1] = k[2];
    dk[2] = 1;
    
    grid = new float[size];
    
    // Fill in most of new grid (all except pad areas)
    BigReal n_sum[3] = {0, 0, 0};
    BigReal p_sum[3] = {0, 0, 0};
    for (int i0 = 0; i0 < k_nopad[0]; i0++) {
	for (int i1 = 0; i1 < k_nopad[1]; i1++) {
	    for (int i2 = 0; i2 < k_nopad[2]; i2++) {
		// Edges are special cases -- take force there to be
		// zero for smooth transition across potential
		// boundary
		
		int ind_nopad = i0*dk_nopad[0] + i1*dk_nopad[1] + i2*dk_nopad[2];
		int j0 = (cont[0]) ? i0 : i0 + border;
		int j1 = (cont[1]) ? i1 : i1 + border;
		int j2 = (cont[2]) ? i2 : i2 + border;
		int ind = j0*dk[0] + j1*dk[1] + j2*dk[2];
		
		if (i0 == 0)			n_sum[0] += grid_nopad[ind_nopad];
		else if (i0 == k_nopad[0]-1)	p_sum[0] += grid_nopad[ind_nopad];
		if (i1 == 0)			n_sum[1] += grid_nopad[ind_nopad];
		else if (i1 == k_nopad[1]-1)	p_sum[1] += grid_nopad[ind_nopad];
		if (i2 == 0)			n_sum[2] += grid_nopad[ind_nopad];
		else if (i2 == k_nopad[2]-1)	p_sum[2] += grid_nopad[ind_nopad];
		
		grid[ind] = grid_nopad[ind_nopad];
	    }
	}
    }
    
    
    BigReal n_avg[3], p_avg[3];
    for (int i0 = 0; i0 < 3; i0++) {
	int i1 = (i0 + 1) % 3;
	int i2 = (i0 + 2) % 3;
	n_avg[i0] = n_sum[i0] / (k_nopad[i1] * k_nopad[i2]);
	p_avg[i0] = p_sum[i0] / (k_nopad[i1] * k_nopad[i2]);
	
	if (cont[i0] && fabs(offset[i0] - (p_avg[i0]-n_avg[i0])) > modThresh) {
	    iout << iWARN << "GRID FORCE POTENTIAL DIFFERENCE IN K" << i0
		 << " DIRECTION IS " << offset[i0] - (p_avg[i0]-n_avg[i0]) << " KCAL/MOL*E\n" << endi;
	}
	DebugM(3, "n_avg[" << i0 << "] = " << n_avg[i0] << "\n" << endi);
	DebugM(3, "p_avg[" << i0 << "] = " << p_avg[i0] << "\n" << endi);
    }
    
    
    Bool twoPadVals = (cont[0] + cont[1] + cont[2] == 2);
    float padVal = 0.0;
    int weight = 0;
    if (!twoPadVals) {
	// Determine pad value (must average)
	if (!cont[0]) {
	    padVal += p_sum[0] + n_sum[0];
	    weight += 2 * k_nopad[1] * k_nopad[2];
	}
	if (!cont[1]) {
	    padVal += p_sum[1] + n_sum[1];
	    weight += 2 * k_nopad[0] * k_nopad[2];
	}
	if (!cont[2]) {
	    padVal += p_sum[2] + n_sum[2];
	    weight += 2 * k_nopad[0] * k_nopad[1];
	}
	padVal /= weight;
    }
    
    for (int i = 0; i < 3; i++) {
	pad_n[i] = (cont[i]) ? 0.0 : (twoPadVals) ? n_avg[i] : padVal;
	pad_p[i] = (cont[i]) ? 0.0 : (twoPadVals) ? p_avg[i] : padVal;
    }
    
    if (cont[0] && cont[1] && cont[2]) {
	// Nothing to do
	return;
    }
    
    // Now fill in rest of new grid
    for (int i0 = 0; i0 < k[0]; i0++) {
	for (int i1 = 0; i1 < k[1]; i1++) {
	    for (int i2 = 0; i2 < k[2]; i2++) {
		if ( (cont[0] || (i0 >= border && i0 < k[0]-border)) &&
		     (cont[1] || (i1 >= border && i1 < k[1]-border)) &&
		     (cont[2] || i2 == border) )
		{
		    i2 += k_nopad[2]-1;
		    continue;
		}
		
		int ind = i0*dk[0] + i1*dk[1] + i2*dk[2];
		Position pos = e * Position(i0, i1, i2);
		int var[3] = {i0, i1, i2};
		
		for (int dir = 0; dir < 3; dir++) {
		    if (cont[dir]) continue;
		    
		    if (var[dir] < border)
			grid[ind] = pad_n[dir];
		    else if (var[dir] >= k[dir]-border)
			grid[ind] = pad_p[dir];
		}
		
		DebugM(2, "grid[" << ind << "; " << i0 << ", " << i1 << ", " << i2 << "] = " << grid[ind] << "\n" << endi);
	    }
	}
    }
    
    // Finally, clean up our garbage
    delete [] grid_nopad;
}

int GridforceGrid::get_box(Box *box, Vector pos) const
{
    Vector p = pos - origin;
    int inds[3];
    int ind;
    Vector g, dg;
    
    DebugM(3, "pos = " << pos << "\n" << "origin = " << origin << "\n" << "p = " << p << "\n" << endi);
    
    g = inv * p;
    for (int i = 0; i < 3; i++) {
	inds[i] = (int)floor(g[i]);
	dg[i] = g[i] - inds[i];
    }
    DebugM(3, "dg = " << dg << "\n" << endi);
    DebugM(3, "ind + dg = " << inds[0]+dg[0] << " " << inds[1]+dg[1] << " " << inds[2]+dg[2] << "\n" << endi);
    
    for (int i = 0; i < 3; i++) {
	if (inds[i] < 0 || inds[i] >= k[i]-1) {
	    if (cont[i]) inds[i] = k[i]-1;
	    else return -1;	// Outside potential and grid is not continuous
	}
    }
    
    ind = inds[0]*dk[0] + inds[1]*dk[1] + inds[2]*dk[2];
    
    for (int i0 = 0; i0 < 8; i0++) {
	int inds2[3];
	int ind2 = 0;
	
	float voff = 0.0;
	int bit = 1;	// bit = 2^i1 in the below loop
	for (int i1 = 0; i1 < 3; i1++) {
	    inds2[i1] = (inds[i1] + ((i0 & bit) ? 1 : 0)) % k[i1];
	    ind2 += inds2[i1] * dk[i1];
	    
	    // Deal with voltage offsets
	    if (cont[i1] && inds[i1] == (k[i1]-1) && inds2[i1] == 0) {
		voff += offset[i1];
		DebugM(3, "offset[" << i1 << "] = " << offset[i1] << "\n" << endi);
	    }
	    
	    bit <<= 1;	// i.e. multiply by 2
	}
	
	DebugM(3, "inds2 = " << inds2[0] << " " << inds2[1] << " " << inds2[2] << "\n" << endi);
	
	// NOTE: leaving everything in terms of unit cell coordinates for now,
	// eventually will multiply by inv tensor when applying the force
	
	// First set variables 'dk_{hi,lo}' (glob notation). The 'hi'
	// ('lo') variable in a given dimension is the number added (subtracted)
	// to go up (down) one grid point in that dimension; both are normally
	// just the corresponding 'dk[i]'. However, if we are sitting on a
	// boundary and we are using a continuous grid, then we want to map the
	// next point off the grid back around to the other side. e.g. point
	// (k[0], i1, k) maps to point (0, i1, k), which would be
	// accomplished by changing 'dk1_hi' to -(k[0]-1)*dk1.
	
	Bool zero_derivs = FALSE;
	int dk_hi[3];
	int dk_lo[3];
	float voffs[3];
	for (int i1 = 0; i1 < 3; i1++) {
	    if (inds2[i1] == 0) {
		if (cont[i1]) {
		    dk_hi[i1] = dk[i1];
		    dk_lo[i1] = -(k[i1]-1) * dk[i1];
		    voffs[i1] = offset[i1];
		}
		else zero_derivs = TRUE;
	    }
	    else if (inds2[i1] == k[i1]-1) {
		if (cont[i1]) {
		    dk_hi[i1] = -(k[i1]-1) * dk[i1];
		    dk_lo[i1] = dk[i1];
		    voffs[i1] = offset[i1];
		}
		else zero_derivs = TRUE;
	    }
	    else {
		dk_hi[i1] = dk[i1];
		dk_lo[i1] = dk[i1];
		voffs[i1] = 0.0;
	    }
	}
	
	DebugM(2, "zero_derivs = " << zero_derivs << "\n" << endi);
	DebugM(2, "dk_hi = " << dk_hi[0] << " " << dk_hi[1] << " " << dk_hi[2] << "\n" << endi);
	DebugM(2, "dk_lo = " << dk_lo[0] << " " << dk_lo[1] << " " << dk_lo[2] << "\n" << endi);
	
	
	// ### Now fill in box ###
	
	box->loc = dg;
	
	// V
	box->b[i0] = grid[ind2] + voff;
	
	if (zero_derivs) {
	    box->b[8+i0] = 0.0;
	    box->b[16+i0] = 0.0;
	    box->b[24+i0] = 0.0;
	    box->b[32+i0] = 0.0;
	    box->b[40+i0] = 0.0;
	    box->b[48+i0] = 0.0;
	    box->b[56+i0] = 0.0;
	} else {
	    box->b[8+i0] = 0.5 * (grid[ind2 + dk_hi[0]] - grid[ind2 - dk_lo[0]] + voffs[0]);	//  dV/dx
	    box->b[16+i0] = 0.5 * (grid[ind2 + dk_hi[1]] - grid[ind2 - dk_lo[1]] + voffs[1]);	//  dV/dy
	    box->b[24+i0] = 0.5 * (grid[ind2 + dk_hi[2]] - grid[ind2 - dk_lo[2]] + voffs[2]);	//  dV/dz
	    box->b[32+i0] = 0.25 * (grid[ind2 + dk_hi[0] + dk_hi[1]] - grid[ind2 - dk_lo[0] + dk_hi[1]]
		    - grid[ind2 + dk_hi[0] - dk_lo[1]] + grid[ind2 - dk_lo[0] - dk_lo[1]]);	//  d2V/dxdy
	    box->b[40+i0] = 0.25 * (grid[ind2 + dk_hi[0] + dk_hi[2]] - grid[ind2 - dk_lo[0] + dk_hi[2]]
		    - grid[ind2 + dk_hi[0] - dk_lo[2]] + grid[ind2 - dk_lo[0] - dk_lo[2]]);	//  d2V/dxdz
	    box->b[48+i0] = 0.25 * (grid[ind2 + dk_hi[1] + dk_hi[2]] - grid[ind2 - dk_lo[1] + dk_hi[2]]
		    - grid[ind2 + dk_hi[1] - dk_lo[2]] + grid[ind2 - dk_lo[1] - dk_lo[2]]);	//  d2V/dydz
	
	    box->b[56+i0] =									// d3V/dxdydz
		0.125 * (grid[ind2 + dk_hi[0] + dk_hi[1] + dk_hi[2]] - grid[ind2 + dk_hi[0]+ dk_hi[1] - dk_lo[2]]
			 - grid[ind2 + dk_hi[0] - dk_lo[1] + dk_hi[2]] - grid[ind2 - dk_lo[0] + dk_hi[1] + dk_hi[2]]
			 + grid[ind2 + dk_hi[0] - dk_lo[1] - dk_lo[2]] + grid[ind2 - dk_lo[0] + dk_hi[1] - dk_lo[2]]
			 + grid[ind2 - dk_lo[0] - dk_lo[1] + dk_hi[2]] - grid[ind2 - dk_lo[0] - dk_lo[1] - dk_lo[2]]);
	}
	
	DebugM(2, "V = " << box->b[i0] << "\n");
	
	DebugM(2, "dV/dx = " << box->b[8+i0] << "\n");
	DebugM(2, "dV/dy = " << box->b[16+i0] << "\n");
	DebugM(2, "dV/dz = " << box->b[24+i0] << "\n");
	
	DebugM(2, "d2V/dxdy = " << box->b[32+i0] << "\n");
	DebugM(2, "d2V/dxdz = " << box->b[40+i0] << "\n");
	DebugM(2, "d2V/dydz = " << box->b[48+i0] << "\n");
	
	DebugM(2, "d3V/dxdydz = " << box->b[56+i0] << "\n" << endi);
    }
    
    return 0;
}
