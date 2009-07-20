/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <iostream>

#include "GridForceGrid.h"
#include "Vector.h"
#include "Node.h"
#include "SimParameters.h"
#include "Molecule.h"
#include "InfoStream.h"
#include "common.h"
#include "ComputeGridForce.h"

#include "MGridforceParams.h"

#define MIN_DEBUG_LEVEL 3

//#define DEBUGM
#include "Debug.h"


GridforceGrid::GridforceGrid(int gridnum)
{
//    iout << iINFO << "GridforceGrid "
//         << gridnum << " created on " << CkMyPe() << "\n" << endi;
//    //    grid = NULL;
    mygridnum = gridnum;
    cont[0]=cont[1]=cont[2]=FALSE;
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
  //    delete [] grid;
}

void GridforceGrid::pack(MOStream *msg) const
{
    DebugM(2, "Packing message\n" << endi);
    
    msg->put(3*sizeof(int), (char*)k);
    msg->put(3*sizeof(int), (char*)k_nopad);
    msg->put(size);
    msg->put(size_nopad);
    msg->put(3*sizeof(int), (char*)dk);
    msg->put(3*sizeof(int), (char*)dk_nopad);
    msg->put(factor);
    
    msg->put(sizeof(Vector), (char*)&origin);
    msg->put(sizeof(Vector), (char*)&center);
    msg->put(sizeof(Tensor), (char*)&e);
    msg->put(sizeof(Tensor), (char*)&inv);
	     
//    msg->put(3*sizeof(float), (char*)pad_p);
//    msg->put(3*sizeof(float), (char*)pad_n);
    msg->put(3*sizeof(Bool), (char*)cont);
    msg->put(3*sizeof(float), (char*)offset);
    msg->put(3*sizeof(float), (char*)gap);
    msg->put(3*sizeof(float), (char*)gapinv);
    msg->put(sizeof(Vector), (char*)&scale);
    
    //    msg->put(size*sizeof(float), (char*)grid);
}

void GridforceGrid::unpack(MIStream *msg)
{
    DebugM(2, "Unpacking message\n" << endi);
//    iout << iINFO << CkMyPe() << " Unpacking message\n" << endi;
    
    msg->get(3*sizeof(int), (char*)k);
    msg->get(3*sizeof(int), (char*)k_nopad);
    msg->get(size);
    msg->get(size_nopad);
    msg->get(3*sizeof(int), (char*)dk);
    msg->get(3*sizeof(int), (char*)dk_nopad);
    msg->get(factor);
    
    msg->get(sizeof(Vector), (char*)&origin);
    msg->get(sizeof(Vector), (char*)&center);
    msg->get(sizeof(Tensor), (char*)&e);
    msg->get(sizeof(Tensor), (char*)&inv);
	     
//    msg->get(3*sizeof(float), (char*)pad_p);
//    msg->get(3*sizeof(float), (char*)pad_n);
    msg->get(3*sizeof(Bool), (char*)cont);
    msg->get(3*sizeof(float), (char*)offset);
    msg->get(3*sizeof(float), (char*)gap);
    msg->get(3*sizeof(float), (char*)gapinv);
    msg->get(sizeof(Vector), (char*)&scale);

    //    if (size) {
    //      delete [] grid;
    //      grid = new float[size];
    //      msg->get(size*sizeof(float), (char*)grid);
    //    }
    
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
//    DebugM(4, "cont = " << cont[0] << " " << cont[1] << " " << cont[2] << "\n");
//     DebugM(2, "offset = " << offset[0] << " " << offset[1] << " " << offset[2] << "\n");
//     DebugM(2, "gap = " << gap[0] << " " << gap[1] << " " << gap[2] << endi << "\n");
}

void GridforceGrid::init1(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams)
// On node 0: Read grid params, send to everyone
{
//  iout << iINFO << "[" << CkMyPe() << "] init1 called\n" << endi;
  poten_fp = Fopen(potfilename, "r");
  if (!poten_fp) {
    NAMD_die("Problem reading grid force potential file");
  }
    
  char line[256];
  do {
    fgets(line, 256, poten_fp);	// Read comment lines
  } while (line[0] == '#');
    
  sscanf(line, "object 1 class gridpositions counts %d %d %d\n",
         &k_nopad[0], &k_nopad[1], &k_nopad[2]);
    
  // Read origin
  fscanf(poten_fp, "origin %lf %lf %lf\n",
         &origin.x, &origin.y, &origin.z);
        
  // Read delta (unit vectors)
  // These are column vectors, so must fill gridfrcE tensor to reflect this
  fscanf(poten_fp, "delta %lf %lf %lf\n", &e.xx, &e.yx, &e.zx);
  fscanf(poten_fp, "delta %lf %lf %lf\n", &e.xy, &e.yy, &e.zy);
  fscanf(poten_fp, "delta %lf %lf %lf\n", &e.xz, &e.yz, &e.zz);
    
  center = origin + e * 0.5 
           * Position(k_nopad[0]-1, k_nopad[1]-1, k_nopad[2]-1);

  BigReal tmp[3];   // Temporary storage
  fscanf(poten_fp, "object 2 class gridconnections counts %lf %lf %lf\n",
         tmp, tmp+1, tmp+2);
  fscanf(poten_fp,
    "object 3 class array type double rank 0 items %lf data follows\n", 
    tmp);
    
  // Calculate inverse tensor
  BigReal det;
  det = e.xx*(e.yy*e.zz - e.yz*e.zy) - e.xy*(e.yx*e.zz - e.yz*e.zx) 
        + e.xz*(e.yx*e.zy - e.yy*e.zx);
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
  size_nopad = k_nopad[0] * k_nopad[1] * k_nopad[2];
  factor=1.0;
  if (mgridParams->gridforceVolts)
  {
    factor /= 0.0434;  // convert V -> kcal/mol*e
  }

  scale = mgridParams->gridforceScale;
  dk_nopad[3];
    
  // Shortcuts for accessing 1-D array with four indices
  dk_nopad[0] = k_nopad[1] * k_nopad[2];
  dk_nopad[1] = k_nopad[2];
  dk_nopad[2] = 1;
    
  Vector Kvec[3];
  Kvec[0] = e * Position(k_nopad[0]-1, 0, 0);
  Kvec[1] = e * Position(0, k_nopad[1]-1, 0);
  Kvec[2] = e * Position(0, 0, k_nopad[2]-1);
  Vector Avec[3];
  Avec[0] = simParams->lattice.a();
  Avec[1] = simParams->lattice.b();
  Avec[2] = simParams->lattice.c();
   
  // Decide whether we're wrapping
  Bool gridforceCont[3];
  gridforceCont[0] = simParams->gridforceContA1;
  gridforceCont[1] = simParams->gridforceContA2;
  gridforceCont[2] = simParams->gridforceContA3;

  for (int i0 = 0; i0 < 3; i0++) {
    if (mgridParams->gridforceCont[i0])
    {
      Bool found = FALSE;
      for (int i1 = 0; i1 < 3; i1++) {
	if (cross(Avec[i0].unit(), Kvec[i1].unit()).length() < 1e-4) {
	  found = TRUE;
	  cont[i1] = TRUE;
	  offset[i1] = mgridParams->gridforceVOffset[i0] * factor;
	  // want in grid-point units (normal = 1)
	  gap[i1] = (inv * (Avec[i0] - Kvec[i1])).length();	
          gapinv[i1] = 1.0/gap[i1];
		    
	  if (gap[i1] < 0) {
	    NAMD_die("Gridforce Grid overlap!");
          }
		    
          DebugM(4, "cont[" << i1 << "] = " << cont[i1] << "\n" << endi);
          DebugM(4, "gap[" << i1 << "] = " << gap[i1] << "\n" << endi);
          DebugM(4, "gapinv[" << i1 << "] = " << gapinv[i1] << "\n" << endi);
        }
      }
	    
      if (!found) {
        NAMD_die("No Gridforce unit vector found parallel to requested continuous grid direction!");
      }
    } else {
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
  DebugM(4, "delta = " << e * delta << " (" << delta << ")\n");
  origin += e * delta;
    
  size = k[0] * k[1] * k[2];
  dk[0] = k[1] * k[2];
  dk[1] = k[2];
  dk[2] = 1;
//  iout << iINFO 
//       << "size=" << size
//       << " dk=" << dk[0] << "," << dk[1] << "," << dk[2]
//       << " size_nopad=" << size_nopad
//       << " dk_nopad=" << dk_nopad[0] << "," 
//           << dk_nopad[1] << "," << dk_nopad[2]
//       << "\n" << endi;
  grid_elems_read = 0;
}

void GridforceGrid::init2()
// On all: Receive grid params... Figure out what grid values are needed
{
//  iout << iINFO << "[" << CkMyPe() << "] init2 called "
//       << " dk=" << dk[0] << "," << dk[1] << "," << dk[2]
//       << " k=" << k[0] << "," << k[1] << "," << k[2]
//       << " gridnum=" << mygridnum
//       << "\n" << endi;

  // Initialize the grid block parameters
  bxsz = (k[0] + GRIDBLOCKX - 1) / GRIDBLOCKX;
  bysz = (k[1] + GRIDBLOCKY - 1) / GRIDBLOCKY;
  bzsz = (k[2] + GRIDBLOCKZ - 1) / GRIDBLOCKZ;
            
  // We have the grid parameters now, so we can allocate the cache
  allocateGridCache();
  allocateGridRequestTables();
  precomputeIndexing();
  
  // Fill in most of new grid (all except pad areas)
  
  n_sum[0] = n_sum[1] = n_sum[2] = 0;
  p_sum[0] = p_sum[1] = p_sum[2] = 0;
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
	//        if (CkMyPe() == 1) {
	//          iout << iINFO << "[" << CkMyPe() << "]"
	//	       << "Need " << j0 << "," << j1 << "," << j2 << ":" << ind
	//	       << "\n" << endi;
	//	}
        if (isGridValCached(ind)) {
	  gridvals_needed[getGridBlockIdx(ind)] = true;
        }
      }
    }
  }
  // Initialize border areas
  if (!cont[0]) {
    for (int i0 = 0; i0 < border; i0++)
      for (int i1 = 0; i1 < k[1]; i1++)
	for (int i2 = 0; i2 < k[2]; i2++) {
	  int ind;
          ind = i0*dk[0] + i1*dk[1] + i2*dk[2];
          if (isGridValCached(ind)) {
  	    addGridValue(ind,0);
	  }
	  ind = (k[0]-1-i0)*dk[0] + i1*dk[1] + i2*dk[2];
          if (isGridValCached(ind)) {
  	    addGridValue(ind,0);
	  }
	}
  }
  if (!cont[1]) {
    for (int i0 = 0; i0 < k[0]; i0++)
      for (int i1 = 0; i1 < border; i1++)
	for (int i2 = 0; i2 < k[2]; i2++) {
	  int ind;
          ind = i0*dk[0] + i1*dk[1] + i2*dk[2];
          if (isGridValCached(ind)) {
  	    addGridValue(ind,0);
	  }
	  ind = i0*dk[0] + (k[1]-1-i1)*dk[1] + i2*dk[2];
          if (isGridValCached(ind)) {
  	    addGridValue(ind,0);
	  }
	}
  }
  if (!cont[2]) {
    for (int i0 = 0; i0 < k[0]; i0++)
      for (int i1 = 0; i1 < k[1]; i1++)
	for (int i2 = 0; i2 < border; i2++) {
	  int ind;
          ind = i0*dk[0] + i1*dk[1] + i2*dk[2];
          if (isGridValCached(ind)) {
  	    addGridValue(ind,0);
	  }
	  ind = i0*dk[0] + i1*dk[1] + (k[2]-1-i2)*dk[2];
          if (isGridValCached(ind)) {
  	    addGridValue(ind,0);
	  }
	}
  }
}

void GridforceGrid::init3(float *gridseg, int *start_indx, int *count)
// On node 0: When everyone has reported in, start reading and sending
// grid values
{
//  iout << iINFO << "[" << CkMyPe() << "] init3 called\n" << endi;
  
  int to_read = GRIDSEGSZ;
  int elems_left = size_nopad - grid_elems_read;
  if (elems_left < to_read) {
    to_read = elems_left;
  }
  *start_indx = grid_elems_read;
  *count = to_read;
//  iout << iINFO << "Sending " << *count << " elements starting at "
//       << *start_indx
//       << " size = " << size_nopad
//       << " left = " << elems_left
//       << "\n" << endi;
  float tmp2;
  int i;
  for (i = 0; i < to_read; i++) {
    int err = fscanf(poten_fp, "%f", &tmp2);
    if (err == EOF || err == 0) {
      NAMD_die("Grid force potential file incorrectly formatted");
    }
    gridseg[i] = tmp2 * factor;
    grid_elems_read++;
  }
  if (grid_elems_read == size_nopad && *count > 0) {
    fclose(poten_fp);
  }

  return;
}

void GridforceGrid::init4(float *gridseg, int start_indx, int count)
// On all: Receive grid vals, store them
{
//  iout << iINFO << "[" << CkMyPe() << "] init4 processing " 
//       << count << " elements starting at " << start_indx
//       << " mygridnum=" << mygridnum
//       << " size=" << size
//       << " dk=" << dk[0] << "," << dk[1] << "," << dk[2]
//       << " size_nopad=" << size_nopad
//       << " dk_nopad=" << dk_nopad[0] << "," 
//           << dk_nopad[1] << "," << dk_nopad[2]
//       << "\n" << endi;
  
  // PARAMETERS which probably ought to be somewhere else
  const BigReal modThresh = 1.0;
  
  int i;
  for(i=0; i < count; i++) {
    int ind_nopad = start_indx + i;
    int i0 = ind_nopad / dk_nopad[0];
    int i1 = (ind_nopad - i0*dk_nopad[0]) / dk_nopad[1];
    int i2 = (ind_nopad - i0*dk_nopad[0] - i1*dk_nopad[1]) / dk_nopad[2];
    int j0 = (cont[0]) ? i0 : i0 + border;
    int j1 = (cont[1]) ? i1 : i1 + border;
    int j2 = (cont[2]) ? i2 : i2 + border;
    int ind = j0 * dk[0] + j1*dk[1] + j2*dk[2];

//    int tmp = i0 * dk_nopad[0] + i1 * dk_nopad[1] + i2 * dk_nopad[2];
//    iout << iINFO << "[" << CkMyPe() << "] init4 processing " 
//       << " ind_nopad=" << ind_nopad
//       << " i0=" << i0 << "," << i1 << "," << i2
//       << " j0=" << j0 << "," << j1 << "," << j2
//       << " ind=" << ind << "," << tmp
//       << "\n" << endi;
    if (i0 == 0)
      n_sum[0] += gridseg[i];
    else if (i0 == k_nopad[0]-1)
      p_sum[0] += gridseg[i];
    if (i1 == 0)
      n_sum[1] += gridseg[i];
    else if (i1 == k_nopad[1]-1)
      p_sum[1] += gridseg[i];
    if (i2 == 0)
      n_sum[2] += gridseg[i];
    else if (i2 == k_nopad[2]-1)
      p_sum[2] += gridseg[i];
      
    if (gridvals_needed[getGridBlockIdx(ind)]) {
      addGridValue(ind,gridseg[i]);
    }
  }
  //  iout << iINFO 
  //     << "[" << CkMyPe() << "] Done w/ msg  "
  //     << grid_cache->size() << "," << gridvals_needed->size()
  //     << "\n" << endi;

  // See if we're all done so we can finish the border regions
  if (start_indx + count == size_nopad) {
//    iout << iINFO 
//         << "[" << CkMyPe() << "] Done reading grid " << mygridnum
//         << "\n" << endi;

    BigReal n_avg[3], p_avg[3];
    int i0;
    for (int i0 = 0; i0 < 3; i0++) {
      int i1 = (i0 + 1) % 3;
      int i2 = (i0 + 2) % 3;
      n_avg[i0] = n_sum[i0] / (k_nopad[i1] * k_nopad[i2]);
      p_avg[i0] = p_sum[i0] / (k_nopad[i1] * k_nopad[i2]);
	
      if (cont[i0] 
          && fabs(offset[i0] - (p_avg[i0]-n_avg[i0])) > modThresh) 
      {
        iout << iWARN << "GRID FORCE POTENTIAL DIFFERENCE IN K" << i0
	     << " DIRECTION IS " 
	     << offset[i0] - (p_avg[i0]-n_avg[i0]) 
	     << " KCAL/MOL*E\n" << endi;
      }
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
      DebugM(4, "pad_n[" << i << "] = " << pad_n[i] << "\n" << endi);
      DebugM(4, "pad_p[" << i << "] = " << pad_p[i] << "\n" << endi);
    }
    
    if (cont[0] && cont[1] && cont[2]) {
      // Nothing to do
      iout << iINFO << "FIX ME\n" << endi;
      // return;
    }
    
    // Now fill in rest of new grid
    for (int i0 = 0; i0 < k[0]; i0++) {
      for (int i1 = 0; i1 < k[1]; i1++) {
        for (int i2 = 0; i2 < k[2]; i2++) {
          if ( (cont[0] || (i0 >= border && i0 < k[0]-border)) 
               && (cont[1] || (i1 >= border && i1 < k[1]-border)) 
               && (cont[2] || i2 == border) )
          {
            i2 += k_nopad[2]-1;
            continue;
          }
		
          int ind = i0*dk[0] + i1*dk[1] + i2*dk[2];
          if (gridvals_needed[getGridBlockIdx(ind)]) {
            Position pos = e * Position(i0, i1, i2);
            int var[3] = {i0, i1, i2};
		
            for (int dir = 0; dir < 3; dir++) {
              if (cont[dir]) 
                continue;
  		    
              if (var[dir] < border)
                addGridValue(ind,pad_n[dir]);
              else if (var[dir] >= k[dir]-border)
                addGridValue(ind,pad_p[dir]);
            }
		
            DebugM(2, "grid[" << ind << "; " << i0 << ", " << i1
              << ", " << i2 << "] = " << get_grid(ind)
              << "\n" << endi);
          }
        }
      }
    }
//    iout << iINFO << "[" << CkMyPe() << "] cleaning up " 
//           << "\n" << endi;
    // Check cache
    //    GridCache::iterator it;
    //    int key0 = 0;
    //    for (it = grid_cache->begin(); it != grid_cache->end(); it++) {
    //      if (key0 != (*it).first) {
    //        iout << iINFO << "[" << CkMyPe() << "] key="
    //	     << key0 << "," << (*it).first
    //	     << "\n" << endi;
    //      }
    //      key0++;
    //    }
    delete [] gridvals_needed;
  }
}
  /*
void GridforceGrid::initialize(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams)
{

    iout << iINFO << "GridforceGrid initialize " << CkMyPe() << "\n" << endi;
    
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
    
    center = origin + e * 0.5 * Position(k_nopad[0]-1, k_nopad[1]-1, k_nopad[2]-1);

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
    if (mgridParams->gridforceVolts)
    {
	factor = 1.0/0.0434;  // convert V -> kcal/mol*e
    } else {
	factor = 1.0;
    }

    scale = mgridParams->gridforceScale;
    
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
    cont[0] = FALSE;
    cont[1] = FALSE;
    cont[2] = FALSE;
    Bool gridforceCont[3];
    gridforceCont[0] = simParams->gridforceContA1;
    gridforceCont[1] = simParams->gridforceContA2;
    gridforceCont[2] = simParams->gridforceContA3;

    for (int i0 = 0; i0 < 3; i0++) {
	if (mgridParams->gridforceCont[i0])
	{
	    Bool found = FALSE;
	    for (int i1 = 0; i1 < 3; i1++) {
		if (cross(Avec[i0].unit(), Kvec[i1].unit()).length() < 1e-4) {
		    found = TRUE;
		    cont[i1] = TRUE;
		    offset[i1] = mgridParams->gridforceVOffset[i0] * factor;
		    gap[i1] = (inv * (Avec[i0] - Kvec[i1])).length();	// want in grid-point units (normal = 1)
		    gapinv[i1] = 1.0/gap[i1];
		    
		    if (gap[i1] < 0) {
			NAMD_die("Gridforce Grid overlap!");
		    }
		    
		    DebugM(4, "cont[" << i1 << "] = " << cont[i1] << "\n" << endi);
		    DebugM(4, "gap[" << i1 << "] = " << gap[i1] << "\n" << endi);
		    DebugM(4, "gapinv[" << i1 << "] = " << gapinv[i1] << "\n" << endi);
		}
	    }
	    
	    if (!found) {
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
    DebugM(4, "delta = " << e * delta << " (" << delta << ")\n");
    origin += e * delta;
    
    size = k[0] * k[1] * k[2];
    dk[0] = k[1] * k[2];
    dk[1] = k[2];
    dk[2] = 1;
    
    grid = new float[size];
    elems_used = new char[size];
    
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
	DebugM(4, "n_avg[" << i0 << "] = " << n_avg[i0] << "\n" << endi);
	DebugM(4, "p_avg[" << i0 << "] = " << p_avg[i0] << "\n" << endi);
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
	DebugM(4, "pad_n[" << i << "] = " << pad_n[i] << "\n" << endi);
	DebugM(4, "pad_p[" << i << "] = " << pad_p[i] << "\n" << endi);
    }
    
    if (cont[0] && cont[1] && cont[2]) {
      // Nothing to do
      // return;
      iout << iINFO << "FIX ME\n" << endi;
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
//    iout << iINFO << "GridforceGrid initialize done " << CkMyPe() << "\n" << endi;
}
  */

void GridforceGrid::allocateGridCache()
{ 
  grid_cache = new float*[bxsz*bysz*bzsz];
  int i,j,k;
  for(i=0;i<bxsz;i++)
    for(j=0; j < bysz; j++)
      for(k=0; k < bzsz; k++) {
        grid_cache[(i*bysz+j)*bzsz + k] = 0;
      }
}

void GridforceGrid::allocateGridRequestTables()
{ 
  const int nodesz = CkNodeSize(CkMyNode());
  
  grid_proc_request = new Bool*[nodesz];
  int rank;
  for(rank=0;rank < nodesz; rank++)
    grid_proc_request[rank] = new Bool[bxsz*bysz*bzsz];

  gridvals_needed = new Bool[bxsz*bysz*bzsz];
  grid_request = new Bool[bxsz*bysz*bzsz];
    
  int i,j,k;
  for(i=0;i<bxsz;i++)
    for(j=0; j < bysz; j++)
      for(k=0; k < bzsz; k++) {
        const int idx = (i*bysz+j)*bzsz + k;
        for(rank=0;rank < nodesz; rank++)
          grid_proc_request[rank][idx] = false;
        gridvals_needed[idx] = false;
        grid_request[idx] = false;
      }
}

void GridforceGrid::allocateGridEntry(int bi)
{
  if (grid_cache[bi] != 0) {
    iout << iWARN << "GridforceGrid trying to allocate already allocated entry " 
         << bi
         << "\n" << endi;
  } else {
    const int blksz = GRIDBLOCKX*GRIDBLOCKY*GRIDBLOCKZ;
    float *block = new float[blksz];
    int i;
    for(i=0; i < blksz; i++)
      block[i] = 0.;
    
    grid_cache[bi] = block;
  }
  return;
}

void GridforceGrid::consolidateGridRequests() 
{
  const int maxblock = bxsz * bysz * bzsz;
  const int nodesz = CkNodeSize(CkMyNode());

  num_requests = 0;
  int block;
  for(block=0; block < maxblock; block++) {
    int rank;
    // Now copy results from proc request to make request
    for(rank=0; rank < nodesz; rank++) {
      const Bool procblock = grid_proc_request[rank][block];
      if (grid_proc_request[rank][block]) {
        if (!grid_request[block]) {
//          CkPrintf("[%d] Block %d requesed\n",rank,block);
          num_requests++;
          grid_request[block] = true;
        }
        grid_proc_request[rank][block] = false;
      }
    }
  }
}

void GridforceGrid::getGridIndices(int *gridStartIndex, 
                                   int *gridIndexList,
		                   int *gridNodeList,
				   int *gridNodeCount)
{
  Molecule *mol = Node::Object()->molecule;
  int gridnum;
  int grid_index=0;

//  CkPrintf("Starting getGridIndices\n");
  int i;
  for(i=0; i < CkNumNodes(); i++)
    gridNodeCount[i] = 0;

  // Count how much total data we need, and store the pointers
  // for the start of each grid's requested data
  for (gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
    const GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
    int req_sz = grid->getNumGridRequests();
    gridStartIndex[gridnum] = grid_index;
    grid_index += req_sz;
  }
  gridStartIndex[mol->numGridforceGrids] = grid_index;

  // Fill the gridIndexList
  grid_index = 0;
  for (gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
    const GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
    const int bxsz = (grid->k[0] + GRIDBLOCKX - 1) / GRIDBLOCKX;
    const int bysz = (grid->k[1] + GRIDBLOCKY - 1) / GRIDBLOCKY;
    const int bzsz = (grid->k[2] + GRIDBLOCKZ - 1) / GRIDBLOCKZ;
    const int maxblock = bxsz * bysz * bzsz;
    int blockx, blocky, blockz;
    int bix, biy, biz;
    for(blockx=0;blockx < bxsz; blockx++)
      for(blocky=0;blocky < bysz; blocky++)
        for(blockz=0;blockz < bzsz; blockz++)  {
          const int block = (blockx * bysz + blocky) * bzsz + blockz;
          Bool thisblock = grid->grid_request[block];
          if (thisblock) {
            gridIndexList[grid_index] = block;
            gridNodeList[grid_index] = grid->getGridValHomeNode(block);
            gridNodeCount[gridNodeList[grid_index]]++;
            grid_index++;
          }
          grid->grid_request[block] = false;
        }
  }
//  CkPrintf("End getGridIndices\n");
}

void GridforceGrid::precomputeIndexing()
{
  // Allocate index tables
  int dim;
  for(dim=0; dim < 3; dim++)  {
    block_index[dim] = new int[k[dim]];
    block_offset[dim] = new int[k[dim]];
    int inds;
    const int bsz[3] = { GRIDBLOCKX, GRIDBLOCKY, GRIDBLOCKZ };
    for(inds=0; inds < k[dim]; inds++) {
      const int b = inds / bsz[dim];
      const int o = inds - b * bsz[dim];
      block_index[dim][inds] = b;
      block_offset[dim][inds] = o;
    }
    
    int i0;
    const int bit = 1 << dim;
    for(i0=0; i0 < 8; i0++)  {
      grid_index_table[i0][dim] = new GridIndices[k[dim]];
      for(inds=0; inds < k[dim]; inds++) {
        GridIndices * const indx = &(grid_index_table[i0][dim][inds]);
        indx->dk_hi = -111;
        indx->dk_lo = -222;
        indx->zero_derivs = -333;
        indx->inds2 = (inds + ((i0 & bit) ? 1 : 0)) % k[dim];

        indx->zero_derivs = FALSE;
        if (indx->inds2 == 0) {
          if (cont[dim]) {
            indx->dk_hi = dk[dim];
            indx->dk_lo = -( k[dim] - 1 ) * dk[dim];
          } else indx->zero_derivs = TRUE;
        } else if (indx->inds2 == k[dim] - 1) {
          if (cont[dim]) {
            indx->dk_hi = -(k[dim]-1) * dk[dim];
            indx->dk_lo = dk[dim];
          } else indx->zero_derivs = TRUE;
        } else {
          indx->dk_hi = dk[dim];
          indx->dk_lo = dk[dim];
        }
      }
    }
  }
//  for(dim=0; dim < 3; dim++)
//    for(int i0=0; i0 < 8; i0++)
//      for(int inds=0; inds < k[dim]; inds++) {
//        GridIndices * const indx = &(grid_index_table[i0][dim][inds]);
//        CkPrintf("[%d,%d,%d] k[dim]=%d hi=%d lo=%d inds2=%d derivs=%d\n",
//                 dim,i0,inds,k[dim],
//                 indx->dk_hi,indx->dk_lo,indx->inds2, indx->zero_derivs);
//      }
    
}

int GridforceGrid::request_box(Vector pos)
{
    Vector p = pos - origin;
    int inds[3];
    int ind;
    Vector g;

    DebugM(3, "pos = " << pos << "\n" << "origin = " << origin << "\n" << "p = " << p << "\n" << endi);
    
    g = inv * p;
    for (int i = 0; i < 3; i++) {
	inds[i] = (int)floor(g[i]);
    }
   
    for (int i = 0; i < 3; i++) {
	if (inds[i] < 0 || inds[i] >= k[i]-1) {
	    if (cont[i])
	      inds[i] = k[i]-1;
	    else {
	      return -1;	// Outside potential and grid is not continuous
	    }
	}
    }
    ind = inds[0]*dk[0] + inds[1]*dk[1] + inds[2]*dk[2];
    
    for (int i0 = 0; i0 < 8; i0++) {
        int ind2 = 0;
        int zero_derivs = FALSE;
        int dk_hi[3], dk_lo[3];
        for(int i1 = 0; i1 < 3; i1++) {
          ind2 += grid_index_table[i0][i1][inds[i1]].inds2 * dk[i1];
          zero_derivs |= grid_index_table[i0][i1][inds[i1]].zero_derivs;
          dk_hi[i1] = grid_index_table[i0][i1][inds[i1]].dk_hi;
          dk_lo[i1] = grid_index_table[i0][i1][inds[i1]].dk_lo;
        }
        const int myrank = CkMyRank();
//        if (ind2 < 0 || ind2 >= bxsz * GRIDBLOCKX * bysz * GRIDBLOCKY * bzsz * GRIDBLOCKZ) {
//          CkPrintf("request_box err %d %d %d %d %d\n",
//                   ind2,bxsz,bysz,bzsz,bxsz*bysz*bzsz);
//          for(int i1=0;i1<3;i1++)  {
//          CkPrintf("i0=%d i1=%d inds[i1]=%d inds2=%d hi=%d lo=%d\n",
//          i0,i1,inds[i1],
//          grid_index_table[i0][i1][inds[i1]].inds2,
//          grid_index_table[i0][i1][inds[i1]].dk_hi,
//          grid_index_table[i0][i1][inds[i1]].dk_lo);
//          }
//        }
	req_grid(myrank,ind2);
	
	if (!zero_derivs) {
	    req_grid(myrank,ind2 + dk_hi[0]);
	    req_grid(myrank,ind2 - dk_lo[0]);	//  dV/dx
	    
            req_grid(myrank,ind2 + dk_hi[1]);
            req_grid(myrank,ind2 - dk_lo[1]);	//  dV/dy
            
	    req_grid(myrank,ind2 + dk_hi[2]);
	    req_grid(myrank,ind2 - dk_lo[2]);	//  dV/dz

	    req_grid(myrank,ind2 + dk_hi[0] + dk_hi[1]);
	    req_grid(myrank,ind2 - dk_lo[0] + dk_hi[1]);
	    req_grid(myrank,ind2 + dk_hi[0] - dk_lo[1]);
	    req_grid(myrank,ind2 - dk_lo[0] - dk_lo[1]);  //  d2V/dxdy

	    req_grid(myrank,ind2 + dk_hi[0] + dk_hi[2]);
	    req_grid(myrank,ind2 - dk_lo[0] + dk_hi[2]);
	    req_grid(myrank,ind2 + dk_hi[0] - dk_lo[2]);
	    req_grid(myrank,ind2 - dk_lo[0] - dk_lo[2]);	//  d2V/dxdz
	    
	    req_grid(myrank,ind2 + dk_hi[1] + dk_hi[2]);
	    req_grid(myrank,ind2 - dk_lo[1] + dk_hi[2]);
	    req_grid(myrank,ind2 + dk_hi[1] - dk_lo[2]);
	    req_grid(myrank,ind2 - dk_lo[1] - dk_lo[2]);	//  d2V/dydz
	
	    // d3V/dxdydz
            req_grid(myrank,ind2 + dk_hi[0] + dk_hi[1] + dk_hi[2]);
            req_grid(myrank,ind2 + dk_hi[0] + dk_hi[1] - dk_lo[2]);
            req_grid(myrank,ind2 + dk_hi[0] - dk_lo[1] + dk_hi[2]);
            req_grid(myrank,ind2 - dk_lo[0] + dk_hi[1] + dk_hi[2]);
            req_grid(myrank,ind2 + dk_hi[0] - dk_lo[1] - dk_lo[2]);
            req_grid(myrank,ind2 - dk_lo[0] + dk_hi[1] - dk_lo[2]);
            req_grid(myrank,ind2 - dk_lo[0] - dk_lo[1] + dk_hi[2]);
            req_grid(myrank,ind2 - dk_lo[0] - dk_lo[1] - dk_lo[2]);
	}
    }
    return 0;
}

int GridforceGrid::get_box(Box *box, Vector pos) const
{
    SimParameters *simParams = Node::Object()->simParameters;
    Vector p = pos - origin;
    int inds[3];
    int ind;
    Vector g, dg;
    Vector gapscale = Vector(1, 1, 1);
    
//    iout << iINFO << "Getting box " << pos << "\n" << endi;
    DebugM(3, "pos = " << pos << "\n" << "origin = " << origin << "\n" << "p = " << p << "\n" << endi);
    
    g = inv * p;
    for (int i = 0; i < 3; i++) {
	inds[i] = (int)floor(g[i]);
	dg[i] = g[i] - inds[i];
    }
    
    for (int i = 0; i < 3; i++) {
	if (inds[i] < 0 || inds[i] >= k[i]-1) {
	    if (cont[i]) inds[i] = k[i]-1;
	    else {
	      return -1;	// Outside potential and grid is not continuous
	    }
	}
	if (cont[i] && inds[i] == k[i]-1) {
	    // Correct for non-unit spacing between continuous grid images
	    gapscale[i] *= gapinv[i];
	    if (g[i] < 0.0) dg[i] = 1.0 + g[i]*gapinv[i]; // = (gap[i] + g[i]) * gapinv[i]
	    else dg[i] = (g[i] - inds[i]) * gapinv[i];
	}
    }
    
    DebugM(3, "gapscale = " << gapscale << "\n" << endi);
    box->scale = Tensor::diagonal(gapscale);
    
    DebugM(3, "dg = " << dg << "\n" << endi);
    DebugM(3, "ind + dg = " << inds[0]+dg[0] << " " << inds[1]+dg[1] << " " << inds[2]+dg[2] << "\n" << endi);
    
    ind = inds[0]*dk[0] + inds[1]*dk[1] + inds[2]*dk[2];
    
    for (int i0 = 0; i0 < 8; i0++) {
        int inds2[3];
        int ind2 = 0;
        int zero_derivs = FALSE;
        int dk_hi[3], dk_lo[3];
        for(int i1 = 0; i1 < 3; i1++) {
          inds2[i1] = grid_index_table[i0][i1][inds[i1]].inds2;
          ind2 += grid_index_table[i0][i1][inds[i1]].inds2 * dk[i1];
          zero_derivs |= grid_index_table[i0][i1][inds[i1]].zero_derivs;
          dk_hi[i1] = grid_index_table[i0][i1][inds[i1]].dk_hi;
          dk_lo[i1] = grid_index_table[i0][i1][inds[i1]].dk_lo;
        }
	
	float voff = 0.0;
	for (int i1 = 0; i1 < 3; i1++) {
	    // Deal with voltage offsets
	    if (cont[i1] && inds[i1] == (k[i1]-1) && inds2[i1] == 0) {
		voff += offset[i1];
		DebugM(3, "offset[" << i1 << "] = " << offset[i1] << "\n" << endi);
	    }
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
	
	float voffs[3];
	float dscales[3] = {0.5, 0.5, 0.5};
	for (int i1 = 0; i1 < 3; i1++) {
	    if (inds2[i1] == 0) {
		if (cont[i1]) {
		    voffs[i1] = offset[i1];
		    dscales[i1] = 1.0/(1.0 + gap[i1]) * 1.0/gapscale[i1];
		}
	    }
	    else if (inds2[i1] == k[i1]-1) {
		if (cont[i1]) {
		    voffs[i1] = offset[i1];
		    dscales[i1] = 1.0/(1.0 + gap[i1]) * 1.0/gapscale[i1];
		}
	    }
	    else {
		voffs[i1] = 0.0;
	    }
	}
	
	DebugM(2, "zero_derivs = " << zero_derivs << "\n" << endi);
	DebugM(2, "dk_hi = " << dk_hi[0] << " " << dk_hi[1] << " " << dk_hi[2] << "\n" << endi);
	DebugM(2, "dk_lo = " << dk_lo[0] << " " << dk_lo[1] << " " << dk_lo[2] << "\n" << endi);
	
	// ### Now fill in box ###
	
	box->loc = dg;
	
	// V
	box->b[i0] = get_grid(ind2) + voff;
	
	if (zero_derivs) {
	    box->b[8+i0] = 0.0;
	    box->b[16+i0] = 0.0;
	    box->b[24+i0] = 0.0;
	    box->b[32+i0] = 0.0;
	    box->b[40+i0] = 0.0;
	    box->b[48+i0] = 0.0;
	    box->b[56+i0] = 0.0;
	} else {
	    box->b[8+i0] = dscales[0] * (get_grid(ind2 + dk_hi[0]) - get_grid(ind2 - dk_lo[0]) + voffs[0]);	//  dV/dx
	    box->b[16+i0] = dscales[1] * (get_grid(ind2 + dk_hi[1]) - get_grid(ind2 - dk_lo[1]) + voffs[1]);	//  dV/dy
	    box->b[24+i0] = dscales[2] * (get_grid(ind2 + dk_hi[2]) - get_grid(ind2 - dk_lo[2]) + voffs[2]);	//  dV/dz
	    box->b[32+i0] = dscales[0] * dscales[1]
		* (get_grid(ind2 + dk_hi[0] + dk_hi[1]) - get_grid(ind2 - dk_lo[0] + dk_hi[1])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[1]) + get_grid(ind2 - dk_lo[0] - dk_lo[1]));	//  d2V/dxdy
	    box->b[40+i0] = dscales[0] * dscales[2]
		* (get_grid(ind2 + dk_hi[0] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] + dk_hi[2])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[2]) + get_grid(ind2 - dk_lo[0] - dk_lo[2]));	//  d2V/dxdz
	    box->b[48+i0] = dscales[1] * dscales[2]
		* (get_grid(ind2 + dk_hi[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[1] + dk_hi[2])
		   - get_grid(ind2 + dk_hi[1] - dk_lo[2]) + get_grid(ind2 - dk_lo[1] - dk_lo[2]));	//  d2V/dydz
	
	    box->b[56+i0] = dscales[0] * dscales[1] * dscales[2]			// d3V/dxdydz
		* (get_grid(ind2 + dk_hi[0] + dk_hi[1] + dk_hi[2]) - get_grid(ind2 + dk_hi[0]+ dk_hi[1] - dk_lo[2])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] + dk_hi[1] + dk_hi[2])
		   + get_grid(ind2 + dk_hi[0] - dk_lo[1] - dk_lo[2]) + get_grid(ind2 - dk_lo[0] + dk_hi[1] - dk_lo[2])
		   + get_grid(ind2 - dk_lo[0] - dk_lo[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] - dk_lo[1] - dk_lo[2]));
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
