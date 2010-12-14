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

#define MIN_DEBUG_LEVEL 2
//#define DEBUGM
#include "Debug.h"


GridforceGrid::GridforceGrid(void)
{
    cont[0] = cont[1] = cont[2] = FALSE;
}

GridforceGrid::~GridforceGrid() {
  //    delete [] grid;
}

void GridforceGrid::pack(MOStream *msg) const
{
    DebugM(2, "Packing message\n" << endi);
    
    msg->put(numSubgrids);
    msg->put(generation);
    
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
    
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->pack(msg);
    }
}

void GridforceGrid::unpack(MIStream *msg)
{
    DebugM(3, "Unpacking message\n" << endi);
//    iout << iINFO << CkMyPe() << " Unpacking message\n" << endi;
    
    msg->get(numSubgrids);
    msg->get(generation);
    
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

    if (numSubgrids) {
	DebugM(3, "Creating subgrids array, size " << numSubgrids << "\n");
	subgrids = new GridforceSubGrid *[numSubgrids];
	for (int i = 0; i < numSubgrids; i++) {
	    subgrids[i] = new GridforceSubGrid(this);
	    subgrids[i]->unpack(msg);
	}
    }
}

void GridforceGrid::readHeader(SimParameters *simParams, MGridforceParams *mgridParams)
{
  char line[256];
  long int poten_offset;
  do {
      poten_offset = ftell(poten_fp);
      fgets(line, 256, poten_fp);	// Read comment lines
      DebugM(4, "Read line: " << line << endi);
  } while (line[0] == '#');
  fseek(poten_fp, poten_offset, SEEK_SET);
  
  // read grid dimensions
  fscanf(poten_fp, "object %*d class gridpositions counts %d %d %d\n",
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

  fscanf(poten_fp, "object %*d class gridconnections counts %*lf %*lf %*lf\n");
  fscanf(poten_fp, "object %*d class array type double rank 0 items %*d data follows\n");
    
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
  
  DebugM(4, "origin = " << origin << "\n" << endi);
  DebugM(4, "e = " << e << "\n" << endi);
  DebugM(4, "inv = " << inv << "\n" << endi);
  
  // Allocate storage for potential
  size_nopad = k_nopad[0] * k_nopad[1] * k_nopad[2];
  factor=1.0;
  if (mgridParams->gridforceVolts)
  {
    factor /= 0.0434;  // convert V -> kcal/mol*e
  }

  scale = mgridParams->gridforceScale;
  grid_elems_read = 0;
}

void GridforceGrid::init2()
{
    // Initialize the grid block parameters
    bxsz = (k[0] + GRIDBLOCKX - 1) / GRIDBLOCKX;
    bysz = (k[1] + GRIDBLOCKY - 1) / GRIDBLOCKY;
    bzsz = (k[2] + GRIDBLOCKZ - 1) / GRIDBLOCKZ;
    
    DebugM(3, "blocksize = " << bxsz << " " << bysz << " " << bzsz << "\n" << endi);
    
    // We have the grid parameters now, so we can allocate the cache
    allocateGridCache();
    allocateGridRequestTables();
    precomputeIndexing();
    
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

void GridforceGrid::addDepositMsg(GridDepositMsg **outmsgs, int &arrayIdx, int startIdx, int grandTotalGrids)
{
    outmsgs[arrayIdx] = new GridDepositMsg;
    outmsgs[arrayIdx]->gridnum = startIdx + arrayIdx;
    outmsgs[arrayIdx]->grid = this;
    outmsgs[arrayIdx]->num_grids = grandTotalGrids;
    arrayIdx++;
    
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->addDepositMsg(outmsgs, arrayIdx, startIdx, grandTotalGrids);
    }
}


void GridforceGrid::readSubgridHierarchy(FILE *poten, int &totalGrids)
{
    DebugM(4, "Beginning of readSubgridHierarchy, generation = " << generation << ", totalGrids = " << totalGrids << "\n" << endi);
    
    int elems, generation_in;
    
    subgrids = new GridforceSubGrid *[numSubgrids];
    
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i] = new GridforceSubGrid(this);
	elems = fscanf(poten_fp, "# namdnugrid subgrid %d generation %d min %d %d %d max %d %d %d subgrids count %d\n",
	       &subgrids[i]->subgridIdx, &generation_in,
	       &subgrids[i]->pmin[0], &subgrids[i]->pmin[1], &subgrids[i]->pmin[2],
	       &subgrids[i]->pmax[0], &subgrids[i]->pmax[1], &subgrids[i]->pmax[2],
	       &subgrids[i]->numSubgrids);
	if (elems < 9) {
	    char msg[256];
	    sprintf(msg, "Problem reading Gridforce potential file! (%d < 9)", elems);
	    NAMD_die(msg);
	}
	
	totalGrids++;
	
	if (subgrids[i]->subgridIdx != (totalGrids - 1)) {
	    char msg[256];
	    sprintf(msg, "Problem reading Gridforce potential file! (%d != %d)", subgrids[i]->subgridIdx, totalGrids - 1);
	    NAMD_die(msg);
	}
	if (subgrids[i]->generation != generation_in) {
	    char msg[256];
	    sprintf(msg, "Problem reading Gridforce potential file! (%d != %d)", subgrids[i]->generation, generation_in);
	    NAMD_die(msg);
	}
	
// 	DebugM(3, "setting maingrid\n");
// 	subgrids[i]->maingrid->subgrids_flat[subgrids[i]->subgridIdx] = subgrids[i];
// 	DebugM(3, "reading subgrid hierarchy\n");
	
	subgrids[i]->readSubgridHierarchy(poten, totalGrids);
    }
}


void GridforceGrid::allocateGridCache()
{ 
  DebugM(3, "allocateGridCache() called for " << CkMyPe() << "\n" << endi);
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
  
  DebugM(4, "allocateGridRequestTables(): nodesz = " << nodesz << "\n" << endi);
  
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
      DebugM(3, "[" << CkMyPe() << "] allocateGridEntry(" << bi << ") called pointer " << (int64)this << " generation = " << generation << "\n" << endi);
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
  gridnum = 0;
  int maingridnum;
  int subgridIdx;
  for (maingridnum = 0; maingridnum < mol->numGridforceGrids; maingridnum++) {
      const GridforceMainGrid *maingrid = mol->get_gridfrc_grid(maingridnum);
      for (subgridIdx = 0; subgridIdx < maingrid->totalGrids; subgridIdx++) {
	  const GridforceGrid *grid = (subgridIdx == 0) ? (GridforceGrid *)maingrid : maingrid->subgrids_flat[subgridIdx-1];
	  int req_sz = grid->getNumGridRequests();
	  gridStartIndex[gridnum] = grid_index;
	  grid_index += req_sz;
	  gridnum++;
      }
  }
  gridStartIndex[gridnum] = grid_index;

  // Fill the gridIndexList
  grid_index = 0;
  gridnum = 0;
  for (maingridnum = 0; maingridnum < mol->numGridforceGrids; maingridnum++) {
      const GridforceMainGrid *maingrid = mol->get_gridfrc_grid(maingridnum);
      for (subgridIdx = 0; subgridIdx < maingrid->totalGrids; subgridIdx++) {
	  const GridforceGrid *grid = (subgridIdx == 0) ? (GridforceGrid *)maingrid : maingrid->subgrids_flat[subgridIdx-1];
	  const int bxsz = (grid->k[0] + GRIDBLOCKX - 1) / GRIDBLOCKX;
	  const int bysz = (grid->k[1] + GRIDBLOCKY - 1) / GRIDBLOCKY;
	  const int bzsz = (grid->k[2] + GRIDBLOCKZ - 1) / GRIDBLOCKZ;
	  const int maxblock = bxsz * bysz * bzsz;
	  int blockx, blocky, blockz;
	  int bix, biy, biz;
	  for(blockx=0;blockx < bxsz; blockx++) {
	      for(blocky=0;blocky < bysz; blocky++) {
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
	  }
	  gridnum++;
      }
  }
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
}

int GridforceGrid::request_box(Vector pos)
{
    int inds[3];
    int ind;
    Vector g, dg;
    Vector gapscale;
    
    int err = get_inds(pos, inds, dg, gapscale);
    if (err) {
	return -1;
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
    
    for (int i = 0; i < numSubgrids; i++) {
	if (inds[0] >= subgrids[i]->pmin[0] && inds[0] <= subgrids[i]->pmax[0] &&
	    inds[1] >= subgrids[i]->pmin[1] && inds[1] <= subgrids[i]->pmax[1] &&
	    inds[2] >= subgrids[i]->pmin[2] && inds[2] <= subgrids[i]->pmax[2])
	{
 	    DebugM(2, "calling request_box() for subgrid " << i << "\n" << endi);
	    subgrids[i]->request_box(pos);
	}
    }
    
    return 0;
}


int GridforceGrid::get_inds(Position pos, int *inds, Vector &dg, Vector &gapscale) const
{
    Vector p = pos - origin;
    Vector g;
    
    g = inv * p;
    
    for (int i = 0; i < 3; i++) {
	inds[i] = (int)floor(g[i]);
	dg[i] = g[i] - inds[i];
    }
    
    for (int i = 0; i < 3; i++) {
	if (inds[i] < 0 || inds[i] >= k[i]-1) {
	    if (cont[i]) inds[i] = k[i]-1;
	    else return -1;	// Outside potential and grid is not continuous
	}
	if (cont[i] && inds[i] == k[i]-1) {
	    // Correct for non-unit spacing between continuous grid images
	    gapscale[i] *= gapinv[i];
	    if (g[i] < 0.0) dg[i] = 1.0 + g[i]*gapinv[i]; // = (gap[i] + g[i]) * gapinv[i]
	    else dg[i] = (g[i] - inds[i]) * gapinv[i];
	}
    }
    
    return 0;
}


int GridforceGrid::compute_VdV(Position pos, float &V, Vector &dV) const
{
    SimParameters *simParams = Node::Object()->simParameters;
    int inds[3];
    int ind;
    Vector g, dg;
    Vector gapscale = Vector(1, 1, 1);
    
    int err = get_inds(pos, inds, dg, gapscale);
    if (err) {
	return -1;
    }
    
    DebugM(1, "gapscale = " << gapscale << "\n" << endi);
    DebugM(1, "dg = " << dg << "\n" << endi);
    DebugM(1, "ind + dg = " << inds[0]+dg[0] << " " << inds[1]+dg[1] << " " << inds[2]+dg[2] << "\n" << endi);
    DebugM(3, "compute_VdV: generation = " << generation << "\n" << endi);
    
    ind = inds[0]*dk[0] + inds[1]*dk[1] + inds[2]*dk[2];
    
    // Pass to subgrid if one exists here
    for (int i = 0; i < numSubgrids; i++) {
	if (((inds[0] >= subgrids[i]->pmin[0] && inds[0] <= subgrids[i]->pmax[0]) || subgrids[i]->cont[0]) &&
	    ((inds[1] >= subgrids[i]->pmin[1] && inds[1] <= subgrids[i]->pmax[1]) || subgrids[i]->cont[1]) &&
	    ((inds[2] >= subgrids[i]->pmin[2] && inds[2] <= subgrids[i]->pmax[2]) || subgrids[i]->cont[2]))
	{
	    return subgrids[i]->compute_VdV(pos, V, dV);
	}
    }
    
    // Compute b
    float b[64];	// Matrix of values at 8 box corners
    compute_b(b, inds, gapscale);
    for (int j = 0; j < 64; j++) DebugM(1, "b[" << j << "] = " << b[j] << "\n" << endi);
    
    // Compute a
    float a[64];
    compute_a(a, b);
    for (int j = 0; j < 64; j++) DebugM(1, "a[" << j << "] = " << a[j] << "\n" << endi);
	    
    // Calculate powers of x, y, z for later use
    // e.g. x[2] = x^2
    float x[4], y[4], z[4];
    x[0] = 1; y[0] = 1; z[0] = 1;
    for (int j = 1; j < 4; j++) {
	x[j] = x[j-1] * dg.x;
	y[j] = y[j-1] * dg.y;
	z[j] = z[j-1] * dg.z;
    }
    
    V = compute_V(a, x, y, z);
    dV = Tensor::diagonal(gapscale) * (compute_dV(a, x, y, z) * inv);
    
    return 0;
}


float GridforceGrid::compute_V(float *a, float *x, float *y, float *z) const
{
    float V = 0.0;
    int ind = 0;
    for (int l = 0; l < 4; l++) {
	for (int k = 0; k < 4; k++) {
	    for (int j = 0; j < 4; j++) {
		V += a[ind] * x[j] * y[k] * z[l];
		ind++;
	    }
	}
    }
    return V;
}


Vector GridforceGrid::compute_dV(float *a, float *x, float *y, float *z) const
{
    Vector dV = 0;
    int ind = 0;
    for (int l = 0; l < 4; l++) {
	for (int k = 0; k < 4; k++) {
	    for (int j = 0; j < 4; j++) {
		if (j > 0) dV.x += a[ind] * j * x[j-1] * y[k]   * z[l];		// dV/dx
		if (k > 0) dV.y += a[ind] * k * x[j]   * y[k-1] * z[l];		// dV/dy
		if (l > 0) dV.z += a[ind] * l * x[j]   * y[k]   * z[l-1];	// dV/dz
		ind++;
	    }
	}
    }
    return dV;
}


Vector GridforceGrid::compute_d2V(float *a, float *x, float *y, float *z) const
{
    Vector d2V = 0;
    int ind = 0;
    for (int l = 0; l < 4; l++) {
	for (int k = 0; k < 4; k++) {
	    for (int j = 0; j < 4; j++) {
		if (j > 0 && k > 0) d2V.x += a[ind] * j * k * x[j-1] * y[k-1] * z[l];	// d2V/dxdy
		if (j > 0 && l > 0) d2V.y += a[ind] * j * l * x[j-1] * y[k]   * z[l-1];	// d2V/dxdz
		if (k > 0 && l > 0) d2V.z += a[ind] * k * l * x[j]   * y[k-1] * z[l-1];	// d2V/dydz
		ind++;
	    }
	}
    }
    return d2V;
}


float GridforceGrid::compute_d3V(float *a, float *x, float *y, float *z) const
{
    float d3V = 0.0;
    int ind = 0;
    for (int l = 0; l < 4; l++) {
	for (int k = 0; k < 4; k++) {
	    for (int j = 0; j < 4; j++) {
		if (j > 0 && k > 0 && l > 0) d3V += a[ind] * j * k * l * x[j-1] * y[k-1] * z[l-1];	// d3V/dxdydz
		ind++;
	    }
	}
    }
    return d3V;
}


void GridforceGrid::compute_a(float *a, float *b) const
{
    // Static sparse 64x64 matrix times vector ... nicer looking way than this?
    a[0] = b[0];
    a[1] = b[8];
    a[2] = -3*b[0] + 3*b[1] - 2*b[8] - b[9];
    a[3] = 2*b[0] - 2*b[1] + b[8] + b[9];
    a[4] = b[16];
    a[5] = b[32];
    a[6] = -3*b[16] + 3*b[17] - 2*b[32] - b[33];
    a[7] = 2*b[16] - 2*b[17] + b[32] + b[33];
    a[8] = -3*b[0] + 3*b[2] - 2*b[16] - b[18];
    a[9] = -3*b[8] + 3*b[10] - 2*b[32] - b[34];
    a[10] = 9*b[0] - 9*b[1] - 9*b[2] + 9*b[3] + 6*b[8] + 3*b[9] - 6*b[10] - 3*b[11]
	+ 6*b[16] - 6*b[17] + 3*b[18] - 3*b[19] + 4*b[32] + 2*b[33] + 2*b[34] + b[35];
    a[11] = -6*b[0] + 6*b[1] + 6*b[2] - 6*b[3] - 3*b[8] - 3*b[9] + 3*b[10] + 3*b[11]
	- 4*b[16] + 4*b[17] - 2*b[18] + 2*b[19] - 2*b[32] - 2*b[33] - b[34] - b[35];
    a[12] = 2*b[0] - 2*b[2] + b[16] + b[18];
    a[13] = 2*b[8] - 2*b[10] + b[32] + b[34];
    a[14] = -6*b[0] + 6*b[1] + 6*b[2] - 6*b[3] - 4*b[8] - 2*b[9] + 4*b[10] + 2*b[11]
	- 3*b[16] + 3*b[17] - 3*b[18] + 3*b[19] - 2*b[32] - b[33] - 2*b[34] - b[35];
    a[15] = 4*b[0] - 4*b[1] - 4*b[2] + 4*b[3] + 2*b[8] + 2*b[9] - 2*b[10] - 2*b[11]
	+ 2*b[16] - 2*b[17] + 2*b[18] - 2*b[19] + b[32] + b[33] + b[34] + b[35];
    a[16] = b[24];
    a[17] = b[40];
    a[18] = -3*b[24] + 3*b[25] - 2*b[40] - b[41];
    a[19] = 2*b[24] - 2*b[25] + b[40] + b[41];
    a[20] = b[48];
    a[21] = b[56];
    a[22] = -3*b[48] + 3*b[49] - 2*b[56] - b[57];
    a[23] = 2*b[48] - 2*b[49] + b[56] + b[57];
    a[24] = -3*b[24] + 3*b[26] - 2*b[48] - b[50];
    a[25] = -3*b[40] + 3*b[42] - 2*b[56] - b[58];
    a[26] = 9*b[24] - 9*b[25] - 9*b[26] + 9*b[27] + 6*b[40] + 3*b[41] - 6*b[42] - 3*b[43]
	+ 6*b[48] - 6*b[49] + 3*b[50] - 3*b[51] + 4*b[56] + 2*b[57] + 2*b[58] + b[59];
    a[27] = -6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 3*b[40] - 3*b[41] + 3*b[42] + 3*b[43]
	- 4*b[48] + 4*b[49] - 2*b[50] + 2*b[51] - 2*b[56] - 2*b[57] - b[58] - b[59];
    a[28] = 2*b[24] - 2*b[26] + b[48] + b[50];
    a[29] = 2*b[40] - 2*b[42] + b[56] + b[58];
    a[30] = -6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 4*b[40] - 2*b[41] + 4*b[42] + 2*b[43]
	- 3*b[48] + 3*b[49] - 3*b[50] + 3*b[51] - 2*b[56] - b[57] - 2*b[58] - b[59];
    a[31] = 4*b[24] - 4*b[25] - 4*b[26] + 4*b[27] + 2*b[40] + 2*b[41] - 2*b[42] - 2*b[43]
	+ 2*b[48] - 2*b[49] + 2*b[50] - 2*b[51] + b[56] + b[57] + b[58] + b[59];
    a[32] = -3*b[0] + 3*b[4] - 2*b[24] - b[28];
    a[33] = -3*b[8] + 3*b[12] - 2*b[40] - b[44];
    a[34] = 9*b[0] - 9*b[1] - 9*b[4] + 9*b[5] + 6*b[8] + 3*b[9] - 6*b[12] - 3*b[13]
	+ 6*b[24] - 6*b[25] + 3*b[28] - 3*b[29] + 4*b[40] + 2*b[41] + 2*b[44] + b[45];
    a[35] = -6*b[0] + 6*b[1] + 6*b[4] - 6*b[5] - 3*b[8] - 3*b[9] + 3*b[12] + 3*b[13]
	- 4*b[24] + 4*b[25] - 2*b[28] + 2*b[29] - 2*b[40] - 2*b[41] - b[44] - b[45];
    a[36] = -3*b[16] + 3*b[20] - 2*b[48] - b[52];
    a[37] = -3*b[32] + 3*b[36] - 2*b[56] - b[60];
    a[38] = 9*b[16] - 9*b[17] - 9*b[20] + 9*b[21] + 6*b[32] + 3*b[33] - 6*b[36] - 3*b[37]
	+ 6*b[48] - 6*b[49] + 3*b[52] - 3*b[53] + 4*b[56] + 2*b[57] + 2*b[60] + b[61];
    a[39] = -6*b[16] + 6*b[17] + 6*b[20] - 6*b[21] - 3*b[32] - 3*b[33] + 3*b[36] + 3*b[37]
	- 4*b[48] + 4*b[49] - 2*b[52] + 2*b[53] - 2*b[56] - 2*b[57] - b[60] - b[61];
    a[40] = 9*b[0] - 9*b[2] - 9*b[4] + 9*b[6] + 6*b[16] + 3*b[18] - 6*b[20] - 3*b[22]
	+ 6*b[24] - 6*b[26] + 3*b[28] - 3*b[30] + 4*b[48] + 2*b[50] + 2*b[52] + b[54];
    a[41] = 9*b[8] - 9*b[10] - 9*b[12] + 9*b[14] + 6*b[32] + 3*b[34] - 6*b[36] - 3*b[38]
	+ 6*b[40] - 6*b[42] + 3*b[44] - 3*b[46] + 4*b[56] + 2*b[58] + 2*b[60] + b[62];
    a[42] = -27*b[0] + 27*b[1] + 27*b[2] - 27*b[3] + 27*b[4] - 27*b[5] - 27*b[6] + 27*b[7]
	- 18*b[8] - 9*b[9] + 18*b[10] + 9*b[11] + 18*b[12] + 9*b[13] - 18*b[14] - 9*b[15]
	- 18*b[16] + 18*b[17] - 9*b[18] + 9*b[19] + 18*b[20] - 18*b[21] + 9*b[22] - 9*b[23]
	- 18*b[24] + 18*b[25] + 18*b[26] - 18*b[27] - 9*b[28] + 9*b[29] + 9*b[30] - 9*b[31]
	- 12*b[32] - 6*b[33] - 6*b[34] - 3*b[35] + 12*b[36] + 6*b[37] + 6*b[38] + 3*b[39]
	- 12*b[40] - 6*b[41] + 12*b[42] + 6*b[43] - 6*b[44] - 3*b[45] + 6*b[46] + 3*b[47]
	- 12*b[48] + 12*b[49] - 6*b[50] + 6*b[51] - 6*b[52] + 6*b[53] - 3*b[54] + 3*b[55]
	- 8*b[56] - 4*b[57] - 4*b[58] - 2*b[59] - 4*b[60] - 2*b[61] - 2*b[62] - b[63];
    a[43] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
	+ 9*b[8] + 9*b[9] - 9*b[10] - 9*b[11] - 9*b[12] - 9*b[13] + 9*b[14] + 9*b[15]
	+ 12*b[16] - 12*b[17] + 6*b[18] - 6*b[19] - 12*b[20] + 12*b[21] - 6*b[22] + 6*b[23]
	+ 12*b[24] - 12*b[25] - 12*b[26] + 12*b[27] + 6*b[28] - 6*b[29] - 6*b[30] + 6*b[31]
	+ 6*b[32] + 6*b[33] + 3*b[34] + 3*b[35] - 6*b[36] - 6*b[37] - 3*b[38] - 3*b[39]
	+ 6*b[40] + 6*b[41] - 6*b[42] - 6*b[43] + 3*b[44] + 3*b[45] - 3*b[46] - 3*b[47]
	+ 8*b[48] - 8*b[49] + 4*b[50] - 4*b[51] + 4*b[52] - 4*b[53] + 2*b[54] - 2*b[55]
	+ 4*b[56] + 4*b[57] + 2*b[58] + 2*b[59] + 2*b[60] + 2*b[61] + b[62] + b[63];
    a[44] = -6*b[0] + 6*b[2] + 6*b[4] - 6*b[6] - 3*b[16] - 3*b[18] + 3*b[20] + 3*b[22]
	- 4*b[24] + 4*b[26] - 2*b[28] + 2*b[30] - 2*b[48] - 2*b[50] - b[52] - b[54];
    a[45] = -6*b[8] + 6*b[10] + 6*b[12] - 6*b[14] - 3*b[32] - 3*b[34] + 3*b[36] + 3*b[38]
	- 4*b[40] + 4*b[42] - 2*b[44] + 2*b[46] - 2*b[56] - 2*b[58] - b[60] - b[62];
    a[46] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
	+ 12*b[8] + 6*b[9] - 12*b[10] - 6*b[11] - 12*b[12] - 6*b[13] + 12*b[14] + 6*b[15]
	+ 9*b[16] - 9*b[17] + 9*b[18] - 9*b[19] - 9*b[20] + 9*b[21] - 9*b[22] + 9*b[23]
	+ 12*b[24] - 12*b[25] - 12*b[26] + 12*b[27] + 6*b[28] - 6*b[29] - 6*b[30] + 6*b[31]
	+ 6*b[32] + 3*b[33] + 6*b[34] + 3*b[35] - 6*b[36] - 3*b[37] - 6*b[38] - 3*b[39]
	+ 8*b[40] + 4*b[41] - 8*b[42] - 4*b[43] + 4*b[44] + 2*b[45] - 4*b[46] - 2*b[47]
	+ 6*b[48] - 6*b[49] + 6*b[50] - 6*b[51] + 3*b[52] - 3*b[53] + 3*b[54] - 3*b[55]
	+ 4*b[56] + 2*b[57] + 4*b[58] + 2*b[59] + 2*b[60] + b[61] + 2*b[62] + b[63];
    a[47] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
	- 6*b[8] - 6*b[9] + 6*b[10] + 6*b[11] + 6*b[12] + 6*b[13] - 6*b[14] - 6*b[15]
	- 6*b[16] + 6*b[17] - 6*b[18] + 6*b[19] + 6*b[20] - 6*b[21] + 6*b[22] - 6*b[23]
	- 8*b[24] + 8*b[25] + 8*b[26] - 8*b[27] - 4*b[28] + 4*b[29] + 4*b[30] - 4*b[31]
	- 3*b[32] - 3*b[33] - 3*b[34] - 3*b[35] + 3*b[36] + 3*b[37] + 3*b[38] + 3*b[39]
	- 4*b[40] - 4*b[41] + 4*b[42] + 4*b[43] - 2*b[44] - 2*b[45] + 2*b[46] + 2*b[47]
	- 4*b[48] + 4*b[49] - 4*b[50] + 4*b[51] - 2*b[52] + 2*b[53] - 2*b[54] + 2*b[55]
	- 2*b[56] - 2*b[57] - 2*b[58] - 2*b[59] - b[60] - b[61] - b[62] - b[63];
    a[48] = 2*b[0] - 2*b[4] + b[24] + b[28];
    a[49] = 2*b[8] - 2*b[12] + b[40] + b[44];
    a[50] = -6*b[0] + 6*b[1] + 6*b[4] - 6*b[5] - 4*b[8] - 2*b[9] + 4*b[12] + 2*b[13]
	- 3*b[24] + 3*b[25] - 3*b[28] + 3*b[29] - 2*b[40] - b[41] - 2*b[44] - b[45];
    a[51] = 4*b[0] - 4*b[1] - 4*b[4] + 4*b[5] + 2*b[8] + 2*b[9] - 2*b[12] - 2*b[13]
	+ 2*b[24] - 2*b[25] + 2*b[28] - 2*b[29] + b[40] + b[41] + b[44] + b[45];
    a[52] = 2*b[16] - 2*b[20] + b[48] + b[52];
    a[53] = 2*b[32] - 2*b[36] + b[56] + b[60];
    a[54] = -6*b[16] + 6*b[17] + 6*b[20] - 6*b[21] - 4*b[32] - 2*b[33] + 4*b[36] + 2*b[37]
	- 3*b[48] + 3*b[49] - 3*b[52] + 3*b[53] - 2*b[56] - b[57] - 2*b[60] - b[61];
    a[55] = 4*b[16] - 4*b[17] - 4*b[20] + 4*b[21] + 2*b[32] + 2*b[33] - 2*b[36] - 2*b[37]
	+ 2*b[48] - 2*b[49] + 2*b[52] - 2*b[53] + b[56] + b[57] + b[60] + b[61];
    a[56] = -6*b[0] + 6*b[2] + 6*b[4] - 6*b[6] - 4*b[16] - 2*b[18] + 4*b[20] + 2*b[22]
	- 3*b[24] + 3*b[26] - 3*b[28] + 3*b[30] - 2*b[48] - b[50] - 2*b[52] - b[54];
    a[57] = -6*b[8] + 6*b[10] + 6*b[12] - 6*b[14] - 4*b[32] - 2*b[34] + 4*b[36] + 2*b[38]
	- 3*b[40] + 3*b[42] - 3*b[44] + 3*b[46] - 2*b[56] - b[58] - 2*b[60] - b[62];
    a[58] = 18*b[0] - 18*b[1] - 18*b[2] + 18*b[3] - 18*b[4] + 18*b[5] + 18*b[6] - 18*b[7]
	+ 12*b[8] + 6*b[9] - 12*b[10] - 6*b[11] - 12*b[12] - 6*b[13] + 12*b[14] + 6*b[15]
	+ 12*b[16] - 12*b[17] + 6*b[18] - 6*b[19] - 12*b[20] + 12*b[21] - 6*b[22] + 6*b[23]
	+ 9*b[24] - 9*b[25] - 9*b[26] + 9*b[27] + 9*b[28] - 9*b[29] - 9*b[30] + 9*b[31]
	+ 8*b[32] + 4*b[33] + 4*b[34] + 2*b[35] - 8*b[36] - 4*b[37] - 4*b[38] - 2*b[39]
	+ 6*b[40] + 3*b[41] - 6*b[42] - 3*b[43] + 6*b[44] + 3*b[45] - 6*b[46] - 3*b[47]
	+ 6*b[48] - 6*b[49] + 3*b[50] - 3*b[51] + 6*b[52] - 6*b[53] + 3*b[54] - 3*b[55]
	+ 4*b[56] + 2*b[57] + 2*b[58] + b[59] + 4*b[60] + 2*b[61] + 2*b[62] + b[63];
    a[59] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
	- 6*b[8] - 6*b[9] + 6*b[10] + 6*b[11] + 6*b[12] + 6*b[13] - 6*b[14] - 6*b[15]
	- 8*b[16] + 8*b[17] - 4*b[18] + 4*b[19] + 8*b[20] - 8*b[21] + 4*b[22] - 4*b[23]
	- 6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 6*b[28] + 6*b[29] + 6*b[30] - 6*b[31]
	- 4*b[32] - 4*b[33] - 2*b[34] - 2*b[35] + 4*b[36] + 4*b[37] + 2*b[38] + 2*b[39]
	- 3*b[40] - 3*b[41] + 3*b[42] + 3*b[43] - 3*b[44] - 3*b[45] + 3*b[46] + 3*b[47]
	- 4*b[48] + 4*b[49] - 2*b[50] + 2*b[51] - 4*b[52] + 4*b[53] - 2*b[54] + 2*b[55]
	- 2*b[56] - 2*b[57] - b[58] - b[59] - 2*b[60] - 2*b[61] - b[62] - b[63];
    a[60] = 4*b[0] - 4*b[2] - 4*b[4] + 4*b[6] + 2*b[16] + 2*b[18] - 2*b[20] - 2*b[22]
	+ 2*b[24] - 2*b[26] + 2*b[28] - 2*b[30] + b[48] + b[50] + b[52] + b[54];
    a[61] = 4*b[8] - 4*b[10] - 4*b[12] + 4*b[14] + 2*b[32] + 2*b[34] - 2*b[36] - 2*b[38]
	+ 2*b[40] - 2*b[42] + 2*b[44] - 2*b[46] + b[56] + b[58] + b[60] + b[62];
    a[62] = -12*b[0] + 12*b[1] + 12*b[2] - 12*b[3] + 12*b[4] - 12*b[5] - 12*b[6] + 12*b[7]
	- 8*b[8] - 4*b[9] + 8*b[10] + 4*b[11] + 8*b[12] + 4*b[13] - 8*b[14] - 4*b[15]
	- 6*b[16] + 6*b[17] - 6*b[18] + 6*b[19] + 6*b[20] - 6*b[21] + 6*b[22] - 6*b[23]
	- 6*b[24] + 6*b[25] + 6*b[26] - 6*b[27] - 6*b[28] + 6*b[29] + 6*b[30] - 6*b[31]
	- 4*b[32] - 2*b[33] - 4*b[34] - 2*b[35] + 4*b[36] + 2*b[37] + 4*b[38] + 2*b[39]
	- 4*b[40] - 2*b[41] + 4*b[42] + 2*b[43] - 4*b[44] - 2*b[45] + 4*b[46] + 2*b[47]
	- 3*b[48] + 3*b[49] - 3*b[50] + 3*b[51] - 3*b[52] + 3*b[53] - 3*b[54] + 3*b[55]
	- 2*b[56] - b[57] - 2*b[58] - b[59] - 2*b[60] - b[61] - 2*b[62] - b[63];
    a[63] = 8*b[0] - 8*b[1] - 8*b[2] + 8*b[3] - 8*b[4] + 8*b[5] + 8*b[6] - 8*b[7]
	+ 4*b[8] + 4*b[9] - 4*b[10] - 4*b[11] - 4*b[12] - 4*b[13] + 4*b[14] + 4*b[15]
	+ 4*b[16] - 4*b[17] + 4*b[18] - 4*b[19] - 4*b[20] + 4*b[21] - 4*b[22] + 4*b[23]
	+ 4*b[24] - 4*b[25] - 4*b[26] + 4*b[27] + 4*b[28] - 4*b[29] - 4*b[30] + 4*b[31]
	+ 2*b[32] + 2*b[33] + 2*b[34] + 2*b[35] - 2*b[36] - 2*b[37] - 2*b[38] - 2*b[39]
	+ 2*b[40] + 2*b[41] - 2*b[42] - 2*b[43] + 2*b[44] + 2*b[45] - 2*b[46] - 2*b[47]
	+ 2*b[48] - 2*b[49] + 2*b[50] - 2*b[51] + 2*b[52] - 2*b[53] + 2*b[54] - 2*b[55]
	+ b[56] + b[57] + b[58] + b[59] + b[60] + b[61] + b[62] + b[63];
}


/*********************/
/* GRIDFORCEMAINGRID */
/*********************/

GridforceMainGrid::GridforceMainGrid(int gridnum)
{
    mygridnum = gridnum;
    generation = 0;
}


void GridforceMainGrid::pack(MOStream *msg) const
{
    DebugM(4, "Packing maingrid\n" << endi);
    
    GridforceGrid::pack(msg);
    
//     msg->put(3*sizeof(float), (char*)pad_p);
//     msg->put(3*sizeof(float), (char*)pad_n);
    msg->put(totalGrids);
    msg->put(mygridnum);
}


void GridforceMainGrid::unpack(MIStream *msg)
{
    DebugM(4, "Unpacking maingrid\n" << endi);
    
    GridforceGrid::unpack(msg);
    
//     msg->get(3*sizeof(float), (char*)pad_p);
//     msg->get(3*sizeof(float), (char*)pad_n);
    msg->get(totalGrids);
    msg->get(mygridnum);
    
    buildSubgridsFlat();
    
    DebugM(4, "size  = " << size << "\n" << endi);
    DebugM(4, "numSubgrids = " << numSubgrids << "\n" << endi);
    DebugM(4, "gapinv = " << gapinv[0] << " " << gapinv[2] << " " << gapinv[2] << " " << "\n" << endi);
    DebugM(4, "generation = " << generation << "\n" << endi);
}


void GridforceMainGrid::buildSubgridsFlat(void)
{
    DebugM(4, "buildSubgridsFlat() called, totalGrids-1 = " << totalGrids-1 << "\n" << endi);
    subgrids_flat = new GridforceSubGrid *[totalGrids-1];
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->addToSubgridsFlat();
    }
    for (int i = 0; i < totalGrids-1; i++) {
	DebugM(4, "subgrids_flat[" << i << "]->numSubgrids = " << subgrids_flat[i]->numSubgrids << "\n" << endi);
    }
    for (int i = 0; i < numSubgrids; i++) {
	DebugM(4, "subgrids[" << i << "]->numSubgrids = " << subgrids[i]->numSubgrids << "\n" << endi);
    }
}


void GridforceMainGrid::init1(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams)
{
    //FILE *poten = Fopen(potfilename, "r");
    poten_fp = Fopen(potfilename, "r");
    if (!poten_fp) {
	NAMD_die("Problem reading grid force potential file");
    }
    
    // Read special comment fields and create subgrid objects
    totalGrids = 1;
    char line[256];
    Bool flag = FALSE;
    numSubgrids = 0;
    float version;
    long int poten_offset;
    do {
	poten_offset = ftell(poten_fp);
	fgets(line, 256, poten_fp);	// Read comment lines
	//flag = sscanf(line, "# maingrid subgrids count %d\n", &numSubgrids);
	flag = sscanf(line, "# namdnugrid version %f\n", &version);
    } while (line[0] == '#' && !flag);
    
    
    if (flag) {
	if (version != 1.0) {
	    NAMD_die("Unsupported version of non-uniform grid file format!");
	}
	fscanf(poten_fp, "# namdnugrid maingrid subgrids count %d\n", &numSubgrids);
	readSubgridHierarchy(poten_fp, totalGrids);
	buildSubgridsFlat();
    } else {
	fseek(poten_fp, poten_offset, SEEK_SET);
    }
    
    readHeader(simParams, mgridParams);
    
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
    
    // Call for each subgrid
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->init1(potfilename, simParams, mgridParams);
    }
}


void GridforceMainGrid::init2()
// On all: Receive grid params... Figure out what grid values are needed
{
    GridforceGrid::init2();
    
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
		if (isGridValCached(ind)) {
		    gridvals_needed[getGridBlockIdx(ind)] = true;
		}
	    }
	}
    }
    // Initialize border areas
    if (!cont[0]) {
	for (int i0 = 0; i0 < border; i0++) {
	    for (int i1 = 0; i1 < k[1]; i1++) {
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
	}
    }
    if (!cont[1]) {
	for (int i0 = 0; i0 < k[0]; i0++) {
	    for (int i1 = 0; i1 < border; i1++) {
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
	}
    }
    if (!cont[2]) {
	for (int i0 = 0; i0 < k[0]; i0++) {
	    for (int i1 = 0; i1 < k[1]; i1++) {
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
    }
}


void GridforceMainGrid::init4(float *gridseg, int start_indx, int count)
// On all: Receive grid vals, store them
{
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

  // See if we're all done so we can finish the border regions
  if (start_indx + count == size_nopad) {
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
    delete [] gridvals_needed;
  }
}


int GridforceMainGrid::buildDepositMsgs(GridDepositMsg **outmsgs, int startIdx, int grandTotalGrids)
{
    int arrayIdx = 0;
    addDepositMsg(outmsgs, arrayIdx, startIdx, grandTotalGrids);
    
    for (int i = 0; i < totalGrids; i++) {
	DebugM(4, "msg->gridnum = " << outmsgs[i]->gridnum << "\n" << endi);
	DebugM(4, "msg->num_grids = " << outmsgs[i]->num_grids << "\n" << endi);
    }
    
    return totalGrids;
}


void GridforceMainGrid::compute_b(float *b, int *inds, Vector gapscale) const
{
    for (int i0 = 0; i0 < 8; i0++) {
	int inds2[3];
	int ind2 = 0;
	int zero_derivs = FALSE;
	int dk_hi[3], dk_lo[3];
	for (int i1 = 0; i1 < 3; i1++) {
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

 	DebugM(1, "inds2 = " << inds2[0] << " " << inds2[1] << " " << inds2[2] << "\n" << endi);
	
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
	
	DebugM(2, "cont = " << cont[0] << " " << cont[1] << " " << cont[2] << "\n" << endi);
	DebugM(2, "zero_derivs = " << zero_derivs << "\n" << endi);
	DebugM(2, "dk_hi = " << dk_hi[0] << " " << dk_hi[1] << " " << dk_hi[2] << "\n" << endi);
	DebugM(2, "dk_lo = " << dk_lo[0] << " " << dk_lo[1] << " " << dk_lo[2] << "\n" << endi);
	DebugM(2, "dscales = " << dscales[0] << " " << dscales[1] << " " << dscales[2] << "\n" << endi);
	
	// V
	b[i0] = get_grid(ind2) + voff;
	
	if (zero_derivs) {
	    DebugM(2, "zero_derivs\n" << endi);
	    b[8+i0] = 0.0;
	    b[16+i0] = 0.0;
	    b[24+i0] = 0.0;
	    b[32+i0] = 0.0;
	    b[40+i0] = 0.0;
	    b[48+i0] = 0.0;
	    b[56+i0] = 0.0;
	} else {
	    b[8+i0] = dscales[0] * (get_grid(ind2 + dk_hi[0]) - get_grid(ind2 - dk_lo[0]) + voffs[0]);	//  dV/dx
	    b[16+i0] = dscales[1] * (get_grid(ind2 + dk_hi[1]) - get_grid(ind2 - dk_lo[1]) + voffs[1]);	//  dV/dy
	    b[24+i0] = dscales[2] * (get_grid(ind2 + dk_hi[2]) - get_grid(ind2 - dk_lo[2]) + voffs[2]);	//  dV/dz
	    b[32+i0] = dscales[0] * dscales[1]
		* (get_grid(ind2 + dk_hi[0] + dk_hi[1]) - get_grid(ind2 - dk_lo[0] + dk_hi[1])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[1]) + get_grid(ind2 - dk_lo[0] - dk_lo[1]));	//  d2V/dxdy
	    b[40+i0] = dscales[0] * dscales[2]
		* (get_grid(ind2 + dk_hi[0] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] + dk_hi[2])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[2]) + get_grid(ind2 - dk_lo[0] - dk_lo[2]));	//  d2V/dxdz
	    b[48+i0] = dscales[1] * dscales[2]
		* (get_grid(ind2 + dk_hi[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[1] + dk_hi[2])
		   - get_grid(ind2 + dk_hi[1] - dk_lo[2]) + get_grid(ind2 - dk_lo[1] - dk_lo[2]));	//  d2V/dydz
	
	    b[56+i0] = dscales[0] * dscales[1] * dscales[2]					// d3V/dxdydz
		* (get_grid(ind2 + dk_hi[0] + dk_hi[1] + dk_hi[2]) - get_grid(ind2 + dk_hi[0]+ dk_hi[1] - dk_lo[2])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] + dk_hi[1] + dk_hi[2])
		   + get_grid(ind2 + dk_hi[0] - dk_lo[1] - dk_lo[2]) + get_grid(ind2 - dk_lo[0] + dk_hi[1] - dk_lo[2])
		   + get_grid(ind2 - dk_lo[0] - dk_lo[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] - dk_lo[1] - dk_lo[2]));
	}
	
	DebugM(2, "V = " << b[i0] << "\n");
	
	DebugM(2, "dV/dx = " << b[8+i0] << "\n");
	DebugM(2, "dV/dy = " << b[16+i0] << "\n");
	DebugM(2, "dV/dz = " << b[24+i0] << "\n");
	
	DebugM(2, "d2V/dxdy = " << b[32+i0] << "\n");
	DebugM(2, "d2V/dxdz = " << b[40+i0] << "\n");
	DebugM(2, "d2V/dydz = " << b[48+i0] << "\n");
	
	DebugM(2, "d3V/dxdydz = " << b[56+i0] << "\n" << endi);
    }
}


/********************/
/* GRIDFORCESUBGRID */
/********************/

GridforceSubGrid::GridforceSubGrid(GridforceGrid *parent_in) {
    parent = parent_in;
    generation = parent->generation + 1;
    GridforceGrid *tmp = parent;
    while (tmp->generation > 0) {
	tmp = ((GridforceSubGrid *)tmp)->parent;
    }
    maingrid = (GridforceMainGrid *)tmp;
    DebugM(4, "generation = " << generation << "\n" << endi);
}


void GridforceSubGrid::pack(MOStream *msg) const
{
    DebugM(4, "Packing subgrid\n" << endi);
    
    GridforceGrid::pack(msg);
    
    msg->put(sizeof(Tensor), (char*)&scale_dV);
    msg->put(sizeof(Tensor), (char*)&scale_d2V);
    msg->put(sizeof(float), (char*)&scale_d3V);
    
    msg->put(3*sizeof(int), (char*)pmin);
    msg->put(3*sizeof(int), (char*)pmax);
    msg->put(subgridIdx);
}


void GridforceSubGrid::unpack(MIStream *msg)
{
    DebugM(4, "Unpacking subgrid\n" << endi);
    
    GridforceGrid::unpack(msg);
    
    msg->get(sizeof(Tensor), (char*)&scale_dV);
    msg->get(sizeof(Tensor), (char*)&scale_d2V);
    msg->get(sizeof(float), (char*)&scale_d3V);
    
    msg->get(3*sizeof(int), (char*)pmin);
    msg->get(3*sizeof(int), (char*)pmax);
    msg->get(subgridIdx);
    
    DebugM(4, "size  = " << size << "\n" << endi);
    DebugM(4, "numSubgrids = " << numSubgrids << "\n" << endi);
    DebugM(4, "gapinv = " << gapinv[0] << " " << gapinv[2] << " " << gapinv[2] << " " << "\n" << endi);
    DebugM(4, "generation = " << generation << "\n" << endi);
}


void GridforceSubGrid::addToSubgridsFlat(void)
{
    DebugM(4, "addToSubgridsFlat() called, subgridIdx = " << subgridIdx << ", maingrid->numSubgrids = " << maingrid->numSubgrids << "\n" << endi);
    maingrid->subgrids_flat[subgridIdx-1] = this;
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->addToSubgridsFlat();
    }
}


void GridforceSubGrid::init1(char *potfilename, SimParameters *simParams, MGridforceParams *mgridParams)
{
    poten_fp = Fopen(potfilename, "r");
    if (!poten_fp) {
	NAMD_die("Problem reading grid force potential file");
    }
    
    // Read through previous grids in file
    for (int i = 0; i < subgridIdx; i++) {
	DebugM(4, "skipping grid " << i << "\n" << endi);
	
	// Skip comment lines
	char line[256];
	long int poten_offset;
	do {
	    poten_offset = ftell(poten_fp);
	    fgets(line, 256, poten_fp);
	} while (line[0] == '#');
	fseek(poten_fp, poten_offset, SEEK_SET);
	
	// Read header, ignoring most
	int tmp;
	
	fscanf(poten_fp, "object %*d class gridpositions counts %*d %*d %*d\n");
	fscanf(poten_fp, "origin %*lf %*lf %*lf\n");
	fscanf(poten_fp, "delta %*lf %*lf %*lf\n");
	fscanf(poten_fp, "delta %*lf %*lf %*lf\n");
	fscanf(poten_fp, "delta %*lf %*lf %*lf\n");
	fscanf(poten_fp, "object %*d class gridconnections counts %*lf %*lf %*lf\n");
	fscanf(poten_fp, "object %*d class array type double rank 0 items %d data follows\n", &tmp);
	
	// Skip data
	DebugM(4, "Skipping " << tmp << " data objects\n" << endi);
	for (int i = 0; i < tmp; i++) {
	    fscanf(poten_fp, "%*f");
	}
	fscanf(poten_fp, "\n");
	
	// Skip 'attribute's
	DebugM(3, "Skipping 'attribute' keywords...\n" << endi);
	char str[256];
	do {
	    poten_offset = ftell(poten_fp);
	    fscanf(poten_fp, "%s", str);
	    fgets(line, 256, poten_fp);
	    DebugM(4, "Read line " << str << " " << line << endi);
	} while (strcmp(str, "attribute") == 0);
	fseek(poten_fp, poten_offset, SEEK_SET);
	
	// Skip 'field' object
	DebugM(3, "Skipping 'field' object\n" << endi);
	fscanf(poten_fp, "object");
	int n;
	n = fscanf(poten_fp, "\"%[^\"]\" class field\n", str);
	if (n == 0) {
	    n = fscanf(poten_fp, "%d class field\n", &tmp);
	}
	
	if (n == 0) {
	    NAMD_die("Error reading gridforce grid! Could not find field object!\n");
	}
	
	// Skip 'component's
	DebugM(3, "Skipping 'component' keywords\n" << endi);
	do {
	    poten_offset = ftell(poten_fp);
	    fscanf(poten_fp, "%s", str);
	    fgets(line, 256, poten_fp);
	} while (strcmp(str, "component") == 0);
	fseek(poten_fp, poten_offset, SEEK_SET);
    }
    
    readHeader(simParams, mgridParams);
    
    for (int i = 0; i < 3; i++) {
	k[i] = k_nopad[i];	// subgrids aren't padded
    }
    
    // Make sure that each subgrid dimension is an integral
    // number of spanned supergrid cells. This is to ensure that no
    // supergrid nodes land in the middle of a subgrid, because in
    // this case forces could not be matched properly.
    for (int i = 0; i < 3; i++) {
	if ((k[i] - 1) % (pmax[i] - pmin[i] + 1) != 0) {
	    iout << (k[i] - 1) << " % " << (pmax[i] - pmin[i] + 1) << " != 0\n" << endi;
	    NAMD_die("Error reading gridforce grid! Subgrid dimensions must be an integral number spanned parent cells!");
	}
    }
    
    
    for (int i = 0; i < 3; i++) {
	if (parent->cont[i]) {
	    cont[i] = (pmin[i] == 0 && pmax[i] == parent->k[i]-2) ? TRUE : FALSE;
	    DebugM(3, "pmin[" << i << "] = " << pmin[i] << " pmax[" << i << "] = " << pmax[i] << " parent->k[" << i << "] = " << parent->k[i] << " cont[" << i << "] = " << cont[i] << "\n" << endi);
	} else {
	    cont[i] = false;
	    if (parent->generation == 0) {
		// Add one to pmin, pmax since parent got an extra gridpoint layer
		pmin[i] += GridforceMainGrid::border;
		pmax[i] += GridforceMainGrid::border;
	    }
	}		
    }
    
    DebugM(4, "pmin = " << pmin[0] << " " << pmin[1] << " " << pmin[2] << "\n" << endi);
    DebugM(4, "pmax = " << pmax[0] << " " << pmax[1] << " " << pmax[2] << "\n" << endi);
    
    Vector origin2 = parent->origin + parent->e * Position(pmin[0], pmin[1], pmin[2]);
    Vector escale, invscale;
    for (int i = 0; i < 3; i++) {
	escale[i] = double(pmax[i] - pmin[i] + 1)/(k[i]-1);
	invscale[i] = 1.0/escale[i];
	if (cont[i]) { pmax[i]++; }
    }
    Tensor e2 = tensorMult(parent->e, Tensor::diagonal(escale));
    
    // Check that lattice parameters agree with min and max numbers
    // from subgrid hierarchy.
    double TOL2 = 1e-4;	// Totally arbitrary
    if (pow(origin2.x-origin.x, 2) > TOL2 ||
	pow(origin2.y-origin.y, 2) > TOL2 ||
	pow(origin2.z-origin.z, 2) > TOL2 ||
	pow(e2.xx-e.xx, 2) > TOL2 ||
	pow(e2.xy-e.xy, 2) > TOL2 ||
	pow(e2.xz-e.xz, 2) > TOL2 ||
	pow(e2.yx-e.yx, 2) > TOL2 ||
	pow(e2.yy-e.yy, 2) > TOL2 ||
	pow(e2.yz-e.yz, 2) > TOL2 ||
	pow(e2.zx-e.zx, 2) > TOL2 ||
	pow(e2.zy-e.zy, 2) > TOL2 ||
	pow(e2.zz-e.zz, 2) > TOL2)
    {
	NAMD_die("Error reading gridforce grid! Subgrid lattice does not match!");
    }
    
    // Overwrite what was read from the header
    origin = origin2;
    e = e2;
    
    inv = tensorMult(Tensor::diagonal(invscale), parent->inv);
    for (int i = 0; i < 3; i++) {
	gap[i] = escale[i] * parent->gap[i];
	gapinv[i] = invscale[i] * parent->gapinv[i];
	offset[i] = parent->offset[i];
    }
    center = origin + e * 0.5 * Position(k[0], k[1], k[2]);
    
    DebugM(4, "origin = " << origin << "\n" << endi);
    DebugM(4, "e = " << e << "\n" << endi);
    DebugM(4, "inv = " << inv << "\n" << endi);
    DebugM(4, "gap = " << gap[0] << " " << gap[2] << " " << gap[2] << " " << "\n" << endi);
    DebugM(4, "gapinv = " << gapinv[0] << " " << gapinv[2] << " " << gapinv[2] << " " << "\n" << endi);
    DebugM(4, "numSubgrids = " << numSubgrids << "\n" << endi);
    DebugM(4, "k = " << k[0] << " " << k[1] << " " << k[2] << "\n" << endi);
    DebugM(4, "escale = " << escale << "\n" << endi);
    DebugM(4, "invscale = " << invscale << "\n" << endi);
    
    /*** Set members ***/
    size = k[0] * k[1] * k[2];
    dk[0] = k[1] * k[2];
    dk[1] = k[2];
    dk[2] = 1;
    
    scale_dV = Tensor::diagonal(escale);
    scale_d2V = Tensor::diagonal(Vector(escale.x*escale.y, escale.x*escale.z, escale.y*escale.z));
    scale_d3V = escale.x * escale.y * escale.z;
    
    DebugM(4, "scale_dV = " << scale_dV << "\n" << endi);
    DebugM(4, "scale_d2V = " << scale_d2V << "\n" << endi);
    DebugM(4, "scale_d3V = " << scale_d3V << "\n" << endi);
    
    // Call for each subgrid
    for (int i = 0; i < numSubgrids; i++) {
	subgrids[i]->init1(potfilename, simParams, mgridParams);
    }
}


void GridforceSubGrid::init2()
// On all: Receive grid params... Figure out what grid values are needed
{
    GridforceGrid::init2();
    
    // Fill in grid
    for (int ind = 0; ind < size; ind++) {
	if (isGridValCached(ind)) {
	    gridvals_needed[getGridBlockIdx(ind)] = true;
	}
    }
}


void GridforceSubGrid::init4(float *gridseg, int start_indx, int count)
{
    int i;
    for(i=0; i < count; i++) {
	int ind = start_indx + i;
	if (gridvals_needed[getGridBlockIdx(ind)]) {
	    addGridValue(ind,gridseg[i]);
	}
    }
}


void GridforceSubGrid::compute_b(float *b, int *inds, Vector gapscale) const
{
    for (int i0 = 0; i0 < 8; i0++) {
	int inds2[3];
	int ind2 = 0;
	int dk_hi[3], dk_lo[3];
	for (int i1 = 0; i1 < 3; i1++) {
	    inds2[i1] = grid_index_table[i0][i1][inds[i1]].inds2;
	    ind2 += grid_index_table[i0][i1][inds[i1]].inds2 * dk[i1];
	    dk_hi[i1] = grid_index_table[i0][i1][inds[i1]].dk_hi;
	    dk_lo[i1] = grid_index_table[i0][i1][inds[i1]].dk_lo;
	}
	
	float voff = 0.0;
	for (int i1 = 0; i1 < 3; i1++) {
	    // Deal with voltage offsets
	    if (inds[i1] == (k[i1]-1) && inds2[i1] == 0 && cont[i1]) {
		voff += offset[i1];
	    }
	}
	
 	DebugM(3, "inds2 = " << inds2[0] << " " << inds2[1] << " " << inds2[2] << "\n" << endi);
	
	float voffs[3];
	float dscales[3] = {0.5, 0.5, 0.5};
	for (int i1 = 0; i1 < 3; i1++) {
	    if (inds2[i1] == 0 && cont[i1]) {
		voffs[i1] = offset[i1];
		dscales[i1] = 1.0/(1.0 + gap[i1]) * 1.0/gapscale[i1];
	    }
	    else if (inds2[i1] == k[i1]-1 && cont[i1]) {
		voffs[i1] = offset[i1];
		dscales[i1] = 1.0/(1.0 + gap[i1]) * 1.0/gapscale[i1];
	    }
	    else {
		voffs[i1] = 0.0;
	    }
	}
	
	bool edge = false;
	for (int i1 = 0; i1 < 3; i1++) {
	    if (!cont[i1] && (inds2[i1] == 0 || inds2[i1] == k[i1]-1)) {
		edge = true;
	    }
	}
	
	if (inds2[2] == 0) {
	    DebugM(3, "cont = " << cont[0] << " " << cont[1] << " " << cont[2] << " dk_hi = " << dk_hi[0] << " " << dk_hi[1] << " " << dk_hi[2] << " dk_lo = " << dk_lo[0] << " " << dk_lo[1] << " " << dk_lo[2] << " dscales = " << dscales[0] << " " << dscales[1] << " " << dscales[2] << "\n" << endi);
	}
	
	if (edge) {
	    DebugM(2, "Edge!\n");
	    
	    // Must get derivatives from parent
	    Position pos = e * Vector(inds2[0], inds2[1], inds2[2]) + origin;	// Gridpoint position in realspace
	    Vector g = parent->inv * (pos - parent->origin);	// Gridpoint position in parent's gridspace
	    Vector dg;
	    int inds3[3];
	    
	    DebugM(2, "g = " << g << "\n" << endi);
	    
	    for (int i = 0; i < 3; i++) {
		inds3[i] = (int)floor(g[i]);
		dg[i] = g[i] - inds3[i];
	    }
    
	    float x[4], y[4], z[4];
	    x[0] = 1; y[0] = 1; z[0] = 1;
	    for (int j = 1; j < 4; j++) {
		x[j] = x[j-1] * dg.x;
		y[j] = y[j-1] * dg.y;
		z[j] = z[j-1] * dg.z;
		DebugM(1, "x[" << j << "] = " << x[j] << "\n" << endi);
		DebugM(1, "y[" << j << "] = " << y[j] << "\n" << endi);
		DebugM(1, "z[" << j << "] = " << z[j] << "\n" << endi);
	    }
	    
	    // Compute parent matrices
	    float b_parent[64];
	    parent->compute_b(b_parent, inds3, gapscale);
	    
	    float a_parent[64];
	    parent->compute_a(a_parent, b_parent);
	    
	    // Compute parent derivatives
	    float V = parent->compute_V(a_parent, x, y, z);
	    Vector dV = scale_dV * parent->compute_dV(a_parent, x, y, z);
	    Vector d2V = scale_d2V * parent->compute_d2V(a_parent, x, y, z);
	    float d3V = scale_d3V * parent->compute_d3V(a_parent, x, y, z);
	    
	    b[i0] = V;
	    b[8+i0] = dV[0];
	    b[16+i0] = dV[1];
	    b[24+i0] = dV[2];
	    b[32+i0] = d2V[0];
	    b[40+i0] = d2V[1];
	    b[48+i0] = d2V[2];
	    b[56+i0] = d3V;
	} else {
	    b[i0] = get_grid(ind2) + voff;	// V
	
	    b[8+i0] = dscales[0] * (get_grid(ind2 + dk_hi[0]) - get_grid(ind2 - dk_lo[0]) + voffs[0]);	//  dV/dx
	    b[16+i0] = dscales[1] * (get_grid(ind2 + dk_hi[1]) - get_grid(ind2 - dk_lo[1]) + voffs[1]);	//  dV/dy
	    b[24+i0] = dscales[2] * (get_grid(ind2 + dk_hi[2]) - get_grid(ind2 - dk_lo[2]) + voffs[2]);	//  dV/dz
	    
	    b[32+i0] = dscales[0] * dscales[1]
		* (get_grid(ind2 + dk_hi[0] + dk_hi[1]) - get_grid(ind2 - dk_lo[0] + dk_hi[1])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[1]) + get_grid(ind2 - dk_lo[0] - dk_lo[1]));	//  d2V/dxdy
	    b[40+i0] = dscales[0] * dscales[2]
		* (get_grid(ind2 + dk_hi[0] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] + dk_hi[2])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[2]) + get_grid(ind2 - dk_lo[0] - dk_lo[2]));	//  d2V/dxdz
	    b[48+i0] = dscales[1] * dscales[2]
		* (get_grid(ind2 + dk_hi[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[1] + dk_hi[2])
		   - get_grid(ind2 + dk_hi[1] - dk_lo[2]) + get_grid(ind2 - dk_lo[1] - dk_lo[2]));	//  d2V/dydz
	
	    b[56+i0] = dscales[0] * dscales[1] * dscales[2]					// d3V/dxdydz
		* (get_grid(ind2 + dk_hi[0] + dk_hi[1] + dk_hi[2]) - get_grid(ind2 + dk_hi[0]+ dk_hi[1] - dk_lo[2])
		   - get_grid(ind2 + dk_hi[0] - dk_lo[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] + dk_hi[1] + dk_hi[2])
		   + get_grid(ind2 + dk_hi[0] - dk_lo[1] - dk_lo[2]) + get_grid(ind2 - dk_lo[0] + dk_hi[1] - dk_lo[2])
		   + get_grid(ind2 - dk_lo[0] - dk_lo[1] + dk_hi[2]) - get_grid(ind2 - dk_lo[0] - dk_lo[1] - dk_lo[2]));
	}
	
	if (inds2[0] == 1 && inds2[1] == 1 && inds2[2] == 0) {
	    DebugM(2, "Sub V = " << b[i0] << "\n");
	
	    DebugM(2, "Sub dV/dx = " << b[8+i0] << "\n");
	    DebugM(2, "Sub dV/dy = " << b[16+i0] << "\n");
	    DebugM(2, "Sub dV/dz = " << b[24+i0] << "\n");
	
	    DebugM(2, "Sub d2V/dxdy = " << b[32+i0] << "\n");
	    DebugM(2, "Sub d2V/dxdz = " << b[40+i0] << "\n");
	    DebugM(2, "Sub d2V/dydz = " << b[48+i0] << "\n");
	
	    DebugM(2, "Sub d3V/dxdydz = " << b[56+i0] << "\n" << endi);
	}
    }
}

