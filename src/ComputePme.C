/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifdef NAMD_FFTW
#include <fftw.h>
#include <rfftw.h>
#endif

#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputePme.h"
#include "ComputePmeMgr.decl.h"
#include "PmeRealSpace.h"
#include "PmeKSpace.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "varsizemsg.h"
#include "Random.h"

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif


class PmeGridMsg : public CMessage_PmeGridMsg {
public:

  int sourceNode;
  Lattice lattice;
  PmeReduction evir;
  int start;
  int len;
  int zlistlen;
  int *zlist;
  char *fgrid;
  float *qgrid;
};

class PmeTransMsg : public CMessage_PmeTransMsg {
public:

  int sourceNode;
  Lattice lattice;
  int x_start;
  int nx;
  double *qgrid;
};


class PmeUntransMsg : public CMessage_PmeUntransMsg {
public:

  int sourceNode;
  PmeReduction evir;
  int y_start;
  int ny;
  double *qgrid;

};


struct LocalPmeInfo {
  int nx, x_start;
  int ny_after_transpose, y_start_after_transpose;
};

class ComputePmeMgr : public BOCclass {
public:
  ComputePmeMgr();
  ~ComputePmeMgr();

  void initialize(CkQdMsg*);

  void sendGrid(void);
  void recvGrid(PmeGridMsg *);
  void gridCalc1(void);
  void sendTrans(void);
  void recvTrans(PmeTransMsg *);
  void gridCalc2(void);
  void sendUntrans(void);
  void recvUntrans(PmeUntransMsg *);
  void gridCalc3(void);
  void sendUngrid(void);
  void recvUngrid(PmeGridMsg *);
  void ungridCalc(void);

  void setCompute(ComputePme *c) { pmeCompute = c; }

private:
  CProxy_ComputePmeMgr pmeProxy;
  ComputePme *pmeCompute;
  PmeGrid myGrid;
  Lattice lattice;
  PmeKSpace *myKSpace;
  double *qgrid;
  double *kgrid;

#ifdef NAMD_FFTW
  fftw_plan forward_plan_x, backward_plan_x;
  rfftwnd_plan forward_plan_yz, backward_plan_yz;
  fftw_complex *work;
#else
  double *work;
#endif

  int fepOn, lesOn, lesFactor, numGrids;

  LocalPmeInfo *localInfo;
  int qgrid_size;
  int qgrid_start;
  int qgrid_len;
  int fgrid_start;
  int fgrid_len;

  int numSources;
  int numGridPes;
  int numTransPes;
  int numDestRecipPes;
  int myGridPe;
  int myTransPe;
  int *gridPeMap;
  int *transPeMap;
  int *recipPeDest;
  int *gridPeOrder;
  int *transPeOrder;
  int grid_count;
  int trans_count;
  int untrans_count;
  int ungrid_count;
  PmeGridMsg **gridmsg_reuse;
  PmeReduction recip_evir;
  PmeReduction recip_evir2;
};

ComputePmeMgr::ComputePmeMgr() : pmeProxy(thisgroup), pmeCompute(0) {
  CpvAccess(BOCclass_group).computePmeMgr = thisgroup;
  myKSpace = 0;
  localInfo = new LocalPmeInfo[CkNumPes()];
  gridPeMap = new int[CkNumPes()];
  transPeMap = new int[CkNumPes()];
  recipPeDest = new int[CkNumPes()];
  gridPeOrder = new int[CkNumPes()];
  transPeOrder = new int[CkNumPes()];
  qgrid = 0;
  kgrid = 0;
  work = 0;
  grid_count = 0;
  trans_count = 0;
  untrans_count = 0;
  ungrid_count = 0;
  gridmsg_reuse= new PmeGridMsg*[CkNumPes()];
}

void ComputePmeMgr::initialize(CkQdMsg *msg) {
  delete msg;

  SimParameters *simParams = Node::Object()->simParameters;

  fepOn = simParams->fepOn;
  numGrids = fepOn ? 2 : 1;
  lesOn = simParams->lesOn;
  if ( lesOn ) {
    lesFactor = simParams->lesFactor;
    numGrids = lesFactor;
  }

  {  // decide how many pes to use for reciprocal sum

    // rules based on work available
    int minslices = 1;
    int dimx = simParams->PMEGridSizeX;
    int nrpx = ( dimx + minslices - 1 ) / minslices;
    int dimy = simParams->PMEGridSizeY;
    int nrpy = ( dimy + minslices - 1 ) / minslices;

    // rules based on processors available
    int nrpp = CkNumPes();
    // if ( nrpp > 32 ) nrpp = 32;  // cap to limit messages
    if ( nrpp < nrpx ) nrpx = nrpp;
    if ( nrpp < nrpy ) nrpy = nrpp;

    // user override
    int nrps = simParams->PMEProcessors;
    if ( nrps > CkNumPes() ) nrps = CkNumPes();
    if ( nrps > 0 ) nrpx = nrps;
    if ( nrps > 0 ) nrpy = nrps;

    // make sure there aren't any totally empty processors
    int bx = ( dimx + nrpx - 1 ) / nrpx;
    nrpx = ( dimx + bx - 1 ) / bx;
    int by = ( dimy + nrpy - 1 ) / nrpy;
    nrpy = ( dimy + by - 1 ) / by;
    if ( bx != ( dimx + nrpx - 1 ) / nrpx )
      NAMD_bug("Error in selecting number of PME processors.");
    if ( by != ( dimy + nrpy - 1 ) / nrpy )
      NAMD_bug("Error in selecting number of PME processors.");

    numGridPes = nrpx;
    numTransPes = nrpy;
  }
  if ( ! CkMyPe() ) {
    iout << iINFO << "PME using " << numGridPes << " and " << numTransPes <<
      " processors for FFT and reciprocal sum.\n" << endi;
  }
  { // generate random orderings for grid and trans messages
    for ( int i = 0; i < numGridPes; ++i ) {
      gridPeOrder[i] = i;
    }
    for ( int i = 0; i < numTransPes; ++i ) {
      transPeOrder[i] = i;
    }
    Random rand(CkMyPe());
    rand.reorder(gridPeOrder,numGridPes);
    rand.reorder(transPeOrder,numTransPes);
  }

  {  // decide which pes to use by bit reversal
    int i;
    int ncpus = CkNumPes();

    // find next highest power of two
    int npow2 = 1;  int nbits = 0;
    while ( npow2 < ncpus ) { npow2 *= 2; nbits += 1; }

    // build bit reversal sequence
    SortableResizeArray<int> seq(ncpus);
    SortableResizeArray<int> seq2(ncpus);
    i = 0;
    for ( int icpu=0; icpu<ncpus; ++icpu ) {
      int ri;
      for ( ri = ncpus; ri >= ncpus; ++i ) {
        ri = 0;
        int pow2 = 1;
        int rpow2 = npow2 / 2;
        for ( int j=0; j<nbits; ++j ) {
          ri += rpow2 * ( ( i / pow2 ) % 2 );
          pow2 *= 2;  rpow2 /= 2;
        }
      }
      seq[icpu] = ri;
      seq2[icpu] = ri;
    }

    // extract and sort PME locations
    for ( i=0; i<numGridPes; ++i ) {
      seq[i] = seq[ncpus - numGridPes + i];
    }
    seq.resize(numGridPes);
    seq.sort();
    if ( ncpus > numTransPes ) {
      seq2.del(0);  // node 0 should be first in list
    }
    seq2.resize(numTransPes);
    seq2.sort();

    myGridPe = -1;
    for ( i=0; i<numGridPes; ++i ) {
      gridPeMap[i] = seq[i];
      if ( gridPeMap[i] == CkMyPe() ) myGridPe = i;
    }
    myTransPe = -1;
    for ( i=0; i<numTransPes; ++i ) {
      transPeMap[i] = seq2[i];
      if ( transPeMap[i] == CkMyPe() ) myTransPe = i;
    }
  }

  if ( ! CkMyPe() ) {
    iout << iINFO << "PME GRID LOCATIONS:";
    int i;
    for ( i=0; i<numGridPes && i<10; ++i ) {
      iout << " " << gridPeMap[i];
    }
    if ( i < numGridPes ) iout << " ...";
    iout << "\n" << endi;
    iout << iINFO << "PME TRANS LOCATIONS:";
    for ( i=0; i<numTransPes && i<10; ++i ) {
      iout << " " << transPeMap[i];
    }
    if ( i < numTransPes ) iout << " ...";
    iout << "\n" << endi;
  }

  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
  myGrid.block1 = ( myGrid.K1 + numGridPes - 1 ) / numGridPes;
  myGrid.block2 = ( myGrid.K2 + numTransPes - 1 ) / numTransPes;

  int nx = 0;
  for ( int pe = 0; pe < numGridPes; ++pe ) {
    localInfo[pe].x_start = nx;
    nx += myGrid.block1;
    if ( nx > myGrid.K1 ) nx = myGrid.K1;
    localInfo[pe].nx = nx - localInfo[pe].x_start;
  }
  int ny = 0;
  for ( int pe = 0; pe < numTransPes; ++pe ) {
    localInfo[pe].y_start_after_transpose = ny;
    ny += myGrid.block2;
    if ( ny > myGrid.K2 ) ny = myGrid.K2;
    localInfo[pe].ny_after_transpose =
			ny - localInfo[pe].y_start_after_transpose;
  }

  {  // decide how many pes this node exchanges charges with

  PatchMap *patchMap = PatchMap::Object();
  Lattice lattice = simParams->lattice;
  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal cutoff = simParams->cutoff;
  BigReal patchdim = simParams->patchDimension;
  int numPatches = patchMap->numPatches();
  int numNodes = CkNumPes();
  int *source_flags = new int[numNodes];
  int node;
  for ( node=0; node<numNodes; ++node ) {
    source_flags[node] = 0;
    recipPeDest[node] = 0;
  }

  // // make sure that we don't get ahead of ourselves on this node
  // if ( CkMyPe() < numPatches && myRecipPe >= 0 ) {
  //   source_flags[CkMyPe()] = 1;
  //   recipPeDest[myRecipPe] = 1;
  // }

  for ( int pid=0; pid < numPatches; ++pid ) {
    int pnode = patchMap->node(pid);
    BigReal minx = patchMap->min_a(pid);
    BigReal maxx = patchMap->max_a(pid);
    BigReal margina = 0.5 * ( patchdim - cutoff ) / sysdima;
    // min1 (max1) is smallest (largest) grid line for this patch
    int min1 = ((int) floor(myGrid.K1 * (minx - margina))) - myGrid.order + 1;
    int max1 = ((int) floor(myGrid.K1 * (maxx + margina)));
    for ( int i=min1; i<=max1; ++i ) {
      int ix = i;
      while ( ix >= myGrid.K1 ) ix -= myGrid.K1;
      while ( ix < 0 ) ix += myGrid.K1;
      // set source_flags[pnode] if this patch sends to our node
      if ( myGridPe >= 0 && ix >= localInfo[myGridPe].x_start &&
           ix < localInfo[myGridPe].x_start + localInfo[myGridPe].nx ) {
        source_flags[pnode] = 1;
      }
      // set dest_flags[] for node that our patch sends to
      if ( pnode == CkMyPe() ) {
        recipPeDest[ix / myGrid.block1] = 1;
      }
    }
  }

  numSources = 0;
  numDestRecipPes = 0;
  for ( node=0; node<numNodes; ++node ) {
    if ( source_flags[node] ) ++numSources;
    if ( recipPeDest[node] ) ++numDestRecipPes;
  }

  delete [] source_flags;

  // CkPrintf("PME on node %d has %d sources and %d destinations\n",
  //           CkMyPe(), numSources, numDestRecipPes);

  }  // decide how many pes this node exchanges charges with (end)

  ungrid_count = numDestRecipPes;

  if ( myGridPe < 0 && myTransPe < 0 ) return;
  // the following only for nodes doing reciprocal sum

  if ( myTransPe >= 0 ) {
  int k2_start = localInfo[myTransPe].y_start_after_transpose;
  int k2_end = k2_start + localInfo[myTransPe].ny_after_transpose;
  myKSpace = new PmeKSpace(myGrid, k2_start, k2_end);
  }

  int local_size = myGrid.block1 * myGrid.K2 * myGrid.dim3;
  int local_size_2 = myGrid.block2 * myGrid.K1 * myGrid.dim3;
  if ( local_size < local_size_2 ) local_size = local_size_2;
  qgrid = new double[local_size*numGrids];
  if ( numGridPes > 1 || numTransPes > 1 ) {
    kgrid = new double[local_size*numGrids];
  } else {
    kgrid = qgrid;
  }
  qgrid_size = local_size;

  if ( myGridPe >= 0 ) {
  qgrid_start = localInfo[myGridPe].x_start * myGrid.K2 * myGrid.dim3;
  qgrid_len = localInfo[myGridPe].nx * myGrid.K2 * myGrid.dim3;
  fgrid_start = localInfo[myGridPe].x_start * myGrid.K2;
  fgrid_len = localInfo[myGridPe].nx * myGrid.K2;
  }

  int n[3]; n[0] = myGrid.K1; n[1] = myGrid.K2; n[2] = myGrid.K3;

#ifdef NAMD_FFTW
  work = new fftw_complex[n[0]];

  if ( ! CkMyPe() ) iout << iINFO << "Optimizing 4 FFT steps.  1..." << endi;
  if ( myGridPe >= 0 ) {
  forward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_REAL_TO_COMPLEX,
	FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
  }
  if ( ! CkMyPe() ) iout << " 2..." << endi;
  if ( myTransPe >= 0 ) {
  forward_plan_x = fftw_create_plan_specific(n[0], FFTW_REAL_TO_COMPLEX,
	FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) kgrid,
	localInfo[myTransPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
  }
  if ( ! CkMyPe() ) iout << " 3..." << endi;
  if ( myTransPe >= 0 ) {
  backward_plan_x = fftw_create_plan_specific(n[0], FFTW_COMPLEX_TO_REAL,
	FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) kgrid,
	localInfo[myTransPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
  }
  if ( ! CkMyPe() ) iout << " 4..." << endi;
  if ( myGridPe >= 0 ) {
  backward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_COMPLEX_TO_REAL,
	FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
  }
  if ( ! CkMyPe() ) iout << "   Done.\n" << endi;
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  if ( myGridPe >= 0 && numSources == 0 )
		NAMD_bug("PME grid elements exist without sources.");
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_size * numGrids * sizeof(double) );
  trans_count = numGridPes;
}

ComputePmeMgr::~ComputePmeMgr() {
  delete myKSpace;
  delete [] localInfo;
  delete [] gridPeMap;
  delete [] transPeMap;
  delete [] recipPeDest;
  delete [] gridPeOrder;
  delete [] transPeOrder;
  delete [] qgrid;
  if ( kgrid != qgrid ) delete [] kgrid;
  delete [] work;
  delete [] gridmsg_reuse;
}

void ComputePmeMgr::sendGrid(void) {
  pmeCompute->sendData(numGridPes,gridPeOrder,recipPeDest,gridPeMap);
}

void ComputePmeMgr::recvGrid(PmeGridMsg *msg) {
  // CkPrintf("recvGrid from %d on Pe(%d)\n",msg->sourceNode,CkMyPe());
  if ( grid_count == 0 ) {
    NAMD_bug("Message order failure in ComputePmeMgr::recvGrid\n");
  }
  if ( grid_count == numSources ) {
    lattice = msg->lattice;
  }

  int zdim = myGrid.dim3;
  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  float *qmsg = msg->qgrid;
  for ( int g=0; g<numGrids; ++g ) {
    char *f = msg->fgrid + fgrid_len * g;
    double *q = qgrid + qgrid_size * g;
    for ( int i=0; i<fgrid_len; ++i ) {
      if ( f[i] ) {
        for ( int k=0; k<zlistlen; ++k ) {
          q[zlist[k]] += *(qmsg++);
        }
      }
      q += zdim;
    }
  }

  gridmsg_reuse[numSources-grid_count] = msg;
  --grid_count;

  if ( grid_count == 0 ) {
#if CHARM_VERSION > 050402
    pmeProxy[CkMyPe()].gridCalc1();
#else
    pmeProxy.gridCalc1(CkMyPe());
#endif
  }
}

void ComputePmeMgr::gridCalc1(void) {
  // CkPrintf("gridCalc1 on Pe(%d)\n",CkMyPe());

#ifdef NAMD_FFTW
  for ( int g=0; g<numGrids; ++g ) {
    rfftwnd_real_to_complex(forward_plan_yz, localInfo[myGridPe].nx,
	qgrid + qgrid_size * g, 1, myGrid.dim2 * myGrid.dim3, 0, 0, 0);
  }
#endif

#if CHARM_VERSION > 050402
  pmeProxy[CkMyPe()].sendTrans();
#else
  pmeProxy.sendTrans(CkMyPe());
#endif
}

void ComputePmeMgr::sendTrans(void) {

  // send data for transpose
  int zdim = myGrid.dim3;
  int nx = localInfo[myGridPe].nx;
  int x_start = localInfo[myGridPe].x_start;
  int slicelen = myGrid.K2 * zdim;
  for (int j=0; j<numTransPes; j++) {
    int pe = transPeOrder[j];  // different order on each node
    LocalPmeInfo &li = localInfo[pe];
    int cpylen = li.ny_after_transpose * zdim;
    PmeTransMsg *newmsg = new (nx * cpylen * numGrids,0) PmeTransMsg;
    newmsg->sourceNode = myGridPe;
    newmsg->lattice = lattice;
    newmsg->x_start = x_start;
    newmsg->nx = nx;
    for ( int g=0; g<numGrids; ++g ) {
      double *q = qgrid + qgrid_size * g + li.y_start_after_transpose * zdim;
      double *qmsg = newmsg->qgrid + nx * cpylen * g;
      for ( int x = 0; x < nx; ++x ) {
        memcpy((void*)qmsg, (void*)q, cpylen*sizeof(double));
        q += slicelen;
        qmsg += cpylen;
      }
    }
#if CHARM_VERSION > 050402
    pmeProxy[transPeMap[pe]].recvTrans(newmsg);
#else
    pmeProxy.recvTrans(newmsg,transPeMap[pe]);
#endif
  }

  untrans_count = numTransPes;
}

void ComputePmeMgr::recvTrans(PmeTransMsg *msg) {
  // CkPrintf("recvTrans on Pe(%d)\n",CkMyPe());
  if ( trans_count == numGridPes ) {
    lattice = msg->lattice;
  }

  int zdim = myGrid.dim3;
  // int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;
  int x_start = msg->x_start;
  int nx = msg->nx;
  for ( int g=0; g<numGrids; ++g ) {
    memcpy((void*)(kgrid + qgrid_size * g + x_start*ny*zdim),
	(void*)(msg->qgrid + nx*ny*zdim*g), nx*ny*zdim*sizeof(double));
  }

  delete msg;
  --trans_count;

  if ( trans_count == 0 ) {
#if CHARM_VERSION > 050402
    pmeProxy[CkMyPe()].gridCalc2();
#else
    pmeProxy.gridCalc2(CkMyPe());
#endif
  }
}

void ComputePmeMgr::gridCalc2(void) {
  // CkPrintf("gridCalc2 on Pe(%d)\n",CkMyPe());

  int zdim = myGrid.dim3;
  // int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;

  for ( int g=0; g<numGrids; ++g ) {
    // finish forward FFT (x dimension)
#ifdef NAMD_FFTW
    fftw(forward_plan_x, ny * zdim / 2, (fftw_complex *)(kgrid+qgrid_size*g),
	ny * zdim / 2, 1, work, 1, 0);
#endif

    // reciprocal space portion of PME
    BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
    recip_evir2[7*g] = myKSpace->compute_energy(kgrid+qgrid_size*g,
			lattice, ewaldcof, &(recip_evir2[7*g+1]));
    // CkPrintf("Ewald reciprocal energy = %f\n", recip_evir2[7*g]);

    // start backward FFT (x dimension)
#ifdef NAMD_FFTW
    fftw(backward_plan_x, ny * zdim / 2, (fftw_complex *)(kgrid+qgrid_size*g),
	ny * zdim / 2, 1, work, 1, 0);
#endif
  }

#if CHARM_VERSION > 050402
  pmeProxy[CkMyPe()].sendUntrans();
#else
  pmeProxy.sendUntrans(CkMyPe());
#endif
}

void ComputePmeMgr::sendUntrans(void) {

  int zdim = myGrid.dim3;
  int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;

  // send data for reverse transpose
  for (int j=0; j<numGridPes; j++) {
    int pe = gridPeOrder[j];  // different order on each node
    LocalPmeInfo &li = localInfo[pe];
    int x_start =li.x_start;
    int nx = li.nx;
    PmeUntransMsg *newmsg = new (nx*ny*zdim*numGrids,0) PmeUntransMsg;
    newmsg->sourceNode = myTransPe;
    if ( j == 0 ) {  // only need these once
      newmsg->evir = recip_evir2;
    } else {
      newmsg->evir = 0.;
    }
    newmsg->y_start = y_start;
    newmsg->ny = ny;
    for ( int g=0; g<numGrids; ++g ) {
      memcpy((void*)(newmsg->qgrid+nx*ny*zdim*g),
		(void*)(kgrid + qgrid_size*g + x_start*ny*zdim),
		nx*ny*zdim*sizeof(double));
    }
#if CHARM_VERSION > 050402
    pmeProxy[gridPeMap[pe]].recvUntrans(newmsg);
#else
    pmeProxy.recvUntrans(newmsg,gridPeMap[pe]);
#endif
  }

  trans_count = numGridPes;
}

void ComputePmeMgr::recvUntrans(PmeUntransMsg *msg) {
  // CkPrintf("recvUntrans on Pe(%d)\n",CkMyPe());
  if ( untrans_count == numTransPes ) {
    recip_evir = 0.;
  }

  recip_evir += msg->evir;

  int zdim = myGrid.dim3;
  // int x_start = localInfo[myGridPe].x_start;
  int nx = localInfo[myGridPe].nx;
  int y_start = msg->y_start;
  int ny = msg->ny;
  int slicelen = myGrid.K2 * zdim;
  int cpylen = ny * zdim;
  for ( int g=0; g<numGrids; ++g ) {
    double *q = qgrid + qgrid_size * g + y_start * zdim;
    double *qmsg = msg->qgrid + nx * cpylen * g;
    for ( int x = 0; x < nx; ++x ) {
      memcpy((void*)q, (void*)qmsg, cpylen*sizeof(double));
      q += slicelen;
      qmsg += cpylen;
    }
  }

  delete msg;
  --untrans_count;

  if ( untrans_count == 0 ) {
#if CHARM_VERSION > 050402
    pmeProxy[CkMyPe()].gridCalc3();
#else
    pmeProxy.gridCalc3(CkMyPe());
#endif
  }
}

void ComputePmeMgr::gridCalc3(void) {
  // CkPrintf("gridCalc3 on Pe(%d)\n",CkMyPe());

  // finish backward FFT
#ifdef NAMD_FFTW
  for ( int g=0; g<numGrids; ++g ) {
    rfftwnd_complex_to_real(backward_plan_yz, localInfo[myGridPe].nx,
	(fftw_complex *) (qgrid + qgrid_size * g),
	1, myGrid.dim2 * myGrid.dim3 / 2, 0, 0, 0);
  }
#endif

#if CHARM_VERSION > 050402
  pmeProxy[CkMyPe()].sendUngrid();
#else
  pmeProxy.sendUngrid(CkMyPe());
#endif
}

void ComputePmeMgr::sendUngrid(void) {

  for ( int j=0; j<numSources; ++j ) {
    // int msglen = qgrid_len;
    PmeGridMsg *newmsg = gridmsg_reuse[j];
    int pe = newmsg->sourceNode;
    if ( j == 0 ) {  // only need these once
      newmsg->evir = recip_evir;
    } else {
      newmsg->evir = 0.;
    }
    int zdim = myGrid.dim3;
    int flen = newmsg->len;
    int fstart = newmsg->start;
    int zlistlen = newmsg->zlistlen;
    int *zlist = newmsg->zlist;
    float *qmsg = newmsg->qgrid;
    for ( int g=0; g<numGrids; ++g ) {
      char *f = newmsg->fgrid + fgrid_len * g;
      double *q = qgrid + qgrid_size * g + (fstart-fgrid_start) * zdim;
      for ( int i=0; i<flen; ++i ) {
        if ( f[i] ) {
          for ( int k=0; k<zlistlen; ++k ) {
            *(qmsg++) = q[zlist[k]];
          }
        }
        q += zdim;
      }
    }
    newmsg->sourceNode = myGridPe;

#if CHARM_VERSION > 050402
    pmeProxy[pe].recvUngrid(newmsg);
#else
    pmeProxy.recvUngrid(newmsg,pe);
#endif
  }
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_size * numGrids * sizeof(double) );
}

void ComputePmeMgr::recvUngrid(PmeGridMsg *msg) {
  // CkPrintf("recvUngrid on Pe(%d)\n",CkMyPe());
  if ( ungrid_count == 0 ) {
    NAMD_bug("Message order failure in ComputePmeMgr::recvUngrid\n");
  }

  pmeCompute->copyResults(msg);
  delete msg;
  --ungrid_count;

  if ( ungrid_count == 0 ) {
#if CHARM_VERSION > 050402
    pmeProxy[CkMyPe()].ungridCalc();
#else
    pmeProxy.ungridCalc(CkMyPe());
#endif
  }
}

void ComputePmeMgr::ungridCalc(void) {
  // CkPrintf("ungridCalc on Pe(%d)\n",CkMyPe());

  pmeCompute->ungridForces();

  ungrid_count = numDestRecipPes;
}


static void scale_coordinates(PmeParticle p[], int N, Lattice lattice, PmeGrid grid) {
  Vector origin = lattice.origin();
  Vector recip1 = lattice.a_r();
  Vector recip2 = lattice.b_r();
  Vector recip3 = lattice.c_r();
  double ox = origin.x;
  double oy = origin.y;
  double oz = origin.z;
  double r1x = recip1.x;
  double r1y = recip1.y;
  double r1z = recip1.z;
  double r2x = recip2.x;
  double r2y = recip2.y;
  double r2z = recip2.z;
  double r3x = recip3.x;
  double r3y = recip3.y;
  double r3z = recip3.z;
  int K1 = grid.K1;
  int K2 = grid.K2;
  int K3 = grid.K3;

  for (int i=0; i<N; i++) {
    double px = p[i].x - ox;
    double py = p[i].y - oy;
    double pz = p[i].z - oz;
    double sx = px*r1x + py*r1y + pz*r1z;
    double sy = px*r2x + py*r2y + pz*r2z;
    double sz = px*r3x + py*r3y + pz*r3z;
    p[i].x = K1 * ( sx - floor(sx) );
    p[i].y = K2 * ( sy - floor(sy) );
    p[i].z = K3 * ( sz - floor(sz) );
    //  Check for rare rounding condition where K * ( 1 - epsilon ) == K
    //  which was observed with g++ on Intel x86 architecture.
    if ( p[i].x == K1 ) p[i].x = 0;
    if ( p[i].y == K2 ) p[i].y = 0;
    if ( p[i].z == K3 ) p[i].z = 0;
  }
}

static void scale_forces(Vector f[], int N, Lattice lattice) {
  Vector recip1 = lattice.a_r();
  Vector recip2 = lattice.b_r();
  Vector recip3 = lattice.c_r();
  double r1x = recip1.x;
  double r1y = recip1.y;
  double r1z = recip1.z;
  double r2x = recip2.x;
  double r2y = recip2.y;
  double r2z = recip2.z;
  double r3x = recip3.x;
  double r3y = recip3.y;
  double r3z = recip3.z;

  for (int i=0; i<N; i++) {
    double f1 = f[i].x;
    double f2 = f[i].y;
    double f3 = f[i].z;
    f[i].x = f1*r1x + f2*r2x + f3*r3x;
    f[i].y = f1*r1y + f2*r2y + f3*r3y;
    f[i].z = f1*r1z + f2*r2z + f3*r3z;
  }
}
 

ComputePme::ComputePme(ComputeID c) :
  ComputeHomePatches(c)
{
  DebugM(4,"ComputePme created.\n");

  CProxy_ComputePmeMgr::ckLocalBranch(
	CpvAccess(BOCclass_group).computePmeMgr)->setCompute(this);

  useAvgPositions = 1;

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  SimParameters *simParams = Node::Object()->simParameters;

  fepOn = simParams->fepOn;
  numGrids = fepOn ? 2 : 1;
  lesOn = simParams->lesOn;
  if ( lesOn ) {
    lesFactor = simParams->lesFactor;
    numGrids = lesFactor;
  }

  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
  qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  fsize = myGrid.K1 * myGrid.dim2;
  q_arr = new double*[fsize*numGrids];
  memset( (void*) q_arr, 0, fsize*numGrids * sizeof(double*) );
  f_arr = new char[fsize*numGrids];
  fz_arr = new char[myGrid.K3];
}

ComputePme::~ComputePme()
{
  for (int i=0; i<fsize*numGrids; ++i) {
    if ( q_arr[i] ) {
      delete [] q_arr[i];
    }
  }
  delete [] q_arr;
  delete [] f_arr;
  delete [] fz_arr;
}

void ComputePme::doWork()
{
  DebugM(4,"Entering ComputePme::doWork().\n");

  ResizeArrayIter<PatchElem> ap(patchList);

  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics )
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&r);
    }
    reduction->submit();
    return;
  }

  // allocate storage
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  Lattice lattice = patchList[0].p->flags.lattice;

  localData = new PmeParticle[numLocalAtoms*(numGrids+(numGrids>1?1:0))];
  localPartition = new char[numLocalAtoms];

  int g;
  for ( g=0; g<numGrids; ++g ) {
    localGridData[g] = localData + numLocalAtoms*(g+1);
  }

  // get positions and charges
  PmeParticle * data_ptr = localData;
  char * part_ptr = localPartition;
  const BigReal coloumb_sqrt = sqrt( COLOUMB * ComputeNonbondedUtil::scaling
				* ComputeNonbondedUtil::dielectric_1 );

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    CompAtom *x = (*ap).positionBox->open();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      data_ptr->x = x[i].position.x;
      data_ptr->y = x[i].position.y;
      data_ptr->z = x[i].position.z;
      data_ptr->cg = coloumb_sqrt * x[i].charge;
      ++data_ptr;
      *part_ptr = x[i].partition;
      ++part_ptr;
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
  }

  // copy to other grids if needed
  if ( fepOn || lesOn ) {
    for ( g=0; g<numGrids; ++g ) {
      PmeParticle *lgd = localGridData[g];
      int nga = 0;
      for(int i=0; i<numLocalAtoms; ++i) {
        if ( localPartition[i] == 0 || localPartition[i] == (g+1) ) {
          lgd[nga++] = localData[i];
        }
      }
      numGridAtoms[g] = nga;
    }
  } else {
    localGridData[0] = localData;
    numGridAtoms[0] = numLocalAtoms;
  }

  evir = 0;
  memset( (void*) fz_arr, 0, myGrid.K3 * sizeof(char) );

  // calculate self energy
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  for ( g=0; g<numGrids; ++g ) {
    BigReal selfEnergy = 0;
    data_ptr = localGridData[g];
    int i;
    for(i=0; i<numGridAtoms[g]; ++i)
    {
      selfEnergy += data_ptr->cg * data_ptr->cg;
      ++data_ptr;
    }
    selfEnergy *= -1. * ewaldcof / SQRT_PI;
    evir[7*g] += selfEnergy;

    double **q = q_arr + g*fsize;
    for (i=0; i<fsize; ++i) {
      if ( q[i] ) {
        memset( (void*) (q[i]), 0, myGrid.dim3 * sizeof(double) );
      }
    }

    char *f = f_arr + g*fsize;
    memset( (void*) f, 0, fsize * sizeof(char) );
    myRealSpace[g] = new PmeRealSpace(myGrid,numGridAtoms[g]);
    scale_coordinates(localGridData[g], numGridAtoms[g], lattice, myGrid);
    myRealSpace[g]->fill_charges(q, f, fz_arr, localGridData[g]);
  }

  CProxy_ComputePmeMgr pmeProxy(CpvAccess(BOCclass_group).computePmeMgr);
#if CHARM_VERSION > 050402
  pmeProxy[CkMyPe()].sendGrid();
#else
  pmeProxy.sendGrid(CkMyPe());
#endif
}

void ComputePme::sendData(int numRecipPes, int *recipPeOrder,
				int *recipPeDest, int *gridPeMap) {

  // iout << "Sending charge grid for " << numLocalAtoms << " atoms to FFT on " << iPE << ".\n" << endi;

  myGrid.block1 = ( myGrid.K1 + numRecipPes - 1 ) / numRecipPes;
  myGrid.block2 = ( myGrid.K2 + numRecipPes - 1 ) / numRecipPes;
  bsize = myGrid.block1 * myGrid.dim2 * myGrid.dim3;

  Lattice lattice = patchList[0].p->flags.lattice;

  resultsRemaining = numRecipPes;

  int errcount = 0;

  CProxy_ComputePmeMgr pmeProxy(CpvAccess(BOCclass_group).computePmeMgr);
  for (int j=0; j<numRecipPes; j++) {
    int pe = recipPeOrder[j];  // different order
    int start = pe * bsize;
    int len = bsize;
    if ( start >= qsize ) { start = 0; len = 0; }
    if ( start + len > qsize ) { len = qsize - start; }
    int zdim = myGrid.dim3;
    int fstart = start / zdim;
    int flen = len / zdim;
    int fcount = 0;
    int i;

    int g;
    for ( g=0; g<numGrids; ++g ) {
      char *f = f_arr + fstart + g*fsize;
      int fcount_g = 0;
      for ( i=0; i<flen; ++i ) {
        fcount_g += ( f[i] ? 1 : 0 );
      }
      fcount += fcount_g;
      if ( ! recipPeDest[pe] ) {
        if ( fcount_g ) {
          ++errcount;
          iout << iERROR << CkMyPe() << " sending to " << gridPeMap[pe] << ":";
          int iz = -1;
          for ( i=0; i<flen; ++i ) {
            if ( f[i] ) {
              int jz = (i+fstart)/myGrid.K2;
              if ( iz != jz ) { iout << " " << jz;  iz = jz; }
            }
          }
          iout << "\n" << endi;
        }
      }
    }

    if ( ! recipPeDest[pe] ) continue;

    int zlistlen = 0;
    for ( i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) ++zlistlen;
    }

    PmeGridMsg *msg = new (zlistlen, flen*numGrids, fcount*zlistlen, 0) PmeGridMsg;
    msg->sourceNode = CkMyPe();
    msg->lattice = lattice;
    msg->start = fstart;
    msg->len = flen;
    msg->zlistlen = zlistlen;
    int *zlist = msg->zlist;
    zlistlen = 0;
    for ( i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) zlist[zlistlen++] = i;
    }
    float *qmsg = msg->qgrid;
    for ( g=0; g<numGrids; ++g ) {
      char *f = f_arr + fstart + g*fsize;
      memcpy((void*)(msg->fgrid+g*flen),(void*)f,flen*sizeof(char));
      double **q = q_arr + fstart + g*fsize;
      for ( i=0; i<flen; ++i ) {
        if ( f[i] ) {
          for ( int k=0; k<zlistlen; ++k ) {
            *(qmsg++) = q[i][zlist[k]];
          }
        }
      }
    }

#if CHARM_VERSION > 050402
    pmeProxy[gridPeMap[pe]].recvGrid(msg);
#else
    pmeProxy.recvGrid(msg,gridPeMap[pe]);
#endif
  }

  if ( errcount ) NAMD_bug("Stray PME grid charges detected.");

  for (int i=0; i<fsize; ++i) {
    if ( q_arr[i] ) {
      memset( (void*) (q_arr[i]), -1, myGrid.dim3 * sizeof(double) );
    }
  }

}

void ComputePme::copyResults(PmeGridMsg *msg) {

  evir += msg->evir;
  int zdim = myGrid.dim3;
  int flen = msg->len;
  int fstart = msg->start;
  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  float *qmsg = msg->qgrid;
  int g;
  for ( g=0; g<numGrids; ++g ) {
    char *f = msg->fgrid + g*flen;
    double **q = q_arr + fstart + g*fsize;
    for ( int i=0; i<flen; ++i ) {
      if ( f[i] ) {
        for ( int k=0; k<zlistlen; ++k ) {
          q[i][zlist[k]] = *(qmsg++);
        }
      }
    }
  }
}

void ComputePme::ungridForces() {

    SimParameters *simParams = Node::Object()->simParameters;

    Vector *localResults = new Vector[numLocalAtoms*(numGrids>1?2:1)];
    Vector *gridResults;
    if ( fepOn || lesOn ) {
      for(int i=0; i<numLocalAtoms; ++i) { localResults[i] = 0.; }
      gridResults = localResults + numLocalAtoms;
    } else {
      gridResults = localResults;
    }

    Lattice lattice = patchList[0].p->flags.lattice;
    int g = 0;
    for ( g=0; g<numGrids; ++g ) {
      myRealSpace[g]->compute_forces(q_arr+g*fsize, localGridData[g], gridResults);
      delete myRealSpace[g];
      scale_forces(gridResults, numGridAtoms[g], lattice);

      if ( fepOn || lesOn ) {
        double scale = 1.;
        if ( fepOn ) {
          if ( g == 0 ) scale = simParams->lambda;
          else if ( g == 1 ) scale = 1. - simParams->lambda;
        } else if ( lesOn ) {
          scale = 1.0 / (double)lesFactor;
        }
        int nga = 0;
        for(int i=0; i<numLocalAtoms; ++i) {
          if ( localPartition[i] == 0 || localPartition[i] == (g+1) ) {
            localResults[i] += gridResults[nga++] * scale;
          }
        }
      }
    }

    delete [] localData;
    delete [] localPartition;

    Vector *results_ptr = localResults;
    ResizeArrayIter<PatchElem> ap(patchList);

    // add in forces
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      Results *r = (*ap).forceBox->open();
      Force *f = r->f[Results::slow];
      int numAtoms = (*ap).p->getNumAtoms();

      for(int i=0; i<numAtoms; ++i)
        {
        f[i].x += results_ptr->x;
        f[i].y += results_ptr->y;
        f[i].z += results_ptr->z;
        ++results_ptr;
      }
  
      (*ap).forceBox->close(&r);
    }

    delete [] localResults;
   
    for ( g=0; g<numGrids; ++g ) {
      double scale = 1.;
      if ( fepOn ) {
        if ( g == 0 ) scale = simParams->lambda;
        else if ( g == 1 ) scale = 1. - simParams->lambda;
      } else if ( lesOn ) {
        scale = 1.0 / (double)lesFactor;
      }
      reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += evir[0+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XX) += evir[1+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XY) += evir[2+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += evir[3+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YX) += evir[2+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YY) += evir[4+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += evir[5+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += evir[3+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += evir[5+7*g] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += evir[6+7*g] * scale;

      double scale2 = 0.;
      if ( fepOn && g == 0 ) scale2 = simParams->lambda2;
      else if ( fepOn && g == 1 ) scale2 = 1. - simParams->lambda2;
      reduction->item(REDUCTION_ELECT_ENERGY_SLOW_F) += evir[0+7*g] * scale2;
    }
    reduction->submit();

}

#include "ComputePmeMgr.def.h"

