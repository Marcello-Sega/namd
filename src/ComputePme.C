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

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif


class PmeGridMsg : public CMessage_PmeGridMsg {
public:

  int sourceNode;
  Lattice lattice;
  double energy;
  double virial[6];
  int start;
  int len;
  char *fgrid;
  double *qgrid;
};

class PmeTransMsg : public CMessage_PmeTransMsg {
public:

  int sourceNode;
  int x_start;
  int nx;
  double *qgrid;
};


class PmeUntransMsg : public CMessage_PmeUntransMsg {
public:

  int sourceNode;
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

#ifdef NAMD_FFTW
  fftw_plan forward_plan_x, backward_plan_x;
  rfftwnd_plan forward_plan_yz, backward_plan_yz;
  fftw_complex *work;
#else
  double *work;
#endif

  LocalPmeInfo *localInfo;
  int qgrid_start;
  int qgrid_len;
  int fgrid_start;
  int fgrid_len;

  int numSources;
  int numRecipPes;
  int numDestRecipPes;
  int firstDestRecipPe;
  int myRecipPe;
  int *recipPeMap;
  int *recipPeDest;
  int grid_count;
  int trans_count;
  int untrans_count;
  int ungrid_count;
  PmeTransMsg **trans_buf;
  int trans_buf_len;
  PmeUntransMsg **untrans_buf;
  int untrans_buf_len;
  PmeGridMsg **gridmsg_reuse;
  double recipEnergy;
  double recip_vir[6];
};

ComputePmeMgr::ComputePmeMgr() : pmeProxy(thisgroup), pmeCompute(0) {
  CpvAccess(BOCclass_group).computePmeMgr = thisgroup;
  myKSpace = 0;
  localInfo = new LocalPmeInfo[CkNumPes()];
  recipPeMap = new int[CkNumPes()];
  recipPeDest = new int[CkNumPes()];
  qgrid = 0;
  work = 0;
  grid_count = 0;
  trans_count = 0;
  untrans_count = 0;
  ungrid_count = 0;
  trans_buf = new PmeTransMsg*[CkNumPes()];
  trans_buf_len = 0;
  untrans_buf = new PmeUntransMsg*[CkNumPes()];
  untrans_buf_len = 0;
  gridmsg_reuse= new PmeGridMsg*[CkNumPes()];
}

void ComputePmeMgr::initialize(CkQdMsg *msg) {
  delete msg;

  SimParameters *simParams = Node::Object()->simParameters;

  {  // decide how many pes to use for reciprocal sum
    int nrp = 1;

    // rules based on work available
    int minslices = 1;
    int dimx = simParams->PMEGridSizeX;
    int nrpx = ( dimx + minslices - 1 ) / minslices;
    if ( nrpx > nrp ) nrp = nrpx;
    int dimy = simParams->PMEGridSizeY;
    int nrpy = ( dimy + minslices - 1 ) / minslices;
    if ( nrpy > nrp ) nrp = nrpy;

    // rules based on processors available
    int nrpp = CkNumPes();
    // if ( nrpp > 32 ) nrpp = 32;  // cap to limit messages
    if ( nrpp < nrp ) nrp = nrpp;

    // user override
    int nrps = simParams->PMEProcessors;
    if ( nrps > CkNumPes() ) nrps = CkNumPes();
    if ( nrps > 0 ) nrp = nrps;

    // make sure there aren't any totally empty processors
    int bx = ( dimx + nrp - 1 ) / nrp;
    int nrpbx = ( dimx + bx - 1 ) / bx;
    int by = ( dimy + nrp - 1 ) / nrp;
    int nrpby = ( dimy + by - 1 ) / by;
    nrp = ( nrpby > nrpbx ? nrpby : nrpbx );
    if ( bx != ( dimx + nrp - 1 ) / nrp )
      NAMD_bug("Error in selecting number of PME processors.");
    if ( by != ( dimy + nrp - 1 ) / nrp )
      NAMD_bug("Error in selecting number of PME processors.");

    numRecipPes = nrp;
  }
  if ( ! CkMyPe() ) {
    iout << iINFO << "PME using " << numRecipPes <<
      " processors for FFT and reciprocal sum.\n" << endi;
  }

  myRecipPe = -1;
  for ( int i=0; i<numRecipPes; ++i ) {
    recipPeMap[i] = CkNumPes() - numRecipPes + i;
    if ( recipPeMap[i] == CkMyPe() ) myRecipPe = i;
  }

  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
  myGrid.block1 = ( myGrid.K1 + numRecipPes - 1 ) / numRecipPes;
  myGrid.block2 = ( myGrid.K2 + numRecipPes - 1 ) / numRecipPes;

  int nx = 0;
  int ny = 0;
  for ( int pe = 0; pe < numRecipPes; ++pe ) {
    localInfo[pe].x_start = nx;
    nx += myGrid.block1;
    if ( nx > myGrid.K1 ) nx = myGrid.K1;
    localInfo[pe].nx = nx - localInfo[pe].x_start;
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

  // make sure that we don't get ahead of ourselves on this node
  if ( CkMyPe() < numPatches && myRecipPe >= 0 ) {
    source_flags[CkMyPe()] = 1;
    recipPeDest[myRecipPe] = 1;
  }

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
      if ( myRecipPe >= 0 && ix >= localInfo[myRecipPe].x_start &&
           ix < localInfo[myRecipPe].x_start + localInfo[myRecipPe].nx ) {
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

  firstDestRecipPe = CkMyPe() % numDestRecipPes;
  int c = 0;
  for ( node=0; node<numNodes ; ++node ) {
    if ( recipPeDest[node] ) {
      if ( c == firstDestRecipPe ) {
        firstDestRecipPe = node;
        break;
      }
      ++c;
    }
  }

  // CkPrintf("PME on node %d has %d sources and %d destinations (first=%d)\n",
  //           CkMyPe(), numSources, numDestRecipPes,firstDestRecipPe);

  }  // decide how many pes this node exchanges charges with (end)

  ungrid_count = numDestRecipPes;

  if ( myRecipPe < 0 ) return;
  // the following only for nodes doing reciprocal sum

  int k2_start = localInfo[myRecipPe].y_start_after_transpose;
  int k2_end = k2_start + localInfo[myRecipPe].ny_after_transpose;
  myKSpace = new PmeKSpace(myGrid, k2_start, k2_end);

  int local_size = myGrid.block1 * myGrid.K2 * myGrid.dim3;
  int local_size_2 = myGrid.block2 * myGrid.K1 * myGrid.dim3;
  if ( local_size < local_size_2 ) local_size = local_size_2;
  qgrid = new double[local_size];

  qgrid_start = localInfo[myRecipPe].x_start * myGrid.K2 * myGrid.dim3;
  qgrid_len = localInfo[myRecipPe].nx * myGrid.K2 * myGrid.dim3;
  fgrid_start = localInfo[myRecipPe].x_start * myGrid.K2;
  fgrid_len = localInfo[myRecipPe].nx * myGrid.K2;

  int n[3]; n[0] = myGrid.K1; n[1] = myGrid.K2; n[2] = myGrid.K3;

#ifdef NAMD_FFTW
  work = new fftw_complex[n[0]];

  if ( ! CkMyPe() ) iout << iINFO << "Optimizing 4 FFT steps.  1..." << endi;
  forward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_REAL_TO_COMPLEX,
	FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
  if ( ! CkMyPe() ) iout << " 2..." << endi;
  forward_plan_x = fftw_create_plan_specific(n[0], FFTW_REAL_TO_COMPLEX,
	FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) qgrid,
	localInfo[myRecipPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
  if ( ! CkMyPe() ) iout << " 3..." << endi;
  backward_plan_x = fftw_create_plan_specific(n[0], FFTW_COMPLEX_TO_REAL,
	FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) qgrid,
	localInfo[myRecipPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
  if ( ! CkMyPe() ) iout << " 4..." << endi;
  backward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_COMPLEX_TO_REAL,
	FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
  if ( ! CkMyPe() ) iout << "   Done.\n" << endi;
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  if ( numSources == 0 ) NAMD_bug("PME grid elements exist without sources.");
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_len * sizeof(double) );
}

ComputePmeMgr::~ComputePmeMgr() {
  delete myKSpace;
  delete [] localInfo;
  delete [] recipPeMap;
  delete [] recipPeDest;
  delete [] qgrid;
  delete [] work;
  delete [] trans_buf;
  delete [] untrans_buf;
  delete [] gridmsg_reuse;
}

void ComputePmeMgr::sendGrid(void) {
  pmeCompute->sendData(numRecipPes,firstDestRecipPe,recipPeDest,recipPeMap);
}

void ComputePmeMgr::recvGrid(PmeGridMsg *msg) {
  // CkPrintf("recvGrid from %d on Pe(%d)\n",msg->sourceNode,CkMyPe());
  if ( grid_count == 0 ) {
    NAMD_bug("Message order failure in ComputePmeMgr::recvGrid\n");
  }
  if ( grid_count == numSources ) {
    lattice = msg->lattice;
  }

  char *f = msg->fgrid;
  int zdim = myGrid.dim3;
  double *q = qgrid;
  double *qmsg = msg->qgrid;
  for ( int i=0; i<fgrid_len; ++i ) {
    if ( f[i] ) {
      for ( int j=0; j<zdim; ++j ) {
        *(q++) += *(qmsg++);
      }
    } else {
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
  rfftwnd_real_to_complex(forward_plan_yz, localInfo[myRecipPe].nx,
	qgrid, 1, myGrid.dim2 * myGrid.dim3, 0, 0, 0);
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
  int nx = localInfo[myRecipPe].nx;
  int x_start = localInfo[myRecipPe].x_start;
  int slicelen = myGrid.K2 * zdim;
  for (int j=0; j<numRecipPes; j++) {
    int pe = ( j + myRecipPe ) % numRecipPes;  // different order on each node
    LocalPmeInfo &li = localInfo[pe];
    int cpylen = li.ny_after_transpose * zdim;
    PmeTransMsg *newmsg = new (nx * cpylen,0) PmeTransMsg;
    newmsg->sourceNode = myRecipPe;
    newmsg->x_start = x_start;
    newmsg->nx = nx;
    double *q = qgrid + li.y_start_after_transpose * zdim;
    double *qmsg = newmsg->qgrid;
    for ( int x = 0; x < nx; ++x ) {
      memcpy((void*)qmsg, (void*)q, cpylen*sizeof(double));
      q += slicelen;
      qmsg += cpylen;
    }
#if CHARM_VERSION > 050402
    pmeProxy[recipPeMap[pe]].recvTrans(newmsg);
#else
    pmeProxy.recvTrans(newmsg,recipPeMap[pe]);
#endif
  }

  trans_count = numRecipPes;
  if ( trans_buf_len ) {
    // CkPrintf("resending %d recvTrans on Pe(%d)\n",trans_buf_len,CkMyPe());
    for ( int m=0; m<trans_buf_len; ++m ) {
#if CHARM_VERSION > 050402
      pmeProxy[CkMyPe()].recvTrans(trans_buf[m]);
#else
      pmeProxy.recvTrans(trans_buf[m],CkMyPe());
#endif
    }
    trans_buf_len = 0;
  }
}

void ComputePmeMgr::recvTrans(PmeTransMsg *msg) {
  // CkPrintf("recvTrans on Pe(%d)\n",CkMyPe());
  if ( trans_count == 0 ) {
    trans_buf[trans_buf_len++] = msg;
    return;
  }

  int zdim = myGrid.dim3;
  // int y_start = localInfo[myRecipPe].y_start_after_transpose;
  int ny = localInfo[myRecipPe].ny_after_transpose;
  int x_start = msg->x_start;
  int nx = msg->nx;
  memcpy((void*)(qgrid + x_start*ny*zdim), (void*)(msg->qgrid),
				nx*ny*zdim*sizeof(double));

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
  // int y_start = localInfo[myRecipPe].y_start_after_transpose;
  int ny = localInfo[myRecipPe].ny_after_transpose;

  // finish forward FFT (x dimension)
#ifdef NAMD_FFTW
  fftw(forward_plan_x, ny * zdim / 2, (fftw_complex *) qgrid,
	ny * zdim / 2, 1, work, 1, 0);
#endif

  // reciprocal space portion of PME
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  recipEnergy = myKSpace->compute_energy(qgrid, lattice, ewaldcof, recip_vir);
  // CkPrintf("Ewald reciprocal energy = %f\n", recipEnergy);

  // start backward FFT (x dimension)
#ifdef NAMD_FFTW
  fftw(backward_plan_x, ny * zdim / 2, (fftw_complex *) qgrid,
	ny * zdim / 2, 1, work, 1, 0);
#endif

#if CHARM_VERSION > 050402
  pmeProxy[CkMyPe()].sendUntrans();
#else
  pmeProxy.sendUntrans(CkMyPe());
#endif
}

void ComputePmeMgr::sendUntrans(void) {

  int zdim = myGrid.dim3;
  int y_start = localInfo[myRecipPe].y_start_after_transpose;
  int ny = localInfo[myRecipPe].ny_after_transpose;

  // send data for reverse transpose
  for (int j=0; j<numRecipPes; j++) {
    int pe = ( j + myRecipPe ) % numRecipPes;  // different order on each node
    LocalPmeInfo &li = localInfo[pe];
    int x_start =li.x_start;
    int nx = li.nx;
    PmeUntransMsg *newmsg = new (nx*ny*zdim,0) PmeUntransMsg;
    newmsg->sourceNode = myRecipPe;
    newmsg->y_start = y_start;
    newmsg->ny = ny;
    memcpy((void*)(newmsg->qgrid), (void*)(qgrid + x_start*ny*zdim),
				nx*ny*zdim*sizeof(double));
#if CHARM_VERSION > 050402
    pmeProxy[recipPeMap[pe]].recvUntrans(newmsg);
#else
    pmeProxy.recvUntrans(newmsg,recipPeMap[pe]);
#endif
  }

  untrans_count = numRecipPes;
  if ( untrans_buf_len ) {
    // CkPrintf("resending %d recvUntrans on Pe(%d)\n",untrans_buf_len,CkMyPe());
    for ( int m=0; m<untrans_buf_len; ++m ) {
#if CHARM_VERSION > 050402
      pmeProxy[CkMyPe()].recvUntrans(untrans_buf[m]);
#else
      pmeProxy.recvUntrans(untrans_buf[m],CkMyPe());
#endif
    }
    untrans_buf_len = 0;
  }
}

void ComputePmeMgr::recvUntrans(PmeUntransMsg *msg) {
  // CkPrintf("recvUntrans on Pe(%d)\n",CkMyPe());
  if ( untrans_count == 0 ) {
    untrans_buf[untrans_buf_len++] = msg;
    return;
  }

  int zdim = myGrid.dim3;
  // int x_start = localInfo[myRecipPe].x_start;
  int nx = localInfo[myRecipPe].nx;
  int y_start = msg->y_start;
  int ny = msg->ny;
  int slicelen = myGrid.K2 * zdim;
  int cpylen = ny * zdim;
  double *q = qgrid + y_start * zdim;
  double *qmsg = msg->qgrid;
  for ( int x = 0; x < nx; ++x ) {
    memcpy((void*)q, (void*)qmsg, cpylen*sizeof(double));
    q += slicelen;
    qmsg += cpylen;
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
  rfftwnd_complex_to_real(backward_plan_yz, localInfo[myRecipPe].nx,
	(fftw_complex *) qgrid, 1, myGrid.dim2 * myGrid.dim3 / 2, 0, 0, 0);
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
      newmsg->energy = recipEnergy;
      for ( int i=0; i<6; ++i ) {
        newmsg->virial[i] = recip_vir[i];
      }
    } else {
      newmsg->energy = 0.;
      for ( int i=0; i<6; ++i ) {
        newmsg->virial[i] = 0.;
      }
    }
    int zdim = myGrid.dim3;
    int flen = newmsg->len;
    int fstart = newmsg->start;
    char *f = newmsg->fgrid;
    double *qmsg = newmsg->qgrid;
    double *q = qgrid + (fstart-fgrid_start) * zdim;
    for ( int i=0; i<flen; ++i ) {
      if ( f[i] ) {
        memcpy((void*)(qmsg),(void*)(q),zdim*sizeof(double));
        qmsg += zdim;
      }
      q += zdim;
    }
    newmsg->sourceNode = myRecipPe;

#if CHARM_VERSION > 050402
    pmeProxy[pe].recvUngrid(newmsg);
#else
    pmeProxy.recvUngrid(newmsg,pe);
#endif
  }
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_len * sizeof(double) );
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
  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
  qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  fsize = myGrid.K1 * myGrid.dim2;
  q_arr = new double*[fsize];
  memset( (void*) q_arr, 0, fsize * sizeof(double*) );
  f_arr = new char[fsize];
}

ComputePme::~ComputePme()
{
  for (int i=0; i<fsize; ++i) {
    if ( q_arr[i] ) {
      delete [] q_arr[i];
    }
  }
  delete [] q_arr;
  delete [] f_arr;
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

  localData = new PmeParticle[numLocalAtoms];

  // get positions and charges
  PmeParticle * data_ptr = localData;
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
      Vector tmp = lattice.delta(x[i].position);
      data_ptr->x = tmp.x;
      data_ptr->y = tmp.y;
      data_ptr->z = tmp.z;
      data_ptr->cg = coloumb_sqrt * x[i].charge;
      ++data_ptr;
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
  }

  int i;
  for ( i=0; i<6; ++i ) { virial[i] = 0; }
  energy = 0;

  // calculate self energy
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  BigReal selfEnergy = 0;
  data_ptr = localData;
  for(i=0; i<numLocalAtoms; ++i)
  {
    selfEnergy += data_ptr->cg * data_ptr->cg;
    ++data_ptr;
  }
  selfEnergy *= -1. * ewaldcof / SQRT_PI;
  energy += selfEnergy;

  for (i=0; i<fsize; ++i) {
    if ( q_arr[i] ) {
      memset( (void*) (q_arr[i]), 0, myGrid.dim3 * sizeof(double) );
    }
  }
  memset( (void*) f_arr, 0, fsize * sizeof(char) );
  myRealSpace = new PmeRealSpace(myGrid,numLocalAtoms);
  scale_coordinates(localData, numLocalAtoms, lattice, myGrid);
  myRealSpace->fill_charges(q_arr, f_arr, localData);

  CProxy_ComputePmeMgr pmeProxy(CpvAccess(BOCclass_group).computePmeMgr);
#if CHARM_VERSION > 050402
  pmeProxy[CkMyPe()].sendGrid();
#else
  pmeProxy.sendGrid(CkMyPe());
#endif
}

void ComputePme::sendData(int numRecipPes, int firstDestRecipPe,
				int *recipPeDest, int *recipPeMap) {

  // iout << "Sending charge grid for " << numLocalAtoms << " atoms to FFT on " << iPE << ".\n" << endi;

  myGrid.block1 = ( myGrid.K1 + numRecipPes - 1 ) / numRecipPes;
  myGrid.block2 = ( myGrid.K2 + numRecipPes - 1 ) / numRecipPes;
  bsize = myGrid.block1 * myGrid.dim2 * myGrid.dim3;

  Lattice lattice = patchList[0].p->flags.lattice;

  resultsRemaining = numRecipPes;

  CProxy_ComputePmeMgr pmeProxy(CpvAccess(BOCclass_group).computePmeMgr);
  for (int j=0; j<numRecipPes; j++) {
    int pe = ( j + firstDestRecipPe ) % numRecipPes;  // different order
    int start = pe * bsize;
    int len = bsize;
    if ( start >= qsize ) { start = 0; len = 0; }
    if ( start + len > qsize ) { len = qsize - start; }
    int zdim = myGrid.dim3;
    int fstart = start / zdim;
    int flen = len / zdim;
    char *f = f_arr + fstart;
    int fcount = 0;
    int i;
    for ( i=0; i<flen; ++i ) {
      fcount += ( f[i] ? 1 : 0 );
    }
    // CkPrintf("count(%d -> %d) = %d\n",CkMyPe(),pe,fcount);

    if ( ! recipPeDest[pe] ) {
      if ( fcount ) NAMD_bug("Stray PME grid charges detected.");
      continue;
    }

    PmeGridMsg *msg = new (flen, fcount*zdim, 0) PmeGridMsg;
    msg->sourceNode = CkMyPe();
    msg->lattice = lattice;
    msg->start = fstart;
    msg->len = flen;
    memcpy((void*)(msg->fgrid),(void*)(f),flen*sizeof(char));

    double **q = q_arr + fstart;
    double *qmsg = msg->qgrid;
    for ( i=0; i<flen; ++i ) {
      if ( f[i] ) {
        memcpy((void*)(qmsg),(void*)(q[i]),zdim*sizeof(double));
        qmsg += zdim;
      }
    }

#if CHARM_VERSION > 050402
    pmeProxy[recipPeMap[pe]].recvGrid(msg);
#else
    pmeProxy.recvGrid(msg,recipPeMap[pe]);
#endif
  }

}

void ComputePme::copyResults(PmeGridMsg *msg) {

  // if ( CkMyPe() == 0 ) {
    energy += msg->energy;
    for ( int i=0; i<6; ++i ) {
      virial[i] += msg->virial[i];
    }
  // }
  int zdim = myGrid.dim3;
  int flen = msg->len;
  int fstart = msg->start;
  char *f = msg->fgrid;
  double *qmsg = msg->qgrid;
  double **q = q_arr + fstart;
  for ( int i=0; i<flen; ++i ) {
    if ( f[i] ) {
      memcpy((void*)(q[i]),(void*)(qmsg),zdim*sizeof(double));
      qmsg += zdim;
    }
  }
}

void ComputePme::ungridForces() {

    Vector *localResults = new Vector[numLocalAtoms];
    myRealSpace->compute_forces(q_arr, localData, localResults);
    delete [] localData;
    delete myRealSpace;
    Lattice lattice = patchList[0].p->flags.lattice;
    scale_forces(localResults, numLocalAtoms, lattice);

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
   
    reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += energy;
    reduction->item(REDUCTION_VIRIAL_SLOW_XX) += (BigReal)(virial[0]);
    reduction->item(REDUCTION_VIRIAL_SLOW_XY) += (BigReal)(virial[1]);
    reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += (BigReal)(virial[2]);
    reduction->item(REDUCTION_VIRIAL_SLOW_YX) += (BigReal)(virial[1]);
    reduction->item(REDUCTION_VIRIAL_SLOW_YY) += (BigReal)(virial[3]);
    reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += (BigReal)(virial[4]);
    reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += (BigReal)(virial[2]);
    reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += (BigReal)(virial[4]);
    reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += (BigReal)(virial[5]);
    reduction->submit();

}

#include "ComputePmeMgr.def.h"

