/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputePme.h"
#include "ComputePmeMsgs.h"
#include <fftw.h>
#include <rfftw.h>
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

#include "PmeCoulomb.h"

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif

class PmeNullMsg : public CMessage_PmeNullMsg {
};

class PmeGridMsg : public CMessage_PmeGridMsg {
public:

  int sourceNode;
  Lattice lattice;
  int start;
  int len;
  double *qgrid;

  VARSIZE_DECL(PmeGridMsg);
};

VARSIZE_MSG(PmeGridMsg,
  VARSIZE_ARRAY(qgrid);
)

class PmeTransMsg : public CMessage_PmeTransMsg {
public:

  int sourceNode;
  int x_start;
  int nx;
  double *qgrid;

  VARSIZE_DECL(PmeTransMsg);
};

VARSIZE_MSG(PmeTransMsg,
  VARSIZE_ARRAY(qgrid);
)

class PmeUntransMsg : public CMessage_PmeUntransMsg {
public:

  int sourceNode;
  double energy;
  double virial[6];
  int y_start;
  int ny;
  double *qgrid;

  VARSIZE_DECL(PmeUntransMsg);
};

VARSIZE_MSG(PmeUntransMsg,
  VARSIZE_ARRAY(qgrid);
)

class PmeUngridMsg : public CMessage_PmeUngridMsg {
public:

  int sourceNode;
  Lattice lattice;
  int start;
  int len;
  double *qgrid;

  VARSIZE_DECL(PmeUngridMsg);
};

VARSIZE_MSG(PmeUngridMsg,
  VARSIZE_ARRAY(qgrid);
)

struct LocalPmeInfo {
  int nx, x_start;
  int ny_after_transpose, y_start_after_transpose;
  int total_size;
};

class ComputePmeMgr : public BOCclass {
public:
  ComputePmeMgr();
  ~ComputePmeMgr();

  void sendGrid(PmeNullMsg *);
  void recvGrid(PmeGridMsg *);
  void gridCalc1(PmeNullMsg *);
  void sendTrans(PmeNullMsg *);
  void recvTrans(PmeTransMsg *);
  void gridCalc2(PmeNullMsg *);
  void sendUntrans(PmeNullMsg *);
  void recvUntrans(PmeUntransMsg *);
  void gridCalc3(PmeNullMsg *);
  void sendUngrid(PmeNullMsg *);
  void recvUngrid(PmeUngridMsg *);
  void ungridCalc(PmeNullMsg *);

  void setCompute(ComputePme *c) { pmeCompute = c; }

private:
  int initialized;
  void initialize();

  CProxy_ComputePmeMgr pmeProxy;
  ComputePme *pmeCompute;
  PmeGrid myGrid;
  Lattice lattice;
  PmeKSpace *myKSpace;
  double *qgrid;

  fftw_plan forward_plan_x, backward_plan_x;
  rfftwnd_plan forward_plan_yz, backward_plan_yz;
  fftw_complex *work;

  LocalPmeInfo *localInfo;
  int qgrid_start;
  int qgrid_len;

  int numSources;
  int grid_count;
  int trans_count;
  int untrans_count;
  int ungrid_count;
  PmeTransMsg **trans_buf;
  int trans_buf_len;
  PmeUntransMsg **untrans_buf;
  int untrans_buf_len;
  double recipEnergy;
  double recip_vir[6];
};

ComputePmeMgr::ComputePmeMgr() : pmeProxy(thisgroup), pmeCompute(0),
						initialized(0) {
  CpvAccess(BOCclass_group).computePmeMgr = thisgroup;
  myKSpace = 0;
  localInfo = 0;
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
}

void ComputePmeMgr::initialize() {
  if ( initialized ) { return; }
  initialized = 1;

  SimParameters *simParams = Node::Object()->simParameters;
  WorkDistrib *workDistrib = Node::Object()->workDistrib;

  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
  myGrid.block1 = ( myGrid.K1 + CkNumPes() - 1 ) / CkNumPes();
  myGrid.block2 = ( myGrid.K2 + CkNumPes() - 1 ) / CkNumPes();

  int n[3]; n[0] = myGrid.K1; n[1] = myGrid.K2; n[2] = myGrid.K3;
  forward_plan_x = fftw_create_plan(n[0], FFTW_REAL_TO_COMPLEX,
					FFTW_MEASURE | FFTW_IN_PLACE);
  forward_plan_yz = rfftwnd_create_plan(2, n+1, FFTW_REAL_TO_COMPLEX,
					FFTW_MEASURE | FFTW_IN_PLACE);
  backward_plan_x = fftw_create_plan(n[0], FFTW_COMPLEX_TO_REAL,
					FFTW_MEASURE | FFTW_IN_PLACE);
  backward_plan_yz = rfftwnd_create_plan(2, n+1, FFTW_COMPLEX_TO_REAL,
					FFTW_MEASURE | FFTW_IN_PLACE);
  work = new fftw_complex[n[0]];

  localInfo = new LocalPmeInfo[CkNumPes()];

  int nx = 0;
  int ny = 0;
  for ( int pe = 0; pe < CkNumPes(); ++pe ) {
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

  int k2_start = localInfo[CkMyPe()].y_start_after_transpose;
  int k2_end = k2_start + localInfo[CkMyPe()].ny_after_transpose;
  myKSpace = new PmeKSpace(myGrid, k2_start, k2_end);

  int local_size = myGrid.block1 * myGrid.K2 * myGrid.dim3;
  int local_size_2 = myGrid.block2 * myGrid.K1 * myGrid.dim3;
  if ( local_size < local_size_2 ) local_size = local_size_2;
  qgrid = new double[local_size];

  qgrid_start = localInfo[CkMyPe()].x_start * myGrid.K2 * myGrid.dim3;
  qgrid_len = localInfo[CkMyPe()].nx * myGrid.K2 * myGrid.dim3;

  numSources = CkNumPes();
  int npatches = (PatchMap::Object())->numPatches();
  if ( numSources > npatches ) numSources = npatches;

  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_len * sizeof(double) );
  ungrid_count = CkNumPes();
}

ComputePmeMgr::~ComputePmeMgr() {
  delete myKSpace;
  delete [] localInfo;
  delete [] qgrid;
  delete [] work;
  delete [] trans_buf;
  delete [] untrans_buf;
}

void ComputePmeMgr::sendGrid(PmeNullMsg *msg) {
  delete msg;
  pmeCompute->sendData();
}

void ComputePmeMgr::recvGrid(PmeGridMsg *msg) {
  initialize();
  // CkPrintf("recvGrid on Pe(%d)\n",CkMyPe());
  if ( grid_count == 0 ) {
    NAMD_bug("Message order failure in ComputePmeMgr::recvGrid\n");
  }
  if ( grid_count == numSources ) {
    lattice = msg->lattice;
  }

  double *q = qgrid + (msg->start - qgrid_start);
  double *qm = msg->qgrid;
  for ( int i=0; i<msg->len; ++i ) {
    *(q++) += *(qm++);
  }
  delete msg;
  --grid_count;

  if ( grid_count == 0 ) {
    pmeProxy.gridCalc1(new PmeNullMsg,CkMyPe());
  }
}

void ComputePmeMgr::gridCalc1(PmeNullMsg *msg) {
  // CkPrintf("gridCalc1 on Pe(%d)\n",CkMyPe());
  delete msg;

  rfftwnd_real_to_complex(forward_plan_yz, localInfo[CkMyPe()].nx,
	qgrid, 1, myGrid.dim2 * myGrid.dim3, 0, 0, 0);

  pmeProxy.sendTrans(new PmeNullMsg, CkMyPe());
}

void ComputePmeMgr::sendTrans(PmeNullMsg *msg) {
  delete msg;

  // send data for transpose
  int zdim = myGrid.dim3;
  int nx = localInfo[CkMyPe()].nx;
  int x_start = localInfo[CkMyPe()].x_start;
  int slicelen = myGrid.K2 * zdim;
  for ( int pe=0; pe<CkNumPes(); ++pe ) {
    LocalPmeInfo &li = localInfo[pe];
    int cpylen = li.ny_after_transpose * zdim;
    int msglen = nx * cpylen;
    PmeTransMsg *newmsg = new (&msglen,0) PmeTransMsg;
    newmsg->sourceNode = CkMyPe();
    newmsg->x_start = x_start;
    newmsg->nx = nx;
    double *q = qgrid + li.y_start_after_transpose * zdim;
    double *qmsg = newmsg->qgrid;
    for ( int x = 0; x < nx; ++x ) {
      memcpy((void*)qmsg, (void*)q, cpylen*sizeof(double));
      q += slicelen;
      qmsg += cpylen;
    }
    pmeProxy.recvTrans(newmsg,pe);
  }

  trans_count = CkNumPes();
  if ( trans_buf_len ) {
    // CkPrintf("resending %d recvTrans on Pe(%d)\n",trans_buf_len,CkMyPe());
    for ( int m=0; m<trans_buf_len; ++m ) {
      pmeProxy.recvTrans(trans_buf[m],CkMyPe());
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
  int y_start = localInfo[CkMyPe()].y_start_after_transpose;
  int ny = localInfo[CkMyPe()].ny_after_transpose;
  int x_start = msg->x_start;
  int nx = msg->nx;
  memcpy((void*)(qgrid + x_start*ny*zdim), (void*)(msg->qgrid),
				nx*ny*zdim*sizeof(double));

  delete msg;
  --trans_count;

  if ( trans_count == 0 ) {
    pmeProxy.gridCalc2(new PmeNullMsg,CkMyPe());
  }
}

void ComputePmeMgr::gridCalc2(PmeNullMsg *msg) {
  // CkPrintf("gridCalc2 on Pe(%d)\n",CkMyPe());
  delete msg;

  int zdim = myGrid.dim3;
  int y_start = localInfo[CkMyPe()].y_start_after_transpose;
  int ny = localInfo[CkMyPe()].ny_after_transpose;

  // finish forward FFT (x dimension)
  fftw(forward_plan_x, ny * zdim / 2, (fftw_complex *) qgrid,
	ny * zdim / 2, 1, work, 1, 0);

  // reciprocal space portion of PME
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  recipEnergy = myKSpace->compute_energy(qgrid, lattice, ewaldcof, recip_vir);
  // CkPrintf("Ewald reciprocal energy = %f\n", recipEnergy);

  // start backward FFT (x dimension)
  fftw(backward_plan_x, ny * zdim / 2, (fftw_complex *) qgrid,
	ny * zdim / 2, 1, work, 1, 0);

  pmeProxy.sendUntrans(new PmeNullMsg, CkMyPe());
}

void ComputePmeMgr::sendUntrans(PmeNullMsg *msg) {
  delete msg;

  int zdim = myGrid.dim3;
  int y_start = localInfo[CkMyPe()].y_start_after_transpose;
  int ny = localInfo[CkMyPe()].ny_after_transpose;

  // send data for reverse transpose
  for ( int pe=0; pe<CkNumPes(); ++pe ) {
    LocalPmeInfo &li = localInfo[pe];
    int x_start =li.x_start;
    int nx = li.nx;
    int msglen = nx*ny*zdim;
    PmeUntransMsg *newmsg = new (&msglen,0) PmeUntransMsg;
    if ( pe == 0 ) {  // only need these once
      newmsg->energy = recipEnergy;
      for ( int i=0; i<6; ++i ) { newmsg->virial[i] = recip_vir[i]; }
    }
    newmsg->sourceNode = CkMyPe();
    newmsg->y_start = y_start;
    newmsg->ny = ny;
    memcpy((void*)(newmsg->qgrid), (void*)(qgrid + x_start*ny*zdim),
				nx*ny*zdim*sizeof(double));
    pmeProxy.recvUntrans(newmsg,pe);
  }

  untrans_count = CkNumPes();
  if ( untrans_buf_len ) {
    // CkPrintf("resending %d recvUntrans on Pe(%d)\n",untrans_buf_len,CkMyPe());
    for ( int m=0; m<untrans_buf_len; ++m ) {
      pmeProxy.recvUntrans(untrans_buf[m],CkMyPe());
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

  if ( CkMyPe() == 0 ) {
    pmeCompute->copyEnergy(msg);
  }

  int zdim = myGrid.dim3;
  int x_start = localInfo[CkMyPe()].x_start;
  int nx = localInfo[CkMyPe()].nx;
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
    pmeProxy.gridCalc3(new PmeNullMsg,CkMyPe());
  }
}

void ComputePmeMgr::gridCalc3(PmeNullMsg *msg) {
  // CkPrintf("gridCalc3 on Pe(%d)\n",CkMyPe());
  delete msg;

  // finish backward FFT
  rfftwnd_complex_to_real(backward_plan_yz, localInfo[CkMyPe()].nx,
	(fftw_complex *) qgrid, 1, myGrid.dim2 * myGrid.dim3 / 2, 0, 0, 0);

  pmeProxy.sendUngrid(new PmeNullMsg, CkMyPe());
}

void ComputePmeMgr::sendUngrid(PmeNullMsg *msg) {
  delete msg;

  for ( int pe=0; pe<numSources; ++pe ) {
    int msglen = qgrid_len;
    PmeUngridMsg *newmsg = new (&msglen,0) PmeUngridMsg;
    newmsg->sourceNode = CkMyPe();
    newmsg->start = qgrid_start;
    newmsg->len = qgrid_len;
    memcpy((void*)(newmsg->qgrid), (void*)qgrid, qgrid_len*sizeof(double));
    pmeProxy.recvUngrid(newmsg,pe);
  }
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_len * sizeof(double) );
}

void ComputePmeMgr::recvUngrid(PmeUngridMsg *msg) {
  // CkPrintf("recvUngrid on Pe(%d)\n",CkMyPe());
  if ( ungrid_count == 0 ) {
    NAMD_bug("Message order failure in ComputePmeMgr::recvUngrid\n");
  }

  pmeCompute->copyResults(msg);
  delete msg;
  --ungrid_count;

  if ( ungrid_count == 0 ) {
    pmeProxy.ungridCalc(new PmeNullMsg,CkMyPe());
  }
}

void ComputePmeMgr::ungridCalc(PmeNullMsg *msg) {
  // CkPrintf("ungridCalc on Pe(%d)\n",CkMyPe());
  delete msg;

  pmeCompute->ungridForces();

  ungrid_count = CkNumPes();
}

#if 0
void ComputePmeMgr::start() {
#ifndef NAMD_FFTW
  NAMD_die("FFTW (http://www.fftw.org/) is required to use PME.");
#else
  // iout << "ComputePmeMgr thread " << thisIndex << " running on " << iPE << "\n" << endi;

  SimParameters *simParams = Node::Object()->simParameters;
  WorkDistrib *workDistrib = Node::Object()->workDistrib;
  PmeGrid myGrid;
  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
  myGrid.block1 = ( myGrid.K1 + CkNumPes() - 1 ) / CkNumPes();
  myGrid.block2 = ( myGrid.K2 + CkNumPes() - 1 ) / CkNumPes();

  int pe = thisIndex;
  int k1_start = pe * myGrid.block1;
  int k1_end = k1_start + myGrid.block1;
  if ( k1_start > myGrid.K1 ) { k1_start = myGrid.K1; }
  if ( k1_end > myGrid.K1 ) { k1_end = myGrid.K1; }
  PmeKSpace myKSpace(myGrid, k1_start, k1_end);

  rfftwnd_mpi_plan forward_plan;
  forward_plan = rfftw3d_mpi_create_plan(AMPI_COMM_WORLD,
	myGrid.K1,myGrid.K2,myGrid.K3,
	FFTW_REAL_TO_COMPLEX, FFTW_MEASURE | FFTW_IN_PLACE);

  rfftwnd_mpi_plan backward_plan;
  backward_plan = rfftw3d_mpi_create_plan(AMPI_COMM_WORLD,
	myGrid.K1,myGrid.K2,myGrid.K3,
	FFTW_COMPLEX_TO_REAL, FFTW_MEASURE | FFTW_IN_PLACE);

  int local_nx, local_x_start, local_ny_after_transpose,
	local_y_start_after_transpose, total_local_size;
  rfftwnd_mpi_local_sizes(forward_plan,
	&local_nx,&local_x_start,&local_ny_after_transpose,
	&local_y_start_after_transpose,&total_local_size);

  double *q_arr = new double[total_local_size];
  double *q_work = new double[total_local_size];

  int qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  int bsize = myGrid.block1 * myGrid.dim2 * myGrid.dim3;
  int lqsize = bsize;
  if ( (pe+1) * bsize > qsize ) lqsize = qsize - pe * bsize;
  if ( pe * bsize >= qsize ) lqsize = 0;
  int lqstart = ( pe * bsize );
  if ( lqstart > qsize ) lqstart = qsize;
  double *q_buf = new double[lqsize];

  int numSources = CkNumPes();
  int npatches = (PatchMap::Object())->numPatches();
  if ( numSources > npatches ) numSources = npatches;

  int i,j;

while ( 1 ) {

  // get data from ComputePme on all nodes
  int first = 1;
  AMPI_Status status;
  for ( j=0; j<numSources; ++j ) {
    //ckTempoRecv(PME_TAG_QARR,(void*)q_buf,lqsize*sizeof(double));
    AMPI_Recv((void*)q_buf, lqsize*sizeof(double), AMPI_BYTE,
	AMPI_ANY_SOURCE, PME_TAG_QARR, AMPI_COMM_WORLD, &status);
    if ( first ) {
      // this is needed so that q_arr isn't zeroed before being packed
      // can be eliminated when compressed format is used in future
      for ( i=0; i<total_local_size; ++i ) q_arr[i] = 0.;
      first = 0;
    }
    // iout << "ComputePmeMgr thread " << thisIndex << " received " << j << " on " << iPE << "\n" << endi;
    for ( i=0; i<lqsize; ++i ) q_arr[i] += q_buf[i];
  }
  // iout << "ComputePmeMgr thread " << thisIndex << " has all data on " << iPE << "\n" << endi;

  // forward transform
  rfftwnd_mpi(forward_plan, 1, q_arr, q_work, FFTW_NORMAL_ORDER);

  // fill in guts of PME here
  Lattice lattice;
  // ckTempoRecv(PME_TAG_LATTICE,(void*)&lattice,sizeof(Lattice));
  AMPI_Recv((void*)&lattice, sizeof(Lattice), AMPI_BYTE,
	AMPI_ANY_SOURCE, PME_TAG_LATTICE, AMPI_COMM_WORLD, &status);
  // iout << "ComputePmeMgr thread " << thisIndex << " received lattice.\n" << endi;
  double recipEnergy;
  double recip_vir[6];
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  recipEnergy = myKSpace.compute_energy(q_arr, lattice, ewaldcof, recip_vir);

  // backward transform
  rfftwnd_mpi(backward_plan, 1, q_arr, q_work, FFTW_NORMAL_ORDER);

  // send results back to nodes
  for ( j=0; j<numSources; ++j ) {
    ComputePmeResultsMsg *msg = new ComputePmeResultsMsg;
    if ( ! j ) {  // only need these once
      msg->energy = recipEnergy;
      for ( i=0; i<6; ++i ) { msg->virial[i] = recip_vir[i]; }
    }
    msg->start = lqstart;
    msg->q_len = lqsize;
    msg->will_delete_array = 0;
    msg->q_arr = q_arr;
    CProxy_ComputeMgr cm(CpvAccess(BOCclass_group).computeMgr);
    cm.recvComputePmeResults(msg, j);
  }

} // while ( 1 )

  delete [] q_buf;
  delete [] q_work;
  delete [] q_arr;
  // iout << "ComputePmeMgr thread " << thisIndex << " finished on " << iPE << "\n" << endi;
#endif  // NAMD_FFTW
}
#endif


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
 

class ComputePmeMaster {
private:
  friend class ComputePme;
  ComputePme *host;
  ComputePmeMaster(ComputePme *);
  ~ComputePmeMaster();
  void recvData(ComputePmeDataMsg *);
  ResizeArray<int> homeNode;
  ResizeArray<int> endForNode;
  int numWorkingPes;
  int numLocalAtoms;
  PmeParticle *localData;
  SubmitReduction *reduction;
  int runcount;
  PmeCoulomb *myPme;
};

ComputePme::ComputePme(ComputeID c, ComputeMgr *m) :
  ComputeHomePatches(c), comm(m)
{
  DebugM(4,"ComputePme created.\n");

  CProxy_ComputePmeMgr::ckLocalBranch(
	CpvAccess(BOCclass_group).computePmeMgr)->setCompute(this);

  useAvgPositions = 1;

  int numWorkingPes = CkNumPes();
  {
    int npatches=(PatchMap::Object())->numPatches();
    if ( numWorkingPes > npatches ) numWorkingPes = npatches;
  }

/*
  masterNode = numWorkingPes - 1;
  if ( CkMyPe() == masterNode ) {
    master = new ComputePmeMaster(this);
    master->numWorkingPes = numWorkingPes;
  }
  else master = 0;
*/
  master = 0;

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  SimParameters *simParams = Node::Object()->simParameters;
  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
  myGrid.block1 = ( myGrid.K1 + CkNumPes() - 1 ) / CkNumPes();
  myGrid.block2 = ( myGrid.K2 + CkNumPes() - 1 ) / CkNumPes();
  qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  bsize = myGrid.block1 * myGrid.dim2 * myGrid.dim3;
  q_arr = new double[qsize];
}

ComputePme::~ComputePme()
{
  delete [] q_arr;

  delete master;
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

  for (i=0; i<qsize; q_arr[i++] = 0.0);
  myRealSpace = new PmeRealSpace(myGrid,numLocalAtoms);
  scale_coordinates(localData, numLocalAtoms, lattice, myGrid);
  myRealSpace->fill_charges(q_arr, localData);

  CProxy_ComputePmeMgr pmeProxy(CpvAccess(BOCclass_group).computePmeMgr);
  pmeProxy.sendGrid(new PmeNullMsg, CkMyPe());
}

void ComputePme::sendData() {

  // iout << "Sending charge grid for " << numLocalAtoms << " atoms to FFT on " << iPE << ".\n" << endi;

  Lattice lattice = patchList[0].p->flags.lattice;

  resultsRemaining = CkNumPes();

  CProxy_ComputePmeMgr pmeProxy(CpvAccess(BOCclass_group).computePmeMgr);
  for (int pe=0; pe<CkNumPes(); pe++) {
    int start = pe * bsize;
    int len = bsize;
    if ( start >= qsize ) { start = 0; len = 0; }
    if ( start + len > qsize ) { len = qsize - start; }
    int msglens[1];
    msglens[0] = len;
    PmeGridMsg *msg = new (msglens,0) PmeGridMsg;
    msg->sourceNode = CkMyPe();
    msg->lattice = lattice;
    msg->start = start;
    msg->len = len;
    memcpy((void*)(msg->qgrid),(void*)(q_arr+start),len*sizeof(double));
    pmeProxy.recvGrid(msg,pe);
  }

}

void ComputePme::recvData(ComputePmeDataMsg *msg)
{
  if ( master ) {
    master->recvData(msg);
  }
  else NAMD_die("ComputePme::master is NULL!");
}

ComputePmeMaster::ComputePmeMaster(ComputePme *h) :
  host(h), numLocalAtoms(0), runcount(0)
{
  DebugM(4,"ComputePmeMaster created.\n");

//  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  Molecule * molecule = Node::Object()->molecule;
  localData = new PmeParticle[molecule->numAtoms];
  SimParameters * simParams = Node::Object()->simParameters;
  PmeGrid grid;
  grid.K1 = simParams->PMEGridSizeX;
  grid.K2 = simParams->PMEGridSizeY;
  grid.K3 = simParams->PMEGridSizeZ;
  grid.order = simParams->PMEInterpOrder;
  myPme = new PmeCoulomb(grid, molecule->numAtoms);
}

ComputePmeMaster::~ComputePmeMaster()
{
  delete reduction;
  delete [] localData;
  delete myPme;
}

void ComputePmeMaster::recvData(ComputePmeDataMsg *msg)
{ 
  DebugM(4,"ComputePmeMaster::recvData() " << msg->numParticles
	<< " particles from node " << msg->node << "\n");

  {
    homeNode.add(msg->node);
    PmeParticle *data_ptr = localData + numLocalAtoms;
    for ( int j = 0; j < msg->numParticles; ++j, ++data_ptr ) {
      *data_ptr = msg->particles[j];
    }
    numLocalAtoms += msg->numParticles;
    endForNode.add(numLocalAtoms);
    delete msg;
  }

  if ( homeNode.size() < numWorkingPes ) return;  // messages outstanding

  DebugM(4,"ComputePmeMaster::recvData() running serial code.\n");

  // single processor version

  Lattice lattice = host->getFlags()->lattice;
  SimParameters * simParams = Node::Object()->simParameters;
  int i;
  Vector *localResults;
  double recip_vir[6];
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  localResults = new Vector[numLocalAtoms];

  // perform calculations
  BigReal electEnergy = 0;

  // calculate self energy
  PmeParticle *data_ptr = localData;
  for(i=0; i<numLocalAtoms; ++i)
  {
    electEnergy += data_ptr->cg * data_ptr->cg;
    ++data_ptr;
  }
  electEnergy *= -1. * ewaldcof / SQRT_PI;

  DebugM(4,"Ewald self energy: " << electEnergy << "\n");

  double recipEnergy;
  DebugM(4,"Calling compute_recip.\n");
  double pme_start_time = 0;
  if ( runcount == 1 ) pme_start_time = CmiTimer();
  // Last argument should be nonzero if this is the last call
  // Compute it using mytime, tsteps, etc.
  // I'll just let the runtime environment take care of it for now.
  recipEnergy = myPme->compute_recip(localData, lattice, ewaldcof, recip_vir, localResults);
  if ( runcount == 1 ) {
    iout << iINFO << "PME reciprocal sum CPU time per evaluation: "
         << (CmiTimer() - pme_start_time) << "\n" << endi;
  }
  electEnergy += recipEnergy;
  DebugM(4,"Returned from PmeRecipCoulomb->calc_recip.\n");
  // No need to reverse sign from new code

  // send out reductions
  DebugM(4,"Timestep : " << host->getFlags()->step << "\n");
  DebugM(4,"Reciprocal sum energy: " << electEnergy << "\n");
  DebugM(4,"Reciprocal sum virial: " << recip_vir[0] << " " <<
	recip_vir[1] << " " << recip_vir[2] << " " << recip_vir[3] << " " <<
	recip_vir[4] << " " << recip_vir[5] << "\n");
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += electEnergy;
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += (BigReal)(recip_vir[0]);
  reduction->item(REDUCTION_VIRIAL_SLOW_XY) += (BigReal)(recip_vir[1]);
  reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += (BigReal)(recip_vir[2]);
  reduction->item(REDUCTION_VIRIAL_SLOW_YX) += (BigReal)(recip_vir[1]);
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += (BigReal)(recip_vir[3]);
  reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += (BigReal)(recip_vir[4]);
  reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += (BigReal)(recip_vir[2]);
  reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += (BigReal)(recip_vir[4]);
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += (BigReal)(recip_vir[5]);
  reduction->submit();

  Vector *results_ptr;
  results_ptr = localResults;

//  numLocalAtoms = 0;
//  for ( i = 0; i < homeNode.size(); ++i ) {
//    ComputePmeResultsMsg *msg = new ComputePmeResultsMsg;
//    msg->node = homeNode[i];
//    msg->numParticles = endForNode[i] - numLocalAtoms;
//    msg->forces = new Vector[msg->numParticles];
//    for ( int j = 0; j < msg->numParticles; ++j, ++results_ptr ) {
//      msg->forces[j] = *results_ptr;
//    }
//    numLocalAtoms = endForNode[i];
//    host->comm->sendComputePmeResults(msg,homeNode[i]);
//  }
  delete [] localResults;

  // reset
  runcount += 1;
  numLocalAtoms = 0;
  homeNode.resize(0);
  endForNode.resize(0);

}

void ComputePme::copyResults(PmeUngridMsg *msg) {
  memcpy((void*)(q_arr + msg->start), (void*)(msg->qgrid),
				msg->len * sizeof(double));
}

void ComputePme::copyEnergy(PmeUntransMsg *msg) {
  energy += msg->energy;
  for ( int i=0; i<6; ++i ) { virial[i] += msg->virial[i]; }
}

void ComputePme::recvResults(ComputePmeResultsMsg *msg)
{
  --resultsRemaining;

  // iout << "ComputePme received results " << msg->start << " to "
// 	<< (msg->start + msg->q_len) << ".\n" << endi;

  int start = msg->start;
  int end = start + msg->q_len;
  double *data = msg->q_arr;
  int i;
  for ( i=start; i<end; ++i ) { q_arr[i] = *(data++); }
  for ( i=0; i<6; ++i ) { virial[i] += msg->virial[i]; }
  energy += msg->energy;
  delete msg;

  if ( ! resultsRemaining ) ungridForces();

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

