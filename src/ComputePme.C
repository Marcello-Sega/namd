/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_NO_TYPE_PREFIX
#include <fftw.h>
#include <rfftw.h>
#else
#include <sfftw.h>
#include <srfftw.h>
#endif
#endif

#include "InfoStream.h"
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
#include "Priorities.h"

// commlib has been observed to cause hangs when starting load balancing
// #define USE_COMM_LIB 1
#ifdef USE_COMM_LIB
#include "EachToManyMulticastStrategy.h"
#endif

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif

char *pencilPMEProcessors;


class PmeGridMsg : public CMessage_PmeGridMsg {
public:

  int sourceNode;
  int sequence;
  Lattice lattice;
  PmeReduction *evir;
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
  int sequence;
  Lattice lattice;
  int x_start;
  int nx;
  float *qgrid;
};


class PmeUntransMsg : public CMessage_PmeUntransMsg {
public:

  int sourceNode;
  PmeReduction *evir;
  int y_start;
  int ny;
  float *qgrid;

};


// use this idiom since messages don't have copy constructors
struct PmePencilInitMsgData {
  PmeGrid grid;
  int xBlocks, yBlocks, zBlocks;
  CProxy_PmeXPencil xPencil;
  CProxy_PmeYPencil yPencil;
  CProxy_PmeZPencil zPencil;
  CProxy_ComputePmeMgr pmeProxy;
};

class PmePencilInitMsg : public CMessage_PmePencilInitMsg {
public:
   PmePencilInitMsg(PmePencilInitMsgData &d) { data = d; }
   PmePencilInitMsgData data;
};


struct LocalPmeInfo {
  int nx, x_start;
  int ny_after_transpose, y_start_after_transpose;
};


//Assigns gridPeMap and transPeMap to the same set of processors.
void generatePmePeList(int *peMap, int numPes){
  // decide which pes to use by bit reversal
  int i;
  int ncpus = CkNumPes();
  
  // find next highest power of two
  int npow2 = 1;  int nbits = 0;
  while ( npow2 < ncpus ) { npow2 *= 2; nbits += 1; }
  
  // build bit reversal sequence
  SortableResizeArray<int> seq(ncpus);
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
  }
  
  // extract and sort PME locations
  for ( i=0; i<numPes; ++i ) {
    seq[i] = seq[ncpus - numPes + i];
  }
  seq.resize(numPes);
  seq.sort();
  
  for ( i=0; i<numPes; ++i ) 
      peMap[i] = seq[i];

  //peMap[0] = 0;
}

//Assigns gridPeMap and transPeMap to different set of processors.
void generatePmePeList2(int *gridPeMap, int numGridPes, int *transPeMap, int numTransPes){
  // decide which pes to use by bit reversal
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
  int firstTransPe = ncpus - numGridPes - numTransPes;
  if ( firstTransPe < 0 ) {
    firstTransPe = 0;
    // 0 should be first in list, skip if possible
    if ( ncpus > numTransPes ) firstTransPe = 1;
  }
  for ( i=0; i<numTransPes; ++i ) {
    seq2[i] = seq2[firstTransPe + i];
  }
  seq2.resize(numTransPes);
  seq2.sort();
  
  for ( i=0; i<numGridPes; ++i ) 
    gridPeMap[i] = seq[i];

  for ( i=0; i<numTransPes; ++i ) 
    transPeMap[i] = seq2[i];
}

#if USE_TOPOMAP 
//Topology aware PME allocation
bool generateBGLORBPmePeList(int *pemap, int numPes, int *block_pes=0, 
			     int nbpes=0);
#endif

class ComputePmeMgr : public BOCclass {
public:
  friend class ComputePme;
  ComputePmeMgr();
  ~ComputePmeMgr();

  void initialize(CkQdMsg*);
  void initialize_pencils(CkQdMsg*);
  void activate_pencils(CkQdMsg*);
  void recvArrays(CProxy_PmeXPencil, CProxy_PmeYPencil, CProxy_PmeZPencil);

  void sendGrid(void);
  void recvGrid(PmeGridMsg *);
  void gridCalc1(void);
  void sendTransBarrier(void);
  void sendTrans(void);
  void recvTrans(PmeTransMsg *);
  void gridCalc2(void);
  void sendUntrans(void);
  void recvUntrans(PmeUntransMsg *);
  void gridCalc3(void);
  void sendUngrid(void);
  void recvUngrid(PmeGridMsg *);
  void ungridCalc(void);

  void setCompute(ComputePme *c) { pmeCompute = c; c->setMgr(this); }

  //Tells if the current processor is a PME processor or not. Called by NamdCentralLB
  int isPmeProcessor(int p);  

#ifdef NAMD_FFTW
  static CmiNodeLock fftw_plan_lock;
#endif

private:
  CProxy_ComputePmeMgr pmeProxy;
  CProxy_ComputePmeMgr pmeProxyDir;
  ComputePme *pmeCompute;
  PmeGrid myGrid;
  Lattice lattice;
  PmeKSpace *myKSpace;
  float *qgrid;
  float *kgrid;

#ifdef NAMD_FFTW
  fftw_plan forward_plan_x, backward_plan_x;
  rfftwnd_plan forward_plan_yz, backward_plan_yz;
  fftw_complex *work;
#else
  float *work;
#endif

  int fepOn, lesOn, lesFactor, pairOn, selfOn, numGrids;

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
  char *isPmeFlag;
  int grid_count;
  int trans_count;
  int untrans_count;
  int ungrid_count;
  PmeGridMsg **gridmsg_reuse;
  PmeReduction recip_evir[PME_MAX_EVALS];
  PmeReduction recip_evir2[PME_MAX_EVALS];

  int sequence;  // used for priorities
  int useBarrier;
  int sendTransBarrier_received;

  int usePencils;
  int xBlocks, yBlocks, zBlocks;
  CProxy_PmeXPencil xPencil;
  CProxy_PmeYPencil yPencil;
  CProxy_PmeZPencil zPencil;
  char *pencilActive;
  int numPencilsActive;
};

#ifdef NAMD_FFTW
  CmiNodeLock ComputePmeMgr::fftw_plan_lock;
#endif

int isPmeProcessor(int p){ 
  return CProxy_ComputePmeMgr::ckLocalBranch(CpvAccess(BOCclass_group).computePmeMgr)->isPmeProcessor(p);
}

int ComputePmeMgr::isPmeProcessor(int p){ 
  return ( usePencils ? pencilPMEProcessors[p] : isPmeFlag[p] );
}

ComputePmeMgr::ComputePmeMgr() : pmeProxy(thisgroup), 
				 pmeProxyDir(thisgroup), pmeCompute(0) {

  CpvAccess(BOCclass_group).computePmeMgr = thisgroup;

#ifdef USE_COMM_LIB
  ComlibDelegateProxy(&pmeProxy);
#endif

#ifdef NAMD_FFTW
  if ( CmiMyRank() == 0 ) {
    fftw_plan_lock = CmiCreateLock();
  }
#endif

  myKSpace = 0;
  localInfo = new LocalPmeInfo[CkNumPes()];
  gridPeMap = new int[CkNumPes()];
  transPeMap = new int[CkNumPes()];
  recipPeDest = new int[CkNumPes()];
  gridPeOrder = new int[CkNumPes()];
  transPeOrder = new int[CkNumPes()];
  isPmeFlag = new char[CkNumPes()];
  kgrid = 0;
  work = 0;
  grid_count = 0;
  trans_count = 0;
  untrans_count = 0;
  ungrid_count = 0;
  gridmsg_reuse= new PmeGridMsg*[CkNumPes()];
  useBarrier = 0;
  sendTransBarrier_received = 0;
  usePencils = 0;
}


void ComputePmeMgr::recvArrays(
	CProxy_PmeXPencil x, CProxy_PmeYPencil y, CProxy_PmeZPencil z) {
  xPencil = x;  yPencil = y;  zPencil = z;
}


void ComputePmeMgr::initialize(CkQdMsg *msg) {
  delete msg;

  SimParameters *simParams = Node::Object()->simParameters;
  PatchMap *patchMap = PatchMap::Object();

  fepOn = simParams->fepOn;
  numGrids = fepOn ? 2 : 1;
  lesOn = simParams->lesOn;
  useBarrier = simParams->PMEBarrier;
  if ( lesOn ) {
    lesFactor = simParams->lesFactor;
    numGrids = lesFactor;
  }
  selfOn = 0;
  pairOn = simParams->pairInteractionOn;
  if ( pairOn ) {
    selfOn = simParams->pairInteractionSelf;
    if ( selfOn ) pairOn = 0;  // make pairOn and selfOn exclusive
    numGrids = selfOn ? 1 : 3;
  }

  if ( numGrids != 1 || simParams->PMEPencils == 0 ) usePencils = 0;
  else if ( simParams->PMEPencils > 0 ) usePencils = 1;
  else {
    int dimx = simParams->PMEGridSizeX;
    int dimy = simParams->PMEGridSizeY;
    int maxslabs = 1 + (dimx - 1) / simParams->PMEMinSlices;
    if ( maxslabs > CkNumPes() ) maxslabs = CkNumPes();
    int maxpencils = ( simParams->PMEGridSizeX * simParams->PMEGridSizeY
		* simParams->PMEGridSizeZ ) / simParams->PMEMinPoints;
    if ( maxpencils > CkNumPes() ) maxpencils = CkNumPes();
    if ( maxpencils > 3 * maxslabs ) usePencils = 1;
    else usePencils = 0;
  }

  if ( usePencils ) {
    if ( simParams->PMEPencils > 1 ) {
      xBlocks = yBlocks = zBlocks = simParams->PMEPencils;
    } else {
      int nb2 = ( simParams->PMEGridSizeX * simParams->PMEGridSizeY
		* simParams->PMEGridSizeZ ) / simParams->PMEMinPoints;
      if ( nb2 > CkNumPes() ) nb2 = CkNumPes();
      int nb = (int) sqrt((float)nb2);
      xBlocks = zBlocks = nb;
      yBlocks = nb2 / nb;
    }

    int dimx = simParams->PMEGridSizeX;
    int bx = 1 + ( dimx - 1 ) / xBlocks;
    xBlocks = 1 + ( dimx - 1 ) / bx;

    int dimy = simParams->PMEGridSizeY;
    int by = 1 + ( dimy - 1 ) / yBlocks;
    yBlocks = 1 + ( dimy - 1 ) / by;

    int dimz = simParams->PMEGridSizeZ / 2 + 1;  // complex
    int bz = 1 + ( dimz - 1 ) / zBlocks;
    zBlocks = 1 + ( dimz - 1 ) / bz;

    if ( ! CkMyPe() ) {
      iout << iINFO << "PME using " << xBlocks << " x " <<
        yBlocks << " x " << zBlocks <<
        " pencil grid for FFT and reciprocal sum.\n" << endi;
    }
  } else { // usePencils

  {  // decide how many pes to use for reciprocal sum

    // rules based on work available
    int minslices = simParams->PMEMinSlices;
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
    int i;
    for ( i = 0; i < numGridPes; ++i ) {
      gridPeOrder[i] = i;
    }
    for ( i = 0; i < numTransPes; ++i ) {
      transPeOrder[i] = i;
    }
    Random rand(CkMyPe());
    rand.reorder(gridPeOrder,numGridPes);
    rand.reorder(transPeOrder,numTransPes);
  }

  int sum_npes = numTransPes + numGridPes;
  int max_npes = (numTransPes > numGridPes)?numTransPes:numGridPes;

#if USE_TOPOMAP
  PatchMap * pmap = PatchMap::Object();
  
  int patch_pes = pmap->numNodesWithPatches();
  TopoManager tmgr;
  if(tmgr.hasMultipleProcsPerNode())
    patch_pes *= 2;
 
  bool done = false;
#ifndef USE_COMM_LIB  
  if(CkNumPes() > 2*sum_npes + patch_pes) {    
    done = generateBGLORBPmePeList(transPeMap, numTransPes);
    done &= generateBGLORBPmePeList(gridPeMap, numGridPes, transPeMap, numTransPes);    
  }
  else 
#endif
    if(CkNumPes() > 2 *max_npes + patch_pes) {
      done = generateBGLORBPmePeList(transPeMap, max_npes);
      gridPeMap = transPeMap;
    }

  if (!done)
#endif
    {
      //generatePmePeList(transPeMap, max_npes);
      //gridPeMap = transPeMap;
      generatePmePeList2(gridPeMap, numGridPes, transPeMap, numTransPes);
    }
  
#ifdef USE_COMM_LIB  
  if(CkMyPe() == 0) {
      ComlibInstanceHandle cinst1 = CkCreateComlibInstance();
      EachToManyMulticastStrategy *strat = new 
          EachToManyMulticastStrategy(USE_DIRECT, numGridPes, 
                                      gridPeMap, numTransPes, transPeMap);
      cinst1.setStrategy(strat);
      
      ComlibInstanceHandle cinst2 = CkCreateComlibInstance();
      strat = new EachToManyMulticastStrategy(USE_DIRECT, numTransPes, transPeMap
                                              , numGridPes, gridPeMap);
      cinst2.setStrategy(strat);
      ComlibDoneCreating();
  }
#endif

  myGridPe = -1;
  int i = 0;
  for ( i=0; i<CkNumPes(); ++i )
    isPmeFlag[i] = 0;
  for ( i=0; i<numGridPes; ++i ) {
    if ( gridPeMap[i] == CkMyPe() ) myGridPe = i;
    isPmeFlag[gridPeMap[i]] |= 1;
  }
  myTransPe = -1;
  for ( i=0; i<numTransPes; ++i ) {
    if ( transPeMap[i] == CkMyPe() ) myTransPe = i;
    isPmeFlag[transPeMap[i]] |= 2;
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

  } // ! usePencils

  myGrid.K1 = simParams->PMEGridSizeX;
  myGrid.K2 = simParams->PMEGridSizeY;
  myGrid.K3 = simParams->PMEGridSizeZ;
  myGrid.order = simParams->PMEInterpOrder;
  myGrid.dim2 = myGrid.K2;
  myGrid.dim3 = 2 * (myGrid.K3/2 + 1);

  if ( ! usePencils ) {
    myGrid.block1 = ( myGrid.K1 + numGridPes - 1 ) / numGridPes;
    myGrid.block2 = ( myGrid.K2 + numTransPes - 1 ) / numTransPes;
    myGrid.block3 = myGrid.dim3 / 2;  // complex
  }

  if ( usePencils ) {
    myGrid.block1 = ( myGrid.K1 + xBlocks - 1 ) / xBlocks;
    myGrid.block2 = ( myGrid.K2 + yBlocks - 1 ) / yBlocks;
    myGrid.block3 = ( myGrid.K3/2 + 1 + zBlocks - 1 ) / zBlocks;  // complex

    if ( CkMyPe() == 0 ) {
      int basepe = 0;  int npe = CkNumPes();
      if ( npe > xBlocks*yBlocks &&
		npe > xBlocks*zBlocks &&
		npe > yBlocks*zBlocks ) {
        // avoid node 0
        ++basepe;
        --npe;
      }

      zPencil = CProxy_PmeZPencil::ckNew();  // (xBlocks,yBlocks,1);
      yPencil = CProxy_PmeYPencil::ckNew();  // (xBlocks,1,zBlocks);
      xPencil = CProxy_PmeXPencil::ckNew();  // (1,yBlocks,zBlocks);
      
#if 1

      PatchMap *pmap = PatchMap::Object();
      int npatches = pmap->numHomePatches();

      int *pmemap = new int [npe];
      memset (pmemap, 0, sizeof (int) * npe);
      
      //Use max of x*y, y*z, z*x 
      int n_pme_pes = xBlocks * yBlocks; 
      if ( n_pme_pes < xBlocks * zBlocks ) n_pme_pes = xBlocks * zBlocks;
      if ( n_pme_pes < yBlocks * zBlocks ) n_pme_pes = yBlocks * zBlocks;
      int n_avail_pes = 0;
      
      //Grab all processors where we can store pme chares
      if (npe > npatches + 2 * n_pme_pes) {
	//Use non patch processors to assign pme chares as we have
	//many processors, check base nodes later
	for (int count = 0; count < npe; count++)
	  if(pmap->numPatchesOnNode(basepe + count) == 0)
	    pmemap[n_avail_pes++] = basepe + count;      
      }
      else {  //Use all processors to assign pme chares
	for (int count = 0; count < npe; count++)
	  pmemap [n_avail_pes ++] = basepe + count;
      }

      double pe = 0.0;
      double stride = 1.0; 
            
      int x,y,z;      
      stride = (1.0 * n_avail_pes) / (n_pme_pes);

      pencilPMEProcessors = new char [CkNumPes()];
      memset (pencilPMEProcessors, 0, sizeof(char) * CkNumPes());

      for (pe=0.0, x = 0; x < xBlocks; x ++)
	for (y = 0; y < yBlocks; y ++) {
	  if (pe >= n_avail_pes) pe = 0.0;    
	  zPencil(x,y,0).insert (pmemap[(int) pe]);
	  pencilPMEProcessors [pmemap[(int) pe]] = 1;
	  pe += stride;
	}
      zPencil.doneInserting();
      
      for (pe=1.0, z = 0; z < zBlocks; z ++)
	for (x = 0; x < xBlocks; x ++) {
	  if (pe >= n_avail_pes) pe = 1.0;
	  yPencil(x,0,z).insert (pmemap[(int) pe]);
	  pencilPMEProcessors [pmemap[(int) pe]] = 1;
	  pe += stride;
	}
      yPencil.doneInserting();
      
      for (pe=0.0, y = 0; y < yBlocks; y ++)	
	for (z = 0; z < zBlocks; z ++) {
	  if (pe >= n_avail_pes) pe = 0.0;
	  xPencil(0,y,z).insert (pmemap[(int) pe]);
	  pencilPMEProcessors [pmemap[(int) pe]] = 1;
	  pe += stride;
	}
      xPencil.doneInserting();      
	
      delete [] pmemap;

#else
      int pe = 0;

      for ( int i=0; i<xBlocks; ++i )
       for ( int j=0; j<yBlocks; ++j )
        zPencil(i,j,0).insert(basepe + pe++ % npe);
      zPencil.doneInserting();

      for ( int i=0; i<xBlocks; ++i )
       for ( int k=0; k<zBlocks; ++k )
        yPencil(i,0,k).insert(basepe + pe++ % npe);
      yPencil.doneInserting();

      for ( int j=0; j<yBlocks; ++j )
       for ( int k=0; k<zBlocks; ++k )
        xPencil(0,j,k).insert(basepe + pe++ % npe);
      xPencil.doneInserting();

#endif

      pmeProxy.recvArrays(xPencil,yPencil,zPencil);
      PmePencilInitMsgData msgdata;
      msgdata.grid = myGrid;
      msgdata.xBlocks = xBlocks;
      msgdata.yBlocks = yBlocks;
      msgdata.zBlocks = zBlocks;
      msgdata.xPencil = xPencil;
      msgdata.yPencil = yPencil;
      msgdata.zPencil = zPencil;
      msgdata.pmeProxy = pmeProxyDir;
      xPencil.init(new PmePencilInitMsg(msgdata));
      yPencil.init(new PmePencilInitMsg(msgdata));
      zPencil.init(new PmePencilInitMsg(msgdata));
    }
    return;  // continue in initialize_pencils() at next startup stage
  }


  int pe;
  int nx = 0;
  for ( pe = 0; pe < numGridPes; ++pe ) {
    localInfo[pe].x_start = nx;
    nx += myGrid.block1;
    if ( nx > myGrid.K1 ) nx = myGrid.K1;
    localInfo[pe].nx = nx - localInfo[pe].x_start;
  }
  int ny = 0;
  for ( pe = 0; pe < numTransPes; ++pe ) {
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

#if 0
  if ( numSources ) {
    iout << iINFO << "PME " << CkMyPe() << " sources:";
    for ( node=0; node<numNodes; ++node ) {
      if ( source_flags[node] ) iout << " " << node;
    }
    iout << "\n" << endi;
  }
#endif

  delete [] source_flags;

  // CkPrintf("PME on node %d has %d sources and %d destinations\n",
  //           CkMyPe(), numSources, numDestRecipPes);

  }  // decide how many pes this node exchanges charges with (end)

  ungrid_count = numDestRecipPes;

  sendTransBarrier_received = 0;

  if ( myGridPe < 0 && myTransPe < 0 ) return;
  // the following only for nodes doing reciprocal sum

  if ( myTransPe >= 0 ) {
      int k2_start = localInfo[myTransPe].y_start_after_transpose;
      int k2_end = k2_start + localInfo[myTransPe].ny_after_transpose;
      myKSpace = new PmeKSpace(myGrid, k2_start, k2_end, 0, myGrid.dim3/2);
  }

  int local_size = myGrid.block1 * myGrid.K2 * myGrid.dim3;
  int local_size_2 = myGrid.block2 * myGrid.K1 * myGrid.dim3;
  if ( local_size < local_size_2 ) local_size = local_size_2;
  qgrid = new float[local_size*numGrids];
  if ( numGridPes > 1 || numTransPes > 1 ) {
    kgrid = new float[local_size*numGrids];
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
  CmiLock(fftw_plan_lock);

  work = new fftw_complex[n[0]];

  if ( ! CkMyPe() ) iout << iINFO << "Optimizing 4 FFT steps.  1..." << endi;
  if ( myGridPe >= 0 ) {
  forward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_REAL_TO_COMPLEX,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
  }
  if ( ! CkMyPe() ) iout << " 2..." << endi;
  if ( myTransPe >= 0 ) {
  forward_plan_x = fftw_create_plan_specific(n[0], FFTW_REAL_TO_COMPLEX,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) kgrid,
	localInfo[myTransPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
  }
  if ( ! CkMyPe() ) iout << " 3..." << endi;
  if ( myTransPe >= 0 ) {
  backward_plan_x = fftw_create_plan_specific(n[0], FFTW_COMPLEX_TO_REAL,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) kgrid,
	localInfo[myTransPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
  }
  if ( ! CkMyPe() ) iout << " 4..." << endi;
  if ( myGridPe >= 0 ) {
  backward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_COMPLEX_TO_REAL,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
  }
  if ( ! CkMyPe() ) iout << "   Done.\n" << endi;

  CmiUnlock(fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  if ( myGridPe >= 0 && numSources == 0 )
		NAMD_bug("PME grid elements exist without sources.");
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_size * numGrids * sizeof(float) );
  trans_count = numGridPes;
}


void ComputePmeMgr::initialize_pencils(CkQdMsg *msg) {
  delete msg;
  if ( ! usePencils ) return;

  SimParameters *simParams = Node::Object()->simParameters;

  PatchMap *patchMap = PatchMap::Object();
  Lattice lattice = simParams->lattice;
  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal cutoff = simParams->cutoff;
  BigReal patchdim = simParams->patchDimension;
  int numPatches = patchMap->numPatches();

  pencilActive = new char[xBlocks*yBlocks];
  for ( int i=0; i<xBlocks; ++i ) {
    for ( int j=0; j<yBlocks; ++j ) {
      pencilActive[i*yBlocks+j] = 0;
    }
  }

  for ( int pid=0; pid < numPatches; ++pid ) {
    int pnode = patchMap->node(pid);
    if ( pnode != CkMyPe() ) continue;

    BigReal minx = patchMap->min_a(pid);
    BigReal maxx = patchMap->max_a(pid);
    BigReal margina = 0.5 * ( patchdim - cutoff ) / sysdima;
    // min1 (max1) is smallest (largest) grid line for this patch
    int min1 = ((int) floor(myGrid.K1 * (minx - margina))) - myGrid.order + 1;
    int max1 = ((int) floor(myGrid.K1 * (maxx + margina)));

    BigReal miny = patchMap->min_b(pid);
    BigReal maxy = patchMap->max_b(pid);
    BigReal marginb = 0.5 * ( patchdim - cutoff ) / sysdimb;
    // min2 (max2) is smallest (largest) grid line for this patch
    int min2 = ((int) floor(myGrid.K2 * (miny - marginb))) - myGrid.order + 1;
    int max2 = ((int) floor(myGrid.K2 * (maxy + marginb)));

    for ( int i=min1; i<=max1; ++i ) {
      int ix = i;
      while ( ix >= myGrid.K1 ) ix -= myGrid.K1;
      while ( ix < 0 ) ix += myGrid.K1;
      for ( int j=min2; j<=max2; ++j ) {
        int jy = j;
        while ( jy >= myGrid.K2 ) jy -= myGrid.K2;
        while ( jy < 0 ) jy += myGrid.K2;
        pencilActive[(ix / myGrid.block1)*yBlocks + (jy / myGrid.block2)] = 1;
      }
    }
  }

  numPencilsActive = 0;
  for ( int i=0; i<xBlocks; ++i ) {
    for ( int j=0; j<yBlocks; ++j ) {
      if ( pencilActive[i*yBlocks+j] ) {
        ++numPencilsActive;
        zPencil(i,j,0).dummyRecvGrid(CkMyPe(),0);
      }
    }
  }
  //if ( numPencilsActive ) {
  //  CkPrintf("node %d sending to %d pencils\n", CkMyPe(), numPencilsActive);
  //}

  ungrid_count = numPencilsActive;
}


void ComputePmeMgr::activate_pencils(CkQdMsg *msg) {
  if ( ! usePencils ) return;
  if ( CkMyPe() == 0 ) zPencil.dummyRecvGrid(CkMyPe(),1);
}


ComputePmeMgr::~ComputePmeMgr() {

#ifdef NAMD_FFTW
  if ( CmiMyRank() == 0 ) {
    CmiDestroyLock(fftw_plan_lock);
  }
#endif

  delete myKSpace;
  delete [] localInfo;
  delete [] gridPeMap;
  delete [] transPeMap;
  delete [] recipPeDest;
  delete [] gridPeOrder;
  delete [] transPeOrder;
  delete [] isPmeFlag;
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
    sequence = msg->sequence;
  }

  int zdim = myGrid.dim3;
  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  float *qmsg = msg->qgrid;
  for ( int g=0; g<numGrids; ++g ) {
    char *f = msg->fgrid + fgrid_len * g;
    float *q = qgrid + qgrid_size * g;
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
    pmeProxyDir[CkMyPe()].gridCalc1();
    if ( useBarrier ) pmeProxyDir[0].sendTransBarrier();
#else
    pmeProxyDir.gridCalc1(CkMyPe());
    if ( useBarrier ) pmeProxyDir.sendTransBarrier(0);
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
  if ( ! useBarrier ) pmeProxyDir[CkMyPe()].sendTrans();
#else
  if ( ! useBarrier ) pmeProxyDir.sendTrans(CkMyPe());
#endif
}

void ComputePmeMgr::sendTransBarrier(void) {
  sendTransBarrier_received += 1;
  // CkPrintf("sendTransBarrier on %d %d\n",myGridPe,numGridPes-sendTransBarrier_received);
  if ( sendTransBarrier_received < numGridPes ) return;
  sendTransBarrier_received = 0;
  for ( int i=0; i<numGridPes; ++i ) {
#if CHARM_VERSION > 050402
    pmeProxyDir[gridPeMap[i]].sendTrans();
#else
    pmeProxyDir.sendTrans(gridPeMap[i]);
#endif
  }
}

void ComputePmeMgr::sendTrans(void) {
  // CkPrintf("sendTrans on %d\n",myTransPe);

  // send data for transpose
  int zdim = myGrid.dim3;
  int nx = localInfo[myGridPe].nx;
  int x_start = localInfo[myGridPe].x_start;
  int slicelen = myGrid.K2 * zdim;

#ifdef USE_COMM_LIB
  ComlibInstanceHandle cinst1 = CkGetComlibInstance(0);
  cinst1.beginIteration();
#endif

#if CMK_VERSION_BLUEGENE
  CmiNetworkProgressAfter (0);
#endif

  for (int j=0; j<numTransPes; j++) {
    int pe = transPeOrder[j];  // different order on each node
    LocalPmeInfo &li = localInfo[pe];
    int cpylen = li.ny_after_transpose * zdim;
    PmeTransMsg *newmsg = new (nx * cpylen * numGrids,
				PRIORITY_SIZE) PmeTransMsg;
    newmsg->sourceNode = myGridPe;
    newmsg->lattice = lattice;
    newmsg->x_start = x_start;
    newmsg->nx = nx;
    for ( int g=0; g<numGrids; ++g ) {
      float *q = qgrid + qgrid_size * g + li.y_start_after_transpose * zdim;
      float *qmsg = newmsg->qgrid + nx * cpylen * g;
      for ( int x = 0; x < nx; ++x ) {
        CmiMemcpy((void*)qmsg, (void*)q, cpylen*sizeof(float));
        q += slicelen;
        qmsg += cpylen;
      }
    }
    newmsg->sequence = sequence;
    SET_PRIORITY(newmsg,sequence,PME_TRANS_PRIORITY)
#if CHARM_VERSION > 050402
    pmeProxy[transPeMap[pe]].recvTrans(newmsg);
#else
    pmeProxy.recvTrans(newmsg,transPeMap[pe]);
#endif
  }
 
  untrans_count = numTransPes;

#ifdef USE_COMM_LIB
  cinst1.endIteration();
#endif  
}

void ComputePmeMgr::recvTrans(PmeTransMsg *msg) {
  // CkPrintf("recvTrans on Pe(%d)\n",CkMyPe());
  if ( trans_count == numGridPes ) {
    lattice = msg->lattice;
    sequence = msg->sequence;
  }

  int zdim = myGrid.dim3;
  // int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;
  int x_start = msg->x_start;
  int nx = msg->nx;
  for ( int g=0; g<numGrids; ++g ) {
    CmiMemcpy((void*)(kgrid + qgrid_size * g + x_start*ny*zdim),
	(void*)(msg->qgrid + nx*ny*zdim*g), nx*ny*zdim*sizeof(float));
  }

  delete msg;
  --trans_count;

  if ( trans_count == 0 ) {
#if CHARM_VERSION > 050402
    pmeProxyDir[CkMyPe()].gridCalc2();
#else
    pmeProxyDir.gridCalc2(CkMyPe());
#endif
  }
}

void ComputePmeMgr::gridCalc2(void) {
  // CkPrintf("gridCalc2 on Pe(%d)\n",CkMyPe());

#if CMK_VERSION_BLUEGENE
  CmiNetworkProgressAfter (0);
#endif

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
    recip_evir2[g][0] = myKSpace->compute_energy(kgrid+qgrid_size*g,
			lattice, ewaldcof, &(recip_evir2[g][1]));
    // CkPrintf("Ewald reciprocal energy = %f\n", recip_evir2[g][0]);

    // start backward FFT (x dimension)
#ifdef NAMD_FFTW
    fftw(backward_plan_x, ny * zdim / 2, (fftw_complex *)(kgrid+qgrid_size*g),
	ny * zdim / 2, 1, work, 1, 0);
#endif
  }

#if CHARM_VERSION > 050402
  pmeProxyDir[CkMyPe()].sendUntrans();
#else
  pmeProxyDir.sendUntrans(CkMyPe());
#endif
}

void ComputePmeMgr::sendUntrans(void) {

  int zdim = myGrid.dim3;
  int y_start = localInfo[myTransPe].y_start_after_transpose;
  int ny = localInfo[myTransPe].ny_after_transpose;

#ifdef USE_COMM_LIB
  ComlibInstanceHandle cinst2 = CkGetComlibInstance(1); 
  cinst2.beginIteration();
#endif  

#if CMK_VERSION_BLUEGENE
  CmiNetworkProgressAfter (0);
#endif

  // send data for reverse transpose
  for (int j=0; j<numGridPes; j++) {
    int pe = gridPeOrder[j];  // different order on each node
    LocalPmeInfo &li = localInfo[pe];
    int x_start =li.x_start;
    int nx = li.nx;
    PmeUntransMsg *newmsg = new (nx*ny*zdim*numGrids,numGrids,
				PRIORITY_SIZE) PmeUntransMsg;
    newmsg->sourceNode = myTransPe;
    newmsg->y_start = y_start;
    newmsg->ny = ny;
    for ( int g=0; g<numGrids; ++g ) {
      if ( j == 0 ) {  // only need these once
        newmsg->evir[g] = recip_evir2[g];
      } else {
        newmsg->evir[g] = 0.;
      }
      CmiMemcpy((void*)(newmsg->qgrid+nx*ny*zdim*g),
		(void*)(kgrid + qgrid_size*g + x_start*ny*zdim),
		nx*ny*zdim*sizeof(float));
    }
    SET_PRIORITY(newmsg,sequence,PME_UNTRANS_PRIORITY)
#if CHARM_VERSION > 050402
    pmeProxy[gridPeMap[pe]].recvUntrans(newmsg);
#else
    pmeProxy.recvUntrans(newmsg,gridPeMap[pe]);
#endif
  }

#ifdef USE_COMM_LIB
  cinst2.endIteration();
#endif  

  trans_count = numGridPes;
}

void ComputePmeMgr::recvUntrans(PmeUntransMsg *msg) {
  // CkPrintf("recvUntrans on Pe(%d)\n",CkMyPe());
  if ( untrans_count == numTransPes ) {
    for ( int g=0; g<numGrids; ++g ) {
      recip_evir[g] = 0.;
    }
  }

#if CMK_VERSION_BLUEGENE
  CmiNetworkProgressAfter (0);
#endif

  int g;
  for ( g=0; g<numGrids; ++g ) {
    recip_evir[g] += msg->evir[g];
  }

  int zdim = myGrid.dim3;
  // int x_start = localInfo[myGridPe].x_start;
  int nx = localInfo[myGridPe].nx;
  int y_start = msg->y_start;
  int ny = msg->ny;
  int slicelen = myGrid.K2 * zdim;
  int cpylen = ny * zdim;
  for ( g=0; g<numGrids; ++g ) {
    float *q = qgrid + qgrid_size * g + y_start * zdim;
    float *qmsg = msg->qgrid + nx * cpylen * g;
    for ( int x = 0; x < nx; ++x ) {
      CmiMemcpy((void*)q, (void*)qmsg, cpylen*sizeof(float));
      q += slicelen;
      qmsg += cpylen;
    }
  }

  delete msg;
  --untrans_count;

  if ( untrans_count == 0 ) {
#if CHARM_VERSION > 050402
    pmeProxyDir[CkMyPe()].gridCalc3();
#else
    pmeProxyDir.gridCalc3(CkMyPe());
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
  pmeProxyDir[CkMyPe()].sendUngrid();
#else
  pmeProxyDir.sendUngrid(CkMyPe());
#endif
}

void ComputePmeMgr::sendUngrid(void) {

  for ( int j=0; j<numSources; ++j ) {
    // int msglen = qgrid_len;
    PmeGridMsg *newmsg = gridmsg_reuse[j];
    int pe = newmsg->sourceNode;
    if ( j == 0 ) {  // only need these once
      for ( int g=0; g<numGrids; ++g ) {
        newmsg->evir[g] = recip_evir[g];
      }
    } else {
      for ( int g=0; g<numGrids; ++g ) {
        newmsg->evir[g] = 0.;
      }
    }
    int zdim = myGrid.dim3;
    int flen = newmsg->len;
    int fstart = newmsg->start;
    int zlistlen = newmsg->zlistlen;
    int *zlist = newmsg->zlist;
    float *qmsg = newmsg->qgrid;
    for ( int g=0; g<numGrids; ++g ) {
      char *f = newmsg->fgrid + fgrid_len * g;
      float *q = qgrid + qgrid_size * g + (fstart-fgrid_start) * zdim;
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

    SET_PRIORITY(newmsg,sequence,PME_UNGRID_PRIORITY)
#if CHARM_VERSION > 050402
    pmeProxyDir[pe].recvUngrid(newmsg);
#else
    pmeProxyDir.recvUngrid(newmsg,pe);
#endif
  }
  grid_count = numSources;
  memset( (void*) qgrid, 0, qgrid_size * numGrids * sizeof(float) );
}

void ComputePmeMgr::recvUngrid(PmeGridMsg *msg) {
  // CkPrintf("recvUngrid on Pe(%d)\n",CkMyPe());
  if ( ungrid_count == 0 ) {
    NAMD_bug("Message order failure in ComputePmeMgr::recvUngrid\n");
  }

  if ( usePencils ) pmeCompute->copyPencils(msg);
  else pmeCompute->copyResults(msg);
  delete msg;
  --ungrid_count;

  if ( ungrid_count == 0 ) {
#if CHARM_VERSION > 050402
    pmeProxyDir[CkMyPe()].ungridCalc();
#else
    pmeProxyDir.ungridCalc(CkMyPe());
#endif
  }
}

void ComputePmeMgr::ungridCalc(void) {
  // CkPrintf("ungridCalc on Pe(%d)\n",CkMyPe());

  pmeCompute->ungridForces();

  ungrid_count = (usePencils ? numPencilsActive : numDestRecipPes );
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
  selfOn = 0;
  pairOn = simParams->pairInteractionOn;
  if ( pairOn ) {
    selfOn = simParams->pairInteractionSelf;
    if ( selfOn ) pairOn = 0;  // make pairOn and selfOn exclusive
    numGrids = selfOn ? 1 : 3;
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

  localData = new PmeParticle[numLocalAtoms*(numGrids+
					((numGrids>1 || selfOn)?1:0))];
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
#ifdef NETWORK_PROGRESS
    CmiNetworkProgress();
#endif

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
#ifdef NETWORK_PROGRESS
      CmiNetworkProgress();
#endif

      PmeParticle *lgd = localGridData[g];
      int nga = 0;
      for(int i=0; i<numLocalAtoms; ++i) {
        if ( localPartition[i] == 0 || localPartition[i] == (g+1) ) {
          lgd[nga++] = localData[i];
        }
      }
      numGridAtoms[g] = nga;
    }
  } else if ( selfOn ) {
    if ( numGrids != 1 ) NAMD_bug("ComputePme::doWork assertion 1 failed");
    g = 0;
    PmeParticle *lgd = localGridData[g];
    int nga = 0;
    for(int i=0; i<numLocalAtoms; ++i) {
      if ( localPartition[i] == 1 ) {
        lgd[nga++] = localData[i];
      }
    }
    numGridAtoms[g] = nga;
  } else if ( pairOn ) {
    if ( numGrids != 3 ) NAMD_bug("ComputePme::doWork assertion 2 failed");
    g = 0;
    PmeParticle *lgd = localGridData[g];
    int nga = 0;
    for(int i=0; i<numLocalAtoms; ++i) {
      if ( localPartition[i] == 1 || localPartition[i] == 2 ) {
        lgd[nga++] = localData[i];
      }
    }
    numGridAtoms[g] = nga;
    for ( g=1; g<3; ++g ) {
      PmeParticle *lgd = localGridData[g];
      int nga = 0;
      for(int i=0; i<numLocalAtoms; ++i) {
        if ( localPartition[i] == g ) {
          lgd[nga++] = localData[i];
        }
      }
      numGridAtoms[g] = nga;
    }
  } else {
    if ( numGrids != 1 ) NAMD_bug("ComputePme::doWork assertion 3 failed");
    localGridData[0] = localData;
    numGridAtoms[0] = numLocalAtoms;
  }

  memset( (void*) fz_arr, 0, myGrid.K3 * sizeof(char) );

  // calculate self energy
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  for ( g=0; g<numGrids; ++g ) {
#ifdef NETWORK_PROGRESS
    CmiNetworkProgress();
#endif

    evir[g] = 0;
    BigReal selfEnergy = 0;
    data_ptr = localGridData[g];
    int i;
    for(i=0; i<numGridAtoms[g]; ++i)
    {
      selfEnergy += data_ptr->cg * data_ptr->cg;
      ++data_ptr;
    }
    selfEnergy *= -1. * ewaldcof / SQRT_PI;
    evir[g][0] += selfEnergy;

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

  if ( myMgr->usePencils ) {
    sendPencils();
  } else {
#if 0
  CProxy_ComputePmeMgr pmeProxy(CpvAccess(BOCclass_group).computePmeMgr);
#if CHARM_VERSION > 050402
  pmeProxy[CkMyPe()].sendGrid();
#else
  pmeProxy.sendGrid(CkMyPe());
#endif
#else
  sendData(myMgr->numGridPes,myMgr->gridPeOrder,
		myMgr->recipPeDest,myMgr->gridPeMap);
#endif
  }
}


void ComputePme::sendPencils() {

  // iout << "Sending charge grid for " << numLocalAtoms << " atoms to FFT on " << iPE << ".\n" << endi;

  int xBlocks = myMgr->xBlocks;
  int yBlocks = myMgr->yBlocks;
  int zBlocks = myMgr->zBlocks;
  myGrid.block1 = ( myGrid.K1 + xBlocks - 1 ) / xBlocks;
  myGrid.block2 = ( myGrid.K2 + yBlocks - 1 ) / yBlocks;
  int K1 = myGrid.K1;
  int K2 = myGrid.K2;
  int dim2 = myGrid.dim2;
  int dim3 = myGrid.dim3;
  int block1 = myGrid.block1;
  int block2 = myGrid.block2;

  Lattice lattice = patchList[0].p->flags.lattice;

  resultsRemaining = myMgr->numPencilsActive;
  const char *pencilActive = myMgr->pencilActive;

  strayChargeErrors = 0;

  for (int ib=0; ib<xBlocks; ++ib) {
   for (int jb=0; jb<yBlocks; ++jb) {
    int ibegin = ib*block1;
    int iend = ibegin + block1;  if ( iend > K1 ) iend = K1;
    int jbegin = jb*block2;
    int jend = jbegin + block2;  if ( jend > K2 ) jend = K2;
    int flen = numGrids * (iend - ibegin) * (jend - jbegin);
    int fcount = 0;

    for ( int g=0; g<numGrids; ++g ) {
      char *f = f_arr + g*fsize;
      int fcount_g = 0;
      for ( int i=ibegin; i<iend; ++i ) {
       for ( int j=jbegin; j<jend; ++j ) {
        fcount_g += ( f[i*dim2+j] ? 1 : 0 );
       }
      }
      fcount += fcount_g;
      if ( ! pencilActive[ib*yBlocks+jb] ) {
        if ( fcount_g ) {
          ++strayChargeErrors;
          iout << iERROR << "Stray PME grid charges detected: "
		<< CkMyPe() << " sending to (x,y)";
          for ( int i=ibegin; i<iend; ++i ) {
           for ( int j=jbegin; j<jend; ++j ) {
            if ( f[i*dim2+j] ) { iout << " (" << i << "," << j << ")"; }
           }
          }
          iout << "\n" << endi;
        }
      }
    }

#ifdef NETWORK_PROGRESS
    CmiNetworkProgress();
#endif

    if ( ! pencilActive[ib*yBlocks+jb] ) continue;

    int zlistlen = 0;
    for ( int i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) ++zlistlen;
    }

    PmeGridMsg *msg = new (fcount*zlistlen, zlistlen, flen,
	numGrids, PRIORITY_SIZE) PmeGridMsg;
    msg->sourceNode = CkMyPe();
    msg->lattice = lattice;
#if 0
    msg->start = fstart;
    msg->len = flen;
#else
    msg->start = -1;   // obsolete?
    msg->len = -1;   // obsolete?
#endif
    msg->zlistlen = zlistlen;
    int *zlist = msg->zlist;
    zlistlen = 0;
    for ( int i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) zlist[zlistlen++] = i;
    }
    char *fmsg = msg->fgrid;
    float *qmsg = msg->qgrid;
    for ( int g=0; g<numGrids; ++g ) {
      char *f = f_arr + g*fsize;
      double **q = q_arr + g*fsize;
      for ( int i=ibegin; i<iend; ++i ) {
       for ( int j=jbegin; j<jend; ++j ) {
        *(fmsg++) = f[i*dim2+j];
        if( f[i*dim2+j] ) {
          for ( int k=0; k<zlistlen; ++k ) {
            *(qmsg++) = q[i*dim2+j][zlist[k]];
          }
        }
       }
      }
    }

    msg->sequence = sequence();
    SET_PRIORITY(msg,sequence(),PME_GRID_PRIORITY)
    myMgr->zPencil(ib,jb,0).recvGrid(msg);
   }
  }

  for (int i=0; i<fsize; ++i) {
    if ( q_arr[i] ) {
      memset( (void*) (q_arr[i]), -1, myGrid.dim3 * sizeof(double) );
    }
  }

}


void ComputePme::copyPencils(PmeGridMsg *msg) {

  int xBlocks = myMgr->xBlocks;
  int yBlocks = myMgr->yBlocks;
  int zBlocks = myMgr->zBlocks;
  myGrid.block1 = ( myGrid.K1 + xBlocks - 1 ) / xBlocks;
  myGrid.block2 = ( myGrid.K2 + yBlocks - 1 ) / yBlocks;
  int K1 = myGrid.K1;
  int K2 = myGrid.K2;
  int dim2 = myGrid.dim2;
  int dim3 = myGrid.dim3;
  int block1 = myGrid.block1;
  int block2 = myGrid.block2;

  // msg->sourceNode = thisIndex.x * initdata.yBlocks + thisIndex.y;
  int ib = msg->sourceNode / yBlocks;
  int jb = msg->sourceNode % yBlocks;

  int ibegin = ib*block1;
  int iend = ibegin + block1;  if ( iend > K1 ) iend = K1;
  int jbegin = jb*block2;
  int jend = jbegin + block2;  if ( jend > K2 ) jend = K2;

  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  float *qmsg = msg->qgrid;
  int g;
  for ( g=0; g<numGrids; ++g ) {
    evir[g] += msg->evir[g];
    char *f = f_arr + g*fsize;
    double **q = q_arr + g*fsize;
    for ( int i=ibegin; i<iend; ++i ) {
     for ( int j=jbegin; j<jend; ++j ) {
      if( f[i*dim2+j] ) {
        for ( int k=0; k<zlistlen; ++k ) {
          q[i*dim2+j][zlist[k]] = *(qmsg++);
        }
      }
     }
    }
  }
}


void ComputePme::sendData(int numRecipPes, int *recipPeOrder,
				int *recipPeDest, int *gridPeMap) {

  // iout << "Sending charge grid for " << numLocalAtoms << " atoms to FFT on " << iPE << ".\n" << endi;

  myGrid.block1 = ( myGrid.K1 + numRecipPes - 1 ) / numRecipPes;
  myGrid.block2 = ( myGrid.K2 + numRecipPes - 1 ) / numRecipPes;
  bsize = myGrid.block1 * myGrid.dim2 * myGrid.dim3;

  Lattice lattice = patchList[0].p->flags.lattice;

  resultsRemaining = numRecipPes;

  strayChargeErrors = 0;

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
          ++strayChargeErrors;
          iout << iERROR << "Stray PME grid charges detected: "
		<< CkMyPe() << " sending to " << gridPeMap[pe] << " for planes";
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

#ifdef NETWORK_PROGRESS
    CmiNetworkProgress();
#endif

    if ( ! recipPeDest[pe] ) continue;

    int zlistlen = 0;
    for ( i=0; i<myGrid.K3; ++i ) {
      if ( fz_arr[i] ) ++zlistlen;
    }

    PmeGridMsg *msg = new (fcount*zlistlen, zlistlen, flen*numGrids,
				numGrids, PRIORITY_SIZE) PmeGridMsg;
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
      CmiMemcpy((void*)(msg->fgrid+g*flen),(void*)f,flen*sizeof(char));
      double **q = q_arr + fstart + g*fsize;
      for ( i=0; i<flen; ++i ) {
        if ( f[i] ) {
          for ( int k=0; k<zlistlen; ++k ) {
            *(qmsg++) = q[i][zlist[k]];
          }
        }
      }
    }

    msg->sequence = sequence();
    SET_PRIORITY(msg,sequence(),PME_GRID_PRIORITY)
#if CHARM_VERSION > 050402
    pmeProxy[gridPeMap[pe]].recvGrid(msg);
#else
    pmeProxy.recvGrid(msg,gridPeMap[pe]);
#endif
  }

  for (int i=0; i<fsize; ++i) {
    if ( q_arr[i] ) {
      memset( (void*) (q_arr[i]), -1, myGrid.dim3 * sizeof(double) );
    }
  }

}

void ComputePme::copyResults(PmeGridMsg *msg) {

  int zdim = myGrid.dim3;
  int flen = msg->len;
  int fstart = msg->start;
  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  float *qmsg = msg->qgrid;
  int g;
  for ( g=0; g<numGrids; ++g ) {
    evir[g] += msg->evir[g];
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

    Vector *localResults = new Vector[numLocalAtoms*
					((numGrids>1 || selfOn)?2:1)];
    Vector *gridResults;
    if ( fepOn || lesOn || selfOn || pairOn ) {
      for(int i=0; i<numLocalAtoms; ++i) { localResults[i] = 0.; }
      gridResults = localResults + numLocalAtoms;
    } else {
      gridResults = localResults;
    }

    Vector pairForce = 0.;
    Lattice lattice = patchList[0].p->flags.lattice;
    int g = 0;
    for ( g=0; g<numGrids; ++g ) {
#ifdef NETWORK_PROGRESS
      CmiNetworkProgress();
#endif

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
      } else if ( selfOn ) {
        PmeParticle *lgd = localGridData[g];
        int nga = 0;
        for(int i=0; i<numLocalAtoms; ++i) {
          if ( localPartition[i] == 1 ) {
            pairForce += gridResults[nga];  // should add up to almost zero
            localResults[i] += gridResults[nga++];
          }
        }
      } else if ( pairOn ) {
        if ( g == 0 ) {
          int nga = 0;
          for(int i=0; i<numLocalAtoms; ++i) {
            if ( localPartition[i] == 1 ) {
              pairForce += gridResults[nga];
            }
            if ( localPartition[i] == 1 || localPartition[i] == 2 ) {
              localResults[i] += gridResults[nga++];
            }
          }
        } else if ( g == 1 ) {
          int nga = 0;
          for(int i=0; i<numLocalAtoms; ++i) {
            if ( localPartition[i] == g ) {
              pairForce -= gridResults[nga];  // should add up to almost zero
              localResults[i] -= gridResults[nga++];
            }
          }
        } else {
          int nga = 0;
          for(int i=0; i<numLocalAtoms; ++i) {
            if ( localPartition[i] == g ) {
              localResults[i] -= gridResults[nga++];
            }
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

      if ( ! strayChargeErrors ) {
        for(int i=0; i<numAtoms; ++i) {
          f[i].x += results_ptr->x;
          f[i].y += results_ptr->y;
          f[i].z += results_ptr->z;
          ++results_ptr;
        }
      }
  
      (*ap).forceBox->close(&r);
    }

    delete [] localResults;
   
    if ( pairOn || selfOn ) {
        ADD_VECTOR_OBJECT(reduction,REDUCTION_PAIR_ELECT_FORCE,pairForce);
    }

    for ( g=0; g<numGrids; ++g ) {
      double scale = 1.;
      if ( fepOn ) {
        if ( g == 0 ) scale = simParams->lambda;
        else if ( g == 1 ) scale = 1. - simParams->lambda;
      } else if ( lesOn ) {
        scale = 1.0 / (double)lesFactor;
      } else if ( pairOn ) {
        scale = ( g == 0 ? 1. : -1. );
      }
      reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += evir[g][0] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XX) += evir[g][1] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XY) += evir[g][2] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += evir[g][3] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YX) += evir[g][2] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YY) += evir[g][4] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += evir[g][5] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += evir[g][3] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += evir[g][5] * scale;
      reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += evir[g][6] * scale;

      double scale2 = 0.;
      if ( fepOn && g == 0 ) scale2 = simParams->lambda2;
      else if ( fepOn && g == 1 ) scale2 = 1. - simParams->lambda2;
      reduction->item(REDUCTION_ELECT_ENERGY_SLOW_F) += evir[g][0] * scale2;
    }
    reduction->item(REDUCTION_STRAY_CHARGE_ERRORS) += strayChargeErrors;
    reduction->submit();

}

#if USE_TOPOMAP 

#define NPRIMES 8
const static unsigned int NAMDPrimes[] = {
  3,
  5,
  7,
  11,
  13,
  17,
  19,
  23,  
  29,
  31,
  37,
  59,
  73,
  93,
  113,
  157,
  307,
  617,
  1217                  //This should b enough for 64K nodes of BGL. 
};

#include "RecBisection.h"

/***-----------------------------------------------------**********
    The Orthogonal Recursive Bisection strategy, which allocates PME
    objects close to the patches they communicate, and at the same
    time spreads them around the grid 
****----------------------------------------------------------****/

bool generateBGLORBPmePeList(int *pemap, int numPes, 
			     int *block_pes, int nbpes) {

  PatchMap *pmap = PatchMap::Object();
  int *pmemap = new int [CkNumPes()];

  if (pemap == NULL)
    return false;

  TopoManager tmgr;

  memset(pmemap, 0, sizeof(int) * CkNumPes());

  for(int count = 0; count < CkNumPes(); count++) {
    if(count < nbpes)
      pmemap[block_pes[count]] = 1;
    
    if(pmap->numPatchesOnNode(count)) {
      pmemap[count] = 1;
      
      //Assumes an XYZT mapping !!
      if(tmgr.hasMultipleProcsPerNode()) {
	pmemap[(count + CkNumPes()/2)% CkNumPes()] = 1;
      }
    }
  }

  if(numPes + nbpes + pmap->numNodesWithPatches() > CkNumPes())
    //NAMD_bug("PME ORB Allocator: Processors Unavailable\n");
    return false;

  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  //first split PME processors into patch groups

  int xsize = 0, ysize = 0, zsize = 0;

  xsize = tmgr.getDimX();
  ysize = tmgr.getDimY();
  zsize = tmgr.getDimZ();
  
  int nx = xsize, ny = ysize, nz = zsize;
  DimensionMap dm;
  
  dm.x = 0;
  dm.y = 1;
  dm.z = 2;
  
  findOptimalDimensions(xsize, ysize, zsize, nx, ny, nz, dm);

  //group size processors have to be allocated to each YZ plane
  int group_size = numPes/nx;
  if(numPes % nx)
    group_size ++;

  int my_prime = NAMDPrimes[0];
  int density = (ny * nz)/group_size + 1;
  int count = 0;
  
  //Choose a suitable prime Number
  for(count = 0; count < NPRIMES; count ++) {
    //Find a prime just greater than the density
    if(density < NAMDPrimes[count]) {
      my_prime = NAMDPrimes[count];
      break;
    }      
  }
  
  if(count == NPRIMES)
    my_prime = NAMDPrimes[NPRIMES-1];

  //int gcount = numPes/2;
  int gcount = 0;
  int npme_pes = 0;
  
  int coord[3];

  for(int x = 0; x < nx; x++) {
    coord[0] = (x + nx/2)%nx;
    
    for(count=0; count < group_size && npme_pes < numPes; count++) {
      int dest = (count + 1) * my_prime;      
      dest = dest % (ny * nz);      
      
      coord[2] = dest / ny;
      coord[1] = dest - coord[2] * ny;
      
      //Locate where in the actual grid the processor is
      int destPe = coord[dm.x] + coord[dm.y] * xsize + 
	coord[dm.z] * xsize* ysize;
      
      if(pmemap[destPe] == 0) {
        pemap[gcount++] = destPe;
        pmemap[destPe] = 1;
	
	if(tmgr.hasMultipleProcsPerNode())
	  pmemap[(destPe + CkNumPes()/2) % CkNumPes()] = 1;	

        npme_pes ++;
      }
      else {
        for(int pos = 1; pos < ny * nz; pos++) {
          
          coord[2] += pos / ny;
          coord[1] += pos % ny;
          
          coord[2] = coord[2] % nz;
          coord[1] = coord[1] % ny;       
          
          int newdest = coord[dm.x] + coord[dm.y] * xsize + 
	    coord[dm.z] * xsize * ysize;
          
          if(pmemap[newdest] == 0) {
            pemap[gcount++] = newdest;
            pmemap[newdest] = 1;
	    
	    if(tmgr.hasMultipleProcsPerNode())
	      pmemap[(newdest + CkNumPes()/2) % CkNumPes()] = 1;	
	    
            npme_pes ++;
            break;
          }
        }
      }      
    }   
    
    if(gcount == numPes)
      gcount = 0;    
    
    if(npme_pes >= numPes)
      break;
  }
  
  delete [] pmemap;
  
  if(npme_pes != numPes)
    //NAMD_bug("ORB PME allocator failed\n");
    return false;

  return true;
}

#endif

template <class T> class PmePencil : public T {
public:
  PmePencil() {
    data = 0;
    work = 0;
  }
  ~PmePencil() {
    delete [] data;
    delete [] work;
  }
  void base_init(PmePencilInitMsg *msg) {
    initdata = msg->data;
  }
  void order_init(int nBlocks) {
    send_order = new int[nBlocks];
    for ( int i=0; i<nBlocks; ++i ) send_order[i] = i;
    Random rand(CkMyPe());
    rand.reorder(send_order,nBlocks);
  }
  PmePencilInitMsgData initdata;
  Lattice lattice;
  PmeReduction evir;
  int sequence;  // used for priorities
  int imsg;  // used in sdag code
  float *data;
  float *work;
  int *send_order;
};

class PmeZPencil : public PmePencil<CBase_PmeZPencil> {
public:
    PmeZPencil_SDAG_CODE;
    PmeZPencil() { __sdag_init(); setMigratable(false); }
    PmeZPencil(CkMigrateMessage *) { __sdag_init();  setMigratable (false); }
    void fft_init();
    void recv_grid(const PmeGridMsg *);
    void forward_fft();
    void send_trans(int dest);
    void recv_untrans(const PmeUntransMsg *);
    void backward_fft();
    void send_ungrid(PmeGridMsg *);
private:
    ResizeArray<PmeGridMsg *> grid_msgs;
#ifdef NAMD_FFTW
    rfftwnd_plan forward_plan, backward_plan;
#endif
    int nx, ny;
};

class PmeYPencil : public PmePencil<CBase_PmeYPencil> {
public:
    PmeYPencil_SDAG_CODE;
    PmeYPencil() { __sdag_init(); setMigratable(false); }
    PmeYPencil(CkMigrateMessage *) { __sdag_init(); }
    void fft_init();
    void recv_trans(const PmeTransMsg *);
    void forward_fft();
    void send_trans(int dest);
    void recv_untrans(const PmeUntransMsg *);
    void backward_fft();
    void send_untrans(int dest);
private:
#ifdef NAMD_FFTW
    fftw_plan forward_plan, backward_plan;
#endif
    int nx, nz;
};

class PmeXPencil : public PmePencil<CBase_PmeXPencil> {
public:
    PmeXPencil_SDAG_CODE;
    PmeXPencil() { __sdag_init();  myKSpace = 0; setMigratable(false); }
    PmeXPencil(CkMigrateMessage *) { __sdag_init(); }
    void fft_init();
    void recv_trans(const PmeTransMsg *);
    void forward_fft();
    void pme_kspace();
    void backward_fft();
    void send_untrans(int dest);
#ifdef NAMD_FFTW
    fftw_plan forward_plan, backward_plan;
#endif
    int ny, nz;
    PmeKSpace *myKSpace;
};

void PmeZPencil::fft_init() {
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int K3 = initdata.grid.K3;
  int dim3 = initdata.grid.dim3;
  int block1 = initdata.grid.block1;
  int block2 = initdata.grid.block2;

  nx = block1;
  if ( (thisIndex.x + 1) * block1 > K1 ) nx = K1 - thisIndex.x * block1;
  ny = block2;
  if ( (thisIndex.y + 1) * block2 > K2 ) ny = K2 - thisIndex.y * block2;

  data = new float[nx*ny*dim3];
  work = new float[dim3];

  order_init(initdata.zBlocks);

#ifdef NAMD_FFTW
  CmiLock(ComputePmeMgr::fftw_plan_lock);

  forward_plan = rfftwnd_create_plan_specific(1, &K3, FFTW_REAL_TO_COMPLEX,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, data, 1, work, 1);
  backward_plan = rfftwnd_create_plan_specific(1, &K3, FFTW_COMPLEX_TO_REAL,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, data, 1, work, 1);

  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif
}

void PmeYPencil::fft_init() {
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int dim2 = initdata.grid.dim2;
  int dim3 = initdata.grid.dim3;
  int block1 = initdata.grid.block1;
  int block3 = initdata.grid.block3;

  nx = block1;
  if ( (thisIndex.x + 1) * block1 > K1 ) nx = K1 - thisIndex.x * block1;
  nz = block3;
  if ( (thisIndex.z+1)*block3 > dim3/2 ) nz = dim3/2 - thisIndex.z*block3;

  data = new float[nx*dim2*nz*2];
  work = new float[2*K2];

  order_init(initdata.yBlocks);

#ifdef NAMD_FFTW
  CmiLock(ComputePmeMgr::fftw_plan_lock);

  forward_plan = fftw_create_plan_specific(K2, FFTW_FORWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	nz, (fftw_complex *) work, 1);
  backward_plan = fftw_create_plan_specific(K2, FFTW_BACKWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	nz, (fftw_complex *) work, 1);

  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif
}

void PmeXPencil::fft_init() {
  CProxy_Node nd(CpvAccess(BOCclass_group).node);
  Node *node = nd.ckLocalBranch();
  SimParameters *simParams = node->simParameters;

  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int dim3 = initdata.grid.dim3;
  int block2 = initdata.grid.block2;
  int block3 = initdata.grid.block3;

  ny = block2;
  if ( (thisIndex.y + 1) * block2 > K2 ) ny = K2 - thisIndex.y * block2;
  nz = block3;
  if ( (thisIndex.z+1)*block3 > dim3/2 ) nz = dim3/2 - thisIndex.z*block3;

  data = new float[K1*block2*block3*2];
  work = new float[2*K1];

  order_init(initdata.xBlocks);

#ifdef NAMD_FFTW
  CmiLock(ComputePmeMgr::fftw_plan_lock);

  forward_plan = fftw_create_plan_specific(K1, FFTW_FORWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	ny*nz, (fftw_complex *) work, 1);
  backward_plan = fftw_create_plan_specific(K1, FFTW_BACKWARD,
	( simParams->FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	| FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) data,
	ny*nz, (fftw_complex *) work, 1);

  CmiUnlock(ComputePmeMgr::fftw_plan_lock);
#else
  NAMD_die("Sorry, FFTW must be compiled in to use PME.");
#endif

  myKSpace = new PmeKSpace(initdata.grid,
		thisIndex.y*block2, thisIndex.y*block2 + ny,
		thisIndex.z*block3, thisIndex.z*block3 + nz);
}

// #define FFTCHECK   // run a grid of integers through the fft
// #define ZEROCHECK  // check for suspicious zeros in fft

void PmeZPencil::recv_grid(const PmeGridMsg *msg) {

  int dim3 = initdata.grid.dim3;
  if ( imsg == 0 ) {
    lattice = msg->lattice;
    sequence = msg->sequence;
    memset(data, 0, sizeof(float) * nx*ny*dim3);
  }

  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  char *fmsg = msg->fgrid;
  float *qmsg = msg->qgrid;
  float *d = data;
  int numGrids = 1;  // pencil FFT doesn't support multiple grids
  for ( int g=0; g<numGrids; ++g ) {
    for ( int i=0; i<nx; ++i ) {
     for ( int j=0; j<ny; ++j, d += dim3 ) {
      if( *(fmsg++) ) {
        for ( int k=0; k<zlistlen; ++k ) {
          d[zlist[k]] += *(qmsg++);
        }
      }
     }
    }
  }
}

void PmeZPencil::forward_fft() {
#ifdef FFTCHECK
  int dim3 = initdata.grid.dim3;
  int K3 = initdata.grid.K3;
  float std_base = 100. * (thisIndex.x+1.) + 10. * (thisIndex.y+1.);
  float *d = data;
  for ( int i=0; i<nx; ++i ) {
   for ( int j=0; j<ny; ++j, d += dim3 ) {
    for ( int k=0; k<dim3; ++k ) {
      d[k] = 10. * (10. * (10. * std_base + i) + j) + k;
    }
   }
  }
#endif
#ifdef NAMD_FFTW
  rfftwnd_real_to_complex(forward_plan, nx*ny,
	data, 1, initdata.grid.dim3, (fftw_complex *) work, 1, 0);
#endif
#ifdef ZEROCHECK
  int dim3 = initdata.grid.dim3;
  int K3 = initdata.grid.K3;
  float *d = data;
  for ( int i=0; i<nx; ++i ) {
   for ( int j=0; j<ny; ++j, d += dim3 ) {
    for ( int k=0; k<dim3; ++k ) {
      if ( d[k] == 0. ) CkPrintf("0 in Z at %d %d %d %d %d %d %d %d %d\n",
	thisIndex.x, thisIndex.y, i, j, k, nx, ny, dim3);
    }
   }
  }
#endif
}

void PmeZPencil::send_trans(int dest) {
  int zBlocks = initdata.zBlocks;
  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;
  for ( int isend=0; isend<zBlocks; ++isend ) {
    int kb = send_order[isend];
    int nz = block3;
    if ( (kb+1)*block3 > dim3/2 ) nz = dim3/2 - kb*block3;
    PmeTransMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) PmeTransMsg;
    msg->lattice = lattice;
    msg->sourceNode = thisIndex.y;
    msg->nx = ny;
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i ) {
     for ( int j=0; j<ny; ++j, d += dim3 ) {
      for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
        *(md++) = d[2*k];
        *(md++) = d[2*k+1];
      }
     }
    }
    msg->sequence = sequence;
    SET_PRIORITY(msg,sequence,PME_TRANS_PRIORITY)
    initdata.yPencil(thisIndex.x,0,kb).recvTrans(msg);
  }
}

void PmeYPencil::recv_trans(const PmeTransMsg *msg) {
  if ( imsg == 0 ) {
    lattice = msg->lattice;
    sequence = msg->sequence;
  }
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  int jb = msg->sourceNode;
  int ny = msg->nx;
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
   for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
    for ( int k=0; k<nz; ++k ) {
#ifdef ZEROCHECK
      if ( (*md) == 0. ) CkPrintf("0 in ZY at %d %d %d %d %d %d %d %d %d\n",
	thisIndex.x, jb, thisIndex.z, i, j, k, nx, ny, nz);
#endif
      d[2*(j*nz+k)] = *(md++);
      d[2*(j*nz+k)+1] = *(md++);
    }
   }
  }
}

void PmeYPencil::forward_fft() {
#ifdef NAMD_FFTW
  for ( int i=0; i<nx; ++i ) {
    fftw(forward_plan, nz,
	((fftw_complex *) data) + i * nz * initdata.grid.K2,
	nz, 1, (fftw_complex *) work, 1, 0);
  }
#endif
}

void PmeYPencil::send_trans(int dest) {
  int yBlocks = initdata.yBlocks;
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  for ( int isend=0; isend<yBlocks; ++isend ) {
    int jb = send_order[isend];
    int ny = block2;
    if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
    PmeTransMsg *msg = new (nx*ny*nz*2,PRIORITY_SIZE) PmeTransMsg;
    msg->lattice = lattice;
    msg->sourceNode = thisIndex.x;
    msg->nx = nx;
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
     for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
        *(md++) = d[2*(j*nz+k)];
        *(md++) = d[2*(j*nz+k)+1];
#ifdef ZEROCHECK
        if ( *(md-2) == 0. ) CkPrintf("send 0 in YX at %d %d %d %d %d %d %d %d %d\n",
	thisIndex.x, jb, thisIndex.z, i, j, k, nx, ny, nz);
#endif
      }
     }
    }
    if ( md != msg->qgrid + nx*ny*nz*2 ) CkPrintf("error in YX at %d %d %d\n",
	thisIndex.x, jb, thisIndex.z);
    msg->sequence = sequence;
    SET_PRIORITY(msg,sequence,PME_TRANS2_PRIORITY)
    initdata.xPencil(0,jb,thisIndex.z).recvTrans(msg);
  }
}

void PmeXPencil::recv_trans(const PmeTransMsg *msg) {
  if ( imsg == 0 ) {
    lattice = msg->lattice;
    sequence = msg->sequence;
  }
  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;
  int ib = msg->sourceNode;
  int nx = msg->nx;
  const float *md = msg->qgrid;
  for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
   float *d = data + i*ny*nz*2;
   for ( int j=0; j<ny; ++j, d += nz*2 ) {
    for ( int k=0; k<nz; ++k ) {
#ifdef ZEROCHECK
      if ( (*md) == 0. ) CkPrintf("0 in YX at %d %d %d %d %d %d %d %d %d\n",
	ib, thisIndex.y, thisIndex.z, i, j, k, nx, ny, nz);
#endif
      d[2*k] = *(md++);
      d[2*k+1] = *(md++);
    }
   }
  }
}

void PmeXPencil::forward_fft() {
#ifdef NAMD_FFTW
  fftw(forward_plan, ny*nz,
	((fftw_complex *) data), ny*nz, 1, (fftw_complex *) work, 1, 0);
#endif
}

void PmeXPencil::pme_kspace() {

  evir = 0.;

#ifdef FFTCHECK
  return;
#endif

  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;

  int numGrids = 1;
  for ( int g=0; g<numGrids; ++g ) {
    evir[0] = myKSpace->compute_energy(data+0*g,
		lattice, ewaldcof, &(evir[1]));
  }

}

void PmeXPencil::backward_fft() {
#ifdef NAMD_FFTW
  fftw(backward_plan, ny*nz,
	((fftw_complex *) data), ny*nz, 1, (fftw_complex *) work, 1, 0);
#endif
}

void PmeXPencil::send_untrans(int dest) {
  int xBlocks = initdata.xBlocks;
  int block1 = initdata.grid.block1;
  int K1 = initdata.grid.K1;
  for ( int isend=0; isend<xBlocks; ++isend ) {
    int ib = send_order[isend];
    int nx = block1;
    if ( (ib+1)*block1 > K1 ) nx = K1 - ib*block1;
    PmeUntransMsg *msg = new (nx*ny*nz*2,(ib==0?1:0),PRIORITY_SIZE) PmeUntransMsg;
    if ( ib == 0 ) msg->evir[0] = evir;
    msg->sourceNode = thisIndex.y;
    msg->ny = ny;
    float *md = msg->qgrid;
    for ( int i=ib*block1; i<(ib*block1+nx); ++i ) {
     float *d = data + i*ny*nz*2;
     for ( int j=0; j<ny; ++j, d += nz*2 ) {
      for ( int k=0; k<nz; ++k ) {
        *(md++) = d[2*k];
        *(md++) = d[2*k+1];
      }
     }
    }
    SET_PRIORITY(msg,sequence,PME_UNTRANS_PRIORITY)
    initdata.yPencil(ib,0,thisIndex.z).recvUntrans(msg);
  }
}

void PmeYPencil::recv_untrans(const PmeUntransMsg *msg) {
  if ( imsg == 0 ) evir = 0.;
  if ( thisIndex.x == 0 ) evir += msg->evir[0];
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  int jb = msg->sourceNode;
  int ny = msg->ny;
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
#if CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif   
    for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
#ifdef ZEROCHECK
	if ( (*md) == 0. ) CkPrintf("0 in XY at %d %d %d %d %d %d %d %d %d\n",
				    thisIndex.x, jb, thisIndex.z, i, j, k, nx, ny, nz);
#endif
	d[2*(j*nz+k)] = *(md++);
	d[2*(j*nz+k)+1] = *(md++);
      }
    }
  }
}

void PmeYPencil::backward_fft() {
#ifdef NAMD_FFTW
  for ( int i=0; i<nx; ++i ) {
#if CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif

    fftw(backward_plan, nz,
	((fftw_complex *) data) + i * nz * initdata.grid.K2,
	nz, 1, (fftw_complex *) work, 1, 0);
  }
#endif
}

void PmeYPencil::send_untrans(int dest) {
  int yBlocks = initdata.yBlocks;
  int block2 = initdata.grid.block2;
  int K2 = initdata.grid.K2;
  for ( int isend=0; isend<yBlocks; ++isend ) {
    int jb = send_order[isend];
    int ny = block2;
    if ( (jb+1)*block2 > K2 ) ny = K2 - jb*block2;
    PmeUntransMsg *msg = new (nx*ny*nz*2,(jb==0?1:0),PRIORITY_SIZE) PmeUntransMsg;
    if ( jb == 0 ) msg->evir[0] = evir;
    msg->sourceNode = thisIndex.z;
    msg->ny = nz;
    float *md = msg->qgrid;
    const float *d = data;
    for ( int i=0; i<nx; ++i, d += K2*nz*2 ) {
     for ( int j=jb*block2; j<(jb*block2+ny); ++j ) {
      for ( int k=0; k<nz; ++k ) {
        *(md++) = d[2*(j*nz+k)];
        *(md++) = d[2*(j*nz+k)+1];
      }
     }
    }
    SET_PRIORITY(msg,sequence,PME_UNTRANS2_PRIORITY)
    initdata.zPencil(thisIndex.x,jb,0).recvUntrans(msg);
  }
}

void PmeZPencil::recv_untrans(const PmeUntransMsg *msg) {
  if ( imsg == 0 ) evir = 0.;
  if ( thisIndex.y == 0 ) evir += msg->evir[0];
  int block3 = initdata.grid.block3;
  int dim3 = initdata.grid.dim3;
  int kb = msg->sourceNode;
  int nz = msg->ny;
  const float *md = msg->qgrid;
  float *d = data;
  for ( int i=0; i<nx; ++i ) {
#if CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif   
    for ( int j=0; j<ny; ++j, d += dim3 ) {
      for ( int k=kb*block3; k<(kb*block3+nz); ++k ) {
#ifdef ZEROCHECK
	if ( (*md) == 0. ) CkPrintf("0 in YZ at %d %d %d %d %d %d %d %d %d\n",
				    thisIndex.x, thisIndex.y, kb, i, j, k, nx, ny, nz);
#endif
	d[2*k] = *(md++);
	d[2*k+1] = *(md++);
      }
    }
  }
}

void PmeZPencil::backward_fft() {
#ifdef NAMD_FFTW
  rfftwnd_complex_to_real(backward_plan, nx*ny,
	    (fftw_complex *) data, 1, initdata.grid.dim3/2, work, 1, 0);
#endif
  
#if CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

#ifdef FFTCHECK
  int dim3 = initdata.grid.dim3;
  int K1 = initdata.grid.K1;
  int K2 = initdata.grid.K2;
  int K3 = initdata.grid.K3;
  float scale = 1. / (1. * K1 * K2 * K3);
  float maxerr = 0.;
  float maxstd = 0.;
  int mi, mj, mk;  mi = mj = mk = -1;
  float std_base = 100. * (thisIndex.x+1.) + 10. * (thisIndex.y+1.);
  const float *d = data;
  for ( int i=0; i<nx; ++i ) {
   for ( int j=0; j<ny; ++j, d += dim3 ) {
    for ( int k=0; k<K3; ++k ) {
      float std = 10. * (10. * (10. * std_base + i) + j) + k;
      float err = scale * d[k] - std;
      if ( fabsf(err) > fabsf(maxerr) ) {
        maxerr = err;
        maxstd = std;
        mi = i;  mj = j;  mk = k;
      }
    }
   }
  }
  CkPrintf("pencil %d %d max error %f at %d %d %d (should be %f)\n",
		thisIndex.x, thisIndex.y, maxerr, mi, mj, mk, maxstd);
#endif
}

void PmeZPencil::send_ungrid(PmeGridMsg *msg) {
  if ( imsg == 0 ) msg->evir[0] = evir; else msg->evir[0] = 0.;

  int pe = msg->sourceNode;
  msg->sourceNode = thisIndex.x * initdata.yBlocks + thisIndex.y;
  int dim3 = initdata.grid.dim3;
  int zlistlen = msg->zlistlen;
  int *zlist = msg->zlist;
  char *fmsg = msg->fgrid;
  float *qmsg = msg->qgrid;
  float *d = data;
  int numGrids = 1;  // pencil FFT doesn't support multiple grids
  for ( int g=0; g<numGrids; ++g ) {
#if CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif    
    for ( int i=0; i<nx; ++i ) {
      for ( int j=0; j<ny; ++j, d += dim3 ) {
	if( *(fmsg++) ) {
	  for ( int k=0; k<zlistlen; ++k ) {
	    *(qmsg++) = d[zlist[k]];
	  }
	}
      }
    }
  }

  SET_PRIORITY(msg,sequence,PME_UNGRID_PRIORITY)
  initdata.pmeProxy[pe].recvUngrid(msg);
}


#include "ComputePmeMgr.def.h"

