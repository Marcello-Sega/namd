/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "PDB.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeMsm.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "Debug.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "varsizemsg.h"
#include <stdio.h>
#include "MsmMap.h"


// use the decomposition of grid cutoff to create more work units
#define MSM_GRID_CUTOFF_DECOMP
//#undef MSM_GRID_CUTOFF_DECOMP

// skip over pairs of blocks that do not actually interact
#define MSM_SKIP_TOO_DISTANT_BLOCKS
//#undef MSM_SKIP_TOO_DISTANT_BLOCKS

// node aware mapping of chare arrays
#define MSM_NODE_MAPPING
#undef MSM_NODE_MAPPING

// report timings for compute routines
#define MSM_TIMING
#undef MSM_TIMING

// report profiling for compute routines
#define MSM_PROFILING
#undef MSM_PROFILING

// use fixed size grid message
#define MSM_FIXED_SIZE_GRID_MSG
//#undef MSM_FIXED_SIZE_GRID_MSG

#ifdef MSM_FIXED_SIZE_GRID_MSG
#define MSM_MAX_BLOCK_SIZE 8
#define MSM_MAX_BLOCK_VOLUME 512
#endif

//
// This is the main message that gets passed between compute chares.
// It is used to bundle blocks of charge (sendUp and send to MsmGridCutoff) 
// and blocks of potential (sendAcross, sendDown, and sendPatch).  
//
// Higher priority has a numerically lower value.  
//
// The priorities are set as follows:
//
//   sendUp priority = level+1
//
//   (send to MsmGridCutoff) and sendAcross priority
//       = nlevels + 2*(nlevels - level) - 1
//
//   sendDown and sendPatch priority
//       = nlevels + 2*(nlevels - level)
//
// This puts the priority on going up the hierarchy before going across 
// and puts the priority on finishing the top levels and down before 
// finishing the lower levels.
//
class GridFloatMsg : public CMessage_GridFloatMsg {
  public:
    int idnum;
    int nlower_i;
    int nlower_j;
    int nlower_k;
    int nextent_i;
    int nextent_j;
    int nextent_k;
    int nelems;
#ifdef MSM_FIXED_SIZE_GRID_MSG
    Float gdata[MSM_MAX_BLOCK_VOLUME];
#else
    Float *gdata;
#endif
    // put a grid into an allocated message to be sent
    void put(const msm::Grid<Float>& g, int id) {
      idnum = id;
      nlower_i = g.lower().i;
      nlower_j = g.lower().j;
      nlower_k = g.lower().k;
      nextent_i = g.extent().i;
      nextent_j = g.extent().j;
      nextent_k = g.extent().k;
      nelems = g.data().len();
      const Float *p = g.data().buffer();
      for (int i = 0;  i < nelems;  i++) {
        gdata[i] = p[i];
      }
    }
    // get the grid from a received message
    void get(msm::Grid<Float>& g, int& id) {
      id = idnum;
      g.set(nlower_i, nextent_i, nlower_j, nextent_j,
          nlower_k, nextent_k);
      ASSERT(g.data().len() == nelems);
      Float *p = g.data().buffer();
      for (int i = 0;  i < nelems;  i++) {
        p[i] = gdata[i];
      }
    }
};


class MsmBlockProxyMsg : public CMessage_MsmBlockProxyMsg {
  public:
    int len;
    char *msmBlockProxyData;
    // put an array into an allocated message to be sent
    void put(const msm::Array<CProxy_MsmBlock>& a) {
      len = a.len();
      memcpy(msmBlockProxyData, a.buffer(), len*sizeof(CProxy_MsmBlock));
#if 0
      for (int i = 0;  i < len;  i++) {
        msmBlock[i] = a[i];
      }
#endif
    }
    // get the array from a received message
    void get(msm::Array<CProxy_MsmBlock>& a) {
      a.resize(len);
      memcpy(a.buffer(), msmBlockProxyData, len*sizeof(CProxy_MsmBlock));
#if 0
      for (int i = 0;  i < len;  i++) {
        a[i] = msmBlock[i];
      }
#endif
    }
};


class MsmGridCutoffProxyMsg : public CMessage_MsmGridCutoffProxyMsg {
  public:
    char *msmGridCutoffProxyData;
    // put proxy into an allocated message to be sent
    void put(const CProxy_MsmGridCutoff *p) {
      memcpy(msmGridCutoffProxyData, p, sizeof(CProxy_MsmGridCutoff));
    }
    // get the proxy from a received message
    void get(CProxy_MsmGridCutoff *p) {
      memcpy(p, msmGridCutoffProxyData, sizeof(CProxy_MsmGridCutoff));
    }
};


class MsmGridCutoffInitMsg : public CMessage_MsmGridCutoffInitMsg {
  public:
    msm::BlockIndex qhBlockIndex;  // charge block index
    msm::BlockSend ehBlockSend;    // potential block sending address
    MsmGridCutoffInitMsg(const msm::BlockIndex& i, const msm::BlockSend& b)
      : qhBlockIndex(i), ehBlockSend(b) { }
};


class MsmTimer : public CBase_MsmTimer {
  public:
    enum { ANTERP=0, INTERP, RESTRICT, PROLONGATE, GRIDCUTOFF, COMM, MAX };

    MsmTimer() {
      for (int i = 0;  i < MAX;  i++)  timing[i] = 0;
    }
    void done(double tm[], int n) {
      for (int i = 0;  i < MAX;  i++)  timing[i] = tm[i];
      print();
    }
    void print() {
      CkPrintf("MSM timings:\n");
      CkPrintf("   anterpolation   %8.6f sec\n", timing[ANTERP]);
      CkPrintf("   interpolation   %8.6f sec\n", timing[INTERP]);
      CkPrintf("   restriction     %8.6f sec\n", timing[RESTRICT]);
      CkPrintf("   prolongation    %8.6f sec\n", timing[PROLONGATE]);
      CkPrintf("   grid cutoff     %8.6f sec\n", timing[GRIDCUTOFF]);
      CkPrintf("   communication   %8.6f sec\n", timing[COMM]);
    }

    double timing[MAX];
};


class MsmProfiler : public CBase_MsmProfiler {
  public:
    enum { MAX = MSM_MAX_BLOCK_SIZE+1 };

    MsmProfiler() {
      for (int i = 0;  i < MAX;  i++)  xloopcnt[i] = 0;
    }
    void done(int lc[], int n) {
      for (int i = 0;  i < MAX;  i++)  xloopcnt[i] = lc[i];
      print();
    }
    void print() {
      int sum = 0;
      for (int i = 0;  i < MAX;  i++)  sum += xloopcnt[i];
      CkPrintf("MSM profiling:\n");
      CkPrintf("   total executions of inner loop:   %d\n", sum);
      for (int i = 0;  i < MAX;  i++) {
        CkPrintf("   executing %d times:   %d  (%5.2f%%)\n",
            i, xloopcnt[i], 100*double(xloopcnt[i])/sum);
      }
    }

    int xloopcnt[MAX];
};


//////////////////////////////////////////////////////////////////////////////
//
//  ComputeMsmMgr
//  chare group containing MSM parameters and constants;
//  one chare object per PE
//

class ComputeMsmMgr : public BOCclass {
  friend struct msm::PatchData;
  friend class MsmBlock;
  friend class MsmGridCutoff;

public:
  ComputeMsmMgr();                    // entry
  ~ComputeMsmMgr();

  void initialize(MsmInitMsg *);      // entry with message
  void recvMsmBlockProxy(MsmBlockProxyMsg *);  // entry with message
  void recvMsmGridCutoffProxy(MsmGridCutoffProxyMsg *);  // entry with message

  void update(CkQdMsg *);             // entry with message

  void compute(msm::Array<int>& patchIDList);
                                      // called by local ComputeMsm object

  void addPotential(GridFloatMsg *);  // entry with message
  void doneCompute();  // called by each local patch

#ifdef MSM_TIMING
  void initTiming() {
    for (int i = 0;  i < MsmTimer::MAX;  i++)  msmTiming[i] = 0;
    cntTiming = 0;
  }
  // every local object being timed should call this during initialization
  void addTiming() {
    numTiming++;
  }
  void doneTiming() {
    if (++cntTiming >= numTiming) {
      CkCallback cb(CkReductionTarget(MsmTimer, done), msmTimer);
      contribute(MsmTimer::MAX*sizeof(double), msmTiming,
          CkReduction::sum_double, cb);
      initTiming();
    }
  }
#endif

#ifdef MSM_PROFILING
  void initProfiling() {
    for (int i = 0;  i < MsmProfiler::MAX;  i++)  xLoopCnt[i] = 0;
    cntProfiling = 0;
  }
  // every local object being profiled should call this during initialization
  void addProfiling() {
    numProfiling++;
  }
  void doneProfiling() {
    if (++cntProfiling >= numProfiling) {
      CkCallback cb(CkReductionTarget(MsmProfiler, done), msmProfiler);
      contribute(MsmProfiler::MAX*sizeof(int), xLoopCnt,
          CkReduction::sum_int, cb);
      initProfiling();  // reset accumulators for next visit
    }
  }
#endif

  void setCompute(ComputeMsm *c) { msmCompute = c;  c->setMgr(this); } // local

  msm::PatchPtrArray& patchPtrArray() { return patchPtr; }

  msm::Map& mapData() { return map; }

  int numLevels() const { return nlevels; }

  // sign(n) = -1 if n < 0,  0 if n == 0,  or  1 if n > 0
  static inline int sign(int n) {
    return (n < 0 ? -1 : (n > 0 ? 1 : 0));
  }
private:
  void setup_hgrid_1d(BigReal len, BigReal& hh, int& nn,
      int& ia, int& ib, int isperiodic);
  void setup_periodic_blocksize(int& bsize, int n);

  CProxy_ComputeMsmMgr msmProxy;
  ComputeMsm *msmCompute;

  msm::Array<CProxy_MsmBlock> msmBlock;

  CProxy_MsmGridCutoff msmGridCutoff;
  int numGridCutoff;  // length of msmGridCutoff chare array

  msm::Map map;

  // find patch by patchID
  // array is length number of patches, initialized to NULL
  // allocate PatchData for only those patches on this PE
  msm::PatchPtrArray patchPtr;

  // allocate subgrid used for receiving message data in addPotential()
  // and sending on to PatchData::addPotential()
  msm::Grid<Float> subgrid;

#ifdef MSM_NODE_MAPPING
  msm::Array<int> blockAssign;
  msm::Array<int> gcutAssign;
  msm::Array<int> nodecnt;
  int blockFlatIndex(int level, int i, int j, int k) {
    int n = 0;
    for (int l = 0;  l < level;  l++) {
      n += map.blockLevel[l].nn();
    }
    return (n + map.blockLevel[level].flatindex(i,j,k));
  }
#endif

#ifdef MSM_TIMING
  CProxy_MsmTimer msmTimer;
  double msmTiming[MsmTimer::MAX];
  int numTiming;  // total number of objects being timed
  int cntTiming;  // count the objects as they provide timing results
  CkCallback *cbTiming;
#endif

#ifdef MSM_PROFILING
  CProxy_MsmProfiler msmProfiler;
  int xLoopCnt[MsmProfiler::MAX];
  int numProfiling;  // total number of objects being profiled
  int cntProfiling;  // count the objects as they provide profiling results
  CkCallback *cbProfiling;
#endif

  Vector c, u, v, w;    // rescaled center and lattice vectors
  Vector ru, rv, rw;    // row vectors to transform to unit space
  int ispu, ispv, ispw; // is periodic along u, v, w?

  Lattice lattice;      // keep local copy of lattice
  ScaledPosition smin;  // keep min values for non-periodic dimensions
  ScaledPosition smax;  // keep max values for non-periodic dimensions
  BigReal gridspacing;  // preferred grid spacing
  BigReal padding;      // padding for non-periodic boundaries
  Vector h;             // finest level grid spacings
  BigReal a;            // cutoff distance
  int nhx, nhy, nhz;    // number of h spacings that cover cell
  int approx;           // ID for approximation
  int split;            // ID for splitting
  int nlevels;          // number of grid levels
  int dispersion;       // calculating dispersion forces?
  BigReal gzero;        // self energy factor from splitting

  Vector sglower;       // lower corner of grid in scaled space
                        // corresponds to index (0,0,0)

  BigReal shx, shy, shz;  // grid spacings in scaled space
  BigReal shx_1, shy_1, shz_1;
  Vector sx_shx;          // row vector to transform interpolated force x
  Vector sy_shy;          // row vector to transform interpolated force y
  Vector sz_shz;          // row vector to transform interpolated force z
  Float srx_x, srx_y, srx_z;  // float version of sx_shx
  Float sry_x, sry_y, sry_z;  // float version of sy_shy
  Float srz_x, srz_y, srz_z;  // float version of sz_shz

  int s_edge;
  int omega;

  enum Approx { CUBIC=0, QUINTIC, QUINTIC2,
    SEPTIC, SEPTIC3, NONIC, NONIC4, NUM_APPROX };

  enum Split { TAYLOR2=0, TAYLOR3, TAYLOR4,
    TAYLOR5, TAYLOR6, TAYLOR7, TAYLOR8,
    TAYLOR2_DISP, TAYLOR3_DISP, TAYLOR4_DISP, TAYLOR5_DISP,
    TAYLOR6_DISP, TAYLOR7_DISP, TAYLOR8_DISP, NUM_SPLIT };

  enum {
    // Approximation formulas with up to degree 9 polynomials.
    MAX_POLY_DEGREE = 9,

    // Max stencil length for polynomial approximation.
    MAX_NSTENCIL_SIZE = (2*MAX_POLY_DEGREE + 1),

    // Max stencil length when skipping zeros
    // (almost half entries are zero for interpolating polynomials).
    MAX_NSTENCIL_SKIP_ZERO = (MAX_POLY_DEGREE + 2)
  };

  // Degree of polynomial basis function Phi.
  static const int PolyDegree[NUM_APPROX];

  // The stencil array lengths below.
  static const int Nstencil[NUM_APPROX];

  // Index offsets from the stencil-centered grid element, to get
  // to the correct contributing grid element.
  static const int IndexOffset[NUM_APPROX][MAX_NSTENCIL_SKIP_ZERO];

  // The grid transfer stencils for the non-factored restriction and
  // prolongation procedures.
  static const Float PhiStencil[NUM_APPROX][MAX_NSTENCIL_SKIP_ZERO];

  // Calculate the smoothing function and its derivative:
  // g(R) and (d/dR)g(R), where R=r/a.
  static int splitting(BigReal& g, BigReal& dg, BigReal r_a, int _split) {
    BigReal s = r_a * r_a;  // s = (r/a)^2, assuming 0 <= s <= 1
    switch (_split) {
      case TAYLOR2:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4));
        break;
      case TAYLOR3:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16)));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16)));
        break;
      case TAYLOR4:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                + (s-1)*(35./128))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32))));
        break;
      case TAYLOR5:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                + (s-1)*(35./128 + (s-1)*(-63./256)))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32 + (s-1)*(-315./256)))));
        break;
      case TAYLOR6:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                + (s-1)*(35./128 + (s-1)*(-63./256
                    + (s-1)*(231./1024))))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32 + (s-1)*(-315./256
                    + (s-1)*(693./512))))));
        break;
      case TAYLOR7:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
            + (s-1)*(35./128 + (s-1)*(-63./256
                + (s-1)*(231./1024 + (s-1)*(-429./2048)))))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32 + (s-1)*(-315./256
                    + (s-1)*(693./512 + (s-1)*(-3003./2048)))))));
        break;
      case TAYLOR8:
        g = 1 + (s-1)*(-1./2 + (s-1)*(3./8 + (s-1)*(-5./16
                + (s-1)*(35./128 + (s-1)*(-63./256
                    + (s-1)*(231./1024 + (s-1)*(-429./2048
                        + (s-1)*(6435./32768))))))));
        dg = (2*r_a)*(-1./2 + (s-1)*(3./4 + (s-1)*(-15./16
                + (s-1)*(35./32 + (s-1)*(-315./256
                    + (s-1)*(693./512 + (s-1)*(-3003./2048
                        + (s-1)*(6435./4096))))))));
        break;
      case TAYLOR2_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6));
        dg = (2*r_a)*(-3 + (s-1)*(12));
        break;
      case TAYLOR3_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10)));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30)));
        break;
      case TAYLOR4_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10 + (s-1)*(15))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30 + (s-1)*(60))));
        break;
      case TAYLOR5_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10
                + (s-1)*(15 + (s-1)*(-21)))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30
                + (s-1)*(60 + (s-1)*(-105)))));
        break;
      case TAYLOR6_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10
                + (s-1)*(15 + (s-1)*(-21 + (s-1)*(28))))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30
                + (s-1)*(60 + (s-1)*(-105 + (s-1)*(168))))));
        break;
      case TAYLOR7_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10
                + (s-1)*(15 + (s-1)*(-21 + (s-1)*(28
                      + (s-1)*(-36)))))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30
                + (s-1)*(60 + (s-1)*(-105 + (s-1)*(168
                      + (s-1)*(-252)))))));
        break;
      case TAYLOR8_DISP:
        g = 1 + (s-1)*(-3 + (s-1)*(6 + (s-1)*(-10
                + (s-1)*(15 + (s-1)*(-21 + (s-1)*(28
                      + (s-1)*(-36 + (s-1)*(45))))))));
        dg = (2*r_a)*(-3 + (s-1)*(12 + (s-1)*(-30
                + (s-1)*(60 + (s-1)*(-105 + (s-1)*(168
                      + (s-1)*(-252 + (s-1)*(360))))))));
        break;
      default:
        return -1;
    }
    return 0;
  }

  void stencil_1d(Float phi[], Float t) {
    switch (approx) {
      case CUBIC:
        phi[0] = 0.5f * (1 - t) * (2 - t) * (2 - t);
        t--;
        phi[1] = (1 - t) * (1 + t - 1.5f * t * t);
        t--;
        phi[2] = (1 + t) * (1 - t - 1.5f * t * t);
        t--;
        phi[3] = 0.5f * (1 + t) * (2 + t) * (2 + t);
        break;
      case QUINTIC:
        phi[0] = (1.f/24) * (1-t) * (2-t) * (3-t) * (3-t) * (4-t);
        t--;
        phi[1] = (1-t) * (2-t) * (3-t) * ((1.f/6)
            + t * (0.375f - (5.f/24)*t));
        t--;
        phi[2] = (1-t*t) * (2-t) * (0.5f + t * (0.25f - (5.f/12)*t));
        t--;
        phi[3] = (1-t*t) * (2+t) * (0.5f - t * (0.25f + (5.f/12)*t));
        t--;
        phi[4] = (1+t) * (2+t) * (3+t) * ((1.f/6)
            - t * (0.375f + (5.f/24)*t));
        t--;
        phi[5] = (1.f/24) * (1+t) * (2+t) * (3+t) * (3+t) * (4+t);
        break;
      case QUINTIC2:
        phi[0] = (1.f/24) * (3-t) * (3-t) * (3-t) * (t-2) * (5*t-8);
        t--;
        phi[1] = (-1.f/24) * (2-t) * (t-1) * (-48+t*(153+t*(-114+t*25)));
        t--;
        phi[2] = (1.f/12) * (1-t) * (12+t*(12+t*(-3+t*(-38+t*25))));
        t--;
        phi[3] = (1.f/12) * (1+t) * (12+t*(-12+t*(-3+t*(38+t*25))));
        t--;
        phi[4] = (-1.f/24) * (2+t) * (t+1) * (48+t*(153+t*(114+t*25)));
        t--;
        phi[5] = (1.f/24) * (3+t) * (3+t) * (3+t) * (t+2) * (5*t+8);
        break;
      case SEPTIC:
        phi[0] = (-1.f/720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-4)*(t-5)*(t-6);
        t--;
        phi[1] = (1.f/720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)*(-6+t*(-20+7*t));
        t--;
        phi[2] = (-1.f/240)*(t*t-1)*(t-2)*(t-3)*(t-4)*(-10+t*(-12+7*t));
        t--;
        phi[3] = (1.f/144)*(t*t-1)*(t*t-4)*(t-3)*(-12+t*(-4+7*t));
        t--;
        phi[4] = (-1.f/144)*(t*t-1)*(t*t-4)*(t+3)*(-12+t*(4+7*t));
        t--;
        phi[5] = (1.f/240)*(t*t-1)*(t+2)*(t+3)*(t+4)*(-10+t*(12+7*t));
        t--;
        phi[6] = (-1.f/720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(-6+t*(20+7*t));
        t--;
        phi[7] = (1.f/720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+4)*(t+5)*(t+6);
        break;
      case SEPTIC3:
        phi[0] = (3632.f/5) + t*((-7456.f/5) + t*((58786.f/45) + t*(-633
                + t*((26383.f/144) + t*((-22807.f/720) + t*((727.f/240)
                      + t*(-89.f/720)))))));
        t--;
        phi[1] = -440 + t*((25949.f/20) + t*((-117131.f/72) + t*((2247.f/2)
                + t*((-66437.f/144) + t*((81109.f/720) + t*((-727.f/48)
                      + t*(623.f/720)))))));
        t--;
        phi[2] = (138.f/5) + t*((-8617.f/60) + t*((12873.f/40) + t*((-791.f/2)
                + t*((4557.f/16) + t*((-9583.f/80) + t*((2181.f/80)
                      + t*(-623.f/240)))))));
        t--;
        phi[3] = 1 + t*t*((-49.f/36) + t*t*((-959.f/144) + t*((2569.f/144)
                + t*((-727.f/48) + t*(623.f/144)))));
        t--;
        phi[4] = 1 + t*t*((-49.f/36) + t*t*((-959.f/144) + t*((-2569.f/144)
                + t*((-727.f/48) + t*(-623.f/144)))));
        t--;
        phi[5] = (138.f/5) + t*((8617.f/60) + t*((12873.f/40) + t*((791.f/2)
                + t*((4557.f/16) + t*((9583.f/80) + t*((2181.f/80)
                      + t*(623.f/240)))))));
        t--;
        phi[6] = -440 + t*((-25949.f/20) + t*((-117131.f/72) + t*((-2247.f/2)
                + t*((-66437.f/144) + t*((-81109.f/720) + t*((-727.f/48)
                      + t*(-623.f/720)))))));
        t--;
        phi[7] = (3632.f/5) + t*((7456.f/5) + t*((58786.f/45) + t*(633
                + t*((26383.f/144) + t*((22807.f/720) + t*((727.f/240)
                      + t*(89.f/720)))))));
        break;
      case NONIC:
        phi[0] = (-1.f/40320)*(t-8)*(t-7)*(t-6)*(t-5)*(t-5)*(t-4)*(t-3)*
          (t-2)*(t-1);
        t--;
        phi[1] = (1.f/40320)*(t-7)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*
          (-8+t*(-35+9*t));
        t--;
        phi[2] = (-1.f/10080)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*
          (-14+t*(-25+9*t));
        t--;
        phi[3] = (1.f/1440)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*
          (-6+t*(-5+3*t));
        t--;
        phi[4] = (-1.f/2880)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*
          (-20+t*(-5+9*t));
        t--;
        phi[5] = (1.f/2880)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*
          (-20+t*(5+9*t));
        t--;
        phi[6] = (-1.f/1440)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*
          (-6+t*(5+3*t));
        t--;
        phi[7] = (1.f/10080)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*
          (-14+t*(25+9*t));
        t--;
        phi[8] = (-1.f/40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*(t+7)*
          (-8+t*(35+9*t));
        t--;
        phi[9] = (1.f/40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+5)*(t+6)*
          (t+7)*(t+8);
        break;
      case NONIC4:
        phi[0] = 439375.f/7+t*(-64188125.f/504+t*(231125375.f/2016
              +t*(-17306975.f/288+t*(7761805.f/384+t*(-2895587.f/640
                    +t*(129391.f/192+t*(-259715.f/4032+t*(28909.f/8064
                          +t*(-3569.f/40320)))))))));
        t--;
        phi[1] = -56375+t*(8314091.f/56+t*(-49901303.f/288+t*(3763529.f/32
                +t*(-19648027.f/384+t*(9469163.f/640+t*(-545977.f/192
                      +t*(156927.f/448+t*(-28909.f/1152
                          +t*(3569.f/4480)))))))));
        t--;
        phi[2] = 68776.f/7+t*(-1038011.f/28+t*(31157515.f/504+t*(-956669.f/16
                +t*(3548009.f/96+t*(-2422263.f/160+t*(197255.f/48
                      +t*(-19959.f/28+t*(144545.f/2016
                          +t*(-3569.f/1120)))))))));
        t--;
        phi[3] = -154+t*(12757.f/12+t*(-230123.f/72+t*(264481.f/48
                +t*(-576499.f/96+t*(686147.f/160+t*(-96277.f/48
                      +t*(14221.f/24+t*(-28909.f/288+t*(3569.f/480)))))))));
        t--;
        phi[4] = 1+t*t*(-205.f/144+t*t*(91.f/192+t*(-6181.f/320
                +t*(6337.f/96+t*(-2745.f/32+t*(28909.f/576
                      +t*(-3569.f/320)))))));
        t--;
        phi[5] = 1+t*t*(-205.f/144+t*t*(91.f/192+t*(6181.f/320
                +t*(6337.f/96+t*(2745.f/32+t*(28909.f/576
                      +t*(3569.f/320)))))));
        t--;
        phi[6] = -154+t*(-12757.f/12+t*(-230123.f/72+t*(-264481.f/48
                +t*(-576499.f/96+t*(-686147.f/160+t*(-96277.f/48
                      +t*(-14221.f/24+t*(-28909.f/288+t*(-3569.f/480)))))))));
        t--;
        phi[7] = 68776.f/7+t*(1038011.f/28+t*(31157515.f/504+t*(956669.f/16
                +t*(3548009.f/96+t*(2422263.f/160+t*(197255.f/48
                      +t*(19959.f/28+t*(144545.f/2016+t*(3569.f/1120)))))))));
        t--;
        phi[8] = -56375+t*(-8314091.f/56+t*(-49901303.f/288+t*(-3763529.f/32
                +t*(-19648027.f/384+t*(-9469163.f/640+t*(-545977.f/192
                      +t*(-156927.f/448+t*(-28909.f/1152
                          +t*(-3569.f/4480)))))))));
        t--;
        phi[9] = 439375.f/7+t*(64188125.f/504+t*(231125375.f/2016
              +t*(17306975.f/288+t*(7761805.f/384+t*(2895587.f/640
                    +t*(129391.f/192+t*(259715.f/4032+t*(28909.f/8064
                          +t*(3569.f/40320)))))))));
        break;
      default:
        NAMD_die("Unknown MSM approximation.");
    } // switch
  }

  void d_stencil_1d(Float dphi[], Float phi[], Float t) {
    switch (approx) {
      case CUBIC:
        phi[0] = 0.5f * (1 - t) * (2 - t) * (2 - t);
        dphi[0] = (1.5f * t - 2) * (2 - t);
        t--;
        phi[1] = (1 - t) * (1 + t - 1.5f * t * t);
        dphi[1] = (-5 + 4.5f * t) * t;
        t--;
        phi[2] = (1 + t) * (1 - t - 1.5f * t * t);
        dphi[2] = (-5 - 4.5f * t) * t;
        t--;
        phi[3] = 0.5f * (1 + t) * (2 + t) * (2 + t);
        dphi[3] = (1.5f * t + 2) * (2 + t);
        break;
      case QUINTIC:
        phi[0] = (1.f/24) * (1-t) * (2-t) * (3-t) * (3-t) * (4-t);
        dphi[0] = ((-1.f/24) * ((3-t) * (3-t) * (14 + t * (-14 + 3*t))
              + 2 * (1-t) * (2-t) * (3-t) * (4-t)));
        t--;
        phi[1] = (1-t) * (2-t) * (3-t) * ((1.f/6)
            + t * (0.375f - (5.f/24)*t));
        dphi[1] = (-((1.f/6) + t * (0.375f - (5.f/24)*t)) *
            (11 + t * (-12 + 3*t)) + (1-t) * (2-t) * (3-t) *
            (0.375f - (5.f/12)*t));
        t--;
        phi[2] = (1-t*t) * (2-t) * (0.5f + t * (0.25f - (5.f/12)*t));
        dphi[2] = (-(0.5f + t * (0.25f - (5.f/12)*t)) * (1 + t * (4 - 3*t))
            + (1-t*t) * (2-t) * (0.25f - (5.f/6)*t));
        t--;
        phi[3] = (1-t*t) * (2+t) * (0.5f - t * (0.25f + (5.f/12)*t));
        dphi[3] = ((0.5f + t * (-0.25f - (5.f/12)*t)) * (1 + t * (-4 - 3*t))
            - (1-t*t) * (2+t) * (0.25f + (5.f/6)*t));
        t--;
        phi[4] = (1+t) * (2+t) * (3+t) * ((1.f/6)
            - t * (0.375f + (5.f/24)*t));
        dphi[4] = (((1.f/6) + t * (-0.375f - (5.f/24)*t)) *
            (11 + t * (12 + 3*t)) - (1+t) * (2+t) * (3+t) *
            (0.375f + (5.f/12)*t));
        t--;
        phi[5] = (1.f/24) * (1+t) * (2+t) * (3+t) * (3+t) * (4+t);
        dphi[5] = ((1.f/24) * ((3+t) * (3+t) * (14 + t * (14 + 3*t))
              + 2 * (1+t) * (2+t) * (3+t) * (4+t)));
        break;
      case QUINTIC2:
        phi[0] = (1.f/24) * (3-t) * (3-t) * (3-t) * (t-2) * (5*t-8);
        dphi[0] = ((1.f/24) * (3-t) * (3-t) * ((3-t)*(5*t-8)
              - 3*(t-2)*(5*t-8) + 5*(t-2)*(3-t)));
        t--;
        phi[1] = (-1.f/24) * (2-t) * (t-1) * (-48+t*(153+t*(-114+t*25)));
        dphi[1] = ((-1.f/24) * ((2-t)*(-48+t*(153+t*(-114+t*25)))
              - (t-1)* (-48+t*(153+t*(-114+t*25)))
              + (2-t)*(t-1)*(153+t*(-228+t*75))));
        t--;
        phi[2] = (1.f/12) * (1-t) * (12+t*(12+t*(-3+t*(-38+t*25))));
        dphi[2] = ((1.f/12) * (-(12+t*(12+t*(-3+t*(-38+t*25))))
              + (1-t)*(12+t*(-6+t*(-114+t*100)))));
        t--;
        phi[3] = (1.f/12) * (1+t) * (12+t*(-12+t*(-3+t*(38+t*25))));
        dphi[3] = ((1.f/12) * ((12+t*(-12+t*(-3+t*(38+t*25))))
              + (1+t)*(-12+t*(-6+t*(114+t*100)))));
        t--;
        phi[4] = (-1.f/24) * (2+t) * (t+1) * (48+t*(153+t*(114+t*25)));
        dphi[4] = ((-1.f/24) * ((2+t)*(48+t*(153+t*(114+t*25)))
              + (t+1)* (48+t*(153+t*(114+t*25)))
              + (2+t)*(t+1)*(153+t*(228+t*75))));
        t--;
        phi[5] = (1.f/24) * (3+t) * (3+t) * (3+t) * (t+2) * (5*t+8);
        dphi[5] = ((1.f/24) * (3+t) * (3+t) * ((3+t)*(5*t+8)
              + 3*(t+2)*(5*t+8) + 5*(t+2)*(3+t)));
        break;
      case SEPTIC:
        phi[0] = (-1.f/720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-4)*(t-5)*(t-6);
        dphi[0] = (-1.f/720)*(t-4)*(-1944+t*(3644+t*(-2512+t*(807
                  +t*(-122+t*7)))));
        t--;
        phi[1] = (1.f/720)*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)*(-6+t*(-20+7*t));
        dphi[1] = (1.f/720)*(756+t*(-9940+t*(17724+t*(-12740+t*(4445
                    +t*(-750+t*49))))));
        t--;
        phi[2] = (-1.f/240)*(t*t-1)*(t-2)*(t-3)*(t-4)*(-10+t*(-12+7*t));
        dphi[2] = (-1.f/240)*(-28+t*(1260+t*(-756+t*(-1260+t*(1365
                    +t*(-450+t*49))))));
        t--;
        phi[3] = (1.f/144)*(t*t-1)*(t*t-4)*(t-3)*(-12+t*(-4+7*t));
        dphi[3] = (1.f/144)*t*(-560+t*(84+t*(644+t*(-175
                  +t*(-150+t*49)))));
        t--;
        phi[4] = (-1.f/144)*(t*t-1)*(t*t-4)*(t+3)*(-12+t*(4+7*t));
        dphi[4] = (-1.f/144)*t*(560+t*(84+t*(-644+t*(-175
                  +t*(150+t*49)))));
        t--;
        phi[5] = (1.f/240)*(t*t-1)*(t+2)*(t+3)*(t+4)*(-10+t*(12+7*t));
        dphi[5] = (1.f/240)*(-28+t*(-1260+t*(-756+t*(1260+t*(1365
                    +t*(450+t*49))))));
        t--;
        phi[6] = (-1.f/720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(-6+t*(20+7*t));
        dphi[6] = (-1.f/720)*(756+t*(9940+t*(17724+t*(12740+t*(4445
                    +t*(750+t*49))))));
        t--;
        phi[7] = (1.f/720)*(t+1)*(t+2)*(t+3)*(t+4)*(t+4)*(t+5)*(t+6);
        dphi[7] = (1.f/720)*(t+4)*(1944+t*(3644+t*(2512+t*(807
                  +t*(122+t*7)))));
        break;
      case SEPTIC3:
        phi[0] = (3632.f/5) + t*((-7456.f/5) + t*((58786.f/45) + t*(-633
                + t*((26383.f/144) + t*((-22807.f/720) + t*((727.f/240)
                      + t*(-89.f/720)))))));
        dphi[0] = ((-7456.f/5) + t*((117572.f/45) + t*(-1899
                + t*((26383.f/36) + t*((-22807.f/144) + t*((727.f/40)
                      + t*(-623.f/720)))))));
        t--;
        phi[1] = -440 + t*((25949.f/20) + t*((-117131.f/72) + t*((2247.f/2)
                + t*((-66437.f/144) + t*((81109.f/720) + t*((-727.f/48)
                      + t*(623.f/720)))))));
        dphi[1] = ((25949.f/20) + t*((-117131.f/36) + t*((6741.f/2)
                + t*((-66437.f/36) + t*((81109.f/144) + t*((-727.f/8)
                      + t*(4361.f/720)))))));
        t--;
        phi[2] = (138.f/5) + t*((-8617.f/60) + t*((12873.f/40) + t*((-791.f/2)
                + t*((4557.f/16) + t*((-9583.f/80) + t*((2181.f/80)
                      + t*(-623.f/240)))))));
        dphi[2] = ((-8617.f/60) + t*((12873.f/20) + t*((-2373.f/2)
                + t*((4557.f/4) + t*((-9583.f/16) + t*((6543.f/40)
                      + t*(-4361.f/240)))))));
        t--;
        phi[3] = 1 + t*t*((-49.f/36) + t*t*((-959.f/144) + t*((2569.f/144)
                + t*((-727.f/48) + t*(623.f/144)))));
        dphi[3] = (t*((-49.f/18) + t*t*((-959.f/36) + t*((12845.f/144)
                  + t*((-727.f/8) + t*(4361.f/144))))));
        t--;
        phi[4] = 1 + t*t*((-49.f/36) + t*t*((-959.f/144) + t*((-2569.f/144)
                + t*((-727.f/48) + t*(-623.f/144)))));
        dphi[4] = (t*((-49.f/18) + t*t*((-959.f/36) + t*((-12845.f/144)
                  + t*((-727.f/8) + t*(-4361.f/144))))));
        t--;
        phi[5] = (138.f/5) + t*((8617.f/60) + t*((12873.f/40) + t*((791.f/2)
                + t*((4557.f/16) + t*((9583.f/80) + t*((2181.f/80)
                      + t*(623.f/240)))))));
        dphi[5] = ((8617.f/60) + t*((12873.f/20) + t*((2373.f/2)
                + t*((4557.f/4) + t*((9583.f/16) + t*((6543.f/40)
                      + t*(4361.f/240)))))));
        t--;
        phi[6] = -440 + t*((-25949.f/20) + t*((-117131.f/72) + t*((-2247.f/2)
                + t*((-66437.f/144) + t*((-81109.f/720) + t*((-727.f/48)
                      + t*(-623.f/720)))))));
        dphi[6] = ((-25949.f/20) + t*((-117131.f/36) + t*((-6741.f/2)
                + t*((-66437.f/36) + t*((-81109.f/144) + t*((-727.f/8)
                      + t*(-4361.f/720)))))));
        t--;
        phi[7] = (3632.f/5) + t*((7456.f/5) + t*((58786.f/45) + t*(633
                + t*((26383.f/144) + t*((22807.f/720) + t*((727.f/240)
                      + t*(89.f/720)))))));
        dphi[7] = ((7456.f/5) + t*((117572.f/45) + t*(1899
                + t*((26383.f/36) + t*((22807.f/144) + t*((727.f/40)
                      + t*(623.f/720)))))));
        break;
      case NONIC:
        phi[0] = (-1.f/40320)*(t-8)*(t-7)*(t-6)*(t-5)*(t-5)*(t-4)*(t-3)*
          (t-2)*(t-1);
        dphi[0] = (-1.f/40320)*(t-5)*(-117648+t*(256552+t*(-221416
                +t*(99340+t*(-25261+t*(3667+t*(-283+t*9)))))));
        t--;
        phi[1] = (1.f/40320)*(t-7)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*
          (-8+t*(-35+9*t));
        dphi[1] = (1.f/40320)*(71856+t*(-795368+t*(1569240+t*(-1357692
                  +t*(634725+t*(-172116+t*(27090+t*(-2296
                          +t*81))))))));
        t--;
        phi[2] = (-1.f/10080)*(t-6)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*
          (-14+t*(-25+9*t));
        dphi[2] = (1.f/10080)*(3384+t*(-69080+t*(55026+t*(62580
                  +t*(-99225+t*(51660+t*(-13104+t*(1640
                          +t*(-81)))))))));
        t--;
        phi[3] = (1.f/1440)*(t-5)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*
          (-6+t*(-5+3*t));
        dphi[3] = (1.f/1440)*(72+t*(-6344+t*(2070+t*(7644
                  +t*(-4725+t*(-828+t*(1260+t*(-328+t*27))))))));
        t--;
        phi[4] = (-1.f/2880)*(t-4)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*
          (-20+t*(-5+9*t));
        dphi[4] = (-1.f/2880)*t*(10792+t*(-972+t*(-12516
                +t*(2205+t*(3924+t*(-882+t*(-328+t*81)))))));
        t--;
        phi[5] = (1.f/2880)*(t-3)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*
          (-20+t*(5+9*t));
        dphi[5] = (1.f/2880)*t*(-10792+t*(-972+t*(12516
                +t*(2205+t*(-3924+t*(-882+t*(328+t*81)))))));
        t--;
        phi[6] = (-1.f/1440)*(t-2)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*
          (-6+t*(5+3*t));
        dphi[6] = (1.f/1440)*(-72+t*(-6344+t*(-2070+t*(7644
                  +t*(4725+t*(-828+t*(-1260+t*(-328+t*(-27)))))))));
        t--;
        phi[7] = (1.f/10080)*(t-1)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*
          (-14+t*(25+9*t));
        dphi[7] = (1.f/10080)*(-3384+t*(-69080+t*(-55026+t*(62580
                  +t*(99225+t*(51660+t*(13104+t*(1640
                          +t*81))))))));
        t--;
        phi[8] = (-1.f/40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+6)*(t+7)*
          (-8+t*(35+9*t));
        dphi[8] = (-1.f/40320)*(71856+t*(795368+t*(1569240
                +t*(1357692+t*(634725+t*(172116+t*(27090+t*(2296
                          +t*81))))))));
        t--;
        phi[9] = (1.f/40320)*(t+1)*(t+2)*(t+3)*(t+4)*(t+5)*(t+5)*(t+6)*
          (t+7)*(t+8);
        dphi[9] = (1.f/40320)*(t+5)*(117648+t*(256552+t*(221416
                +t*(99340+t*(25261+t*(3667+t*(283+t*9)))))));
        break;
      case NONIC4:
        phi[0] = 439375.f/7+t*(-64188125.f/504+t*(231125375.f/2016
              +t*(-17306975.f/288+t*(7761805.f/384+t*(-2895587.f/640
                    +t*(129391.f/192+t*(-259715.f/4032+t*(28909.f/8064
                          +t*(-3569.f/40320)))))))));
        dphi[0] = (-64188125.f/504+t*(231125375.f/1008
              +t*(-17306975.f/96+t*(7761805.f/96+t*(-2895587.f/128
                    +t*(129391.f/32+t*(-259715.f/576+t*(28909.f/1008
                          +t*(-3569.f/4480)))))))));
        t--;
        phi[1] = -56375+t*(8314091.f/56+t*(-49901303.f/288+t*(3763529.f/32
                +t*(-19648027.f/384+t*(9469163.f/640+t*(-545977.f/192
                      +t*(156927.f/448+t*(-28909.f/1152
                          +t*(3569.f/4480)))))))));
        dphi[1] = (8314091.f/56+t*(-49901303.f/144+t*(11290587.f/32
                +t*(-19648027.f/96+t*(9469163.f/128+t*(-545977.f/32
                      +t*(156927.f/64+t*(-28909.f/144
                          +t*(32121.f/4480)))))))));
        t--;
        phi[2] = 68776.f/7+t*(-1038011.f/28+t*(31157515.f/504+t*(-956669.f/16
                +t*(3548009.f/96+t*(-2422263.f/160+t*(197255.f/48
                      +t*(-19959.f/28+t*(144545.f/2016
                          +t*(-3569.f/1120)))))))));
        dphi[2] = (-1038011.f/28+t*(31157515.f/252+t*(-2870007.f/16
                +t*(3548009.f/24+t*(-2422263.f/32+t*(197255.f/8
                      +t*(-19959.f/4+t*(144545.f/252
                          +t*(-32121.f/1120)))))))));
        t--;
        phi[3] = -154+t*(12757.f/12+t*(-230123.f/72+t*(264481.f/48
                +t*(-576499.f/96+t*(686147.f/160+t*(-96277.f/48
                      +t*(14221.f/24+t*(-28909.f/288+t*(3569.f/480)))))))));
        dphi[3] = (12757.f/12+t*(-230123.f/36+t*(264481.f/16
                +t*(-576499.f/24+t*(686147.f/32+t*(-96277.f/8
                      +t*(99547.f/24+t*(-28909.f/36
                          +t*(10707.f/160)))))))));
        t--;
        phi[4] = 1+t*t*(-205.f/144+t*t*(91.f/192+t*(-6181.f/320
                +t*(6337.f/96+t*(-2745.f/32+t*(28909.f/576
                      +t*(-3569.f/320)))))));
        dphi[4] = t*(-205.f/72+t*t*(91.f/48+t*(-6181.f/64
                +t*(6337.f/16+t*(-19215.f/32+t*(28909.f/72
                      +t*(-32121.f/320)))))));
        t--;
        phi[5] = 1+t*t*(-205.f/144+t*t*(91.f/192+t*(6181.f/320
                +t*(6337.f/96+t*(2745.f/32+t*(28909.f/576
                      +t*(3569.f/320)))))));
        dphi[5] = t*(-205.f/72+t*t*(91.f/48+t*(6181.f/64
                +t*(6337.f/16+t*(19215.f/32+t*(28909.f/72
                      +t*(32121.f/320)))))));
        t--;
        phi[6] = -154+t*(-12757.f/12+t*(-230123.f/72+t*(-264481.f/48
                +t*(-576499.f/96+t*(-686147.f/160+t*(-96277.f/48
                      +t*(-14221.f/24+t*(-28909.f/288+t*(-3569.f/480)))))))));
        dphi[6] = (-12757.f/12+t*(-230123.f/36+t*(-264481.f/16
                +t*(-576499.f/24+t*(-686147.f/32+t*(-96277.f/8
                      +t*(-99547.f/24+t*(-28909.f/36
                          +t*(-10707.f/160)))))))));
        t--;
        phi[7] = 68776.f/7+t*(1038011.f/28+t*(31157515.f/504+t*(956669.f/16
                +t*(3548009.f/96+t*(2422263.f/160+t*(197255.f/48
                      +t*(19959.f/28+t*(144545.f/2016+t*(3569.f/1120)))))))));
        dphi[7] = (1038011.f/28+t*(31157515.f/252+t*(2870007.f/16
                +t*(3548009.f/24+t*(2422263.f/32+t*(197255.f/8
                      +t*(19959.f/4+t*(144545.f/252
                          +t*(32121.f/1120)))))))));
        t--;
        phi[8] = -56375+t*(-8314091.f/56+t*(-49901303.f/288+t*(-3763529.f/32
                +t*(-19648027.f/384+t*(-9469163.f/640+t*(-545977.f/192
                      +t*(-156927.f/448+t*(-28909.f/1152
                          +t*(-3569.f/4480)))))))));
        dphi[8] = (-8314091.f/56+t*(-49901303.f/144+t*(-11290587.f/32
                +t*(-19648027.f/96+t*(-9469163.f/128+t*(-545977.f/32
                      +t*(-156927.f/64+t*(-28909.f/144
                          +t*(-32121.f/4480)))))))));
        t--;
        phi[9] = 439375.f/7+t*(64188125.f/504+t*(231125375.f/2016
              +t*(17306975.f/288+t*(7761805.f/384+t*(2895587.f/640
                    +t*(129391.f/192+t*(259715.f/4032+t*(28909.f/8064
                          +t*(3569.f/40320)))))))));
        dphi[9] = (64188125.f/504+t*(231125375.f/1008
              +t*(17306975.f/96+t*(7761805.f/96+t*(2895587.f/128
                    +t*(129391.f/32+t*(259715.f/576+t*(28909.f/1008
                          +t*(3569.f/4480)))))))));
        break;
      default:
        NAMD_die("Unknown MSM approximation.");
    } // switch
  }

}; // ComputeMsmMgr


// Degree of polynomial basis function Phi.
const int ComputeMsmMgr::PolyDegree[NUM_APPROX] = {
  3, 5, 5, 7, 7, 9, 9,
};

// The stencil array lengths below.
const int ComputeMsmMgr::Nstencil[NUM_APPROX] = {
  5, 7, 7, 9, 9, 11, 11,
};

// Index offsets from the stencil-centered grid element, to get
// to the correct contributing grid element.
const int ComputeMsmMgr::IndexOffset[NUM_APPROX][MAX_NSTENCIL_SKIP_ZERO] = {
  // cubic
  {-3, -1, 0, 1, 3},

  // quintic C1
  {-5, -3, -1, 0, 1, 3, 5},

  // quintic C2  (same as quintic C1)
  {-5, -3, -1, 0, 1, 3, 5},

  // septic C1
  {-7, -5, -3, -1, 0, 1, 3, 5, 7},

  // septic C3  (same as septic C1)
  {-7, -5, -3, -1, 0, 1, 3, 5, 7},

  // nonic C1
  {-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9},

  // nonic C4  (same as nonic C1)
  {-9, -7, -5, -3, -1, 0, 1, 3, 5, 7, 9},
};

// The grid transfer stencils for the non-factored restriction and
// prolongation procedures.
const Float ComputeMsmMgr::PhiStencil[NUM_APPROX][MAX_NSTENCIL_SKIP_ZERO] = {
  // cubic
  {-1.f/16, 9.f/16, 1, 9.f/16, -1.f/16},

  // quintic C1
  {3.f/256, -25.f/256, 75.f/128, 1, 75.f/128, -25.f/256, 3.f/256},

  // quintic C2  (same as quintic C1)
  {3.f/256, -25.f/256, 75.f/128, 1, 75.f/128, -25.f/256, 3.f/256},

  // septic C1
  { -5.f/2048, 49.f/2048, -245.f/2048, 1225.f/2048, 1, 1225.f/2048,
    -245.f/2048, 49.f/2048, -5.f/2048 },

  // septic C3  (same as septic C3)
  { -5.f/2048, 49.f/2048, -245.f/2048, 1225.f/2048, 1, 1225.f/2048,
    -245.f/2048, 49.f/2048, -5.f/2048 },

  // nonic C1
  { 35.f/65536, -405.f/65536, 567.f/16384, -2205.f/16384, 
    19845.f/32768, 1, 19845.f/32768, -2205.f/16384, 567.f/16384, 
    -405.f/65536, 35.f/65536 },

  // nonic C4  (same as nonic C1)
  { 35.f/65536, -405.f/65536, 567.f/16384, -2205.f/16384, 
    19845.f/32768, 1, 19845.f/32768, -2205.f/16384, 567.f/16384, 
    -405.f/65536, 35.f/65536 },
};


namespace msm {

  //
  // PatchData
  //
 
  struct PatchData {
    ComputeMsmMgr *mgr;
    Map *map;
    PatchDiagram *pd;
    AtomCoordArray coord;
    ForceArray force;
    Grid<Float> qh;
    Grid<Float> eh;
    Grid<Float> subgrid;
    BigReal energy;
    BigReal virial[3][3];
    int cntRecvs;
    int patchID;

    AtomCoordArray& coordArray() { return coord; }
    ForceArray& forceArray() { return force; }

    PatchData(ComputeMsmMgr *pmgr, int pid);
    void init(int natoms);
    void anterpolation();
    void sendCharge();
    void addPotential(const Grid<Float>& epart);
    void interpolation();
  };

} // namespace msm


/////////////////
//
// MsmGridCutoff

class MsmGridCutoff : public CBase_MsmGridCutoff {
  public:
    CProxy_ComputeMsmMgr mgrProxy;
    ComputeMsmMgr *mgrLocal;     // for quick access to data
    msm::Map *map;
    msm::BlockIndex qhblockIndex;  // source of charges
    msm::BlockSend ehblockSend;    // destination for potentials
    msm::Grid<Float> qh;
    msm::Grid<Float> eh;

    MsmGridCutoff() {
      mgrProxy = CProxy_ComputeMsmMgr(CkpvAccess(BOCclass_group).computeMsmMgr);
      mgrLocal = CProxy_ComputeMsmMgr::ckLocalBranch(
          CkpvAccess(BOCclass_group).computeMsmMgr);
      map = &(mgrLocal->mapData());
#ifdef MSM_TIMING
      mgrLocal->addTiming();
#endif
#ifdef MSM_PROFILING
      mgrLocal->addProfiling();
#endif
    }

    MsmGridCutoff(CkMigrateMessage *m) { }

    void initialize(MsmGridCutoffInitMsg *bmsg) {
      qhblockIndex = bmsg->qhBlockIndex;
      ehblockSend = bmsg->ehBlockSend;
      delete bmsg;
#ifdef DEBUG_MSM_GRID
      printf("MsmGridCutoff[%d]:  initialize()"
          " send to level=%d block=(%d,%d,%d)\n",
          thisIndex, ehblockSend.nblock_wrap.level,
          ehblockSend.nblock_wrap.n.i,
          ehblockSend.nblock_wrap.n.j,
          ehblockSend.nblock_wrap.n.k);
#endif
    }

#define CLIP_POTENTIAL_GRID
#undef CLIP_POTENTIAL_GRID

    void compute(GridFloatMsg *gmsg) {
#ifdef DEBUG_MSM_GRID
      printf("MsmGridCutoff %d:  compute()\n", thisIndex);
#endif
#ifdef MSM_TIMING
      double startTime, stopTime;
      startTime = CkWallTimer();
#endif
      //
      // receive block of charges
      //
      int priority = mgrLocal->nlevels
        + 2*(mgrLocal->nlevels - ehblockSend.nblock_wrap.level) - 1;
      int pid;
      // qh is resized only the first time, memory allocation persists
      gmsg->get(qh, pid);
      delete gmsg;
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif

      //
      // grid cutoff calculation
      // this charge block -> this potential block
      //

#ifdef MSM_TIMING
      startTime = stopTime;
#endif
      // resets indexing on block
      // eh is resized only the first time, memory allocation persists
      eh.init(ehblockSend.nrange);
      eh.reset(0);
      // destination block index for potentials
      const msm::BlockIndex& bindex = ehblockSend.nblock_wrap;
      // need grid of weights for this level
      msm::Grid<Float>& gc = map->gc[bindex.level];
      // index range of weights
      int gia = gc.ia();
      int gib = gc.ib();
      int gja = gc.ja();
      int gjb = gc.jb();
      int gka = gc.ka();
      int gkb = gc.kb();
      int gni = gc.ni();
      int gnj = gc.nj();
      // index range of charge grid
      int qia = qh.ia();
      int qib = qh.ib();
      int qja = qh.ja();
      int qjb = qh.jb();
      int qka = qh.ka();
      int qkb = qh.kb();
      int qni = qh.ni();
      int qnj = qh.nj();
      // index range of potentials
      int ia = eh.ia();
      int ib = eh.ib();
      int ja = eh.ja();
      int jb = eh.jb();
      int ka = eh.ka();
      int kb = eh.kb();
#ifdef CLIP_POTENTIAL_GRID
      int eni = eh.ni();
      int enj = eh.nj();
#endif
      // access buffers directly
      const Float *gcbuffer = gc.data().buffer();
      const Float *qhbuffer = qh.data().buffer();
      Float *ehbuffer = eh.data().buffer();

#ifdef CLIP_POTENTIAL_GRID
      // clip potential grid
      int eia = (qia + gia < ia ? ia : qia + gia);
      int eib = (qib + gib > ib ? ib : qib + gib);
      int eja = (qja + gja < ja ? ja : qja + gja);
      int ejb = (qjb + gjb > jb ? jb : qjb + gjb);
      int eka = (qka + gka < ka ? ka : qka + gka);
      int ekb = (qkb + gkb > kb ? kb : qkb + gkb);
#else
      int eia = ia;
      int eib = ib;
      int eja = ja;
      int ejb = jb;
      int eka = ka;
      int ekb = kb;

      int index = 0;
#endif

      // loop over potentials
      for (int k = eka;  k <= ekb;  k++) {
        // clip charges to weights along k
        int mka = ( qka >= gka + k ? qka : gka + k );
        int mkb = ( qkb <= gkb + k ? qkb : gkb + k );
#ifdef CLIP_POTENTIAL_GRID
        int ekoff = (k - eka) * enj;
#endif

        for (int j = eja;  j <= ejb;  j++) {
          // clip charges to weights along j
          int mja = ( qja >= gja + j ? qja : gja + j );
          int mjb = ( qjb <= gjb + j ? qjb : gjb + j );
#ifdef CLIP_POTENTIAL_GRID
          int ejkoff = (ekoff + (j - eja)) * eni;
#endif

          for (int i = eia;  i <= eib;  i++) {
            // clip charges to weights along i
            int mia = ( qia >= gia + i ? qia : gia + i );
            int mib = ( qib <= gib + i ? qib : gib + i );
#ifdef CLIP_POTENTIAL_GRID
            int eijkoff = ejkoff + (i - eia);
#endif

            // accumulate sum to this eh point
            Float ehsum = 0;

#if 0
            // loop over charge grid
            for (int qk = mka;  qk <= mkb;  qk++) {
              int qkoff = (qk - qka) * qnj;
              int gkoff = ((qk-k) - gka) * gnj;

              for (int qj = mja;  qj <= mjb;  qj++) {
                int qjkoff = (qkoff + qj - qja) * qni;
                int gjkoff = (gkoff + (qj-j) - gja) * gni;

// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                for (int qi = mia;  qi <= mib;  qi++) {
                  int qijkoff = qjkoff + qi - qia;
                  int gijkoff = gjkoff + (qi-i) - gia;

                  ehsum += gcbuffer[gijkoff] * qhbuffer[qijkoff];
                }
              }
            } // end loop over charge grid
#else

#if 0
            // loop over charge grid
            int nn = mib - mia + 1;
            for (int qk = mka;  qk <= mkb;  qk++) {
              int qkoff = (qk - qka) * qnj;
              int gkoff = ((qk-k) - gka) * gnj;

              for (int qj = mja;  qj <= mjb;  qj++) {
                int qjkoff = (qkoff + qj - qja) * qni;
                int gjkoff = (gkoff + (qj-j) - gja) * gni;

                const Float *qbuf = qhbuffer + (qjkoff - qia + mia);
                const Float *gbuf = gcbuffer + (gjkoff - i - gia + mia);
#ifdef MSM_PROFILING
                mgrLocal->xLoopCnt[nn]++;
#endif
// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                for (int ii = 0;  ii < nn;  ii++) {
                  ehsum += gbuf[ii] * qbuf[ii];
                }
              }
            } // end loop over charge grid
#else
            // loop over charge grid
            int nn = mib - mia + 1;
            if (nn == 8) {  // hard coded inner loop = 8
              int qnji = qnj * qni;
              int qkoff = -qka*qnji - qja*qni - qia + mia;
              int gnji = gnj * gni;
              int gkoff = (-k-gka)*gnji + (-j-gja)*gni - i - gia + mia;

              for (int qk = mka;  qk <= mkb;  qk++) {
                int qjkoff = qkoff + qk*qnji;
                int gjkoff = gkoff + qk*gnji;

                for (int qj = mja;  qj <= mjb;  qj++) {
                  const Float *qbuf = qhbuffer + (qjkoff + qj*qni);
                  const Float *gbuf = gcbuffer + (gjkoff + qj*gni);
#ifdef MSM_PROFILING
                  mgrLocal->xLoopCnt[nn]++;
#endif
// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                  for (int ii = 0;  ii < 8;  ii++) {
                    ehsum += gbuf[ii] * qbuf[ii];
                  }
                }
              } // end loop over charge grid
            }
            else {  // variable length inner loop < 8
              int qnji = qnj * qni;
              int qkoff = -qka*qnji - qja*qni - qia + mia;
              int gnji = gnj * gni;
              int gkoff = (-k-gka)*gnji + (-j-gja)*gni - i - gia + mia;

              for (int qk = mka;  qk <= mkb;  qk++) {
                int qjkoff = qkoff + qk*qnji;
                int gjkoff = gkoff + qk*gnji;

                for (int qj = mja;  qj <= mjb;  qj++) {
                  const Float *qbuf = qhbuffer + (qjkoff + qj*qni);
                  const Float *gbuf = gcbuffer + (gjkoff + qj*gni);
#ifdef MSM_PROFILING
                  mgrLocal->xLoopCnt[nn]++;
#endif
// help the vectorizer make reasonable decisions
#if defined(__INTEL_COMPILER)
#pragma vector always 
#endif
                  for (int ii = 0;  ii < nn;  ii++) {
                    ehsum += gbuf[ii] * qbuf[ii];
                  }
                }
              } // end loop over charge grid
            }
#endif

#endif

#ifdef CLIP_POTENTIAL_GRID
            ehbuffer[eijkoff] = ehsum;
#else
            ehbuffer[index] = ehsum;
            index++;
#endif
          }
        }
      } // end loop over potentials
#ifdef MSM_PROFILING
      mgrLocal->doneProfiling();
#endif

#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::GRIDCUTOFF] += stopTime - startTime;
#endif

      //
      // send block of potentials
      //

#ifdef MSM_TIMING
      startTime = stopTime;
#endif
      // shift grid index range to its true (wrapped) values
      eh.updateLower( ehblockSend.nrange_wrap.lower() );
      // place eh into message
      int nelems = eh.data().len();
#ifdef MSM_FIXED_SIZE_GRID_MSG
      GridFloatMsg *gm = new(sizeof(int)) GridFloatMsg;
#else
      GridFloatMsg *gm = new(nelems, sizeof(int)) GridFloatMsg;
#endif
      *(int *)CkPriorityPtr(gm) = priority;
      CkSetQueueing(gm, CK_QUEUEING_IFIFO);
      gm->put(eh, bindex.level);
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
      // lookup in ComputeMsmMgr proxy array by level
      mgrLocal->msmBlock[bindex.level](
          bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);
#ifdef MSM_TIMING
      mgrLocal->doneTiming();
#endif
    } // compute()

};

// MsmGridCutoff
//
/////////////////


/////////////////
//
// MsmBlock

class MsmBlock : public CBase_MsmBlock {
  public:
    CProxy_ComputeMsmMgr mgrProxy;
    ComputeMsmMgr *mgrLocal;  // for quick access to data
    msm::Map *map;
    msm::BlockDiagram *bd;
    msm::Grid<Float> qh;
    msm::Grid<Float> eh;
#ifndef MSM_GRID_CUTOFF_DECOMP
    msm::Grid<Float> ehCutoff;
#endif
    msm::Grid<Float> qhRestricted;
    msm::Grid<Float> ehProlongated;
    int cntRecvsCharge;
    int cntRecvsPotential;
    msm::BlockIndex blockIndex;
 
    //msm::Grid<Float> qpart, epart;
    msm::Grid<Float> subgrid;

    MsmBlock(int level);
    MsmBlock(CkMigrateMessage *m) { }

    void init();

    void addCharge(GridFloatMsg *);  // entry

    void restriction();
    void sendUpCharge();
    void gridCutoff();
#ifndef MSM_GRID_CUTOFF_DECOMP
    void sendAcrossPotential();
#endif

    void addPotential(GridFloatMsg *);  // entry

    void prolongation();
    void sendDownPotential();
    void sendPatch();
};

MsmBlock::MsmBlock(int level) {
  blockIndex.level = level;
  blockIndex.n = msm::Ivec(thisIndex.x, thisIndex.y, thisIndex.z);
  mgrProxy = CProxy_ComputeMsmMgr(CkpvAccess(BOCclass_group).computeMsmMgr);
  mgrLocal = CProxy_ComputeMsmMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMgr);
  map = &(mgrLocal->mapData());
  bd = &(map->blockLevel[blockIndex.level](blockIndex.n));
  qh.init( bd->nrange );
  eh.init( bd->nrange );
#ifndef MSM_GRID_CUTOFF_DECOMP
  ehCutoff.init( bd->nrangeCutoff );
#endif
  qhRestricted.init( bd->nrangeRestricted );
  ehProlongated.init( bd->nrangeProlongated );
#ifdef DEBUG_MSM_GRID
  printf("MsmBlock level=%d, n=%d %d %d:  constructor\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif
#ifdef MSM_TIMING
  mgrLocal->addTiming();
#endif
  init();
}

void MsmBlock::init() {
  qh.reset(0);
  eh.reset(0);
#ifndef MSM_GRID_CUTOFF_DECOMP
  ehCutoff.reset(0);
#endif
  qhRestricted.reset(0);
  ehProlongated.reset(0);
  cntRecvsCharge = 0;
  cntRecvsPotential = 0;
}

void MsmBlock::addCharge(GridFloatMsg *gm)
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  //msm::Grid<BigReal> qpart;
  int pid;
  gm->get(subgrid, pid);
  delete gm;
  qh += subgrid;
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
  if (++cntRecvsCharge == bd->numRecvsCharge) {
    int nlevels = mgrLocal->numLevels();
    if (blockIndex.level < nlevels-1) {
      restriction();
    }
    gridCutoff();
  }
}

void MsmBlock::restriction()
{
#ifdef DEBUG_MSM_GRID
  printf("MsmBlock level=%d, id=%d %d %d:  restriction\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif

#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  // stencil data for approximating charge on restricted grid
  const int approx = mgrLocal->approx;
  const int nstencil = ComputeMsmMgr::Nstencil[approx];
  const int *offset = ComputeMsmMgr::IndexOffset[approx];
  const Float *phi = ComputeMsmMgr::PhiStencil[approx];

  // index range for h grid charges
  int ia1 = qh.ia();
  int ib1 = qh.ib();
  int ja1 = qh.ja();
  int jb1 = qh.jb();
  int ka1 = qh.ka();
  int kb1 = qh.kb();

  // index range for restricted (2h) grid charges
  int ia2 = qhRestricted.ia();
  int ib2 = qhRestricted.ib();
  int ja2 = qhRestricted.ja();
  int jb2 = qhRestricted.jb();
  int ka2 = qhRestricted.ka();
  int kb2 = qhRestricted.kb();

  // loop over restricted (2h) grid
  for (int k2 = ka2;  k2 <= kb2;  k2++) {
    int k1 = 2 * k2;
    for (int j2 = ja2;  j2 <= jb2;  j2++) {
      int j1 = 2 * j2;
      for (int i2 = ia2;  i2 <= ib2;  i2++) {
        int i1 = 2 * i2;

        // loop over stencils on h grid
        Float q2hsum = 0;  // sum charge to restricted (2h) grid point
        for (int k = 0;  k < nstencil;  k++) {
          int kn = k1 + offset[k];
          if      (kn < ka1) continue;
          else if (kn > kb1) break;

          for (int j = 0;  j < nstencil;  j++) {
            int jn = j1 + offset[j];
            if      (jn < ja1) continue;
            else if (jn > jb1) break;

            for (int i = 0;  i < nstencil;  i++) {
              int in = i1 + offset[i];
              if      (in < ia1) continue;
              else if (in > ib1) break;

              q2hsum += qh(in,jn,kn) * phi[i] * phi[j] * phi[k];
            }
          }
        } // end loop over stencils on h grid

        qhRestricted(i2,j2,k2) = q2hsum;
      }
    }
  } // end loop over restricted (2h) grid
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::RESTRICT] += stopTime - startTime;
#endif

  sendUpCharge();
}

void MsmBlock::sendUpCharge()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
#endif
  int lnext = blockIndex.level + 1;
  // buffer portions of grid to send to Blocks on next level
  // allocate the largest buffer space we'll need
  //msm::Grid<BigReal> subgrid;
  //subgrid.resize(map->bsx[lnext] * map->bsy[lnext] * map->bsz[lnext]);
  for (int n = 0;  n < bd->sendUp.len();  n++) {
#ifdef MSM_TIMING
    startTime = CkWallTimer();
#endif
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendUp[n].nrange );
    // extract the values from the larger grid into the subgrid
    qhRestricted.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendUp[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendUp[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    // SET MESSAGE PRIORITY
    int nelems = subgrid.data().len();
#ifdef MSM_FIXED_SIZE_GRID_MSG
    GridFloatMsg *gm = new(sizeof(int)) GridFloatMsg;
#else
    GridFloatMsg *gm = new(nelems, sizeof(int)) GridFloatMsg;
#endif
    *(int *)CkPriorityPtr(gm) = lnext;
    CkSetQueueing(gm, CK_QUEUEING_IFIFO);
    gm->put(subgrid, bindex.level);
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addCharge(gm);
  } // for
}

void MsmBlock::gridCutoff()
{
#ifdef DEBUG_MSM_GRID
  printf("MsmBlock level=%d, id=%d %d %d:  grid cutoff\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif
#ifndef MSM_GRID_CUTOFF_DECOMP
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  // need grid of weights for this level
  msm::Grid<Float>& gc = map->gc[blockIndex.level];
  // index range of weights
  int gia = gc.ia();
  int gib = gc.ib();
  int gja = gc.ja();
  int gjb = gc.jb();
  int gka = gc.ka();
  int gkb = gc.kb();
  // index range of charge grid
  int qia = qh.ia();
  int qib = qh.ib();
  int qja = qh.ja();
  int qjb = qh.jb();
  int qka = qh.ka();
  int qkb = qh.kb();
  // index range of potentials
  int ia = ehCutoff.ia();
  int ib = ehCutoff.ib();
  int ja = ehCutoff.ja();
  int jb = ehCutoff.jb();
  int ka = ehCutoff.ka();
  int kb = ehCutoff.kb();
  // loop over potentials
  for (int k = ka;  k <= kb;  k++) {
    for (int j = ja;  j <= jb;  j++) {
      for (int i = ia;  i <= ib;  i++) {
        // clip charges to weights
        int mia = ( qia >= gia + i ? qia : gia + i );
        int mib = ( qib <= gib + i ? qib : gib + i );
        int mja = ( qja >= gja + j ? qja : gja + j );
        int mjb = ( qjb <= gjb + j ? qjb : gjb + j );
        int mka = ( qka >= gka + k ? qka : gka + k );
        int mkb = ( qkb <= gkb + k ? qkb : gkb + k );
        // accumulate sum to this eh point
        Float ehsum = 0;
        // loop over smaller charge grid
        for (int qk = mka;  qk <= mkb;  qk++) {
          for (int qj = mja;  qj <= mjb;  qj++) {
            for (int qi = mia;  qi <= mib;  qi++) {
              ehsum += gc(qi-i, qj-j, qk-k) * qh(qi,qj,qk);
            }
          }
        } // end loop over smaller charge grid
        ehCutoff(i,j,k) = ehsum;
      }
    }
  } // end loop over potentials
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::GRIDCUTOFF] += stopTime - startTime;
#endif

  sendAcrossPotential();

#else

  // send charge block to MsmGridCutoff compute objects
#ifdef MSM_TIMING
  double startTime, stopTime;
#endif
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level)-1;
  int nelems = qh.data().len();
  int len = bd->indexGridCutoff.len();
  for (int n = 0;  n < len;  n++) {
#ifdef MSM_TIMING
    startTime = CkWallTimer();
#endif
    int index = bd->indexGridCutoff[n];
#ifdef MSM_FIXED_SIZE_GRID_MSG
    GridFloatMsg *gm = new(sizeof(int)) GridFloatMsg;
#else
    GridFloatMsg *gm = new(nelems, sizeof(int)) GridFloatMsg;
#endif
    *(int *)CkPriorityPtr(gm) = priority;
    CkSetQueueing(gm, CK_QUEUEING_IFIFO);
    gm->put(qh, blockIndex.level);
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    mgrLocal->msmGridCutoff[index].compute(gm);
  }
#endif
}

#ifndef MSM_GRID_CUTOFF_DECOMP
void MsmBlock::sendAcrossPotential()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
#endif
  int lnext = blockIndex.level;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level)-1;
  // buffer portions of grid to send to Blocks on this level
  // allocate the largest buffer space we'll need
  //msm::Grid<BigReal> subgrid;
  //subgrid.resize(map->bsx[lnext] * map->bsy[lnext] * map->bsz[lnext]);
  for (int n = 0;  n < bd->sendAcross.len();  n++) {
#ifdef MSM_TIMING
    startTime = CkWallTimer();
#endif
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendAcross[n].nrange );
    // extract the values from the larger grid into the subgrid
    ehCutoff.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendAcross[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendAcross[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    int nelems = subgrid.data().len();
#ifdef MSM_FIXED_SIZE_GRID_MSG
    GridFloatMsg *gm = new(sizeof(int)) GridFloatMsg;
#else
    GridFloatMsg *gm = new(nelems, sizeof(int)) GridFloatMsg;
#endif
    *(int *)CkPriorityPtr(gm) = priority;
    CkSetQueueing(gm, CK_QUEUEING_IFIFO);
    gm->put(subgrid, bindex.level);
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);
  } // for
}
#endif

void MsmBlock::addPotential(GridFloatMsg *gm)
{
#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  //msm::Grid<Float> epart;
  int pid;
  gm->get(subgrid, pid);
  delete gm;
  eh += subgrid;
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
  if (++cntRecvsPotential == bd->numRecvsPotential) {
    if (blockIndex.level > 0) {
      prolongation();
    }
    else {
      sendPatch();
    }
  }
}

void MsmBlock::prolongation()
{
#ifdef DEBUG_MSM_GRID
  printf("MsmBlock level=%d, id=%d %d %d:  prolongation\n",
      blockIndex.level, blockIndex.n.i, blockIndex.n.j, blockIndex.n.k);
#endif

#ifdef MSM_TIMING
  double startTime, stopTime;
  startTime = CkWallTimer();
#endif
  // stencil data for approximating potential on prolongated grid
  const int approx = mgrLocal->approx;
  const int nstencil = ComputeMsmMgr::Nstencil[approx];
  const int *offset = ComputeMsmMgr::IndexOffset[approx];
  const Float *phi = ComputeMsmMgr::PhiStencil[approx];

  // index range for prolongated h grid potentials
  int ia1 = ehProlongated.ia();
  int ib1 = ehProlongated.ib();
  int ja1 = ehProlongated.ja();
  int jb1 = ehProlongated.jb();
  int ka1 = ehProlongated.ka();
  int kb1 = ehProlongated.kb();

  // index range for 2h grid potentials
  int ia2 = eh.ia();
  int ib2 = eh.ib();
  int ja2 = eh.ja();
  int jb2 = eh.jb();
  int ka2 = eh.ka();
  int kb2 = eh.kb();

  // loop over 2h grid
  for (int k2 = ka2;  k2 <= kb2;  k2++) {
    int k1 = 2 * k2;
    for (int j2 = ja2;  j2 <= jb2;  j2++) {
      int j1 = 2 * j2;
      for (int i2 = ia2;  i2 <= ib2;  i2++) {
        int i1 = 2 * i2;

        // loop over stencils on prolongated h grid
        for (int k = 0;  k < nstencil;  k++) {
          int kn = k1 + offset[k];
          if      (kn < ka1) continue;
          else if (kn > kb1) break;

          for (int j = 0;  j < nstencil;  j++) {
            int jn = j1 + offset[j];
            if      (jn < ja1) continue;
            else if (jn > jb1) break;

            for (int i = 0;  i < nstencil;  i++) {
              int in = i1 + offset[i];
              if      (in < ia1) continue;
              else if (in > ib1) break;

              ehProlongated(in,jn,kn) +=
                eh(i2,j2,k2) * phi[i] * phi[j] * phi[k];
            }
          }
        } // end loop over stencils on prolongated h grid

      }
    }
  } // end loop over 2h grid
#ifdef MSM_TIMING
  stopTime = CkWallTimer();
  mgrLocal->msmTiming[MsmTimer::PROLONGATE] += stopTime - startTime;
#endif

  sendDownPotential();
}

void MsmBlock::sendDownPotential()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
#endif
  int lnext = blockIndex.level - 1;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level);
  // buffer portions of grid to send to Blocks on next level
  // allocate the largest buffer space we'll need
  //msm::Grid<BigReal> subgrid;
  //subgrid.resize(map->bsx[lnext] * map->bsy[lnext] * map->bsz[lnext]);
  for (int n = 0;  n < bd->sendDown.len();  n++) {
#ifdef MSM_TIMING
    startTime = CkWallTimer();
#endif
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendDown[n].nrange );
    // extract the values from the larger grid into the subgrid
    ehProlongated.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendDown[n].nrange_wrap.lower() );
    // add the subgrid charges into the block
    msm::BlockIndex& bindex = bd->sendDown[n].nblock_wrap;
    ASSERT(bindex.level == lnext);
    // place subgrid into message
    int nelems = subgrid.data().len();
#ifdef MSM_FIXED_SIZE_GRID_MSG
    GridFloatMsg *gm = new(sizeof(int)) GridFloatMsg;
#else
    GridFloatMsg *gm = new(nelems, sizeof(int)) GridFloatMsg;
#endif
    *(int *)CkPriorityPtr(gm) = priority;
    CkSetQueueing(gm, CK_QUEUEING_IFIFO);
    gm->put(subgrid, bindex.level);
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    // lookup in ComputeMsmMgr proxy array by level
    mgrLocal->msmBlock[lnext](
        bindex.n.i, bindex.n.j, bindex.n.k).addPotential(gm);
  } // for
#ifdef MSM_TIMING
  mgrLocal->doneTiming();
#endif
  init();  // reinitialize for next computation
}

void MsmBlock::sendPatch()
{
#ifdef MSM_TIMING
  double startTime, stopTime;
#endif
  int lnext = blockIndex.level;
  int priority = mgrLocal->nlevels + 2*(mgrLocal->nlevels - blockIndex.level);
  ASSERT(lnext == 0);
  // buffer portions of grid to send to Blocks on next level
  // allocate the largest buffer space we'll need
  //msm::Grid<Float> subgrid;
  //subgrid.resize(map->bsx[lnext] * map->bsy[lnext] * map->bsz[lnext]);
  for (int n = 0;  n < bd->sendPatch.len();  n++) {
#ifdef MSM_TIMING
    startTime = CkWallTimer();
#endif
    // initialize the proper subgrid indexing range
    subgrid.init( bd->sendPatch[n].nrange );
    // extract the values from the larger grid into the subgrid
    eh.extract(subgrid);
    // translate the subgrid indexing range to match the MSM block
    subgrid.updateLower( bd->sendPatch[n].nrange_unwrap.lower() );
    // add the subgrid charges into the block, need its patch ID
    int pid = bd->sendPatch[n].patchID;
    // place subgrid into message
    int nelems = subgrid.data().len();
#ifdef MSM_FIXED_SIZE_GRID_MSG
    GridFloatMsg *gm = new(sizeof(int)) GridFloatMsg;
#else
    GridFloatMsg *gm = new(nelems, sizeof(int)) GridFloatMsg;
#endif
    *(int *)CkPriorityPtr(gm) = priority;
    CkSetQueueing(gm, CK_QUEUEING_IFIFO);
    gm->put(subgrid, pid);
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgrLocal->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    // lookup which PE has this patch
    PatchMap *pm = PatchMap::Object();
    int pe = pm->node(pid);
    mgrProxy[pe].addPotential(gm);
  }
#ifdef MSM_TIMING
  mgrLocal->doneTiming();
#endif
  init();  // reinitialize for next computation
}

// MsmBlock
//
//////////////////


ComputeMsmMgr::ComputeMsmMgr() :
  msmProxy(thisgroup), msmCompute(0)
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  (constructor) PE %d\n", CkMyPe());
#endif
  CkpvAccess(BOCclass_group).computeMsmMgr = thisgroup;

#ifdef MSM_TIMING
  if (CkMyPe() == 0) {
    msmTimer = CProxy_MsmTimer::ckNew();
  }
  initTiming();
#endif
#ifdef MSM_PROFILING
  if (CkMyPe() == 0) {
    msmProfiler = CProxy_MsmProfiler::ckNew();
  }
  initProfiling();
#endif
}

ComputeMsmMgr::~ComputeMsmMgr()
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  (destructor) PE %d\n", CkMyPe());
#endif
  // free memory?
}


//
// Given basis vector length "len" (and with user grid spacing)
// If using periodic boundary conditions along this basis vector,
// h is calculated to be close to desired grid spacing such that 
// nn = 2^k or 3*2^k.  For non-periodic boundaries, we can set h
// to the desired grid spacing, and set ia and ib to pad 1/2 the 
// interpolating stencil width.  
//
void ComputeMsmMgr::setup_hgrid_1d(BigReal len, BigReal& hh, int& nn,
    int& ia, int& ib, int isperiodic)
{
  ASSERT(gridspacing > 0);
  if (isperiodic) {
    const BigReal hmin = (4./5) * gridspacing;
    const BigReal hmax = 1.5 * hmin;
    hh = len;
    nn = 1;  // start with one grid point across length
    while (hh >= hmax) {
      hh *= 0.5;  // halve spacing and double grid points
      nn <<= 1;
    }
    if (hh < hmin) {
      if (nn < 4) {
        NAMD_die("Basis vector is too short or MSM grid spacing is too large");
      }
      hh *= (4./3);  // scale hh by 4/3 and nn by 3/4
      nn >>= 2;
      nn *= 3;
    }
    // now we have:  hmin <= h < hmax,
    // where nn is a power of 2 times no more than one power of 3
    ia = 0;
    ib = nn-1;
  }
  else {
    hh = gridspacing;
    // Instead of "nn = (int) ceil(len / hh);"
    // len is divisible by hh, up to roundoff error, so round to closest nn
    nn = (int) floor(len/hh + 0.5);
    ia = -s_edge;
    ib = nn + s_edge;
  }
}


// make sure that block sizes divide evenly into periodic dimensions
// call only for periodic dimensions
void ComputeMsmMgr::setup_periodic_blocksize(int& bsize, int n)
{
  if (n % bsize != 0) {
    // n is either 2^k or 3*2^k
    int newbsize = 1;
    if (n % 3 == 0) newbsize = 3;
    while (newbsize < bsize && newbsize < n) newbsize *= 2;
    if (bsize < newbsize) newbsize /= 2;
    if (n % newbsize != 0) {
      NAMD_die("MSM grid size for periodic dimensions must be "
          "a power of 2 times at most one power of 3");
    }
    bsize = newbsize;
  }
  return;
}


void ComputeMsmMgr::initialize(MsmInitMsg *msg)
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  initialize() PE %d\n", CkMyPe());
#endif

  smin = msg->smin;
  smax = msg->smax;
  delete msg;

#ifdef DEBUG_MSM_VERBOSE
  printf("smin = %g %g %g  smax = %g %g %g\n",
      smin.x, smin.y, smin.z, smax.x, smax.y, smax.z);
#endif

  SimParameters *simParams = Node::Object()->simParameters;

  // get required sim params, check validity
  lattice = simParams->lattice;

  a = simParams->cutoff;

  gridspacing = simParams->MSMGridSpacing;
  if (gridspacing <= 0 || gridspacing >= a) {
    NAMD_die("MSM: illegal grid spacing requested (MSMGridSpacing)");
  }

  padding = simParams->MSMPadding;
  if (padding < 0) {
    NAMD_die("MSM: illegal padding requested (MSMPadding)");
  }

  approx = simParams->MSMApprox;
  if (approx < 0 || approx >= NUM_APPROX) {
    NAMD_die("MSM: unknown approximation requested (MSMApprox)");
  }
  split = simParams->MSMSplit;
  if (split < 0 || split >= NUM_SPLIT) {
    NAMD_die("MSM: unknown splitting requested (MSMSplit)");
  }
  if (CkMyPe() == 0) {
    const char *approx_str, *split_str;
    switch (approx) {
      case CUBIC:    approx_str = "C1 cubic";   break;
      case QUINTIC:  approx_str = "C1 quintic"; break;
      case QUINTIC2: approx_str = "C2 quintic"; break;
      case SEPTIC:   approx_str = "C1 septic";  break;
      case SEPTIC3:  approx_str = "C3 septic";  break;
      case NONIC:    approx_str = "C1 nonic";   break;
      case NONIC4:   approx_str = "C4 nonic";   break;
      default:       approx_str = "unknown";    break;
    }
    switch (split) {
      case TAYLOR2:  split_str = "C2 Taylor";   break;
      case TAYLOR3:  split_str = "C3 Taylor";   break;
      case TAYLOR4:  split_str = "C4 Taylor";   break;
      case TAYLOR5:  split_str = "C5 Taylor";   break;
      case TAYLOR6:  split_str = "C6 Taylor";   break;
      case TAYLOR7:  split_str = "C7 Taylor";   break;
      case TAYLOR8:  split_str = "C8 Taylor";   break;
      default:       split_str = "unknown";     break;
    }
    iout << iINFO << "MSM using "
                  << approx_str << " interpolation\n";
    iout << iINFO << "MSM using "
                  << split_str << " splitting function\n";
  }

  // set maximum number of levels (default 0 adapts levels to system)
  nlevels = simParams->MSMLevels;

  // XXX dispersion unused for now
  dispersion = 0;
  if ( ! dispersion && split >= TAYLOR2_DISP) {
    NAMD_die("MSM: requested splitting for long-range dispersion "
        "(not implemented)");
  }

  // set block sizes for grid decomposition
  int bsx = simParams->MSMBlockSizeX;
  int bsy = simParams->MSMBlockSizeY;
  int bsz = simParams->MSMBlockSizeZ;
  if (bsx <= 0 || bsy <= 0 || bsz <= 0) {
    NAMD_die("MSM: invalid block size requested (MSMBlockSize[XYZ])");
  }
#ifdef MSM_FIXED_SIZE_GRID_MSG
  else if (bsx * bsy * bsz > MSM_MAX_BLOCK_VOLUME) {
    NAMD_die("MSM: requested block size (MSMBlockSize[XYZ]) too big");
  }
#endif
  if (CkMyPe() == 0) {
    iout << iINFO << "MSM block size decomposition along X is "
                  << bsx << " grid points\n";
    iout << iINFO << "MSM block size decomposition along Y is "
                  << bsy << " grid points\n";
    iout << iINFO << "MSM block size decomposition along Z is "
                  << bsz << " grid points\n";
  }

  s_edge = (PolyDegree[approx] - 1) / 2;  // stencil edge size
  omega = 2 * PolyDegree[approx];         // smallest non-periodic grid length

  BigReal xlen, ylen, zlen;
  Vector sgmin, sgmax;  // grid min and max, in scaled coordinates
  int ispx = lattice.a_p();
  int ispy = lattice.b_p();
  int ispz = lattice.c_p();
  int ispany = (ispx || ispy || ispz);  // is there any periodicity?

  if (ispx) {  // periodic along basis vector
    xlen = lattice.a().length();
    sgmax.x = 0.5;
    sgmin.x = -0.5;
  }
  else {  // non-periodic
    sgmax.x = smax.x + padding;  // pad the edges
    sgmin.x = smin.x - padding;
    ASSERT(gridspacing > 0);
    // restrict center to be on a grid point
    BigReal mupper = ceil(sgmax.x / (2*gridspacing));
    BigReal mlower = floor(sgmin.x / (2*gridspacing));
    sgmax.x = 2*gridspacing*mupper;
    sgmin.x = 2*gridspacing*mlower;
    xlen = sgmax.x - sgmin.x;
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("xlen = %g   sgmin.x = %g   sgmax.x = %g\n", xlen, sgmin.x, sgmax.x);
#endif

  if (ispy) {  // periodic along basis vector
    ylen = lattice.b().length();
    sgmax.y = 0.5;
    sgmin.y = -0.5;
  }
  else {  // non-periodic
    sgmax.y = smax.y + padding;  // pad the edges
    sgmin.y = smin.y - padding;
    ASSERT(gridspacing > 0);
    // restrict center to be on a grid point
    BigReal mupper = ceil(sgmax.y / (2*gridspacing));
    BigReal mlower = floor(sgmin.y / (2*gridspacing));
    sgmax.y = 2*gridspacing*mupper;
    sgmin.y = 2*gridspacing*mlower;
    ylen = sgmax.y - sgmin.y;
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("ylen = %g   sgmin.y = %g   sgmax.y = %g\n", ylen, sgmin.y, sgmax.y);
#endif

  if (ispz) {  // periodic along basis vector
    zlen = lattice.c().length();
    sgmax.z = 0.5;
    sgmin.z = -0.5;
  }
  else {  // non-periodic
    sgmax.z = smax.z + padding;  // pad the edges
    sgmin.z = smin.z - padding;
    ASSERT(gridspacing > 0);
    // restrict center to be on a grid point
    BigReal mupper = ceil(sgmax.z / (2*gridspacing));
    BigReal mlower = floor(sgmin.z / (2*gridspacing));
    sgmax.z = 2*gridspacing*mupper;
    sgmin.z = 2*gridspacing*mlower;
    zlen = sgmax.z - sgmin.z;
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("zlen = %g   sgmin.z = %g   sgmax.z = %g\n", zlen, sgmin.z, sgmax.z);
#endif
  sglower = sgmin;

  BigReal hxlen, hylen, hzlen;
  int ia, ib, ja, jb, ka, kb;
  setup_hgrid_1d(xlen, hxlen, nhx, ia, ib, ispx);
  setup_hgrid_1d(ylen, hylen, nhy, ja, jb, ispy);
  setup_hgrid_1d(zlen, hzlen, nhz, ka, kb, ispz);
  if (CkMyPe() == 0) {
    if (ispx || ispy || ispz) {
      iout << iINFO << "MSM grid spacing along X is "<< hxlen << " A\n";
      iout << iINFO << "MSM grid spacing along Y is "<< hylen << " A\n";
      iout << iINFO << "MSM grid spacing along Z is "<< hzlen << " A\n";
    }
    else {
      iout << iINFO << "MSM grid spacing is " << gridspacing << " A\n";
    }
    if ( ! ispx || ! ispy || ! ispz ) {
      iout << iINFO<<"MSM non-periodic padding is "<< padding << " A\n";
    }
  }

  int ni = ib - ia + 1;
  int nj = jb - ja + 1;
  int nk = kb - ka + 1;
  int n;

#if 0
  // reserve temp space for factored grid transfer operation
  n = (nk > omega ? nk : omega);  // row along z-dimension
  lzd.resize(n);
  n *= (nj > omega ? nj : omega);  // plane along yz-dimensions
  lyzd.resize(n);
#endif

  int lastnelems = 1;
  int smallestnbox = 1;
  int isclamped = 0;
  int maxlevels = nlevels;  // user-defined number of levels

#ifdef DEBUG_MSM_VERBOSE
  printf("maxlevels = %d\n", maxlevels);
#endif
  if (nlevels <= 0) {  // instead we set number of levels
    n = ni;
    if (n < nj) n = nj;
    if (n < nk) n = nk;
    for (maxlevels = 1;  n > 0;  n >>= 1)  maxlevels++;
    if (ispany == 0) {  // no periodicity
      // use rule of thumb 3/4 diameter of grid cutoff sphere
      int ngci = (int) ceil(3*a / hxlen) - 1;
      int ngcj = (int) ceil(3*a / hylen) - 1;
      int ngck = (int) ceil(3*a / hzlen) - 1;
      int omega3 = omega * omega * omega;
      int nhalf = (int) sqrt(ni * nj * nk);
      lastnelems = (nhalf > omega3 ? nhalf : omega3);
      smallestnbox = ngci * ngcj * ngck;  // smaller grids don't reduce work
      isclamped = 1;
    }
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("maxlevels = %d\n", maxlevels);
#endif

  // allocate space for storing grid dimensions for each level
  map.gridrange.resize(maxlevels);

  // set periodicity flags
  map.ispx = ispx;
  map.ispy = ispy;
  map.ispz = ispz;

  int level = 0;
  int done = 0;
  int alldone = 0;
  do {
    map.gridrange[level].setbounds(ia, ib, ja, jb, ka, kb);

    // Msm index?

    if (++level == nlevels)  done |= 0x07;  // user limit on levels

    if (isclamped) {
      int nelems = ni * nj * nk;
      if (nelems <= lastnelems)    done |= 0x07;
      if (nelems <= smallestnbox)  done |= 0x07;
    }

    alldone = (done == 0x07);  // make sure all dimensions are done

    if (ispx) {
      ni >>= 1;
      ib = ni-1;
      if (ni & 1)        done |= 0x07;  // == 3 or 1
      else if (ni == 2)  done |= 0x01;  // can do one more
    }
    else {
      ia = -((-ia+1)/2) - s_edge;
      ib = (ib+1)/2 + s_edge;
      ni = ib - ia + 1;
      if (ni <= omega)   done |= 0x01;  // can do more restrictions
    }

    if (ispy) {
      nj >>= 1;
      jb = nj-1;
      if (nj & 1)        done |= 0x07;  // == 3 or 1
      else if (nj == 2)  done |= 0x02;  // can do one more
    }
    else {
      ja = -((-ja+1)/2) - s_edge;
      jb = (jb+1)/2 + s_edge;
      nj = jb - ja + 1;
      if (nj <= omega)   done |= 0x02;  // can do more restrictions
    }

    if (ispz) {
      nk >>= 1;
      kb = nk-1;
      if (nk & 1)        done |= 0x07;  // == 3 or 1
      else if (nk == 2)  done |= 0x04;  // can do one more
    }
    else {
      ka = -((-ka+1)/2) - s_edge;
      kb = (kb+1)/2 + s_edge;
      nk = kb - ka + 1;
      if (nk <= omega)   done |= 0x04;  // can do more restrictions
    }
  } while ( ! alldone );
  nlevels = level;

  // for periodic boundaries, don't visit top level (all 0)
  // toplevel visited only for all nonperiodic boundaries
  int toplevel = (ispany ? nlevels : nlevels - 1);

  // resize down to the actual number of levels (does not change alloc)
  map.gridrange.resize(nlevels);
  map.gc.resize(nlevels);

  // print out some information about MSM
  if (CkMyPe() == 0) {
    iout << iINFO << "MSM using " << nlevels << " levels\n";
    for (n = 0;  n < nlevels;  n++) {
      char s[100];
      snprintf(s, sizeof(s), "    level %d:  "
          "[%d..%d] x [%d..%d] x [%d..%d]\n", n,
          map.gridrange[n].ia(), map.gridrange[n].ib(),
          map.gridrange[n].ja(), map.gridrange[n].jb(),
          map.gridrange[n].ka(), map.gridrange[n].kb());
      iout << iINFO << s;
    }
  }

  // find grid spacing basis vectors
  Vector hu = hxlen * lattice.a().unit();
  Vector hv = hylen * lattice.b().unit();
  Vector hw = hzlen * lattice.c().unit();

  const Vector& ru = lattice.a_r();
  const Vector& rv = lattice.b_r();
  const Vector& rw = lattice.c_r();

  // determine grid spacings in scaled space
  shx = ru * hu;
  shy = rv * hv;
  shz = rw * hw;
  shx_1 = 1 / shx;
  shy_1 = 1 / shy;
  shz_1 = 1 / shz;

  // row vectors to transform interpolated force back to real space
  sx_shx = shx_1 * Vector(ru.x, rv.x, rw.x);
  sy_shy = shy_1 * Vector(ru.y, rv.y, rw.y);
  sz_shz = shz_1 * Vector(ru.z, rv.z, rw.z);
  srx_x = Float(sx_shx.x);
  srx_y = Float(sx_shx.y);
  srx_z = Float(sx_shx.z);
  sry_x = Float(sy_shy.x);
  sry_y = Float(sy_shy.y);
  sry_z = Float(sy_shy.z);
  srz_x = Float(sz_shz.x);
  srz_y = Float(sz_shz.y);
  srz_z = Float(sz_shz.z);

  Vector pu = cross(hv, hw);
  BigReal s = (hu * pu) / (pu * pu);
  pu *= s;  // pu is orthogonal projection of hu onto hv CROSS hw

  Vector pv = cross(hw, hu);
  s = (hv * pv) / (pv * pv);
  pv *= s;  // pv is orthogonal projection of hv onto hw CROSS hu

  Vector pw = cross(hu, hv);
  s = (hw * pw) / (pw * pw);
  pw *= s;  // pw is orthogonal projection of hw onto hu CROSS hv

  // radii for parallelepiped of weights enclosing grid cutoff sphere
  ni = (int) ceil(2*a / pu.length() ) - 1;
  nj = (int) ceil(2*a / pv.length() ) - 1;
  nk = (int) ceil(2*a / pw.length() ) - 1;

  Float scaling = 1;
  Float scaling_factor = 0.5f;
  BigReal a_1 = 1/a;
  BigReal a_p = a_1;
  if (dispersion) {
    a_p = a_p * a_p * a_p;   // = 1/a^3
    a_p = a_p * a_p;         // = 1/a^6
    scaling_factor = 1.f/64;  // = 1/2^6
  }
  int i, j, k;
  for (level = 0;  level < toplevel;  level++) {
    map.gc[level].setbounds(-ni, ni, -nj, nj, -nk, nk);

    for (k = -nk;  k <= nk;  k++) {
      for (j = -nj;  j <= nj;  j++) {
        for (i = -ni;  i <= ni;  i++) {
          if (level == 0) {
            BigReal s, t, gs=0, gt=0, g=0, d=0;
            s = (i*hu + j*hv + k*hw).length() * a_1;
            t = 0.5 * s;
            if (t >= 1) {
              g = 0;
            }
            else {
              splitting(gt, d, t, split);
              if (s >= 1) {
                if (dispersion) {
                  BigReal s_1 = 1/s;
                  gs = s_1 * s_1 * s_1;  // = 1/s^3
                  gs = gs * gs;  // = 1/s^6
                }
                else {
                  gs = 1/s;
                }
              }
              else {
                splitting(gs, d, s, split);
              }
              g = (gs - scaling_factor * gt) * a_p;

              // Msm index?

            }
            map.gc[0](i,j,k) = Float(g);
          } // if level 0
          else {
            map.gc[level](i,j,k) = scaling * map.gc[0](i,j,k);
          }

        } // for i
      } // for j
    } // for k
    scaling *= scaling_factor;

  } // for level

  if (toplevel < nlevels) {
    // nonperiodic along all basis vector directions
    // calculate top level weights where all grid points
    // interact with each other
    ni = map.gridrange[toplevel].ni();
    nj = map.gridrange[toplevel].nj();
    nk = map.gridrange[toplevel].nk();
    map.gc[toplevel].setbounds(-ni, ni, -nj, nj, -nk, nk);

    // Msm index?

    for (k = -nk;  k <= nk;  k++) {
      for (j = -nj;  j <= nj;  j++) {
        for (i = -ni;  i <= ni;  i++) {
          BigReal s, gs, d;
          s = (i*hu + j*hv + k*hw).length() * a_1;
          if (s >= 1) {
            if (dispersion) {
              BigReal s_1 = 1/s;
              gs = s_1 * s_1 * s_1;  // = 1/s^3
              gs = gs * gs;  // = 1/s^6
            }
            else {
              gs = 1/s;
            }
          }
          else {
            splitting(gs, d, s, split);
          }
          map.gc[toplevel](i,j,k) = scaling * Float(gs * a_p);
        } // for i
      } // for j
    } // for k
  } // if toplevel

  // calculate self energy factor for splitting
  BigReal gs=0, d=0;
  splitting(gs, d, 0, split);
  gzero = gs * a_p;

  // allocate map for patches
  PatchMap *pm = PatchMap::Object();
  int numpatches = pm->numPatches();
  map.patchList.resize(numpatches);
#ifdef DEBUG_MSM_VERBOSE
  printf("numPatches = %d\n", numpatches);
#endif

  // allocate map for blocks for each grid level
  map.blockLevel.resize(nlevels);
  map.bsx.resize(nlevels);
  map.bsy.resize(nlevels);
  map.bsz.resize(nlevels);
  for (level = 0;  level < nlevels;  level++) {
    msm::IndexRange& g = map.gridrange[level];
    msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
    int gia = g.ia();
    int gni = g.ni();
    int gja = g.ja();
    int gnj = g.nj();
    int gka = g.ka();
    int gnk = g.nk();
    map.bsx[level] = bsx;
    map.bsy[level] = bsy;
    map.bsz[level] = bsz;
    if (level < nlevels-1 &&
        (map.bsx[level] < gni ||
         map.bsy[level] < gnj ||
         map.bsz[level] < gnk)) {
      // make sure that block sizes divide evenly into periodic dimensions
      if (ispx) setup_periodic_blocksize(map.bsx[level], gni);
      if (ispy) setup_periodic_blocksize(map.bsy[level], gnj);
      if (ispz) setup_periodic_blocksize(map.bsz[level], gnk);
      // subdivide grid into multiple blocks
      //   == ceil(gni / bsx), etc.
      int bni = (gni / map.bsx[level]) + (gni % map.bsx[level] != 0);
      int bnj = (gnj / map.bsy[level]) + (gnj % map.bsy[level] != 0);
      int bnk = (gnk / map.bsz[level]) + (gnk % map.bsz[level] != 0);
      b.set(0, bni, 0, bnj, 0, bnk);
      for (k = 0;  k < bnk;  k++) {
        for (j = 0;  j < bnj;  j++) {
          for (i = 0;  i < bni;  i++) {
            b(i,j,k).reset();
            int ia = gia + i*map.bsx[level];
            int ib = ia + map.bsx[level] - 1;
            int ja = gja + j*map.bsy[level];
            int jb = ja + map.bsy[level] - 1;
            int ka = gka + k*map.bsz[level];
            int kb = ka + map.bsz[level] - 1;
            if (ib >= gia + gni) ib = gia + gni - 1;
            if (jb >= gja + gnj) jb = gja + gnj - 1;
            if (kb >= gka + gnk) kb = gka + gnk - 1;
            b(i,j,k).nrange.setbounds(ia, ib, ja, jb, ka, kb);
          }
        }
      }
    }
    else {
      // make entire grid into single block
      b.set(0, 1, 0, 1, 0, 1);
      b(0,0,0).reset();
      b(0,0,0).nrange.set(gia, gni, gja, gnj, gka, gnk);
      // set correct block dimensions
      map.bsx[level] = gni;
      map.bsy[level] = gnj;
      map.bsz[level] = gnk;
    }
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("Done allocating map for grid levels\n");
  printf("Grid level decomposition:\n");
  for (level = 0;  level < nlevels;  level++) {
    msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
    int bia = b.ia();
    int bib = b.ib();
    int bja = b.ja();
    int bjb = b.jb();
    int bka = b.ka();
    int bkb = b.kb();
    for (k = bka;  k <= bkb;  k++) {
      for (j = bja;  j <= bjb;  j++) {
        for (i = bia;  i <= bib;  i++) {
          int ia = b(i,j,k).nrange.ia();
          int ib = b(i,j,k).nrange.ib();
          int ja = b(i,j,k).nrange.ja();
          int jb = b(i,j,k).nrange.jb();
          int ka = b(i,j,k).nrange.ka();
          int kb = b(i,j,k).nrange.kb();
          printf("level=%d  id=%d %d %d  [%d..%d] x [%d..%d] x [%d..%d]\n",
              level, i, j, k, ia, ib, ja, jb, ka, kb);
        }
      }
    }
  }
#endif

  // initialize grid of PatchDiagram
  // a = cutoff
  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal sysdimc = lattice.c_r().unit() * lattice.c();
  BigReal patchdim = simParams->patchDimension;
  BigReal xmargin = 0.5 * (patchdim - a) / sysdima;
  BigReal ymargin = 0.5 * (patchdim - a) / sysdimb;
  BigReal zmargin = 0.5 * (patchdim - a) / sysdimc;
  int pid;
  for (pid = 0;  pid < numpatches;  pid++) {
    // shortcut reference to this patch diagram
    msm::PatchDiagram& p = map.patchList[pid];
    p.reset();
    // find extent of patch including margin
    BigReal xmin = pm->min_a(pid) - xmargin;
    BigReal xmax = pm->max_a(pid) + xmargin;
    BigReal ymin = pm->min_b(pid) - ymargin;
    BigReal ymax = pm->max_b(pid) + ymargin;
    BigReal zmin = pm->min_c(pid) - zmargin;
    BigReal zmax = pm->max_c(pid) + zmargin;
    // find grid point covering of patch plus outer edge stencil
    int ia = int(floor((xmin - sglower.x) * shx_1)) - s_edge;
    int ib = int(floor((xmax - sglower.x) * shx_1)) + 1 + s_edge;
    int ja = int(floor((ymin - sglower.y) * shy_1)) - s_edge;
    int jb = int(floor((ymax - sglower.y) * shy_1)) + 1 + s_edge;
    int ka = int(floor((zmin - sglower.z) * shz_1)) - s_edge;
    int kb = int(floor((zmax - sglower.z) * shz_1)) + 1 + s_edge;
    // set the index range for this patch's surrounding grid points
    p.nrange.setbounds(ia,ib,ja,jb,ka,kb);
    // find lower and upper blocks of MSM h-grid
    msm::BlockIndex blower = map.blockOfGridIndex(msm::Ivec(ia,ja,ka),0);
    msm::BlockIndex bupper = map.blockOfGridIndex(msm::Ivec(ib,jb,kb),0);
    int maxarrlen = (bupper.n.i - blower.n.i + 1) *
      (bupper.n.j - blower.n.j + 1) * (bupper.n.k - blower.n.k + 1);
    p.send.setmax(maxarrlen);  // allocate space for send array
    // loop over the blocks
    for (int kk = blower.n.k;  kk <= bupper.n.k;  kk++) {
      for (int jj = blower.n.j;  jj <= bupper.n.j;  jj++) {
        for (int ii = blower.n.i;  ii <= bupper.n.i;  ii++) {
          // determine actual block and range to send to
          msm::BlockSend bs;
          bs.nblock.n = msm::Ivec(ii,jj,kk);
          bs.nblock.level = 0;
          bs.nrange = map.clipBlockToIndexRange(bs.nblock, p.nrange);
          map.wrapBlockSend(bs);  // determine wrapping to true block index
          p.send.append(bs);  // append this block to the send array
          // increment counter for receive block
          map.blockLevel[0](bs.nblock_wrap.n).numRecvsCharge++;
          // initialize patch send back from this block
          msm::PatchSend ps;
          ps.nrange = bs.nrange_wrap;
          ps.nrange_unwrap = bs.nrange;
          ps.patchID = pid;
          map.blockLevel[0](bs.nblock_wrap.n).sendPatch.append(ps);
          // increment number of receives back to this patch
          p.numRecvs++;
        }
      }
    }
    // number of receives should be same as number of sends
    ASSERT(p.numRecvs == p.send.len() );
  }
#ifdef DEBUG_MSM_VERBOSE
  printf("Done allocating map for patches\n");
  printf("Patch level decomposition:\n");
  for (pid = 0;  pid < numpatches;  pid++) {
    msm::PatchDiagram& p = map.patchList[pid];
    int ia = p.nrange.ia();
    int ib = p.nrange.ib();
    int ja = p.nrange.ja();
    int jb = p.nrange.jb();
    int ka = p.nrange.ka();
    int kb = p.nrange.kb();
    printf("patch id=%d  [%d..%d] x [%d..%d] x [%d..%d]\n",
        pid, ia, ib, ja, jb, ka, kb);
  }
#endif

  // initialize grid of BlockDiagram for each level
  int polydeg = PolyDegree[approx];
  numGridCutoff = 0;
  for (level = 0;  level < nlevels;  level++) {
    msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
    int bni = b.ni();
    int bnj = b.nj();
    int bnk = b.nk();
#ifdef MSM_SKIP_TOO_DISTANT_BLOCKS
    int bsx = map.bsx[level];
    int bsy = map.bsy[level];
    int bsz = map.bsz[level];
#endif
    for (k = 0;  k < bnk;  k++) {
      for (j = 0;  j < bnj;  j++) {
        for (i = 0;  i < bni;  i++) {
          // Grid cutoff calculation, sendAcross
          int ia = b(i,j,k).nrange.ia() + map.gc[level].ia();
          int ib = b(i,j,k).nrange.ib() + map.gc[level].ib();
          int ja = b(i,j,k).nrange.ja() + map.gc[level].ja();
          int jb = b(i,j,k).nrange.jb() + map.gc[level].jb();
          int ka = b(i,j,k).nrange.ka() + map.gc[level].ka();
          int kb = b(i,j,k).nrange.kb() + map.gc[level].kb();
          msm::Ivec na = map.clipIndexToLevel(msm::Ivec(ia,ja,ka), level);
          msm::Ivec nb = map.clipIndexToLevel(msm::Ivec(ib,jb,kb), level);
          b(i,j,k).nrangeCutoff.setbounds(na.i, nb.i, na.j, nb.j, na.k, nb.k);
          // determine sendAcross blocks
          msm::BlockIndex blower = map.blockOfGridIndex(na, level);
          msm::BlockIndex bupper = map.blockOfGridIndex(nb, level);
          int maxarrlen = (bupper.n.i - blower.n.i + 1) *
            (bupper.n.j - blower.n.j + 1) * (bupper.n.k - blower.n.k + 1);
          b(i,j,k).sendAcross.setmax(maxarrlen);  // allocate send array
          b(i,j,k).indexGridCutoff.setmax(maxarrlen);  // alloc indexing
          // loop over sendAcross blocks
          int ii, jj, kk;
          for (kk = blower.n.k;  kk <= bupper.n.k;  kk++) {
            for (jj = blower.n.j;  jj <= bupper.n.j;  jj++) {
              for (ii = blower.n.i;  ii <= bupper.n.i;  ii++) {
#ifdef MSM_SKIP_TOO_DISTANT_BLOCKS
                // make sure that block (ii,jj,kk) interacts with (i,j,k)
                int si = sign(ii-i);
                int sj = sign(jj-j);
                int sk = sign(kk-k);
                int di = (ii-i)*bsx + si*(1-bsx);
                int dj = (jj-j)*bsy + sj*(1-bsy);
                int dk = (kk-k)*bsz + sk*(1-bsz);
                Vector d = di*hu + dj*hv + dk*hw;
                if (d.length2() >= 4*a*a) continue;
#endif
                // determine actual block and range to send to
                msm::BlockSend bs;
                bs.nblock.n = msm::Ivec(ii,jj,kk);
                bs.nblock.level = level;
                bs.nrange = map.clipBlockToIndexRange(bs.nblock,
                    b(i,j,k).nrangeCutoff);
                map.wrapBlockSend(bs);  // wrap to true block index
                b(i,j,k).sendAcross.append(bs);
                b(i,j,k).indexGridCutoff.append(numGridCutoff);
                numGridCutoff++;  // one MsmGridCutoff for each send across
                // increment counter for receive block
                b(bs.nblock_wrap.n).numRecvsPotential++;
              }
            }
          } // end loop over sendAcross blocks

          // Restriction, sendUp
          if (level < nlevels-1) {
            int ia2, ib2, ja2, jb2, ka2, kb2;
            ia = b(i,j,k).nrange.ia();
            ib = b(i,j,k).nrange.ib();
            ja = b(i,j,k).nrange.ja();
            jb = b(i,j,k).nrange.jb();
            ka = b(i,j,k).nrange.ka();
            kb = b(i,j,k).nrange.kb();
            // determine expansion of h-grid onto 2h-grid
            if ( ia==ib && ((ia & 1)==0) ) {
              ia2 = ib2 = ia / 2;
            }
            else {
              ia2 = (ia / 2) - ((polydeg+1) / 2) + 1;
              ib2 = ((ib+1) / 2) + ((polydeg+1) / 2) - 1;
            }
            if ( ja==jb && ((ja & 1)==0) ) {
              ja2 = jb2 = ja / 2;
            }
            else {
              ja2 = (ja / 2) - ((polydeg+1) / 2) + 1;
              jb2 = ((jb+1) / 2) + ((polydeg+1) / 2) - 1;
            }
            if ( ka==kb && ((ka & 1)==0) ) {
              ka2 = kb2 = ka / 2;
            }
            else {
              ka2 = (ka / 2) - ((polydeg+1) / 2) + 1;
              kb2 = ((kb+1) / 2) + ((polydeg+1) / 2) - 1;
            }
            // clip to boundaries of 2h-grid
            msm::Ivec na2, nb2;
            na2 = map.clipIndexToLevel(msm::Ivec(ia2,ja2,ka2), level+1);
            nb2 = map.clipIndexToLevel(msm::Ivec(ib2,jb2,kb2), level+1);
            b(i,j,k).nrangeRestricted.setbounds(na2.i, nb2.i, na2.j, nb2.j,
                na2.k, nb2.k);
            // determine sendUp blocks
            msm::BlockIndex blower = map.blockOfGridIndex(na2, level+1);
            msm::BlockIndex bupper = map.blockOfGridIndex(nb2, level+1);
            int maxarrlen = (bupper.n.i - blower.n.i + 1) *
              (bupper.n.j - blower.n.j + 1) * (bupper.n.k - blower.n.k + 1);
            b(i,j,k).sendUp.setmax(maxarrlen);  // allocate send array
            // loop over sendUp blocks
            int ii, jj, kk;
            for (kk = blower.n.k;  kk <= bupper.n.k;  kk++) {
              for (jj = blower.n.j;  jj <= bupper.n.j;  jj++) {
                for (ii = blower.n.i;  ii <= bupper.n.i;  ii++) {
                  // determine actual block and range to send to
                  msm::BlockSend bs;
                  bs.nblock.n = msm::Ivec(ii,jj,kk);
                  bs.nblock.level = level+1;
                  bs.nrange = map.clipBlockToIndexRange(bs.nblock,
                      b(i,j,k).nrangeRestricted);
                  map.wrapBlockSend(bs);  // wrap to true block index
                  b(i,j,k).sendUp.append(bs);
                  // increment counter for receive block
                  map.blockLevel[level+1](bs.nblock_wrap.n).numRecvsCharge++;
                }
              }
            } // end loop over sendUp blocks

          } // end if restriction

          // Prolongation, sendDown
          if (level > 0) {
            int ia2 = b(i,j,k).nrange.ia();
            int ib2 = b(i,j,k).nrange.ib();
            int ja2 = b(i,j,k).nrange.ja();
            int jb2 = b(i,j,k).nrange.jb();
            int ka2 = b(i,j,k).nrange.ka();
            int kb2 = b(i,j,k).nrange.kb();
            // determine expansion of 2h-grid onto h-grid
            ia = 2*ia2 - polydeg;
            ib = 2*ib2 + polydeg;
            ja = 2*ja2 - polydeg;
            jb = 2*jb2 + polydeg;
            ka = 2*ka2 - polydeg;
            kb = 2*kb2 + polydeg;
            // clip to boundaries of h-grid
            msm::Ivec na, nb;
            na = map.clipIndexToLevel(msm::Ivec(ia,ja,ka), level-1);
            nb = map.clipIndexToLevel(msm::Ivec(ib,jb,kb), level-1);
            b(i,j,k).nrangeProlongated.setbounds(na.i, nb.i, na.j, nb.j,
                na.k, nb.k);
            // determine sendDown blocks
            msm::BlockIndex blower = map.blockOfGridIndex(na, level-1);
            msm::BlockIndex bupper = map.blockOfGridIndex(nb, level-1);
            int maxarrlen = (bupper.n.i - blower.n.i + 1) *
              (bupper.n.j - blower.n.j + 1) * (bupper.n.k - blower.n.k + 1);
            b(i,j,k).sendDown.setmax(maxarrlen);  // allocate send array
            // loop over sendUp blocks
            int ii, jj, kk;
            for (kk = blower.n.k;  kk <= bupper.n.k;  kk++) {
              for (jj = blower.n.j;  jj <= bupper.n.j;  jj++) {
                for (ii = blower.n.i;  ii <= bupper.n.i;  ii++) {
                  // determine actual block and range to send to
                  msm::BlockSend bs;
                  bs.nblock.n = msm::Ivec(ii,jj,kk);
                  bs.nblock.level = level-1;
                  bs.nrange = map.clipBlockToIndexRange(bs.nblock,
                      b(i,j,k).nrangeProlongated);
                  map.wrapBlockSend(bs);  // wrap to true block index
                  b(i,j,k).sendDown.append(bs);
                  // increment counter for receive block
                  map.blockLevel[level-1](bs.nblock_wrap.n).numRecvsPotential++;
                }
              }
            } // end loop over sendDown blocks

          } // end if prolongation

        }
      }
    }
  }
  // end of Map setup

  // allocate chare arrays

  if (1) {
    PatchMap *pm = PatchMap::Object();
    patchPtr.resize( pm->numPatches() );
    for (int i = 0;  i < pm->numPatches();  i++) {
      patchPtr[i] = NULL;
    }
#ifdef DEBUG_MSM_VERBOSE
    printf("Allocating patchPtr array length %d\n", pm->numPatches());
#endif
    if (CkMyPe() == 0) {
      iout << iINFO << "MSM has " << pm->numPatches()
                    << " interpolation / anterpolation objects"
                    << " (one per patch)\n";
    }
  }

#ifdef MSM_NODE_MAPPING
  if (1) {
    // Node aware initial assignment of chares
    //
    // Create map object for each 3D chare array of MsmBlock and the
    // 1D chare array of MsmGridCutoff.  Design map to equally distribute
    // blocks across nodes, assigned to node PEs in round robin manner.  
    // Attempt to reduce internode communication bandwidth by assigning 
    // each MsmGridCutoff element to either its source node or its 
    // destination node, again assigned to node PEs in round robin manner.
#if 0
    int numNodes = 32;
    int numPes = 1024;
#else
    int numNodes = CkNumNodes();
    int numPes = CkNumPes();
#endif
    int numBlocks = 0;  // find total number of blocks
    for (level = 0;  level < nlevels;  level++) {
      numBlocks += map.blockLevel[level].nn();
    }
    blockAssign.resize(numBlocks);
    gcutAssign.resize(numGridCutoff);
    nodecnt.resize(numNodes);

    // assign blocks to nodes
    // the following algorithm divides as evenly as possible the 
    // blocks across the nodes
    int r = numBlocks % numNodes;
    if (r == 0) {
      int q = numBlocks / numNodes;
      for (n = 0;  n < numNodes;  n++) {
        int moffset = n * q;
        for (int m = 0;  m < q;  m++) {
          blockAssign[moffset + m] = n;
        }
        nodecnt[n] = q;
      }
#if 0
      if (CkMyPe() == 0) {
        CkPrintf("%d objects each assigned to %d nodes\n", q, numNodes);
        CkPrintf("%d  =?  %d\n", numNodes * q, numBlocks);
      }
#endif
    }
    else {  // numNodes does not evenly divide numBlocks
      r = numBlocks % (numNodes-1);
      int q = numBlocks / (numNodes-1);
      int qp, ncnt;
      if (r <= q) {
        ncnt = q - r;  // # of procs to assign q-1 objects
        qp = q - 1;
      }
      else {
        ncnt = r - q;  // # of procs to assign q+1 objects
        qp = q + 1;
      }
      // (numNodes - ncnt) gets q objs, ncnt gets qp objects
#if 0
      if (CkMyPe() == 0) {
        CkPrintf("%d objects to %d nodes\n", q, numNodes-ncnt);
        CkPrintf("%d objects to %d nodes\n", qp, ncnt);
        CkPrintf("%d  =?  %d\n", (numNodes-ncnt)*q + ncnt*qp, numBlocks);
      }
#endif
      for (n = 0;  n < numNodes - ncnt;  n++) {
        int moffset = n * q;
        for (int m = 0;  m < q;  m++) {
          blockAssign[moffset + m] = n;
        }
        nodecnt[n] = q;
      }
      for ( ;  n < numNodes;  n++) {
        int moffset = (numNodes - ncnt)*q + (n - (numNodes - ncnt))*qp;
        for (int m = 0;  m < qp;  m++) {
          blockAssign[moffset + m] = n;
        }
        nodecnt[n] = qp;
      }
    }

    // assign grid cutoff objects to nodes (gcutAssign)
    // choose whichever of source or destination node has less work
    int gindex = 0;  // index for grid cutoff array
    for (level = 0;  level < nlevels;  level++) { // for all levels
      msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
      int bni = b.ni();
      int bnj = b.nj();
      int bnk = b.nk();
      for (k = 0;  k < bnk;  k++) { // for all blocks
        for (j = 0;  j < bnj;  j++) {
          for (i = 0;  i < bni;  i++) {
            int isrc = blockFlatIndex(level, i, j, k);
            int nsrc = blockAssign[isrc];  // node block isrc is assigned
            int numSendAcross = b(i,j,k).sendAcross.len();
            ASSERT( numSendAcross == b(i,j,k).indexGridCutoff.len() );
            for (n = 0;  n < numSendAcross;  n++) {
              msm::BlockIndex &bn = b(i,j,k).sendAcross[n].nblock_wrap;
              int idest = blockFlatIndex(level, bn.n.i, bn.n.j, bn.n.k);
              int ndest = blockAssign[idest];  // node block idest is assigned
              // assign this grid cutoff work to least subscribed node
              if (nodecnt[nsrc] <= nodecnt[ndest]) {
                gcutAssign[gindex] = nsrc;
                nodecnt[nsrc]++;
              }
              else {
                gcutAssign[gindex] = ndest;
                nodecnt[ndest]++;
              }
              gindex++;
            } // end for numSendAcross
          }
        }
      } // end for all blocks
    } // end for all levels

    // now change the node assignments into PE assignments
    // use round robin assignment to PEs within each node
    int ppn = numPes / numNodes;  // num PEs per node
    // reset nodecnt - this array will now store PE offset for that node
    for (n = 0;  n < numNodes;  n++)  nodecnt[n] = 0;
    for (n = 0;  n < numBlocks;  n++) {
      int node = blockAssign[n];
      blockAssign[n] = node * ppn + nodecnt[node];  // PE number within node
      nodecnt[node]++;  // increment to next PE
      if (nodecnt[node] >= ppn)  nodecnt[node] = 0;  // with wrap around
    }
    for (n = 0;  n < numGridCutoff;  n++) {
      int node = gcutAssign[n];
      gcutAssign[n] = node * ppn + nodecnt[node];  // PE number within node
      nodecnt[node]++;  // increment to next PE
      if (nodecnt[node] >= ppn)  nodecnt[node] = 0;  // with wrap around
    }

    // print mapping
#if 0
    if (CkMyPe() == 0) {
      for (n = 0;  n < numBlocks;  n++) {
        CkPrintf("block %d:   node=%d  pe=%d\n",
            n, blockAssign[n]/ppn, blockAssign[n]);
      }
      for (n = 0;  n < numGridCutoff;  n++) {
        CkPrintf("grid cutoff %d:   node=%d  pe=%d\n",
            n, gcutAssign[n]/ppn, gcutAssign[n]);
      }
    }
#endif

  } // end node aware initial assignment of chares
#endif // MSM_NODE_MAPPING

  if (CkMyPe() == 0) {
    // on PE 0, create 3D chare array of MsmBlock for each level;
    // broadcast this array of proxies to the rest of the group
    msmBlock.resize(nlevels);
    for (level = 0;  level < nlevels;  level++) {
      int ni = map.blockLevel[level].ni();
      int nj = map.blockLevel[level].nj();
      int nk = map.blockLevel[level].nk();
      msmBlock[level] = CProxy_MsmBlock::ckNew(level, ni, nj, nk);
#ifdef DEBUG_MSM_VERBOSE
      printf("Create MsmBlock[%d] 3D chare array ( %d x %d x %d )\n",
          level, ni, nj, nk);
#endif
      char msg[128];
      int nijk = ni * nj * nk;
      sprintf(msg, "MSM grid level %d decomposed into %d block%s"
          " ( %d x %d x %d )\n",
          level, nijk, (nijk==1 ? "" : "s"), ni, nj, nk);
      iout << iINFO << msg;
    }
    MsmBlockProxyMsg *msg =
      new(nlevels*sizeof(CProxy_MsmBlock), 0) MsmBlockProxyMsg;
    msg->put(msmBlock);
    msmProxy.recvMsmBlockProxy(msg);  // broadcast

#ifdef MSM_GRID_CUTOFF_DECOMP
    // on PE 0, create 1D chare array of MsmGridCutoff
    // broadcast this array proxy to the rest of the group
    msmGridCutoff = CProxy_MsmGridCutoff::ckNew(numGridCutoff);
    MsmGridCutoffProxyMsg *gcmsg =
      new(sizeof(CProxy_MsmGridCutoff), 0) MsmGridCutoffProxyMsg;
    gcmsg->put(&msmGridCutoff);
    msmProxy.recvMsmGridCutoffProxy(gcmsg);

    // XXX PE 0 initializes each MsmGridCutoff
    // one-to-many
    // for M length chare array, better for each PE to initialize M/P?
    for (level = 0;  level < nlevels;  level++) { // for all levels
      msm::Grid<msm::BlockDiagram>& b = map.blockLevel[level];
      int bni = b.ni();
      int bnj = b.nj();
      int bnk = b.nk();
      for (k = 0;  k < bnk;  k++) { // for all blocks
        for (j = 0;  j < bnj;  j++) {
          for (i = 0;  i < bni;  i++) {
            int numSendAcross = b(i,j,k).sendAcross.len();
            ASSERT( numSendAcross == b(i,j,k).indexGridCutoff.len() );
            for (n = 0;  n < numSendAcross;  n++) {
              msm::BlockIndex bi = msm::BlockIndex(level, msm::Ivec(i,j,k));
              msm::BlockSend &bs = b(i,j,k).sendAcross[n];
              int index = b(i,j,k).indexGridCutoff[n];
              MsmGridCutoffInitMsg *bsmsg = new MsmGridCutoffInitMsg(bi, bs);
              msmGridCutoff[index].initialize(bsmsg);
            } // traverse sendAcross, indexGridCutoff arrays

          }
        }
      } // end for all blocks

    } // end for all levels

    iout << iINFO << "MSM grid cutoff calculation decomposed into "
      << numGridCutoff << " work objects\n";
#endif
    iout << endi;
  }

#ifdef DEBUG_MSM_VERBOSE
  printf("end of initialization\n");
#endif
}

void ComputeMsmMgr::recvMsmBlockProxy(MsmBlockProxyMsg *msg)
{
  msg->get(msmBlock);
  delete(msg);
}

void ComputeMsmMgr::recvMsmGridCutoffProxy(MsmGridCutoffProxyMsg *msg)
{
  msg->get(&msmGridCutoff);
  delete(msg);
}

void ComputeMsmMgr::update(CkQdMsg *msg)
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  update() PE %d\n", CkMyPe());
#endif
  delete msg;

  // XXX how do update for constant pressure simulation?
}


void ComputeMsmMgr::compute(msm::Array<int>& patchIDList)
{
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsmMgr:  compute() PE=%d\n", CkMyPe());
#endif

  int n; 
  for (n = 0;  n < patchIDList.len();  n++) {
    int patchID = patchIDList[n];
    if (patchPtr[patchID] == NULL) {
      char msg[100];
      snprintf(msg, sizeof(msg),
          "Expected MSM data for patch %d does not exist on PE %d",
          patchID, CkMyPe());
      NAMD_die(msg);
    }
    patchPtr[patchID]->anterpolation();
    // all else should follow from here
  }
  return;
}


void ComputeMsmMgr::addPotential(GridFloatMsg *gm)
{
  //msm::Grid<BigReal> subgrid;
  int pid;
  gm->get(subgrid, pid);
  delete gm;
  if (patchPtr[pid] == NULL) {
    char msg[100];
    snprintf(msg, sizeof(msg), "Expecting patch %d to exist on PE %d",
        pid, CkMyPe());
    NAMD_die(msg);
  }
  patchPtr[pid]->addPotential(subgrid);
}


void ComputeMsmMgr::doneCompute()
{
  msmCompute->saveResults();
}


//////////////////////////////////////////////////////////////////////////////
//
//  ComputeMsm
//  MSM compute objects, starts and finishes calculation;
//  there is up to one compute object per PE
//

ComputeMsm::ComputeMsm(ComputeID c) : ComputeHomePatches(c)
{
  CProxy_ComputeMsmMgr::ckLocalBranch(
      CkpvAccess(BOCclass_group).computeMsmMgr)->setCompute(this);
  SimParameters *simParams = Node::Object()->simParameters;
  qscaling = sqrtf(COULOMB / simParams->dielectric);
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsm:  (constructor) PE=%d\n", CkMyPe());
#endif
}

ComputeMsm::~ComputeMsm()
{
  // free memory
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsm:  (destructor) PE=%d\n", CkMyPe());
#endif
}

void ComputeMsm::doWork()
{
#if 0
#ifdef MSM_TIMING
  myMgr->initTiming();
#endif
#ifdef MSM_PROFILING
  myMgr->initProfiling();
#endif
#endif

  // patchList is inherited from ComputeHomePatches
  ResizeArrayIter<PatchElem> ap(patchList);
  numLocalPatches = patchList.size();
  cntLocalPatches = 0;
  ASSERT(cntLocalPatches < numLocalPatches);

  // for each patch do stuff
#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsm:  doWork() PE=%d\n", CkMyPe());
#endif

#ifdef DEBUG_MSM_VERBOSE
  printf("patchList size = %d\n", patchList.size() );
#endif

  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics ) {
    for (ap = ap.begin();  ap != ap.end();  ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&r);
    }
    reduction->submit();
    return;
  }
  msm::Map& map = myMgr->mapData();
  // This is the patchPtr array for MSM; any local patch will be set up
  // with a non-NULL pointer to its supporting data structure.
  msm::PatchPtrArray& patchPtr = myMgr->patchPtrArray();
  // also store just a list of IDs for the local patches
  msm::Array<int> patchIDList(numLocalPatches);
  patchIDList.resize(0);  // to use append on pre-allocated array buffer
  int cnt=0, n;
  for (ap = ap.begin();  ap != ap.end();  ap++) {
    CompAtom *x = (*ap).positionBox->open();
    CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();
    int patchID = (*ap).patchID;
    patchIDList.append(patchID);
    if (patchPtr[patchID] == NULL) {
      // create PatchData if it doesn't exist for this patchID
      patchPtr[patchID] = new msm::PatchData(myMgr, patchID);
#ifdef DEBUG_MSM_VERBOSE
      printf("Creating new PatchData:  patchID=%d  PE=%d\n",
          patchID, CkMyPe());
#endif
    }
    msm::PatchData& patch = *(patchPtr[patchID]);
    patch.init(numAtoms);
    msm::AtomCoordArray& coord = patch.coordArray();
    ASSERT(coord.len() == numAtoms);
    for (n = 0;  n < numAtoms;  n++) {
      coord[n].position = x[n].position;
      coord[n].charge = qscaling * x[n].charge;
      coord[n].id = xExt[n].id;
    }
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).avgPositionBox->close(&x);
    }
    else {
      (*ap).positionBox->close(&x);
    }
  }

  myMgr->compute(patchIDList);
}

void ComputeMsm::saveResults()
{
  if (++cntLocalPatches != numLocalPatches) return;

  // NAMD patches
  ResizeArrayIter<PatchElem> ap(patchList);
#ifdef DEBUG_MSM
  for (ap = ap.begin();  ap != ap.end();  ap++) {
    int patchID = (*ap).patchID;
    ASSERT(myMgr->patchPtrArray()[patchID]->cntRecvs ==
        myMgr->mapData().patchList[patchID].numRecvs);
  }
#endif

  // get results from ComputeMsmMgr
  msm::PatchPtrArray& patchPtr = myMgr->patchPtrArray();

#ifdef DEBUG_MSM_VERBOSE
  printf("ComputeMsm:  saveResults() PE=%d\n", CkMyPe());
#endif
  // store force updates
  // submit reductions

  // add in forces
  int cnt=0, n;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();
    int patchID = (*ap).patchID;
    if (patchPtr[patchID] == NULL) {
      char msg[100];
      snprintf(msg, sizeof(msg), "Expecting patch %d to exist on PE %d",
          patchID, CkMyPe());
      NAMD_die(msg);
    }
    msm::PatchData& patch = *(patchPtr[patchID]);
    ASSERT(numAtoms == patch.force.len() );
    for (n = 0;  n < numAtoms;  n++) {
      f[n] += patch.force[n];
    }
    (*ap).forceBox->close(&r);

    reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += patch.energy;
    reduction->item(REDUCTION_VIRIAL_SLOW_XX) += patch.virial[0][0];
    reduction->item(REDUCTION_VIRIAL_SLOW_XY) += patch.virial[0][1];
    reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += patch.virial[0][2];
    reduction->item(REDUCTION_VIRIAL_SLOW_YX) += patch.virial[1][0];
    reduction->item(REDUCTION_VIRIAL_SLOW_YY) += patch.virial[1][1];
    reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += patch.virial[1][2];
    reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += patch.virial[2][0];
    reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += patch.virial[2][1];
    reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += patch.virial[2][2];
  }
  reduction->submit();
}

// method definitions for PatchData
namespace msm {

  PatchData::PatchData(ComputeMsmMgr *pmgr, int pid) {
    mgr = pmgr;
    map = &(mgr->mapData());
    patchID = pid;
    //PatchMap *pm = PatchMap::Object();
    pd = &(map->patchList[pid]);
    qh.init(pd->nrange);
    eh.init(pd->nrange);
    subgrid.resize(map->bsx[0] * map->bsy[0] * map->bsz[0]);
#ifdef MSM_TIMING
    mgr->addTiming();
#endif
  }

  void PatchData::init(int natoms) {
    coord.resize(natoms);
    force.resize(natoms);
    cntRecvs = 0;
    energy = 0;
    memset(virial, 0, 3*3*sizeof(BigReal));
    for (int i = 0;  i < natoms;  i++)  force[i] = 0;
    qh.reset(0);
    eh.reset(0);
  }

  void PatchData::anterpolation() {
#ifdef DEBUG_MSM_GRID
    printf("patchID %d:  anterpolation\n", patchID);
#endif

#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
    Float xphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float yphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float zphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];

    const Double rs_edge = Double( mgr->s_edge );
    const int s_size = ComputeMsmMgr::PolyDegree[mgr->approx] + 1;

    const int ia = qh.ia();
    const int ib = qh.ib();
    const int ja = qh.ja();
    const int jb = qh.jb();
    const int ka = qh.ka();
    const int kb = qh.kb();
    const int ni = qh.ni();
    const int nj = qh.nj();
    Float *qhbuffer = qh.data().buffer();

    // loop over atoms
    for (int n = 0;  n < coord.len();  n++) {
      Float q = coord[n].charge;
      if (0==q) continue;

      ScaledPosition s = mgr->lattice.scale(coord[n].position);

      BigReal sx_hx = (s.x - mgr->sglower.x) * mgr->shx_1;
      BigReal sy_hy = (s.y - mgr->sglower.y) * mgr->shy_1;
      BigReal sz_hz = (s.z - mgr->sglower.z) * mgr->shz_1;

      BigReal xlo = floor(sx_hx) - rs_edge;
      BigReal ylo = floor(sy_hy) - rs_edge;
      BigReal zlo = floor(sz_hz) - rs_edge;

      // calculate Phi stencils along each dimension
      Float xdelta = Float(sx_hx - xlo);
      mgr->stencil_1d(xphi, xdelta);
      Float ydelta = Float(sy_hy - ylo);
      mgr->stencil_1d(yphi, ydelta);
      Float zdelta = Float(sz_hz - zlo);
      mgr->stencil_1d(zphi, zdelta);

      int ilo = int(xlo);
      int jlo = int(ylo);
      int klo = int(zlo);

      // test to see if stencil is within edges of grid
      int iswithin = ( ia <= ilo && (ilo+(s_size-1)) <= ib &&
                       ja <= jlo && (jlo+(s_size-1)) <= jb &&
                       ka <= klo && (klo+(s_size-1)) <= kb );

      if ( ! iswithin ) {
        char msg[100];
        snprintf(msg, sizeof(msg), "Atom %d is outside of the MSM grid.",
            coord[n].id);
        NAMD_die(msg);
      }

      // determine charge on cube of grid points around atom
      for (int k = 0;  k < s_size;  k++) {
        int koff = ((k+klo) - ka) * nj;
        Float ck = zphi[k] * q;
        for (int j = 0;  j < s_size;  j++) {
          int jkoff = (koff + (j+jlo) - ja) * ni;
          Float cjk = yphi[j] * ck;
          for (int i = 0;  i < s_size;  i++) {
            int ijkoff = jkoff + (i+ilo) - ia;
            qhbuffer[ijkoff] += xphi[i] * cjk;
          }
        }
      }

    } // end loop over atoms
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::ANTERP] += stopTime - startTime;
#endif

    sendCharge();
  }

  void PatchData::sendCharge() {
#ifdef MSM_TIMING
    double startTime, stopTime;
#endif
    int priority = 1;
    // buffer portions of grid to send to Blocks on level 0
    // allocate the largest buffer space we'll need
    //Grid<BigReal> subgrid;
    //subgrid.resize(map->bsx[0] * map->bsy[0] * map->bsz[0]);
    for (int n = 0;  n < pd->send.len();  n++) {
#ifdef MSM_TIMING
      startTime = CkWallTimer();
#endif
      // initialize the proper subgrid indexing range
      subgrid.init( pd->send[n].nrange );
      // extract the values from the larger grid into the subgrid
      qh.extract(subgrid);
      // translate the subgrid indexing range to match the MSM block
      subgrid.updateLower( pd->send[n].nrange_wrap.lower() );
      // add the subgrid charges into the block
      BlockIndex& bindex = pd->send[n].nblock_wrap;
      // place subgrid into message
      int nelems = subgrid.data().len();
#ifdef MSM_FIXED_SIZE_GRID_MSG
      GridFloatMsg *gm = new(sizeof(int)) GridFloatMsg;
#else
      GridFloatMsg *gm = new(nelems, sizeof(int)) GridFloatMsg;
#endif
      *(int *)CkPriorityPtr(gm) = priority;
      CkSetQueueing(gm, CK_QUEUEING_IFIFO);
      gm->put(subgrid, bindex.level);
#ifdef MSM_TIMING
      stopTime = CkWallTimer();
      mgr->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
      mgr->msmBlock[bindex.level](
          bindex.n.i, bindex.n.j, bindex.n.k).addCharge(gm);
    }
  }

  void PatchData::addPotential(const Grid<Float>& epart) {
#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
    eh += epart;
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::COMM] += stopTime - startTime;
#endif
    if (++cntRecvs == pd->numRecvs) {
      interpolation();
    }
  }

  void PatchData::interpolation() {
#ifdef DEBUG_MSM_GRID
    printf("patchID %d:  interpolation\n", patchID);
#endif

#ifdef MSM_TIMING
    double startTime, stopTime;
    startTime = CkWallTimer();
#endif
    BigReal energy_self = 0;

    Float xphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float yphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float zphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float dxphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float dyphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];
    Float dzphi[ComputeMsmMgr::MAX_POLY_DEGREE+1];

    const Double rs_edge = Double( mgr->s_edge );
    const int s_size = ComputeMsmMgr::PolyDegree[mgr->approx] + 1;

    const int ia = eh.ia();
    const int ib = eh.ib();
    const int ja = eh.ja();
    const int jb = eh.jb();
    const int ka = eh.ka();
    const int kb = eh.kb();
    const int ni = eh.ni();
    const int nj = eh.nj();
    Float *ehbuffer = eh.data().buffer();

    // loop over atoms
    for (int n = 0;  n < coord.len();  n++) {
      Float q = coord[n].charge;
      if (0==q) continue;

      ScaledPosition s = mgr->lattice.scale(coord[n].position);

      BigReal sx_hx = (s.x - mgr->sglower.x) * mgr->shx_1;
      BigReal sy_hy = (s.y - mgr->sglower.y) * mgr->shy_1;
      BigReal sz_hz = (s.z - mgr->sglower.z) * mgr->shz_1;

      BigReal xlo = floor(sx_hx) - rs_edge;
      BigReal ylo = floor(sy_hy) - rs_edge;
      BigReal zlo = floor(sz_hz) - rs_edge;

      // calculate Phi stencils along each dimension
      Float xdelta = Float(sx_hx - xlo);
      mgr->d_stencil_1d(dxphi, xphi, xdelta);
      Float ydelta = Float(sy_hy - ylo);
      mgr->d_stencil_1d(dyphi, yphi, ydelta);
      Float zdelta = Float(sz_hz - zlo);
      mgr->d_stencil_1d(dzphi, zphi, zdelta);

      int ilo = int(xlo);
      int jlo = int(ylo);
      int klo = int(zlo);

#if 0
      // XXX don't need to test twice!

      // test to see if stencil is within edges of grid
      int iswithin = ( ia <= ilo && (ilo+(s_size-1)) <= ib &&
                       ja <= jlo && (jlo+(s_size-1)) <= jb &&
                       ka <= klo && (klo+(s_size-1)) <= kb );

      if ( ! iswithin ) {
        char msg[100];
        snprintf(msg, sizeof(msg), "Atom %d is outside of the MSM grid.",
            coord[n].id);
        NAMD_die(msg);
      }
#endif

      // determine force on atom from surrounding potential grid points
      //Force f = 0;
      //BigReal e = 0;
      Float fx=0, fy=0, fz=0, e=0;
      for (int k = 0;  k < s_size;  k++) {
        int koff = ((k+klo) - ka) * nj;
        for (int j = 0;  j < s_size;  j++) {
          int jkoff = (koff + (j+jlo) - ja) * ni;
          Float cx = yphi[j] * zphi[k];
          Float cy = dyphi[j] * zphi[k];
          Float cz = yphi[j] * dzphi[k];
          for (int i = 0;  i < s_size;  i++) {
            int ijkoff = jkoff + (i+ilo) - ia;
            Float ec = ehbuffer[ijkoff];
            fx += ec * dxphi[i] * cx;
            fy += ec * xphi[i] * cy;
            fz += ec * xphi[i] * cz;
            e += ec * xphi[i] * cx;
          }
        }
      }

      force[n].x -= q * (mgr->srx_x * fx + mgr->srx_y * fy + mgr->srx_z * fz);
      force[n].y -= q * (mgr->sry_x * fx + mgr->sry_y * fy + mgr->sry_z * fz);
      force[n].z -= q * (mgr->srz_x * fx + mgr->srz_y * fy + mgr->srz_z * fz);
      energy += q * e;
      energy_self += q * q;

    } // end loop over atoms

    energy_self *= mgr->gzero;
    energy -= energy_self;
    energy *= 0.5;
#ifdef MSM_TIMING
    stopTime = CkWallTimer();
    mgr->msmTiming[MsmTimer::INTERP] += stopTime - startTime;
#endif

#ifdef MSM_TIMING
    mgr->doneTiming();
#endif
    mgr->doneCompute();
  }

} // namespace msm


#include "ComputeMsmMgr.def.h"
