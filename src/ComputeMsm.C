
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
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "varsizemsg.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "pup_stl.h"
#include "MsmMacros.h"

void MsmData::pup(PUP::er &p)
{
  p|ispx, p|ispy, p|ispz;
}

void MsmData::print()
{
  // print something here
}


//////////////////////////////////////////////////////////////////////////////
//
//  ComputeMsmMgr
//  chare group containing MSM parameters and constants;
//  one chare object per PE
//

class ComputeMsmMgr : public BOCclass {
public:
  ComputeMsmMgr();                    // entry
  ~ComputeMsmMgr();

  void initialize(MsmInitMsg *);      // entry with message
  void update(CkQdMsg *);             // entry with message

  void setCompute(ComputeMsm *c) { msmCompute = c; }  // local

private:
  int setup_hgrid_1d(double len, double& hh, int& nn,
      int& ia, int& ib, int isperiodic);

  CProxy_ComputeMsmMgr msmProxy;
  ComputeMsm *msmCompute;

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
  Vector rlower;        // lower corner of grid in unit space
  Vector rh;            // grid spacings in unit space
  Vector rh_1;          // inverse of unit space grid spacings
  Vector rx_rhx;        // row vector to transform interpolated force x
  Vector ry_rhy;        // row vector to transform interpolated force y
  Vector rz_rhz;        // row vector to transform interpolated force z

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
  static const BigReal PhiStencil[NUM_APPROX][MAX_NSTENCIL_SKIP_ZERO];

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
const BigReal ComputeMsmMgr::PhiStencil[NUM_APPROX][MAX_NSTENCIL_SKIP_ZERO] = {
  // cubic
  {-1./16, 9./16, 1, 9./16, -1./16},

  // quintic C1
  {3./256, -25./256, 75./128, 1, 75./128, -25./256, 3./256},

  // quintic C2  (same as quintic C1)
  {3./256, -25./256, 75./128, 1, 75./128, -25./256, 3./256},

  // septic C1
  { -5./2048, 49./2048, -245./2048, 1225./2048, 1, 1225./2048,
    -245./2048, 49./2048, -5./2048 },

  // septic C3  (same as septic C3)
  { -5./2048, 49./2048, -245./2048, 1225./2048, 1, 1225./2048,
    -245./2048, 49./2048, -5./2048 },

  // nonic C1
  { 35./65536, -405./65536, 567./16384, -2205./16384, 
    19845./32768, 1, 19845./32768, -2205./16384, 567./16384, 
    -405./65536, 35./65536 },

  // nonic C4  (same as nonic C1)
  { 35./65536, -405./65536, 567./16384, -2205./16384, 
    19845./32768, 1, 19845./32768, -2205./16384, 567./16384, 
    -405./65536, 35./65536 },
};


ComputeMsmMgr::ComputeMsmMgr() :
  msmProxy(thisgroup), msmCompute(0)
{
  //printf("ComputeMsmMgr:  (constructor) PE %d\n", CkMyPe());
  CkpvAccess(BOCclass_group).computeMsmMgr = thisgroup;
}

ComputeMsmMgr::~ComputeMsmMgr()
{
  //printf("ComputeMsmMgr:  (destructor) PE %d\n", CkMyPe());
  // free memory
}

int ComputeMsmMgr::setup_hgrid_1d(double len, double& hh, int& nn,
    int& ia, int& ib, int isperiodic)
{
  if (isperiodic) {
    const double hmin = (4./5) * gridspacing;
    const double hmax = 1.5 * hmin;
    hh = len;
    nn = 1;  // start with one grid point across length
    while (hh >= hmax) {
      hh *= 0.5;  // halve spacing and double grid points
      nn <<= 1;
    }
    if (hh < hmin) {
      if (nn < 4) {
        return -1;   // either len is too small or gridspacing is too large
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
    nn = (int) ceil(len / hh);
    ia = -s_edge;
    ib = nn + s_edge;
  }
  return 0;
}

void ComputeMsmMgr::initialize(MsmInitMsg *msg)
{
  printf("ComputeMsmMgr:  initialize() PE %d\n", CkMyPe());

  smin = msg->smin;
  smax = msg->smax;
  delete msg;

//  if (CkMyPe() != 0) return;  // initialize() is called on all PEs
                              // but want only PE 0 to do the initialization

  SimParameters *simParams = Node::Object()->simParameters;

  // get required sim params
  lattice = simParams->lattice;
  a = simParams->cutoff;
  gridspacing = simParams->MSMGridSpacing;
  padding = simParams->MSMPadding;
  approx = simParams->MSMApprox;   // XXX needs to be more user friendly
  split = simParams->MSMSplit;     // XXX needs to be more user friendly
  nlevels = simParams->MSMLevels;
  dispersion = 0;                  // XXX unused for now
  rlower = Vector(-0.5);

  s_edge = (PolyDegree[approx] - 1) / 2;  // stencil edge size
  omega = 2 * PolyDegree[approx];         // smallest non-periodic grid length

  BigReal alen = lattice.a().length();
  BigReal blen = lattice.b().length();
  BigReal clen = lattice.c().length();

  int ispx = lattice.a_p();
  int ispy = lattice.b_p();
  int ispz = lattice.c_p();
  if ( ! ispx ) {
    alen = smax.x - smin.x + 2*padding;
    rlower.x = smin.x - padding;
  }
  if ( ! ispy ) {
    blen = smax.y - smin.y + 2*padding;
    rlower.y = smin.y - padding;
  }
  if ( ! ispz ) {
    clen = smax.z - smin.z + 2*padding;
    rlower.z = smin.z - padding;
  }

  int ia, ib, ja, jb, ka, kb;
  int rc = 0;

  if ((rc=setup_hgrid_1d(alen, h.x, nhx, ia, ib, ispx)) ||
      (rc=setup_hgrid_1d(blen, h.y, nhy, ja, jb, ispy)) ||
      (rc=setup_hgrid_1d(clen, h.z, nhz, ka, kb, ispz))) {
    NAMD_die("die");
  }

  // calculate self energy factor for splitting
  BigReal gs=0, d=0, s=0;
  splitting(gs, d, s, split);
  gzero = gs / a;  // XXX changes to a^6 for dispersion

  // determine grid hierarchy
  // msm map, broadcast
  // allocate chare arrays
}

void ComputeMsmMgr::update(CkQdMsg *msg)
{
  printf("ComputeMsmMgr:  update() PE %d\n", CkMyPe());
  delete msg;

//  if (CkMyPe() != 0) return;  // update() is called on all PEs
                              // but want only PE 0 to do the update

  // find grid spacing basis vectors in real space
  Vector hu = h.x * lattice.a().unit();
  Vector hv = h.y * lattice.b().unit();
  Vector hw = h.z * lattice.c().unit();

  Vector pu = lattice.scale(hu + lattice.origin());
  Vector pv = lattice.scale(hv + lattice.origin());
  Vector pw = lattice.scale(hw + lattice.origin());
  rh = Vector(pu.x, pv.y, pw.z);  // keep grid spacings in unit space
  rh_1 = Vector(1/pu.x, 1/pv.y, 1/pw.z);

  // row vectors to transform interpolated force back to real space
  Vector ru = lattice.a_r();
  Vector rv = lattice.b_r();
  Vector rw = lattice.c_r();
  rx_rhx = rh_1.x * Vector(ru.x, rv.x, rw.x);
  ry_rhy = rh_1.y * Vector(ru.y, rv.y, rw.y);
  rz_rhz = rh_1.z * Vector(ru.z, rv.z, rw.z);

  BigReal s;
  pu = cross(hv, hw);
  s = (hu * pu) / (pu * pu);
  pu *= s;  // pu is orthogonal projection of hu onto hv CROSS hw

  pv = cross(hw, hu);
  s = (hv * pv) / (pv * pv);
  pv *= s;  // pv is orthogonal projection of hv onto hw CROSS hu

  pw = cross(hu, hv);
  s = (hw * pw) / (pw * pw);
  pw *= s;  // pw is orthogonal projection of hw onto hu CROSS hv

  // radii for parallelepiped of weights enclosing grid cutoff sphere
  int ni = (int) ceil(2*a / pu.length()) - 1;
  int nj = (int) ceil(2*a / pv.length()) - 1;
  int nk = (int) ceil(2*a / pw.length()) - 1;


  // calculate constants for grid cutoff, broadcast
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
  qscaling = sqrt(COULOMB / simParams->dielectric);
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  printf("ComputeMsm:  (constructor) PE=%d\n", CkMyPe());
}

ComputeMsm::~ComputeMsm()
{
  // free memory
  printf("ComputeMsm:  (destructor) PE=%d\n", CkMyPe());
}

void ComputeMsm::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);

  // for each patch do stuff
  printf("ComputeMsm:  doWork() PE=%d\n", CkMyPe());

  // Skip computations if nothing to do.
  if ( 1 /* ! patchList[0].p->flags.doFullElectrostatics */ ) {
    for (ap = ap.begin();  ap != ap.end();  ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&r);
      reduction->submit();
    }
    return;
  }
}

void ComputeMsm::saveResults(int n, const MsmForce force[], double self_energy)
{
  printf("ComputeMsm:  saveResults() PE=%d\n", CkMyPe());
  // store force updates
  // submit reductions
  BigReal sum = 0;
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += 0.5*sum;
  reduction->submit();
}


#include "ComputeMsmMgr.def.h"

