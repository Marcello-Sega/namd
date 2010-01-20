/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeAniso.h"
#include "Molecule.h"
#include "Parameters.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"

#define CALCULATE_ANISO

// static initialization
int AnisoElem::pressureProfileSlabs = 0;
int AnisoElem::pressureProfileAtomTypes = 1;
BigReal AnisoElem::pressureProfileThickness = 0;
BigReal AnisoElem::pressureProfileMin = 0;

void AnisoElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Aniso** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in AnisoElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numAnisos;
  *byatom = mol->anisosByAtom;
  *structarray = mol->anisos;
#endif
}

void AnisoElem::getParameterPointers(Parameters *p, const AnisoValue **v) {
  *v = NULL;  // parameters are stored in the structure
}

void AnisoElem::computeForce(BigReal *reduction, 
                                BigReal *pressureProfileData)
{
  DebugM(3, "::computeForce() localIndex = " << localIndex[0] << " "
               << localIndex[1] << " " << localIndex[2] << " "
               << localIndex[3] << std::endl);

#ifdef CALCULATE_ANISO
  // used some comments from Ed Harder's implementation in CHARMM

  const BigReal kpar0  = 2*value->k11;  // force constants
  const BigReal kperp0 = 2*value->k22;
  const BigReal kiso0  = 2*value->k33;

  const Position & ri = p[0]->x[localIndex[0]].position;    // atom I
  const Position & rj = p[0]->x[localIndex[0]+1].position;  // atom I's Drude
  const Position & rl = p[1]->x[localIndex[1]].position;    // atom L
  const Position & rm = p[2]->x[localIndex[2]].position;    // atom M
  const Position & rn = p[3]->x[localIndex[3]].position;    // atom N

  // calculate parallel and perpendicular displacement vectors
  const Lattice & lattice = p[0]->p->lattice;
  Vector u1 = lattice.delta(ri,rl);  // shortest vector image:  ri - rl
  Vector u2 = lattice.delta(rm,rn);  // shortest vector image:  rm - rn

  BigReal u1_invlen = u1.rlength();  // need reciprocal lengths of u1, u2
  BigReal u2_invlen = u2.rlength();

  u1 *= u1_invlen;  // normalize u1, u2
  u2 *= u2_invlen;

  Vector dr = rj - ri;  // Drude displacement vector

  BigReal dpar  = dr * u1;  // parallel displacement
  BigReal dperp = dr * u2;  // perpendicular displacement

  // aniso spring energy
  // kpar reduces response along carbonyl vector
  // kperp reduces response perp to bond vector
  //   (reg in and out of plane response)
  BigReal eaniso;
  eaniso = 0.5*kpar0*dpar*dpar + 0.5*kperp0*dperp*dperp + 0.5*kiso0*(dr*dr);

  // calculate force vectors in one direction only

  // force into j
  Vector fj = -kiso0 * dr;
  fj -= kpar0 * dpar * u1;
  fj -= kperp0 * dperp * u2;

  // force from l
  Vector fi_l = kpar0 * dpar * dpar * u1_invlen * u1;
  fi_l -= kpar0 * dpar * u1_invlen * dr;

  // force into m
  Vector fm = kperp0 * dperp * dperp * u2_invlen * u2;
  fm -= kperp0 * dperp * u2_invlen * dr;

  // accumulate forces
  p[0]->f[localIndex[0]] += fi_l - fj;
  p[0]->f[localIndex[0]+1] += fj;
  p[1]->f[localIndex[1]] -= fi_l;
  p[2]->f[localIndex[2]] += fm;
  p[3]->f[localIndex[3]] -= fm;

  // update potential
  reduction[anisoEnergyIndex] += eaniso;

  reduction[virialIndex_XX] += fj.x * dr.x + fi_l.x * u1.x + fm.x * u2.x;
  reduction[virialIndex_XY] += fj.x * dr.y + fi_l.x * u1.y + fm.x * u2.y;
  reduction[virialIndex_XZ] += fj.x * dr.z + fi_l.x * u1.z + fm.x * u2.z;
  reduction[virialIndex_YX] += fj.y * dr.x + fi_l.y * u1.x + fm.y * u2.x;
  reduction[virialIndex_YY] += fj.y * dr.y + fi_l.y * u1.y + fm.y * u2.y;
  reduction[virialIndex_YZ] += fj.y * dr.z + fi_l.y * u1.z + fm.y * u2.z;
  reduction[virialIndex_ZX] += fj.z * dr.x + fi_l.z * u1.x + fm.z * u2.x;
  reduction[virialIndex_ZY] += fj.z * dr.y + fi_l.z * u1.y + fm.z * u2.y;
  reduction[virialIndex_ZZ] += fj.z * dr.z + fi_l.z * u1.z + fm.z * u2.z;

  // update pressure profile data
  if (pressureProfileData) {
    BigReal zi = p[0]->x[localIndex[0]].position.z;
    BigReal zj = p[0]->x[localIndex[0]+1].position.z;
    BigReal zl = p[1]->x[localIndex[1]].position.z;
    BigReal zm = p[2]->x[localIndex[2]].position.z;
    BigReal zn = p[3]->x[localIndex[3]].position.z;
    int ni = (int)floor((zi-pressureProfileMin)/pressureProfileThickness);
    int nj = (int)floor((zj-pressureProfileMin)/pressureProfileThickness);
    int nl = (int)floor((zl-pressureProfileMin)/pressureProfileThickness);
    int nm = (int)floor((zm-pressureProfileMin)/pressureProfileThickness);
    int nn = (int)floor((zn-pressureProfileMin)/pressureProfileThickness);
    pp_clamp(ni, pressureProfileSlabs);
    pp_clamp(nj, pressureProfileSlabs);
    pp_clamp(nl, pressureProfileSlabs);
    pp_clamp(nm, pressureProfileSlabs);
    pp_clamp(nn, pressureProfileSlabs);
    int pi = p[0]->x[localIndex[0]].partition;
    int pj = p[0]->x[localIndex[0]+1].partition;
    int pl = p[1]->x[localIndex[1]].partition;
    int pm = p[2]->x[localIndex[2]].partition;
    int pn = p[3]->x[localIndex[3]].partition;
    int pt = pressureProfileAtomTypes;
    pp_reduction(pressureProfileSlabs, nj, ni,
        pj, pi, pt, fj.x * dr.x, fj.y * dr.y, fj.z * dr.z,
        pressureProfileData);
    pp_reduction(pressureProfileSlabs, ni, nl,
        pi, pl, pt, fi_l.x * u1.x, fi_l.y * u1.y, fi_l.z * u1.z,
        pressureProfileData);
    pp_reduction(pressureProfileSlabs, nm, nn,
        pm, pn, pt, fm.x * u2.x, fm.y * u2.y, fm.z * u2.z,
        pressureProfileData);
  }
#endif
}


void AnisoElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_ANISO_ENERGY) += data[anisoEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

