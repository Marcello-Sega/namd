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
  // calculation and (most) comments from Ed Harder's implementation in CHARMM
#if 0
  fprintf(stderr, "AnisoElem::computeForce() -  localIndex[] = %d %d %d %d\n",
      localIndex[0], localIndex[1], localIndex[2], localIndex[3]);
  fprintf(stderr, "     id = %d %d %d %d\n",
      p[0]->xExt[localIndex[0]].id,
      p[1]->xExt[localIndex[1]].id,
      p[2]->xExt[localIndex[2]].id,
      p[3]->xExt[localIndex[3]].id);
#endif

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

  Vector fi = kiso0 * dr;   // iso spring force
  Vector fj = -kiso0 * dr;

  // par/perp spring forces
  fi.x += kpar0*dpar*(u1.x - (dr.x + u1.x*dpar)*u1_invlen) + kperp0*dperp*u2.x;
  fi.y += kpar0*dpar*(u1.y - (dr.y + u1.y*dpar)*u1_invlen) + kperp0*dperp*u2.y;
  fi.z += kpar0*dpar*(u1.z - (dr.z + u1.z*dpar)*u1_invlen) + kperp0*dperp*u2.z;

  fj.x -= kpar0*dpar*u1.x + kperp0*dperp*u2.x;
  fj.y -= kpar0*dpar*u1.y + kperp0*dperp*u2.y;
  fj.z -= kpar0*dpar*u1.z + kperp0*dperp*u2.z;

  Vector fl, fm, fn;

  fl.x = kpar0*dpar*(dr.x - u1.x*dpar)*u1_invlen;
  fl.y = kpar0*dpar*(dr.y - u1.y*dpar)*u1_invlen;
  fl.z = kpar0*dpar*(dr.z - u1.z*dpar)*u1_invlen;

  fm.x = -kperp0*dperp*(dr.x - u2.x*dperp)*u2_invlen;
  fm.y = -kperp0*dperp*(dr.y - u2.y*dperp)*u2_invlen;
  fm.z = -kperp0*dperp*(dr.z - u2.z*dperp)*u2_invlen;

  fn.x = kperp0*dperp*(dr.x - u2.x*dperp)*u2_invlen;
  fn.y = kperp0*dperp*(dr.y - u2.y*dperp)*u2_invlen;
  fn.z = kperp0*dperp*(dr.z - u2.z*dperp)*u2_invlen;

  // accumulate forces
  p[0]->f[localIndex[0]] += fi;
  p[0]->f[localIndex[0]+1] += fj;
  p[1]->f[localIndex[1]] += fl;
  p[2]->f[localIndex[2]] += fm;
  p[3]->f[localIndex[3]] += fn;

  // update potential
  reduction[anisoEnergyIndex] += eaniso;

#if 0
  // update virial
  // XXX won't the virial get messed up whenever an atom wraps around
  // periodic boundaries?
  Tensor v = outer(fi,ri);
  v += outer(fj,rj);
  v += outer(fl,rl);
  v += outer(fm,rm);
  v += outer(fn,rn);

  reduction[virialIndex_XX] = v.xx;
  reduction[virialIndex_XY] = v.xy;
  reduction[virialIndex_XZ] = v.xz;
  reduction[virialIndex_YX] = v.yx;
  reduction[virialIndex_YY] = v.yy;
  reduction[virialIndex_YZ] = v.yz;
  reduction[virialIndex_ZX] = v.zx;
  reduction[virialIndex_ZY] = v.zy;
  reduction[virialIndex_ZZ] = v.zz;
#endif

  //fprintf(stderr, "     eaniso = %g\n", eaniso);
#endif
}


void AnisoElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  reduction->item(REDUCTION_ANISO_ENERGY) += data[anisoEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NORMAL,data,virialIndex);
}

