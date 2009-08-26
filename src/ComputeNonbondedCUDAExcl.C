/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "ComputeNonbondedCUDAExcl.h"
#include "Molecule.h"
#include "Parameters.h"
#include "LJTable.h"
#include "Node.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "PressureProfile.h"
#include "Debug.h"


// static initialization
int ExclElem::pressureProfileSlabs = 0;
int ExclElem::pressureProfileAtomTypes = 1;
BigReal ExclElem::pressureProfileThickness = 0;
BigReal ExclElem::pressureProfileMin = 0;

void ExclElem::getMoleculePointers
    (Molecule* mol, int* count, int32*** byatom, Exclusion** structarray)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Should not be called in ExclElem::getMoleculePointers in memory optimized version!");
#else
  *count = mol->numExclusions;
  *byatom = mol->exclusionsByAtom;
  *structarray = mol->exclusions;
#endif
}

void ExclElem::getParameterPointers(Parameters *p, const int **v) {
  *v = 0;
}

void ExclElem::computeForce(BigReal *reduction, 
                            BigReal *pressureProfileData)
{
    const CompAtom &p_i = p[0]->x[localIndex[0]];
    const CompAtom &p_j = p[1]->x[localIndex[1]];

    // compute vectors between atoms and their distances
    const Lattice & lattice = p[0]->p->lattice;
    const Vector r12 = lattice.delta(p_i.position, p_j.position);
    BigReal r2 = r12.length2();

    if ( r2 > cutoff2 ) return;

    r2 += r2_delta;

    union { double f; int64 i; } r2i;
    r2i.f = r2;
    const int r2_delta_expc = 64 * (r2_delta_exp - 1023);
    int table_i = (r2i.i >> (32+14)) + r2_delta_expc;  // table_i >= 0

    BigReal diffa = r2 - r2_table[table_i];
    // XXX for short table_four is set to table_short, check for MTS
    const BigReal* const table_four_i = table_noshort + 16*table_i;

    BigReal fast_a, fast_b, fast_c, fast_d;

  if ( modified ) {

    const LJTable::TableEntry * lj_pars =
            ljTable->table_row(p_i.vdwType) + 2 * p_j.vdwType;

    // modified - normal = correction
    const BigReal A = scaling * ( (lj_pars+1)->A - lj_pars->A );
    const BigReal B = scaling * ( (lj_pars+1)->B - lj_pars->B );
    const BigReal kqq = (scale14 - 1.0) *
            COLOUMB * p_i.charge * p_j.charge * scaling * dielectric_1;

    // iout << " A " << A << " B " << B << " kqq " << kqq << "\n" << endi;

    BigReal vdw_d = A * table_four_i[0] - B * table_four_i[2];
    BigReal vdw_c = A * table_four_i[1] - B * table_four_i[3];
    BigReal vdw_b = A * table_four_i[4] - B * table_four_i[6];
    BigReal vdw_a = A * table_four_i[5] - B * table_four_i[7];

    fast_d = vdw_d + kqq * table_four_i[12];
    fast_c = vdw_c + kqq * table_four_i[13];
    fast_b = vdw_b + kqq * table_four_i[14];
    fast_a = vdw_a + kqq * table_four_i[15];  // not used!

  } else {  // full exclusion

    const BigReal kqq = 
            COLOUMB * p_i.charge * p_j.charge * scaling * dielectric_1;

    fast_d = kqq * ( table_four_i[8]  - table_four_i[12] );
    fast_c = kqq * ( table_four_i[9]  - table_four_i[13] );
    fast_b = kqq * ( table_four_i[10] - table_four_i[14] );
    fast_a = kqq * ( table_four_i[11] - table_four_i[15] );  // not used!

  }

    register BigReal fast_dir =
                  (diffa * fast_d + fast_c) * diffa + fast_b;

  const Force f12 = fast_dir * r12;

  //  Now add the forces to each force vector
  p[0]->r->f[Results::nbond][localIndex[0]] += f12;
  p[1]->r->f[Results::nbond][localIndex[1]] -= f12;

  // reduction[nonbondedEnergyIndex] += energy;
  reduction[virialIndex_XX] += f12.x * r12.x;
  // reduction[virialIndex_XY] += f12.x * r12.y;
  // reduction[virialIndex_XZ] += f12.x * r12.z;
  // reduction[virialIndex_YX] += f12.y * r12.x;
  reduction[virialIndex_YY] += f12.y * r12.y;
  // reduction[virialIndex_YZ] += f12.y * r12.z;
  // reduction[virialIndex_ZX] += f12.z * r12.x;
  // reduction[virialIndex_ZY] += f12.z * r12.y;
  reduction[virialIndex_ZZ] += f12.z * r12.z;

}

void ExclElem::submitReductionData(BigReal *data, SubmitReduction *reduction)
{
  // reduction->item(REDUCTION_BOND_ENERGY) += data[bondEnergyIndex];
  ADD_TENSOR(reduction,REDUCTION_VIRIAL_NBOND,data,virialIndex);
}

