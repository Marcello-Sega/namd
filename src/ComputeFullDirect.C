/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#include "Namd.h"
#include "Node.h"
#include "PatchMap.h"
#include "AtomMap.h"
#include "ComputeFullDirect.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

// Only works on one processor, without periodic boundary conditions.  -JCP

ComputeFullDirect::ComputeFullDirect(ComputeID c) : ComputeHomePatches(c)
{
  reduction->Register(REDUCTION_ELECT_ENERGY);
}

ComputeFullDirect::~ComputeFullDirect()
{
  reduction->unRegister(REDUCTION_ELECT_ENERGY);
}


void ComputeFullDirect::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);
  register int i,j;

  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics )
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      Position *x = (*ap).positionBox->open();
      AtomProperties *a = (*ap).atomBox->open();
      Force *f = (*ap).forceBox->open();
      reduction->submit(fake_seq, REDUCTION_ELECT_ENERGY, 0.);
      ++fake_seq;
      (*ap).positionBox->close(&x);
      (*ap).atomBox->close(&a);
      (*ap).forceBox->close(&f);
    }
    return;
  }

  // allocate storage
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  localPositions = new Position[numLocalAtoms];	// freed at end of this method
  localCharges = new BigReal[numLocalAtoms];	// freed at end of this method
  localForces = new Force[numLocalAtoms];	// freed at end of this method

  // get positions and charges
  j = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Position *x = (*ap).positionBox->open();
    AtomProperties *a = (*ap).atomBox->open();
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i, ++j)
    {
      localPositions[j] = x[i];
      localCharges[j] = a[i].charge;
    }

    (*ap).positionBox->close(&x);
    (*ap).atomBox->close(&a);
  } 

  // perform local calculations
  BigReal electEnergy = 0;
  const BigReal coloumb = COLOUMB * ComputeNonbondedUtil::dielectric_1;
  for(i=0; i<numLocalAtoms; ++i)
  {
    register BigReal kq_i = coloumb * localCharges[i];
    register BigReal p_i_x = localPositions[i].x;
    register BigReal p_i_y = localPositions[i].y;
    register BigReal p_i_z = localPositions[i].z;
    register BigReal f_i_x = 0.;
    register BigReal f_i_y = 0.;
    register BigReal f_i_z = 0.;
    for(j=i+1; j<numLocalAtoms; ++j)
    {
      register Position *p_j = localPositions + j;
      register BigReal p_ij_x = p_i_x - p_j->x;
      register BigReal p_ij_y = p_i_y - p_j->y;
      register BigReal p_ij_z = p_i_z - p_j->z;

      register BigReal r_1;
      r_1 = 1./sqrt(p_ij_x * p_ij_x + p_ij_y * p_ij_y + p_ij_z * p_ij_z);
      register BigReal f = localCharges[j] * kq_i * r_1;
      electEnergy += f;
      f *= r_1*r_1;
      p_ij_x *= f;
      p_ij_y *= f;
      p_ij_z *= f;
      f_i_x += p_ij_x;
      f_i_y += p_ij_y;
      f_i_z += p_ij_z;
      localForces[j] -= Vector(p_ij_x,p_ij_y,p_ij_z);
    }
    localForces[i] += Vector(f_i_x,f_i_y,f_i_z);
  }

  // send out reductions
  DebugM(4,"Full-electrostatics energy: " << electEnergy << "\n");
  reduction->submit(fake_seq, REDUCTION_ELECT_ENERGY, electEnergy);
  ++fake_seq;

  // add in forces
  j = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Force *f = (*ap).forceBox->open();
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i, ++j)
    {
      f[i] += localForces[j];
    }

    (*ap).forceBox->close(&f);
  }

  // free storage
  delete [] localPositions;	// allocated at beginning of this method
  delete [] localCharges;	// allocated at beginning of this method
  delete [] localForces;	// allocated at beginning of this method
}

