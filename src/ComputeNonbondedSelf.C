/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedSelf.C
 *
 ***************************************************************************/

#include "ComputeNonbondedSelf.h"
#include "ReductionMgr.h"
#include "Patch.h"
#include "LdbCoordinator.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

ComputeNonbondedSelf::ComputeNonbondedSelf(ComputeID c, PatchID pid,
		int minPartition, int maxPartition, int numPartitions)
  : ComputePatch(c,pid),
    minPart(minPartition), maxPart(maxPartition), numParts(numPartitions)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

void ComputeNonbondedSelf::initialize() {
  ComputePatch::initialize();
  avgPositionBox = patch->registerAvgPositionPickup(cid);
}

ComputeNonbondedSelf::~ComputeNonbondedSelf()
{
  delete reduction;
  if (avgPositionBox != NULL) {
    patch->unregisterAvgPositionPickup(cid,&avgPositionBox);
  }
}


void ComputeNonbondedSelf::doForce(Position* p,
                               Results* r,
                               AtomProperties* a)
{
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called
  LdbCoordinator::Object()->startWork(cid,0); // Timestep not used

  DebugM(2,"doForce() called.\n");
  DebugM(1,numAtoms << " patch 1 atoms\n");
  DebugM(3, "NUMATOMSxNUMATOMS = " << numAtoms*numAtoms << "\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

  if ( patch->flags.doNonbonded )
  {
    nonbonded params;
    params.p[0] = p;
    params.p[1] = p;
    params.ff[0] = r->f[Results::nbond];
    params.ff[1] = r->f[Results::nbond];
    params.a[0] = a;
    params.a[1] = a;
    params.numAtoms[0] = numAtoms-1;
    params.numAtoms[1] = numAtoms;
    params.reduction = reductionData;

    params.minPart = minPart;
    params.maxPart = maxPart;
    params.numParts = numParts;

    if ( patch->flags.doFullElectrostatics )
    {
      params.fullf[0] = r->f[Results::slow];
      params.fullf[1] = r->f[Results::slow];
      if ( patch->flags.doMolly ) {
        calcSelf(&params);
        Position *p_avg = avgPositionBox->open();
        params.p[0] = p_avg;
        params.p[1] = p_avg;
        calcSlowSelf(&params);
        avgPositionBox->close(&p_avg);
      } else {
        calcFullSelf(&params);
      }
    }
    else
      calcSelf(&params);
  }

  submitReductionData(reductionData,reduction);
  reduction->submit();
  // Inform load balancer
  LdbCoordinator::Object()->endWork(cid,0); // Timestep not used
}

