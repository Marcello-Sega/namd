/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeNonbondedPair.h"
#include "ReductionMgr.h"
#include "Patch.h"
#include "LdbCoordinator.h"
#include "PatchMap.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

ComputeNonbondedPair::ComputeNonbondedPair(ComputeID c, PatchID pid[], int trans[])
  : ComputePatchPair(c,pid,trans)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

void ComputeNonbondedPair::initialize() {
  ComputePatchPair::initialize();
  for (int i=0; i<2; i++) {
    avgPositionBox[i] = patch[i]->registerAvgPositionPickup(cid);
  }
}

ComputeNonbondedPair::~ComputeNonbondedPair()
{
  delete reduction;
  for (int i=0; i<2; i++) {
    if (avgPositionBox[i] != NULL) {
      patch[i]->unregisterAvgPositionPickup(cid,&avgPositionBox[i]);
    }
  }
}

int ComputeNonbondedPair::noWork() {

  // return 0;  // for testing
  if ( numAtoms[0] && numAtoms[1] && patch[0]->flags.doNonbonded )
  {
    return 0;  // work to do, enqueue as usual
  } else
  {
    // Inform load balancer
    LdbCoordinator::Object()->startWork(cid,0); // Timestep not used
    // fake out patches and reduction system

    BigReal reductionData[reductionDataSize];
    int i;
    for ( i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

    CompAtom* p[2];
    CompAtom* p_avg[2];
    Results* r[2];

    // Open up positionBox, forceBox, and atomBox
    for (i=0; i<2; i++) {
      p[i] = positionBox[i]->open();
      r[i] = forceBox[i]->open();
      if ( patch[0]->flags.doMolly ) p_avg[i] = avgPositionBox[i]->open();
    }

    // Close up boxes
    for (i=0; i<2; i++) {
      positionBox[i]->close(&p[i]);
      forceBox[i]->close(&r[i]);
      if ( patch[0]->flags.doMolly ) avgPositionBox[i]->close(&p_avg[i]);
    }

    submitReductionData(reductionData,reduction);
    reduction->submit();

    // Inform load balancer
    LdbCoordinator::Object()->endWork(cid,0); // Timestep not used

    return 1;  // no work to do, do not enqueue
  }
}


void ComputeNonbondedPair::doForce(CompAtom* p[2],
                               Results* r[2])
{
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called
  LdbCoordinator::Object()->startWork(cid,0); // Timestep not used

  DebugM(2,"doForce() called.\n");
  DebugM(2, numAtoms[0] << " patch #1 atoms and " <<
	numAtoms[1] << " patch #2 atoms\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;

  if ( numAtoms[0] && numAtoms[1] )
  {
    nonbonded params;
    params.reduction = reductionData;

    // swap to place more atoms in inner loop (second patch)
    if ( numAtoms[0] > numAtoms[1] )
    {
      params.p[0] = p[1];
      params.p[1] = p[0];
      params.ff[0] = r[1]->f[Results::nbond];
      params.ff[1] = r[0]->f[Results::nbond];
      params.numAtoms[0] = numAtoms[1];
      params.numAtoms[1] = numAtoms[0];
      DebugM(3, "NUMATOMSxNUMATOMS = " << numAtoms[0]*numAtoms[1] << "\n" );
      if ( patch[0]->flags.doFullElectrostatics )
      {
	params.fullf[0] = r[1]->f[Results::slow];
	params.fullf[1] = r[0]->f[Results::slow];
	if ( patch[0]->flags.doMolly ) {
          calcPair(&params);
	  CompAtom *p_avg[2];
	  p_avg[0] = avgPositionBox[0]->open();
	  p_avg[1] = avgPositionBox[1]->open();
	  params.p[0] = p_avg[1];
	  params.p[1] = p_avg[0];
	  calcSlowPair(&params);
	  avgPositionBox[0]->close(&p_avg[0]);
	  avgPositionBox[1]->close(&p_avg[1]);
	} else {
	  calcFullPair(&params);
	}
      }
      else
        calcPair(&params);
    }
    else
    {
      params.p[0] = p[0];
      params.p[1] = p[1];
      params.numAtoms[0] = numAtoms[0];
      params.numAtoms[1] = numAtoms[1];
      params.ff[0] = r[0]->f[Results::nbond];
      params.ff[1] = r[1]->f[Results::nbond];
      if ( patch[0]->flags.doFullElectrostatics )
      {
        params.fullf[0] = r[0]->f[Results::slow];
        params.fullf[1] = r[1]->f[Results::slow];
	if ( patch[0]->flags.doMolly ) {
          calcPair(&params);
	  CompAtom *p_avg[2];
	  p_avg[0] = avgPositionBox[0]->open();
	  p_avg[1] = avgPositionBox[1]->open();
	  params.p[0] = p_avg[0];
	  params.p[1] = p_avg[1];
	  calcSlowPair(&params);
	  avgPositionBox[0]->close(&p_avg[0]);
	  avgPositionBox[1]->close(&p_avg[1]);
	} else {
	  calcFullPair(&params);
	}
      }
      else
        calcPair(&params);
    }
  }

  submitReductionData(reductionData,reduction);
  reduction->submit();

  // Inform load balancer
  LdbCoordinator::Object()->endWork(cid,0); // Timestep not used
}

