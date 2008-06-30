/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeNonbondedSelf.h"
#include "ReductionMgr.h"
#include "Patch.h"
#include "LdbCoordinator.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

ComputeNonbondedSelf::ComputeNonbondedSelf(ComputeID c, PatchID pid,
		ComputeNonbondedWorkArrays* _workArrays,
		int minPartition, int maxPartition, int numPartitions)
  : ComputePatch(c,pid), workArrays(_workArrays),
    minPart(minPartition), maxPart(maxPartition), numParts(numPartitions)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  if (pressureProfileOn) {
    pressureProfileReduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_PPROF_NONBONDED);
    int n = pressureProfileAtomTypes;
    pressureProfileData = new BigReal[3*n*n*pressureProfileSlabs];
  } else {
    pressureProfileReduction = NULL;
    pressureProfileData = NULL;
  }
  pairlistsValid = 0;
  pairlistTolerance = 0.;
}

void ComputeNonbondedSelf::initialize() {
  ComputePatch::initialize();
  avgPositionBox = patch->registerAvgPositionPickup(cid);
}

ComputeNonbondedSelf::~ComputeNonbondedSelf()
{
  delete reduction;
  delete pressureProfileReduction;
  delete [] pressureProfileData;
  if (avgPositionBox != NULL) {
    patch->unregisterAvgPositionPickup(cid,&avgPositionBox);
  }
}

#ifdef MEM_OPT_VERSION
void ComputeNonbondedSelf::doForce(CompAtom* p, CompAtomExt* pExt, Results* r)
#else
void ComputeNonbondedSelf::doForce(CompAtom* p, Results* r)
#endif
{
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called
  LdbCoordinator::Object()->startWork(cid,0); // Timestep not used

#ifdef TRACE_COMPUTE_OBJECTS
  double traceObjStartTime = CmiWallTimer();
#endif

  DebugM(2,"doForce() called.\n");
  DebugM(1,numAtoms << " patch 1 atoms\n");
  DebugM(3, "NUMATOMSxNUMATOMS = " << numAtoms*numAtoms << "\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;
  if (pressureProfileOn) {
    int n = pressureProfileAtomTypes;
    memset(pressureProfileData, 0, 3*n*n*pressureProfileSlabs*sizeof(BigReal));
    // adjust lattice dimensions to allow constant pressure
    const Lattice &lattice = patch->lattice;
    pressureProfileThickness = lattice.c().z / pressureProfileSlabs;
    pressureProfileMin = lattice.origin().z - 0.5*lattice.c().z;
  }
  if ( patch->flags.doNonbonded )
  {
    plint maxa = (plint)(-1);
    if ( numAtoms > maxa ) {
      char estr[1024];
      sprintf(estr,"patch has %d atoms, maximum allowed is %d",numAtoms,maxa);
      NAMD_die(estr); 
    }

    int doEnergy = patch->flags.doEnergy;
    nonbonded params;
    params.offset = 0.;
    params.p[0] = p;
    params.p[1] = p;
#ifdef MEM_OPT_VERSION
    params.pExt[0] = pExt;
    params.pExt[1] = pExt;
#endif
    params.ff[0] = r->f[Results::nbond];
    params.ff[1] = r->f[Results::nbond];
    params.numAtoms[0] = numAtoms;
    params.numAtoms[1] = numAtoms;

    // DMK - Atom Separation (water vs. non-water)
    #if NAMD_SeparateWaters != 0
      params.numWaterAtoms[0] = numWaterAtoms;
      params.numWaterAtoms[1] = numWaterAtoms;
    #endif

    params.reduction = reductionData;
    params.pressureProfileReduction = pressureProfileData;

    params.minPart = minPart;
    params.maxPart = maxPart;
    params.numParts = numParts;

    params.workArrays = workArrays;

    params.pairlists = &pairlists;
    params.savePairlists = 0;
    params.usePairlists = 0;
    if ( patch->flags.savePairlists ) {
      params.savePairlists = 1;
      params.usePairlists = 1;
    } else if ( patch->flags.usePairlists ) {
      if ( ! pairlistsValid ||
           ( 2. * patch->flags.maxAtomMovement > pairlistTolerance ) ) {
        reductionData[pairlistWarningIndex] += 1;
      } else { 
        params.usePairlists = 1;
      }
    }
    if ( ! params.usePairlists ) {
      pairlistsValid = 0;
    }
    params.plcutoff = cutoff;
    params.groupplcutoff = cutoff + 2. * patch->flags.maxGroupRadius;
    if ( params.savePairlists ) {
      pairlistsValid = 1;
      pairlistTolerance = 2. * patch->flags.pairlistTolerance;
      params.plcutoff += pairlistTolerance;
      params.groupplcutoff += pairlistTolerance;
    }

    if ( patch->flags.doFullElectrostatics )
    {
      params.fullf[0] = r->f[Results::slow];
      params.fullf[1] = r->f[Results::slow];
      if ( patch->flags.doMolly ) {
        if ( doEnergy ) calcSelfEnergy(&params);
	else calcSelf(&params);
        CompAtom *p_avg = avgPositionBox->open();
        params.p[0] = p_avg;
        params.p[1] = p_avg;
        if ( doEnergy ) calcSlowSelfEnergy(&params);
	else calcSlowSelf(&params);
        avgPositionBox->close(&p_avg);
      } else if ( patch->flags.maxForceMerged == Results::slow ) {
        if ( doEnergy ) calcMergeSelfEnergy(&params);
	else calcMergeSelf(&params);
      } else {
        if ( doEnergy ) calcFullSelfEnergy(&params);
	else calcFullSelf(&params);
      }
    }
    else
      if ( doEnergy ) calcSelfEnergy(&params);
      else calcSelf(&params);
  }

  submitReductionData(reductionData,reduction);
  if (pressureProfileOn)
    submitPressureProfileData(pressureProfileData, pressureProfileReduction);

#ifdef TRACE_COMPUTE_OBJECTS
    traceUserBracketEvent(TRACE_COMPOBJ_IDOFFSET+cid, traceObjStartTime, CmiWallTimer());
#endif

  // Inform load balancer
  LdbCoordinator::Object()->endWork(cid,0); // Timestep not used

  reduction->submit();
  if (pressureProfileOn)
    pressureProfileReduction->submit();
}

