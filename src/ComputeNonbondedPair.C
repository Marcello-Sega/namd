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

ComputeNonbondedPair::ComputeNonbondedPair(ComputeID c, PatchID pid[], int trans[],
		ComputeNonbondedWorkArrays* _workArrays,
		int minPartition, int maxPartition, int numPartitions)
  : ComputePatchPair(c,pid,trans), workArrays(_workArrays),
    minPart(minPartition), maxPart(maxPartition), numParts(numPartitions)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  if (pressureProfileOn) {
    int n = pressureProfileAtomTypes; 
    pressureProfileData = new BigReal[3*n*n*pressureProfileSlabs];
    pressureProfileReduction = ReductionMgr::Object()->willSubmit(
	REDUCTIONS_PPROF_NONBONDED, 3*pressureProfileSlabs*((n*(n+1))/2));
  } else {
    pressureProfileReduction = NULL;
    pressureProfileData = NULL;
  }
  pairlistsValid = 0;
  pairlistTolerance = 0.;
}

void ComputeNonbondedPair::initialize() {
  ComputePatchPair::initialize();
  for (int i=0; i<2; i++) {
    avgPositionBox[i] = patch[i]->registerAvgPositionPickup(cid,trans[i]);
  }
#ifdef NAMD_CUDA
  register_cuda_compute_pair(cid, patchID, trans);
#endif
}

ComputeNonbondedPair::~ComputeNonbondedPair()
{
  delete reduction;
  delete pressureProfileReduction;
  delete [] pressureProfileData;
  for (int i=0; i<2; i++) {
    if (avgPositionBox[i] != NULL) {
      patch[i]->unregisterAvgPositionPickup(cid,&avgPositionBox[i]);
    }
  }
}

int ComputeNonbondedPair::noWork() {

  // return 0;  // for testing
  if ( numAtoms[0] && numAtoms[1] && patch[0]->flags.doNonbonded
#ifdef NAMD_CUDA
	&& patch[0]->flags.doEnergy
#endif
 )
  {
    return 0;  // work to do, enqueue as usual
  } else
  {
    // Inform load balancer
#ifndef NAMD_CUDA
    LdbCoordinator::Object()->startWork(cid,0); // Timestep not used
#endif
    // fake out patches and reduction system

    BigReal reductionData[reductionDataSize];
    int i;
    for ( i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;
    if (pressureProfileOn) {
      int n = pressureProfileAtomTypes;
      memset(pressureProfileData, 0, 3*n*n*pressureProfileSlabs*sizeof(BigReal));
    }
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
    if (pressureProfileOn)
      submitPressureProfileData(pressureProfileData, pressureProfileReduction);

    // Inform load balancer
#ifndef NAMD_CUDA
    LdbCoordinator::Object()->endWork(cid,0); // Timestep not used
#endif

    reduction->submit();
    if (pressureProfileOn) 
      pressureProfileReduction->submit();

    return 1;  // no work to do, do not enqueue
  }
}

#ifdef MEM_OPT_VERSION
void ComputeNonbondedPair::doForce(CompAtom* p[2], CompAtomExt* pExt[2], Results* r[2])
#else
void ComputeNonbondedPair::doForce(CompAtom* p[2], Results* r[2])
#endif
{
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called
#ifndef NAMD_CUDA
  LdbCoordinator::Object()->startWork(cid,0); // Timestep not used
#endif

#ifdef TRACE_COMPUTE_OBJECTS
    double traceObjStartTime = CmiWallTimer();
#endif

  DebugM(2,"doForce() called.\n");
  DebugM(2, numAtoms[0] << " patch #1 atoms and " <<
	numAtoms[1] << " patch #2 atoms\n");

  BigReal reductionData[reductionDataSize];
  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;
  if (pressureProfileOn) {
    int n = pressureProfileAtomTypes; 
    memset(pressureProfileData, 0, 3*n*n*pressureProfileSlabs*sizeof(BigReal));
    // adjust lattice dimensions to allow constant pressure
    const Lattice &lattice = patch[0]->lattice;
    pressureProfileThickness = lattice.c().z / pressureProfileSlabs;
    pressureProfileMin = lattice.origin().z - 0.5*lattice.c().z;
  }
  if ( numAtoms[0] && numAtoms[1] )
  {
    nonbonded params;
    params.reduction = reductionData;
    params.pressureProfileReduction = pressureProfileData;

    params.minPart = minPart;
    params.maxPart = maxPart;
    params.numParts = numParts;

    params.workArrays = workArrays;

    params.pairlists = &pairlists;
    params.savePairlists = 0;
    params.usePairlists = 0;
    if ( patch[0]->flags.savePairlists ) {
      params.savePairlists = 1;
      params.usePairlists = 1;
    } else if ( patch[0]->flags.usePairlists && patch[1]->flags.usePairlists ) {
      if ( ! pairlistsValid ||
           ( patch[0]->flags.maxAtomMovement +
             patch[1]->flags.maxAtomMovement > pairlistTolerance ) ) {
        reductionData[pairlistWarningIndex] += 1;
      } else {
        params.usePairlists = 1;
      }
    }
    if ( ! params.usePairlists ) {
      pairlistsValid = 0;
    }
    params.plcutoff = cutoff;
    params.groupplcutoff = cutoff +
	patch[0]->flags.maxGroupRadius + patch[1]->flags.maxGroupRadius;
    if ( params.savePairlists ) {
      pairlistsValid = 1;
      pairlistTolerance = patch[0]->flags.pairlistTolerance +
                          patch[1]->flags.pairlistTolerance;
      params.plcutoff += pairlistTolerance;
      params.groupplcutoff += pairlistTolerance;
    }

    // swap to place more atoms in inner loop (second patch)
    int a = 0;  int b = 1;
    if ( numAtoms[0] > numAtoms[1] ) { a = 1; b = 0; }

    const Lattice &lattice = patch[0]->lattice;
    params.offset = lattice.offset(trans[a]) - lattice.offset(trans[b]);

    // Atom Sorting : If we are sorting the atoms along the line connecting
    //   the patch centers, then calculate a normalized vector pointing from
    //   patch a to patch b (i.e. outer loop patch to inner loop patch).
    #if NAMD_ComputeNonbonded_SortAtoms != 0

      // Center of patch a (outer-loop; i-loop) and patch b (inner-loop; j/k-loop)
      PatchMap* patchMap = PatchMap::Object();
      ScaledPosition p_a_center, p_b_center;
      p_a_center.x = patchMap->min_a(patchID[a]);
      p_a_center.y = patchMap->min_b(patchID[a]);
      p_a_center.z = patchMap->min_c(patchID[a]);
      p_b_center.x = patchMap->min_a(patchID[b]);
      p_b_center.y = patchMap->min_b(patchID[b]);
      p_b_center.z = patchMap->min_c(patchID[b]);
      p_a_center.x += patchMap->max_a(patchID[a]);
      p_a_center.y += patchMap->max_b(patchID[a]);
      p_a_center.z += patchMap->max_c(patchID[a]);
      p_b_center.x += patchMap->max_a(patchID[b]);
      p_b_center.y += patchMap->max_b(patchID[b]);
      p_b_center.z += patchMap->max_c(patchID[b]);
      p_a_center *= (BigReal)0.5;
      p_b_center *= (BigReal)0.5;
      p_a_center = lattice.unscale(p_a_center);
      p_b_center = lattice.unscale(p_b_center);

      // Adjust patch a's center by the offset
      p_a_center.x += params.offset.x;
      p_a_center.y += params.offset.y;
      p_a_center.z += params.offset.z;

      // Calculate and fill in the projected line vector
      params.projLineVec = p_b_center - p_a_center;
      params.projLineVec /= params.projLineVec.length(); // Normalize the vector

    #endif

    int doEnergy = patch[0]->flags.doEnergy;
      params.p[0] = p[a];
      params.p[1] = p[b];
    #ifdef MEM_OPT_VERSION      
      params.pExt[0] = pExt[a]; 
      params.pExt[1] = pExt[b];
    #endif
      params.ff[0] = r[a]->f[Results::nbond];
      params.ff[1] = r[b]->f[Results::nbond];
#ifdef NAMD_CUDA
      params.ff[0] = r[a]->f[Results::slow];
      params.ff[1] = r[b]->f[Results::slow];
#endif
      params.numAtoms[0] = numAtoms[a];
      params.numAtoms[1] = numAtoms[b];

      // DMK - Atom Separation (water vs. non-water)
      #if NAMD_SeparateWaters != 0
        params.numWaterAtoms[0] = numWaterAtoms[a];
        params.numWaterAtoms[1] = numWaterAtoms[b];
      #endif

      if ( patch[0]->flags.doFullElectrostatics )
      {
	params.fullf[0] = r[a]->f[Results::slow];
	params.fullf[1] = r[b]->f[Results::slow];
	if ( patch[0]->flags.doMolly ) {
          if ( doEnergy ) calcPairEnergy(&params);
	  else calcPair(&params);
	  CompAtom *p_avg[2];
	  p_avg[0] = avgPositionBox[0]->open();
	  p_avg[1] = avgPositionBox[1]->open();
	  params.p[0] = p_avg[a];
	  params.p[1] = p_avg[b];
	  if ( doEnergy ) calcSlowPairEnergy(&params);
	  else calcSlowPair(&params);
	  avgPositionBox[0]->close(&p_avg[0]);
	  avgPositionBox[1]->close(&p_avg[1]);
        } else if ( patch[0]->flags.maxForceMerged == Results::slow ) {
          if ( doEnergy ) calcMergePairEnergy(&params);
	  else calcMergePair(&params);
	} else {
	  if ( doEnergy ) calcFullPairEnergy(&params);
	  else calcFullPair(&params);
	}
      }
      else
        if ( doEnergy ) calcPairEnergy(&params);
        else calcPair(&params);
  }

  submitReductionData(reductionData,reduction);
  if (pressureProfileOn)
    submitPressureProfileData(pressureProfileData, pressureProfileReduction);

#ifdef TRACE_COMPUTE_OBJECTS
    traceUserBracketEvent(TRACE_COMPOBJ_IDOFFSET+cid, traceObjStartTime, CmiWallTimer());
#endif

  // Inform load balancer
#ifndef NAMD_CUDA
  LdbCoordinator::Object()->endWork(cid,0); // Timestep not used
#endif

  reduction->submit();
  if (pressureProfileOn)
    pressureProfileReduction->submit();
}

