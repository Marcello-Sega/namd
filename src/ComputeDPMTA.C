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
#include "ComputeDPMTA.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

ComputeDPMTA::ComputeDPMTA(ComputeID c) : ComputeHomePatches(c)
{
  ;
}

ComputeDPMTA::~ComputeDPMTA()
{
  ;
}


void ComputeDPMTA::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).x = (*ap).positionBox->open();
    (*ap).a = (*ap).atomBox->open();



    (*ap).positionBox->close(&(*ap).x);
    (*ap).atomBox->close(&(*ap).a);
  } 

  ++fake_seq;

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    (*ap).f = (*ap).forceBox->open();



    (*ap).forceBox->close(&(*ap).f);
  }
}

