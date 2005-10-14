/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Patch.h"
#include "PatchMap.h"
#include "Compute.h"

#include "AtomMap.h"
#include "ComputeMap.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "ResizeArrayPrimIter.h"

#include "Sync.h"

typedef ResizeArrayPrimIter<ComputeID> ComputeIDListIter;

//#define  DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

Patch::Patch(PatchID pd) :
   lattice(flags.lattice),
   patchID(pd), numAtoms(0), numFixedAtoms(0),
   positionPtr(0), avgPositionPtr(0),
   positionBox(this,&Patch::positionBoxClosed),
   avgPositionBox(this,&Patch::avgPositionBoxClosed),
   forceBox(this,&Patch::forceBoxClosed),
   boxesOpen(0), _hasNewAtoms(0)
{
  lattice = Node::Object()->simParameters->lattice;
}

PositionBox<Patch>* Patch::registerPositionPickup(ComputeID cid, int trans)
{
   //DebugM(4, "registerPositionPickupa("<<patchID<<") from " << cid << "\n");
   if (positionComputeList.add(cid) < 0)
   {
     DebugM(7, "registerPositionPickup() failed for cid " << cid << endl);
     return NULL;
   }
   return positionBox.checkOut(trans);
}

void Patch::unregisterPositionPickup(ComputeID cid, PositionBox<Patch> **const box)
{
   DebugM(4, "UnregisterPositionPickup from " << cid << "\n");
   positionComputeList.del(cid);
   positionBox.checkIn(*box);
   *box = 0;
}

PositionBox<Patch>* Patch::registerAvgPositionPickup(ComputeID cid, int trans)
{
   //DebugM(4, "registerAvgPositionPickup("<<patchID<<") from " << cid << "\n");
   return avgPositionBox.checkOut(trans);
}

void Patch::unregisterAvgPositionPickup(ComputeID cid, PositionBox<Patch> **const box)
{
   DebugM(4, "UnregisterAvgPositionPickup from " << cid << "\n");
   avgPositionBox.checkIn(*box);
   *box = 0;
}

Box<Patch,Results>* Patch::registerForceDeposit(ComputeID cid)
{
   if (forceComputeList.add(cid) < 0)
   {
     DebugM(7, "registerForceDeposit() failed for cid " << cid << endl);
     DebugM(7, "  size of forceCompueList " << forceComputeList.size() << endl);
     return NULL;
   }
   return forceBox.checkOut();
}

void Patch::unregisterForceDeposit(ComputeID cid, Box<Patch,Results> **const box)
{
   DebugM(4, "unregisterForceDeposit() computeID("<<cid<<")"<<endl);
   forceComputeList.del(cid);
   forceBox.checkIn(*box);
   *box = 0;
}

void Patch::positionBoxClosed(void)
{
   positionPtr = 0;
   this->boxClosed(0);
}

void Patch::forceBoxClosed(void)
{
   DebugM(4, "patchID("<<patchID<<") forceBoxClosed! call\n");
   for (int j = 0; j < Results::maxNumForces; ++j )
   {
     results.f[j] = 0;
   }
   this->boxClosed(1);
}

void Patch::avgPositionBoxClosed(void)
{
   avgPositionPtr = 0;
   this->boxClosed(3);
}

// void Patch::boxClosed(int box) is virtual

void Patch::positionsReady(int doneMigration)
{
   DebugM(4,"Patch::positionsReady() - patchID(" << patchID <<")"<<endl );
   ComputeMap *computeMap = ComputeMap::Object();

   if ( doneMigration ) AtomMap::Object()->registerIDs(patchID,p.begin(),p.end());

   boxesOpen = 2;
   if ( flags.doMolly ) boxesOpen++;
   _hasNewAtoms = (doneMigration != 0);

   // Give all position pickup boxes access to positions
   positionPtr = p.begin();
   positionBox.open(positionPtr,numAtoms,&lattice);
   if ( flags.doMolly ) {
     avgPositionPtr = p_avg.begin();
     avgPositionBox.open(avgPositionPtr,numAtoms,&lattice);
   }

   // Give all force deposit boxes access to forces
   Force *forcePtr;
   for ( int j = 0; j < Results::maxNumForces; ++j )
   {
      f[j].resize(numAtoms);
      forcePtr = f[j].begin();
      for(register int i=0; i<numAtoms; i++) forcePtr[i] = 0.;
      results.f[j] = forcePtr;
   }
   forceBox.open(&results);

   // Iterate over compute objects that need to be informed we are ready
   ComputeIDListIter cid(positionComputeList);
   // gzheng
   if (useSync) {
     if (Sync::Object()->holdComputes(patchID, cid, doneMigration))
       return;
   }

   int compute_count = 0;
   int seq = flags.sequence;
   for(cid = cid.begin(); cid != cid.end(); cid++)
   {
         compute_count++;
	 computeMap->compute(*cid)->patchReady(patchID,doneMigration,seq);
   }
   if (compute_count == 0 && PatchMap::Object()->node(patchID) != CkMyPe()) {
       iout << iINFO << "PATCH_COUNT: Patch " << patchID 
	    << " on PE " << CkMyPe() <<" home patch " 
	    << PatchMap::Object()->node(patchID)
	    << " does not have any computes\n" 
	    << endi;
   }
}


