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
   avgPositionPtrBegin(0), avgPositionPtrEnd(0),
   positionBox(this,&Patch::positionBoxClosed),
   avgPositionBox(this,&Patch::avgPositionBoxClosed),
   forceBox(this,&Patch::forceBoxClosed),
   boxesOpen(0), _hasNewAtoms(0)

   // DMK - Atom Separation (water vs. non-water)
   #if NAMD_SeparateWaters != 0
     ,numWaterAtoms(-1)
   #endif
{
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
    positionPtrBegin = 0;
    positionPtrEnd = 0;
#endif
  lattice = Node::Object()->simParameters->lattice;
}

Box<Patch,CompAtom>* Patch::registerPositionPickup(ComputeID cid, int trans)
{
   //DebugM(4, "registerPositionPickupa("<<patchID<<") from " << cid << "\n");
   if (positionComputeList.add(cid) < 0)
   {
     DebugM(7, "registerPositionPickup() failed for cid " << cid << std::endl);
     return NULL;
   }
   return positionBox.checkOut();
}

void Patch::unregisterPositionPickup(ComputeID cid, Box<Patch,CompAtom> **const box)
{
   DebugM(4, "UnregisterPositionPickup from " << cid << "\n");
   positionComputeList.del(cid);
   positionBox.checkIn(*box);
   *box = 0;
}

Box<Patch,CompAtom>* Patch::registerAvgPositionPickup(ComputeID cid, int trans)
{
   //DebugM(4, "registerAvgPositionPickup("<<patchID<<") from " << cid << "\n");
   return avgPositionBox.checkOut();
}

void Patch::unregisterAvgPositionPickup(ComputeID cid, Box<Patch,CompAtom> **const box)
{
   DebugM(4, "UnregisterAvgPositionPickup from " << cid << "\n");
   avgPositionBox.checkIn(*box);
   *box = 0;
}

Box<Patch,Results>* Patch::registerForceDeposit(ComputeID cid)
{
   if (forceComputeList.add(cid) < 0)
   {
     DebugM(7, "registerForceDeposit() failed for cid " << cid << std::endl);
     DebugM(7, "  size of forceCompueList " << forceComputeList.size() << std::endl);
     return NULL;
   }
   return forceBox.checkOut();
}

void Patch::unregisterForceDeposit(ComputeID cid, Box<Patch,Results> **const box)
{
   DebugM(4, "unregisterForceDeposit() computeID("<<cid<<")"<<std::endl);
   forceComputeList.del(cid);
   forceBox.checkIn(*box);
   *box = 0;
}

void Patch::positionBoxClosed(void)
{
   //positionPtrBegin = 0;
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
   avgPositionPtrBegin = 0;
   this->boxClosed(3);
}

// void Patch::boxClosed(int box) is virtual

void Patch::positionsReady(int doneMigration)
{
   DebugM(4,"Patch::positionsReady() - patchID(" << patchID <<")"<<std::endl );
   ComputeMap *computeMap = ComputeMap::Object();

   if ( doneMigration ){
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
       AtomMap::Object()->registerIDs(patchID,positionPtrBegin,positionPtrEnd);       
#else
       AtomMap::Object()->registerIDs(patchID,p.begin(),p.end());
#endif
   }

   boxesOpen = 2;
   if ( flags.doMolly ) boxesOpen++;
   _hasNewAtoms = (doneMigration != 0);

#if CMK_VERSION_BLUEGENE
   CmiNetworkProgressAfter (0);
#endif

   // Give all position pickup boxes access to positions
   //positionPtrBegin = p.begin();
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
   positionBox.open(positionPtrBegin);
#else
   positionBox.open(p.begin());
#endif
   if ( flags.doMolly ) {
     //avgPositionPtrBegin = p_avg.begin();
     avgPositionBox.open(avgPositionPtrBegin);
   }

#if CMK_VERSION_BLUEGENE
   CmiNetworkProgressAfter (0);
#endif
   
   // Give all force deposit boxes access to forces
   Force *forcePtr;
   for ( int j = 0; j < Results::maxNumForces; ++j )
    {
      f[j].resize(numAtoms);
      forcePtr = f[j].begin();
      memset (forcePtr, 0, sizeof (Force) * numAtoms);
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


