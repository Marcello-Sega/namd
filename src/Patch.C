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

Patch::~Patch() {
  delete atomMapper;
}

Patch::Patch(PatchID pd) :
   lattice(flags.lattice),
   patchID(pd), numAtoms(0), numFixedAtoms(0),
   avgPositionPtrBegin(0), avgPositionPtrEnd(0),
   velocityPtrBegin(0), velocityPtrEnd(0),	// BEGIN LA, END LA
   positionBox(this,&Patch::positionBoxClosed,pd,0),
   avgPositionBox(this,&Patch::avgPositionBoxClosed,pd,3),
   velocityBox(this,&Patch::velocityBoxClosed,pd,4), // BEGIN LA, END LA
   psiSumBox(this,&Patch::psiSumBoxClosed,pd,5), // begin gbis
   intRadBox(this,&Patch::intRadBoxClosed,pd,6),
   bornRadBox(this,&Patch::bornRadBoxClosed,pd,7),
   dEdaSumBox(this,&Patch::dEdaSumBoxClosed,pd,8),
   dHdrPrefixBox(this,&Patch::dHdrPrefixBoxClosed,pd,9), //end gbis
   forceBox(this,&Patch::forceBoxClosed,pd,1),
   boxesOpen(0), _hasNewAtoms(0)

   // DMK - Atom Separation (water vs. non-water)
   #if NAMD_SeparateWaters != 0
     ,numWaterAtoms(-1)
   #endif
{
  //CkPrintf("GBIS: PatchCreated\n");
#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
    positionPtrBegin = 0;
    positionPtrEnd = 0;
#endif

	nChild = 0;
	child = NULL;
#ifdef NODEAWARE_PROXY_SPANNINGTREE
  #ifdef USE_NODEPATCHMGR
  nodeChildren = NULL;
  numNodeChild = 0;
  #endif
#endif

  lattice = Node::Object()->simParameters->lattice;
  atomMapper = new AtomMapper(pd);
}

Box<Patch,CompAtom>* Patch::registerPositionPickup(ComputeID cid)
{
   //DebugM(4, "registerPositionPickupa("<<patchID<<") from " << cid << "\n");
   if (positionComputeList.add(cid) < 0)
   {
     DebugM(7, "registerPositionPickup() failed for cid " << cid << std::endl);
     return NULL;
   }
   return positionBox.checkOut(cid);
}

void Patch::unregisterPositionPickup(ComputeID cid, Box<Patch,CompAtom> **const box)
{
   DebugM(4, "UnregisterPositionPickup from " << cid << "\n");
   positionComputeList.del(cid);
   positionBox.checkIn(*box);
   *box = 0;
}

Box<Patch,CompAtom>* Patch::registerAvgPositionPickup(ComputeID cid)
{
   //DebugM(4, "registerAvgPositionPickup("<<patchID<<") from " << cid << "\n");
   return avgPositionBox.checkOut(cid);
}

void Patch::unregisterAvgPositionPickup(ComputeID cid, Box<Patch,CompAtom> **const box)
{
   DebugM(4, "UnregisterAvgPositionPickup from " << cid << "\n");
   avgPositionBox.checkIn(*box);
   *box = 0;
}

// BEGIN LA
Box<Patch,CompAtom>* Patch::registerVelocityPickup(ComputeID cid)
{
   //DebugM(4, "registerVelocityPickup("<<patchID<<") from " << cid << "\n");
   return velocityBox.checkOut(cid);
}

void Patch::unregisterVelocityPickup(ComputeID cid, Box<Patch,CompAtom> **const box)
{
   DebugM(4, "UnregisterVelocityPickup from " << cid << "\n");
   velocityBox.checkIn(*box);
   *box = 0;
}
// END LA

//begin gbis
//deposit, not pickup
Box<Patch,BigReal>* Patch::registerPsiSumDeposit(ComputeID cid) {

  if (psiSumComputeList.add(cid) < 0) {
    DebugM(7, "registerPsiSumDeposit() failed for cid " << cid << std::endl);
    DebugM(7, "  size of psiSumCompueList " << psiSumComputeList.size() << std::endl);
     return NULL;
  }
  return psiSumBox.checkOut(cid);
}

void Patch::unregisterPsiSumDeposit(ComputeID cid,Box<Patch,BigReal> **const box) {
  psiSumComputeList.del(cid);
  psiSumBox.checkIn(*box);
  *box = 0;
}
Box<Patch,Real>* Patch::registerIntRadPickup(ComputeID cid) {
  return intRadBox.checkOut(cid);
}
void Patch::unregisterIntRadPickup(ComputeID cid,Box<Patch,Real> **const box) {
  intRadBox.checkIn(*box);
  *box = 0;
}

Box<Patch,BigReal>* Patch::registerBornRadPickup(ComputeID cid) {
  return bornRadBox.checkOut(cid);
}
void Patch::unregisterBornRadPickup(ComputeID cid,Box<Patch,BigReal> **const box) {
  bornRadBox.checkIn(*box);
  *box = 0;
}

Box<Patch,BigReal>* Patch::registerDEdaSumDeposit(ComputeID cid) {
  if (dEdaSumComputeList.add(cid) < 0) {
    DebugM(7, "registerDEdaSumDeposit() failed for cid " << cid << std::endl);
    DebugM(7, "  size of dEdaSumCompueList " << dEdaSumComputeList.size() << std::endl);
     return NULL;
  }
  return dEdaSumBox.checkOut(cid);
}
void Patch::unregisterDEdaSumDeposit(ComputeID cid,Box<Patch,BigReal> **const box){
  dEdaSumComputeList.del(cid);
  dEdaSumBox.checkIn(*box);
  *box = 0;
}

Box<Patch,BigReal>* Patch::registerDHdrPrefixPickup(ComputeID cid)
{
  return dHdrPrefixBox.checkOut(cid);
}
void Patch::unregisterDHdrPrefixPickup(ComputeID cid,Box<Patch,BigReal> **const box) {
  dHdrPrefixBox.checkIn(*box);
  *box = 0;
}
//end gbis

Box<Patch,Results>* Patch::registerForceDeposit(ComputeID cid)
{
   if (forceComputeList.add(cid) < 0)
   {
     DebugM(7, "registerForceDeposit() failed for cid " << cid << std::endl);
     DebugM(7, "  size of forceCompueList " << forceComputeList.size() << std::endl);
     return NULL;
   }
   return forceBox.checkOut(cid);
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

// BEGIN LA
void Patch::velocityBoxClosed(void)
{
   DebugM(4, "patchID("<<patchID<<") velocityBoxClosed! call\n");
   velocityPtrBegin = 0;
   this->boxClosed(4);	// ?? Don't know about number
}
// END LA

// void Patch::boxClosed(int box) is virtual

// begin gbis
void Patch::psiSumBoxClosed(void) {
  this->boxClosed(5);
}
void Patch::intRadBoxClosed(void) {
   //dHdrPrefixPtr = 0;
   this->boxClosed(6);
}
void Patch::bornRadBoxClosed(void) {
   //bornRadPtr = 0;
   this->boxClosed(7);
}
void Patch::dEdaSumBoxClosed(void) {
   //dEdaSumPtr = 0;
   this->boxClosed(8);
}
void Patch::dHdrPrefixBoxClosed(void) {
   //dHdrPrefixPtr = 0;
   this->boxClosed(9);
}
// end gbis

void Patch::positionsReady(int doneMigration)
{
   DebugM(4,"Patch::positionsReady() - patchID(" << patchID <<")"<<std::endl );
   ComputeMap *computeMap = ComputeMap::Object();

   if ( doneMigration ){
// #ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
//       AtomMap::Object()->registerIDs(patchID,positionPtrBegin,positionPtrEnd);       
// #else
       atomMapper->registerIDsCompAtomExt(pExt.begin(),pExt.end());
// #endif
   }

   boxesOpen = 2;
   if ( flags.doMolly ) boxesOpen++;
   // BEGIN LA
   if (flags.doLoweAndersen) {
       DebugM(4, "Patch::positionsReady, flags.doMolly = " << flags.doMolly << "\n");
       boxesOpen++;
   }
   // END LA
   _hasNewAtoms = (doneMigration != 0);

#if CMK_BLUEGENEL
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
   
   // BEGIN LA
   if (flags.doLoweAndersen) {
       velocityBox.open(velocityPtrBegin);
   }
   // END LA
   // begin gbis
    if (flags.doGBIS) {
      boxesOpen += 5;
      //intRad should already be taken care of
      intRadBox.open(intRad.begin());
      psiSum.resize(numAtoms);//resize array
      psiSum.setall(0);
      psiSumBox.open(psiSum.begin());
      psiFin.resize(numAtoms);//has no box
      psiFin.setall(0);
      bornRad.resize(numAtoms);
      bornRad.setall(0);
      bornRadBox.open(bornRad.begin());
      dEdaSum.resize(numAtoms);//resize array
      dEdaSum.setall(0);
      dEdaSumBox.open(dEdaSum.begin());
      dHdrPrefix.resize(numAtoms);
      dHdrPrefix.setall(0);
      dHdrPrefixBox.open(dHdrPrefix.begin());
    }
   // end gbis

#if CMK_BLUEGENEL
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

// begin gbis

void Patch::gbisP2Ready() {
 ComputeMap *computeMap = ComputeMap::Object();
 ComputeIDListIter cid(positionComputeList);

  int compute_count = 0;
  int seq = flags.sequence;
  for(cid = cid.begin(); cid != cid.end(); cid++) {
    if ( computeMap->compute(*cid)->type() == computeNonbondedSelfType ||
         computeMap->compute(*cid)->type() == computeNonbondedPairType) {
      compute_count++;
      computeMap->compute(*cid)->gbisP2PatchReady(patchID,seq);
    }
  }
}

void Patch::gbisP3Ready() {

  ComputeMap *computeMap = ComputeMap::Object();
  ComputeIDListIter cid(positionComputeList);

  int compute_count = 0;
  int seq = flags.sequence;
  for(cid = cid.begin(); cid != cid.end(); cid++) {
    if ( computeMap->compute(*cid)->type() == computeNonbondedSelfType ||
         computeMap->compute(*cid)->type() == computeNonbondedPairType) {
      compute_count++;
      computeMap->compute(*cid)->gbisP3PatchReady(patchID,seq);
    }
  }
}

//end gbis
