/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PATCH_H
#define PATCH_H

#include "NamdTypes.h"
#include "OwnerBox.h"
#include "Box.h"
#include "UniqueSortedArray.h"
#include "Lattice.h"
#include "PatchTypes.h"

typedef SortedArray<ComputeID> ComputeIDList;

class Compute;
class Sequencer;
class PatchMap;
class AtomMapper;

// This the base class of homepatches and proxy patches. It maintains
// common functions of these patches. These include managing dependences
// between compute (force) objects and the patch and updating atom map.

class Patch
{
  public:

     Patch(PatchID pd);
     int hasNewAtoms() { return _hasNewAtoms; }
     virtual ~Patch();

     // methods for use by Compute objects
     Box<Patch,CompAtom>* registerPositionPickup(ComputeID cid);
     void unregisterPositionPickup(ComputeID cid,
				   Box<Patch,CompAtom>**const box);
     Box<Patch,CompAtom>* registerAvgPositionPickup(ComputeID cid);
     void unregisterAvgPositionPickup(ComputeID cid,
				   Box<Patch,CompAtom>**const box);
     // BEGIN LA
     Box<Patch,CompAtom>* registerVelocityPickup(ComputeID cid);
     void unregisterVelocityPickup(ComputeID cid,
                                  Box<Patch,CompAtom>**const box);
     // END LA

    //begin gbis
    Box<Patch,Real>* registerIntRadPickup(ComputeID cid);
    void unregisterIntRadPickup(ComputeID cid, Box<Patch,Real>**const box);

    Box<Patch,BigReal>* registerPsiSumDeposit(ComputeID cid);
    void unregisterPsiSumDeposit(ComputeID cid, Box<Patch,BigReal>**const box);

    Box<Patch,BigReal>* registerBornRadPickup(ComputeID cid);
    void unregisterBornRadPickup(ComputeID cid, Box<Patch,BigReal>**const box);

    Box<Patch,BigReal>* registerDEdaSumDeposit(ComputeID cid);
    void unregisterDEdaSumDeposit(ComputeID cid,Box<Patch,BigReal> **const box);

    Box<Patch,BigReal>* registerDHdrPrefixPickup(ComputeID cid);
    void unregisterDHdrPrefixPickup(ComputeID cid, Box<Patch,BigReal>**const box);
     //end gbis

     Box<Patch,Results>* registerForceDeposit(ComputeID cid);
     void unregisterForceDeposit(ComputeID cid, Box<Patch,Results> **const box);

     // methods for use by Sequencer or ProxyManager
     // void positionsReady(void) { positionsReady(0); }
     void positionsReady(int n=0);

     // methods for Box callbacks
     void positionBoxClosed(void);
     void forceBoxClosed(void);
     void avgPositionBoxClosed(void);
     // BEGIN LA
     void velocityBoxClosed(void);
     // END LA

     //begin gbis
     void intRadBoxClosed(void);// intrinsic radii
     void psiSumBoxClosed(void);// sum screening 
     void bornRadBoxClosed(void);// born radius
     void dEdaSumBoxClosed(void);// sum dEda contributions
     void dHdrPrefixBoxClosed(void);//dHdr prefix
     void gbisP2Ready();
     void gbisP3Ready();
     //end gbis

     int getNumAtoms() { return numAtoms; }

     // DMK - Atom Separation (water vs. non-water)
     #if NAMD_SeparateWaters != 0
       int getNumWaterAtoms() { return numWaterAtoms; }
     #endif

     int getNumFixedAtoms() { return numFixedAtoms; }  // not updated
     void setNumFixedAtoms(int numFixed) { numFixedAtoms=numFixed; }  // not updated
     PatchID getPatchID() { return patchID; }
     int getNumComputes() { return positionComputeList.size(); }

     CompAtomExt* getCompAtomExtInfo() { return pExt.begin(); }

     Lattice &lattice;
     Flags flags;

  protected:

     const PatchID patchID;
     int           numAtoms;
     int           numFixedAtoms;
     CompAtomList  p;
     CompAtomList  p_avg;
     // BEGIN LA
     CompAtomList  v;
     // END LA

     AtomMapper *atomMapper;

     // begin gbis
     RealList intRad;
     BigRealList psiSum;
     BigRealList psiFin;
     BigRealList bornRad;
     BigRealList dHdrPrefix;
     BigRealList dEdaSum;
     // end gbis

     // DMK - Atom Separation (water vs. non-water)
     #if NAMD_SeparateWaters != 0
       int numWaterAtoms;  // Set numWaters to the number of water atoms at
                           //   the lead of the atoms list.  If numWaters is
                           //   set to -1, this should indicate that
                           //   atoms has not been separated yet.
     #endif

     CompAtomExtList pExt;

#ifdef REMOVE_PROXYDATAMSG_EXTRACOPY
     //1. Those fields are declared for reusing position info
     //inside the ProxyDataMsg msg at every step so that the
     //extra copy is avoided.
     //Regarding the CompAtomExt list inside the msg of ProxyAllMsg type
     //we cannot avoid the copy in the current scheme because this information
     //will be lost as the msg will be deleted at the next timestep. But the
     //overhead is amortized among the steps that atoms don't migrate
     //2. positionPtrBegin is better to be made 32-byte aligned so we could
     // have better cache performance in the force calculation part. This
     // is especially needed for BG/L machine.
     // --Chao Mei
     CompAtom      *positionPtrBegin;
     CompAtom      *positionPtrEnd;     
#endif
     CompAtom      *avgPositionPtrBegin;
     CompAtom      *avgPositionPtrEnd;

     // BEGIN LA
     CompAtom      *velocityPtrBegin;
     CompAtom      *velocityPtrEnd;
     // END LA

     ForceList     f[Results::maxNumForces];
     Results	   results;

     OwnerBox<Patch,CompAtom> positionBox;
     ComputeIDList              positionComputeList;
     OwnerBox<Patch,CompAtom> avgPositionBox;
     ComputeIDList              avgPositionComputeList;
     // BEGIN LA
     OwnerBox<Patch,CompAtom> velocityBox;
     ComputeIDList              velocityComputeList;
     // END LA

     //begin gbis
     OwnerBox<Patch,Real> intRadBox;
     ComputeIDList           intRadComputeList;
     OwnerBox<Patch,BigReal> psiSumBox;
     ComputeIDList           psiSumComputeList;
     OwnerBox<Patch,BigReal> bornRadBox;
     ComputeIDList           bornRadComputeList;
     OwnerBox<Patch,BigReal> dEdaSumBox;
     ComputeIDList           dEdaSumComputeList;
     OwnerBox<Patch,BigReal> dHdrPrefixBox;
     ComputeIDList           dHdrPrefixComputeList;
     //end gbis

     OwnerBox<Patch,Results>    forceBox;
     ComputeIDList              forceComputeList;

     virtual void boxClosed(int /* box */) = 0;
     int boxesOpen;

     int _hasNewAtoms;

#ifdef NODEAWARE_PROXY_SPANNINGTREE    
    //its own children in proxy tree
    #ifdef USE_NODEPATCHMGR
    //the immediate children (in terms of node id) also cotains two parts
    //as the above variable shows
    //If this patch has proxies residing on the same node, then the last entry
    //of "nodeChildren" stores this node id. It is same with that variable
    //in ProxyPatch      
    //If this proxy resides on node A, then the last entry
    //of "nodeChildren" has to be A
    //It is same with that variable in HomePatch
    int *nodeChildren;
    int numNodeChild;
    #endif
#endif
	int *child;
	int nChild;

  private:
  

};


#endif

