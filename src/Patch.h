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


// This the base class of homepatches and proxy patches. It maintains
// common functions of these patches. These include managing dependences
// between compute (force) objects and the patch and updating atom map.

class Patch
{
  public:

     Patch(PatchID pd);
     int hasNewAtoms() { return _hasNewAtoms; }
     virtual ~Patch(void) { };

     // methods for use by Compute objects
     Box<Patch,CompAtom>* registerPositionPickup(ComputeID cid, int trans = 13);
     void unregisterPositionPickup(ComputeID cid,
				   Box<Patch,CompAtom>**const box);
     Box<Patch,CompAtom>* registerAvgPositionPickup(ComputeID cid, int trans = 13);
     void unregisterAvgPositionPickup(ComputeID cid,
				   Box<Patch,CompAtom>**const box);
     Box<Patch,Results>* registerForceDeposit(ComputeID cid);
     void unregisterForceDeposit(ComputeID cid, Box<Patch,Results> **const box);

     // methods for use by Sequencer or ProxyManager
     // void positionsReady(void) { positionsReady(0); }
     void positionsReady(int n=0);

     // methods for Box callbacks
     void positionBoxClosed(void);
     void forceBoxClosed(void);
     void avgPositionBoxClosed(void);

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

     ForceList     f[Results::maxNumForces];
     Results	   results;

     OwnerBox<Patch,CompAtom> positionBox;
     ComputeIDList              positionComputeList;
     OwnerBox<Patch,CompAtom> avgPositionBox;
     ComputeIDList              avgPositionComputeList;
     OwnerBox<Patch,Results>    forceBox;
     ComputeIDList              forceComputeList;

     virtual void boxClosed(int /* box */) = 0;
     int boxesOpen;

     int _hasNewAtoms;

#ifdef NODEAWARE_PROXY_SPANNINGTREE    
    //its own children in proxy tree
    int *children;
    int numChild;
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


  private:
  

};


#endif

