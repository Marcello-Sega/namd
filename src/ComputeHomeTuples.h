//-*-c++-*-
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

#ifndef COMPUTEHOMETUPLES_H
#define COMPUTEHOMETUPLES_H

#include "NamdTypes.h"
#include "common.h"
#include "structures.h"
#include "Compute.h"
#include "HomePatch.h"

#include "Box.h"
#include "OwnerBox.h"
#include "PositionBox.h"
#include "PositionOwnerBox.h"
#include "UniqueSet.h"

#include "Namd.h"
#include "Node.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeHomeTuples.h"
#include "PatchMgr.h"
#include "HomePatchList.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "Inform.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"

class TuplePatchElem {
  public:
    PatchID patchID;
    Patch *p;
    PositionBox<Patch> *positionBox;
    Box<Patch,Results> *forceBox;
    Box<Patch,AtomProperties> *atomBox;
    Position *x;
    Results *r;
    Force *f;
    AtomProperties *a;

    int hash() const { return patchID; }

  TuplePatchElem(PatchID p = -1) {
    patchID = p;
    p = NULL;
    positionBox = NULL;
    forceBox = NULL;
    atomBox = NULL;
    x = NULL;
    r = NULL;
    f = NULL;
    a = NULL;
  }

  TuplePatchElem(Patch *p, ComputeID cid) {
    patchID = p->getPatchID();
    this->p = p;
    positionBox = p->registerPositionPickup(cid);
    forceBox = p->registerForceDeposit(cid);
    atomBox = p->registerAtomPickup(cid);
    x = NULL;
    r = NULL;
    f = NULL;
    a = NULL;
  }
    
  ~TuplePatchElem() {};

  int operator==(const TuplePatchElem &a) const {
    return (a.patchID == patchID);
  }

  int operator<(const TuplePatchElem &a) const {
    return (patchID < a.patchID);
  }
};

typedef UniqueSet<TuplePatchElem> TuplePatchList;
typedef UniqueSetIter<TuplePatchElem> TuplePatchListIter;

class AtomMap;
class ReductionMgr;

template <class T, class S> class ComputeHomeTuples : public Compute {

  private:
  
    virtual void loadTuples(void) {
      int numTuples;
      int **tuplesByAtom;
      S *tupleStructs;
    
      T::getMoleculePointers(node->molecule,
		    &numTuples, &tuplesByAtom, &tupleStructs);
    
      char *tupleFlag = new char[numTuples];
      int *tupleStack = new int[numTuples];
      int *nextTuple = tupleStack;
    
      memset((void *)tupleFlag, 0, numTuples*sizeof(char));
    
      // cycle through each patch and gather all tuples
      TuplePatchListIter ai(tuplePatchList);
    
      for ( ai = ai.begin(); ai != ai.end(); ai++ )
      {
        Patch *patch = (*ai).p;
        AtomIDList atomID = patch->getAtomIDList();
        int numAtoms = patch->getNumAtoms();
    
        // cycle through each atom in the patch and load up tuples
        for (int i=0; i < numAtoms; i++)
        {
           /* get list of all tuples for the atom */
           register int *tuples = tuplesByAtom[atomID[i]];
    
           /* cycle through each tuple */
           register int t;
           while((t = *tuples) != -1) {
	     if (!tupleFlag[t])
             {
	       *nextTuple = t;
	       nextTuple++;
               tupleFlag[t] = 1;
             }
	     tuples++;
           }
        }
      }
    
      delete [] tupleFlag;
    
      tupleList.clear();
    
      int *curTuple;
      LocalID aid[T::size];
      for (curTuple = tupleStack; curTuple != nextTuple; curTuple++ ) {
        register int al = *curTuple;
        T t(&tupleStructs[al]);
        register int i;
        aid[0] = atomMap->localID(t.atomID[0]);
        int homepatch = aid[0].pid;
        for (i=1; i < T::size; i++) {
	    aid[i] = atomMap->localID(t.atomID[i]);
	    homepatch = patchMap->downstream(homepatch,aid[i].pid);
        }
        if ( homepatch != notUsed && patchMap->node(homepatch) == CMyPe() ) {
          for (i=0; i < T::size; i++) {
	    t.p[i] = tuplePatchList.find(TuplePatchElem(aid[i].pid));
	    t.localIndex[i] = aid[i].index;
          }
          tupleList.load(t);
        }
      }
      delete [] tupleStack;
    }

    int doLoadTuples;
  
  protected:
  
    UniqueSet<T> tupleList;
    TuplePatchList tuplePatchList;
  
    PatchMap *patchMap;
    AtomMap *atomMap;
    ReductionMgr *reduction;
  
  public:
  
    ComputeHomeTuples(ComputeID c) : Compute(c) {
      patchMap = PatchMap::Object();
      atomMap = AtomMap::Object();
      reduction = ReductionMgr::Object();
      T::registerReductionData(reduction);
      doLoadTuples = false;
    }

    virtual ~ComputeHomeTuples() {
      T::unregisterReductionData(reduction);
    }

    //======================================================================
    // initialize() - Method is invoked only the first time
    // atom maps, patchmaps etc are ready and we are about to start computations
    //======================================================================
    void initialize(void) {
    
      // Gather all HomePatches
      HomePatchList *a = patchMap->homePatchList();
      ResizeArrayIter<HomePatchElem> ai(*a);
    
      // Start with empty list
      tuplePatchList.clear();
    
      for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
        tuplePatchList.add(TuplePatchElem((*ai).patch, cid));
      }
    
      // Gather all proxy patches (neighbors, that is)
      PatchID neighbors[PatchMap::MaxOneOrTwoAway];
    
      for ( ai = ai.begin(); ai != ai.end(); ai++ ) {
        int numNeighbors = patchMap->upstreamNeighbors((*ai).pid,neighbors);
        for ( int i = 0; i < numNeighbors; ++i ) {
          if ( patchMap->node(neighbors[i]) != CMyPe() &&
	       ! tuplePatchList.find(TuplePatchElem(neighbors[i])) ) {
            Patch *patch = patchMap->patch(neighbors[i]);
	    tuplePatchList.add(TuplePatchElem(patch, cid));
          }
        }
      }
      setNumPatches(tuplePatchList.size());
      doLoadTuples = true;

      basePriority = 0;
    }

    //======================================================================
    // atomUpdate() - Method is invoked after anytime that atoms have been
    // changed in patches used by this Compute object.
    //======================================================================
    void atomUpdate(void) {
      doLoadTuples = true;
    }

//-------------------------------------------------------------------
// Routine which is called by enqueued work msg.  It wraps
// actualy Force computation with the apparatus needed
// to get access to atom positions, return forces etc.
//-------------------------------------------------------------------
    void doWork(void) {
      if ( doLoadTuples ) {
        loadTuples();
        doLoadTuples = false;
      }
    
      // Open Boxes - register tFat we are using Positions
      // and will be depositing Forces.
      UniqueSetIter<TuplePatchElem> ap(tuplePatchList);
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->x = ap->positionBox->open();
        ap->a = ap->atomBox->open();
        ap->r = ap->forceBox->open();
        ap->f = ap->r->f[Results::normal];
      } 
    
      BigReal reductionData[T::reductionDataSize];
      for ( int i = 0; i < T::reductionDataSize; ++i ) reductionData[i] = 0;
    
      // take triplet and pass with tuple info to force eval
      UniqueSetIter<T> al(tupleList);
      for (al = al.begin(); al != al.end(); al++ ) {
        al->computeForce(reductionData);
      }
    
      T::submitReductionData(reductionData,reduction,ap.begin()->p->flags.seq);
    
      // Close boxes - i.e. signal we are done with Positions and
      // AtomProperties and that we are depositing Forces
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->positionBox->close(&(ap->x));
        ap->atomBox->close(&(ap->a));
        ap->forceBox->close(&(ap->r));
      }
    }
};


#endif

