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
    PositionBox<Patch> *avgPositionBox;
    Box<Patch,Results> *forceBox;
    Box<Patch,AtomProperties> *atomBox;
    Position *x;
    Position *x_avg;
    Results *r;
    Force *f;
    AtomProperties *a;

    int hash() const { return patchID; }

  TuplePatchElem(PatchID pid = -1) {
    patchID = pid;
    p = NULL;
    positionBox = NULL;
    avgPositionBox = NULL;
    forceBox = NULL;
    atomBox = NULL;
    x = NULL;
    x_avg = NULL;
    r = NULL;
    f = NULL;
    a = NULL;
  }

  TuplePatchElem(Patch *p_param, ComputeID cid) {
    patchID = p_param->getPatchID();
    p = p_param;
    positionBox = p_param->registerPositionPickup(cid);
    avgPositionBox = p_param->registerAvgPositionPickup(cid);
    forceBox = p_param->registerForceDeposit(cid);
    atomBox = p_param->registerAtomPickup(cid);
    x = NULL;
    x_avg = NULL;
    r = NULL;
    f = NULL;
    a = NULL;
  }
    
  ~TuplePatchElem() {};

  int operator==(const TuplePatchElem &elem) const {
    return (elem.patchID == patchID);
  }

  int operator<(const TuplePatchElem &elem) const {
    return (patchID < elem.patchID);
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

      tupleList.clear();

      LocalID aid[T::size];

      // cycle through each patch and gather all tuples
      TuplePatchListIter ai(tuplePatchList);
    
      for ( ai = ai.begin(); ai != ai.end(); ai++ )
      {
        Patch *patch = (*ai).p;
        AtomIDList atomID = patch->getAtomIDList();
        int numAtoms = patch->getNumAtoms();
    
        // cycle through each atom in the patch and load up tuples
        for (int j=0; j < numAtoms; j++)
        {
           /* get list of all tuples for the atom */
           int *curTuple = tuplesByAtom[atomID[j]];
    
           /* cycle through each tuple */
           for( ; *curTuple != -1; ++curTuple) {
             T t(&tupleStructs[*curTuple]);
             register int i;
             aid[0] = atomMap->localID(t.atomID[0]);
             int homepatch = aid[0].pid;
             for (i=1; i < T::size; i++) {
	         aid[i] = atomMap->localID(t.atomID[i]);
	         homepatch = patchMap->downstream(homepatch,aid[i].pid);
             }
             if ( homepatch != notUsed && patchMap->node(homepatch) == CkMyPe() ) {
               for (i=0; i < T::size; i++) {
	         t.p[i] = tuplePatchList.find(TuplePatchElem(aid[i].pid));
	         t.localIndex[i] = aid[i].index;
               }
               tupleList.load(t);
             }
           }
        }
      }
    }

    int doLoadTuples;
  
  protected:
  
    UniqueSet<T> tupleList;
    TuplePatchList tuplePatchList;
  
    PatchMap *patchMap;
    AtomMap *atomMap;
    SubmitReduction *reduction;
  
  public:
  
    ComputeHomeTuples(ComputeID c) : Compute(c) {
      patchMap = PatchMap::Object();
      atomMap = AtomMap::Object();
      reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
      doLoadTuples = false;
    }

    virtual ~ComputeHomeTuples() {
      delete reduction;
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
          if ( patchMap->node(neighbors[i]) != CkMyPe() &&
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
    
      // Open Boxes - register that we are using Positions
      // and will be depositing Forces.
      UniqueSetIter<TuplePatchElem> ap(tuplePatchList);
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->x = ap->positionBox->open();
        if ( ap->p->flags.doMolly ) ap->x_avg = ap->avgPositionBox->open();
        ap->a = ap->atomBox->open();
        ap->r = ap->forceBox->open();
        ap->f = ap->r->f[Results::normal];
      } 
    
      BigReal reductionData[T::reductionDataSize];
      for ( int i = 0; i < T::reductionDataSize; ++i ) reductionData[i] = 0;
      int tupleCount = 0;
    
      // take triplet and pass with tuple info to force eval
      UniqueSetIter<T> al(tupleList);
      for (al = al.begin(); al != al.end(); al++ ) {
        al->computeForce(reductionData);
        tupleCount += 1;
      }
    
      T::submitReductionData(reductionData,reduction);
      reduction->item(T::reductionChecksumLabel) += (BigReal)tupleCount;
      reduction->submit();
    
      // Close boxes - i.e. signal we are done with Positions and
      // AtomProperties and that we are depositing Forces
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->positionBox->close(&(ap->x));
        if ( ap->p->flags.doMolly ) ap->avgPositionBox->close(&(ap->x_avg));
        ap->atomBox->close(&(ap->a));
        ap->forceBox->close(&(ap->r));
      }
    }
};


#endif

