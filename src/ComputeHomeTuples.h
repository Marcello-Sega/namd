/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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

#include "Node.h"
#include "SimParameters.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeHomeTuples.h"
#include "PatchMgr.h"
#include "HomePatchList.h"
#include "Molecule.h"
#include "Parameters.h"
#include "ReductionMgr.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"

class TuplePatchElem {
  public:
    PatchID patchID;
    Patch *p;
    PositionBox<Patch> *positionBox;
    PositionBox<Patch> *avgPositionBox;
    Box<Patch,Results> *forceBox;
    CompAtom *x;
    CompAtom *x_avg;
    Results *r;
    Force *f;

    int hash() const { return patchID; }

  TuplePatchElem(PatchID pid = -1) {
    patchID = pid;
    p = NULL;
    positionBox = NULL;
    avgPositionBox = NULL;
    forceBox = NULL;
    x = NULL;
    x_avg = NULL;
    r = NULL;
    f = NULL;
  }

  TuplePatchElem(Patch *p_param, ComputeID cid) {
    patchID = p_param->getPatchID();
    p = p_param;
    positionBox = p_param->registerPositionPickup(cid);
    avgPositionBox = p_param->registerAvgPositionPickup(cid);
    forceBox = p_param->registerForceDeposit(cid);
    x = NULL;
    x_avg = NULL;
    r = NULL;
    f = NULL;
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

template <class T, class S, class P> class ComputeHomeTuples : public Compute {

  protected:
  
    virtual void loadTuples(void) {
      int numTuples;
      int32 **tuplesByAtom;
      /* const (need to propagate const) */ S *tupleStructs;
      const P *tupleValues;
    
      T::getMoleculePointers(node->molecule,
		    &numTuples, &tuplesByAtom, &tupleStructs);
      T::getParameterPointers(node->parameters, &tupleValues);

      tupleList.clear();

      LocalID aid[T::size];

      Real invLesFactor = node->simParameters->lesOn ? 
                          1.0/node->simParameters->lesFactor :
                          1.0;

      // cycle through each patch and gather all tuples
      TuplePatchListIter ai(tuplePatchList);
    
      for ( ai = ai.begin(); ai != ai.end(); ai++ )
      {
        CompAtom *atom = (*ai).x;
        Patch *patch = (*ai).p;
        int numAtoms = patch->getNumAtoms();
    
        // cycle through each atom in the patch and load up tuples
        for (int j=0; j < numAtoms; j++)
        {
           /* get list of all tuples for the atom */
           int32 *curTuple = tuplesByAtom[atom[j].id];
    
           /* cycle through each tuple */
           for( ; *curTuple != -1; ++curTuple) {
             T t(&tupleStructs[*curTuple],tupleValues);
             
             register int i;
             aid[0] = atomMap->localID(t.atomID[0]);
             int homepatch = aid[0].pid;
             int samepatch = 1;
             int has_les = node->molecule->get_fep_type(t.atomID[0]);
             for (i=1; i < T::size; i++) {
	         aid[i] = atomMap->localID(t.atomID[i]);
	         samepatch = samepatch && ( homepatch == aid[i].pid );
                 has_les |= node->molecule->get_fep_type(t.atomID[i]);
             }
             if ( samepatch ) continue;
             t.scale = has_les ? invLesFactor : 1;
             for (i=1; i < T::size; i++) {
	         homepatch = patchMap->downstream(homepatch,aid[i].pid);
             }
             if ( homepatch != notUsed && isBasePatch[homepatch] ) {
               for (i=0; i < T::size; i++) {
	         TuplePatchElem *p;
	         t.p[i] = p = tuplePatchList.find(TuplePatchElem(aid[i].pid));
	         if ( ! p ) {
	           iout << iWARN << "Tuple " << *curTuple << " with atoms ";
	           int erri;
	           for( erri = 0; erri < T::size; erri++ ) {
	             iout << t.atomID[erri] << "(" <<  aid[erri].pid << ") ";
	           }
	           iout << "missing patch " << aid[i].pid << "\n" << endi;
	           
	           NAMD_die("Patch needed for tuple is missing.\n");
	         }
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
    char *isBasePatch;
  
    ComputeHomeTuples(ComputeID c) : Compute(c) {
      patchMap = PatchMap::Object();
      atomMap = AtomMap::Object();
      reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
      doLoadTuples = false;
      isBasePatch = 0;
    }

    ComputeHomeTuples(ComputeID c, PatchIDList pids) : Compute(c) {
      patchMap = PatchMap::Object();
      atomMap = AtomMap::Object();
      reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
      doLoadTuples = false;
      int nPatches = patchMap->numPatches();
      isBasePatch = new char[nPatches];
      int i;
      for (i=0; i<nPatches; ++i) { isBasePatch[i] = 0; }
      for (i=0; i<pids.size(); ++i) { isBasePatch[pids[i]] = 1; }
    }

  public:
  
    virtual ~ComputeHomeTuples() {
      delete reduction;
      delete [] isBasePatch;
    }

    //======================================================================
    // initialize() - Method is invoked only the first time
    // atom maps, patchmaps etc are ready and we are about to start computations
    //======================================================================
    virtual void initialize(void) {
    
      // Start with empty list
      tuplePatchList.clear();
    
      int nPatches = patchMap->numPatches();
      int pid;
      for (pid=0; pid<nPatches; ++pid) {
        if ( isBasePatch[pid] ) {
          Patch *patch = patchMap->patch(pid);
	  tuplePatchList.add(TuplePatchElem(patch, cid));
        }
      }
    
      // Gather all proxy patches (neighbors, that is)
      PatchID neighbors[PatchMap::MaxOneOrTwoAway];
    
      for (pid=0; pid<nPatches; ++pid) if ( isBasePatch[pid] ) {
        int numNeighbors = patchMap->upstreamNeighbors(pid,neighbors);
        for ( int i = 0; i < numNeighbors; ++i ) {
          if ( ! tuplePatchList.find(TuplePatchElem(neighbors[i])) ) {
            Patch *patch = patchMap->patch(neighbors[i]);
	    tuplePatchList.add(TuplePatchElem(patch, cid));
          }
        }
      }
      setNumPatches(tuplePatchList.size());
      doLoadTuples = true;

      basePriority = 1;
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
    virtual void doWork(void) {

      // Open Boxes - register that we are using Positions
      // and will be depositing Forces.
      UniqueSetIter<TuplePatchElem> ap(tuplePatchList);
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->x = ap->positionBox->open();
        if ( ap->p->flags.doMolly ) ap->x_avg = ap->avgPositionBox->open();
        ap->r = ap->forceBox->open();
        ap->f = ap->r->f[Results::normal];
      } 
    
      if ( doLoadTuples ) {
        loadTuples();
        doLoadTuples = false;
      }
    
      BigReal reductionData[T::reductionDataSize];
      for ( int i = 0; i < T::reductionDataSize; ++i ) reductionData[i] = 0;
      int tupleCount = 0;
    
      // take triplet and pass with tuple info to force eval
      UniqueSetIter<T> al(tupleList);
      if ( Node::Object()->simParameters->commOnly ) {
        for (al = al.begin(); al != al.end(); al++ ) {
          tupleCount += 1;
        }
      } else {
        for (al = al.begin(); al != al.end(); al++ ) {
          al->computeForce(reductionData);
          tupleCount += 1;
        }
      }
    
      T::submitReductionData(reductionData,reduction);
      reduction->item(T::reductionChecksumLabel) += (BigReal)tupleCount;
      reduction->submit();
    
      // Close boxes - i.e. signal we are done with Positions and
      // AtomProperties and that we are depositing Forces
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->positionBox->close(&(ap->x));
        if ( ap->p->flags.doMolly ) ap->avgPositionBox->close(&(ap->x_avg));
        ap->forceBox->close(&(ap->r));
      }
    }
};


#endif

