/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTESELFTUPLES_H
#define COMPUTESELFTUPLES_H

#include "ComputeHomeTuples.h"
#include "LdbCoordinator.h"
    
template <class T, class S, class P> class ComputeSelfTuples :
	public ComputeHomeTuples<T,S,P> {

  private:
  
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
      // There should be only one!
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
             if ( samepatch ) {
               t.scale = has_les ? invLesFactor : 1;
	       TuplePatchElem *p;
	       p = tuplePatchList.find(TuplePatchElem(homepatch));
               for (i=0; i < T::size; i++) {
	         t.p[i] = p;
	         t.localIndex[i] = aid[i].index;
               }
               tupleList.load(t);
             }
           }
        }
      }
    }

    PatchID patchID;

  public:

    ComputeSelfTuples(ComputeID c, PatchID p) : ComputeHomeTuples<T,S,P>(c) {
      patchID = p;
    }

    virtual ~ComputeSelfTuples() {
      UniqueSetIter<TuplePatchElem> ap(tuplePatchList);
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->p->unregisterPositionPickup(cid,&(ap->positionBox));
        ap->p->unregisterAvgPositionPickup(cid,&(ap->avgPositionBox));
        ap->p->unregisterForceDeposit(cid,&(ap->forceBox));
      }
    }


    //======================================================================
    // initialize() - Method is invoked only the first time
    // atom maps, patchmaps etc are ready and we are about to start computations
    //======================================================================
    virtual void initialize(void) {
    
      // Start with empty list
      tuplePatchList.clear();
    
      tuplePatchList.add(TuplePatchElem(ComputeHomeTuples<T,S,P>::patchMap->patch(patchID), cid));
    
      setNumPatches(tuplePatchList.size());
      doLoadTuples = true;

      int myNode = CkMyPe();
      if ( PatchMap::Object()->node(patchID) != myNode )
      {
        basePriority = 64 + patchID % 64;
      }
      else
      {
        basePriority = 2 * 64 + (patchID % 64);
      }
    }

    void doWork(void) {
      LdbCoordinator::Object()->startWork(cid,0); // Timestep not used

      ComputeHomeTuples<T,S,P>::doWork();

      LdbCoordinator::Object()->endWork(cid,0); // Timestep not used
    }

#if 0
//-------------------------------------------------------------------
// Routine which is called by enqueued work msg.  It wraps
// actualy Force computation with the apparatus needed
// to get access to atom positions, return forces etc.
//-------------------------------------------------------------------
    void doWork(void) {

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
#endif

};


#endif

