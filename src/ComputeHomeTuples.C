/***************************************************************************/
/*           (C) Copyright 1996,1997 The Board of Trustees of the          */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: General template class for computing Bonded Forces 
 *              like Angles, Dihedrals etc..
 *
 ***************************************************************************/

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
#include "Templates/UniqueSet.h"
#include "Templates/UniqueSetIter.h"

//#define DEBUGM
#undef MIN_DEBUG_LEVEL
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

template <class T, class S>
ComputeHomeTuples<T,S>::ComputeHomeTuples(ComputeID c) : Compute(c) {
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();
  reduction = ReductionMgr::Object();
  T::registerReductionData(reduction);
  doLoadTuples = false;
}

template <class T, class S>
ComputeHomeTuples<T,S>::~ComputeHomeTuples()
{
  T::unregisterReductionData(reduction);
}


//===========================================================================
// initialize() - Method is invoked only the first time
// atom maps, patchmaps etc are ready and we are about to start computations
//===========================================================================
template <class T, class S>
void ComputeHomeTuples<T,S>::initialize() {

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
    for ( int i = 0; i < numNeighbors; ++i )
    {
      if ( patchMap->node(neighbors[i]) != CMyPe() &&
	   ! tuplePatchList.find(TuplePatchElem(neighbors[i])) )
      {
        Patch *patch = patchMap->patch(neighbors[i]);
	tuplePatchList.add(TuplePatchElem(patch, cid));
      }
    }
  }

  setNumPatches(tuplePatchList.size());

  doLoadTuples = true;
}

//===========================================================================
// atomUpdate() - Method is invoked after anytime that atoms have been
// changed in patches used by this Compute object.
//===========================================================================
template <class T, class S>
void ComputeHomeTuples<T,S>::atomUpdate() {
  doLoadTuples = true;
}

template <class T, class S>
void ComputeHomeTuples<T,S>::loadTuples() {

  int numTuples;
  int **tuplesByAtom;
  S *tupleStructs;

  T::getMoleculePointers(node->molecule,
		&numTuples, &tuplesByAtom, &tupleStructs);

  char *tupleFlag = new char[numTuples];
  int *tupleStack = new int[numTuples];
  int *nextTuple = tupleStack;

//  register char *cmax = tupleFlag + numTuples;
//  for (register char *c = tupleFlag; c < cmax; *c++ = 0);

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

#if 0
  tupleList.clear();

  if ( node->simParameters->fixedAtomsOn ) {
    Molecule *molecule = node->molecule;
    for (register int i=0; i<numTuples; i++) {
      if (tupleFlag[i]) {
        register int j;
	T t(&tupleStructs[i]);
        register int all_fixed = molecule->is_atom_fixed(t.atomID[0]);
        for ( j = 1; j < T::size && all_fixed; j++ ) {
	  all_fixed = molecule->is_atom_fixed(t.atomID[j]);
        }
        if ( ! all_fixed ) tupleList.load(t);
      }
    }
  } else {
    for (register int i=0; i<numTuples; i++) {
      if (tupleFlag[i]) tupleList.load(T(&tupleStructs[i]));
    }
  }

  delete [] tupleFlag;


  // Resolve all atoms in tupleList to correct PatchList element and index
  // and eliminate tuples we aren't responsible for
  UniqueSetIter<T> al(tupleList);
  UniqueSet<T> tupleList2;


  tupleList.clear();

  int al;
  LocalID aid[T::size];
  for (al = 0; al < numTuples; al++ ) {
//    if (!tupleFlag[al])
//      continue;
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
	/*
        if ( ! (al->p)[i] ) {
 	  iout << iERROR << "ComputeHomeTuples couldn't find patch " 
 	    << aid[i].pid << " for atom " << al->atomID[i] 
 	    << ", aborting.\n" << endi;
 	  Namd::die();
        }
	*/
	t.localIndex[i] = aid[i].index;
      }
      tupleList.load(t);
    }
  }

#endif

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
	/*
        if ( ! (al->p)[i] ) {
 	  iout << iERROR << "ComputeHomeTuples couldn't find patch " 
 	    << aid[i].pid << " for atom " << al->atomID[i] 
 	    << ", aborting.\n" << endi;
 	  Namd::die();
        }
	*/
	t.localIndex[i] = aid[i].index;
      }
      tupleList.load(t);
    }
  }
  delete [] tupleStack;
}


//-------------------------------------------------------------------
// Routine which is called by enqueued work msg.  It wraps
// actualy Force computation with the apparatus needed
// to get access to atom positions, return forces etc.
//-------------------------------------------------------------------
template <class T, class S>
void ComputeHomeTuples<T,S>::doWork() {
  if ( doLoadTuples ) {
    loadTuples();
    doLoadTuples = false;
  }

  DebugM(1, "ComputeHomeTuples::doWork() -- started " << endl );

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
  DebugM(1, "ComputeHomeTuples::doWork() -- done" << endl);
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: ComputeHomeTuples.C,v $
 *      $Author: brunner $  $Locker:  $             $State: Exp $
 *      $Revision: 1.1021 $     $Date: 1997/11/13 00:21:01 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeHomeTuples.C,v $
 * Revision 1.1021  1997/11/13 00:21:01  brunner
 * Revised the tuple-generation procedure.  It now uses an array of tuples
 * to be checked, and flags to prevent duplicate tuples.
 *
 * Revision 1.1020  1997/10/17 20:26:31  jim
 * Sped up tuples, but slowed down fixed atom code.
 * There is a much more efficient way of doing fixed atoms that can
 * be included later.
 *
 * Revision 1.1019  1997/10/17 17:16:46  jim
 * Switched from hash tables to checklists, eliminated special exclusion code.
 *
 * Revision 1.1018  1997/10/06 00:12:28  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 * Revision 1.1017  1997/10/02 22:01:20  jim
 * Moved loadTuples() out of recvProxyAll entry point and into enqueueWork.
 *
 * Revision 1.1016  1997/09/30 16:57:42  jim
 * Fixed bug dealing with atoms on unknown patches.
 *
 * Revision 1.1015  1997/09/28 22:36:49  jim
 * Modified tuple-based computations to not duplicate calculations and
 * only require "upstream" proxies.
 *
 * Revision 1.1014  1997/09/28 10:19:05  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1013  1997/09/22 03:36:00  jim
 * Sped up simulations involving fixed atoms.
 *
 * Revision 1.1012  1997/08/26 16:26:12  jim
 * Revamped prioritites for petter performance and easier changes.
 *
 * Revision 1.1011  1997/04/06 22:44:58  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1010  1997/04/04 23:34:16  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1009  1997/03/18 21:35:25  jim
 * Eliminated fake_seq.  Reductions now use Patch::flags.seq.
 *
 * Revision 1.1008  1997/03/13 22:39:35  jim
 * Fixed some bugs in multiple-force return / full electrostatics.
 *
 * Revision 1.1007  1997/03/12 22:06:35  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1006  1997/03/11 23:46:26  ari
 * Improved ComputeNonbondedExcl loadTuples() by overloading the default
 * template method from ComputeHomeTuples and used the checklist suggested
 * by Jim.  Good performance gain.
 *
 * Revision 1.1005  1997/03/10 17:40:05  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1004  1997/03/04 22:38:16  ari
 * Reworked ResizeArray - much more rational.  Overall about
 * same performance as before (sometimes a little better).
 * Needed tricks to make ResizeArray(Raw) work fast.
 * Clean up of code.
 *
 *
 ***************************************************************************/
