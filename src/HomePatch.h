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

#ifndef HOMEPATCH_H
#define HOMEPATCH_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"
#include "Patch.h"
#include "Templates/UniqueSortedArray.h"


class LocalAtomID {
public:
    AtomID atomID;
    int index;
    LocalAtomID(AtomID a, int i) : atomID(a), index(i) {};
    LocalAtomID() {};
    ~LocalAtomID() {};
    int operator < (const LocalAtomID &a) const {
	return (atomID < a.atomID);
    }
    int operator== (const LocalAtomID &a) const {
       return (atomID == a.atomID);
    }
};


class LocalIndex : public UniqueSortedArray<LocalAtomID> { };

class HomePatch : public Patch {
   friend PatchMgr;

   private:
      PatchID       patchID;
      PositionList  p;
      PositionList  pBegin;
      VelocityList  v; 
      ForceList     f;
      ForceList     f_short;
      ForceList     f_long;
      AtomIDList    atomIDList;
      LocalIndex    localIndex;

      //Sequencer  *sequencer;
      // keeps track of remaining proxies that i am waiting to updat my data
      //int        proxyCounter; 
      // declaration for 
      // local computations i.e., computations that needs the local data only
      // for example 
      //       -- harmonic constraints
      //       -- shake 
      //       etc

      inline int beginCycle() { return 0; }

   public:

      HomePatch();
      HomePatch(PatchID pd, AtomIDList al, 
	PositionList pl, VelocityList vl) : 
	  patchID(pd), atomIDList(al), p(pl), pBegin(pl), v(vl) {
	if (atomIDList.size() != p.size() || atomIDList.size() != v.size()) {
	  CPrintf("HomePatch::HomePatch(...) : Whoa - got different # of\
	  atom coordinates, id's or velocities!\n");
	}
	AtomIDListIter a(atomIDList);
	int i = 0;
	for ( a = a.begin(); a != a.end(); a++ ) {
	  LocalAtomID la(*a, i++);
	  localIndex.load(la);
	}
	localIndex.sort();
	localIndex.uniq();
      }
	
      ~HomePatch();


      // void use_sequencer(Sequencer *sequencerPtr) {sequencer=sequencerPtr;}

      // 
      /*
      void prepare_for_next_cycle();
      void prepare_for_next_step();
      */
      
      // calculations that depends on remote data
      /*
      void compute_f_long();
      void compute_f_short();
      */


      // local calculations
      /*
      void update_f_at_cycle_begin();
      void update_f_at_step(int);
      void advance_x();
      void update_f_at_cycle_end();
      void update_v();
      void output();
      void f_short_done();
      void f_long_done();
      void atom_redist_data(int, int *);
      */

      // pack and unpack functions

      /*
      void dispose(char *&);
      void updateData(char *&);
      void packInitData(char *&);
      void packData(char *&);
      */
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: HomePatch.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/08/29 00:50:42 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatch.h,v $
 * Revision 1.2  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.11  1996/07/15 22:57:05  gursoy
 * *** empty log message ***
 *
 * Revision 1.10  1996/07/15 21:18:01  gursoy
 * *** empty log message ***
 *
 * Revision 1.9  1996/07/10 16:19:52  gursoy
 * *** empty log message ***
 *
 * Revision 1.8  1996/07/08 21:32:38  gursoy
 * ,
 *
 * Revision 1.7  1996/06/12 16:34:46  brunner
 * *** empty log message ***
 *
 * Revision 1.6  1996/06/11 22:36:35  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/06/11 20:07:22  gursoy
 * *** empty log message ***
 *
 * Revision 1.4  1996/06/10 22:04:14  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/10 20:31:57  gursoy
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/10 20:30:00  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/05/30 21:31:36  gursoy
 * Initial revision
 *
 ***************************************************************************/
