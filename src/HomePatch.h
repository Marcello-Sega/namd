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

#include "Sequencer.h"

#include "HomePatchTypes.h"
#include "main.h"


class HomePatch : public Patch {
   friend PatchMgr;
   private: // for PatchMgr to use only!!
      HomePatch(PatchID, AtomIDList, PositionList, VelocityList);

   public:

      ~HomePatch();

      void registerProxy(RegisterProxyMsg *);
      void unregisterProxy(UnregisterProxyMsg *);
      void useSequencer(Sequencer *sequencerPtr) {sequencer=sequencerPtr;}
      void runSequencer(int numberOfCycles = 0)
		{ sequencer->run(numberOfCycles); }


   protected:
      virtual void boxClosed(int);

   private:

      PositionList  pInit;
      VelocityList  v; 
      ForceList     f_short;
      ForceList     f_long;

      ProxyList	    proxy;

      Sequencer  *sequencer;



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
 *	$Revision: 1.9 $	$Date: 1996/12/05 23:45:09 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatch.h,v $
 * Revision 1.9  1996/12/05 23:45:09  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/12/05 01:44:16  ari
 * started toward proxy management
 *
 * Revision 1.7  1996/12/01 22:46:11  jim
 * switched to use simParams for number of cycles
 *
 * Revision 1.6  1996/11/30 00:35:51  jim
 * implemented boxClosed(), useSequencer(), runSequencer()
 *
 * Revision 1.5  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/04 21:14:47  jim
 * Moved functionality to Patch
 *
 * Revision 1.3  1996/09/03 22:54:25  ari
 * *** empty log message ***
 *
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
