//-*-c++-*-
/***************************************************************************/
/*              (C) Copyright 1996,1997 The Board of Trustees of the       */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Used with single instance owned by main()
 *
 ***************************************************************************/

#ifndef _NAMD_H
#define _NAMD_H

#include "charm++.h"
#include "converse.h"

#include <stdlib.h>

#include "NamdState.h"
#include "Node.h"	// for determining time to compute 1ns
#include "SimParameters.h" // for timesteps

class Node;

class Namd {
public:
  // Constructor for startup node
  Namd(void);             

  // Shut down slaves and clean up
  ~Namd();                

  // read in various input files by invoking
  // proper classes (parameters, molecule etc)
  void startup(char *);   

  // last call of system
  static void namdDone(void) { 
    CkPrintf("==========================================\n");
    Real CPUtime = CmiCpuTimer()-cmiCpuFirstStart;
    Real Walltime = CmiWallTimer()-cmiWallFirstStart;
    CkPrintf("WallClock : %f  CPUTime : %f \n",Walltime,CPUtime);

    // femtoseconds
    SimParameters *params = Node::Object()->simParameters;
    BigReal fs = (params->N - params->firstTimestep) * params->dt;
    // scale to nanoseconds
    BigReal ns = fs / 1000000.0;
    BigReal days = 1.0 / (24.0 * 60.0 * 60.0);
    CkPrintf("Days per ns:  WallClock : %lf  CPUTime : %lf \n",
	days*Walltime/ns, days*CPUtime/ns);
    CkExit(); 
  }

  // Emergency bailout 
  static void die() { CmiAbort("NAMD ABORTING DUE TO Namd::die().\n"); }
  static void startTimer() {
	cmiWallStart = CmiWallTimer(); 
	cmiCpuStart = CmiCpuTimer();
	if ( ! cmiFirstStart ) {
		cmiFirstStart = 1;
		cmiWallFirstStart = cmiWallStart;
		cmiCpuFirstStart = cmiCpuStart;
	}
  }

  static float cmiWallStart;
  static float cmiCpuStart;

private:
  static int cmiFirstStart;
  static float cmiWallFirstStart;
  static float cmiCpuFirstStart;

  Node *node;
  int nodeGroup;
  int workDistribGroup;
  int patchMgrGroup;

  NamdState namdState;

};

#endif /* _NAMD_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Namd.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1012 $	$Date: 1999/05/11 23:56:36 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Namd.h,v $
 * Revision 1.1012  1999/05/11 23:56:36  brunner
 * Changes for new charm version
 *
 * Revision 1.1011  1998/05/25 21:10:36  jim
 * Started using CmiAbort().
 *
 * Revision 1.1010  1998/05/22 00:33:59  jim
 * Fixed final timing when load balancing is used.
 *
 * Revision 1.1009  1998/03/03 23:05:17  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1008  1997/08/18 23:01:17  jim
 * Properly reports "Days per ns" with nonzero firstTimeStep.
 *
 * Revision 1.1007  1997/07/08 15:48:10  milind
 * Made namd2 to work with Origin2000: Again...
 *
 * Revision 1.1006  1997/05/28 19:30:56  nealk
 * Added output to show number of days to compute 1 ns.
 *
 * Revision 1.1005  1997/03/18 18:09:05  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1004  1997/03/06 22:06:05  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1003  1997/03/04 22:37:11  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1002  1997/02/13 23:17:18  ari
 * Fixed a final bug in AtomMigration - numatoms in ComputePatchPair.C not
 * set correctly in atomUpdate()
 *
 * Revision 1.1001  1997/02/13 16:17:14  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1000  1997/02/06 15:58:47  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:56  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:26  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/24 22:00:33  jim
 * Changes for periodic boundary conditions.
 *
 * Revision 1.777  1997/01/17 19:36:30  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.5  1996/12/05 21:37:53  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/16 21:19:34  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 04:39:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 01:54:16  ari
 * Initial revision
 *
 * Revision 1.4  1996/08/03 20:08:09  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/12 16:36:18  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/04 21:42:12  gursoy
 * minor changes
 *
 * Revision 1.1  1996/05/30 20:16:09  brunner
 * Initial revision
 *
 * Revision 1.1  1996/05/30 20:11:09  brunner
 * Initial revision
 *
 ***************************************************************************/
