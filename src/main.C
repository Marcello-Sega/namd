/***************************************************************************/
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/main.C,v 1.1004 1997/10/01 16:47:04 milind Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"

// Needed for namd.1.X components
#include "Namd.h"
#include "Communicate.h"
#include "Inform.h"

// Needed for namd.1.X components
Communicate *comm = NULL; // Should be initialized in Node::node()

Inform namdErr("ERROR");
Inform namdWarn("Warning");
Inform namdInfo("Info");
Inform namdDebug("** DEBUG **");

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

class main : public chare_object
{
public:
  main(int argc, char **argv)
  {
    // Namd object is only on Pe(0)
    Namd *namd = new Namd;

    if (argc >= 2) {
	namd->startup(argv[argc-1]);
    }
    else {
       CPrintf("main::main() no arguments, exiting\n");
    }
    DebugM(1, "main() - leaving - Charm should queue up messages now!\n");
  }
};

#include "main.bot.h"
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: main.C,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1004 $	$Date: 1997/10/01 16:47:04 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.C,v $
 * Revision 1.1004  1997/10/01 16:47:04  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.1003  1997/03/19 11:55:03  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1002  1997/03/04 22:37:20  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1001  1997/02/26 16:53:21  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1000  1997/02/06 15:59:35  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:37  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:37:13  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.11  1996/12/27 22:21:54  nealk
 * Added some debugging code.
 *
 * Revision 1.10  1996/12/11 00:04:23  milind
 * *** empty log message ***
 *
 * Revision 1.9  1996/12/10 00:13:12  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/12/06 19:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/11/07 23:08:38  ari
 * *** empty log message ***
 *
 * Revision 1.6  1996/11/05 16:59:58  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/08/16 04:39:46  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/16 01:54:59  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/15 20:32:14  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/06 20:38:38  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/02 19:20:13  gursoy
 * Initial revision
 *
 *
 ***************************************************************************/

