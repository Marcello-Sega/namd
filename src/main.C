/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/main.C,v 1.778 1997/01/28 00:31:37 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"

// Needed for namd.1.X components
#include "Namd.h"
#include "CommunicateConverse.h"
#include "Inform.h"

// Needed for namd.1.X components
Communicate *comm = NULL;
Inform namdErr("ERROR");
Inform namdWarn("Warning");
Inform namdInfo("Info");
Inform namdDebug("** DEBUG **");

class main : public chare_object
{
public:
  main(int argc, char **argv)
  {
    Namd *namd = new Namd;
    CPrintf("main::main() - Namd::Namd() should have been invoked\n");


    if (argc >= 2)
	namd->startup(argv[argc-1]);
    else
       CPrintf("main::main() no arguments, exiting\n");
    CPrintf("main() - leaving - Charm should queue up messages now!\n");
  }
};

#include "main.bot.h"
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: main.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.778 $	$Date: 1997/01/28 00:31:37 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.C,v $
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

