/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/main.C,v 1.4 1996/08/16 01:54:59 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"

#include "Namd.h"
#include "Communicate.h"
#include "Inform.h"

Communicate *comm;
Inform namdErr("ERROR");
Inform namdWarn("Warning");
Inform namdInfo("Info");
Inform namdDebug("** DEBUG **");


class main : public chare_object
{
public:
  main(int argc, char **argv)
  {
    Namd namd;

    comm = new Communicate();
    if (argc == 2)
	namd.startup(argv[1]);
    else
       CPrintf("main::main() no arguments, exiting\n");

    CharmExit();
  }
};

#include "main.bot.h"
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: main.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1996/08/16 01:54:59 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.C,v $
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

