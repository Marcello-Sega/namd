/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: main.C,v $
 *	$Author: gursoy $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/02 19:20:13 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.C,v $
 * Revision 1.1  1996/08/02 19:20:13  gursoy
 * Initial revision
 *
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/main.C,v 1.1 1996/08/02 19:20:13 gursoy Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"

class main : public chare_object
{
public:
  main(int argc, char **argv)
  {
    if (argc != 2)
    {
      CPrintf("Usage: %s conf-file\n",argv[0]);
      CharmExit();
    }
    else
    {
      CPrintf("Running :namd2 %s\n",argv[1]);
      CharmExit();
    }
  }
};

#include "main.bot.h"





