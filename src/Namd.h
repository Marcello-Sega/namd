//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

#ifndef _NAMD_H
#define _NAMD_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdState.h"

class Node;

class Namd
{

public:
    // Constructor for startup node
  Namd(void);             
    // Shut down slaves and clean up
  ~Namd();                
    // read in various input files by invoking
    // proper classes (parameters, molecule etc)
  void startup(char *);   
  static void die() { CharmExit(); }

private:
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
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:58:47 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Namd.h,v $
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
