/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

#ifndef _NAMD_H
#define _NAMD_H

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
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:36:30 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Namd.h,v $
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
