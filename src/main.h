/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef MAIN_H
#define MAIN_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

class InitMsg : public comm_object
{
  int x;
};

class MovePatchDoneMsg : public comm_object
{
  int x;
};

class PatchMgrInitMsg : public comm_object
{
public:
  int workDistribGroup;
};

class NodeInitMsg : public comm_object
{
public:
  int workDistribGroup;
  int patchMgrGroup;
};


class Compute;

class LocalWorkMsg : public comm_object
{
public:
  Compute *compute;
};

#endif /* MAIN_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: main.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.7 $	$Date: 1996/09/15 20:43:11 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.h,v $
 * Revision 1.7  1996/09/15 20:43:11  jim
 * fixed missing semicolon
 *
 * Revision 1.6  1996/09/10 04:45:06  ari
 * added LocalWorkMsg
 *
 * Revision 1.5  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/22 16:15:07  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/19 22:05:31  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 20:52:30  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/02 19:20:13  gursoy
 * Initial revision
 *
 *
 ***************************************************************************/
