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
 *	$RCSfile: main.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/08/19 22:05:31 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.h,v $
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

#ifndef MAIN_H
#define MAIN_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

class InitMsg : public comm_object
{
  int x;
};

class NodeInitMsg : public comm_object
{
public:
  int workDistribGroup;
  int patchMgrGroup;
};

#endif /* MAIN_H */



