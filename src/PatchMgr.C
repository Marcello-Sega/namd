/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************/
/* DESCRIPTION:                                                            */
/*								           */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/PatchMgr.C,v 1.1 1996/08/19 22:07:49 ari Exp $";

#include "PatchMgr.top.h"

#include "main.top.h"
#include "main.h"

#include "NamdTypes.h"
#include "Compute.h"
#include "PatchMgr.h"
#include "HomePatch.h"


PatchMgr::PatchMgr(InitMsg *msg)
{
}

PatchMgr::~PatchMgr()
{
}

void PatchMgr::createPatch(PatchID pid, PositionList& p, VelocityList& v)
{
    
}

void PatchMgr::movePatch(PatchIDList& pid, NodeID nodeID)
{

}




#include "PatchMgr.bot.h"

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PatchMgr.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/19 22:07:49 $
 *
 * REVISION HISTORY:
 *
 * $Log: PatchMgr.C,v $
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 ***************************************************************************/
