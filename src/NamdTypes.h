/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdTypes.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/08/16 21:42:58 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdTypes.h,v $
 * Revision 1.2  1996/08/16 21:42:58  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 21:00:37  brunner
 * Initial revision
 *
 ***************************************************************************/

#ifndef NAMDTYPES_H

#define NAMDTYPES_H

typedef double Position;

typedef int PatchID;
typedef int ComputeID;

enum ComputeType
{
  electForceType,
  bondForceType,
  angleForceType,
  dihedralForceType,
  improperForceType
};

enum Boolean
{
  false=0,
  true=1
};

#endif /* NAMDTYPES_H */

