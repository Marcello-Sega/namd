/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: Common operations for ComputeNonbonded classes
 *
 ***************************************************************************/

#ifndef COMPUTENONBONDEDUTIL_H
#define COMPUTENONBONDEDUTIL_H

#include "NamdTypes.h"

class ComputeNonbondedUtil {

public:

  static void select(void);
  static void (*calcPair)(Position*[2],Force*[2],AtomProperties*[2],int[2]);
  static void (*calcSelf)(Position*,Force*,AtomProperties*,int);
  static void (*calcExcl)(const Position &, const Position &,
		Force &, Force &,
		const AtomProperties &, const AtomProperties &,
		int);

#define DECLARATION
#undef DEFINITON

#define NBPAIR
#undef NBSELF
#undef NBEXCL

#define MODIFY14
#define SWITCHING
#include "ComputeNonbondedBase.h"

#define MODIFY14
#undef SWITCHING
#include "ComputeNonbondedBase.h"

#undef MODIFY14
#define SWITCHING
#include "ComputeNonbondedBase.h"

#undef MODIFY14
#undef SWITCHING
#include "ComputeNonbondedBase.h"


#undef NBPAIR
#define NBSELF
#undef NBEXCL

#define MODIFY14
#define SWITCHING
#include "ComputeNonbondedBase.h"

#define MODIFY14
#undef SWITCHING
#include "ComputeNonbondedBase.h"

#undef MODIFY14
#define SWITCHING
#include "ComputeNonbondedBase.h"

#undef MODIFY14
#undef SWITCHING
#include "ComputeNonbondedBase.h"


#undef NBPAIR
#undef NBSELF
#define NBEXCL

#define MODIFY14
#define SWITCHING
#include "ComputeNonbondedBase.h"

#define MODIFY14
#undef SWITCHING
#include "ComputeNonbondedBase.h"

#undef MODIFY14
#define SWITCHING
#include "ComputeNonbondedBase.h"

#undef MODIFY14
#undef SWITCHING
#include "ComputeNonbondedBase.h"


};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedUtil.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 1996/12/03 21:05:09 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedUtil.h,v $
 * Revision 1.5  1996/12/03 21:05:09  jim
 * added support for exclusion correction computes
 *
 * Revision 1.4  1996/11/21 01:00:15  jim
 * made methods public, got rid of friends
 *
 * Revision 1.3  1996/11/21 00:00:40  jim
 * added select(), calcPair, and calcSelf
 *
 * Revision 1.2  1996/11/20 23:16:39  jim
 * first compiling version of generic nonbonded function
 *
 * Revision 1.1  1996/10/31 22:35:04  jim
 * Initial revision
 *
 *
 ***************************************************************************/

