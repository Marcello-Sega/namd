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

#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "Node.h"

void ComputeNonbondedUtil::select(void)
{
#undef DECLARATION
#undef DEFINITION

  if ( Node::Object()->simParameters->exclude == SCALED14 )
#define MODIFY14
  {
    if ( Node::Object()->simParameters->switchingActive )
#define SWITCHING
    {
#define NBPAIR
#undef NBSELF
#undef NBEXCL
#include "ComputeNonbondedHack.h"
      calcPair = NAME;
#undef NBPAIR
#define NBSELF
#undef NBEXCL
#include "ComputeNonbondedHack.h"
      calcSelf = NAME;
    }
    else
#undef SWITCHING
    {
#define NBPAIR
#undef NBSELF
#undef NBEXCL
#include "ComputeNonbondedHack.h"
      calcPair = NAME;
#undef NBPAIR
#define NBSELF
#undef NBEXCL
#include "ComputeNonbondedHack.h"
      calcSelf = NAME;
    }
  }
  else
#undef MODIFY14
  {
    if ( Node::Object()->simParameters->switchingActive )
#define SWITCHING
    {
#define NBPAIR
#undef NBSELF
#undef NBEXCL
#include "ComputeNonbondedHack.h"
      calcPair = NAME;
#undef NBPAIR
#define NBSELF
#undef NBEXCL
#include "ComputeNonbondedHack.h"
      calcSelf = NAME;
    }
    else
#undef SWITCHING
    {
#define NBPAIR
#undef NBSELF
#undef NBEXCL
#include "ComputeNonbondedHack.h"
      calcPair = NAME;
#undef NBPAIR
#define NBSELF
#undef NBEXCL
#include "ComputeNonbondedHack.h"
      calcSelf = NAME;
    }
  }
}

#undef DECLARATION
#define DEFINITION

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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedUtil.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/11/21 00:00:40 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedUtil.C,v $
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

