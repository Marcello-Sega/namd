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
#include "Molecule.h"
#include "LJTable.h"
#include "ReductionMgr.h"


void ComputeNonbondedUtil::registerReductionData(ReductionMgr *reduction)
{
  reduction->Register(REDUCTION_ELECT_ENERGY);
  reduction->Register(REDUCTION_LJ_ENERGY);
}

void ComputeNonbondedUtil::submitReductionData(BigReal *data, ReductionMgr *reduction, int seq)
{
  reduction->submit(seq, REDUCTION_ELECT_ENERGY, data[electEnergyIndex]);
  reduction->submit(seq, REDUCTION_LJ_ENERGY, data[vdwEnergyIndex]);
}

void ComputeNonbondedUtil::unregisterReductionData(ReductionMgr *reduction)
{
  reduction->unRegister(REDUCTION_ELECT_ENERGY);
  reduction->unRegister(REDUCTION_LJ_ENERGY);
}


void ComputeNonbondedUtil::select(void)
{

  cutoff = Node::Object()->simParameters->cutoff;
  cutoff2 = cutoff*cutoff;
  dielectric_1 = 1/Node::Object()->simParameters->dielectric;
  ljTable = LJTable::Instance();
  mol = Node::Object()->molecule;
  scale14 = Node::Object()->simParameters->scale14;
  switchOn = Node::Object()->simParameters->switchingDist;
  switchOn2 = switchOn*switchOn;
  c0 = 1/(cutoff2-switchOn2);
  c1 = c0*c0*c0;
  c3 = c1 * 4;
  c5 = 1/cutoff2;
  c6 = -4 * c5;

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
#undef NBPAIR
#undef NBSELF
#define NBEXCL
#include "ComputeNonbondedHack.h"
      calcExcl = NAME;
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
#undef NBPAIR
#undef NBSELF
#define NBEXCL
#include "ComputeNonbondedHack.h"
      calcExcl = NAME;
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
#undef NBPAIR
#undef NBSELF
#define NBEXCL
#include "ComputeNonbondedHack.h"
      calcExcl = NAME;
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
#undef NBPAIR
#undef NBSELF
#define NBEXCL
#include "ComputeNonbondedHack.h"
      calcExcl = NAME;
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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedUtil.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:58:13 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedUtil.C,v $
 * Revision 1.1000  1997/02/06 15:58:13  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:25  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:35:59  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.6  1997/01/16 20:00:23  jim
 * Added reduction calls to ComputeNonbondedSelf and ...Pair.
 * Also moved some code from ...Excl to ...Util.
 *
 * Revision 1.5  1996/12/04 17:16:32  jim
 * ComputeNonbondedUtil::select() now caches simulation parameters
 *
 * Revision 1.4  1996/12/03 21:05:09  jim
 * added support for exclusion correction computes
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

