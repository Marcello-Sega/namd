/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#ifndef NAMD_RESTRICT
#define restrict
#endif

#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "Node.h"
#include "Molecule.h"
#include "LJTable.h"
#include "ReductionMgr.h"

#include "PressureProfile.h"

// clear all
// define interaction type (pair or self)
#define NBPAIR	1
#define NBSELF	2

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define FEPFLAG
#define CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#undef CALCENERGY
#undef FEPFLAG

#define LESFLAG

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#define CALCENERGY
#include "ComputeNonbondedBase.h"
#undef CALCENERGY
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#undef LESFLAG

#define INTFLAG
#define CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#undef  NBTYPE

#undef CALCENERGY
#undef INTFLAG

