/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedBase.h support hacks.
 *
 ***************************************************************************/

// Several special cases are defined:
//   NBPAIR, NBSELF, NBEXCL switch environment (mutually exclusive)
//   MODIFY14 modified 1-4 parameters?
//   SWITCHING switching function?

#undef DECL
#undef NODECL
#ifdef DECLARATION
#define DECL(X) X
#define NODECL(X)
#else
#define DECL(X)
#define NODECL(X) X
#endif

#undef CLASS

#undef PAIR
#ifdef NBPAIR
#define CLASS ComputeNonbondedPair
#define PAIR(X) X
#else
#define PAIR(X)
#endif

#undef SELF
#ifdef NBSELF
#define CLASS ComputeNonbondedSelf
#define SELF(X) X
#else
#define SELF(X)
#endif

#undef EXCL
#ifdef NBEXCL
#define CLASS ComputeNonbondedExcl
#define EXCL(X) X
#else
#define EXCL(X)
#endif

#define NAME M14NAME(calc)

#undef M14NAME
#undef M14
#undef NOM14
#ifdef MODIFY14
#define M14NAME(X) SWNAME( X ## _m14 )
#define M14(X) X
#define NOM14(X)
#else
#define M14NAME(X) SWNAME( X )
#define M14(X)
#define NOM14(X) X
#endif

#define LAST(X) X

#undef SWNAME
#undef SW
#undef NOSW
#ifdef SWITCHING
#define SWNAME(X) LAST( X ## _sw )
#define SW(X) X
#define NOSW(X)
#else
#define SWNAME(X) LAST( X )
#define SW(X)
#define NOSW(X) X
#endif

#undef PLEN
#undef I_SUB
#undef I_LOWER
#undef I_UPPER
#undef J_SUB
#undef J_LOWER
#undef J_UPPER
#if defined NBPAIR
#define PLEN [2]
#define I_SUB 0][i
#define I_LOWER 0
#define I_UPPER numAtoms[0]
#define J_SUB 1][j
#define J_LOWER 0
#define J_UPPER numAtoms[1]
#elif defined NBSELF
#define PLEN
#define I_SUB [i]
#define I_LOWER 0
#define I_UPPER numAtoms
#define J_SUB [j]
#define J_LOWER i + 1
#define J_UPPER numAtoms
#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedHack.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1996/11/08 02:37:12 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedHack.h,v $
 * Revision 1.4  1996/11/08 02:37:12  jim
 * split into two files to hide some preprocessor code from user
 *
 *
 ***************************************************************************/

