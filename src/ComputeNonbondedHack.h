//-*-c++-*-
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
//   FULLELECT full electrostatics calculation?

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
#undef NOEXCL
#ifdef NBEXCL
#define CLASS ComputeNonbondedExcl
#define EXCL(X) X
#define NOEXCL(X)
#else
#define EXCL(X)
#define NOEXCL(X) X
#endif

#define NAME CLASSNAME(calc)

#undef CLASSNAME
#ifdef NBPAIR
#define CLASSNAME(X) FULLELECTNAME( X ## _pair )
#endif
#ifdef NBSELF
#define CLASSNAME(X) FULLELECTNAME( X ## _self )
#endif
#ifdef NBEXCL
#define CLASSNAME(X) FULLELECTNAME( X ## _excl )
#endif

#undef FULLELECTNAME
#undef FULL
#undef NOFULL
#ifdef FULLELECT
#define FULLELECTNAME(X) SPLITTINGNAME( X ## _fullelect )
#define FULL(X) X
#define NOFULL(X)
#else
#define FULLELECTNAME(X) SPLITTINGNAME( X )
#define FULL(X)
#define NOFULL(X) X
#endif

#undef SPLITTINGNAME

#undef SHIFTING
#define SHIFTING(X)
#ifdef NOSPLIT
#define SPLITTINGNAME(X) LAST( X )
#undef SHIFTING
#define SHIFTING(X) X
#endif

#undef XPLORSPLITTING
#define XPLORSPLITTING(X)
#ifdef SPLIT_XPLOR
#define SPLITTINGNAME(X) LAST( X ## _xplor )
#undef XPLORSPLITTING
#define XPLORSPLITTING(X) X
#endif

#undef C1SPLITTING
#define C1SPLITTING(X)
#ifdef SPLIT_C1
#define SPLITTINGNAME(X) LAST( X ## _c1 )
#undef C1SPLITTING
#define C1SPLITTING(X) X
#endif

#define LAST(X) X

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
#define I_SUB i
#define I_LOWER 0
#define I_UPPER numAtoms - 1
#define J_SUB j
#define J_LOWER i + 1
#define J_UPPER numAtoms
#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedHack.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1004 $	$Date: 1997/04/09 16:15:30 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedHack.h,v $
 * Revision 1.1004  1997/04/09 16:15:30  jim
 * Fixed ABR error in purify.
 *
 * Revision 1.1003  1997/03/14 06:44:55  jim
 * First working versions of full electrostatics splitting functions.
 *
 * Revision 1.1002  1997/02/28 04:47:03  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1001  1997/02/21 20:45:12  jim
 * Eliminated multiple function for switching and modified 1-4 interactions.
 * Now assumes a switching function, but parameters are such that nothing
 * happens, same for modified 1-4.  Slight penalty for rare simulations
 * in which these features are not used, but otherwise no loss and
 * simplifies code.
 *
 * Revision 1.1000  1997/02/06 15:58:09  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:21  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:06  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:35:55  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.6  1996/12/03 21:05:09  jim
 * added support for exclusion correction computes
 *
 * Revision 1.5  1996/11/20 23:16:39  jim
 * first compiling version of generic nonbonded function
 *
 * Revision 1.4  1996/11/08 02:37:12  jim
 * split into two files to hide some preprocessor code from user
 *
 *
 ***************************************************************************/

