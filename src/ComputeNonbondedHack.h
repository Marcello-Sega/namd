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
#define CLASSNAME(X) M14NAME( X ## _pair )
#endif
#ifdef NBSELF
#define CLASSNAME(X) M14NAME( X ## _self )
#endif
#ifdef NBEXCL
#define CLASSNAME(X) M14NAME( X ## _excl )
#endif

#undef M14FLAG
#undef M14NAME
#undef M14
#undef NOM14
#ifdef MODIFY14
#define M14FLAG 1
#define M14NAME(X) SWNAME( X ## _m14 )
#define M14(X) X
#define NOM14(X)
#else
#define M14FLAG 0
#define M14NAME(X) SWNAME( X )
#define M14(X)
#define NOM14(X) X
#endif

#define LAST(X) X

#undef SWFLAG
#undef SWNAME
#undef SW
#undef NOSW
#ifdef SWITCHING
#define SWFLAG 1
#define SWNAME(X) LAST( X ## _sw )
#define SW(X) X
#define NOSW(X)
#else
#define SWFLAG 0
#define SWNAME(X) LAST( X )
#define SW(X)
#define NOSW(X) X
#endif

#define NBINDEX 2 * ( M14FLAG ) + SWFLAG
#define NBLENGTH 2 * 2

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
#define I_UPPER numAtoms
#define J_SUB j
#define J_LOWER i + 1
#define J_UPPER numAtoms
#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedHack.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:58:09 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedHack.h,v $
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

