/***************************************************************************/
/*        (C) Copyright 1996,1997 The Board of Trustees of the             */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *  Inline functions for ComputeAngles.
 *
 ***************************************************************************/

#ifndef COMPUTEANGLE_INL
#define COMPUTEANGLE_INL

#include "ComputeAngles.h"

inline AngleElem::AngleElem()
  {
    atomID[0] = -1;
    atomID[1] = -1;
    atomID[2] = -1;
    p[0] = NULL;
    p[1] = NULL;
    p[2] = NULL;
  }

inline AngleElem::AngleElem(const Angle *a)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    angleType = a->angle_type;
  }

inline AngleElem::AngleElem(AtomID atom0, AtomID atom1, AtomID atom2)
  {
    if (atom0 > atom2) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom2; atom2 = atom0; atom0 = tmp; 
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
    atomID[2] = atom2;
  }

inline int AngleElem::operator==(const AngleElem &a) const
  {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1] &&
        a.atomID[2] == atomID[2]);
  }

inline int AngleElem::operator<(const AngleElem &a) const
  {
    return (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1] ||
            (atomID[1] == a.atomID[1] &&
             atomID[2] < a.atomID[2]) )));
  }

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeAngles.inl,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1997/03/20 23:53:28 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeAngles.inl,v $
 * Revision 1.1003  1997/03/20 23:53:28  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1002  1997/03/16 22:31:19  jim
 * Removed include Debug.h which was screwing up debugging in .C file.
 *
 * Revision 1.1001  1997/02/21 17:38:31  nealk
 * Got DPMTA to work!
 *
 * Revision 1.1000  1997/02/06 15:57:45  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:03  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:35:36  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1997/01/14 15:29:47  nealk
 * Initial revision
 *
 * Revision 1.7  1996/11/22 01:45:28  jim
 * restructured, fixed bugs, now seems to work
 *
 * Revision 1.6  1996/11/19 06:58:37  jim
 * first compiling templated version, needed ugly void* hack
 *
 * Revision 1.5  1996/11/18 21:28:48  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/11/04 20:06:17  nealk
 * Now it compiles :-)
 *
 * Revision 1.3  1996/11/04 19:29:02  nealk
 * Added angleForce() to system, but it is untested.
 *
 * Revision 1.2  1996/11/04 16:55:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/11/01 21:20:45  ari
 * Initial revision
 *
 * Revision 1.3  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.4  1996/07/16 01:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.3  96/07/16  01:10:26  01:10:26  ari (Aritomo Shinozaki)
 * Fixed comments, added methods
 * 
 * Revision 1.2  1996/06/25 21:10:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/24 14:12:26  gursoy
 * Initial revision
 *
 ***************************************************************************/

