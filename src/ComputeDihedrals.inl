/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef COMPUTEDIHEDRALS_INL
#define COMPUTEDIHEDRALS_INL

#include "ComputeDihedrals.h"

inline DihedralElem::DihedralElem()
  {
	atomID[0] = -1;
	atomID[1] = -1;
	atomID[2] = -1;
	atomID[3] = -1;
	p[0] = NULL;
	p[1] = NULL;
	p[2] = NULL;
	p[3] = NULL;
  }

inline DihedralElem::DihedralElem(const Dihedral *a)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    atomID[3] = a->atom4;
    dihedralType = a->dihedral_type;
  }

inline DihedralElem::DihedralElem(AtomID atom0, AtomID atom1,
				  AtomID atom2, AtomID atom3)
  {
    if (atom0 > atom3) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom3; atom3 = atom0; atom0 = tmp; 
      tmp = atom1; atom1 = atom2; atom2 = tmp;
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
    atomID[2] = atom2;
    atomID[3] = atom3;
  }

inline int DihedralElem::operator==(const DihedralElem &a) const
  {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1] &&
        a.atomID[2] == atomID[2] && a.atomID[3] == atomID[3]);
  }

inline int DihedralElem::operator<(const DihedralElem &a) const
  {
    return  (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1] ||
            (atomID[1] == a.atomID[1] &&
            (atomID[2] < a.atomID[2] ||
            (atomID[2] == a.atomID[2] &&
             atomID[3] < a.atomID[3] 
	     ))))));
  }

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeDihedrals.inl,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1997/01/14 15:29:47 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDihedrals.inl,v $
 * Revision 1.1  1997/01/14 15:29:47  nealk
 * Initial revision
 *
 * Revision 1.1  1996/11/26 16:33:35  nealk
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

