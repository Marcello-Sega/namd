/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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

inline DihedralElem::DihedralElem(const Dihedral *a, const DihedralValue *v)
  {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    atomID[3] = a->atom4;
    value = &v[a->dihedral_type];
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

