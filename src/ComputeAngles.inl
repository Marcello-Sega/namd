/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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

