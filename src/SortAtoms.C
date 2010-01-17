/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "SortAtoms.h"
#include "NamdTypes.h"
#include <algorithm>

// #include "charm++.h"


struct sortop_base {
  const CompAtom * const a;
  sortop_base(const CompAtom* atoms) : a(atoms) { }
};

struct sortop_x : public sortop_base {
  sortop_x(const CompAtom* atoms) : sortop_base(atoms) { }
  bool operator() (int i, int j) const {
    return ( a[i].position.x < a[j].position.x );
  }
};

struct sortop_y : public sortop_base {
  sortop_y(const CompAtom* atoms) : sortop_base(atoms) { }
  bool operator() (int i, int j) const {
    return ( a[i].position.y < a[j].position.y );
  }
};

struct sortop_z : public sortop_base {
  sortop_z(const CompAtom* atoms) : sortop_base(atoms) { }
  bool operator() (int i, int j) const {
    return ( a[i].position.z < a[j].position.z );
  }
};


static void partition(int *order, const CompAtom *atoms, int begin, int end) {

  //  Applies orthogonal recursive bisection with splittings limited
  //  to multiples of 32 for warps.  Partitions containing a multiple
  //  of 16 but not a multiple of 32 are sorted but not split.

  // must have at least 3 atoms for sorting to matter 
  if ( end - begin < 3 ) return;

  // must be a multiple of 16 between begin and end for sorting to matter
  if ( (begin+16)/16 > (end-1)/16 ) return;

  BigReal xmin, ymin, zmin, xmax, ymax, zmax;
  {
    const Position &pos = atoms[order[begin]].position;
    xmin = pos.x;
    ymin = pos.y;
    zmin = pos.z;
    xmax = pos.x;
    ymax = pos.y;
    zmax = pos.z;
  }
  for ( int i=begin+1; i<end; ++i ) {
    const Position &pos = atoms[order[i]].position;
    if ( pos.x < xmin ) { xmin = pos.x; }
    if ( pos.y < ymin ) { ymin = pos.y; }
    if ( pos.z < zmin ) { zmin = pos.z; }
    if ( pos.x > xmax ) { xmax = pos.x; }
    if ( pos.y > ymax ) { ymax = pos.y; }
    if ( pos.z > zmax ) { zmax = pos.z; }
  }
  xmax -= xmin;
  ymax -= ymin;
  zmax -= zmin;
  if ( xmax >= ymax && xmax >= zmax ) {
    std::sort(order+begin, order+end, sortop_x(atoms));
  } else if ( ymax >= xmax && ymax >= zmax ) {
    std::sort(order+begin, order+end, sortop_y(atoms));
  } else {
    std::sort(order+begin, order+end, sortop_z(atoms));
  }

  // must be a multiple of 32 between begin and end to split at
  if ( (begin+32)/32 > (end-1)/32 ) {
    // CkPrintf("%d %d %f %f %f\n", begin, end, xmax, ymax, zmax);
    return;
  }

  // find a multiple of 32 near the median
  int split = ((begin + end + 32) / 64) * 32;

  // recursively partition before and after split
  partition(order, atoms, begin, split);
  partition(order, atoms, split, end);

}

void sortAtomsForCUDA(int *order, const CompAtom *atoms, int nfree, int n) {

  // partition free atoms
  // CkPrintf("%d %d\n", 0, nfree);
  partition(order, atoms, 0, nfree);

  // partition fixed atoms
  // CkPrintf("%d %d\n", nfree, n);
  partition(order, atoms, nfree, n);

}


