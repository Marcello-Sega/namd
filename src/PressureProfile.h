/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "common.h"

inline void pp_reduction(BigReal thickness, BigReal min, 
                        int nslabs, BigReal z1, BigReal z2, 
                        int n1, int n2,
                        BigReal vxx, BigReal vyy, BigReal vzz,
                        BigReal *reduction) {

  const BigReal zcell = nslabs*thickness;
  // wrap n and z around
  if (n1 < 0) {
    n1 += nslabs;
    z1 += zcell;
  } else if (n1 >= nslabs) {
    n1 -= nslabs;
    z1 -= zcell;
  }
  if (n2 < 0) {
    n2 += nslabs;
    z2 += zcell;
  } else if (n2 >= nslabs) {
    n2 -= nslabs;
    z2 -= zcell;
  }
  if (n1 == n2) {
    reduction[3*n1] += vxx;
    reduction[3*n1+1] += vyy;
    reduction[3*n1+2] += vzz; 
    return;
  }
  // WOLOG set things up so that n1 > n2.
  if (n1 < n2) {
    // swap
    int tmp = n1; n1 = n2; n2 = tmp;
    BigReal btmp = z1; z1 = z2; z2 = btmp;
  }
  int within_cell = (2*(n1-n2) < nslabs);
  BigReal dz1;  BigReal dz2; BigReal idelta;
  if (within_cell) {
    // dz1 is distance from z1 to bottom of slab; 
    // dz2 is distance from z2 to top of slab.
    idelta = 1.0/(z1-z2);  // must be positive
    dz1 = z1 - (min+n1*thickness);
    dz2 = (min+(n2+1)*thickness) - z2;
  } else {
    // the other side
    idelta = 1.0/(z2-z1+nslabs*thickness);  // must be positive
    dz2 = z2 - (min+n2*thickness);
    dz1 = (min+(n1+1)*thickness) - z1;
  }

// sanity checks
if (dz1*idelta > 1 || dz2*idelta > 1 || idelta < 0) {
CkPrintf("Warning: z1=%f, z2=%f, n1=%d, n2=%d, dz1=%f, dz2=%f, idelta=%f\n",
  (float)z1, (float)z2, n1, n2, (float)dz1, (float)dz2, (float)idelta);
static int first = 1;
if (first) {
  CkPrintf("min=%f, thickness=%f\n", min, thickness);
  first = 0;
}
}

  // Add contributions to the slabs in which the particles are found
  vxx *= idelta;
  vyy *= idelta;
  vzz *= idelta;
  reduction[3*n1] += dz1 * vxx;
  reduction[3*n1+1] += dz1 * vyy;
  reduction[3*n1+2] += dz1 * vzz;
  reduction[3*n2] += dz2 * vxx;
  reduction[3*n2+1] += dz2 * vyy;
  reduction[3*n2+2] += dz2 * vzz;
  if (within_cell) {
    for (int islab = n2+1; islab < n1; islab++) {
      reduction[3*islab] += thickness * vxx;
      reduction[3*islab+1] += thickness * vyy;
      reduction[3*islab+2] += thickness * vzz;
    }
  } else {
    int islab;
    for (islab = 0; islab < n2; islab++) {
      reduction[3*islab] += thickness * vxx;
      reduction[3*islab+1] += thickness * vyy;
            reduction[3*islab+2] += thickness * vzz;
    }
    for (islab = n1+1; islab < nslabs; islab++) {
      reduction[3*islab] += thickness * vxx;
      reduction[3*islab+1] += thickness * vyy;
      reduction[3*islab+2] += thickness * vzz;
    }
  }
}

