/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PME_COULOMB_H__
#define PME_COULOMB_H__

#include "PmeBase.h"
#include "PmeFFT.h"
#include "PmeRealSpace.h"
#include "PmeKSpace.h"

class PmeCoulomb {

public:
 
  PmeCoulomb(PmeGrid grid, int natoms);
  ~PmeCoulomb();

  double compute_recip(PmeParticle particles[], Lattice lattice,
                       double ewaldcof, double virial[], Vector forces[]);

private:
  PmeFFT *myFFT;
  PmeRealSpace *myRealSpace;
  PmeKSpace *myKSpace;

  // The 3-d grid
  double *q_arr; 

  PmeGrid myGrid;   // Grid dimensions
  const int N;            // Number of atoms

  void scale_coordinates(PmeParticle particles[], Lattice lattice);

  void scale_forces(Vector forces[], Lattice lattice);

};

#endif
