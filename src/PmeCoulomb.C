/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "PmeCoulomb.h"
#include "common.h"

PmeCoulomb::PmeCoulomb(PmeGrid grid, int natoms) 
  : myGrid(grid), N(natoms) {
  int qsize;

  myFFT = new PmeFFT(myGrid.K1, myGrid.K2, myGrid.K3);
  myFFT->getdims(&(myGrid.dim2), &(myGrid.dim3));
  myRealSpace = new PmeRealSpace(myGrid, natoms);
  myKSpace = new PmeKSpace(myGrid);

  qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  q_arr = new double[qsize];
  for (int i=0; i<qsize; q_arr[i++] = 0.0);
}

PmeCoulomb::~PmeCoulomb() {
  delete [] q_arr;
}

double PmeCoulomb::compute_recip(PmeParticle p[], Lattice lattice, 
                        double ewaldcof, double virial[], Vector f[]) {
  double energy;
  int i, qsize;

  qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  for (i=0; i<qsize; q_arr[i++] = 0.0);

  scale_coordinates(p, lattice);  // move to nodes
  myRealSpace->fill_charges(q_arr, p);

  myFFT->forward(q_arr);

  energy = myKSpace->compute_energy(q_arr, lattice, ewaldcof, virial);

  myFFT->backward(q_arr);

  myRealSpace->compute_forces(q_arr, p, f);
  scale_forces(f, lattice);  // move to nodes

  return energy;
}

void PmeCoulomb::scale_coordinates(PmeParticle p[], Lattice lattice) { 
  Vector origin = lattice.origin();
  Vector recip1 = lattice.a_r();
  Vector recip2 = lattice.b_r();
  Vector recip3 = lattice.c_r();
  double ox = origin.x;
  double oy = origin.y;
  double oz = origin.z;
  double r1x = recip1.x;
  double r1y = recip1.y;
  double r1z = recip1.z;
  double r2x = recip2.x;
  double r2y = recip2.y;
  double r2z = recip2.z;
  double r3x = recip3.x;
  double r3y = recip3.y;
  double r3z = recip3.z;
  int K1 = myGrid.K1;
  int K2 = myGrid.K2;
  int K3 = myGrid.K3;

  for (int i=0; i<N; i++) {
    double px = p[i].x - ox;
    double py = p[i].y - oy;
    double pz = p[i].z - oz;
    double sx = px*r1x + py*r1y + pz*r1z;
    double sy = px*r2x + py*r2y + pz*r2z;
    double sz = px*r3x + py*r3y + pz*r3z;
    p[i].x = K1 * ( sx - floor(sx) );
    p[i].y = K2 * ( sy - floor(sy) );
    p[i].z = K3 * ( sz - floor(sz) );
  }
}

void PmeCoulomb::scale_forces(Vector f[], Lattice lattice) {
  Vector recip1 = lattice.a_r();
  Vector recip2 = lattice.b_r();
  Vector recip3 = lattice.c_r();
  double r1x = recip1.x;
  double r1y = recip1.y;
  double r1z = recip1.z;
  double r2x = recip2.x;
  double r2y = recip2.y;
  double r2z = recip2.z;
  double r3x = recip3.x;
  double r3y = recip3.y;
  double r3z = recip3.z;
  
  for (int i=0; i<N; i++) {
    double f1 = f[i].x;
    double f2 = f[i].y;
    double f3 = f[i].z;
    f[i].x = f1*r1x + f2*r2x + f3*r3x;
    f[i].y = f1*r1y + f2*r2y + f3*r3y;
    f[i].z = f1*r1z + f2*r2z + f3*r3z;
  }
}

