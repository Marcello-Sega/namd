
#include "PmeCoulomb.h"

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

double PmeCoulomb::compute_recip(PmeParticle p[], PmeBox *box, 
                        double virial[], PmeVector f[]) {
  double energy;
  int i, qsize;

  qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
  for (i=0; i<qsize; q_arr[i++] = 0.0);

  scale_coordinates(p, box);
  myRealSpace->fill_charges(q_arr, p);

  myFFT->forward(q_arr);

  energy = myKSpace->compute_energy(q_arr, box, virial);

  myFFT->backward(q_arr);

  myRealSpace->compute_forces(q_arr, p, f);
  scale_forces(f, box);

  return energy;
}

void PmeCoulomb::scale_coordinates(PmeParticle p[], PmeBox *box) { 
  int i;
  double recipx, recipy, recipz, fr1, fr2, fr3; 
  int K1, K2, K3;
  recipx=box->recipx;
  recipy=box->recipy;
  recipz=box->recipz;
  K1 = myGrid.K1;
  K2 = myGrid.K2;
  K3 = myGrid.K3;

  for (i=0; i<N; i++) {
    fr1=p[i].x*recipx + 0.5;
    fr2=p[i].y*recipy + 0.5;
    fr3=p[i].z*recipz + 0.5;
    p[i].x = K1*(fr1-(int)(fr1+0.5) + 0.5);
    p[i].y = K2*(fr2-(int)(fr2+0.5) + 0.5);
    p[i].z = K3*(fr3-(int)(fr3+0.5) + 0.5);
  }
}

void PmeCoulomb::scale_forces(PmeVector f[], PmeBox *box) {
  int i;
  double rxx, ryy, rzz;
  double f1, f2, f3;
  rxx=box->recipx;
  ryy=box->recipy;
  rzz=box->recipz;
  
  for (i=0; i<N; i++) {
    f1=rxx*f[i].x;
    f2=             ryy*f[i].y;
    f3=                           rzz*f[i].z;
    f[i].x = f1;
    f[i].y = f2;
    f[i].z = f3;
  }
}

