
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

  double compute_recip(PmeParticle particles[], PmeBox *box, double virial[],
                       PmeVector forces[]);

private:
  PmeFFT *myFFT;
  PmeRealSpace *myRealSpace;
  PmeKSpace *myKSpace;

  // The 3-d grid
  double *q_arr; 

  PmeGrid myGrid;   // Grid dimensions
  const int N;            // Number of atoms

  void scale_coordinates(PmeParticle particles[], PmeBox *box);

  void scale_forces(PmeVector forces[], PmeBox *box);
};

#endif
