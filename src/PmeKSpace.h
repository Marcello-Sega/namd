
#ifndef PME_K_SPACE_H__
#define PME_K_SPACE_H__

#include "PmeBase.h"
#include "Lattice.h"

class PmeKSpace {

public:
  PmeKSpace(PmeGrid grid);
  ~PmeKSpace();

  double compute_energy(double q_arr[], Lattice lattice, double ewald,
                        double virial[]);
  
private:
  // b-spline moduli
  double *bm1, *bm2, *bm3; 
  double *exp1, *exp2, *exp3;
  double i_pi_volume, piob;

  const PmeGrid myGrid;

  void init_exp(double *xp, int K, double recip);
};

#endif
