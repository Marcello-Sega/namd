
#ifndef PME_REAL_SPACE_H__
#define PME_REAL_SPACE_H__

#include "PmeBase.h"

class PmeRealSpace {
  
public:
  PmeRealSpace(PmeGrid grid, int natoms);
  ~PmeRealSpace();

  void fill_charges(double *q_arr, PmeParticle p[]); 
  void compute_forces(const double *q_arr, const PmeParticle p[], 
                      PmeVector f[]);

private:
  void fill_b_spline(PmeParticle p[]);

  const int N;
  const PmeGrid myGrid;
  double *M, *dM;
};


#endif

