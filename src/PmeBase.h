
#ifndef PME_BASE_H__
#define PME_BASE_H__

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct PmeGrid {
  int K1, K2, K3;
  int dim2, dim3;
  int order;
};

struct PmeParticle {
  double x, y, z;
  double cg;
};

void compute_b_spline(double frac[3], double *M, double *dM, int order);

#endif
