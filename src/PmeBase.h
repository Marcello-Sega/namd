
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

struct PmeBox {
  double volume;
  double recipx, recipy, recipz;
  double ewald;
};

struct PmeParticle {
  double x, y, z;
  double cg;
};

// Need this to match DPME's definition for now.
#ifndef DPME_DEFINES_PMEVECTOR
struct PmeVector {
  double x, y, z;
};
#endif

void compute_b_spline(double frac[3], double *M, double *dM, int order);

int find_ewaldcof(double *cutoff, double *dtol, double *ewaldcof);
#endif
