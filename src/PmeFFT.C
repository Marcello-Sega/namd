
#include "PmeFFT.h"
#include "common.h"


PmeFFT::PmeFFT(int K1, int K2, int K3) {
#ifdef NAMD_FFTW
  forward_plan = rfftw3d_create_plan(K1, K2, K3, FFTW_REAL_TO_COMPLEX,
		FFTW_MEASURE | FFTW_IN_PLACE);
  backward_plan = rfftw3d_create_plan(K1, K2, K3, FFTW_COMPLEX_TO_REAL,
		FFTW_MEASURE | FFTW_IN_PLACE);
  dim2 = K2;
  dim3 = K3 + (K3 & 1 ? 1 : 2);
  return;
#endif
  NAMD_die("No FFT available for PmeFFT.  Try \"useDPME on\" or recompile.\n");
}

PmeFFT::~PmeFFT() {
#ifdef NAMD_FFTW
  rfftwnd_destroy_plan(forward_plan);
  rfftwnd_destroy_plan(backward_plan);
#endif
}

void PmeFFT::forward(double *q) {
#ifdef NAMD_FFTW
  rfftwnd_one_real_to_complex(forward_plan, q, NULL); 
#endif
}

void PmeFFT::backward(double *q) {
#ifdef NAMD_FFTW
    rfftwnd_one_complex_to_real(backward_plan, (fftw_complex *)q, NULL);
#endif
}


