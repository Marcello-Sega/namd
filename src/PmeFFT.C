/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "PmeFFT.h"
#include "common.h"
#include "InfoStream.h"
#include "pub3dfft.h"

PmeFFT::PmeFFT(int K1, int K2, int K3) {
  k1 = K1; k2 = K2; k3 = K3;
  dim2 = K2;
  dim3 = 2 * (K3/2 + 1);
#if defined(NAMD_FFTW)
  iout << iINFO << "PME using FFTW (http://www.fftw.org/).\n" << endi;
  forward_plan = rfftw3d_create_plan(K1, K2, K3, FFTW_REAL_TO_COMPLEX,
		FFTW_MEASURE | FFTW_IN_PLACE);
  backward_plan = rfftw3d_create_plan(K1, K2, K3, FFTW_COMPLEX_TO_REAL,
		FFTW_MEASURE | FFTW_IN_PLACE);
#elif defined(NAMD_SGI_COMPLIB_FFT)
  iout << iINFO << "PME using FFT from SGI complib.\n" << endi;
  sgi_complib_fft_coeff = dzfft3dui(K3, K2, K1, NULL);
#else
  iout << iINFO << "PME using PUB3DFFT, for better performance recompile with FFTW.\n" << endi;
  if ( K3 & 1 ) {
    NAMD_die("Third dimension of FFT grid must be even for PUB3DFFT.");
  }
  int maxdim = K1 > K2 ? K1 : K2; if ( (K3/2+1) > maxdim ) maxdim = K3/2+1;
  ntable = 4*maxdim + 15;
  table = new double[3*ntable];
  work = new double[2*maxdim];
  pubd3di(K3,K2,K1,table,ntable);
#endif
}

PmeFFT::~PmeFFT() {
#if defined(NAMD_FFTW)
  rfftwnd_destroy_plan(forward_plan);
  rfftwnd_destroy_plan(backward_plan);
#elif defined(NAMD_SGI_COMPLIB_FFT)
  // this should free sgi_complib_fft_coeff somehow
#else
  delete [] table;
  delete [] work;
#endif
}

void PmeFFT::forward(double *q) {
#if defined(NAMD_FFTW)
  rfftwnd_one_real_to_complex(forward_plan, q, NULL); 
#elif defined(NAMD_SGI_COMPLIB_FFT)
  dzfft3du(-1 ,k3, k2, k1, q, dim3, dim2, sgi_complib_fft_coeff);
#else
  pubdz3d(-1,k3,k2,k1,q,dim3,dim2,table,ntable,work);
#endif
}

void PmeFFT::backward(double *q) {
#if defined(NAMD_FFTW)
    rfftwnd_one_complex_to_real(backward_plan, (fftw_complex *)q, NULL);
#elif defined(NAMD_SGI_COMPLIB_FFT)
  zdfft3du(1 ,k3, k2, k1, q, dim3, dim2, sgi_complib_fft_coeff);
#else
  pubzd3d(1,k3,k2,k1,q,dim3,dim2,table,ntable,work);
#endif
}


