/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */


/*************************************************/
/* This program calls functions from pubfft/libpubfft.a which are
   compiled from fortran files, that's why the function pubz3di_ and 
   pubz3d_ may have underscores that may  be taken out. If your 
   platform has a specific 3DFFT routine , eg. KSR, you can use it
*/
#include "dpme2.h"

#ifdef NAMD_FFTW
#include "fftw.h"

fftwnd_plan forward_plan;
fftwnd_plan backward_plan;

#endif

#ifdef NAMD_SGI_COMPLIB_FFT
#include <fft.h>

zomplex *sgi_complib_fft_coeff;

#endif

/* this file has been modified to exclude the SGI & CRAY */
/* optimizatoin stuff, by default we are using the Pubff stuff */
/* this is necessary for the f2c to work. ayt 6/29/94 */
/* ___________________________________________________________ */
/* Subroutine */ int get_fftdims( int *nfft1,  int *nfft2,  int *
	nfft3,  int *nfftdim1,  int *nfftdim2,  int *nfftdim3, 
	 int *nfftable,  int *nffwork,  int *sizfftab,  int *
	sizffwrk)
{
    /* System generated locals */
     int i_1; /* ok */

    /* Local variables */
      int n, nfftmax;

/* Computing MAX */
    i_1 = max(*nfft2,*nfft1);
    nfftmax = max(*nfft3,i_1);
    *nfftdim1 = *nfft1;
#ifndef NAMD_FFTW
    n = *nfft1 / 2;
    if (*nfft1 == n << 1) {
	*nfftdim1 = *nfft1 + 1;
    }
#endif
    *nfftdim2 = *nfft2;
#ifndef NAMD_FFTW
    n = *nfft2 / 2;
    if (*nfft2 == n << 1) {
	*nfftdim2 = *nfft2 + 1;
    }
#endif
    *nfftdim3 = *nfft3;
#ifndef NAMD_FFTW
    n = *nfft3 / 2;
    if (*nfft3 == n << 1) {
	*nfftdim3 = *nfft3 + 1;
    }
#endif
/*     the pubfft is active now */
    *nfftable = (nfftmax << 2) + 15;
    *nffwork = nfftmax;
    *sizfftab = *nfftable * 3;
    *sizffwrk = nfftmax << 1;
    return 0;
} /* get_fftdims */

/* --------------------------------------------------------------- */
/* Subroutine */ int fft_setup(double * DUMMY1 /* array */, double *fftable, 
	double * DUMMY2 /* ffwork */,  int *nfft1,  int *nfft2,  int *nfft3, 
	 int * DUMMY3 /* nfftdim1 */,  int * DUMMY4 /* nfftdim2 */,  int * DUMMY5 /* nfftdim3 */,
	 int *nfftable,  int * DUMMY6 /* nffwork */)
{

    /* Local variables */
    extern /* Subroutine */ int pubz3di( int *,  int *,  int *, 
	    double *,  int *);

#ifdef NAMD_FFTW

  forward_plan = fftw3d_create_plan(*nfft3, *nfft2, *nfft1,
	1, FFTW_MEASURE | FFTW_IN_PLACE);
  backward_plan = fftw3d_create_plan(*nfft3, *nfft2, *nfft1,
	-1, FFTW_MEASURE | FFTW_IN_PLACE);

#else
#ifdef NAMD_SGI_COMPLIB_FFT

  sgi_complib_fft_coeff = zfft3di(*nfft1, *nfft2, *nfft3, NULL);

#else

#if VERBOSE       
    printf("using public domain fft code...\n");
#endif
    pubz3di(nfft1, nfft2, nfft3, fftable, nfftable);

#endif
#endif

    return 0;
} /* fft_setup */

/* ----------------------------------------------------------- */
/* Subroutine */ int fft_forward(doublecomplex *array, double *fftable, 
	doublecomplex *ffwork,  int *nfft1,  int *nfft2,  int *nfft3, 
	 int *nfftdim1,  int *nfftdim2,  int * DUMMY1 /* nfftdim3 */,
	 int *nfftable,  int *nffwork)
{
      int isign;
    extern /* Subroutine */ int pubz3d( int *,  int *,  int *, 
	     int *, doublecomplex *,  int *,  int *, double *, 
	     int *, doublecomplex *,  int *);

#ifdef NAMD_FFTW

  fftwnd_one(forward_plan,(fftw_complex*)array,(fftw_complex*)array);

#else
#ifdef NAMD_SGI_COMPLIB_FFT

  zfft3d(1, *nfft1, *nfft2, *nfft3, (zomplex*)array,
          *nfftdim1, *nfftdim2, sgi_complib_fft_coeff);

#else

    /* Function Body */
    isign = 1;
    pubz3d(&isign, nfft1, nfft2, nfft3, array, nfftdim1, nfftdim2,
	    fftable, nfftable, ffwork, nffwork);

#endif
#endif

    return 0;
} /* fft_forward */

/* ----------------------------------------------------------- */
/* Subroutine */ int fft_back(doublecomplex *array, double *fftable, 
	 doublecomplex *ffwork,  int *nfft1,  int *nfft2,  int *nfft3, 
	 int *nfftdim1,  int *nfftdim2,  int * DUMMY1 /* nfftdim3 */,
	 int *nfftable,  int *nffwork)
{
  int isign;
  extern int pubz3d( int *,  int *,  int *, 
		    int *, doublecomplex *,  int *,  int *, double *, 
		    int *, doublecomplex *,  int *);
  
#ifdef NAMD_FFTW

  fftwnd_one(backward_plan,(fftw_complex*)array,(fftw_complex*)array);

#else
#ifdef NAMD_SGI_COMPLIB_FFT

  zfft3d(-1, *nfft1, *nfft2, *nfft3, (zomplex*)array,
          *nfftdim1, *nfftdim2, sgi_complib_fft_coeff);

#else

  isign = -1;
  pubz3d(&isign, nfft1, nfft2, nfft3, array, nfftdim1, nfftdim2,
	 fftable, nfftable, ffwork, nffwork);

#endif
#endif

  return 0;
} /* fft_back */

