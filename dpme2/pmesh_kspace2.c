/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */

#include "dpme2.h"

/* Subroutine */ 
int pmesh_kspace_get_sizes( int *nfft1,  int *nfft2, 
     int *nfft3,  int *numatoms,  int *order,  int *sizfftab, 
     int *sizffwrk,  int *siztheta,  int *siz_q,  int *
    sizheap,  int *sizstack)
{

    /* Local variables */
    int nfftable;
  extern /* Subroutine */ int get_fftdims( int *,  int *,  int *,
	 int *,  int *,  int *,  int *,  int *,  int *,
	      int *);
    int nffwork, nfftdim1, nfftdim2, nfftdim3;
   
/* INPUT */
/*      nfft1,nfft2,nfft3,numatoms,order */
/*      nfft1,nfft2,nfft3 are the dimensions of the charge grid array */
/*      numatoms is number of atoms */
/*      order is the order of B-spline interpolation */
/* OUTPUT */
/*      sizfftab,sizffwrk,siztheta,siz_Q */
/*      sizfftab is permanent 3d fft table storage */
/*      sizffwrk is temporary 3d fft work storage */
/*      siztheta is size of arrays theta1-3 dtheta1-3 */
/*      sizheap is total size of permanent storage */
/*      sizstack is total size of temporary storage */
/* This routine computes the above output parameters needed for */
/* heap or stack allocation. */
    get_fftdims(nfft1, nfft2, nfft3, &nfftdim1, &nfftdim2, &nfftdim3, &
	    nfftable, &nffwork, sizfftab, sizffwrk);
    *siztheta = *numatoms * *order;
    *siz_q = (nfftdim1 << 1) * nfftdim2 * nfftdim3;
    *sizheap = *nfft1 + *nfft2 + *nfft3 + *sizfftab;
    *sizstack = *siz_q + *siztheta * 6 + *sizffwrk + *numatoms * 3;
    printf("total HEAP storage needed = %d \n",(*sizheap));
    printf("total STACK storage needed = %d \n",(*sizstack));
    return 0;
} /* pmesh_kspace_get_sizes */

/* ---------------------------------------------------- */
int pmesh_kspace_setup(double *bsp_mod1, double *
		       bsp_mod2, double *bsp_mod3, double *fftable, double *ffwork,  
		       int *nfft1,  int *nfft2,  int *nfft3,  int *order,  
		       int *sizfftab,  int *sizffwrk)
{
  int nfftable, sfft, sffw;
  double *dummy;
  extern int fft_setup(double *, double *, 
		       double *,  int *,  int *,  int *,  int *,  int 
		       *,  int *,  int *,  int *), 
  get_fftdims( int *, 
	      int *,  int *,  int *,  int *,  int *,  int *, 
	      int *,  int *,  int *);
  int nffwork;
  extern int load_bsp_moduli(double *, double *, 
			     double *,  int *,  int *,  int *,  int *);
  int nfftdim1, nfftdim2, nfftdim3;
  
    /* Parameter adjustments */
      --bsp_mod1;
      --bsp_mod2;
      --bsp_mod3;
      --fftable;
      /* --ffwork; */


    /* Function Body */
/*  see DO_PMESH_KSPACE for explanation of arguments */
  get_fftdims(nfft1, nfft2, nfft3, &nfftdim1, &nfftdim2, &nfftdim3, &
	      nfftable, &nffwork, &sfft, &sffw);

  load_bsp_moduli(&bsp_mod1[1], &bsp_mod2[1], &bsp_mod3[1], nfft1, 
		  nfft2, nfft3, order);

  fft_setup(dummy, &fftable[1],ffwork, nfft1, nfft2, nfft3,
	    &nfftdim1, &nfftdim2, &nfftdim3, &nfftable, &nffwork);

  return 0;
} /* pmesh_kspace_setup */

/* ---------------------------------------------------- */
#if NOTDEFINED
this is not working with DPME2
int do_pmesh_kspace2( int *numatoms,
		    PmeParticlePtr ParticlePtr , double *recip, 
		    double *volume, double *ewald_coeff,  int *order, 
		    int *nfft1,  int *nfft2,  int *nfft3, double *eer, 
		    PmeVectorPtr rfparticle, double *virial, int *sizfftab,  
		    int *sizffwrk,  int *siztheta,  int *siz_q, double *bsp_mod1, 
		    double *bsp_mod2, double *bsp_mod3, double *fftable, double *q, 
		    double *ffwork, double *theta1, double *theta2, double *theta3, 
		    double *dtheta1, double *dtheta2, double *dtheta3, 
		    double *fr1, double *fr2, double *fr3, double *time)
{
  
  int nfftable, sfft, sffw;
  int nfftdim1, nfftdim2, nfftdim3, nffwork;
  struct timeval time1, time2;
  struct timezone tzp;
  double walltime; /* ayt */
  int i;
 /* Parameter adjustments */
/*********************************/
    recip -= 4;
    --virial;
    --time;
    --q;
    --ffwork;
      
    --bsp_mod1;
    --bsp_mod2;
    --bsp_mod3;
    --fftable;
   
    --theta1;
    --theta2;
    --theta3;
    --dtheta1;
    --dtheta2;
    --dtheta3;
    --fr1;
    --fr2;
    --fr3;
  

    /* Function Body */
/* INPUT */
/*       numatoms:  number of atoms */
/*       x,y,z:   atomic coords */
/*       charge  atomic charges */
/*      recip: 3x3 array of reciprocal unit cell vectors (stored as columns)*/
/*       volume: the volume of the unit cell */
/*       ewald_coeff:   ewald convergence parameter */
/*       order: the order of Bspline interpolation. E.g. cubic is order 4 */
/*          fifth degree is order 6 etc. The order must be an even number */
/*          and at least 4. */
/*       nfft1,nfft2,nfft3: the dimensions of the charge grid array */
/* OUTPUT */
/*       eer:  ewald reciprocal or k-space  energy */
/*       rfparticle(x,y,z): forces incremented by k-space sum */
/*       virial:  virial due to k-space sum (valid for atomic scaling; */
/*                rigid molecule virial needs a correction term not */
/*                computed here */
/*       time: used to profile the different component routines */
/* SIZES of some arrays */
/* HEAP STORAGE:  These arrays need to be preserved throughout simulation */
/* STACK STORAGE: These arrays can be tossed after leaving this routine */
/*******************************************************************************/

/*  get some  int array dimensions */
    get_fftdims(nfft1, nfft2, nfft3, &nfftdim1, &nfftdim2, &nfftdim3, &
	    nfftable, &nffwork, &sfft, &sffw);
  /* ayt, 4/96 this subroutine is called elsewhere    
     get_scaled_fractionals(numatoms, ParticlePtr, &recip[4], nfft1, 
	nfft2, nfft3, &fr1[1], &fr2[1], &fr3[1]);
   */
#if TIMEME
    gettimeofday(&(time1),&tzp);  /* eqv to second */
#endif

    get_bspline_coeffs(numatoms, &fr1[1], &fr2[1], &fr3[1], order, &theta1[
	    1], &theta2[1], &theta3[1], &dtheta1[1], &dtheta2[1], &dtheta3[1]);

#if TIMEME
    gettimeofday(&(time2),&tzp);  /* eqv to second */
    time[1] +=  ( (double)time2.tv_sec
                        - (double)time1.tv_sec 
                        + ((double)time2.tv_usec
                        - (double)time1.tv_usec) / 1000000) ;
    time1 = time2;
#endif
    
    fill_charge_grid(numatoms, ParticlePtr, &theta1[1], &theta2[1], &theta3[
	    1], &fr1[1], &fr2[1], &fr3[1], order, nfft1, nfft2, nfft3, &

	    nfftdim1, &nfftdim2, &nfftdim3, &q[1]);

#if TIMEME
    gettimeofday(&(time2),&tzp);  /* eqv to second */
    time[2] +=  ( (double)time2.tv_sec
                        - (double)time1.tv_sec 
                        + ((double)time2.tv_usec
                        - (double)time1.tv_usec) / 1000000) ;
    time1 = time2;
#endif

    fft_back(&q[1], &fftable[1], &ffwork[1], nfft1, nfft2, nfft3, &nfftdim1,
	     &nfftdim2, &nfftdim3, &nfftable, &nffwork);
#if TIMEME
    gettimeofday(&(time2),&tzp);  /* eqv to second */
    time[3] +=  ( (double)time2.tv_sec
                        - (double)time1.tv_sec 
                        + ((double)time2.tv_usec
                        - (double)time1.tv_usec) / 1000000) ;
     time1 = time2;
#endif

    scalar_sum(&q[1], ewald_coeff, volume, &recip[4], &bsp_mod1[1], &
	    bsp_mod2[1], &bsp_mod3[1], nfft1, nfft2, nfft3, &nfftdim1, &
	    nfftdim2, &nfftdim3, eer, &virial[1]);

#if TIMEME
    gettimeofday(&(time2),&tzp);  /* eqv to second */
    time[4] +=  ( (double)time2.tv_sec
                        - (double)time1.tv_sec 
                        + ((double)time2.tv_usec
                        - (double)time1.tv_usec) / 1000000) ;
    time1 = time2;
#endif

    fft_forward(&q[1], &fftable[1], &ffwork[1], nfft1, nfft2, nfft3, &
	    nfftdim1, &nfftdim2, &nfftdim3, &nfftable, &nffwork);

#if TIMEME
    gettimeofday(&(time2),&tzp);  /* eqv to second */
    time[3] +=  ( (double)time2.tv_sec
                        - (double)time1.tv_sec 
                        + ((double)time2.tv_usec
                        - (double)time1.tv_usec) / 1000000) ;
    time1 = time2;
#endif

  grad_sum(numatoms, ParticlePtr, &recip[4], &theta1[1], &theta2[1], &
	    theta3[1], &dtheta1[1], &dtheta2[1], &dtheta3[1], rfparticle,
	    &fr1[1], &fr2[1], &fr3[1], order, nfft1, nfft2, nfft3, &
	    nfftdim1, &nfftdim2, &nfftdim3, &q[1]);
 
#if TIMEME
  gettimeofday(&(time2),&tzp);  /* eqv to second */
  time[5] +=  ( (double)time2.tv_sec
                        - (double)time1.tv_sec 
                        + ((double)time2.tv_usec
                        - (double)time1.tv_usec) / 1000000) ;
#endif
    /* printf(" from the  do_pmesh_kspace() the value of eer = %lf \n",*eer);*/
    return 0;

} /* do_pmesh_kspace */
#endif /* undefined */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ayt 4/96 this routine is different from its first version in that it :
   1- does Not require adjusting the arrays 
   2- accepts Pme2Particle data structure 
   *3- it adds a 0.5 to w to account for the scaling of paritcle coordinates
       that is accociated wiht LJS method. Note that you need to add
       different corrections if the box is not a cube/rectangle and recip[] of-diagnols
       are not 0.
*/
/* INPUT: */
/*      numatoms: number of atoms */
/*      x,y,z: arrays of cartesian coords */
/*      recip: the 3x3 array of reciprocal vectors stored as columns */
/* OUTPUT: */
/*     fr1,fr2,fr3 the scaled and shifted fractional coords */

/* ayt 4/96 added 1/2 (ie prd2/box) to offset LJS method  scaling (note above)*/
/* ayt 12/96 added atompnt and list to properly access the nlocal atoms */
/* ayt 3/97 deleted atompnt,list,tmpcg for the master worker version */

int get_scaled_fractionals2( int *numatoms,
			    Pme2Particle *ParticlePtr, double *recip,  int *nfft1, 
			    int *nfft2,  int *nfft3, double *fr1, double *fr2, 
			    double *fr3)
{
    
  /* Local variables */
  int n,ii;
  double w;
    
#if 0
  printf("get_scaled_frac... ORIG \n");
  for (n = 0; n < (*numatoms); ++n) {
    w = ParticlePtr[n].x * recip[0] +  ParticlePtr[n].y * recip[1] +  ParticlePtr[n].z * recip[2] + 0.5;
    fr1[n] = *nfft1 * (w - Nint(w) + .5);
    w = ParticlePtr[n].x * recip[3] +  ParticlePtr[n].y * recip[4] +  ParticlePtr[n].z * recip[5] + 0.5;
    fr2[n] = *nfft2 * (w - Nint(w) + .5);
    w = ParticlePtr[n].x * recip[6] + ParticlePtr[n].y * recip[7] + ParticlePtr[n].z * recip[8] + 0.5;
    fr3[n] = *nfft3 * (w - Nint(w) + .5);
  }
#endif
  
#if 1

#if DPME_DEBUG 
  printf("get_scaled_frac... MODFD  \n");
#endif

  --ParticlePtr;    /* if u pass Myparticle[1] in resip_sum_setup() */

  n= 1; /* ayt 3/97 to comply with master worker model */
  
  for (ii= 0; ii < (*numatoms); ++ii,++n) {
  
    w = ParticlePtr[n].x * recip[0] +  
      ParticlePtr[n].y * recip[1] +  ParticlePtr[n].z * recip[2] + 0.5;
   
    /* used to be fr123[n], but this will create a prob. in packing later */
    fr1[ii] = *nfft1 * (w - Nint(w) + .5);
    
    w = ParticlePtr[n].x * recip[3] + 
      ParticlePtr[n].y * recip[4] +  ParticlePtr[n].z * recip[5] + 0.5;
    
    fr2[ii] = *nfft2 * (w - Nint(w) + .5);
    
    w = ParticlePtr[n].x * recip[6] + 
      ParticlePtr[n].y * recip[7] + ParticlePtr[n].z * recip[8] + 0.5;
    
    fr3[ii] = *nfft3 * (w - Nint(w) + .5);

        
  }
#endif

  return 0;
} /* get_scaled_fractionals2 */

/* --------------------------------------------------------------- */
/* this routine loads the moduli of the inverse DFT of the B splines */
/* bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions, */
/* Order is the order of the B spline approx. */

int load_bsp_moduli(double *bsp_mod1, double *
		    bsp_mod2, double *bsp_mod3,  int *nfft1,  int *nfft2, 
		    int *nfft3,  int *order)
{
  /* System generated locals */
  int i_1; /* ok */
  
  /* Local variables */
  int nmax;
  int i; 
  double w;
  
  double *array, *darray,*bsp_arr;
  
  
  extern /* Subroutine */ int dftmod(double *, double *,  int *)  ;
  extern /* Subroutine */ int fill_bspline(double *,  int *, 
					   double *, double *); 
  /* Parameter adjustments */
  --bsp_mod1;
  --bsp_mod2;
  --bsp_mod3;

  i_1 = max(*nfft2,*nfft1);
  nmax = max(*nfft3,i_1);

  /* alloc memory */
  array= dvector_pme(*order +1 );
  darray= dvector_pme(*order +1 );
  bsp_arr= dvector_pme(nmax +1);

  w = 0.;
  fill_bspline(&w, order, array, darray);
  
  for (i = 1; i <=nmax; ++i) {
    bsp_arr[i - 1] = 0.;
  }
  i_1 = *order + 1;
  for (i = 2; i <= i_1; ++i) {
    bsp_arr[i - 1] = array[i - 2];
  }
  dftmod(&bsp_mod1[1], bsp_arr, nfft1);
  dftmod(&bsp_mod2[1], bsp_arr, nfft2);
  dftmod(&bsp_mod3[1], bsp_arr, nfft3);

  free(array),free(darray),free(bsp_arr);
  return 0;
} /* load_bsp_moduli */

/*------------------------------------------------------------------------*/
/* Subroutine */ int dftmod(double *bsp_mod, double *bsp_arr, 
	 int *nfft)
{
    
    /* Builtin functions */
    double cos(double), sin(double);

    /* Local variables */
     double tiny;
      int j, k;
     double twopi, arg, sum1, sum2;

    /* Parameter adjustments */
    --bsp_mod;
    --bsp_arr;

    /* Function Body */
/* Computes the modulus of the discrete fourier transform of bsp_arr, */
/*  storing it into bsp_mod */
    twopi =  2.0 * 3.14159265358979323846;
    tiny = 1.e-7;
    
    for (k = 1; k <=(*nfft); ++k) {
	sum1 = 0.;
	sum2 = 0.;
	
	for (j = 1; j <= (*nfft); ++j) {
	    arg = twopi * (k - 1) * (j - 1) / *nfft;
	    sum1 += bsp_arr[j] * cos(arg);
	    sum2 += bsp_arr[j] * sin(arg);
/* L250: */
	}
	bsp_mod[k] = sum1*sum1 + sum2*sum2;
/* L300: */
    }
   
    for (k = 1; k <=(*nfft); ++k) {
	if (bsp_mod[k] < tiny) {
	    bsp_mod[k] = (bsp_mod[k - 1] + bsp_mod[k + 1]) * .5;
	}
/* L400: */
    }
    return 0;
} /* dftmod_ */
/******************************************************************/
/* alloc memory */
double *dvector_pme(int nh)
{
  double *v;
  v=(double *)malloc(((nh)*sizeof(double)));
  if (!v) { 
    fprintf(stderr,"ERROR in allocating a pme vector,exiting!!!!\n");
    exit(0);
  }
  return v;
}
/******************************************************************/
