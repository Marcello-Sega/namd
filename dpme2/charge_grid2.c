/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */

#include "dpme2.h"

int fill_charge_grid2(int *numatoms, int nlocal,
		      Pme2Particle *ParticlePtr,
		      double *theta1, double *theta2, double *theta3, 
	double *fr1, double *fr2, double *fr3,  int *order, 
	int *nfft1,  int *nfft2,  int *nfft3,  int *nfftdim1, 
	int *nfftdim2, int *nfftdim3, double *q)
{
/* --------------------------------------------------------------------- */
/* INPUT: */
/*      numatoms:  number of atoms */
/*      nlocal: the number of atoms this processor maintains */
/*      charge: the array of atomic charges */
/*      theta1,theta2,theta3: the spline coeff arrays */
/*      fr1,fr2,fr3 the scaled and shifted fractional coords */
/*      nfft1,nfft2,nfft3: the charge grid dimensions */
/*      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims */
/*      order: the order of spline interpolation */
/* OUTPUT: */
/*      Q the charge grid */
/* --------------------------------------------------------------------- */
  int theta1_dim1, theta1_offset, theta2_dim1, theta2_offset, 
      theta3_dim1, theta3_offset, q_dim2, q_dim3, q_offset;

    /* Local variables */
    double prod;
    int ntot, i, j, k, n, i0, j0, k0;
    extern /* Subroutine */ int clearq(double *, int *);
    int ith1, ith2, ith3;
    int nn, ncg ;

    /* Parameter adjustments */
    theta1_dim1 = *order;
    theta1_offset = theta1_dim1 + 1;
    theta1 -= theta1_offset;
    theta2_dim1 = *order;
    theta2_offset = theta2_dim1 + 1;
    theta2 -= theta2_offset;
    theta3_dim1 = *order;
    theta3_offset = theta3_dim1 + 1;
    theta3 -= theta3_offset;

    /* ayt 4/96 changes --fr1;    --fr2;    --fr3; */

    q_dim2 = *nfftdim1;
    q_dim3 = *nfftdim2;
    q_offset = (q_dim2 * (q_dim3 + 1) + 1 << 1) + 1;
    q -= q_offset;

    ntot = (*nfftdim1 * 2) * *nfftdim2 * *nfftdim3;
   
    clearq(&q[q_offset], &ntot);


#if 0
    /* PART 1 , loop over my local particles */
    nn=0; /* use in indexing arrays created w/dvecotor ie fr123 */
   
    for (n = 1; n <= (nlocal); n++,++nn) {
      k0 = (int)(fr3[nn]) - *order;
      for (ith3 = 1; ith3 <= (*order); ++ith3) {
	++k0;
	k = k0 + 1 + (*nfft3 - Nsign(*nfft3, k0)) / 2;
	j0 = (int) fr2[nn] - *order;
	for (ith2 = 1; ith2 <= (*order); ++ith2) {
	  ++j0;
	  j = j0 + 1 + (*nfft2 - Nsign(*nfft2, j0)) / 2;
	  prod = theta2[ith2 + n * theta2_dim1] *
	    theta3[ith3 + n * theta3_dim1] * ParticlePtr[n].cg;
	  i0 = (int)fr1[nn] - *order;
	  for (ith1 = 1; ith1 <=(*order); ++ith1) {
	    ++i0;
	    i = i0 + 1 + (*nfft1 - Nsign(*nfft1, i0)) / 2;
	    /* the *2 below replaced a shift oprtr "<<1"  ayt */
	    q[(i + (j + k * q_dim3) * q_dim2 << 1) + 1] += 
	      theta1[ith1 + n * theta1_dim1] * prod;
	    
	  }
	}
      }
    }
    
    /* part 2 loop over the rest of atoms. deleted ayt 3/97 */
  
#endif

#if 1

#if DPME_DEBUG
    printf("fill_charge_grid modified ...\n");
#endif
   /* PART 1 , loop over my local particles */
    nn=1; /* ayt 3/97. use in indexing arrays created w/dvecotor ie fr123 */
   
    for (n = 1; n <= (nlocal); n++,nn++) {
      k0 = (int)(fr3[n-1]) - *order;
      for (ith3 = 1; ith3 <= (*order); ++ith3) {
	++k0;
	k = k0 + 1 + (*nfft3 - Nsign(*nfft3, k0)) / 2;
	j0 = (int) fr2[n-1] - *order;
	for (ith2 = 1; ith2 <= (*order); ++ith2) {
	  ++j0;
	  j = j0 + 1 + (*nfft2 - Nsign(*nfft2, j0)) / 2;
	  prod = theta2[ith2 + n * theta2_dim1] *
	    theta3[ith3 + n * theta3_dim1] * ParticlePtr[nn].cg;
	  i0 = (int)fr1[n-1] - *order;
	  for (ith1 = 1; ith1 <=(*order); ++ith1) {
	    ++i0;
	    i = i0 + 1 + (*nfft1 - Nsign(*nfft1, i0)) / 2;
	    /* the *2 below replaced a shift oprtr "<<1"  ayt */
	    q[(i + (j + k * q_dim3) * q_dim2 << 1) + 1] += 
	      theta1[ith1 + n * theta1_dim1] * prod;
	    
	  }
	}
      }
    }
    
    /* part 2 loop over the rest of atoms, is deleted ayt 3/97  */
   
#endif    
    return 0;
  } /* fill_charge_grid2 */

/* ----------------------------------------------------------- */
/* Subroutine */ int clearq(double *q, int *ntot)
{
  /* Local variables */
    int i;

    /* Parameter adjustments */
    --q;

    /* Function Body */
    for (i = 1; i <= (*ntot); ++i) {
	q[i] = 0.;
/* L10: */
    }
    return 0;
} /* clearq_ */

/* ----------------------------------------------------------- */
/* Subroutine */ int scalar_sum(double *q, double *ewaldcof, 
	double *volume, double *recip, double *bsp_mod1, 
	double *bsp_mod2, double *bsp_mod3, int *nfft1, 
	int *nfft2, int *nfft3, int *nfftdim1, int *nfftdim2, 
	int *nfftdim3, double *eer, double *vir)
{
    /* System generated locals */
     int q_dim2, q_dim3, q_offset;
    double d_1, d_2;  /* ok */

    /* Builtin functions */
    double exp(double);

    /* Local variables */
    double mhat1, mhat2, mhat3;
    int k;
    double denom, eterm;
    int k1;
    double vterm;
    int k2, k3, m1, m2, m3;
    double struc2, pi, energy;
    int indtop, nf1, nf2, nf3;
    double fac;
    int nff, ind, jnd;
    double msq;

    /* Parameter adjustments */
    q_dim2 = *nfftdim1;
    q_dim3 = *nfftdim2;
    q_offset = (q_dim2 * (q_dim3 + 1) + 1 << 1) + 1;
    q -= q_offset;
    /* recip -= 4; */ /* ayt 5/96 */
    --bsp_mod1;
    --bsp_mod2;
    --bsp_mod3;
    /* --vir;  ayt 5/96 */

    /* Function Body */
    indtop = *nfft1 * *nfft2 * *nfft3;
    pi = 3.14159265358979323846;
    fac = pi*pi / (*ewaldcof * *ewaldcof);
    nff = *nfft1 * *nfft2;
    nf1 = *nfft1 / 2;
    if (nf1 << 1 < *nfft1) {
	++nf1;
    }
    nf2 = *nfft2 / 2;
    if (nf2 << 1 < *nfft2) {
	++nf2;
    }
    nf3 = *nfft3 / 2;
    if (nf3 << 1 < *nfft3) {
	++nf3;
    }
    energy = 0.;
#if VIRIAL
    for (k = 1; k <= 6; ++k) 
      vir[k] = 0.;
#endif

    for (ind = 1; ind <= (indtop - 1); ++ind) {
      /* get k1,k2,k3 from the relationship:*/
      /* ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1 */
	k3 = ind / nff + 1;
	jnd = ind - (k3 - 1) * nff;
	k2 = jnd / *nfft1 + 1;
	k1 = jnd - (k2 - 1) * *nfft1 + 1;
	
	m1 = k1 - 1; 
	m2 = k2 - 1; 
	m3 = k3 - 1;

	if (k1 > nf1) 
	  m1 -=  *nfft1;
	if (k2 > nf2)
	  m2 -= *nfft2;
	if (k3 > nf3) 
	  m3 -= *nfft3;
	
	mhat1 = recip[0] * m1 + recip[3] * m2 + recip[6] * m3;
	mhat2 = recip[1] * m1 + recip[4] * m2 + recip[7] * m3;
	mhat3 = recip[2] * m1 + recip[5] * m2 + recip[8] * m3;
	msq = mhat1 * mhat1 + mhat2 * mhat2 + mhat3 * mhat3;
	denom = pi * *volume * bsp_mod1[k1] * bsp_mod2[k2] * bsp_mod3[
		k3] * msq;
	eterm = exp(-fac * msq) / denom;
	
/* Computing 2nd power */
         d_1 = q[(k1 + (k2 + k3 * q_dim3) * q_dim2 << 1) + 1]; 
        /* printf(" $$$$$$$$$at charge_grid.c d1=%f \n",d_1);
	 * d_1 = q[(k1 + (k2 + k3 * q_dim3) * q_dim2 * 2) + 1];
	 * printf("$$$$$d1=%f \n",d_1); 
	 */
/* Computing 2nd power */
        d_2 = q[(k1 + (k2 + k3 * q_dim3) * q_dim2 << 1) + 2];
	/* d_2 = q[(k1 + (k2 + k3 * q_dim3) * q_dim2 * 2) + 2]; ayt */
	struc2 = d_1 * d_1 + d_2 * d_2;
	energy += eterm * struc2;
        /* printf("eterm %f %f %f %f \n",eterm,struc2,eterm*struc2,energy);  */
	/* printf(" %f  %f \n",d_1,struc2); */
#if VIRIAL
	vterm = (fac * msq + 1.) * 2. / msq;
	vir[1] += eterm * struc2 * (vterm * mhat1 * mhat1 - 1.);
	vir[2] += eterm * struc2 * (vterm * mhat1 * mhat2);
	vir[3] += eterm * struc2 * (vterm * mhat1 * mhat3);
	vir[4] += eterm * struc2 * (vterm * mhat2 * mhat2 - 1.);
	vir[5] += eterm * struc2 * (vterm * mhat2 * mhat3);
	vir[6] += eterm * struc2 * (vterm * mhat3 * mhat3 - 1.);
#endif
	q[(k1 + (k2 + k3 * q_dim3) * q_dim2 << 1) + 1] = eterm * q[(k1 + (k2 
		+ k3 * q_dim3) * q_dim2 << 1) + 1];
	q[(k1 + (k2 + k3 * q_dim3) * q_dim2 << 1) + 2] = eterm * q[(k1 + (k2 
		+ k3 * q_dim3) * q_dim2 << 1) + 2];
	
      } /* end for */
 
    *eer = energy * .5;
   
#if VIRIAL
    for (k = 1; k <= 6; ++k) 
	vir[k] *= .5;
#endif

    return 0;
} /* scalar_sum */

/* ----------------------------------------------------------- */
/* Subroutine */ int grad_sum( int *nlocal,
			      Pme2Particle *ParticlePtr, 
	double *recip, double *theta1, double *theta2, double 
	*theta3, double *dtheta1, double *dtheta2, double *
	dtheta3, PmeVector *rfparticle, double *
	fr1, double *fr2, double *fr3,  int *order,  int *nfft1,
	  int *nfft2,  int *nfft3,  int *nfftdim1,  int *nfftdim2,
	  int *nfftdim3, double *q)
{
    /* System generated locals */
    int theta1_dim1, theta2_dim1, theta3_dim1, dtheta1_dim1,
	    dtheta2_dim1, dtheta3_dim1, q_dim2, q_dim3, q_offset;

    /* Local variables */
    double term;
    int i, j, k, n,nn;
    double f1, f2;
    int i0, j0, k0;
    double f3;
    int ith1, ith2, ith3;
    int index1, index2,index3,ifr;

    /* Parameter adjustments */
   
/* july 4 96    recip -= 4; */
    theta1_dim1 = *order;
    theta1 -= *order +1 ;
    theta2_dim1 = *order;
    theta2 -= *order +1 ;
    theta3_dim1 = *order;
    theta3 -=  *order +1 ;

    dtheta1_dim1 = *order;
    dtheta1 -=  *order +1 ;
    dtheta2_dim1 = *order;
    dtheta2 -=  *order +1 ;
    dtheta3_dim1 = *order;
    dtheta3 -=  *order +1 ;
  
  

    q_dim2 = *nfftdim1;
    q_dim3 = *nfftdim2;
    q_offset = (q_dim2 * (q_dim3 + 1) + 1 << 1) + 1;
    q -= q_offset;

#if 0
    /* this is NOT CORRECT because it doesn't use list[] atom pointers */

    /* do so because I start n-atoms loop counter from 1*/
    /* u can do fr(n-1) if u wish and not do this */
    --fr1; 
    --fr2;
    --fr3; 
    
    /* each processor calc recip-force on its own set */
    for (n = 1; n <= (*nlocal); ++n) {
      f1 = 0.;
      f2 = 0.;
      f3 = 0.;
      k0 = ( int) fr3[n] - *order;
      for (ith3 = 1; ith3 <= ( *order); ++ith3) {
	++k0;
	k = k0 + 1 + (*nfft3 - Nsign(*nfft3, k0)) / 2;
	j0 = ( int) fr2[n] - *order;
	
	for (ith2 = 1; ith2 <=(*order); ++ith2) {
	  ++j0;
	  j = j0 + 1 + (*nfft2 - Nsign(*nfft2, j0)) / 2;
	  i0 = ( int) fr1[n] - *order;

	  for (ith1 = 1; ith1 <= (*order); ++ith1) {
	    ++i0;
	    i = i0 + 1 + (*nfft1 - Nsign(*nfft1, i0)) / 2;
	    term = ParticlePtr[n].cg * q[(i + (j + k * q_dim3) * q_dim2 << 1) 
					 + 1];
	    /* force is negative of grad */
	    f1 -= *nfft1 * term * dtheta1[ith1 + n * dtheta1_dim1] * 
	      theta2[ith2 + n * theta2_dim1] * theta3[ith3 + n *
						      theta3_dim1];
	    f2 -= *nfft2 * term * theta1[ith1 + n * theta1_dim1] * 
	      dtheta2[ith2 + n * dtheta2_dim1] * theta3[ith3 + 
							n * theta3_dim1];
	    f3 -= *nfft3 * term * theta1[ith1 + n * theta1_dim1] * 
	      theta2[ith2 + n * theta2_dim1] * dtheta3[ith3 + n 
						       * dtheta3_dim1];
	  }
	}
      }
      /* this was n-1 now corrected dec 96 */
      rfparticle[n].x += recip[0] * f1 + recip[3] * f2 + recip[6] * f3;
      rfparticle[n].y += recip[1] * f1 + recip[4] * f2 + recip[7] * f3;
      rfparticle[n].z += recip[2] * f1 + recip[5] * f2 + recip[8] * f3;
      }
#endif

#if 1

#if DPME_DEBUG
    printf("grad_sum MODIFIED..\n");
#endif

    nn=1; /* ayt 3/97 made this 1 to comply with master/worker model */

    ifr=0; /* index the fr123 arrays = n-1 */
    /* each processor calc recip-force on its own set */
    for (n = 1; n <= (*nlocal); ++n, ++ifr, ++nn) {
      f1 = 0.;
      f2 = 0.;
      f3 = 0.;
      k0 = ( int) fr3[ifr] - *order;
      for (ith3 = 1; ith3 <= ( *order); ++ith3) {
	++k0;
	k = k0 + 1 + (*nfft3 - Nsign(*nfft3, k0)) / 2;
	j0 = ( int) fr2[ifr] - *order;
	index3= ith3 + n * (*order) ;

	for (ith2 = 1; ith2 <=(*order); ++ith2) {
	  ++j0;
	  j = j0 + 1 + (*nfft2 - Nsign(*nfft2, j0)) / 2;
	  i0 = ( int) fr1[ifr] - *order;
	  index2 = ith2 + n * (*order) ;

	  for (ith1 = 1; ith1 <= (*order); ++ith1) {
	    ++i0;
	    i = i0 + 1 + (*nfft1 - Nsign(*nfft1, i0)) / 2;
	    term = ParticlePtr[nn].cg * q[(i + (j + k * q_dim3) * q_dim2 << 1) 
					 + 1];
	    index1=  ith1 + n * (*order) ;
	    /* force is negative of grad */
	    f1 -= *nfft1 * term * dtheta1[index1] * theta2[index2] * theta3[index3];
	    f2 -= *nfft2 * term * theta1[index1] * dtheta2[index2] * theta3[index3];
	    f3 -= *nfft3 * term * theta1[index1] * theta2[index2] * dtheta3[index3];
	  }
	}
      }

      /* this was n-1 now corrected dec 96 */
      rfparticle[nn].x += recip[0] * f1 + recip[3] * f2 + recip[6] * f3;
      rfparticle[nn].y += recip[1] * f1 + recip[4] * f2 + recip[7] * f3;
      rfparticle[nn].z += recip[2] * f1 + recip[5] * f2 + recip[8] * f3;
      /* printf("n=%d nn=%d recipF1 (%3.4f,%3.4f,%3.4f) \n",n,nn,
       * rfparticle[nn].x, rfparticle[nn].y,  rfparticle[nn].z);
       */
    }
#endif /* modfied grad sum */
    return 0;
  } /* grad_sum */










