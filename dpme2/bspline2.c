/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */

/* --------------------------------------------------------------------- */
/* Subroutine */ int get_bspline_coeffs2( int *numatoms, double *fr1, 
	double *fr2, double *fr3,  int *order, double *theta1, 
	double *theta2, double *theta3, double *dtheta1, 
	double *dtheta2, double *dtheta3)
{
/* ayt 4/96 I have corrected the fr123 indexes but not the theta/dtheta123
 * this means that u have to declare fr's with dvector2 and theta/dtheta
 * with dvector()
 */ 
/* ---------------------------------------------------------------- */
/* INPUT: */
/*      numatoms: number of atoms */
/*      fr1,fr2,fr3 the scaled and shifted fractional coords */
/*      order: the order of spline interpolation */
/* OUTPUT */
/*      theta1,theta2,theta3: the spline coeff arrays */
/*      dtheta1,dtheta2,dtheta3: the 1st deriv of spline coeff arrays */
/* --------------------------------------------------------------------- 
*/
    /* System generated locals */
     int theta1_dim1, theta1_offset, theta2_dim1, theta2_offset, 
	    theta3_dim1, theta3_offset, dtheta1_dim1, dtheta1_offset, 
	    dtheta2_dim1, dtheta2_offset, dtheta3_dim1, dtheta3_offset;

    /* Local variables */
    extern /* Subroutine */ int fill_bspline(double *,  int *, 
	    double *, double *);
     int n,nn;
     double w;
     /* took out the --fr123 stuff 4/96 */
     theta1_dim1 = *order;
     theta1_offset = theta1_dim1 + 1;
     theta1 -= theta1_offset;
     theta2_dim1 = *order;
     theta2_offset = theta2_dim1 + 1;
     theta2 -= theta2_offset;
     theta3_dim1 = *order;
     theta3_offset = theta3_dim1 + 1;
     theta3 -= theta3_offset;
     dtheta1_dim1 = *order;
     dtheta1_offset = dtheta1_dim1 + 1;
     dtheta1 -= dtheta1_offset;
     dtheta2_dim1 = *order;
     dtheta2_offset = dtheta2_dim1 + 1;
     dtheta2 -= dtheta2_offset;
     dtheta3_dim1 = *order;
     dtheta3_offset = dtheta3_dim1 + 1;
     dtheta3 -= dtheta3_offset;

     nn=0;
     for (n = 1; n <= (*numatoms); ++n,++nn) {
       w = fr1[nn] - (int) fr1[nn];

       fill_bspline(&w, order, &theta1[n * theta1_dim1 + 1], 
		    &dtheta1[n * dtheta1_dim1 + 1]);
       w = fr2[nn] - (int) fr2[nn];
       fill_bspline(&w, order, &theta2[n * theta2_dim1 + 1], 
		    &dtheta2[n *  dtheta2_dim1 + 1]);
       w = fr3[nn] - (int) fr3[nn];
       fill_bspline(&w, order, &theta3[n * theta3_dim1 + 1], 
		    &dtheta3[n * dtheta3_dim1 + 1]);
     }
     return 0;
} /* get_bspline_coeffs */

/* --------------------------------------------------- */
/* Subroutine */ int fill_bspline(double *w,   int *order, double 
	*array, double *darray)
{
  /* Local variables */
    extern /* Subroutine */ int diff(double *, double *,  int *), 
	    init(double *, double *,  int *), one_pass(
	    double *, double *,  int *);
    int k;

    /* Parameter adjustments */
    --array;
    --darray;

    /* Function Body */
/* ---------- use standard B-spline recursions: see doc file */
/* do linear case */
    init(&array[1], w, order);
/* compute standard b-spline recursion */
    
    for (k = 3; k <= ( *order - 1); ++k) {
	one_pass(&array[1], w, &k);
/* L10: */
    }
/* perform standard b-spline differentiation */
    diff(&array[1], &darray[1], order);
/* one more recursion */
    one_pass(&array[1], w, order);
    return 0;
} /* fill_bspline */

/* --------------------------------------------------- */
/* Subroutine */ int init(double *c, double *x,  int *order)
{
    /* Parameter adjustments */
    --c;

    /* Function Body */
    c[*order] = 0.;
    c[2] = *x;
    c[1] = 1. - *x;
    return 0;
} /* init_ */

/* ------------------------------------- */
/* Subroutine */ int one_pass(double *c, double *x,  int *k)
{
  /* Local variables */
    int j;
    double div;

    /* Parameter adjustments */
    --c;

    /* Function Body */
    div = 1. / (*k - 1);
    c[*k] = div * *x * c[*k - 1];
    
    for (j = 1; j <= ( *k - 2); ++j) {
	c[*k - j] = div * ((*x + j) * c[*k - j - 1] + (*k - j - *x) * c[*k - j]);
    }
    c[1] = div * (1 - *x) * c[1];
    return 0;
} /* one_pass */

/* ------------------------------------- */
/* Subroutine */ int diff(double *c, double *d,  int *order)
{
  /* Local variables */
     int j;

    /* Parameter adjustments */
    --c;
    --d;

    /* Function Body */
    d[1] = -c[1];
    
    for (j = 2; j <= (*order); ++j) {
	d[j] = c[j - 1] - c[j];
/* L10: */
    }
    return 0;
} /* diff_ */



