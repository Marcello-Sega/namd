#include "f2c.h"
/****************************************************************************
**/

/* 	3D (slow) Fourier Transform */
/*   this 1d->3d code is brute force approach */
/*   the 1d code is a double precision version of fftpack from netlib */
/*   due to Paul N Swartztrauber at NCAR Boulder Coloraso */

/****************************************************************************
**/
/* Subroutine */ int pubz3di(int *n1, int *n2, int *n3, 
	double *table, int *ntable)
{
    /* System generated locals */
    int table_dim1, table_offset;

    /* Local variables */
    extern /* Subroutine */ int cffti(int *, double *);

    /* Parameter adjustments */
    table_dim1 = *ntable;
    table_offset = table_dim1 + 1;
    table -= table_offset;

    /* Function Body */
/* ntable should be 4*max(n1,n2,n3) +15 */
    cffti(n1, &table[table_dim1 + 1]);
    cffti(n2, &table[(table_dim1 << 1) + 1]);
    cffti(n3, &table[table_dim1 * 3 + 1]);
    return 0;
} /* pubz3di_ */

/****************************************************************************
**/
/* Subroutine */ int pubz3d(int *isign, int *n1, int *n2, 
	int *n3, doublecomplex *w, int *ld1, int *ld2, double 
	*table, int *ntable, doublecomplex *work, int * DUMMY1 /* nwork */)
{
    /* System generated locals */
    int w_dim1, w_dim2, w_offset, table_dim1, table_offset, i_1, i_2, i_3,
	     i_4, i_5;

    /* Local variables */
    static int i, j, k;
    extern /* Subroutine */ int cfftb(int *, doublecomplex *, double 
	    *), cfftf(int *, doublecomplex *, double *);

    /* Parameter adjustments */
    w_dim1 = *ld1;
    w_dim2 = *ld2;
    w_offset = w_dim1 * (w_dim2 + 1) + 1;
    w -= w_offset;
    table_dim1 = *ntable;
    table_offset = table_dim1 + 1;
    table -= table_offset;
    --work;

    /* Function Body */
/* ntable should be 4*max(n1,n2,n3) +15 */
/* nwork should be max(n1,n2,n3) */

/*   transform along X  first ... */

    i_1 = *n3;
    for (k = 1; k <= i_1; ++k) {
	i_2 = *n2;
	for (j = 1; j <= i_2; ++j) {
	    i_3 = *n1;
	    for (i = 1; i <= i_3; ++i) {
		i_4 = i;
		i_5 = i + (j + k * w_dim2) * w_dim1;
		work[i_4].r = w[i_5].r, work[i_4].i = w[i_5].i;
/* L70: */
	    }
	    if (*isign == -1) {
		cfftf(n1, &work[1], &table[table_dim1 + 1]);
	    }
	    if (*isign == 1) {
		cfftb(n1, &work[1], &table[table_dim1 + 1]);
	    }
	    i_3 = *n1;
	    for (i = 1; i <= i_3; ++i) {
		i_4 = i + (j + k * w_dim2) * w_dim1;
		i_5 = i;
		w[i_4].r = work[i_5].r, w[i_4].i = work[i_5].i;
/* L80: */
	    }
/* L90: */
	}
/* L100: */
    }

/*   transform along Y then ... */

    i_1 = *n3;
    for (k = 1; k <= i_1; ++k) {
	i_2 = *n1;
	for (i = 1; i <= i_2; ++i) {
	    i_3 = *n2;
	    for (j = 1; j <= i_3; ++j) {
		i_4 = j;
		i_5 = i + (j + k * w_dim2) * w_dim1;
		work[i_4].r = w[i_5].r, work[i_4].i = w[i_5].i;
/* L170: */
	    }
	    if (*isign == -1) {
		cfftf(n2, &work[1], &table[(table_dim1 << 1) + 1]);
	    }
	    if (*isign == 1) {
		cfftb(n2, &work[1], &table[(table_dim1 << 1) + 1]);
	    }
	    i_3 = *n2;
	    for (j = 1; j <= i_3; ++j) {
		i_4 = i + (j + k * w_dim2) * w_dim1;
		i_5 = j;
		w[i_4].r = work[i_5].r, w[i_4].i = work[i_5].i;
/* L180: */
	    }
/* L190: */
	}
/* L200: */
    }

/*   transform along Z finally ... */

    i_1 = *n1;
    for (i = 1; i <= i_1; ++i) {
	i_2 = *n2;
	for (j = 1; j <= i_2; ++j) {
	    i_3 = *n3;
	    for (k = 1; k <= i_3; ++k) {
		i_4 = k;
		i_5 = i + (j + k * w_dim2) * w_dim1;
		work[i_4].r = w[i_5].r, work[i_4].i = w[i_5].i;
/* L270: */
	    }
	    if (*isign == -1) {
		cfftf(n3, &work[1], &table[table_dim1 * 3 + 1]);
	    }
	    if (*isign == 1) {
		cfftb(n3, &work[1], &table[table_dim1 * 3 + 1]);
	    }
	    i_3 = *n3;
	    for (k = 1; k <= i_3; ++k) {
		i_4 = i + (j + k * w_dim2) * w_dim1;
		i_5 = k;
		w[i_4].r = work[i_5].r, w[i_4].i = work[i_5].i;
/* L280: */
	    }
/* L290: */
	}
/* L300: */
    }
    return 0;
} /* pubz3d_ */

