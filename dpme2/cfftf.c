
/* Subroutine */ int cfftf(int *n, double *c, double *wsave)
{
    extern /* Subroutine */ int cfftf1(int *, double *, double *,
	     double *, double *);
    static int iw1, iw2;

    /* Parameter adjustments */
    --c;
    --wsave;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    iw1 = *n + *n + 1;
    iw2 = iw1 + *n + *n;
    cfftf1(n, &c[1], &wsave[1], &wsave[iw1], &wsave[iw2]);
    return 0;
} /* cfftf_ */

