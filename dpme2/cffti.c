/* Subroutine */ int cffti(int *n, double *wsave)
{
    extern /* Subroutine */ int cffti1(int *, double *, double *)
	    ;
    static int iw1, iw2;

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    iw1 = *n + *n + 1;
    iw2 = iw1 + *n + *n;
    cffti1(n, &wsave[iw1], &wsave[iw2]);
    return 0;
} /* cffti_ */

