
/* Subroutine */ int cfftb1(int *n, double *c, double *ch, 
	double *wa, int *ifac)
{
    /* System generated locals */
    int i_1;

    /* Local variables */
    static int idot, i;
    extern /* Subroutine */ int passb(int *, int *, int *, 
	    int *, int *, double *, double *, double *, 
	    double *, double *, double *);
    static int k1, l1, l2, n2;
    extern /* Subroutine */ int passb2(int *, int *, double *, 
	    double *, double *), passb3(int *, int *, 
	    double *, double *, double *, double *), passb4_(
	    int *, int *, double *, double *, double *, 
	    double *, double *), passb5(int *, int *, 
	    double *, double *, double *, double *, 
	    double *, double *);
    static int na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1;

    /* Parameter adjustments */
    --c;
    --ch;
    --wa;
    --ifac;

    /* Function Body */
    nf = ifac[2];
    na = 0;
    l1 = 1;
    iw = 1;
    i_1 = nf;
    for (k1 = 1; k1 <= i_1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	if (na != 0) {
	    goto L101;
	}
	passb4(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	passb4(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	passb2(&idot, &l1, &c[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	passb2(&idot, &l1, &ch[1], &c[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + idot;
	if (na != 0) {
	    goto L107;
	}
	passb3(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	passb3(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	ix4 = ix3 + idot;
	if (na != 0) {
	    goto L110;
	}
	passb5(&idot, &l1, &c[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
	goto L111;
L110:
	passb5(&idot, &l1, &ch[1], &c[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	passb(&nac, &idot, &ip, &l1, &idl1, &c[1], &c[1], &c[1], &ch[1], &ch[
		1], &wa[iw]);
	goto L114;
L113:
	passb(&nac, &idot, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c[1], &
		c[1], &wa[iw]);
L114:
	if (nac != 0) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * idot;
/* L116: */
    }
    if (na == 0) {
	return 0;
    }
    n2 = *n + *n;
    i_1 = n2;
    for (i = 1; i <= i_1; ++i) {
	c[i] = ch[i];
/* L117: */
    }
    return 0;
} /* cfftb1_ */

