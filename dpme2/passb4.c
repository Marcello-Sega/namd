/* Subroutine */ int passb4(int *ido, int *l1, double *cc, 
	double *ch, double *wa1, double *wa2, double *wa3)
{
    /* System generated locals */
    int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i_1, i_2;

    /* Local variables */
    static int i, k;
    static double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, 
	    tr2, tr3, tr4;

    /* Parameter adjustments */
    cc_dim1 = *ido;
    cc_offset = cc_dim1 * 5 + 1;
    cc -= cc_offset;
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = ch_dim1 * (ch_dim2 + 1) + 1;
    ch -= ch_offset;
    --wa1;
    --wa2;
    --wa3;

    /* Function Body */
    if (*ido != 2) {
	goto L102;
    }
    i_1 = *l1;
    for (k = 1; k <= i_1; ++k) {
	ti1 = cc[((k << 2) + 1) * cc_dim1 + 2] - cc[((k << 2) + 3) * cc_dim1 
		+ 2];
	ti2 = cc[((k << 2) + 1) * cc_dim1 + 2] + cc[((k << 2) + 3) * cc_dim1 
		+ 2];
	tr4 = cc[((k << 2) + 4) * cc_dim1 + 2] - cc[((k << 2) + 2) * cc_dim1 
		+ 2];
	ti3 = cc[((k << 2) + 2) * cc_dim1 + 2] + cc[((k << 2) + 4) * cc_dim1 
		+ 2];
	tr1 = cc[((k << 2) + 1) * cc_dim1 + 1] - cc[((k << 2) + 3) * cc_dim1 
		+ 1];
	tr2 = cc[((k << 2) + 1) * cc_dim1 + 1] + cc[((k << 2) + 3) * cc_dim1 
		+ 1];
	ti4 = cc[((k << 2) + 2) * cc_dim1 + 1] - cc[((k << 2) + 4) * cc_dim1 
		+ 1];
	tr3 = cc[((k << 2) + 2) * cc_dim1 + 1] + cc[((k << 2) + 4) * cc_dim1 
		+ 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = tr2 + tr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = tr2 - tr3;
	ch[(k + ch_dim2) * ch_dim1 + 2] = ti2 + ti3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ti2 - ti3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = tr1 + tr4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = tr1 - tr4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ti1 + ti4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ti1 - ti4;
/* L101: */
    }
    return 0;
L102:
    i_1 = *l1;
    for (k = 1; k <= i_1; ++k) {
	i_2 = *ido;
	for (i = 2; i <= i_2; i += 2) {
	    ti1 = cc[i + ((k << 2) + 1) * cc_dim1] - cc[i + ((k << 2) + 3) * 
		    cc_dim1];
	    ti2 = cc[i + ((k << 2) + 1) * cc_dim1] + cc[i + ((k << 2) + 3) * 
		    cc_dim1];
	    ti3 = cc[i + ((k << 2) + 2) * cc_dim1] + cc[i + ((k << 2) + 4) * 
		    cc_dim1];
	    tr4 = cc[i + ((k << 2) + 4) * cc_dim1] - cc[i + ((k << 2) + 2) * 
		    cc_dim1];
	    tr1 = cc[i - 1 + ((k << 2) + 1) * cc_dim1] - cc[i - 1 + ((k << 2) 
		    + 3) * cc_dim1];
	    tr2 = cc[i - 1 + ((k << 2) + 1) * cc_dim1] + cc[i - 1 + ((k << 2) 
		    + 3) * cc_dim1];
	    ti4 = cc[i - 1 + ((k << 2) + 2) * cc_dim1] - cc[i - 1 + ((k << 2) 
		    + 4) * cc_dim1];
	    tr3 = cc[i - 1 + ((k << 2) + 2) * cc_dim1] + cc[i - 1 + ((k << 2) 
		    + 4) * cc_dim1];
	    ch[i - 1 + (k + ch_dim2) * ch_dim1] = tr2 + tr3;
	    cr3 = tr2 - tr3;
	    ch[i + (k + ch_dim2) * ch_dim1] = ti2 + ti3;
	    ci3 = ti2 - ti3;
	    cr2 = tr1 + tr4;
	    cr4 = tr1 - tr4;
	    ci2 = ti1 + ti4;
	    ci4 = ti1 - ti4;
	    ch[i - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * cr2 - 
		    wa1[i] * ci2;
	    ch[i + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i - 1] * ci2 + wa1[i]
		     * cr2;
	    ch[i - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * cr3 - wa2[
		    i] * ci3;
	    ch[i + (k + ch_dim2 * 3) * ch_dim1] = wa2[i - 1] * ci3 + wa2[i] * 
		    cr3;
	    ch[i - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * cr4 - 
		    wa3[i] * ci4;
	    ch[i + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i - 1] * ci4 + wa3[i]
		     * cr4;
/* L103: */
	}
/* L104: */
    }
    return 0;
} /* passb4_ */

