/* use this if you have fortran PME and no erfc() call */
double  erfcfun(double *x)
{
  double erfc(double);
  return (erfc(*x));
	
}
