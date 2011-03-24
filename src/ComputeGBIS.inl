/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef FS_MAX
#define FS_MAX 1.728 /*maximum screened radius*/
#endif

#ifndef COULOMB
#define COULOMB 332.0636 /* ke [kcal*Ang/e^2] common.h */
#endif

#ifndef TA
#define TA 0.333333333333333 // 1/3
#define TB 0.4               // 2/5
#define TC 0.428571428571428 // 3/7
#define TD 0.444444444444444 // 4/9
#define TE 0.454545454545454 // 5/11
#define DA 1.333333333333333 // 4* 1/3
#define DB 2.4               // 6* 2/5
#define DC 3.428571428571428 // 8* 3/7
#define DD 4.444444444444444 // 10*4/9
#define DE 5.454545454545454 // 12*5/11
#endif

inline BigReal FastTanH( BigReal x ) {
  BigReal a = 2*x+0.02;;
  a *= (6 + a*(3 + a));
  return (a)/(a+12);
}

/******************************************************************************
 * use mass to determine element
 * return Bondi radius
 * returns the radius of i, depends on j though
 ******************************************************************************/
inline Real MassToRadius(Mass mi) {//, Mass mj) {
return
  (mi <   2.50 ) ? 1.20 : // AtmNum = 1;  Elem =  H ; Mass =   1.00
  (mi <   5.47 ) ? 1.40 : // AtmNum = 2;  Elem =  He; Mass =   4.00
  (mi <   7.98 ) ? 1.82 : // AtmNum = 3;  Elem =  Li; Mass =   6.94
  (mi <   9.91 ) ? 2.13 : // AtmNum = 4;  Elem =  Be; Mass =   9.01
  (mi <  11.41 ) ? 2.13 : // AtmNum = 5;  Elem =  B ; Mass =  10.81
  (mi <  13.01 ) ? 1.70 : // AtmNum = 6;  Elem =  C ; Mass =  12.01
  (mi <  15.00 ) ? 1.55 : // AtmNum = 7;  Elem =  N ; Mass =  14.00
  (mi <  17.49 ) ? 1.50 : // AtmNum = 8;  Elem =  O ; Mass =  15.99
  (mi <  19.58 ) ? 1.50 : // AtmNum = 9;  Elem =  F ; Mass =  18.99
  (mi <  21.58 ) ? 1.54 : // AtmNum = 10; Elem =  Ne; Mass =  20.17
  (mi <  23.64 ) ? 2.27 : // AtmNum = 11; Elem =  Na; Mass =  22.98
  (mi <  25.64 ) ? 1.73 : // AtmNum = 12; Elem =  Mg; Mass =  24.30
  (mi <  27.53 ) ? 2.51 : // AtmNum = 13; Elem =  Al; Mass =  26.98
  (mi <  29.53 ) ? 2.10 : // AtmNum = 14; Elem =  Si; Mass =  28.08
  (mi <  31.52 ) ? 1.85 : // AtmNum = 15; Elem =  P ; Mass =  30.97
  (mi <  33.76 ) ? 1.80 : // AtmNum = 16; Elem =  S ; Mass =  32.06
  (mi <  37.28 ) ? 1.70 : // AtmNum = 17; Elem =  Cl; Mass =  35.45
  (mi <  39.29 ) ? 2.75 : // AtmNum = 19; Elem =  K ; Mass =  39.10
  (mi <  49.09 ) ? 1.88 : // AtmNum = 18; Elem =  Ar; Mass =  39.48
  (mi <  61.12 ) ? 1.63 : // AtmNum = 28; Elem =  Ni; Mass =  58.69
  (mi <  64.46 ) ? 1.40 : // AtmNum = 29; Elem =  Cu; Mass =  63.54
  (mi <  67.55 ) ? 1.39 : // AtmNum = 30; Elem =  Zn; Mass =  65.38
  (mi <  71.18 ) ? 1.87 : // AtmNum = 31; Elem =  Ga; Mass =  69.72
  (mi <  73.78 ) ? 2.19 : // AtmNum = 32; Elem =  Ge; Mass =  72.64
  (mi <  76.94 ) ? 1.85 : // AtmNum = 33; Elem =  As; Mass =  74.92
  (mi <  79.43 ) ? 1.90 : // AtmNum = 34; Elem =  Se; Mass =  78.96
  (mi <  81.85 ) ? 1.85 : // AtmNum = 35; Elem =  Br; Mass =  79.90
  (mi <  95.11 ) ? 2.02 : // AtmNum = 36; Elem =  Kr; Mass =  83.79
  (mi < 107.14 ) ? 1.63 : // AtmNum = 46; Elem =  Pd; Mass = 106.42
  (mi < 110.14 ) ? 1.72 : // AtmNum = 47; Elem =  Ag; Mass = 107.86
  (mi < 113.61 ) ? 1.58 : // AtmNum = 48; Elem =  Cd; Mass = 112.41
  (mi < 116.76 ) ? 1.93 : // AtmNum = 49; Elem =  In; Mass = 114.81
  (mi < 120.24 ) ? 2.17 : // AtmNum = 50; Elem =  Sn; Mass = 118.71
  (mi < 124.33 ) ? 2.09 : // AtmNum = 51; Elem =  Sb; Mass = 121.76
  (mi < 127.25 ) ? 1.98 : // AtmNum = 53; Elem =  I ; Mass = 126.90
  (mi < 129.45 ) ? 2.06 : // AtmNum = 52; Elem =  Te; Mass = 127.60
  (mi < 163.19 ) ? 2.16 : // AtmNum = 54; Elem =  Xe; Mass = 131.29
  (mi < 196.02 ) ? 1.75 : // AtmNum = 78; Elem =  Pt; Mass = 195.08
  (mi < 198.78 ) ? 1.66 : // AtmNum = 79; Elem =  Au; Mass = 196.96
  (mi < 202.49 ) ? 1.55 : // AtmNum = 80; Elem =  Hg; Mass = 200.59
  (mi < 205.79 ) ? 1.96 : // AtmNum = 81; Elem =  Tl; Mass = 204.38
  (mi < 222.61 ) ? 2.02 : // AtmNum = 82; Elem =  Pb; Mass = 207.20
  (mi < 119.01 ) ? 1.86 : // AtmNum = 92; Elem =  U ; Mass = 238.02
                   1.50 ;   // Unknown
}

/******************************************************************************
 * Screen radii
 * use masses to determine elements
 * use elements to lookup Sij
 * to scale the coulomb radius
 * from Hawkins, Cramer, Truhlar; 1996
 * mi is descreened atom - calculating it's alpha (outer loop index)
 * mj is descreening atom - contributor (inner loop index)
 ******************************************************************************/
inline Real MassToScreen(Mass mi) {//, Mass mj) {
    return
      (mi <   1.500) ? 0.85 : //H
      (mi <  12.500) ? 0.72 : //C
      (mi <  14.500) ? 0.79 : //N
      (mi <  16.500) ? 0.85 : //O
      (mi <  19.500) ? 0.88 : //F
      (mi <  31.500) ? 0.86 : //P
      (mi <  32.500) ? 0.96 : //S
                        0.8 ; //all others
}


/******************************************************************************
  Piecewise screening functions Hij dHij/drij
  r   distance
  r2  square distance
  ri  inverse distance
  rc  cutoff
  r0  descreened atom radius
  rs  descreening atom radius
  h   return value
  dh  return value
 ******************************************************************************/
inline void h0 ( BigReal r, BigReal r2, BigReal ri,//0*5.3%
    Real rc, BigReal r0, BigReal rs, BigReal & h ) {
  h = 0;
}
inline void dh0 ( BigReal r, BigReal r2, BigReal ri,//0*5.3%
    Real rc, BigReal r0, BigReal rs, BigReal & dh ) {
  dh = 0;
}

inline void h1 ( BigReal r, BigReal r2, BigReal ri, //18.4%
    BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {

  BigReal rci = 1.0/rc;
  BigReal rmrs = r-rs;
  BigReal rmrsi = 1.0/rmrs;
  BigReal rmrs2 = rmrs*rmrs;
  BigReal rs2 = rs*rs;
  BigReal logr = log(rmrs*rci);
  BigReal rci2 = rci*rci;
  h = 0.125*ri*(1 + 2*r*rmrsi + rci2*(r2 - 4*rc*r - rs2) + 2*logr);
}
inline void dh1 ( BigReal r, BigReal r2, BigReal ri, //18.4%
    BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {

  BigReal rci = 1.0/rc;
  BigReal rmrs = r-rs;// 4 times
  BigReal rmrsi = 1.0/rmrs;
  BigReal rmrs2 = rmrs*rmrs;
  BigReal rs2 = rs*rs;
  BigReal logr = log(rmrs*rci);
  BigReal rci2 = rci*rci;
  dh = ri*ri*(-0.25*logr - (rc*rc - rmrs2)*(rs2 + r2)*0.125*rci2*rmrsi*rmrsi);
}

inline void h2 ( BigReal r, BigReal r2, BigReal ri,// 74.5%
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {

    BigReal k = rs*ri; k*=k;//k=(rs/r)^2
    h = rs*ri*ri*k*(TA+k*(TB+k*(TC+k*(TD+k*TE))));
}
inline void dh2 ( BigReal r, BigReal r2, BigReal ri,// 74.5%
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {

    BigReal k = rs*ri; k*=k;//k=(rs/r)^2
    dh = -rs*ri*ri*ri*k*(DA+k*(DB+k*(DC+k*(DD+k*DE))));
}

inline void h3 ( BigReal r, BigReal r2, BigReal ri,// 1.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {
    BigReal r2mrs2i = 1.0/(r2-rs*rs);
    h = 0.5 * ( rs*r2mrs2i + 0.5 * log((r-rs)/(r+rs))*ri );
}
inline void dh3 ( BigReal r, BigReal r2, BigReal ri,// 1.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {
    BigReal rs2 = rs*rs;
    BigReal r2mrs2i = 1.0/(r2-rs2);
    dh = -0.25*ri*(2*(r2+rs2)*rs*r2mrs2i*r2mrs2i + ri*log((r-rs)/(r+rs)));
}

inline void h4 ( BigReal r, BigReal r2, BigReal ri,// 0.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {
    BigReal ri2 = ri*ri;
    BigReal r02 = r0*r0;
    BigReal rs2 = rs*rs;
    BigReal r0i = 1.0/r0;
    BigReal rspri = 1.0/(r+rs);
    BigReal logr = log(r0*rspri);
    BigReal r02mrs2 = r02-rs2;
    BigReal rilogr = ri*logr;
    h = 0.25*( r0i*(2-0.5*(r0i*ri*(r2 + r02 - rs2))) - rspri + rilogr );
}
inline void dh4 ( BigReal r, BigReal r2, BigReal ri,// 0.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {
    BigReal ri2 = ri*ri;
    BigReal r02 = r0*r0;
    BigReal rs2 = rs*rs;
    BigReal r0i = 1.0/r0;
    BigReal rspri = 1.0/(r+rs);
    BigReal logr = log(r0*rspri);
    BigReal r02mrs2 = r02-rs2;
    BigReal rilogr = ri*logr;
    dh = 0.25*( (-0.5+(r2*r02mrs2 - 2*r*rs*rs2+rs2*r02mrs2)
        * 0.5*ri2*rspri*rspri)*r0i*r0i - ri*rilogr );
}

inline void h5 ( BigReal r, BigReal r2, BigReal ri,// 0%, r<0.7Ang
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {
    BigReal rs2 = rs*rs;
    BigReal r2mrs2i = 1/(r2-rs2);
    BigReal rsr2mrs2i = rs*r2mrs2i;
    BigReal rprs = r+rs;
    BigReal rmrs = r-rs;
    BigReal logr = 0.5*ri*log(-rmrs/rprs);
    h = 0.5*( rsr2mrs2i + 2/r0 + logr );
}
inline void dh5 ( BigReal r, BigReal r2, BigReal ri,// 0%, r<0.7Ang
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {
    BigReal rs2 = rs*rs;
    BigReal r2mrs2i = 1/(r2-rs2);
    BigReal rsr2mrs2i = rs*r2mrs2i;
    BigReal rprs = r+rs;
    BigReal rmrs = r-rs;
    BigReal logr = 0.5*ri*log(-rmrs/rprs);
    dh = -0.5*ri*((rs2+r2)*rsr2mrs2i*r2mrs2i+logr );
}

inline void h6 ( BigReal r, BigReal r2, BigReal ri,//0%, one atom within other
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {
  h = 0;
}
inline void dh6 ( BigReal r, BigReal r2, BigReal ri,//0%, one atom within other
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {
  dh = 0;
}
#if 0
inline void CalcH ( BigReal r, BigReal r2, BigReal ri,
BigReal rc, BigReal r0, BigReal rs, BigReal & h, int & d) {
  if( r <= rc - rs && r > 4*rs ) {//II
    h2(r,r2,ri,rc,r0,rs,h); d = 2;
  } else if (r <= rc + rs && r > rc - rs) {//I
    h1(r,r2,ri,rc,r0,rs,h); d = 1;
  } else if (r > rc + rs) {//0
    h0(r,r2,ri,rc,r0,rs,h); d = 0;
  } else if( r <= 4*rs && r > r0 + rs ) {//III
    h3(r,r2,ri,rc,r0,rs,h); d = 3;
  } else if ( r <= r0 + rs && r > (r0>rs?r0-rs:rs-r0) ) {//IV
    h4(r,r2,ri,rc,r0,rs,h); d = 4;
  } else if (r0 < rs ) {//V
    h5(r,r2,ri,rc,r0,rs,h); d = 5;
  } else {//VI
    h6(r,r2,ri,rc,r0,rs,h); d = 6;
  }
}
inline void CalcDH ( BigReal r, BigReal r2, BigReal ri,
BigReal rc, BigReal r0, BigReal rs, BigReal & dh, int & d) {
  if( r <= rc - rs && r > 4*rs ) {//II
    dh2(r,r2,ri,rc,r0,rs,dh); d = 2;
  } else if (r <= rc + rs && r > rc - rs) {//I
    dh1(r,r2,ri,rc,r0,rs,dh); d = 1;
  } else if (r > rc + rs) {//0
    dh0(r,r2,ri,rc,r0,rs,dh); d = 0;
  } else if( r <= 4*rs && r > r0 + rs ) {//III
    dh3(r,r2,ri,rc,r0,rs,dh); d = 3;
  } else if ( r <= r0 + rs && r > (r0>rs?r0-rs:rs-r0) ) {//IV
    dh4(r,r2,ri,rc,r0,rs,dh); d = 4;
  } else if (r0 < rs ) {//V
    dh5(r,r2,ri,rc,r0,rs,dh); d = 5;
  } else {//VI
    dh6(r,r2,ri,rc,r0,rs,dh); d = 6;
  }
}
#else
inline void CalcH ( BigReal r, BigReal r2, BigReal ri,
BigReal rc, BigReal r0, BigReal rs, BigReal & h, int & d) {

/*
r - distance
rc - alpha cutoff
rs - screened radius
*/

if (r > 4*rs) { //change this to 1/4 r > rs
  if( r < rc - rs) {//II            68%
    h2(r,r2,ri,rc,r0,rs,h); d = 2;
  } else if (r < rc + rs) {//I      23%
    h1(r,r2,ri,rc,r0,rs,h); d = 1;
  } else /*if (r > rc + rs)*/ {//0  7%
    h0(r,r2,ri,rc,r0,rs,h); d = 0;
  } 
} else {
  if( r > r0 + rs ) {//III          1%
    h3(r,r2,ri,rc,r0,rs,h); d = 3;
  } else if ( r > (r0>rs?r0-rs:rs-r0) ) {//IV 0%
    h4(r,r2,ri,rc,r0,rs,h); d = 4;
  } else if (r0 < rs ) {//V         0%
    h5(r,r2,ri,rc,r0,rs,h); d = 5;
  } else {//VI                      0%
    h6(r,r2,ri,rc,r0,rs,h); d = 6;
  }
}
}
inline void CalcDH ( BigReal r, BigReal r2, BigReal ri,
BigReal rc, BigReal r0, BigReal rs, BigReal & dh, int & d) {
if (r > 4*rs) {
  if( r < rc - rs) {//II
    dh2(r,r2,ri,rc,r0,rs,dh); d = 2;
  } else if (r < rc + rs) {//I
    dh1(r,r2,ri,rc,r0,rs,dh); d = 1;
  } else /*if (r > rc + rs)*/ {//0
    dh0(r,r2,ri,rc,r0,rs,dh); d = 0;
  }
} else {
  if( r > r0 + rs ) {//III
    dh3(r,r2,ri,rc,r0,rs,dh); d = 3;
  } else if (r > (r0>rs?r0-rs:rs-r0) ) {//IV
    dh4(r,r2,ri,rc,r0,rs,dh); d = 4;
  } else if (r0 < rs ) {//V
    dh5(r,r2,ri,rc,r0,rs,dh); d = 5;
  } else {//VI
    dh6(r,r2,ri,rc,r0,rs,dh); d = 6;
  }
}
}
#endif
inline void CalcHPair (
  BigReal r,//distance
  BigReal r2,//distance squared
  BigReal ri,//inverse distance
  BigReal rc,//cutoff
  BigReal ri0,
  BigReal rjs,
  BigReal rj0,
  BigReal ris,
  int & dij,//domain 1
  int & dji,//domain 2
  BigReal & hij,//output
  BigReal & hji//output
) {
  CalcH(r,r2,ri,rc,ri0,rjs,hij,dij);//hij
  CalcH(r,r2,ri,rc,rj0,ris,hji,dji);//hji
}
inline void CalcDHPair (
  BigReal r,//distance
  BigReal r2,
  BigReal ri,
  BigReal rc,//cutoff
  BigReal ri0,
  BigReal rjs,
  BigReal rj0,
  BigReal ris,
  int & dij,//domain 1
  int & dji,//domain 2
  BigReal & dhij,
  BigReal & dhji
) {
//                  swapped
  CalcDH(r,r2,ri,rc,ri0,rjs,dhij,dij);//hij
  CalcDH(r,r2,ri,rc,rj0,ris,dhji,dji);//hji
}

/*
 * Calculate GB Energy, GB dEdr force
 * also output intermediate values used in dEda
 */
inline void Calc_dEdr_Pair(//no longer does i==j
  const BigReal & r,
  const BigReal & r2,
  const BigReal & qiqj,
  const BigReal & ai,
  const BigReal & aj,
  const BigReal & kappa,
  const BigReal & epsilon_p_i,
  const BigReal & epsilon_s_i,
  BigReal & aiaj,
  BigReal & expr2aiaj4,
  BigReal & fij,
  BigReal & f_i,
  BigReal & expkappa,
  BigReal & Dij,
  BigReal & gbE,   //return
  BigReal & ddrGbE //return
) {
  //allocate local variables
  BigReal aiaj4,ddrDij,ddrf_i,ddrfij;

  //calculate GB energy
  aiaj = ai*aj;
  aiaj4 = 4*aiaj;
  //printf("exp(%e)\n",(-r2/aiaj4));
  expr2aiaj4 = exp(-r2/aiaj4);
  fij = sqrt(r2+aiaj*expr2aiaj4);
  f_i = 1/fij;
  expkappa = (kappa > 0) ? exp(-kappa*fij) : 1.0;
  Dij = epsilon_p_i - expkappa*epsilon_s_i;
  //gbE = -COULOMB*qiqj*Dij*f_i;
  gbE = qiqj*Dij*f_i;

  //calculate energy derivatives
  ddrfij = r*f_i*(1 - 0.25*expr2aiaj4);
  ddrf_i = -ddrfij*f_i*f_i;
  ddrDij = kappa*expkappa*ddrfij*epsilon_s_i;
  //ddrGbE = -COULOMB*qiqj*(ddrDij*f_i+Dij*ddrf_i);
  ddrGbE = qiqj*(ddrDij*f_i+Dij*ddrf_i);
}

/*
 * Calculate summation element of dEda array
 * must calculate dEdr previously to retreive intermediate values
 */
inline void Calc_dEda_Pair(
  const BigReal & r2,
  const BigReal & ai,
  const BigReal & aj,
  const BigReal & qiqj,
  const BigReal & kappa,
  const BigReal & aiaj,
  const BigReal & expkappa,
  const BigReal & expr2aiaj4,
  const BigReal & fij,
  const BigReal & f_i,
  const BigReal & Dij,
  const BigReal & epsilon_s_i,
  BigReal & dEdai,//return
  BigReal & dEdaj //return
) {

  //BigReal tmp_dEda = -0.5*COULOMB*qiqj*f_i*f_i
  BigReal tmp_dEda = 0.5*qiqj*f_i*f_i
                      *(kappa*epsilon_s_i*expkappa-Dij*f_i)
                      *(aiaj+0.25*r2)*expr2aiaj4;//0
  dEdai = tmp_dEda/ai;
  dEdaj = tmp_dEda/aj;
}

/*
 * Calculate Coulomb and GB interaction and dEda element
 * for a pair of atoms
 */
inline void Phase2_Pair(//doesn't do self energies

//input values
  const BigReal & r,
  const BigReal & r2,
  const BigReal & r_i,
  const BigReal & qiqj,
  const BigReal & ai,
  const BigReal & aj,
  const BigReal & epsilon_p_i,
  const BigReal & epsilon_s_i,
  const BigReal & kappa,
  const int & doFullElect,

//return values
  BigReal & gbEij,
  BigReal & ddrGbEij,
  BigReal & dEdai,
  BigReal & dEdaj
) {

  //calculate GB energy and force
  BigReal aiaj,expr2aiaj4,fij,f_i,expkappa,Dij;
  Calc_dEdr_Pair(r,r2,qiqj,ai,aj,
      kappa,epsilon_p_i,epsilon_s_i,
      aiaj,expr2aiaj4,fij,f_i,expkappa,
      Dij,gbEij,ddrGbEij);

  //calculate dEda
  if (doFullElect) {
    Calc_dEda_Pair(r2,ai,aj,qiqj,kappa,
              aiaj,expkappa,expr2aiaj4,
              fij,f_i,Dij,epsilon_s_i,dEdai,dEdaj);
  } else {
    dEdai = 0.0;
    dEdaj = 0.0;
  }
}

inline void init_gbisTable (
BigReal **tablePtr,
BigReal kappa,
BigReal maxX,
int numEntriesPerX
) {
  BigReal *table = *tablePtr;//dereference
  BigReal minX = 0;
  int numPts = (maxX-minX) * numEntriesPerX;
  int numVals = 3;
  /*
  table = (BigReal*) malloc(numVals*numPts*sizeof(BigReal));
  for (int i = 0; i < numPts; i++) {
    BigReal x = (1.0*i) / numEntriesPerX + minX; 
      bornRadJ = gbisParams->bornRad[1][j];
      aiaj = bornRadI * bornRadJ;
      aiaj4 = 4.0 * aiaj;
      expr2aiaj4 = exp(-r2/aiaj4);
      fij = sqrt(r2+aiaj*expr2aiaj4);
      f_i = 1.0/fij;
      expkappa = (kappa > 0.0) ? exp(-kappa*fij) : 1.0;
  }
*/
}
