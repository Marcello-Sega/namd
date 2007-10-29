/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Settle.h"
#include <string.h>
#include <math.h>
#include <charm++.h> // for CkPrintf


//
// XXX static and global variables are unsafe for shared memory builds.
// The global and static vars should be eliminated.
// Unfortunately, the routines that use these below are actually
// in use in NAMD.
//

static BigReal mOrmT, mHrmT;
static BigReal ra, rb, rc;
static BigReal rra;
static int settle_first_time = 1;


int settle1isinitted(void) {
  return !settle_first_time;
}

// Initialize various properties of the waters
// settle1() assumes all waters are identical, 
// and will generate bad results if they are not.
int settle1init(BigReal pmO, BigReal pmH, BigReal hhdist, BigReal ohdist) {
  if (settle_first_time) {
    settle_first_time = 0;
    BigReal rmT = 1.0 / (pmO+pmH+pmH);
    mOrmT = pmO * rmT;
    mHrmT = pmH * rmT;
    BigReal t1 = 0.5*pmO/pmH;
    rc = 0.5*hhdist;
    ra = sqrt(ohdist*ohdist-rc*rc)/(1.0+t1);
    rb = t1*ra;
    rra = 1.0 / ra;
  }

  return 0;
}


int settle1(const Vector *ref, Vector *pos, Vector *vel, BigReal invdt) {
  // vectors in the plane of the original positions
  Vector b0 = ref[1]-ref[0];
  Vector c0 = ref[2]-ref[0];
  
  // new center of mass
  Vector d0 = pos[0]*mOrmT + ((pos[1] + pos[2])*mHrmT);
 
  Vector a1 = pos[0] - d0;
  Vector b1 = pos[1] - d0;
  Vector c1 = pos[2] - d0;
  
  // Vectors describing transformation from original coordinate system to
  // the 'primed' coordinate system as in the diagram.  
  Vector n0 = cross(b0, c0);
  Vector n1 = cross(a1, n0); 
  Vector n2 = cross(n0, n1); 
  n0 = n0.unit();
  n1 = n1.unit();
  n2 = n2.unit();

  // this is wasteful, as b0.z is never referenced again
  b0 = Vector(n1*b0, n2*b0, n0*b0);
  // this is wasteful, as c0.z is never referenced again
  c0 = Vector(n1*c0, n2*c0, n0*c0);
 
  // this is wasteful, as a1.x and a1.y are never referenced again
  a1 = Vector(n1*a1, n2*a1, n0*a1);
  b1 = Vector(n1*b1, n2*b1, n0*b1);
  c1 = Vector(n1*c1, n2*c1, n0*c1);

  // now we can compute positions of canonical water 
  BigReal sinphi = a1.z * rra;
  BigReal tmp = 1.0-sinphi*sinphi;
  BigReal cosphi = sqrt(tmp);
  BigReal sinpsi = (b1.z - c1.z)/(2.0*rc*cosphi);
  tmp = 1.0-sinpsi*sinpsi;
  BigReal cospsi = sqrt(tmp);

  BigReal rbphi = -rb*cosphi;
  BigReal tmp1 = rc*sinpsi*sinphi;
  BigReal tmp2 = rc*sinpsi*cosphi;
 
  Vector a2(0, ra*cosphi, ra*sinphi);
  Vector b2(-rc*cospsi, rbphi - tmp1, -rb*sinphi + tmp2);
  Vector c2( rc*cosphi, rbphi + tmp1, -rb*sinphi - tmp2);

  // there are no a0 terms because we've already subtracted the term off 
  // when we first defined b0 and c0.
  BigReal alpha = b2.x*(b0.x - c0.x) + b0.y*b2.y + c0.y*c2.y;
  BigReal beta  = b2.x*(c0.y - b0.y) + b0.x*b2.y + c0.x*c2.y;
  BigReal gama  = b0.x*b1.y - b1.x*b0.y + c0.x*c1.y - c1.x*c0.y;
 
  BigReal a2b2 = alpha*alpha + beta*beta;
  BigReal sintheta = (alpha*gama - beta*sqrt(a2b2 - gama*gama))/a2b2;
  BigReal costheta = sqrt(1.0 - sintheta*sintheta);
  
#if 0
  Vector a3( -a2.y*sintheta, 
              a2.y*costheta,
              a2.z);
  Vector b3(b2.x*costheta - b2.y*sintheta,
              b2.x*sintheta + b2.y*costheta,
              b2.z);
  Vector c3(c2.x*costheta - c2.y*sintheta,
              c2.x*sintheta + c2.y*costheta,
              c2.z);
  
#else
  Vector a3( -a2.y*sintheta, 
              a2.y*costheta,
              a1.z);
  Vector b3(b2.x*costheta - b2.y*sintheta,
              b2.x*sintheta + b2.y*costheta,
              b1.z);
  Vector c3(-b2.x*costheta - c2.y*sintheta,
            -b2.x*sintheta + c2.y*costheta,
              c1.z);

#endif

  // undo the transformation; generate new normal vectors from the transpose.
  Vector m1(n1.x, n2.x, n0.x);
  Vector m2(n1.y, n2.y, n0.y);
  Vector m0(n1.z, n2.z, n0.z);

  pos[0] = Vector(a3*m1, a3*m2, a3*m0) + d0;
  pos[1] = Vector(b3*m1, b3*m2, b3*m0) + d0;
  pos[2] = Vector(c3*m1, c3*m2, c3*m0) + d0;

  // dt can be negative during startup!
  if (invdt != 0) {
    vel[0] = (pos[0]-ref[0])*invdt;
    vel[1] = (pos[1]-ref[1])*invdt;
    vel[2] = (pos[2]-ref[2])*invdt;
  }

  return 0;
}



static int settlev(const Vector *pos, BigReal ma, BigReal mb, Vector *vel,
				   BigReal dt, Tensor *virial) {
  
  Vector rAB = pos[1]-pos[0];
  Vector rBC = pos[2]-pos[1];
  Vector rCA = pos[0]-pos[2];
 
  Vector AB = rAB.unit();
  Vector BC = rBC.unit();
  Vector CA = rCA.unit();
  
  BigReal cosA = -AB * CA;
  BigReal cosB = -BC * AB;
  BigReal cosC = -CA * BC;

  BigReal vab = (vel[1]-vel[0])*AB;
  BigReal vbc = (vel[2]-vel[1])*BC;
  BigReal vca = (vel[0]-vel[2])*CA;

  BigReal mab = ma+mb;
  
  BigReal d = (2*mab*mab + 2*ma*mb*cosA*cosB*cosC - 2*mb*mb*cosA*cosA
               - ma*mab*(cosB*cosB + cosC*cosC))*0.5/mb;

  BigReal tab = (vab*(2*mab - ma*cosC*cosC) +
                vbc*(mb*cosC*cosA - mab*cosB) +
                vca*(ma*cosB*cosC - 2*mb*cosA))*ma/d;
            
  BigReal tbc = (vbc*(mab*mab - mb*mb*cosA*cosA) +
                vca*ma*(mb*cosA*cosB - mab*cosC) +
                vab*ma*(mb*cosC*cosA - mab*cosB))/d;
  
  BigReal tca = (vca*(2*mab - ma*cosB*cosB) +
                vab*(ma*cosB*cosC - 2*mb*cosA) +
                vbc*(mb*cosA*cosB - mab*cosC))*ma/d;
 
  Vector ga = tab*AB - tca*CA;
  Vector gb = tbc*BC - tab*AB;
  Vector gc = tca*CA - tbc*BC;
#if 0
  if (virial) {
    *virial += 0.5*outer(tab, rAB)/dt;
    *virial += 0.5*outer(tbc, rBC)/dt;
    *virial += 0.5*outer(tca, rCA)/dt;
  }
#endif
  vel[0] += (0.5/ma)*ga;
  vel[1] += (0.5/mb)*gb;
  vel[2] += (0.5/mb)*gc;

  return 0;
}


int settle2(BigReal mO, BigReal mH, const Vector *pos,
                  Vector *vel, BigReal dt, Tensor *virial) {

  settlev(pos, mO, mH, vel, dt, virial);
  return 0;
}

