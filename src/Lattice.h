//-*-c++-*-
#ifndef LATTICE_H
#define LATTICE_H

#include "NamdTypes.h"
#include <math.h>

#define rint(X) floor((X)+0.5)

typedef Vector ScaledPosition;

class Lattice
{
public:
  Lattice(void) : a1(0), a2(0), a3(0) {};

  // maps a transformation triplet onto a single integer
  static int index(int i=0, int j=0, int k=0)
  {
    return 9 * (k+1) + 3 * (j+1) + (i+1);
  }

  // sets lattice basis vectors and origin (fixed center)
  // sets lattice basis vectors and origin (fixed center)
  void set(Vector A, Vector B, Vector C, Position Origin)
  {
    a1 = A.x;  b1 = ( a1 ? 1. / a1 : 0 );
    a2 = B.y;  b2 = ( a2 ? 1. / a2 : 0 );
    a3 = C.z;  b3 = ( a3 ? 1. / a3 : 0 );
    o = Origin;
  }

  // rescale lattice dimensions by factor, origin doesn't move
  void rescale(Vector factor)
  {
    a1 *= factor.x;  b1 = ( a1 ? 1. / a1 : 0 );
    a2 *= factor.y;  b2 = ( a2 ? 1. / a2 : 0 );
    a3 *= factor.z;  b3 = ( a3 ? 1. / a3 : 0 );
  }

  // rescale a position, keeping origin constant, assume 3D
  void rescale(Position &p, Vector factor) const
  {
    p -= o;
    p.x *= factor.x;
    p.y *= factor.y;
    p.z *= factor.z;
    p += o;
  }

  // transform scaled position to unscaled position
  Position unscale(ScaledPosition s) const
  {
    return Vector
    (
	( a1 ? ( o.x + a1 * s.x ) : s.x ),
	( a2 ? ( o.y + a2 * s.y ) : s.y ),
	( a3 ? ( o.z + a3 * s.z ) : s.z )
    );
  }

  // transform unscaled position to scaled position
  ScaledPosition scale(Position p) const
  {
    return Vector
    (
	( a1 ? ( b1 * ( p.x - o.x ) ) : p.x ),
	( a2 ? ( b2 * ( p.y - o.y ) ) : p.y ),
	( a3 ? ( b3 * ( p.z - o.z ) ) : p.z )
    );
  }

  // transforms a position nearest to a SCALED reference position
  Position nearest(Position data, ScaledPosition ref) const
  {
    BigReal tmp;
    return Vector
    (
	( a1 ? ( tmp=b1*(data.x-o.x)-ref.x, o.x+a1*(ref.x+tmp-rint(tmp)) ) : data.x ),
	( a2 ? ( tmp=b2*(data.y-o.y)-ref.y, o.y+a2*(ref.y+tmp-rint(tmp)) ) : data.y ),
	( a3 ? ( tmp=b3*(data.z-o.z)-ref.z, o.z+a3*(ref.z+tmp-rint(tmp)) ) : data.z )
    );
  }

  // transforms a position nearest to a SCALED reference position
  // adds transform for later reversal
  Position nearest(Position data, ScaledPosition ref, Transform *t) const
  {
    BigReal tmp;
    BigReal rit;
    Vector v = data;
    if ( a1 ) {
      tmp=b1*(data.x-o.x)-ref.x;
      rit = rint(tmp);
      v.x = o.x+a1*(ref.x+tmp-rit);
      t->i -= (int) rit;
    }
    if ( a2 ) {
      tmp=b2*(data.y-o.y)-ref.y;
      rit = rint(tmp);
      v.y = o.y+a2*(ref.y+tmp-rit);
      t->j -= (int) rit;
    }
    if ( a3 ) {
      tmp=b3*(data.z-o.z)-ref.z;
      rit = rint(tmp);
      v.z = o.z+a3*(ref.z+tmp-rit);
      t->k -= (int) rit;
    }
    return v;
  }

  // reverses cumulative transformations for output
  Position reverse_transform(Position data, const Transform &t) const
  {
    return ( data - Vector(t.i*a1, t.j*a2, t.k*a3) );
  }

  // calculates shortest vector from p2 to p1 (equivalent to p1 - p2)
  Vector delta(Position p1, Position p2) const
  {
    Vector result = p1 - p2;
    result.x -= a1 * rint( b1 * result.x );
    result.y -= a2 * rint( b2 * result.y );
    result.z -= a3 * rint( b3 * result.z );
    return result;
  }

  // calculates shortest vector from origin to p1 (equivalent to p1 - o)
  Vector delta(Position p1) const
  {
    Vector result = p1 - o;
    result.x -= a1 * rint( b1 * result.x );
    result.y -= a2 * rint( b2 * result.y );
    result.z -= a3 * rint( b3 * result.z );
    return result;
  }

  Position* create(Position *d, int n, int i) const
  {
    Position *dt;
    if ( i != 13 )
    {
      dt = new Position[n];
      Vector shift( (i%3-1) * a1 , ((i/3)%3-1) * a2 , (i/9-1) * a3 );
      for( int j = 0; j < n; ++j )
        dt[j] = d[j] + shift;
    }
    else
    {
      dt = d;
    }
    return dt;
  }

  void destroy(Position **d, int i) const
  {
    if ( i != 13 ) delete [] *d;
    *d = NULL;
  }

  BigReal a() const { return a1; }
  BigReal b() const { return a2; }
  BigReal c() const { return a3; }

  Vector dimension() const
  {
    return Vector(a1,a2,a3);
  }

  Vector origin() const
  {
    return o;
  }

  BigReal volume(void) const
  {
    return ( a1 * a2 * a3 );
  }

private:
  BigReal a1,a2,a3; // real lattice vectors (eventually)
  BigReal b1,b2,b3; // reciprocal lattice vectors (eventually)
  Vector o; // origin (fixed center of cell)

};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1010 $	$Date: 1999/01/06 19:19:20 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Lattice.h,v $
 * Revision 1.1010  1999/01/06 19:19:20  jim
 * Broadcast and Sequencers understand anisotropic volume rescaling factors.
 *
 * Revision 1.1009  1998/08/21 01:15:04  jim
 * Eliminated warnings.
 *
 * Revision 1.1008  1998/08/11 16:30:28  jim
 * Modified output from periodic boundary simulations to return atoms to
 * internally consistent coordinates.  We store the transformations which
 * were performed and undo them at the end.  It might be better to do this
 * by always keeping the original coordinates and only doing the transform
 * for the nonbonded terms but this works for now.
 *
 * Revision 1.1007  1998/04/06 16:34:08  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1006  1998/03/30 21:01:17  jim
 * Added nearest-image support for periodic boundary conditions to full direct.
 *
 * Revision 1.1005  1998/03/26 23:28:29  jim
 * Small changes for KCC port.  Altered use of strstream in ComputeFreeEnergy.
 *
 * Revision 1.1004  1997/03/27 08:04:18  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1003  1997/03/21 23:05:36  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1002  1997/03/19 11:54:24  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
