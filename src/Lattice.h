//-*-c++-*-
#ifndef LATTICE_H
#define LATTICE_H

#include "NamdTypes.h"
#include <math.h>

class Lattice
{
public:
  Lattice(void) : anz(0), bnz(0), cnz(0) {};

  static int index(int i=0, int j=0, int k=0)
  {
    return 9 * (k+1) + 3 * (j+1) + (i+1);
  }

  void set(Vector aa, Vector bb, Vector cc)
  {
    a = aa; b = bb; c = cc;
    anz = ( a.length2() > 0 );
    bnz = ( b.length2() > 0 );
    cnz = ( c.length2() > 0 );
  }

  Vector nearest(Vector data, Vector ref) const
  {
    Vector result = data;
    if ( anz )
    {
      result.x = ref.x + drem(data.x-ref.x,a.x);
    }
    if ( bnz )
    {
      result.y = ref.y + drem(data.y-ref.y,b.y);
    }
    if ( cnz )
    {
      result.z = ref.z + drem(data.z-ref.z,c.z);
    }
    return result;
  }

  Vector delta(Vector v1, Vector v2) const
  {
    Vector result = v1 - v2;
    if ( anz )
    {
      result.x = drem(result.x,a.x);
    }
    if ( bnz )
    {
      result.y = drem(result.y,b.y);
    }
    if ( cnz )
    {
      result.z = drem(result.z,c.z);
    }
    return result;
  }

  Position* create(Position *d, int n, int i) const
  {
    Position *dt;
    if ( i != 13 )
    {
      dt = new Position[n];
      Vector shift = (i/9-1) * c + ((i/3)%3-1) * b + (i%3-1) * a;
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

  Vector dimension()
  {
  Vector rc;
  rc.x = a.x;
  rc.y = b.y;
  rc.z = c.z;
  return (rc);
  }

  BigReal volume(void) const
  {
    return ( a.x * b.y * c.z );
  }

  void rescale(BigReal factor)
  {
    a *= factor;
    b *= factor;
    c *= factor;
  }

private:
  Vector a,b,c;
  int anz, bnz, cnz;

};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1997/03/21 23:05:36 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Lattice.h,v $
 * Revision 1.1003  1997/03/21 23:05:36  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1002  1997/03/19 11:54:24  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
