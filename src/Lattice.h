//-*-c++-*-
#ifndef LATTICE_H
#define LATTICE_H

#include "NamdTypes.h"

class Lattice
{
public:
  Lattice(void) {};

  static int index(int i=0, int j=0, int k=0)
  {
    return 9 * (k+1) + 3 * (j+1) + (i+1);
  }

  void set(Vector aa, Vector bb, Vector cc)
  {
    a = aa; b = bb; c = cc;
  }

  Position* create(Position *d, int n, int i)
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

  void destroy(Position **d, int i)
  {
    if ( i != 13 ) delete [] *d;
    *d = NULL;
  }

private:
  Vector a,b,c;

};

#endif
