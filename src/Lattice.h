#ifndef LATTICE_H
#define LATTICE_H

import "NamdTypes.h"

class Lattice
{
public:
  Lattice(void) : dims(0) {};

  Position* create(Position *d, int n, int i)
  {
    Position *dt;
    if ( i != 13 )
    {
      dt = new Position[n];
      Vector shift = (i/9-1) * a + ((i/3)%3-1) * b + (i%3-1) * c;
      for( int j = 0; j < n; ++j )
        dt[j] = d[j] + shift;
    }
    else
    {
      dt = d;
    }
    return dt;
  }

  void create(Position **d, int n, int i)
  {
    if ( i != 13 ) delete [] *d;
    *d = NULL;
  }

private:
  int dims;
  Vector a,b,c;

};


class LatticeTransform
{
public:
  LatticeTransform(int i=0, int j=0, int k=0)
  {
    index = 9 * i + 3 * j + k;
  }
  int index(void) { return index; }

private:
  int index;
};

#endif
