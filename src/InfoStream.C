/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Defines a new stream: iout, for "i"nforming consoles.
*/

#include "InfoStream.h"
#include <charm++.h>
#include "Vector.h"
#include "Tensor.h"

/* output using CkPrintf() (end by inform) */
void infostream::endi() {
  *this << ends;
  CkPrintf("%s",iBuffer);
  (*this).seekp(0);   // clear buffer
}

ostream& iPE(ostream& s) {
  return s << "Pe(" << CkMyPe() << ')';
}

ostream& operator<<(ostream& strm, const Vector &v1) {
       strm << v1.x << " " << v1.y << " " << v1.z;
       return strm;
}

infostream& operator<<(infostream& strm, const Vector &v1) {
       strm << v1.x << " " << v1.y << " " << v1.z;
       return strm;
}


ostream& operator<<(ostream& strm, const Tensor &t1) {
       strm << t1.xx << " " << t1.xy << " " << t1.xz << " "
            << t1.yx << " " << t1.yy << " " << t1.yz << " "
            << t1.zx << " " << t1.zy << " " << t1.zz;
       return strm;
}

infostream& operator<<(infostream& strm, const Tensor &t1) {
       strm << t1.xx << " " << t1.xy << " " << t1.xz << " "
            << t1.yx << " " << t1.yy << " " << t1.yz << " "
            << t1.zx << " " << t1.zy << " " << t1.zz;
       return strm;
}


infostream iout;

