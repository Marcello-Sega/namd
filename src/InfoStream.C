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
#include <stdio.h>

#ifndef NO_STRSTREAM_H
infostream::infostream() : ostrstream(iBuffer,sizeof(iBuffer)) {;}
#else
infostream::infostream() : ostringstream() {}
#endif

infostream::~infostream() {;}

/* output using CkPrintf() (end by inform) */
void infostream::endi() {
  *this << ends;
#ifndef NO_STRSTREAM_H
  CkPrintf("%s",iBuffer);
  fflush(stdout);  // since CkPrintf doesn't always flush
  (*this).seekp(0);   // clear buffer
#else
  CkPrintf("%s",str().c_str());
  fflush(stdout);  // since CkPrintf doesn't always flush
  (*this).seekp(0, ios_base::beg);
#endif
}

infostream& endi(infostream& s)  { s.endi(); return s; }

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

/* define how to use the remaining << args */
/** infostream<<ostream (hot to handle inherited modifiers) **/
infostream& infostream::operator<<(ostream& (*f)(ostream&)) { f(*this); return(*this); }
/** infostream<<infostream (how to handle class modifiers) **/
infostream& infostream::operator<<(infostream& (*f)(infostream&)) { return f(*this); }

#define LOCALMOD(type) infostream& infostream::operator<<(type x) \
		{ (ostream&)(*this) << x; return(*this); }
/** << characters **/
LOCALMOD(char)
LOCALMOD(unsigned char)
LOCALMOD(const char *)
/** << integers **/
LOCALMOD(int)
LOCALMOD(long)
LOCALMOD(short)
LOCALMOD(unsigned int)
LOCALMOD(unsigned long)
LOCALMOD(unsigned short)
/** << floats **/
LOCALMOD(float)
LOCALMOD(double)
/** << pointers **/
LOCALMOD(void *)
LOCALMOD(streambuf *)
#undef LOCALMOD

/** common messages **/
/** iINFO, iWARN, iERROR, iDEBUG provide initial headings. **/
/** iINFOF, iWARNF, iERRORF, iDEBUGF provide initial headings with file name
    and line numbers. **/
ostream& iINFO (ostream& s)  { return s << "Info: "; }
ostream& iWARN (ostream& s)  { return s << "Warning: "; }
ostream& iERROR(ostream& s)  { return s << "ERROR: "; }
ostream& iDEBUG(ostream& s)  { return s << "DEBUG: "; }

infostream iout;

