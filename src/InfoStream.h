//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *  Defines a new stream: iout, for "i"nforming consoles.
 ***************************************************************************/
#ifndef INFOSTREAM_H
#define INFOSTREAM_H

#include <iostream.h>	// for cout
#include <strstream.h>	// for ostrstream
#include <stdio.h>	// for CkPrintf
#include "charm++.h"	// for CkPrintf

class infostream : public ostrstream
{
  private:
  char iBuffer[16384];

  public:
  infostream() : ostrstream(iBuffer,sizeof(iBuffer)) {;}
  ~infostream() {;}

  /* output using CkPrintf() (end by inform) */
  void endi()
    {
    *this << ends;
    CkPrintf("%s",iBuffer);
    (*this).seekp(0);	// clear buffer
    }

  /* output to stdout (end by console) */
  void endc()
    {
    *this << ends;
    cout << iBuffer << ends;
    (*this).seekp(0);	// clear buffer
    }

  /* define how to use the remaining << args */
  /** infostream<<ostream (hot to handle inherited modifiers) **/
  infostream& operator<<(ostream& (*f)(ostream&)) { f(*this); return(*this); }
  /** infostream<<infostream (how to handle class modifiers) **/
  infostream& operator<<(infostream& (*f)(infostream&)) { return f(*this); }

  #define LOCALMOD(type) infostream& operator<<(type x) \
		{ (ostream&)(*this) << x; return(*this); }
  /** << characters **/
  LOCALMOD(char);
  LOCALMOD(unsigned char);
  LOCALMOD(const char *);
  /** << integers **/
  LOCALMOD(int);
  LOCALMOD(long);
  LOCALMOD(short);
  LOCALMOD(unsigned int);
  LOCALMOD(unsigned long);
  LOCALMOD(unsigned short);
  /** << floats **/
  LOCALMOD(float);
  LOCALMOD(double);
  /** << pointers **/
  LOCALMOD(void *);
  LOCALMOD(streambuf *);
  #undef LOCALMOD
};

/** modifiers **/
inline infostream& endc(infostream& s)  { s.endc(); return s; }
inline infostream& endi(infostream& s)  { s.endi(); return s; }

/** common messages **/
/** iINFO, iWARN, iERROR, iDEBUG provide initial headings. **/
/** iINFOF, iWARNF, iERRORF, iDEBUGF provide initial headings with file name
    and line numbers. **/
inline ostream& iINFO (ostream& s)  { return s << "Info: "; }
inline ostream& iWARN (ostream& s)  { return s << "Warning: "; }
inline ostream& iERROR(ostream& s)  { return s << "ERROR: "; }
inline ostream& iDEBUG(ostream& s)  { return s << "DEBUG: "; }
inline ostream& iPE(ostream& s)     { return s << "Pe(" << CkMyPe() << ')'; }

#define iFILE __FILE__<<'('<<__LINE__<<"): "
#define iINFOF  iINFO << iFILE
#define iWARNF  iWARN << iFILE
#define iERRORF  iERROR << iFILE
#define iDEBUGF  iDEBUG << iFILE


extern infostream iout;

#endif /* INFOSTREAM_H */

