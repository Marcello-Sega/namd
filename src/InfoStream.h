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

#include <stdio.h>
#include <iostream.h>	// for cout
#include <strstream.h>	// for ostrstream
#include "ckdefs.h"	// for CPrintf

class infostream : public ostrstream
  {
  private:
  char Buffer[1024];

  public:
  infostream() : ostrstream(Buffer,sizeof(Buffer)) {;}

  /* output using CPrintf() (end by inform) */
  infostream& endi()
    {
    *this << ends;
    CPrintf("%s",Buffer);
    this->seekp(0);	// clear buffer
    return *this;
    }

  /* output to stdout (end by console) */
  infostream& endc()
    {
    *this << ends;
    cout << Buffer;
    this->seekp(0);	// clear buffer
    return *this;
    }

  };

extern infostream iout;

/** now we have a stream.  Let's tell the stream when/where to output **/
/** iout << endi;  (send output to information console -- CPrintf) **/
/** iout << endc;  (send output to host console -- cout) **/
inline infostream& endi(infostream& s)   { return s.endi(); }
inline infostream& endc(infostream& s)   { return s.endc(); }
/** common messages **/
inline infostream& iINFO(infostream& s)  { s << "Info: "; return s; }
inline infostream& iWARN(infostream& s)  { s << "Warning: "; return s; }
inline infostream& iERROR(infostream& s) { s << "ERROR: "; return s; }

#endif

