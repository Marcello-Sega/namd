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
#include "ckdefs.h"	// for CPrintf

class infostream : public ostrstream
  {
  private:
  char Buffer[1024];

  public:
  infostream() : ostrstream(Buffer,sizeof(Buffer)) {;}

  /* output using CPrintf() (end by inform) */
  void endi()
    {
    *this << ends;
    CPrintf("%s",Buffer);
    this->seekp(0);
    return;
    }

  /* output to stdout (end by console) */
  void endc()
    {
    *this << ends;
    cout << Buffer;
    this->seekp(0);
    return;
    }

  };

extern infostream iout;

/** now we have a stream.  Let's tell the stream when/where to output **/
/** iout << endi;  (send output to information console -- CPrintf) **/
/** iout << endc;  (send output to host console -- cout) **/
char * endi(infostream& s) { s.endi(); return ""; }
char * endc(infostream& s) { s.endc(); return ""; }
#define endi endi(iout)
#define endc endc(iout)

/** define some basic messages **/
#define iWARN   "Warning"
#define iERROR  "ERROR"
#define iInfo   "Info"

#endif
