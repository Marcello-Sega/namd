/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MEMUSAGE_H
#define MEMUSAGE_H

long memusage();

class memusageinit {
public:
  memusageinit();
private:
  static int initialized;
  static long sbrkval;
  static long memusage_sbrk();
  friend long memusage();
};

static memusageinit memusageinitobject;

#endif

