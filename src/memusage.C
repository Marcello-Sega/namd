/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "memusage.h"
#ifndef WIN32
#include <unistd.h>
#else
int sbrk(int) { return 0; }
#endif

int memusageinit::initialized;
long memusageinit::sbrkval;

memusageinit::memusageinit() {
  if ( initialized == 0 ) {
    sbrkval = (long) sbrk(0);
    initialized = 1;
  }
}

long memusage() {
  long newval = (long) sbrk(0);
  return ( newval - memusageinit::sbrkval );
}

