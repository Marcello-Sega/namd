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

long memusageinit::memusage_sbrk() {
  long newval = (long) sbrk(0);
  return ( newval - memusageinit::sbrkval );
}

#ifdef WIN32
#define _NO_MALLOC_H
#endif

#ifndef _NO_MALLOC_H

#include <malloc.h>

long memusage_mallinfo() {

  struct mallinfo mi = mallinfo();

  long memtotal = mi.usmblks + mi.uordblks + mi.hblkhd;

  return memtotal;

}

#else // ifndef _NO_MALLOC_H

long memusage_mallinfo() { return 0; }

#endif

long memusage() {

  long memtotal = memusage_mallinfo();

  if ( ! memtotal ) memtotal = memusageinit::memusage_sbrk();

  return memtotal;

}

