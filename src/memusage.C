/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "memusage.h"
#include "converse.h"
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
#define MEMUSAGE_USE_SBRK
#endif

#if defined(_NO_MALLOC_H) && !defined(MAC_OSX)
#define MEMUSAGE_USE_SBRK
#endif

#ifndef MEMUSAGE_USE_SBRK

#ifdef MAC_OSX
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

long memusage_mallinfo() {

#ifndef MAC_OSX 
  struct mallinfo mi = mallinfo();

  long memtotal = mi.usmblks + mi.uordblks + mi.hblkhd;
#else
  struct mstats ms = mstats();

  long memtotal = ms.bytes_used;
#endif

  return memtotal;

}

#else // ifndef MEMUSAGE_USE_SBRK

long memusage_mallinfo() { return 0; }

#endif

long memusage() {

  long memtotal = 0;

#if CHARM_VERSION > 50911
  if (CmiMemoryIs(CMI_MEMORY_IS_GNU) ) memtotal = CmiMemoryUsage();
#endif

  if ( ! memtotal ) memtotal = memusage_mallinfo();

  if ( ! memtotal ) memtotal = memusageinit::memusage_sbrk();


  return memtotal;

}

