/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "memusage.h"
#include "converse.h"
#ifndef WIN32
#include <stdio.h>
#include <sys/types.h>
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

#ifdef _NO_MALLOC_H
#define NO_MALLINFO
#endif

#ifdef MEMUSAGE_USE_SBRK
#define NO_MALLINFO
#define NO_PS
#endif


#ifdef MEMUSAGE_USE_MSTATS

#include <malloc/malloc.h>

long memusage_mstats() {
  struct mstats ms = mstats();
  long memtotal = ms.bytes_used;
  return memtotal;
}

#else
inline long memusage_mstats() { return 0; }
#endif


#ifndef NO_MALLINFO

#include <malloc.h>

long memusage_mallinfo() {
  struct mallinfo mi = mallinfo();
  long memtotal = mi.usmblks + mi.uordblks + mi.hblkhd;
  return memtotal;
}

#else
inline long memusage_mallinfo() { return 0; }
#endif


#ifndef NO_PS

inline long memusage_ps() {
  char pscmd[100];
  sprintf(pscmd, "/bin/ps -o vsz= -p %d", getpid());
  long vsz = 0;
  FILE *p = popen(pscmd, "r");
  if ( p ) {
    fscanf(p, "%ld", &vsz);
    pclose(p);
  }
  return ( vsz * 1024 );
}

#else
inline long memusage_ps() { return 0; }
#endif


long memusage() {

  long memtotal = 0;

#if CHARM_VERSION > 50911
  if (CmiMemoryIs(CMI_MEMORY_IS_GNU) ) memtotal = CmiMemoryUsage();
#endif

  if ( ! memtotal ) memtotal = memusage_mstats();

  if ( ! memtotal ) memtotal = memusage_mallinfo();

  if ( ! memtotal ) memtotal = memusageinit::memusage_sbrk();

  if ( ! memtotal ) memtotal = memusage_ps();

  return memtotal;

}

