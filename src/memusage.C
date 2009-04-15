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
unsigned long memusageinit::sbrkval;

memusageinit::memusageinit() {
  if ( initialized == 0 ) {
    sbrkval = (unsigned long) sbrk(0);
    initialized = 1;
  }
}

unsigned long memusageinit::memusage_sbrk() {
  unsigned long newval = (unsigned long) sbrk(0);
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

unsigned long memusage_mstats() {
  struct mstats ms = mstats();
  unsigned long memtotal = ms.bytes_used;
  return memtotal;
}

#else
inline unsigned long memusage_mstats() { return 0; }
#endif


#ifndef NO_MALLINFO

#include <malloc.h>

unsigned long memusage_mallinfo() {
  struct mallinfo mi = mallinfo();
  // unsigned long memtotal = mi.usmblks + mi.uordblks + mi.hblkhd;
  unsigned long memtotal = (unsigned int) mi.uordblks;
  unsigned long memtotal2 = (unsigned int) mi.usmblks;
  memtotal2 += (unsigned int) mi.hblkhd;
  if ( memtotal2 > memtotal ) memtotal = memtotal2;

  // printf("mallinfo %d %d %d\n", mi.usmblks, mi.uordblks, mi.hblkhd);

  return memtotal;
}

#else
inline unsigned long memusage_mallinfo() { return 0; }
#endif


#ifndef NO_PS

inline unsigned long memusage_ps() {
  char pscmd[100];
  sprintf(pscmd, "/bin/ps -o vsz= -p %d", getpid());
  unsigned long vsz = 0;
  FILE *p = popen(pscmd, "r");
  if ( p ) {
    fscanf(p, "%ld", &vsz);
    pclose(p);
  }
  return ( vsz * (unsigned long) 1024 );
}

#else
inline unsigned long memusage_ps() { return 0; }
#endif


inline unsigned long memusage_proc_self_stat() {

  static int failed_once = 0;
  if ( failed_once ) return 0;  // no point in retrying

  FILE *f = fopen("/proc/self/stat","r");
  if ( ! f ) { failed_once = 1; return 0; }
  for ( int i=0; i<22; ++i ) fscanf(f,"%*s");
  unsigned long vsz = 0;  // should remain 0 on failure
  fscanf(f,"%lu",&vsz);
  fclose(f);
  if ( ! vsz ) failed_once = 1;
  // printf("/proc/self/stat reports %d MB\n", vsz/(1024*1024));
  return vsz;

}


unsigned long memusage(const char **source) {

  unsigned long memtotal = 0;
  const char* s = "ERROR";

  if (CmiMemoryIs(CMI_MEMORY_IS_GNU) ) {
    memtotal = CmiMemoryUsage();  s = "CmiMemoryUsage";
  }

  if ( ! memtotal ) {
    memtotal = memusage_proc_self_stat();  s = "/proc/self/stat";
  }

  if ( ! memtotal ) { memtotal = memusage_mstats(); s = "mstats"; }

  if ( ! memtotal ) { memtotal = memusage_mallinfo(); s = "mallinfo"; }

  if ( ! memtotal ) { memtotal = memusageinit::memusage_sbrk(); s = "sbrk"; }

  if ( ! memtotal ) { memtotal = memusage_ps(); s = "ps"; }

  if ( ! memtotal ) s = "nothing";

  if ( source ) *source = s;

  return memtotal;

}

