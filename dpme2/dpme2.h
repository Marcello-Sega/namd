/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */

#ifndef DPME2_H
#define DPME2_H

#include <stdlib.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/times.h>
#include "stdio.h"
#include "math.h"
#include <sys/resource.h>
#include <unistd.h>
#include "malloc.h"

/* below is only for the HP cluster at UIUC */
#if 0
#include <sys/types.h>
#define _INCLUDE_POSIX_SOURCE
#define _INCLUDE_HPUX_SOURCE
#endif

#include "dpme2_pvm.h"
#include "dpme2def.h"
#include "protype.h"
#include "pvm3.h"


/* ONLY MACROS ARE DEFINED HERE */

#define Nabs(dumx) ((dumx) >= 0 ? (dumx) : -(dumx))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)
#define abs(x) ((x) >= 0 ? (x) : -(x))

/* Nint is eqvlnt to rint i.e. round x to the nearest integer */
/*#define Nint(dumx)  ( ((dumx) >=0.0) ? (int)((dumx) + 0.5) : (int)((dumx) - 0.5) )*/
/* use  this if rint is defined on your machine */
#define  Nint(dumx) rint(dumx)

/* below is to be used in charge_grid.c */
#define Nsign(i,j)  ((j) >=0.0 ?  Nabs((i)): -Nabs((i)) )

#endif
