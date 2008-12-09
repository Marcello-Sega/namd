/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2006 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: fortread.h,v $
 *      $Author: jim $       $Locker:  $             $State: Exp $
 *      $Revision: 1.1 $       $Date: 2008/12/09 19:46:22 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Unformatted and formatted fortran file reading routines for
 *   use in various plugins to simplify dealing with fortran i/o quirks.
 ***************************************************************************/

#ifndef FORTREAD_H
#define FORTREAD_H

#include <stdlib.h>
#include <stdio.h>

#include "endianswap.h"


/*  Unformatted reads.
 *
 *   Each function reads the next record from the file (provided it contains
 *   no more than n elements), optionally swapping its contents before
 *   writing it into dest. 
 *   Returns the number of elements on success, 0 on failure.
 *
 *   TODO: These should perhaps rewind the file to the beginning of the
 *   record on failure.
 */

/* Only works with aligned four-byte quantities, swap4_aligned() will */
/* cause a bus error on sume platforms if used on unaligned data      */
static int fortread_4(void *dest, int n, int swap, FILE *fd) {
  int dataBegin, dataEnd, count;

  if (fread(&dataBegin, sizeof(int), 1, fd) != 1) return 0;
  if (swap) swap4_aligned(&dataBegin, 1);
  if ((dataBegin <= 0) || (n < dataBegin/4)) return 0;

  count = fread(dest, 4, dataBegin/4, fd);
  if (count != dataBegin/4) return 0;
  if (swap) swap4_aligned(dest, dataBegin/4);

  if (fread(&dataEnd, sizeof(int), 1, fd) != 1) return 0;
  if (swap) swap4_aligned(&dataBegin, 1);
  if (dataEnd != dataBegin) return 0;

  return count;
}

/* Formatted reads:
 *
 * copy at most 'len' characters from source to target. 
 *
 * leading white space is skipped over but counts towards 'len'.
 * the copy stops at first whitspace or a '\0'.
 * unlike strncpy(3) the result will always be \0 terminated.
 *
 * intended for copying (short) strings from formatted fortran  
 * i/o files that must not contain whitespace (e.g. residue names, 
 * atom name/types etc. in .pdb, .psf and alike.). 
 */
static void strnwscpy(char *target, const char *source, const int len) {
  int i, c;

  for (i=0, c=0; i<len; ++i) {
    if (*source == '\0' || (c > 0 && *source == ' ')) {
      break;
    }

    if (*source == ' ') { 
      source++;
    } else {
      *target++ = *source++;
      c++;
    }
  }
  *target = '\0';
}

#endif

