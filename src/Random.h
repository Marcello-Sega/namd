/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
 * Copyright (c) 1993 Martin Birgmeier
 * All rights reserved.
 *
 * You may redistribute unmodified or modified versions of this source
 * code provided that the above copyright notice and this and the
 * following conditions are retained.
 *
 * This software is provided ``as is'', and comes with no warranties
 * of any kind. I shall in no event be liable for anything that happens
 * to anyone/anything when using this software.
 */

#ifndef RANDOM_H
#define RANDOM_H

#include <math.h>
#include "common.h"

#define	RAND48_SEED_0	(0x330e)
#define	RAND48_SEED_1	(0xabcd)
#define	RAND48_SEED_2	(0x1234)
#define	RAND48_MULT_0	(0xe66d)
#define	RAND48_MULT_1	(0xdeec)
#define	RAND48_MULT_2	(0x0005)
#define	RAND48_ADD_0	(0x000b)
#define	RAND48_ADD_1	(0x0000)
#define	RAND48_ADD_2	(0x0000)

class Random {

private:

  double second_gaussian;
  unsigned short second_gaussian_waiting;
  unsigned short rand48_seed_0;
  unsigned short rand48_seed_1;
  unsigned short rand48_seed_2;
  unsigned short rand48_mult_0;
  unsigned short rand48_mult_1;
  unsigned short rand48_mult_2;
  unsigned short rand48_add_0;
  unsigned short rand48_add_1;
  unsigned short rand48_add_2;

public:

  // default constructor
  Random(void) {
    second_gaussian = 0;
    second_gaussian_waiting = 0;
    rand48_seed_0 = RAND48_SEED_0;
    rand48_seed_1 = RAND48_SEED_1;
    rand48_seed_2 = RAND48_SEED_2;
    rand48_mult_0 = RAND48_MULT_0;
    rand48_mult_1 = RAND48_MULT_1;
    rand48_mult_2 = RAND48_MULT_2;
    rand48_add_0  = RAND48_ADD_0;
    rand48_add_1  = RAND48_ADD_1;
    rand48_add_2  = RAND48_ADD_2;
  }

  // constructor with seed
  Random(long seed) {
    second_gaussian = 0;
    second_gaussian_waiting = 0;
    rand48_seed_0 = RAND48_SEED_0;
    rand48_seed_1 = (unsigned short) seed;
    rand48_seed_2 = (unsigned short) (seed >> 16);
    rand48_mult_0 = RAND48_MULT_0;
    rand48_mult_1 = RAND48_MULT_1;
    rand48_mult_2 = RAND48_MULT_2;
    rand48_add_0  = RAND48_ADD_0;
    rand48_add_1  = RAND48_ADD_1;
    rand48_add_2  = RAND48_ADD_2;
  }

  // reinitialize with seed
  void init(long seed) {
    second_gaussian = 0;
    second_gaussian_waiting = 0;
    rand48_seed_0 = RAND48_SEED_0;
    rand48_seed_1 = (unsigned short) seed;
    rand48_seed_2 = (unsigned short) (seed >> 16);
    rand48_mult_0 = RAND48_MULT_0;
    rand48_mult_1 = RAND48_MULT_1;
    rand48_mult_2 = RAND48_MULT_2;
    rand48_add_0  = RAND48_ADD_0;
    rand48_add_1  = RAND48_ADD_1;
    rand48_add_2  = RAND48_ADD_2;
  }

  // advance generator by one
  void skip(void) {
    unsigned long accu;
    unsigned short temp_0;
    unsigned short temp_1;

    accu = (unsigned long) rand48_mult_0 * (unsigned long) rand48_seed_0 +
           (unsigned long) rand48_add_0;
    temp_0 = (unsigned short) accu;	/* lower 16 bits */
    accu >>= sizeof(unsigned short) * 8;
    accu += (unsigned long) rand48_mult_0 * (unsigned long) rand48_seed_1 +
            (unsigned long) rand48_mult_1 * (unsigned long) rand48_seed_0 +
            (unsigned long) rand48_add_1;
    temp_1 = (unsigned short) accu;	/* middle 16 bits */
    accu >>= sizeof(unsigned short) * 8;
    accu += rand48_mult_0 * rand48_seed_2 + rand48_mult_1 * rand48_seed_1 +
            rand48_mult_2 * rand48_seed_0 + rand48_add_2;
    rand48_seed_0 = temp_0;
    rand48_seed_1 = temp_1;
    rand48_seed_2 = (unsigned short) accu;
  }

  // split into numStreams different steams and take stream iStream
  void split(int iStream, int numStreams) {

    int i;

    // make sure that numStreams is odd to ensure maximum period
    numStreams |= 1;

    // iterate to get to the correct stream
    for ( i = 0; i < iStream; ++i ) skip();

    // save seed and add so we can use skip() for our calculations
    unsigned short save_seed_0 = rand48_seed_0;
    unsigned short save_seed_1 = rand48_seed_1;
    unsigned short save_seed_2 = rand48_seed_2;

    // calculate c *= ( 1 + a + ... + a^(numStreams-1) )
    rand48_seed_0 = rand48_add_0;
    rand48_seed_1 = rand48_add_1;
    rand48_seed_2 = rand48_add_2;
    for ( i = 1; i < numStreams; ++i ) skip();
    unsigned short new_add_0 = rand48_seed_0;
    unsigned short new_add_1 = rand48_seed_1;
    unsigned short new_add_2 = rand48_seed_2;

    // calculate a = a^numStreams
    rand48_seed_0 = rand48_mult_0;
    rand48_seed_1 = rand48_mult_1;
    rand48_seed_2 = rand48_mult_2;
    rand48_add_0  = 0; rand48_add_1  = 0; rand48_add_2  = 0;
    for ( i = 1; i < numStreams; ++i ) skip();
    rand48_mult_0 = rand48_seed_0;
    rand48_mult_1 = rand48_seed_1;
    rand48_mult_2 = rand48_seed_2;

    rand48_add_0  = new_add_0;
    rand48_add_1  = new_add_1;
    rand48_add_2  = new_add_2;
    rand48_seed_0 = save_seed_0;
    rand48_seed_1 = save_seed_1;
    rand48_seed_2 = save_seed_2;

    second_gaussian = 0;
    second_gaussian_waiting = 0;
  }

  // return a number uniformly distributed between 0 and 1
  BigReal uniform(void) {
    skip();
    return ldexp((double) rand48_seed_0, -48) +
           ldexp((double) rand48_seed_1, -32) +
           ldexp((double) rand48_seed_2, -16);
  }

  // return a number from a standard gaussian distribution
  BigReal gaussian(void) {
    BigReal fac, r, v1, v2;

    if (second_gaussian_waiting) {
      second_gaussian_waiting = 0;
      return second_gaussian;
    } else {
      r = 2.;                 // r >= 1.523e-8 ensures abs result < 6
      while (r >=1. || r < 1.523e-8) { // make sure we are within unit circle
        v1 = 2.0 * uniform() - 1.0;
        v2 = 2.0 * uniform() - 1.0;
        r = v1*v1 + v2*v2;
      }
      fac = sqrt(-2.0 * log(r)/r);
      // now make the Box-Muller transformation to get two normally
      // distributed random numbers. Save one and return the other.
      second_gaussian_waiting = 1;
      second_gaussian = v1 * fac;
      return v2 * fac;
    }
  }

  // return a vector of gaussian random numbers
  Vector gaussian_vector(void) {
    return Vector( gaussian(), gaussian(), gaussian() );
  }

  // return a random long
  long integer(void) {
    skip();
    return ((long) rand48_seed_2 << 15) + ((long) rand48_seed_1 >> 1);
  }

};

#endif  // RANDOM_H

