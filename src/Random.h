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
  }

  // return a random double
  double next_double(void) {
    skip();
    return ldexp((double) rand48_seed_0, -48) +
           ldexp((double) rand48_seed_1, -32) +
           ldexp((double) rand48_seed_2, -16);
  }

  // return a random long
  long next_long(void) {
    skip();
    return ((long) rand48_seed_2 << 15) + ((long) rand48_seed_1 >> 1);
  }

};

#endif  // RANDOM_H

