/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "converse.h"
#include "common.h"
#include "Namd.h"
#include "BackEnd.h"


int main(int argc, char **argv) {
  BackEnd::init(argc,argv);
  if ( argc >= 2 ) {
    Namd *namd = new Namd;
    namd->startup(argv[argc-1]);
    CsdScheduler(-1);  // process messages
  } else {
    NAMD_die("No simulation config file specified on command line.");
  }
  BackEnd::exit();
}

