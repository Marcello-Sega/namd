/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "converse.h"
#include "ProcessorPrivate.h"
#include "common.h"
#include "Namd.h"

extern void _initCharm(int, char**);

// called on all procs by namd_init()
void slave_init(int argc, char **argv)
{
  ProcessorPrivateInit();
  _initCharm(argc, argv);
  CsdScheduler(-1);
}

// called on all procs by front end
void namd_init(int argc, char **argv) {
  ConverseInit(argc, argv, slave_init, 1, 1);
  if ( CmiMyPe() ) {
    slave_init(argc, argv);  // for procs that call main
    ConverseExit();  // should never return
  } else {
    ProcessorPrivateInit();
    _initCharm(argc, argv);  // message main Chare
  }
}

// called on proc 0 by front end
void namd_run(char *config) {
  Namd *namd = new Namd;
  namd->startup(config);
  CsdScheduler(-1);  // process messages
}

// called on proc 0 by front end
void namd_exit(void) {
  ConverseExit();
}

int main(int argc, char **argv) {
  namd_init(argc,argv);
  if ( argc >= 2 ) {
    namd_run(argv[argc-1]);
  } else {
    NAMD_die("No simulation config file specified on command line.");
  }
  namd_exit();
}

