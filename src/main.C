/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"

#include "main.decl.h"
#include "main.h"

// Needed for namd.1.X components
#include "Namd.h"
#include "Communicate.h"
#include "Inform.h"

Inform namdErr("ERROR");
Inform namdWarn("Warning");
Inform namdInfo("Info");
Inform namdDebug("** DEBUG **");

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

class main : public Chare
{
public:
  main(CkArgMsg *msg)
  {
    int argc = msg->argc;
    char **argv = msg->argv;

    // print banner
    iout << iINFO << "NAMD 2.1b4"
#ifdef NAMD_FFTW
         << " (includes FFTW, do not distribute)"
#endif
         << "\n"
#if 0
         << iWARN << "\n"
         << iWARN << "          ***  UNRELEASED EXPERIMENTAL VERSION  ***\n"
         << iWARN << "\n"
#endif
         << iINFO << "Please complete the registration form at\n"
         << iINFO << "http://www.ks.uiuc.edu/Research/namd/download.html\n"
         << iINFO << "and send feedback or bug reports to namd@ks.uiuc.edu\n"
         << endi;

    // Namd object is only on Pe(0)
    Namd *namd = new Namd;

    if (argc >= 2) {
	namd->startup(argv[argc-1]);
    }
    else {
       NAMD_die("No simulation config file specified on command line.");
    }
    DebugM(1, "main() - leaving - Charm should queue up messages now!\n");
  }
};

#include "main.def.h"

