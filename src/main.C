/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"
#include "InfoStream.h"

#include "main.decl.h"
#include "main.h"

class main : public Chare
{
public:
  main(CkArgMsg *)
  {
    // print banner
    iout << iINFO << "NAMD 2.1"
#ifdef NAMD_FFTW
         << " (includes FFTW, do not distribute)"
#endif
         << "\n"
#if 1
         << iWARN << "\n"
         << iWARN << "          ***  UNRELEASED EXPERIMENTAL VERSION  ***\n"
         << iWARN << "\n"
#endif
         << iINFO << "Please complete the registration form at\n"
         << iINFO << "http://www.ks.uiuc.edu/Research/namd/download.html\n"
         << iINFO << "and send feedback or bug reports to namd@ks.uiuc.edu\n"
         << endi;
  }
};

#include "main.def.h"

