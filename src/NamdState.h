//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

#ifndef _NAMDSTATE_H
#define _NAMDSTATE_H

#include "Controller.h"

class Molecule;
class SimParameters;
class Parameters;
class ConfigList;
class PDB;

// Everything needed to specify a simulation is in this object
// For the moment it is really only a structure.  Eventually
// I hope it encapsulates a namd state.  
class NamdState {
  friend class Namd; 
  private:
    Molecule *molecule;
    Parameters *parameters;
    SimParameters *simParameters;
    ConfigList *configList;
    PDB *pdb;
    Controller *controller;
  public:
    NamdState(void);
    ~NamdState() {}
    int configFileInit(char *);
    int status();
    void useController(Controller *controllerPtr) {controller=controllerPtr;}
    void runController(int numberOfCycles = 0)
		{ controller->run(numberOfCycles); }
};

#endif /* _NAMDSTATE_H */
 
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdState.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/02/11 18:51:49 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdState.h,v $
 * Revision 1.1001  1997/02/11 18:51:49  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1000  1997/02/06 15:58:48  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:57  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:27  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:32  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.2  1997/01/09 20:48:10  jim
 * added Controller code
 *
 * Revision 1.1  1996/08/16 04:39:46  ari
 * Initial revision
 *
 ***************************************************************************/

