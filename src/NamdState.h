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
  friend Namd; 
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
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1997/01/09 20:48:10 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdState.h,v $
 * Revision 1.2  1997/01/09 20:48:10  jim
 * added Controller code
 *
 * Revision 1.1  1996/08/16 04:39:46  ari
 * Initial revision
 *
 ***************************************************************************/

