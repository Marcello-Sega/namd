//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef _NAMDSTATE_H
#define _NAMDSTATE_H

#include "Lattice.h"

class Molecule;
class SimParameters;
class Parameters;
class ConfigList;
class PDB;
class Controller;
class SMDData;

// Everything needed to specify a simulation is in this object
// For the moment it is really only a structure.  Eventually
// I hope it encapsulates a namd state.  
class NamdState {
  friend class Namd; 
  friend class Node;
  friend class Controller;
  friend class TestController;
  private:
    Molecule *molecule;
    Parameters *parameters;
    SimParameters *simParameters;
    ConfigList *configList;
    PDB *pdb;
    Controller *controller;
    Lattice lattice;
    SMDData *smdData;
  public:
    NamdState(void);
    ~NamdState() {}
    int configFileInit(char *);
    int status();
    void useController(Controller *controllerPtr);
    void runController(void);
};

#endif /* _NAMDSTATE_H */
 
