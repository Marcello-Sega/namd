/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

#ifndef _NAMDSTATE_H
#define _NAMDSTATE_H

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
  public:
    NamdState(void);
    ~NamdState() {}
    int configFileInit(char *);
    int status();
};

#endif /* _NAMDSTATE_H */
 
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdState.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/16 04:39:46 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdState.h,v $
 * Revision 1.1  1996/08/16 04:39:46  ari
 * Initial revision
 *
 ***************************************************************************/

