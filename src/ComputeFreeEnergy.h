/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#ifndef COMPUTEFREEENERGY_H
#define COMPUTEFREEENERGY_H

#include "ComputeGlobalMaster.h"
#include "NamdTypes.h"
#include "iostream.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;
class Molecule;
class SimParameters;

class ComputeFreeEnergy : public ComputeGlobalMaster {
private:
  friend class ComputeGlobal;
  ComputeFreeEnergy(ComputeGlobal *);
  ~ComputeFreeEnergy();
  virtual void initialize();
  virtual void calculate();
  virtual void user_initialize();
  virtual void user_calculate();
  ComputeGlobalConfigMsg *configMsg;
  ComputeGlobalResultsMsg *resultsMsg;
  Molecule *molecule;
  SimParameters *simParams;
  strstream config;
protected:
  // These all return -1 on error.
  int getAtomID(const char *segid, int resid, const char *aname);
  int requestAtom(int atomid);
  int getPosition(int atomid, Position &position);
  int addForce(int atomid, Force force);
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeFreeEnergy.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1998/02/13 22:02:40 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeFreeEnergy.h,v $
 * Revision 1.4  1998/02/13 22:02:40  jim
 * Added script reading from config file and used streams in free energy.
 *
 * Revision 1.3  1998/02/11 17:49:03  jim
 * Added filename parameter to user_initialize().
 *
 * Revision 1.2  1998/02/11 07:31:35  jim
 * Finished interface to free energy perturbation code, including method
 * for determining atomid from segnamde, resid, and atomname.
 *
 * Revision 1.1  1998/02/10 06:45:09  jim
 * Added class ComputeFreeEnergy.
 *
 *
 *
 *
 ***************************************************************************/

