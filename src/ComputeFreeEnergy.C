/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#include "Namd.h"
#include "Node.h"
#include "ComputeFreeEnergy.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.top.h"
#include <stdio.h>
#ifdef GCC
#include <fstream.h>
#endif

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


void ComputeFreeEnergy::user_initialize()
{
  iout << iDEBUG << "Initializing free energy.\n"; 
  iout << iDEBUG << "***********************************\n"; 
  iout << config.str(); 
  iout << iDEBUG << "***********************************\n" << endi; 
}


void ComputeFreeEnergy::user_calculate()
{
  ;
}


int ComputeFreeEnergy::
	getAtomID(const char *segid, int resid, const char *aname)
{
  return molecule->get_atom_from_name(segid,resid,aname);
}

int ComputeFreeEnergy::requestAtom(int atomid)
{
  if ( ! configMsg ) return -1;  // failure
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1;  // failure
  configMsg->aid.add(atomid);
  return 0;  // success
}

int ComputeFreeEnergy::getPosition(int atomid, Position &position)
{
  AtomIDList::iterator a_i = aid.begin();
  AtomIDList::iterator a_e = aid.end();
  PositionList::iterator p_i = p.begin();
  for ( ; a_i != a_e; ++a_i, ++p_i ) {
    if ( *a_i == atomid ) {
      position = *p_i;
      return 0;  // success
    }
  }
  return -1;  // failure
}

int ComputeFreeEnergy::addForce(int atomid, Force force)
{
  if ( ! resultsMsg ) return -1;  // failure
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1;  // failure
  resultsMsg->aid.add(atomid);
  resultsMsg->f.add(force);
  return 0;  // success
}


ComputeFreeEnergy::ComputeFreeEnergy(ComputeGlobal *h) : ComputeGlobalMaster(h) {
  DebugM(3,"Constructing ComputeFreeEnergy\n");
  molecule = Node::Object()->molecule;
  simParams = Node::Object()->simParameters;
  configMsg = 0;  resultsMsg = 0;
}

ComputeFreeEnergy::~ComputeFreeEnergy() {
  DebugM(3,"Destructing ComputeFreeEnergy\n");
}


void ComputeFreeEnergy::initialize() {
  DebugM(4,"Initializing master\n");

  configMsg = new (MsgIndex(ComputeGlobalConfigMsg)) ComputeGlobalConfigMsg;

  // Get our script
  StringList *script = Node::Object()->configList->find("freeEnergyConfig");

  for ( ; script; script = script->next) {
    if ( script->data[0] == '{' ) {  // this is a flag, no } at end
      config << script->data + 1;    // so skip it at beginning
    } else {
      ifstream infile(script->data);
      if ( infile ) infile.get(*config.rdbuf(),'\0');
      if ( ! infile ) {
	char errmsg[100];
	sprintf(errmsg,"Error trying to read file %s!\n",script->data);
	NAMD_die(errmsg);
      }
    }
  }
  config.flush();

  iout << iDEBUG << "Free energy perturbation - initialize()\n" << endi; 
  user_initialize();

  // Send config to clients
  host->comm->sendComputeGlobalConfig(configMsg);
  configMsg = 0;
}


void ComputeFreeEnergy::calculate() {
  DebugM(4,"Calculating forces on master\n");

  resultsMsg = new (MsgIndex(ComputeGlobalResultsMsg)) ComputeGlobalResultsMsg;

  iout << iDEBUG << "Free energy perturbation - calculate()\n" << endi; 
  user_calculate();

  // Send results to clients
  DebugM(3,"Sending results (" << resultsMsg->aid.size() << " forces) on master\n");
  if ( resultsMsg->reconfig ) {
    DebugM(4,"Sending new configuration (" <<
			resultsMsg->newaid.size() << " atoms) on master\n");
  }
  host->comm->sendComputeGlobalResults(resultsMsg);
  resultsMsg = 0;
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1998/03/09 17:06:02 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeFreeEnergy.C,v $
 * Revision 1.6  1998/03/09 17:06:02  milind
 * Included fstream.h for GCC.
 *
 * Revision 1.5  1998/02/14 09:55:22  jim
 * Final changes to allow inline reading of { } delimited input.
 * Strings read this way begin with a { but do not end with a }.
 * This was done to allow inlines to be readily distinguishable.
 *
 * Revision 1.4  1998/02/13 22:02:39  jim
 * Added script reading from config file and used streams in free energy.
 *
 * Revision 1.3  1998/02/11 17:49:03  jim
 * Added filename parameter to user_initialize().
 *
 * Revision 1.2  1998/02/11 07:31:34  jim
 * Finished interface to free energy perturbation code, including method
 * for determining atomid from segnamde, resid, and atomname.
 *
 * Revision 1.1  1998/02/10 06:45:08  jim
 * Added class ComputeFreeEnergy.
 *
 *
 *
 *
 ***************************************************************************/
