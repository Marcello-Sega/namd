/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#include <string.h>
#include <strstream.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"
#include "Vector.h"
#include "FreeEnergyVector.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"
#include "ComputeFreeEnergy.h"
#include "FreeEnergyParse.h"

#include "Namd.h"
#include "Node.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.top.h"
#include <stdio.h>
#include <fstream.h>

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


void ComputeFreeEnergy::update() {
//-----------------------------------------------------------------
// get lambdas from LambdaManager, inform RestraintManager.
// calculate gradients for each center-of-mass of each restraint,
// and apply the forces to the atoms involved
//-----------------------------------------------------------------
  double  LambdaKf, LambdaRef;
  double  Sum_dU_dLambdas;

  if (m_LambdaManager.GetLambdas(LambdaKf, LambdaRef)) {

    // stuff that's done every time step
    m_RestraintManager.SetLambdas(LambdaKf, LambdaRef);
    m_RestraintManager.UpdateCOMs(*this);
    m_RestraintManager.AddForces(*this);
    if (m_LambdaManager.IsTimeToClearAccumulator()) {
      m_LambdaManager.ZeroAccumulator();
    }
    Sum_dU_dLambdas = m_RestraintManager.Sum_dU_dLambdas();
    m_LambdaManager.Accumulate(Sum_dU_dLambdas);

    // for integrating all the MCTI averages
    if (m_LambdaManager.IsEndOf_MCTI_Step()) {
      m_LambdaManager.Integrate_MCTI();
    }

    // stuff that's done when it's time to print
    if (m_LambdaManager.IsFirstStep()) {
      m_LambdaManager.PrintLambdaHeader(simParams->dt);
    }
    if (m_LambdaManager.IsTimeToPrint()) {
      m_LambdaManager.PrintHeader(simParams->dt);
      if (m_LambdaManager.IsTimeToPrint_dU_dLambda()) {
        m_RestraintManager.Print_dU_dLambda_Info();
        if (m_RestraintManager.ThereIsAForcingRestraint()) {
          m_LambdaManager.Print_dU_dLambda_Summary(Sum_dU_dLambdas);
        }
      }
      else {
        m_LambdaManager.PrintSomeSpaces();
      }
      m_RestraintManager.PrintEnergyInfo();
      m_RestraintManager.PrintRestraintInfo();
      if (m_LambdaManager.IsEndOf_MCTI()) {
        m_LambdaManager.Print_MCTI_Integration();
      }
    }
  }
}


void ComputeFreeEnergy::user_initialize() {
//-----------------------------------------------------------------
// get char* from the input stream.  read all the input
//-----------------------------------------------------------------
  char* Str = config;

  iout << iDEBUG << "Initializing free energy.\n"; 
  iout << iDEBUG << "***********************************\n"; 
  iout << Str;
  iout << iDEBUG << "***********************************\n" << endi; 

  ReadInput(Str, m_RestraintManager, m_LambdaManager, *this, simParams->dt);

  // exit if there aren't enough steps to complete all pmf & mcti blocks
  int Total = m_LambdaManager.GetTotalNumSteps();
  if (Total > simParams->N) {
    iout << "FreeEnergy: Not enough steps to complete pfm & mcti blocks" << endl;
    iout << "FreeEnergy:   Num Steps Needed =    " << Total << endl;
    iout << "FreeEnergy:   Num Steps Requested = " << simParams->N << endl << endi;
    NAMD_die("FreeEnergy: Fatal Run-Time Error");
  }
}


void ComputeFreeEnergy::user_calculate() {
//-----------------------------------------------------------------
// this is what's executed every time-step
//-----------------------------------------------------------------
  m_LambdaManager.IncCurrStep();
  update();
}


int ComputeFreeEnergy::
	getAtomID(const char *segid, int resid, const char *aname)
{
  return molecule->get_atom_from_name(segid,resid,aname);
}

int ComputeFreeEnergy::
	getNumAtoms(const char* segid, int resid) // 0 on error
{
  return molecule->get_residue_size(segid,resid);
}

int ComputeFreeEnergy::
	getAtomID(const char *segid, int resid, int index)
{
  return molecule->get_atom_from_index_in_residue(segid,resid,index);
}

double ComputeFreeEnergy::getMass(int atomid)
{
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1.;  // failure
  return molecule->atommass(atomid);
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
  configMsg = 0;  resultsMsg = 0;  config = 0;
}

ComputeFreeEnergy::~ComputeFreeEnergy() {
  DebugM(3,"Destructing ComputeFreeEnergy\n");
  delete config;
}


void ComputeFreeEnergy::initialize() {
  DebugM(4,"Initializing master\n");

  configMsg = new (MsgIndex(ComputeGlobalConfigMsg)) ComputeGlobalConfigMsg;

  // Get our script
  StringList *script = Node::Object()->configList->find("freeEnergyConfig");

  ostrstream oconfig;

  for ( ; script; script = script->next) {
    if ( script->data[0] == '{' ) {  // this is a flag, no } at end
      oconfig << script->data + 1;    // so skip it at beginning
    } else {
      ifstream infile(script->data);
      if ( infile ) infile.get(*oconfig.rdbuf(),'\0');
      if ( ! infile ) {
	char errmsg[100];
	sprintf(errmsg,"Error trying to read file %s!\n",script->data);
	NAMD_die(errmsg);
      }
    }
  }
  oconfig.flush();
  char *configstr = oconfig.str();
  config = configstr;

  iout << iDEBUG << "Free energy perturbation - initialize()\n" << endi; 
  user_initialize();

  // Send config to clients
  host->comm->sendComputeGlobalConfig(configMsg);
  configMsg = 0;
}


void ComputeFreeEnergy::calculate() {
  DebugM(4,"Calculating forces on master\n");

  resultsMsg = new (MsgIndex(ComputeGlobalResultsMsg)) ComputeGlobalResultsMsg;

//  iout << iDEBUG << "Free energy perturbation - calculate()\n" << endi; 
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
 *	$Revision: 1.13 $	$Date: 1998/09/20 16:34:55 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeFreeEnergy.C,v $
 * Revision 1.13  1998/09/20 16:34:55  hurwitz
 * make sure Lambda control objects start and stop on just the right step.
 * made output shorter and more readable (compile with _VERBOSE_PMF for old output)
 * : ----------------------------------------------------------------------
 *
 * Revision 1.12  1998/06/05 22:54:38  hurwitz
 * accumulate dU/dLambda for free energy calculation
 *
 * Revision 1.11  1998/05/25 21:55:04  jim
 * Eliminated compile errors in KCC by avoiding use of istrstream.
 *
 * Revision 1.10  1998/05/22 19:08:29  hurwitz
 * Do NAMD_die if there aren't enough steps to complete all pmf & mcti blocks
 *
 * Revision 1.9  1998/05/21 22:37:31  hurwitz
 * initial check in of code for fixed and forcing restraints
 * -Dave Hurwitz
 *
 * Revision 1.8  1998/04/30 04:53:22  jim
 * Added forces from MDComm and other improvements to ComputeGlobal.
 *
 * Revision 1.7  1998/03/26 23:28:26  jim
 * Small changes for KCC port.  Altered use of strstream in ComputeFreeEnergy.
 *
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
