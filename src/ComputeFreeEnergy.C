/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include <string.h>
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
#include "SimParameters.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
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
// read all the input from config
//-----------------------------------------------------------------

  iout << iINFO << "  FREE ENERGY PERTURBATION CONFIG\n"; 
  iout << iINFO << "***********************************\n"; 
  iout << config;
  iout << iINFO << "***********************************\n" << endi; 

  ReadInput(config, m_RestraintManager, m_LambdaManager, *this, simParams->dt);

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

  configMsg = new ComputeGlobalConfigMsg;

  // Get our script
  StringList *script = Node::Object()->configList->find("freeEnergyConfig");

  config = new char[1];
  config[0] = '\0';

  for ( ; script; script = script->next) {
    if ( script->data[0] == '{' ) {
      // this is a flag, no '}' at end so skip it at beginning
      size_t add_len = strlen(script->data + 1);
      size_t config_len = 0;
      config_len = strlen(config);
      char *new_config = new char[config_len + add_len + 2];
      strcpy(new_config,config);
      strcat(new_config,script->data + 1);
      strcat(new_config,"\n");  // just to be safe
      delete [] config;
      config = new_config;
    } else {
      FILE *infile = fopen(script->data,"r");
      if ( ! infile ) {
	char errmsg[256];
	sprintf(errmsg,"Error trying to read file %s!\n",script->data);
	NAMD_die(errmsg);
      }
      fseek(infile,0,SEEK_END);
      size_t add_len = ftell(infile);
      size_t config_len = 0;
      config_len = strlen(config);
      char *new_config = new char[config_len + add_len + 3];
      strcpy(new_config,config);
      delete [] config;
      config = new_config;
      new_config += config_len;
      rewind(infile);
      fread(new_config,sizeof(char),add_len,infile);
      new_config += add_len;
      new_config[0] = '\n';
      new_config[1] = '\0';
      fclose(infile);
    }
  }

  // iout << iDEBUG << "Free energy perturbation - initialize()\n" << endi; 
  user_initialize();

  // Send config to clients
  host->comm->sendComputeGlobalConfig(configMsg);
  configMsg = 0;
}


void ComputeFreeEnergy::calculate() {
  DebugM(4,"Calculating forces on master\n");

  resultsMsg = new ComputeGlobalResultsMsg;
  resultsMsg->gforce.resize(gmass.size());

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

