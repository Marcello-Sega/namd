/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeGlobal.h"
#include "ComputeMisc.h"
#include "ComputeGlobalMsgs.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "SimParameters.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include <stdio.h>

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


ComputeGlobalEasy::ComputeGlobalEasy(ComputeGlobal *h, const char *cn)
  : ComputeGlobalMaster(h) {

  molecule = Node::Object()->molecule;
  simParams = Node::Object()->simParameters;
  configMsg = 0;  resultsMsg = 0;

  int len = strlen(cn);
  configName = new char[len+1];
  strcpy(configName,cn);

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
}

ComputeGlobalEasy::~ComputeGlobalEasy() {
  delete [] configName;
  delete reduction;
}

int ComputeGlobalEasy::
        getAtomID(const char *segid, int resid, const char *aname)
{
  return molecule->get_atom_from_name(segid,resid,aname);
}

int ComputeGlobalEasy::
        getNumAtoms(const char* segid, int resid) // 0 on error
{
  return molecule->get_residue_size(segid,resid);
}

int ComputeGlobalEasy::
        getAtomID(const char *segid, int resid, int index)
{
  return molecule->get_atom_from_index_in_residue(segid,resid,index);
}

double ComputeGlobalEasy::getMass(int atomid)
{
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1.;  // failure
  return molecule->atommass(atomid);
}


int ComputeGlobalEasy::requestAtom(int atomid)
{
  if ( ! configMsg ) return -1;  // failure
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1;  // failure
  configMsg->aid.add(atomid);
  return 0;  // success
}

int ComputeGlobalEasy::getPosition(int atomid, Position &position)
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

int ComputeGlobalEasy::addForce(int atomid, Force force)
{
  if ( ! resultsMsg ) return -1;  // failure
  if ( atomid < 0 || atomid >= molecule->numAtoms ) return -1;  // failure
  resultsMsg->aid.add(atomid);
  resultsMsg->f.add(force);
  return 0;  // success
}

void ComputeGlobalEasy::addEnergy(BigReal e) {
  energy += e;
}

void ComputeGlobalEasy::initialize() {
  DebugM(4,"Initializing master\n");

  configMsg = new ComputeGlobalConfigMsg;

  // Get the script for subclasses
  StringList *script = Node::Object()->configList->find(configName);

  char *config = new char[1];
  config[0] = '\0';

  for ( ; script; script = script->next) {
    if ( script->data[0] == '{' ) {
      script->data[strlen(script->data)-1] = 0;
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

  // Initialize subclass
  easy_init(config);

  delete [] config;

  // Send config to clients
  host->comm->sendComputeGlobalConfig(configMsg);
  configMsg = 0;
}

void ComputeGlobalEasy::calculate() {

  resultsMsg = new ComputeGlobalResultsMsg;
  resultsMsg->gforce.resize(gmass.size());
  energy = 0.0;

  // Build results here
  easy_calc();

  // Send results to clients
  host->comm->sendComputeGlobalResults(resultsMsg);
  resultsMsg = 0;
  reduction->item(REDUCTION_MISC_ENERGY) += energy;
  reduction->submit();
}

void ComputeGlobalEasy::easy_init(const char *) {
  CkPrintf("Default ComputeGlobalEasy::easy_init() called.\n");
}

void ComputeGlobalEasy::easy_calc() {
  CkPrintf("Default ComputeGlobalEasy::easy_calc() called.\n");
}

