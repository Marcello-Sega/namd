/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Namd() - launches BOC's Node, WorkDistrib, PatchMgr, ProxyMgr
   ReductionMgr, CollectionMgr
*/

#include "unistd.h"

#include "charm++.h"

#include "Namd.h"
#include "NamdState.h"
#include "Node.h"

float Namd::cmiWallStart;
float Namd::cmiCpuStart;
int Namd::cmiFirstStart;
float Namd::cmiWallFirstStart;
float Namd::cmiCpuFirstStart;

// Namd(void ) is the constructor for the startup node.  It needs to
// read in file data,
Namd::Namd(void)
{
  namdState = new NamdState;
}


// ~Namd(void) just needs to tell all the slave nodes to die.
Namd::~Namd(void)
{
  delete namdState;
  CkPrintf("Namd::~Namd() called\n");
}


// startup(char *) 
void Namd::startup(char *confFile)
{
  namdState->configFileInit(confFile);

#ifdef NAMD_TCL
  Node::Object()->enableStartupCont(this);
}

void Namd::startupCont(void) {
  namdState->configFileInitCont();
#endif

  if (namdState->status()) {
    CkPrintf("Namd::startup() - could not initialize namdState.\n");
    CkExit();
  }

  // Give our node[PE = 0] pointers to the data objects, so it can use them,
  // or send them on as messages elsewhere.
  Node::Object()->saveMolDataPointers(namdState);

  Node::messageStartUp(); // tell all nodes to startup
}

// last call of system
void Namd::namdDone(void) {
    Real CPUtime = CmiCpuTimer()-cmiCpuFirstStart;
    Real Walltime = CmiWallTimer()-cmiWallFirstStart;
    CkPrintf("==========================================\n"
	"WallClock : %f  CPUTime : %f \n",Walltime,CPUtime);
    CkExit();
}

