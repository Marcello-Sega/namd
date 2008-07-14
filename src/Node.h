/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Master BOC.  coordinates startup, close down of each PE
   Also owns pointers to common objects needed by system		
   Many utility static methods are owned by Node.
*/

#ifndef _NODE_H
#define _NODE_H

#include "charm++.h"

#include "main.h"

#ifdef CHARMIZE_NAMD
#include "AtomsDisInfo.h"
#endif

#ifdef SOLARIS
extern "C" int gethostname( char *name, int namelen);
#endif

#include "ProcessorPrivate.h"
#include "Node.decl.h"

class PatchMap;
class AtomMap;
class ProxyMgr;
class ComputeMap;
class PatchMgr;
class Molecule;
class Parameters;
class SimParameters;
class ConfigList;
class PDB;
class WorkDistrib;
class PatchMgr;
class ComputeMgr;
class Communicate;
class Namd;
class NamdState;
class Output;
class LdbCoordinator;
class ScriptTcl;
class IMDOutput;
class Vector;

// Message to send our per processor BOC's list of groupIDs of
// all other BOC's
class GroupInitMsg : public CMessage_GroupInitMsg
{
public:
  BOCgroup group;
};

class AllCharmArrsMsg : public CMessage_AllCharmArrsMsg
{
#ifdef CHARMIZE_NAMD
public:
  CProxy_AtomsDisInfo atomsDis;
#endif
};

#define MAX_SCRIPT_PARAM_SIZE 128
class ScriptParamMsg : public CMessage_ScriptParamMsg {
public:
  char param[MAX_SCRIPT_PARAM_SIZE];
  char value[MAX_SCRIPT_PARAM_SIZE];
};

class Node : public BOCclass
{
public:

  Node(GroupInitMsg *msg);
  ~Node(void);

  // Singleton Access method
  inline static Node *Object() {return CkpvAccess(Node_instance);}

  // Run for the number of steps specified in the sim_parameters
  static void messageRun();
  void run();                  

  // Change parameters in mid-run
  void enableScriptBarrier();  
  void scriptBarrier(CkQdMsg *);  
  void scriptParam(ScriptParamMsg *);

  void reloadCharges(const char *filename);
  void reloadCharges(float charge[], int n);

  void sendEnableExitScheduler(void);
  void recvEnableExitScheduler(CkQdMsg *);
  void enableExitScheduler(void);
  void exitScheduler(CkQdMsg *);

  void sendEnableEarlyExit(void);
  void recvEnableEarlyExit(CkQdMsg *);
  void enableEarlyExit(void);
  void earlyExit(CkQdMsg *);

  // Charm Entry point - Read in system data, get all ready to simulate
  static void messageStartUp();
  void startup();  
  void startUp(CkQdMsg *);  

  // Charm Entry point - synchronize on BOC creation and startup
  static void messageBOCCheckIn();
  void BOCCheckIn();
  void awaitBOCCheckIn();

  // Utility for storing away simulation data for Node
  void saveMolDataPointers(NamdState *);

  //#ifdef CHARMIZE_NAMD  
  //It is an entry method which has to be declared
  void sendCharmArrProxies(AllCharmArrsMsg *msg);
  //#endif

  // Made public for pmeAid;
  WorkDistrib *workDistrib;

  // Made public in order to access the ComputeGlobal on the node
  ComputeMgr *computeMgr;
  
  // NAMD 1.X molecule database objects - must be public for now
  Molecule *molecule;
  Parameters *parameters;
  SimParameters *simParameters;
  ConfigList *configList;
  PDB *pdb;
  NamdState *state;
  Output *output;
  IMDOutput *imd;
  Vector *coords;  // Only exists during measure from Tcl

  #ifdef CHARMIZE_NAMD
  CProxy_AtomsDisInfo atomDisArr;
  #endif


  // Remove these calls?
  int myid() { return CkMyPe(); }
  int numNodes() { return CkNumPes(); }

  void setScript(ScriptTcl *s) { script = s; }
  ScriptTcl *getScript(void) { return script; }

protected:
  // Map Databases - they have a singleton this access method ::Object()
  AtomMap    *atomMap;
  PatchMap   *patchMap;
  ComputeMap *computeMap;
  LdbCoordinator *ldbCoordinator;

private:
  #ifdef CHARMIZE_NAMD  
  void populateAtomDisArrs(int startupPhase);
  #endif
  
  void namdOneCommInit();
  void namdOneRecv();
  void namdOneSend();
  void threadInit();
  void buildSequencers();

  PatchMgr *patchMgr;
  ProxyMgr *proxyMgr;
  Namd *namd;
  ScriptTcl *script;

  // Startup phase
  int startupPhase;  
};

#endif /* _NODE_H */

