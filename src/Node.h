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
class SMDData;
class SMDDataMsg;
class ScriptTcl;
class IMDOutput;

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
  inline static Node *Object() {return CpvAccess(Node_instance);}

  void startupCont(CkQdMsg *);
#ifdef NAMD_TCL
  void enableStartupCont(Namd *);
  ScriptTcl *getScript(void) { return script; }
#endif

  // Run for the number of steps specified in the sim_parameters
  static void messageRun();
  void run();                  

  // Change parameters in mid-run
  void enableScriptBarrier();  
  void scriptBarrier(CkQdMsg *);  
  void scriptParam(ScriptParamMsg *);

  // End of run
  void enableHaltBarrier();  
  void haltBarrier(CkQdMsg *);  

  // Deal with quiescence
  void quiescence(CkQdMsg *);

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

  // Init the socket connect for imd
  void IMDinit(void *);

  // Deal with SMD data message
  void sendSMDData(SMDDataMsg *);
  void recvSMDData(SMDDataMsg *);

  // NAMD 1.X molecule database objects - must be public for now
  Molecule *molecule;
  Parameters *parameters;
  SimParameters *simParameters;
  ConfigList *configList;
  PDB *pdb;
  NamdState *state;
  Output *output;
  SMDData *smdData;
  IMDOutput *imd;

  // Remove these calls?
  int myid() { return CkMyPe(); }
  int numNodes() { return CkNumPes(); }

protected:
  // Map Databases - they have a singleton this access method ::Object()
  AtomMap    *atomMap;
  PatchMap   *patchMap;
  ComputeMap *computeMap;
  LdbCoordinator *ldbCoordinator;

private:
  void namdOneCommInit();
  void namdOneRecv();
  void namdOneSend();
  void threadInit();
  void buildSequencers();

  WorkDistrib *workDistrib;
  PatchMgr *patchMgr;
  ComputeMgr *computeMgr;
  ProxyMgr *proxyMgr;
  ScriptTcl *script;
#ifdef NAMD_TCL
  Namd *namd;
#endif

  // Countdown for Node::startup barrier
  int numNodeStartup;

  // Countdown for Node::homeDone termination 
  int numHomePatchesRunning;

  // Countdown for Node::nodeDone termination
  int numNodesRunning;

  // Startup phase
  int startupPhase;
};

#endif /* _NODE_H */

