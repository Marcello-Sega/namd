//-*-c++-*-
/***************************************************************************/
/*              (C) Copyright 1996,1997 The Board of Trustees of the       */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Master BOC.  coordinates startup, close down of each PE
 *		Also owns pointers to common objects needed by system		
 *		Many utility static methods are owned by Node.
 ***************************************************************************/

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
class NamdState;
class Output;
class LdbCoordinator;
class SMDData;
class SMDDataMsg;

class Node : public BOCclass
{
public:

  Node(GroupInitMsg *msg);
  ~Node(void);

  // Singleton Access method
  inline static Node *Object() {return CpvAccess(Node_instance);}

  // Run for the number of steps specified in the sim_parameters
  static void messageRun();
  void run();                  

  // Signal HomePatch and Node is done
  static void messageHomeDone();
  void homeDone();
  void nodeDone();

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

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1012 $	$Date: 1999/05/11 23:56:39 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.h,v $
 * Revision 1.1012  1999/05/11 23:56:39  brunner
 * Changes for new charm version
 *
 * Revision 1.1011  1998/03/03 23:05:19  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1010  1998/01/05 20:32:24  sergei
 * added public member SMDData smdData
 * added functions recvSMDData(), sendSMDData()
 *
 * Revision 1.1009  1997/11/07 20:17:42  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1008  1997/08/12 22:18:42  milind
 * Made NAMD2 to link on Solaris machines.
 *
 * Revision 1.1007  1997/03/27 20:25:50  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1006  1997/03/19 05:50:12  jim
 * Added ComputeSphericalBC, cleaned up make dependencies.
 *
 * Revision 1.1005  1997/03/18 18:09:09  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1004  1997/03/14 21:40:13  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1003  1997/03/04 22:37:15  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1002  1997/02/13 16:17:17  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1001  1997/02/11 22:56:14  jim
 * Added dcd file writing.
 *
 * Revision 1.1000  1997/02/06 15:58:52  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:00  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:28  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:35  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.18  1997/01/09 20:48:10  jim
 * added Controller code
 *
 * Revision 1.17  1996/12/20 22:53:27  jim
 * fixing parallel bugs, going home for Christmass
 *
 * Revision 1.16  1996/12/19 02:26:30  jim
 * Node::startup2 is now triggered by quiescence
 *
 * Revision 1.15  1996/12/13 08:52:37  jim
 * staged startup implemented
 *
 * Revision 1.14  1996/12/12 20:14:50  milind
 * *** empty log message ***
 *
 * Revision 1.13  1996/12/06 19:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.12  1996/11/30 00:44:24  jim
 * added sequencer use, ComputeMgr use, and quiescence detection
 *
 * Revision 1.11  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.10  1996/10/29 23:35:27  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/09/03 22:54:25  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/08/23 22:03:52  brunner
 * Made WorkdDistrib, PatchMgr public members
 *
 * Revision 1.6  1996/08/21 23:58:25  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/08/19 22:05:31  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/19 17:57:47  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/16 21:56:17  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 21:42:58  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 21:19:34  ari
 * Initial revision
 *
 *
 ***************************************************************************/
