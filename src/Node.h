//-*-c++-*-
//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************/
/* DESCRIPTION:							           */	
/*									   */
/***************************************************************************/

#ifndef _NODE_H
#define _NODE_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.h"
#include "PatchMap.h"
#include "AtomMap.h"
#include "ComputeMap.h"
#include "PatchMgr.h"
#include "ProxyMgr.h"

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

class Node : public BOCclass
{
public:

  Node(GroupInitMsg *msg);
  ~Node(void);

  // Singleton Access method
  inline static Node *Object() {return _instance;}

  // Run for the number of steps specified in the sim_parameters
  static void messageRun();
  void run(RunMsg *);                  

  // Deal with quiescence
  void quiescence(QuiescenceMessage *);

  // Charm Entry point - Read in system data, get all ready to simulate
  static void messageStartup();
  void startup(InitMsg *);  
  void startup1(InitMsg *);  
  void messageStartup2(QuiescenceMessage *);
  void startup2(InitMsg *);  
  void messageStartup3(QuiescenceMessage *);
  void startup3(InitMsg *);  
  void messageStartup4(QuiescenceMessage *);
  void startup4(InitMsg *);  

  static void messageStartupDone();
  void startupDone(DoneMsg *);

  // Charm Entry point - synchronize on BOC creation and startup
  static void messageBOCCheckIn();
  void BOCCheckIn(InitMsg *msg);
  void awaitBOCCheckIn();

  // Utility for storing away simulation data for Node
  void saveMolDataPointers(Molecule *, Parameters *,
			   SimParameters *, ConfigList *,
			   PDB *, NamdState *);


  // NAMD 1.X molecule database objects - must be public for now
  Molecule *molecule;
  Parameters *parameters;
  SimParameters *simParameters;
  ConfigList *configList;
  PDB *pdb;
  NamdState *state;
  Output *output;

  // Remove these calls?
  int myid() { return CMyPe(); }
  int numNodes() { return CNumPes(); }

protected:
  // Map Databases - they have a singleton this access method ::Object()
  AtomMap    *atomMap;
  PatchMap   *patchMap;
  ComputeMap *computeMap;

private:
  static Node *_instance;

  WorkDistrib *workDistrib;
  PatchMgr *patchMgr;
  ComputeMgr *computeMgr;
  ProxyMgr *proxyMgr;

  // Countdown for Node::startup barrier
  int numNodeStartup;

};

#endif /* _NODE_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/02/11 22:56:14 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.h,v $
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
