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

class Molecule;
class Parameters;
class SimParameters;
class ConfigList;
class PDB;
class WorkDistrib;
class PatchMgr;
class ComputeMgr;

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

  static void messageStartupDone();
  void startupDone(DoneMsg *);

  // Charm Entry point - synchronize on BOC creation and startup
  static void messageBOCCheckIn();
  void BOCCheckIn(InitMsg *msg);
  void awaitBOCCheckIn();

  // Utility for storing away simulation data for Node
  void saveMolDataPointers(Molecule *, Parameters *,
			   SimParameters *, ConfigList *,
			   PDB *);


  // NAMD 1.X molecule database objects - must be public for now
  Molecule *molecule;
  Parameters *parameters;
  SimParameters *simParameters;
  ConfigList *configList;
  PDB *pdb;

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

  // Countdown for Node::startup barrier
  int numNodeStartup;

};

#endif /* _NODE_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.12 $	$Date: 1996/11/30 00:44:24 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.h,v $
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
