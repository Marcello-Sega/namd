//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef COMPUTEMGR_H
#define COMPUTEMGR_H

#include "charm++.h"
#include "main.h"
#include "new.h"

#include "NamdTypes.h"
#include "BOCgroup.h"

#include "ResizeArray.h"

class Compute;
class ComputeMap;
class CkQdMsg;

class ComputeGlobal;
class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;

class ComputeDPME;
class ComputeDPMEDataMsg;
class ComputeDPMEResultsMsg;

class ComputePme;
class ComputePmeDataMsg;
class ComputePmeResultsMsg;

class ComputeMgr : public BOCclass
{
public:

  ComputeMgr();
  ~ComputeMgr();
  void createComputes(ComputeMap *map);
  void updateComputes(int,int);
  void updateComputes2(CkQdMsg *);
  void updateComputes3();
  void updateLocalComputes();
  void updateLocalComputes2(CkQdMsg *);
  void updateLocalComputes3();
  void doneUpdateLocalComputes();

  void sendComputeGlobalConfig(ComputeGlobalConfigMsg *);
  void recvComputeGlobalConfig(ComputeGlobalConfigMsg *);
  void sendComputeGlobalData(ComputeGlobalDataMsg *);
  void recvComputeGlobalData(ComputeGlobalDataMsg *);
  void sendComputeGlobalResults(ComputeGlobalResultsMsg *);
  void recvComputeGlobalResults(ComputeGlobalResultsMsg *);

  void sendComputeDPMEData(ComputeDPMEDataMsg *);
  void recvComputeDPMEData(ComputeDPMEDataMsg *);
  void sendComputeDPMEResults(ComputeDPMEResultsMsg *, int);
  void recvComputeDPMEResults(ComputeDPMEResultsMsg *);

  void sendComputePmeData(ComputePmeDataMsg *);
  void recvComputePmeData(ComputePmeDataMsg *);
  void sendComputePmeResults(ComputePmeResultsMsg *, int);
  void recvComputePmeResults(ComputePmeResultsMsg *);

private:
  void createCompute(ComputeID, ComputeMap *);
  int numNonbondedSelf;
  int numNonbondedPair;

  ComputeGlobal *computeGlobalObject;
  ComputeDPME *computeDPMEObject;
  ComputePme *computePmeObject;

  int updateComputesCount;
  int updateComputesReturnEP;
  int updateComputesReturnChareID;

  int *computeFlag;

  class ComputeElem {
  public:
    ComputeID   cid;
    Compute *c;

    // new matched with delete
    void * operator new(size_t size) { return ::operator new(size); }	
    // placement new for use within containers
    void * operator new(size_t /* size */, void * ptr) { return ptr; }

    // delete matched with new
    void operator delete(void* ptr) { ::operator delete(ptr); }		
    // error on IRIX!  void operator delete(void*, void*) { ; }		

    int operator<(ComputeElem e) { return (cid < e.cid); }
    int operator==(ComputeElem e) { return (cid == e.cid); }

    ComputeElem(ComputeID id=-1, Compute *compute=NULL) : 
      cid(id), c(compute) {};
    ~ComputeElem() { };
    ComputeElem& operator=(const ComputeElem &e) { 
      cid = e.cid; c = e.c;  // Do not delete c!  This op only used to shuffle
                             // we delete the c here only when the Compute is 
		             // moved off!
      return(*this);
    };
  };

  int numComputes;

  typedef ResizeArray<ComputeElem> ComputeList;
  typedef ResizeArray<int> ComputeIndex;

  // global patch number to local patch table conversion table
  ComputeIndex computeIndex;

  // an array of compute pointers residing on this node
  ComputeList computeList;

  int workDistribGroup;
};

#endif /* COMPUTEMGR_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMgr.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1016 $	$Date: 1999/06/08 14:52:06 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMgr.h,v $
 * Revision 1.1016  1999/06/08 14:52:06  jim
 * Incorporated Justin's faster PME code along side DPME.
 *
 * Revision 1.1015  1999/05/11 23:56:26  brunner
 * Changes for new charm version
 *
 * Revision 1.1014  1999/03/17 17:59:22  jim
 * Eliminated compiler warnings and errors.
 *
 * Revision 1.1013  1998/04/10 04:16:01  jim
 * Finished incorporating DPME.
 *
 * Revision 1.1012  1998/03/26 23:28:27  jim
 * Small changes for KCC port.  Altered use of strstream in ComputeFreeEnergy.
 *
 * Revision 1.1011  1998/03/03 23:05:11  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1010  1997/12/26 23:10:46  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1009  1997/12/19 23:48:50  jim
 * Added Tcl interface for calculating forces.
 *
 * Revision 1.1008  1997/04/17 19:19:36  brunner
 * Put in new AlgSeven.C
 *
 * Revision 1.1007  1997/04/10 09:13:52  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of cout and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1006  1997/04/08 07:08:20  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1005  1997/03/20 23:53:42  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1004  1997/03/19 05:49:57  jim
 * Added ComputeSphericalBC, cleaned up make dependencies.
 *
 * Revision 1.1003  1997/03/04 22:37:09  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1002  1997/02/14 19:07:31  nealk
 * Added new/delete comments.
 * Played with DPMTA.
 *
 * Revision 1.1001  1997/02/07 17:39:37  ari
 * More debugging for atomMigration.
 * Using -w on CC got us some minor fixes
 * using purify got us a major memory problem due to bad sizing of dummy force
 *
 * Revision 1.1000  1997/02/06 15:58:05  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:18  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:03  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:35:51  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.3  1996/12/16 23:16:56  jim
 * eliminated warning about multiple new's
 *
 * Revision 1.2  1996/11/30 00:34:05  jim
 * added some includes, now uses InitMsg
 *
 * Revision 1.1  1996/11/27 20:19:59  jim
 * Initial revision
 *
 *
 ***************************************************************************/
