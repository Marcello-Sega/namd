/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   HomePatch is the key distributed source/sink of Atom data
   including positions, velocities and forces applied
*/

#ifndef HOMEPATCH_H
#define HOMEPATCH_H

#include "charm++.h"

#include "NamdTypes.h"
#include "Patch.h"
#include "PatchMap.h"

#include "MigrateAtomsMsg.h"
#include "HomePatchTypes.h"
#include "main.h"
#include "common.h"
#include "Migration.h"

class RegisterProxyMsg;
class UnregisterProxyMsg;
class ProxyResultMsg;
class ProxyCombinedResultMsg;
class Sequencer;

class HomePatch : public Patch {
  friend class PatchMgr;
  friend class Sequencer;

private: 
  // for PatchMgr to use only
  HomePatch(PatchID, FullAtomList);
  ScaledPosition min, max, center;

public:
  ~HomePatch();

  // Message from ProxyPatch (via ProxyMgr) which registers its existence
  void registerProxy(RegisterProxyMsg *);
  // opposite of above
  void unregisterProxy(UnregisterProxyMsg *);

  // ProxyPatch sends Forces back to here (via ProxyMgr)
  void receiveResults(ProxyResultMsg *msg);
  void receiveResults(ProxyCombinedResultMsg *msg);

  // AtomMigration messages passes from neighbor HomePatches to here.
  void depositMigration(MigrateAtomsMsg *);

  // Bind a Sequencer to this HomePatch
  void useSequencer(Sequencer *sequencerPtr);
  // start simulation over this Patch of atoms
  void runSequencer(void);
  
  //--------------------------------------------------------------------
  // methods for Sequencer to use
  //

  // Signal HomePatch that positions stored are to be now to be used
  void positionsReady(int doMigration=0);

  // methods to implement integration
  void saveForce(const int ftag = Results::normal);
  void addForceToMomentum(const BigReal, const int ftag = Results::normal,
				const int useSaved = 0);
  void addVelocityToPosition(const BigReal);

  // methods for rigidBonds
  int rattle1(const BigReal);
  void rattle2(const BigReal, Tensor *virial);

  // methods for mollified impluse (MOLLY)
  void mollyAverage();
  void mollyMollify(Tensor *virial);
//  Bool average(Vector qtilde[],const Vector q[],BigReal lambda[],const int n,const int m, const BigReal imass[], const BigReal length2[], const int ial[], const int ilb[], const Vector qji[], const BigReal tolf, const int ntrial);
//  void mollify(Vector qtilde[],const Vector q0[],const BigReal lambda[], Vector force[],const int n, const int m, const BigReal imass[],const int ial[],const int ibl[],const Vector refab[]); 

  // methods for CONTRA, etc
  void checkpoint(void);
  void revert(void);

  // methods for QM (ExtForces replacement)
  void replaceForces(ExtForce *f);

  // patch-wise calculations
  BigReal calcKineticEnergy();
  Vector calcMomentum();
  Vector calcAngularMomentum();

  // load-balancing trigger
  void submitLoadStats(int timestep);

  // for ComputeHomePatches
  FullAtomList &getAtomList() { return (atom); }

  // build spanning tree for proxy nodes
  void buildSpanningTree(void);
protected:
  virtual void boxClosed(int);

  // Internal Atom Migration methods and data
  void doGroupSizeCheck();
  void doMarginCheck();
  void doAtomMigration();
  int inMigration;
  int numMlBuf;
  MigrateAtomsMsg *msgbuf[PatchMap::MaxOneAway];
  
private:
  // Store of Atom-wise variables
  FullAtomList  atom;
  ForceList f_saved[Results::maxNumForces];
  ExtForce *replacementForces;

  // checkpointed state
  FullAtomList  checkpoint_atom;
  Lattice  checkpoint_lattice;

  // MOLLY data
  ResizeArray<BigReal> molly_lambda;
  
  // List of Proxies
  ProxyList     proxy;
  
  Sequencer  *sequencer;

  // Needed for initialization
  int patchMapRead;
  void readPatchMap();

  // Atom Migration internals
  int allMigrationIn;
  int migrationSuspended;
  int patchMigrationCounter;
  int numNeighbors;
  MigrationInfo realInfo[PatchMap::MaxOneAway];
  MigrationInfo *mInfo[3][3][3];

  int child[PROXY_SPAN_DIM];	// spanning tree of proxies
  int nChild;
};

#endif

