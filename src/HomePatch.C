/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   HomePatch owns the actual atoms of a Patch of space
   Proxy(s) get messages via ProxyMgr from HomePatch(es)
   to update lists of atoms and their coordinates
   HomePatch(es) also have a Sequencer bound to them

   superclass: 	Patch		
*/

#include <math.h>
#include "charm++.h"

#include "SimParameters.h"
#include "HomePatch.h"
#include "AtomMap.h"
#include "Node.h"
#include "PatchMap.inl"
#include "main.h"
#include "ProxyMgr.decl.h"
#include "ProxyMgr.h"
#include "Migration.h"
#include "Molecule.h"
#include "PatchMgr.h"
#include "Sequencer.h"
#include "LdbCoordinator.h"

#include "Sync.h"

#define TINY 1.0e-20;
#define MAXHGS 10
#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

typedef int HGArrayInt[MAXHGS];
typedef BigReal HGArrayBigReal[MAXHGS];
typedef Vector HGArrayVector[MAXHGS];
typedef BigReal HGMatrixBigReal[MAXHGS][MAXHGS];
typedef Vector HGMatrixVector[MAXHGS][MAXHGS];

int average(CompAtom *qtilde,const HGArrayVector &q,BigReal *lambda,const int n,const int m, const HGArrayBigReal &imass, const HGArrayBigReal &length2, const HGArrayInt &ial, const HGArrayInt &ibl, const HGArrayVector &refab, const BigReal tolf, const int ntrial);

void mollify(CompAtom *qtilde,const HGArrayVector &q0,const BigReal *lambda, HGArrayVector &force,const int n, const int m, const HGArrayBigReal &imass,const HGArrayInt &ial,const HGArrayInt &ibl,const HGArrayVector &refab);


HomePatch::HomePatch(PatchID pd, FullAtomList al) : Patch(pd), atom(al)
{ 
  min.x = PatchMap::Object()->min_a(patchID);
  min.y = PatchMap::Object()->min_b(patchID);
  min.z = PatchMap::Object()->min_c(patchID);
  max.x = PatchMap::Object()->max_a(patchID);
  max.y = PatchMap::Object()->max_b(patchID);
  max.z = PatchMap::Object()->max_c(patchID);
  center = 0.5*(min+max);

  migrationSuspended = false;
  allMigrationIn = false;
  patchMapRead = 0; // We delay read of PatchMap data
		    // to make sure it is really valid
  inMigration = false;
  numMlBuf = 0;
  flags.sequence = -1;

  numAtoms = atom.size();
  replacementForces = 0;

  nChild = 0;	// number of proxy spanning tree children
}

// Bind a Sequencer to this HomePatch
void HomePatch::useSequencer(Sequencer *sequencerPtr)
{ sequencer=sequencerPtr; }

// start simulation over this Patch of atoms
void HomePatch::runSequencer(void)
{ sequencer->run(); }

void HomePatch::readPatchMap() {
  // iout << "Patch " << patchID << " has " << proxy.size() << " proxies.\n" << endi;
  PatchMap *p = PatchMap::Object();
  PatchID nnPatchID[PatchMap::MaxOneAway];

  patchMigrationCounter = numNeighbors 
    = PatchMap::Object()->oneAwayNeighbors(patchID, nnPatchID);
  DebugM( 1, "NumNeighbors for pid " <<patchID<<" is "<< numNeighbors << "\n");
  int n;
  for (n=0; n<numNeighbors; n++) {
    realInfo[n].destNodeID = p->node(realInfo[n].destPatchID = nnPatchID[n]);
     DebugM( 1, " nnPatchID=" <<nnPatchID[n]<<" nnNodeID="<< realInfo[n].destNodeID<< "\n");
    realInfo[n].mList.resize(0);
  }

  // Make mapping from the 3x3x3 cube of pointers to real migration info
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
      {
	int pid =  p->pid(p->index_a(patchID)+i-1, 
	    p->index_b(patchID)+j-1, p->index_c(patchID)+k-1);
	if (pid < 0) {
	   DebugM(5, "ERROR, for patchID " << patchID <<" I got neigh pid = " << pid << "\n");
	}
	if (pid == patchID) {
	  mInfo[i][j][k] = NULL;
	}
	else {
	  // Does not work as expected for periodic with only two patches.
	  // Also need to check which image we want, but OK for now.  -JCP
	  for (n = 0; n<numNeighbors; n++) {
	    if (pid == realInfo[n].destPatchID) {
	      mInfo[i][j][k] = &realInfo[n];
	      break;
	    }
	  }
	  if (n == numNeighbors) { // disaster! 
	    DebugM(4,"BAD News, I could not find PID " << pid << "\n");
	  }
	}
      }

  DebugM(1,"Patch("<<patchID<<") # of neighbors = " << numNeighbors << "\n");
}

HomePatch::~HomePatch()
{
}


void HomePatch::boxClosed(int)
{
  if ( ! --boxesOpen )
  {
    if ( replacementForces ) {
      for ( int i = 0; i < numAtoms; ++i ) {
        if ( replacementForces[i].replace ) {
          for ( int j = 0; j < Results::maxNumForces; ++j ) { f[j][i] = 0; }
          f[Results::normal][i] = replacementForces[i].force;
        }
      }
      replacementForces = 0;
    }
    DebugM(1,patchID << ": " << CthSelf() << " awakening sequencer "
	<< sequencer->thread << "(" << patchID << ") @" << CmiTimer() << "\n");
    // only awaken suspended threads.  Then say it is suspended.
    sequencer->awaken();
    return;
  }
  else
  {
    DebugM(1,patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void HomePatch::registerProxy(RegisterProxyMsg *msg) {
  DebugM(4, "registerProxy("<<patchID<<") - adding node " <<msg->node<<"\n");
  proxy.add(ProxyListElem(msg->node,forceBox.checkOut()));
  delete msg;
}

void HomePatch::unregisterProxy(UnregisterProxyMsg *msg) {
  int n = msg->node;
  ProxyListElem *pe = proxy.begin();
  for ( ; pe->node != n; ++pe );
  forceBox.checkIn(pe->forceBox);
  proxy.del(pe - proxy.begin());
  delete msg;
}

void HomePatch::buildSpanningTree(void)
{
  nChild = 0;
  if (proxy.size() == 0) return;
  NodeIDList tree;
  tree.resize(proxy.size()+1);
  tree[0] = CkMyPe();
  int nPatches = PatchMap::Object()->numPatches();
  int s=1, e=proxy.size()+1;
  ProxyListIter pli(proxy);
  for ( pli = pli.begin(); pli != pli.end(); ++pli )
  {
    if ( pli->node < nPatches ) {
      tree[--e] = pli->node;
    } else {
      tree[s++] = pli->node;
    }
  }
  
//CkPrintf("home: %d:(%d) %d %d %d %d %d\n", patchID, tree.size(),tree[0],tree[1],tree[2],tree[3],tree[4]);
  for (int i=1; i<=PROXY_SPAN_DIM; i++) {
    if (tree.size() <= i) break;
    child[i-1] = tree[i];
    nChild++;
  }
  ProxySpanningTreeMsg *msg = new ProxySpanningTreeMsg;
  msg->patch = patchID;
  msg->node = CkMyPe();
  msg->tree = tree;
  ProxyMgr::Object()->sendSpanningTree(msg);

}

void HomePatch::receiveResults(ProxyResultMsg *msg)
{
  DebugM(4, "patchID("<<patchID<<") receiveRes() nodeID("<<msg->node<<")\n");
  int n = msg->node;
  ProxyListElem *pe = proxy.begin();
  for ( ; pe->node != n; ++pe );
  Results *r = pe->forceBox->open();
  for ( int k = 0; k < Results::maxNumForces; ++k )
  {
    Force *f = r->f[k];
    register ForceList::iterator f_i, f_e;
    f_i = msg->forceList[k].begin();
    f_e = msg->forceList[k].end();
    for ( ; f_i != f_e; ++f_i, ++f ) *f += *f_i;
  }
  pe->forceBox->close(&r);
  delete msg;
}

void HomePatch::receiveResults(ProxyCombinedResultMsg *msg)
{
  DebugM(4, "patchID("<<patchID<<") receiveRes() nodeID("<<msg->node<<")\n");
//CkPrintf("[%d] Homepatch: %d receiveResults from %d nodes\n", CkMyPe(), patchID, n);
  for (int i=0; i<msg->nodes.size(); i++) {
    int node = msg->nodes[i];
    ProxyListElem *pe = proxy.begin();
    for ( ; pe->node != node; ++pe );
    Results *r = pe->forceBox->open();
    if (i==0) {
      for ( int k = 0; k < Results::maxNumForces; ++k )
      {
        Force *f = r->f[k];
        register ForceList::iterator f_i, f_e;
        f_i = msg->forceList[k].begin();
        f_e = msg->forceList[k].end();
        for ( ; f_i != f_e; ++f_i, ++f ) *f += *f_i;
      }
    }
    pe->forceBox->close(&r);
  }
  delete msg;
}

void HomePatch::positionsReady(int doMigration)
{
  flags.sequence++;

  if (!patchMapRead) { readPatchMap(); }
      
  if (numNeighbors) {
    if (doMigration) {
      doAtomMigration();
    } else {
      doMarginCheck();
    }
  }
  doMigration = (doMigration && numNeighbors) || ! patchMapRead;

  // Workaround for oversize groups
  doGroupSizeCheck();

  // Copy information needed by computes and proxys to Patch::p.
  p.resize(numAtoms);
  CompAtom *p_i = p.begin();
  FullAtom *a_i = atom.begin();
  int i; int n = numAtoms;
  for ( i=0; i<n; ++i ) { p_i[i] = a_i[i]; }

  if (flags.doMolly) mollyAverage();

  // Must Add Proxy Changes when migration completed!
  ProxyListIter pli(proxy);
  int *pids;
  int npid;
  if (proxySendSpanning == 0) {
    npid = proxy.size();
    pids = new int[npid];
    int *pidi = pids;
    int *pide = pids + proxy.size();
    int nPatches = PatchMap::Object()->numPatches();
    for ( pli = pli.begin(); pli != pli.end(); ++pli )
    {
      if ( pli->node < nPatches ) {
        *(--pide) = pli->node;
      } else {
        *(pidi++) = pli->node;
      }
    }
  }
  else {
    npid = nChild;
    pids = new int[PROXY_SPAN_DIM];
    for (int i=0; i<nChild; i++) pids[i] = child[i];
  }
  if (npid) {
    int seq = flags.sequence;
    int priority = 64 + (seq % 256) * 256 + (patchID % 64);
    if (doMigration) {
        ProxyAllMsg *allmsg = new (sizeof(int)*8) ProxyAllMsg;
        CkSetQueueing(allmsg, CK_QUEUEING_IFIFO);
        *((int*) CkPriorityPtr(allmsg)) = priority;
        allmsg->patch = patchID;
        allmsg->flags = flags;
        allmsg->positionList = p;
        if (flags.doMolly) allmsg->avgPositionList = p_avg;
        ProxyMgr::Object()->sendProxyAll(allmsg,npid,pids);
    } else {
        ProxyDataMsg *nmsg = new (sizeof(int)*8) ProxyDataMsg;
        CkSetQueueing(nmsg, CK_QUEUEING_IFIFO);
        *((int*) CkPriorityPtr(nmsg)) = priority;
        nmsg->patch = patchID;
        nmsg->flags = flags;
        nmsg->positionList = p;
        if (flags.doMolly) nmsg->avgPositionList = p_avg;
        ProxyMgr::Object()->sendProxyData(nmsg,npid,pids);
    }   
  }
  delete [] pids;
  DebugM(4, "patchID("<<patchID<<") doing positions Ready\n");
  Patch::positionsReady(doMigration);

  patchMapRead = 1;

  // gzheng
  if (useSync) Sync::Object()->PatchReady();
}

void HomePatch::replaceForces(ExtForce *f)
{
  replacementForces = f;
}

void HomePatch::saveForce(const int ftag)
{
  f_saved[ftag].resize(numAtoms);
  for ( int i = 0; i < numAtoms; ++i )
  {
    f_saved[ftag][i] = f[ftag][i];
  }
}

void HomePatch::addForceToMomentum(const BigReal timestep, const int ftag,
							const int useSaved)
{
  const BigReal dt = timestep / TIMEFACTOR;
  if ( useSaved ) {
    for ( int i = 0; i < numAtoms; ++i )
    {
      atom[i].velocity += f_saved[ftag][i] * ( dt / atom[i].mass );
      if ( atom[i].atomFixed ) atom[i].velocity = 0;
    }
  } else {
    for ( int i = 0; i < numAtoms; ++i )
    {
      atom[i].velocity += f[ftag][i] * ( dt / atom[i].mass );
      if ( atom[i].atomFixed ) atom[i].velocity = 0;
    }
  }
}

void HomePatch::addVelocityToPosition(const BigReal timestep)
{
  const BigReal dt = timestep / TIMEFACTOR;
  for ( int i = 0; i < numAtoms; ++i )
  {
    if ( ! ( atom[i].atomFixed ) ) atom[i].position += atom[i].velocity * dt;
  }
}

//  RATTLE algorithm from Allen & Tildesley
int HomePatch::rattle1(const BigReal timestep)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal dt = timestep / TIMEFACTOR;
  BigReal tol2 = 2.0 * simParams->rigidTol;
  int maxiter = simParams->rigidIter;
  int dieOnError = simParams->rigidDie;
  int i, iter;
  BigReal dsq[10], tmp;
  int ial[10], ibl[10];
  Vector ref[10];  // reference position
  Vector refab[10];  // reference vector
  Vector pos[10];  // new position
  Vector vel[10];  // new velocity
  BigReal rmass[10];  // 1 / mass
  int fixed[10];  // is atom fixed?
  
  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      pos[i] = atom[ig+i].position;
      vel[i] = atom[ig+i].velocity;
      rmass[i] = 1. / atom[ig+i].mass;
      fixed[i] = ( atom[ig+i].atomFixed );
      // undo addVelocityToPosition to get proper reference coordinates
      if ( fixed[i] ) rmass[i] = 0.; else ref[i] -= vel[i] * dt;
    }
    int icnt = 0;
    if ( ( tmp = mol->rigid_bond_length(atom[ig].id) ) > 0 ) {  // for water
      if ( hgs != 3 ) {
        NAMD_bug("Hydrogen group error caught in rattle1().");
      }
      if ( !(fixed[1] && fixed[2]) ) {
	dsq[icnt] = tmp * tmp;  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
      }
    }
    for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = mol->rigid_bond_length(atom[ig+i].id) ) > 0 ) {
	if ( !(fixed[0] && fixed[i]) ) {
	  dsq[icnt] = tmp * tmp;  ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
	}
      }
    }
    if ( icnt == 0 ) continue;  // no constraints
    for ( i = 0; i < icnt; ++i ) {
      refab[i] = ref[ial[i]] - ref[ibl[i]];
    }
    int done;
    int consFailure;
    for ( iter = 0; iter < maxiter; ++iter ) {
      done = 1;
      consFailure = 0;
      for ( i = 0; i < icnt; ++i ) {
	int a = ial[i];  int b = ibl[i];
	Vector pab = pos[a] - pos[b];
	BigReal pabsq = pab.x*pab.x + pab.y*pab.y + pab.z*pab.z;
	BigReal rabsq = dsq[i];
	BigReal diffsq = rabsq - pabsq;
	if ( fabs(diffsq) > (rabsq * tol2) ) {
	  Vector &rab = refab[i];
	  BigReal rpab = rab.x*pab.x + rab.y*pab.y + rab.z*pab.z;
	  if ( rpab < ( rabsq * 1.0e-6 ) ) {
	    done = 0;
	    consFailure = 1;
	    continue;
	  }
	  BigReal rma = rmass[a];
	  BigReal rmb = rmass[b];
	  BigReal gab = diffsq / ( 2.0 * ( rma + rmb ) * rpab );
	  Vector dp = rab * gab;
	  pos[a] += rma * dp;
	  pos[b] -= rmb * dp;
	  if ( dt != 0. ) {
	    dp /= dt;
	    vel[a] += rma * dp;
	    vel[b] -= rmb * dp;
	  }
	  done = 0;
	}
      }
      if ( done ) break;
    }
    if ( consFailure ) {
      if ( dieOnError ) {
	iout << iERROR << "Constraint failure in RATTLE algorithm for atom "
			<< (atom[ig].id + 1) << "!\n" << endi;
	return -1;  // triggers early exit
      } else {
	iout << iWARN << "Constraint failure in RATTLE algorithm for atom "
			<< (atom[ig].id + 1) << "!\n" << endi;
      }
    } else if ( ! done ) {
      if ( dieOnError ) {
	iout << iERROR << "Exceeded RATTLE iteration limit for atom "
			<< (atom[ig].id + 1) << "!\n" << endi;
	return -1;  // triggers early exit
      } else {
	iout << iWARN << "Exceeded RATTLE iteration limit for atom "
			<< (atom[ig].id + 1) << "!\n" << endi;
      }
    }
    // store data back to patch
    for ( i = 0; i < hgs; ++i ) {
      atom[ig+i].position = pos[i];
      atom[ig+i].velocity = vel[i];
    }
  }
  return 0;

}

//  RATTLE algorithm from Allen & Tildesley
void HomePatch::rattle2(const BigReal timestep, Tensor *virial)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal dt = timestep / TIMEFACTOR;
  Tensor wc;  // constraint virial
  BigReal tol = simParams->rigidTol;
  int maxiter = simParams->rigidIter;
  int dieOnError = simParams->rigidDie;
  int i, iter;
  BigReal dsqi[10], tmp;
  int ial[10], ibl[10];
  Vector ref[10];  // reference position
  Vector refab[10];  // reference vector
  Vector vel[10];  // new velocity
  BigReal rmass[10];  // 1 / mass
  BigReal redmass[10];  // reduced mass
  int fixed[10];  // is atom fixed?

  //  CkPrintf("In rattle2!\n");
  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    //    CkPrintf("ig=%d\n",ig);
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = atom[ig+i].position;
      vel[i] = atom[ig+i].velocity;
      rmass[i] = 1. / atom[ig+i].mass;
      fixed[i] = ( atom[ig+i].atomFixed );
      if ( fixed[i] ) rmass[i] = 0.;
    }
    int icnt = 0;
    if ( ( tmp = mol->rigid_bond_length(atom[ig].id) ) > 0 ) {  // for water
      if ( hgs != 3 ) {
        NAMD_bug("Hydrogen group error caught in rattle2().");
      }
      if ( !(fixed[1] && fixed[2]) ) {
	redmass[icnt] = 1. / (rmass[1] + rmass[2]);
	dsqi[icnt] = 1. / (tmp * tmp);  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
      }
    }
    //    CkPrintf("Loop 2\n");
    for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = mol->rigid_bond_length(atom[ig+i].id) ) > 0 ) {
        if ( !(fixed[0] && fixed[i]) ) {
	  redmass[icnt] = 1. / (rmass[0] + rmass[i]);
	  dsqi[icnt] = 1. / (tmp * tmp);  ial[icnt] = 0;
	  ibl[icnt] = i;  ++icnt;
	}
      }
    }
    if ( icnt == 0 ) continue;  // no constraints
    //    CkPrintf("Loop 3\n");
    for ( i = 0; i < icnt; ++i ) {
      refab[i] = ref[ial[i]] - ref[ibl[i]];
    }
    //    CkPrintf("Loop 4\n");
    int done;
    for ( iter = 0; iter < maxiter; ++iter ) {
      done = 1;
      for ( i = 0; i < icnt; ++i ) {
	int a = ial[i];  int b = ibl[i];
	Vector vab = vel[a] - vel[b];
	Vector &rab = refab[i];
	BigReal rabsqi = dsqi[i];
	BigReal rvab = rab.x*vab.x + rab.y*vab.y + rab.z*vab.z;
	if ( (fabs(rvab) * dt * rabsqi) > tol ) {
	  Vector dp = rab * (-rvab * redmass[i] * rabsqi);
	  wc += outer(dp,rab);
	  vel[a] += rmass[a] * dp;
	  vel[b] -= rmass[b] * dp;
	  done = 0;
	}
      }
      if ( done ) break;
    }
    if ( ! done ) {
      if ( dieOnError ) {
	NAMD_die("Exceeded maximum number of iterations in rattle2().");
      } else {
	iout << iWARN <<
	  "Exceeded maximum number of iterations in rattle2().\n" << endi;
      }
    }
    // store data back to patch
    for ( i = 0; i < hgs; ++i ) {
      atom[ig+i].velocity = vel[i];
    }
  }
  //  CkPrintf("Leaving rattle2!\n");
  // check that there isn't a constant needed here!
  *virial += wc / ( 0.5 * dt );

}


//  MOLLY algorithm part 1
void HomePatch::mollyAverage()
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  BigReal tol = simParams->mollyTol;
  int maxiter = simParams->mollyIter;
  int i, iter;
  HGArrayBigReal dsq;
  BigReal tmp;
  HGArrayInt ial, ibl;
  HGArrayVector ref;  // reference position
  HGArrayVector refab;  // reference vector
  HGArrayBigReal rmass;  // 1 / mass
  BigReal *lambda;  // Lagrange multipliers
  CompAtom *avg;  // averaged position
  int numLambdas = 0;
  HGArrayInt fixed;  // is atom fixed?

  //  iout<<iINFO << "mollyAverage: "<<endl<<endi;
  p_avg.resize(numAtoms);
  for ( i=0; i<numAtoms; ++i ) p_avg[i] = p[i];

  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
	for ( i = 0; i < hgs; ++i ) {
	  ref[i] = atom[ig+i].position;
	  rmass[i] = 1. / atom[ig+i].mass;
	  fixed[i] = ( atom[ig+i].atomFixed );
	  if ( fixed[i] ) rmass[i] = 0.;
	}
	avg = &(p_avg[ig]);
	int icnt = 0;

	if ( ( tmp = mol->rigid_bond_length(atom[ig].id) ) ) {  // for water
	  if ( hgs != 3 ) {
	    NAMD_die("Hydrogen group error caught in mollyAverage().  It's a bug!\n");
	  }
	  if ( !(fixed[1] && fixed[2]) ) {
	    dsq[icnt] = tmp * tmp;  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
	  }
	}
	for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
	  if ( ( tmp = mol->rigid_bond_length(atom[ig+i].id) ) ) {
	    if ( !(fixed[0] && fixed[i]) ) {
	      dsq[icnt] =  tmp * tmp;  ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
	    }
	  }
	}
	if ( icnt == 0 ) continue;  // no constraints
	numLambdas += icnt;
	molly_lambda.resize(numLambdas);
	lambda = &(molly_lambda[numLambdas - icnt]);
	for ( i = 0; i < icnt; ++i ) {
	  refab[i] = ref[ial[i]] - ref[ibl[i]];
	}
	//	iout<<iINFO<<"hgs="<<hgs<<" m="<<icnt<<endl<<endi;
	iter=average(avg,ref,lambda,hgs,icnt,rmass,dsq,ial,ibl,refab,tol,maxiter);
	if ( iter == maxiter ) {
	  iout << iWARN << "Exceeded maximum number of iterations in mollyAverage().\n"<<endi;
	}
  }
}


//  MOLLY algorithm part 2
void HomePatch::mollyMollify(Tensor *virial)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  Tensor wc;  // constraint virial
  int i;
  HGArrayInt ial, ibl;
  HGArrayVector ref;  // reference position
  CompAtom *avg;  // averaged position
  HGArrayVector refab;  // reference vector
  HGArrayVector force;  // new force
  HGArrayBigReal rmass;  // 1 / mass
  BigReal *lambda;  // Lagrange multipliers
  int numLambdas = 0;
  HGArrayInt fixed;  // is atom fixed?

  for ( int ig = 0; ig < numAtoms; ig += atom[ig].hydrogenGroupSize ) {
    int hgs = atom[ig].hydrogenGroupSize;
    if (hgs == 1 ) continue;  // only one atom in group
	for ( i = 0; i < hgs; ++i ) {
	  ref[i] = atom[ig+i].position;
	  force[i] = f[Results::slow][ig+i];
	  rmass[i] = 1. / atom[ig+i].mass;
	  fixed[i] = ( atom[ig+i].atomFixed );
	  if ( fixed[i] ) rmass[i] = 0.;
	}
	int icnt = 0;
	// c-ji I'm only going to mollify water for now
	if ( ( mol->rigid_bond_length(atom[ig].id) ) ) {  // for water
	  if ( hgs != 3 ) {
	    NAMD_die("Hydrogen group error caught in mollyMollify().  It's a bug!\n");
	  }
	  if ( !(fixed[1] && fixed[2]) ) {
	    ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
	  }
	}
	for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
	  if ( ( mol->rigid_bond_length(atom[ig+i].id) ) ) {
	    if ( !(fixed[0] && fixed[i]) ) {
	      ial[icnt] = 0;  ibl[icnt] = i;  ++icnt;
	    }
	  }
	}

	if ( icnt == 0 ) continue;  // no constraints
	lambda = &(molly_lambda[numLambdas]);
	numLambdas += icnt;
	for ( i = 0; i < icnt; ++i ) {
	  refab[i] = ref[ial[i]] - ref[ibl[i]];
	}
	avg = &(p_avg[ig]);
	mollify(avg,ref,lambda,force,hgs,icnt,rmass,ial,ibl,refab);
	// store data back to patch
	for ( i = 0; i < hgs; ++i ) {
	  wc += outer(force[i]-f[Results::slow][ig+i],ref[i]);
	  f[Results::slow][ig+i] = force[i];
	}
  }
  // check that there isn't a constant needed here!
  *virial += wc;
  p_avg.resize(0);
}

void HomePatch::checkpoint(void) {
  FullAtomList tmp_a(&atom); checkpoint_atom = tmp_a;
  checkpoint_lattice = lattice;
}

void HomePatch::revert(void) {
  FullAtomList tmp_a(&checkpoint_atom); atom = tmp_a;
  numAtoms = atom.size();
  lattice = checkpoint_lattice;
}

BigReal HomePatch::calcKineticEnergy()
{
  BigReal total = 0;
  for ( int i = 0; i < numAtoms; ++i )
    {
      total += 0.5 * atom[i].mass * atom[i].velocity * atom[i].velocity;
    }
  return total;
}

Vector HomePatch::calcMomentum()
{
  Vector total(0,0,0);
  for ( int i = 0; i < numAtoms; ++i )
  {
     total += atom[i].mass * atom[i].velocity;
  }
  return total;
}

Vector HomePatch::calcAngularMomentum()
{
  Vector total(0,0,0);
  for ( int i = 0; i < numAtoms; ++i )
  {
     total += cross(atom[i].mass,atom[i].position,atom[i].velocity); // m r % v
  }
  return total;
}

void HomePatch::submitLoadStats(int timestep)
{
  LdbCoordinator::Object()->patchLoad(patchID,numAtoms,timestep);
}

void HomePatch::doGroupSizeCheck()
{
  if ( ! flags.doNonbonded ) return;

  SimParameters *simParams = Node::Object()->simParameters;
  BigReal hgcut = 0.5 * simParams->hgroupCutoff;  hgcut *= hgcut;

  FullAtomList::iterator p_i = atom.begin();
  FullAtomList::iterator p_e = atom.end();

  while ( p_i != p_e ) {
    int hgs = p_i->hydrogenGroupSize;
    p_i->nonbondedGroupIsAtom = 0;
    BigReal x = p_i->position.x;
    BigReal y = p_i->position.y;
    BigReal z = p_i->position.z;
    ++p_i;
    int oversize = 0;
    // limit spatial extent
    for ( int i = 1; i < hgs; ++i ) {
      p_i->nonbondedGroupIsAtom = 0;
      BigReal dx = p_i->position.x - x;
      BigReal dy = p_i->position.y - y;
      BigReal dz = p_i->position.z - z;
      BigReal r2 = dx * dx + dy * dy + dz * dz;
      ++p_i;
      if ( r2 > hgcut ) oversize = 1;
    }
    // also limit to at most 4 atoms per group
    if ( oversize || hgs > 4 ) {
      p_i -= hgs;
      for ( int i = 0; i < hgs; ++i ) {
        p_i->nonbondedGroupIsAtom = 1;
        ++p_i;
      }
    }
  }
}

void HomePatch::doMarginCheck()
{
  SimParameters *simParams = Node::Object()->simParameters;

  BigReal sysdima = lattice.a_r().unit() * lattice.a();
  BigReal sysdimb = lattice.b_r().unit() * lattice.b();
  BigReal sysdimc = lattice.c_r().unit() * lattice.c();

  BigReal cutoff = simParams->cutoff;

  BigReal margina = 0.5 * ( max.x - min.x - cutoff / sysdima );
  BigReal marginb = 0.5 * ( max.y - min.y - cutoff / sysdimb );
  BigReal marginc = 0.5 * ( max.z - min.z - cutoff / sysdimc );

  BigReal minx = min.x - margina;
  BigReal miny = min.y - marginb;
  BigReal minz = min.z - marginc;
  BigReal maxx = max.x + margina;
  BigReal maxy = max.y + marginb;
  BigReal maxz = max.z + marginc;

  int xdev, ydev, zdev;
  int problemCount = 0;

  FullAtomList::iterator p_i = atom.begin();
  FullAtomList::iterator p_e = atom.end();
  for ( ; p_i != p_e; ++p_i ) {

    ScaledPosition s = lattice.scale(p_i->position);

    // check if atom is within bounds
    if (s.x < minx) xdev = 0;
    else if (maxx <= s.x) xdev = 2; 
    else xdev = 1;

    if (s.y < miny) ydev = 0;
    else if (maxy <= s.y) ydev = 2; 
    else ydev = 1;

    if (s.z < minz) zdev = 0;
    else if (maxz <= s.z) zdev = 2; 
    else zdev = 1;

    if (mInfo[xdev][ydev][zdev]) { // somewhere else to be
	++problemCount;
    }

  }

  if ( problemCount ) {
      iout << iERROR <<
	"Found " << problemCount << " margin violations!\n" << endi;
  } 

}


void
HomePatch::doAtomMigration()
{
  int i;

  for (i=0; i<numNeighbors; i++) {
    realInfo[i].mList.resize(0);
  }

  // Purge the AtomMap
  AtomMap::Object()->unregisterIDs(patchID,p.begin(),p.end());

  // Determine atoms that need to migrate

  BigReal minx = min.x;
  BigReal miny = min.y;
  BigReal minz = min.z;
  BigReal maxx = max.x;
  BigReal maxy = max.y;
  BigReal maxz = max.z;

  int xdev, ydev, zdev;
  int delnum = 0;

  FullAtomList::iterator atom_i = atom.begin();
  FullAtomList::iterator atom_e = atom.end();
  while ( atom_i != atom_e ) {
    if ( atom_i->hydrogenGroupSize ) {

      ScaledPosition s = lattice.scale(atom_i->position);

      // check if atom is within bounds
      if (s.x < minx) xdev = 0;
      else if (maxx <= s.x) xdev = 2;
      else xdev = 1;

      if (s.y < miny) ydev = 0;
      else if (maxy <= s.y) ydev = 2;
      else ydev = 1;

      if (s.z < minz) zdev = 0;
      else if (maxz <= s.z) zdev = 2;
      else zdev = 1;

    }

    if (mInfo[xdev][ydev][zdev]) { // process atom for migration
                                    // Don't migrate if destination is myself

      // See if we have a migration list already
      MigrationList &mCur = mInfo[xdev][ydev][zdev]->mList;
      DebugM(3,"Migrating atom " << atomIDList_i << " from patch "
		<< patchID << " with position " << p_i << "\n");
      mCur.add(*atom_i);

      ++delnum;

    } else {

      if ( delnum ) { *(atom_i-delnum) = *atom_i; }

    }

    ++atom_i;

  }

  int delpos = numAtoms - delnum;
  DebugM(4,"numAtoms " << numAtoms << " deleted " << delnum << "\n");
  atom.del(delpos,delnum);

  numAtoms = atom.size();

  PatchMgr::Object()->sendMigrationMsgs(patchID, realInfo, numNeighbors);

  // signal depositMigration() that we are inMigration mode
  inMigration = true;

  // Drain the migration message buffer
  for (i=0; i<numMlBuf; i++) {
     DebugM(1, "Draining migration buffer ("<<i<<","<<numMlBuf<<")\n");
     depositMigration(msgbuf[i]);
  }
  numMlBuf = 0;
     
  if (!allMigrationIn) {
    DebugM(3,"All Migrations NOT in, we are suspending patch "<<patchID<<"\n");
    migrationSuspended = true;
    sequencer->suspend();
    migrationSuspended = false;
  }
  allMigrationIn = false;

  inMigration = false;
}

void 
HomePatch::depositMigration(MigrateAtomsMsg *msg)
{

  if (!inMigration) { // We have to buffer changes due to migration
		      // until our patch is in migration mode
    msgbuf[numMlBuf++] = msg;
    return;
  } 

  {
    MigrationList &migrationList = msg->migrationList;
    MigrationList::iterator mi;
    for (mi = migrationList.begin(); mi != migrationList.end(); mi++) {
      DebugM(1,"Migrating atom " << mi->atomID << " to patch "
		<< patchID << " with position " << mi->pos << "\n"); 
      mi->position = lattice.nearest(mi->position,center,&(mi->transform));
      atom.add(*mi);
    }
  }

  numAtoms = atom.size();
  delete msg;

  DebugM(3,"Counter on " << patchID << " = " << patchMigrationCounter << "\n");
  if (!--patchMigrationCounter) {
    DebugM(3,"All Migrations are in for patch "<<patchID<<"\n");
    allMigrationIn = true;
    patchMigrationCounter = numNeighbors;
    if (migrationSuspended) {
      DebugM(3,"patch "<<patchID<<" is being awakened\n");
      migrationSuspended = false;
      sequencer->awaken();
      return;
    }
    else {
       DebugM(3,"patch "<<patchID<<" was not suspended\n");
    }
  }
}

inline void lubksb(HGMatrixBigReal &a, int n, HGArrayInt &indx,
                                              HGArrayBigReal &b)
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii >= 0)
			for (j=ii;j<i;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


inline void ludcmp(HGMatrixBigReal &a, int n, HGArrayInt &indx, BigReal *d)
{

	int i,imax,j,k;
	double big,dum,sum,temp;
	HGArrayBigReal vv;
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) NAMD_die("Singular matrix in routine ludcmp\n");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}


inline void G_q(const HGArrayVector &refab, HGMatrixVector &gqij,
     const int n, const int m, const HGArrayInt &ial, const HGArrayInt &ibl) {
  int i; 
  // step through the rows of the matrix
  for(i=0;i<m;i++) {
    gqij[i][ial[i]]=2.0*refab[i];
    gqij[i][ibl[i]]=-gqij[i][ial[i]];
  }
};


// c-ji code for MOLLY 7-31-99
int average(CompAtom *qtilde,const HGArrayVector &q,BigReal *lambda,const int n,const int m, const HGArrayBigReal &imass, const HGArrayBigReal &length2, const HGArrayInt &ial, const HGArrayInt &ibl, const HGArrayVector &refab, const BigReal tolf, const int ntrial) {
  //  input:  n = length of hydrogen group to be averaged (shaked)
  //          q[n] = vector array of original positions
  //          m = number of constraints
  //          imass[n] = inverse mass for each atom
  //          length2[m] = square of reference bond length for each constraint
  //          ial[m] = atom a in each constraint 
  //          ibl[m] = atom b in each constraint 
  //          refab[m] = vector of q_ial(i) - q_ibl(i) for each constraint
  //          tolf = function error tolerance for Newton's iteration
  //          ntrial = max number of Newton's iterations
  //  output: lambda[m] = double array of lagrange multipliers (used by mollify)
  //          qtilde[n] = vector array of averaged (shaked) positions

  int k,k1,i,j;
  BigReal errx,errf,d,tolx;

  HGArrayInt indx;
  HGArrayBigReal p;
  HGArrayBigReal fvec;
  HGMatrixBigReal fjac;
  HGArrayVector avgab;
  HGMatrixVector grhs;
  HGMatrixVector auxrhs;
  HGMatrixVector glhs;

  //  iout <<iINFO << "average: n="<<n<<" m="<<m<<endl<<endi;
  tolx=tolf; 
  
  // initialize lambda, globalGrhs

  for (i=0;i<m;i++) {
    lambda[i]=0.0;
  }

  // define grhs, auxrhs for all iterations
  // grhs= g_x(q)
  //
  G_q(refab,grhs,n,m,ial,ibl);
  for (k=1;k<=ntrial;k++) {
    //    usrfun(qtilde,q0,lambda,fvec,fjac,n,water); 
    HGArrayBigReal gij;
    // this used to be main loop of usrfun
    // compute qtilde given q0, lambda, IMASSes
    {
      BigReal multiplier;
      HGArrayVector tmp;
      for (i=0;i<m;i++) {
	multiplier = lambda[i];
	// auxrhs = M^{-1}grhs^{T}
	for (j=0;j<n;j++) {
	  auxrhs[i][j]=multiplier*imass[j]*grhs[i][j];
	}
      }
      for (j=0;j<n;j++) {
	//      tmp[j]=0.0;      
	for (i=0;i<m;i++) {
	  tmp[j]+=auxrhs[i][j];
	}
      }
 
      for (j=0;j<n;j++) {
	qtilde[j].position=q[j]+tmp[j];
      }
      //      delete [] tmp;
    }
  
    for ( i = 0; i < m; i++ ) {
      avgab[i] = qtilde[ial[i]].position - qtilde[ibl[i]].position;
    }

    //  iout<<iINFO << "Calling Jac" << endl<<endi;
    //  Jac(qtilde, q0, fjac,n,water);
    {
      //  Vector glhs[3*n+3];

      HGMatrixVector grhs2;

      G_q(avgab,glhs,n,m,ial,ibl);
#ifdef DEBUG0
      iout<<iINFO << "G_q:" << endl<<endi;
      for (i=0;i<m;i++) {
	iout<<iINFO << glhs[i*n+0] << " " << glhs[i*n+1] << " " << glhs[i*n+2] << endl<<endi;
      }
#endif
      //      G_q(refab,grhs2,m,ial,ibl);
      // update with the masses
      for (j=0; j<n; j++) { // number of atoms
	for (i=0; i<m;i++) { // number of constraints
	  grhs2[i][j] = grhs[i][j]*imass[j];
	}
      }

      // G_q(qtilde) * M^-1 G_q'(q0) =
      // G_q(qtilde) * grhs'
      for (i=0;i<m;i++) { // number of constraints
	for (j=0;j<m;j++) { // number of constraints
	  fjac[i][j] = 0; 
	  for (k1=0;k1<n;k1++) {
	    fjac[i][j] += glhs[i][k1]*grhs2[j][k1]; 
	  }
	}
      }
#ifdef DEBUG0  
      iout<<iINFO << "glhs" <<endi;
      for(i=0;i<9;i++) {
	iout<<iINFO << glhs[i] << ","<<endi;
      }
      iout<<iINFO << endl<<endi;
      for(i=0;i<9;i++) {
	iout<<iINFO << grhs2[i] << ","<<endi;
      }
      iout<<iINFO << endl<<endi;
#endif
      //      delete[] grhs2;
    }
    // end of Jac calculation
#ifdef DEBUG0
    iout<<iINFO << "Jac" << endl<<endi;
    for (i=0;i<m;i++) 
      for (j=0;j<m;j++)
	iout<<iINFO << fjac[i][j] << " "<<endi;
    iout<< endl<<endi;
#endif
    // calculate constraints in gij for n constraints this being a water
    //  G(qtilde, gij, n, water);
    for (i=0;i<m;i++) {
      gij[i]=avgab[i]*avgab[i]-length2[i];
    }
#ifdef DEBUG0
    iout<<iINFO << "G" << endl<<endi;
    iout<<iINFO << "( "<<endi;
    for(i=0;i<m-1;i++) {
      iout<<iINFO << gij[i] << ", "<<endi;
    }
    iout<<iINFO << gij[m-1] << ")" << endl<<endi;
#endif
    // fill the return vector
    for(i=0;i<m;i++) {
      fvec[i] = gij[i];
    }
    // free up the constraints
    //    delete[] gij;
    // continue Newton's iteration    
    errf=0.0;
    for (i=0;i<m;i++) errf += fabs(fvec[i]);
#ifdef DEBUG0
    iout<<iINFO << "errf: " << errf << endl<<endi;
#endif
    if (errf <= tolf) {
      break;
    }
    for (i=0;i<m;i++) p[i] = -fvec[i];
    //    iout<<iINFO << "Doing dcmp in average " << endl<<endi;
    ludcmp(fjac,m,indx,&d);
    lubksb(fjac,m,indx,p);

    errx=0.0;
    for (i=0;i<m;i++) {
      errx += fabs(p[i]);
    }
    for (i=0;i<m;i++)  
      lambda[i] += p[i];

#ifdef DEBUG0
    iout<<iINFO << "lambda:" << lambda[0] 
	 << " " << lambda[1] << " " << lambda[2] << endl<<endi;
    iout<<iINFO << "errx: " << errx << endl<<endi;
#endif
    if (errx <= tolx) break;
#ifdef DEBUG0
    iout<<iINFO << "Qtilde:" << endl<<endi;
    iout<<iINFO << qtilde[0].position << " " << qtilde[1].position << " " << qtilde[2].position << endl<<endi; 
#endif
  }
#ifdef DEBUG
  iout<<iINFO << "LAMBDA:" << lambda[0] << " " << lambda[1] << " " << lambda[2] << endl<<endi;
#endif

  return k; // 
}

void mollify(CompAtom *qtilde,const HGArrayVector &q0,const BigReal *lambda, HGArrayVector &force,const int n, const int m, const HGArrayBigReal &imass,const HGArrayInt &ial,const HGArrayInt &ibl,const HGArrayVector &refab) {
  int i,j,k;
  BigReal d;
  HGMatrixBigReal fjac;
  Vector zero(0.0,0.0,0.0);
  
  HGArrayVector tmpforce;
  HGArrayVector tmpforce2;
  HGArrayVector y;
  HGMatrixVector grhs;
  HGMatrixVector glhs;
  HGArrayBigReal aux;
  HGArrayInt indx;

  for(i=0;i<n;i++) {
    tmpforce[i]=imass[i]*force[i];
  }

  HGMatrixVector grhs2;
  HGArrayVector avgab;

  for ( i = 0; i < m; i++ ) {
	avgab[i] = qtilde[ial[i]].position - qtilde[ibl[i]].position;
  }

  G_q(avgab,glhs,n,m,ial,ibl);
  G_q(refab,grhs,n,m,ial,ibl);
  // update with the masses
  for (j=0; j<n; j++) { // number of atoms
	for (i=0; i<m;i++) { // number of constraints
	  grhs2[i][j] = grhs[i][j]*imass[j];
	}
  }

  // G_q(qtilde) * M^-1 G_q'(q0) =
  // G_q(qtilde) * grhs'
  for (i=0;i<m;i++) { // number of constraints
	for (j=0;j<m;j++) { // number of constraints
	  fjac[j][i] = 0; 
	  for (k=0;k<n;k++) {
	    fjac[j][i] += glhs[i][k]*grhs2[j][k]; 
	  }
	}
  }

  // aux=gqij*tmpforce
  //  globalGrhs::computeGlobalGrhs(q0,n,water);
  //  G_q(refab,grhs,m,ial,ibl);
  for(i=0;i<m;i++) {
    aux[i]=0.0;
    for(j=0;j<n;j++) {
      aux[i]+=grhs[i][j]*tmpforce[j];
    }
  }

  ludcmp(fjac,m,indx,&d);
  lubksb(fjac,m,indx,aux);

  for(j=0;j<n;j++) {
    y[j] = zero;
    for(i=0;i<m;i++) {
      y[j] += aux[i]*glhs[i][j];
    }
  }
  for(i=0;i<n;i++) {
    y[i]=force[i]-y[i];
  }
    
  // gqq12*y
  for(i=0;i<n;i++) {
    tmpforce2[i]=imass[i]*y[i];
  }

  // here we assume that tmpforce is initialized to zero.
  for (i=0;i<n;i++) {
    tmpforce[i]=zero;
  }
  
  for (j=0;j<m;j++) {
    Vector tmpf = 2.0*lambda[j]*(tmpforce2[ial[j]]-tmpforce2[ibl[j]]);
    tmpforce[ial[j]] += tmpf;
    tmpforce[ibl[j]] -= tmpf;
  }
  // c-ji the other bug for 2 constraint water was this line (2-4-99)
  //  for(i=0;i<m;i++) {
  for(i=0;i<n;i++) {
    force[i]=tmpforce[i]+y[i];
  }

}

