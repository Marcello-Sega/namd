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

int average(Vector *qtilde,const HGArrayVector &q,BigReal *lambda,const int n,const int m, const HGArrayBigReal &imass, const HGArrayBigReal &length2, const HGArrayInt &ial, const HGArrayInt &ibl, const HGArrayVector &refab, const BigReal tolf, const int ntrial);

void mollify(Vector *qtilde,const HGArrayVector &q0,const BigReal *lambda, HGArrayVector &force,const int n, const int m, const HGArrayBigReal &imass,const HGArrayInt &ial,const HGArrayInt &ibl,const HGArrayVector &refab);


HomePatch::HomePatch(PatchID pd, AtomIDList al, TransformList tl,
      PositionList pl, VelocityList vl) : Patch(pd,al,pl), v(vl), t(tl),
      p_checkpoint(&pl)
{ 
  DebugM(4, "HomePatch("<<pd<<") at " << this << "\n");
  if (atomIDList.size() != v.size()) {
    CkPrintf("HomePatch::HomePatch(...) : size mismatch-Velocities and IDs!\n");
  }
  AtomMap::Object()->registerIDs(pd,al);  
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
}

// Bind a Sequencer to this HomePatch
void HomePatch::useSequencer(Sequencer *sequencerPtr)
{ sequencer=sequencerPtr; }

// start simulation over this Patch of atoms
void HomePatch::runSequencer(void)
{ sequencer->run(); }

void HomePatch::readPatchMap() {
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
  ProxyAtomsMsg *nmsg = new ProxyAtomsMsg;
  nmsg->patch = patchID;
  nmsg->atomIDList = AtomIDList(&atomIDList); // copy the array
  ProxyMgr::Object()->sendProxyAtoms(nmsg,msg->node);
  delete msg;
}

void HomePatch::unregisterProxy(UnregisterProxyMsg *msg) {
  int i = proxy.findIndex(ProxyListElem(msg->node));
  forceBox.checkIn(proxy[i].forceBox);
  proxy.del(i);
  delete msg;
}

void HomePatch::receiveResults(ProxyResultMsg *msg)
{
  DebugM(4, "patchID("<<patchID<<") receiveRes() nodeID("<<msg->node<<")\n");
  int i = proxy.findIndex(ProxyListElem(msg->node));
  Results *r = proxy[i].forceBox->open();
  for ( int k = 0; k < Results::maxNumForces; ++k )
  {
    Force *f = r->f[k];
    register ForceList::iterator f_i, f_e;
    f_i = msg->forceList[k].begin();
    f_e = msg->forceList[k].end();
    for ( ; f_i != f_e; ++f_i, ++f ) *f += *f_i;
  }
  proxy[i].forceBox->close(&r);
  delete msg;
}


void HomePatch::positionsReady(int doMigration)
{
  if (!patchMapRead) {
    readPatchMap();
    patchMapRead = 1;
  }
      
  doMigration = doMigration && numNeighbors;

  if (doMigration) {
    doAtomMigration();
  } else {
    doMarginCheck();
  }

  // Workaround for oversize groups now handled by Patch.
  // doGroupSizeCheck();

  if (flags.doMolly) mollyAverage();

  // Must Add Proxy Changes when migration completed!
  ProxyListIter pli(proxy);
  for ( pli = pli.begin(); pli != pli.end(); ++pli )
  {
    if (doMigration) {
      ProxyAllMsg *allmsg = new ProxyAllMsg;
      allmsg->patch = patchID;
      allmsg->flags = flags;
      allmsg->positionList = p;
      if (flags.doMolly) allmsg->avgPositionList = p_avg;
      allmsg->atomIDList = atomIDList;
      DebugM(1, "atomIDList.size() = " << atomIDList.size() << " p.size() = " << p.size() << "\n" );
      ProxyMgr::Object()->sendProxyAll(allmsg,pli->node);
    } else {
      ProxyDataMsg *nmsg = new ProxyDataMsg;
      nmsg->patch = patchID;
      nmsg->flags = flags;
      nmsg->positionList = p;
      if (flags.doMolly) nmsg->avgPositionList = p_avg;
      ProxyMgr::Object()->sendProxyData(nmsg,pli->node);
    }   
  }
  DebugM(4, "patchID("<<patchID<<") doing positions Ready\n");
  Patch::positionsReady(doMigration);

}


void HomePatch::addForceToMomentum(const BigReal timestep, const int ftag)
{
  const BigReal dt = timestep / TIMEFACTOR;
  //if (v.check() == NULL || f.check() == NULL || a.check() == NULL) {
      //DebugM(5, "NULL found! v="<<v<<" f="<<f<<" a="<<a<<"\n");
  //}
  for ( int i = 0; i < numAtoms; ++i )
  {
    v[i] += f[ftag][i] * ( dt / a[i].mass );
    if ( a[i].flags & ATOM_FIXED ) v[i] = 0;
  }
}

void HomePatch::addVelocityToPosition(const BigReal timestep)
{
  const BigReal dt = timestep / TIMEFACTOR;
  for ( int i = 0; i < numAtoms; ++i )
  {
    if ( ! ( a[i].flags & ATOM_FIXED ) ) p[i] += v[i] * dt;
  }
}

//  RATTLE algorithm from Allen & Tildesley
void HomePatch::rattle1(const BigReal timestep)
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
  
  for ( int ig = 0; ig < numAtoms; ig += a[ig].hydrogenGroupSize ) {
    int hgs = a[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = p[ig+i];
      pos[i] = p[ig+i];
      vel[i] = v[ig+i];
      rmass[i] = 1. / a[ig+i].mass;
      fixed[i] = ( a[ig+i].flags & ATOM_FIXED );
      // undo addVelocityToPosition to get proper reference coordinates
      if ( fixed[i] ) rmass[i] = 0.; else ref[i] -= vel[i] * dt;
    }
    int icnt = 0;
    if ( ( tmp = mol->rigid_bond_length(a[ig].id) ) > 0 ) {  // for water
      if ( hgs != 3 ) {
        NAMD_bug("Hydrogen group error caught in rattle1().");
      }
      if ( !(fixed[1] && fixed[2]) ) {
	dsq[icnt] = tmp * tmp;  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
      }
    }
    for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
      if ( ( tmp = mol->rigid_bond_length(a[ig+i].id) ) > 0 ) {
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
	NAMD_die("Constraint failure in RATTLE algorithm!");
      } else {
	iout << iWARN <<
	  "Constraint failure in RATTLE algorithm!\n" << endi;
      }
    } else if ( ! done ) {
      if ( dieOnError ) {
	NAMD_die("Exceeded maximum number of iterations in rattle1().");
      } else {
	iout << iWARN <<
	  "Exceeded maximum number of iterations in rattle1().\n" << endi;
      }
    }
    // store data back to patch
    for ( i = 0; i < hgs; ++i ) {
      p[ig+i] = pos[i];
      v[ig+i] = vel[i];
    }
  }

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
  for ( int ig = 0; ig < numAtoms; ig += a[ig].hydrogenGroupSize ) {
    //    CkPrintf("ig=%d\n",ig);
    int hgs = a[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
    // cache data in local arrays and integrate positions normally
    for ( i = 0; i < hgs; ++i ) {
      ref[i] = p[ig+i];
      vel[i] = v[ig+i];
      rmass[i] = 1. / a[ig+i].mass;
      fixed[i] = ( a[ig+i].flags & ATOM_FIXED );
      if ( fixed[i] ) rmass[i] = 0.;
    }
    int icnt = 0;
    if ( ( tmp = mol->rigid_bond_length(a[ig].id) ) > 0 ) {  // for water
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
      if ( ( tmp = mol->rigid_bond_length(a[ig+i].id) ) > 0 ) {
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
      v[ig+i] = vel[i];
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
  Vector *avg;  // averaged position
  int numLambdas = 0;
  HGArrayInt fixed;  // is atom fixed?

  //  iout<<iINFO << "mollyAverage: "<<endl<<endi;
  p_avg.resize(numAtoms);
  for ( i=0; i<numAtoms; ++i ) p_avg[i] = p[i];

  for ( int ig = 0; ig < numAtoms; ig += a[ig].hydrogenGroupSize ) {
    int hgs = a[ig].hydrogenGroupSize;
    if ( hgs == 1 ) continue;  // only one atom in group
	for ( i = 0; i < hgs; ++i ) {
	  ref[i] = p[ig+i];
	  rmass[i] = 1. / a[ig+i].mass;
	  fixed[i] = ( a[ig+i].flags & ATOM_FIXED );
	  if ( fixed[i] ) rmass[i] = 0.;
	}
	avg = &(p_avg[ig]);
	int icnt = 0;

	if ( ( tmp = mol->rigid_bond_length(a[ig].id) ) ) {  // for water
	  if ( hgs != 3 ) {
	    NAMD_die("Hydrogen group error caught in mollyAverage().  It's a bug!\n");
	  }
	  if ( !(fixed[1] && fixed[2]) ) {
	    dsq[icnt] = tmp * tmp;  ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
	  }
	}
	for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
	  if ( ( tmp = mol->rigid_bond_length(a[ig+i].id) ) ) {
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
  Vector *avg;  // averaged position
  HGArrayVector refab;  // reference vector
  HGArrayVector force;  // new force
  HGArrayBigReal rmass;  // 1 / mass
  BigReal *lambda;  // Lagrange multipliers
  int numLambdas = 0;
  HGArrayInt fixed;  // is atom fixed?

  for ( int ig = 0; ig < numAtoms; ig += a[ig].hydrogenGroupSize ) {
    int hgs = a[ig].hydrogenGroupSize;
    if (hgs == 1 ) continue;  // only one atom in group
	for ( i = 0; i < hgs; ++i ) {
	  ref[i] = p[ig+i];
	  force[i] = f[Results::slow][ig+i];
	  rmass[i] = 1. / a[ig+i].mass;
	  fixed[i] = ( a[ig+i].flags & ATOM_FIXED );
	  if ( fixed[i] ) rmass[i] = 0.;
	}
	int icnt = 0;
	// c-ji I'm only going to mollify water for now
	if ( ( mol->rigid_bond_length(a[ig].id) ) ) {  // for water
	  if ( hgs != 3 ) {
	    NAMD_die("Hydrogen group error caught in mollyMollify().  It's a bug!\n");
	  }
	  if ( !(fixed[1] && fixed[2]) ) {
	    ial[icnt] = 1;  ibl[icnt] = 2;  ++icnt;
	  }
	}
	for ( i = 1; i < hgs; ++i ) {  // normal bonds to mother atom
	  if ( ( mol->rigid_bond_length(a[ig+i].id) ) ) {
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
  for ( int i = 0; i < numAtoms; ++i ) {
    p_checkpoint[i] = p[i];
  }
}

void HomePatch::revert(void) {
  for ( int i = 0; i < numAtoms; ++i ) {
    p[i] = p_checkpoint[i];
  }
}

BigReal HomePatch::calcKineticEnergy()
{
  BigReal total = 0;
  for ( int i = 0; i < numAtoms; ++i )
    {
      total += 0.5 * a[i].mass * v[i] * v[i];
    }
  return total;
}

Vector HomePatch::calcMomentum()
{
  Vector total;
  for ( int i = 0; i < numAtoms; ++i )
  {
     total += a[i].mass * v[i];
  }
  return total;
}

Vector HomePatch::calcAngularMomentum()
{
  Vector total;
  for ( int i = 0; i < numAtoms; ++i )
  {
     total += cross(a[i].mass,p[i],v[i]); // m r % v
  }
  return total;
}

void HomePatch::submitLoadStats(int timestep)
{
  LdbCoordinator::Object()->patchLoad(patchID,numAtoms,timestep);
}


void HomePatch::doGroupSizeCheck()
{
  SimParameters *simParams = Node::Object()->simParameters;
  if ( simParams->splitPatch != SPLIT_PATCH_HYDROGEN ) return;
  BigReal hgcut = 0.5 * simParams->hgroupCutoff;  hgcut *= hgcut;

  AtomPropertiesList::iterator a_i = a.begin();
  PositionList::iterator p_i = p.begin();
  PositionList::iterator p_e = p.end();
  BigReal x=0,y=0,z=0;
  int parentID=0;

  while ( p_i != p_e ) {
    if ( a_i->hydrogenGroupSize ) {  // group parent
      x = p_i->x; y = p_i->y; z = p_i->z; parentID = a_i->id;
    } else {  // child atom
      BigReal dx = p_i->x - x;
      BigReal dy = p_i->y - y;
      BigReal dz = p_i->z - z;
      BigReal r2 = dx * dx + dy * dy + dz * dz;
      if ( r2 > hgcut ) {
        iout << iERROR << "H-group size violation between atoms " <<
		parentID + 1 << " and " << a_i->id + 1 <<
		" at distance " << sqrt(r2) << "!\n" <<
		"You may wish to increase the hgroupCutoff parameter " <<
		"from its current value of " << simParams->hgroupCutoff <<
		".\n" << endi;
      }
    }
    ++p_i; ++ a_i;
  }
}


void HomePatch::doMarginCheck()
{
  SimParameters *simParams = Node::Object()->simParameters;
  Position Min = lattice.unscale(min);
  Position Max = lattice.unscale(max);
  BigReal cutoff = simParams->cutoff;
  BigReal marginx = 0.5 * ( Max.x - Min.x - cutoff );
  BigReal marginy = 0.5 * ( Max.y - Min.y - cutoff );
  BigReal marginz = 0.5 * ( Max.z - Min.z - cutoff );
  BigReal minx = Min.x - marginx;
  BigReal miny = Min.y - marginy;
  BigReal minz = Min.z - marginz;
  BigReal maxx = Max.x + marginx;
  BigReal maxy = Max.y + marginy;
  BigReal maxz = Max.z + marginz;
  int xdev, ydev, zdev;
  int problemCount = 0;

  PositionList::iterator p_i = p.begin();
  PositionList::iterator p_e = p.end();
  for ( ; p_i != p_e; ++p_i ) {

    // check if atom should is within bounds
    if (p_i->x < minx) xdev = 0;
    else if (maxx <= p_i->x) xdev = 2; 
    else xdev = 1;

    if (p_i->y < miny) ydev = 0;
    else if (maxy <= p_i->y) ydev = 2; 
    else ydev = 1;

    if (p_i->z < minz) zdev = 0;
    else if (maxz <= p_i->z) zdev = 2; 
    else zdev = 1;

    if (mInfo[xdev][ydev][zdev]) { // process atom for migration
                                   // Don't migrate if destination is myself
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
  int i,j;
  int xdev=1, ydev=1, zdev=1;

  // Drain the migration message buffer
  //for (i=0; i<numMlBuf; i++) {
  //   DebugM(3, "Draining migration buffer ("<<i<<","<<numMlBuf<<")\n");
  //   depositMigration(srcID[i], mlBuf[i]);
  //}
  //numMlBuf = 0;
     
  // realInfo points to migration lists for neighbors we actually have. 
  //    element of mInfo[3][3][3] points to an element of realInfo
  for (i=0; i<numNeighbors; i++) {
    realInfo[i].mList.resize(0);
  }

  // Purge the AtomMap
  AtomMap::Object()->unregisterIDs(patchID,atomIDList);

  // Determine atoms that need to migrate
  AtomIDList::iterator atomIDList_i = atomIDList.begin();
  AtomIDList::iterator atomIDList_e = atomIDList.end();
  AtomPropertiesList::iterator a_i = a.begin();
  TransformList::iterator t_i = t.begin();
  PositionList::iterator p_i = p.begin();
  PositionList::iterator p_c_i = p_checkpoint.begin();
  VelocityList::iterator v_i = v.begin();
  ForceList::iterator f_i[Results::maxNumForces];
  for ( j = 0; j < Results::maxNumForces; ++j ) f_i[j] = f[j].begin();
  int delnum = 0;
  Position Min = lattice.unscale(min);
  Position Max = lattice.unscale(max);

  while ( atomIDList_i != atomIDList_e )
  {
      DebugM(4, atomIDList_e - atomIDList_i << " iterations remaining\n");

      if ( a_i->hydrogenGroupSize )
	  {
	  // check if atom should is within bounds
	  if (p_i->x < Min.x) xdev = 0;
	  else if (Max.x <= p_i->x) xdev = 2; 
	  else xdev = 1;

	  if (p_i->y < Min.y) ydev = 0;
	  else if (Max.y <= p_i->y) ydev = 2; 
	  else ydev = 1;

	  if (p_i->z < Min.z) zdev = 0;
	  else if (Max.z <= p_i->z) zdev = 2; 
	  else zdev = 1;
	  }

     if (mInfo[xdev][ydev][zdev]) { // process atom for migration
                                    // Don't migrate if destination is myself

       // See if we have a migration list already
       MigrationList &mCur = mInfo[xdev][ydev][zdev]->mList;
       Force force[Results::maxNumForces];
       for ( j = 0; j < Results::maxNumForces; ++j ) force[j] = *(f_i[j]);
       DebugM(3,"Migrating atom " << atomIDList_i << " from patch "
		<< patchID << " with position " << p_i << "\n");
       mCur.add(MigrationElem(*atomIDList_i, *a_i, *t_i, *p_i, *p_c_i, *v_i, force));

       ++delnum;

     } else {

       if ( delnum ) {
         *(atomIDList_i-delnum) = *atomIDList_i;
         *(a_i-delnum) = *a_i;
         *(t_i-delnum) = *t_i;
         *(p_i-delnum) = *p_i;
         *(p_c_i-delnum) = *p_c_i;
         *(v_i-delnum) = *v_i;
         for ( j = 0; j < Results::maxNumForces; ++j ) {
           *(f_i[j]-delnum) = *(f_i[j]);
         }
       }

     }

     ++atomIDList_i;
     ++a_i;
     ++t_i;
     ++p_i;
     ++p_c_i;
     ++v_i;
     for ( j = 0; j < Results::maxNumForces; ++j ) ++(f_i[j]);

  }

  int delpos = numAtoms - delnum;
  DebugM(4,"numAtoms " << numAtoms << " deleted " << delnum << "\n");
  atomIDList.del(delpos,delnum);
  a.del(delpos,delnum);
  t.del(delpos,delnum);
  p.del(delpos,delnum);
  p_checkpoint.del(delpos,delnum);
  v.del(delpos,delnum);
  for ( j = 0; j < Results::maxNumForces; ++j ) f[j].del(delpos,delnum);

  numAtoms = atomIDList.size();

  PatchMgr::Object()->sendMigrationMsgs(patchID, realInfo, numNeighbors);
  // for (i=0; i < numNeighbors; i++) {
  //   PatchMgr::Object()->sendMigrationMsg(patchID, realInfo[i]);
  // }

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
  // indexAtoms();  NEVER USED -JCP

  // reload the AtomMap
  AtomMap::Object()->registerIDs(patchID,atomIDList);

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
      atomIDList.add(mi->atomID);
      a.add(mi->atomProp);
      p.add(lattice.nearest(mi->pos,center,&(mi->trans)));
      // JCP FIX THIS!!!
      p_checkpoint.add(mi->pos_checkpoint);
      t.add(mi->trans);
      v.add(mi->vel);
      for ( int j = 0; j < Results::maxNumForces; ++j )
        f[j].add(mi-> force[j]);
    }
  }

  numAtoms = atomIDList.size();
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
int average(Vector *qtilde,const HGArrayVector &q,BigReal *lambda,const int n,const int m, const HGArrayBigReal &imass, const HGArrayBigReal &length2, const HGArrayInt &ial, const HGArrayInt &ibl, const HGArrayVector &refab, const BigReal tolf, const int ntrial) {
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
	qtilde[j]=q[j]+tmp[j];
      }
      //      delete [] tmp;
    }
  
    for ( i = 0; i < m; i++ ) {
      avgab[i] = qtilde[ial[i]] - qtilde[ibl[i]];
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
    iout<<iINFO << qtilde[0] << " " << qtilde[1] << " " << qtilde[2] << endl<<endi; 
#endif
  }
#ifdef DEBUG
  iout<<iINFO << "LAMBDA:" << lambda[0] << " " << lambda[1] << " " << lambda[2] << endl<<endi;
#endif

  return k; // 
}

void mollify(Vector *qtilde,const HGArrayVector &q0,const BigReal *lambda, HGArrayVector &force,const int n, const int m, const HGArrayBigReal &imass,const HGArrayInt &ial,const HGArrayInt &ibl,const HGArrayVector &refab) {
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
	avgab[i] = qtilde[ial[i]] - qtilde[ibl[i]];
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

