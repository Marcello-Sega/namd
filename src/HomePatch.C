/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: HomePatch owns the actual atoms of a Patch of space
 *		Proxy(s) get messages via ProxyMgr from HomePatch(es)
 *		to update lists of atoms and their coordinates
 *              HomePatch(es) also have a Sequencer bound to them
 *
 * superclass: 	Patch		
 ***************************************************************************/

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


// avoid dissappearence of ident?
char HomePatch::ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/HomePatch.C,v 1.1059 1999/08/31 15:43:31 jim Exp $";

HomePatch::HomePatch(PatchID pd, AtomIDList al, TransformList tl,
      PositionList pl, VelocityList vl) : Patch(pd,al,pl), v(vl), t(tl)
{ 
  DebugM(4, "HomePatch("<<pd<<") at " << this << "\n");
  if (atomIDList.size() != v.size()) {
    CkPrintf("HomePatch::HomePatch(...) : size mismatch-Velocities and IDs!\n");
  }
  AtomMap::Object()->registerIDs(pd,al);  
  min.x = PatchMap::Object()->minX(patchID);
  min.y = PatchMap::Object()->minY(patchID);
  min.z = PatchMap::Object()->minZ(patchID);
  max.x = PatchMap::Object()->maxX(patchID);
  max.y = PatchMap::Object()->maxY(patchID);
  max.z = PatchMap::Object()->maxZ(patchID);
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
    realInfo[n].mList = NULL;
  }

  // Make mapping from the 3x3x3 cube of pointers to real migration info
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
      {
	int pid =  p->pid(p->xIndex(patchID)+i-1, 
	    p->yIndex(patchID)+j-1, p->zIndex(patchID)+k-1);
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
  nmsg->atomIDList = atomIDList;
  nmsg->prepack();
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
        NAMD_die("Hydrogen group error caught in rattle1().  It's a bug!\n");
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
    for ( iter = 0; iter < maxiter; ++iter ) {
      int done = 1;
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
	    NAMD_die("Constraint failure in RATTLE algorithm!\n");
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
    if ( iter == maxiter ) {
      NAMD_die("Exceeded maximum number of iterations in rattle1().\n");
    }
    // store data back to patch
    for ( i = 0; i < hgs; ++i ) {
      p[ig+i] = pos[i];
      v[ig+i] = vel[i];
    }
  }

}

//  RATTLE algorithm from Allen & Tildesley
void HomePatch::rattle2(const BigReal timestep, Vector *virial)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  const BigReal dt = timestep / TIMEFACTOR;
  Vector wc(0.,0.,0.);  // constraint virial
  BigReal tol = simParams->rigidTol;
  int maxiter = simParams->rigidIter;
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
        NAMD_die("Hydrogen group error caught in rattle2().  It's a bug!\n");
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
    for ( iter = 0; iter < maxiter; ++iter ) {
      int done = 1;
      for ( i = 0; i < icnt; ++i ) {
	int a = ial[i];  int b = ibl[i];
	Vector vab = vel[a] - vel[b];
	Vector &rab = refab[i];
	BigReal rabsqi = dsqi[i];
	BigReal rvab = rab.x*vab.x + rab.y*vab.y + rab.z*vab.z;
	if ( (fabs(rvab) * dt * rabsqi) > tol ) {
	  Vector dp = rab * (-rvab * redmass[i] * rabsqi);
	  wc.x += dp.x * rab.x;
	  wc.y += dp.y * rab.y;
	  wc.z += dp.z * rab.z;
	  vel[a] += rmass[a] * dp;
	  vel[b] -= rmass[b] * dp;
	  done = 0;
	}
      }
      if ( done ) break;
    }
    if ( iter == maxiter ) {
      NAMD_die("Exceeded maximum number of iterations in rattle2().\n");
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
void HomePatch::mollyMollify(Vector *virial)
{
  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  Vector wc(0.,0.,0.);  // constraint virial
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
	  f[Results::slow][ig+i] = force[i];
	}
  }
  // check that there isn't a constant needed here!
  *virial += wc;
  p_avg.resize(0);
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
  MigrationList *mCur;

  // Drain the migration message buffer
  //for (i=0; i<numMlBuf; i++) {
  //   DebugM(3, "Draining migration buffer ("<<i<<","<<numMlBuf<<")\n");
  //   depositMigration(srcID[i], mlBuf[i]);
  //}
  //numMlBuf = 0;
     
  // realInfo points to migration lists for neighbors we actually have. 
  //    element of mInfo[3][3][3] points to an element of realInfo
  for (i=0; i<numNeighbors; i++) {
    realInfo[i].mList = NULL;
  }

  // Purge the AtomMap
  AtomMap::Object()->unregisterIDs(patchID,atomIDList);

  // Determine atoms that need to migrate
  AtomIDList::iterator atomIDList_i = atomIDList.begin();
  AtomIDList::iterator atomIDList_e = atomIDList.end();
  AtomPropertiesList::iterator a_i = a.begin();
  TransformList::iterator t_i = t.begin();
  PositionList::iterator p_i = p.begin();
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
       if (NULL == (mCur = mInfo[xdev][ydev][zdev]->mList)) {
	 // new: all mList pointers are actually in realInfo[].mList
	 //one of the following below sendMigrationMsg() does delete
	 // [1] in pack of MigrateAtomsMsg (PatchMgr.C)
	 // [2] in recvMigrateAtoms (PatchMgr.C)
	 mCur = mInfo[xdev][ydev][zdev]->mList = new MigrationList;
       }
       Force force[Results::maxNumForces];
       for ( j = 0; j < Results::maxNumForces; ++j ) force[j] = *(f_i[j]);
       DebugM(3,"Migrating atom " << atomIDList_i << " from patch "
		<< patchID << " with position " << p_i << "\n");
       mCur->add(MigrationElem(*atomIDList_i, *a_i, *t_i, *p_i, *v_i, force));

       ++delnum;

     } else {

       if ( delnum ) {
         *(atomIDList_i-delnum) = *atomIDList_i;
         *(a_i-delnum) = *a_i;
         *(t_i-delnum) = *t_i;
         *(p_i-delnum) = *p_i;
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
     ++v_i;
     for ( j = 0; j < Results::maxNumForces; ++j ) ++(f_i[j]);

  }

  int delpos = numAtoms - delnum;
  DebugM(4,"numAtoms " << numAtoms << " deleted " << delnum << "\n");
  atomIDList.del(delpos,delnum);
  a.del(delpos,delnum);
  t.del(delpos,delnum);
  p.del(delpos,delnum);
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
  // PatchID srcPatchID;
  MigrationList *migrationList;

  if (!inMigration) { // We have to buffer changes due to migration
		      // until our patch is in migration mode
    msgbuf[numMlBuf++] = msg;
    return;
  } 

  // srcPatchID = msg->srcPatchID;
  migrationList = msg->migrationList;

  if (migrationList) {
    MigrationListIter mi(*migrationList);
    for (mi = mi.begin(); mi != mi.end(); mi++) {
      DebugM(1,"Migrating atom " << mi->atomID << " to patch "
		<< patchID << " with position " << mi->pos << "\n"); 
      atomIDList.add(mi->atomID);
      a.add(mi->atomProp);
      p.add(lattice.nearest(mi->pos,center,&(mi->trans)));
      t.add(mi->trans);
      v.add(mi->vel);
      for ( int j = 0; j < Results::maxNumForces; ++j )
        f[j].add(mi-> force[j]);
    }
    delete migrationList;
    migrationList = NULL;
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


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: HomePatch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1059 $	$Date: 1999/08/31 15:43:31 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatch.C,v $
 * Revision 1.1059  1999/08/31 15:43:31  jim
 * Cleaned up MOLLY code.
 *
 * Revision 1.1058  1999/08/27 16:40:19  jim
 * Eliminated memory allocation and 1-based indexing from MOLLY code.
 *
 * Revision 1.1057  1999/08/25 18:38:21  jim
 * Fixed bug in divide-by-zero bug fix.
 *
 * Revision 1.1056  1999/08/25 16:50:58  brunner
 * Fixed a divide-by-zero error with fixed atoms in rattle
 *
 * Revision 1.1055  1999/08/20 19:11:11  jim
 * Added MOLLY - mollified impluse method.
 *
 * Revision 1.1054  1999/07/22 21:43:06  jim
 * Fixed virial calculation during RATTLE.
 *
 * Revision 1.1053  1999/07/08 21:26:48  jim
 * Eliminated compiler warnings.
 *
 * Revision 1.1052  1999/05/11 23:56:31  brunner
 * Changes for new charm version
 *
 * Revision 1.1051  1999/04/28 22:50:38  jim
 * Fixed tolerance and increased speed for SHAKE/RATTLE.
 *
 * Revision 1.1050  1999/03/19 23:03:01  jim
 * Fixed bugs in constant pressure code.
 *
 * Revision 1.1049  1999/03/17 17:59:24  jim
 * Eliminated compiler warnings and errors.
 *
 * Revision 1.1048  1998/10/24 19:57:33  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1047  1998/08/11 16:30:27  jim
 * Modified output from periodic boundary simulations to return atoms to
 * internally consistent coordinates.  We store the transformations which
 * were performed and undo them at the end.  It might be better to do this
 * by always keeping the original coordinates and only doing the transform
 * for the nonbonded terms but this works for now.
 *
 * Revision 1.1046  1998/04/14 05:58:24  jim
 * Added automatic correction if hgroupCutoff is too small.  No more warnings.
 * However, performance wil degrade if many groups are below cutoff size.
 *
 * Revision 1.1045  1998/03/26 23:28:29  jim
 * Small changes for KCC port.  Altered use of strstream in ComputeFreeEnergy.
 *
 * Revision 1.1044  1998/03/03 23:05:13  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1043  1998/02/19 01:21:13  jim
 * Small usability changes.
 *
 * Revision 1.1042  1998/02/17 06:39:19  jim
 * SHAKE/RATTLE (rigidBonds) appears to work!!!  Still needs langevin,
 * proper startup, and degree of freedom tracking.
 *
 * Revision 1.1041  1998/01/14 00:40:57  jim
 * Added hydrogen group size checking (vs. hgroupcutoff parameter).
 *
 * Revision 1.1040  1998/01/13 23:10:57  jim
 * Added margin checking - prelude to automatic migration.
 *
 * Revision 1.1039  1997/12/22 21:29:24  jim
 * Proxies no longer send empty arrays back to HomePatch.  Requires some new
 * flags to be set correctly in Sequencer in order to work.  These are:
 *   maxForceMerged - this and faster are added into Results::normal array
 *   maxForceUsed - all forces slower than this are discarded (assumed zero)
 * Generally maxForceMerged doesn't change but maxForceUsed depends on timestep.
 *
 * Revision 1.1038  1997/11/14 04:56:45  jim
 * Added STL-style iterators, eliminated bad algorithm in doAtomMigration.
 *
 * Revision 1.1037  1997/10/06 00:12:31  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 * Revision 1.1036  1997/09/28 10:19:07  milind
 * Fixed priorities, ReductionMgr etc.
 *
 * Revision 1.1035  1997/09/19 08:55:31  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1034  1997/09/19 05:17:43  jim
 * Cleaned up and tweaked hydrogen-group based temporary pairlist
 * generation for roughly a 6% performance improvement.
 *
 * Revision 1.1033  1997/08/26 16:26:15  jim
 * Revamped prioritites for petter performance and easier changes.
 *
 * Revision 1.1032  1997/08/22 20:12:03  milind
 * Turned on Priorities.
 *
 * Revision 1.1031  1997/07/08 15:48:08  milind
 * Made namd2 to work with Origin2000: Again...
 *
 * Revision 1.1030  1997/04/21 00:58:33  jim
 * Fixed hang on patch migration in systems with only one patch.
 *
 * Revision 1.1029  1997/04/10 22:29:11  jim
 * First steps towards combining atom migration messages.
 *
 * Revision 1.1028  1997/04/10 09:13:57  ari
 * Final debugging for compute migration / proxy creation for load balancing.
 * Lots of debug code added, mostly turned off now.
 * Fixed bug in PositionBox when Patch had no dependencies.
 * Eliminated use of iout<<iINFO and misuse of iout in numerous places.
 *                                            Ari & Jim
 *
 * Revision 1.1027  1997/04/08 07:08:35  ari
 * Modification for dynamic loadbalancing - moving computes
 * Still bug in new computes or usage of proxies/homepatches.
 * Works if ldbStrategy is none as before.
 *
 * Revision 1.1026  1997/04/06 22:45:04  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1025  1997/03/31 16:12:51  nealk
 * Atoms now can migrate by hydrogen groups.
 *
 * Revision 1.1024  1997/03/27 20:25:44  brunner
 * Changes for LdbCoordinator, the load balance control BOC
 *
 * Revision 1.1023  1997/03/27 08:04:16  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1022  1997/03/18 18:09:01  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1021  1997/03/12 22:06:40  jim
 * First step towards multiple force returns and multiple time stepping.
 *
 * Revision 1.1020  1997/03/10 17:40:11  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1019  1997/03/06 22:06:02  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1018  1997/02/28 04:47:08  jim
 * Full electrostatics now works with fulldirect on one node.
 *
 * Revision 1.1017  1997/02/26 21:39:27  jim
 * Fixed migration with periodic boundary conditions to correctly
 * re-center migrated atoms on their new home patch.
 *
 * Revision 1.1016  1997/02/26 16:53:09  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1015  1997/02/17 23:46:59  ari
 * Added files for cleaning up atom migration code
 *
 * Revision 1.1014  1997/02/13 23:17:17  ari
 * Fixed a final bug in AtomMigration - numatoms in ComputePatchPair.C not
 * set correctly in atomUpdate()
 *
 * Revision 1.1013  1997/02/13 17:06:22  jim
 * Turned off debugging.
 *
 * Revision 1.1012  1997/02/13 16:17:13  ari
 * Intermediate debuging commit - working to fix deep bug in migration?
 *
 * Revision 1.1011  1997/02/13 04:43:09  jim
 * Fixed initial hanging (bug in PatchMap, but it still shouldn't have
 * happened) and saved migration messages in the buffer from being
 * deleted, but migration still dies (even on one node).
 *
 * Revision 1.1010  1997/02/11 18:51:46  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1009  1997/02/11 16:31:48  nealk
 * Using a flag to determine suspend/awaken in sequencer.
 * Migration turned off (for now).
 *
 * Revision 1.1008  1997/02/10 19:44:29  nealk
 * More debugging code.  Corrected HomePatch/Sequencer awaken/suspend bug.
 *
 * Revision 1.1007  1997/02/10 08:17:30  jim
 * Commented out sending and allocation of unused data.
 *
 * Revision 1.1006  1997/02/07 17:39:38  ari
 * More debugging for atomMigration.
 * Using -w on CC got us some minor fixes
 * using purify got us a major memory problem due to bad sizing of dummy force
 *
 * Revision 1.1005  1997/02/07 07:51:42  jim
 * pInit aliases p, leading to problems.  I commented out any add or del
 * involving it.  Now seems to integrate correctly.
 *
 * Revision 1.1004  1997/02/07 05:42:30  ari
 * Some bug fixing - atom migration on one node works
 * Atom migration on multiple nodes gets SIGSEGV
 *
 * Revision 1.1003  1997/02/06 23:25:07  jim
 * Fixed bugs.
 *
 * Revision 1.1002  1997/02/06 21:20:50  jim
 * Fixed a couple of atom migration bugs.
 *
 * Revision 1.1001  1997/02/06 18:05:28  nealk
 * Modified (added some, turned off others) debug statements.
 *
 * Revision 1.1000  1997/02/06 15:58:26  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:11  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/05 22:18:13  ari
 * Added migration code - Currently the framework is
 * there with compiling code.  This version does
 * crash shortly after migration is complete.
 * Migration appears to complete, but Patches do
 * not appear to be left in a correct state.
 *
 * Revision 1.778  1997/01/28 00:30:35  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:14  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/21 23:04:44  ari
 * Basic framework for atom migration placed into code.  - Non
 * functional since it is not called.  Works currently without
 * atom migration.
 *
 * Revision 1.777  1997/01/17 19:36:10  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.19  1997/01/15 17:59:11  jim
 * now divide by TIMEFACTOR
 *
 * Revision 1.18  1997/01/15 17:09:43  ari
 * Commented out new Atom migration code
 *
 * Revision 1.17  1997/01/10 22:38:37  jim
 * kinetic energy reporting
 *
 * Revision 1.16  1996/12/17 23:58:02  jim
 * proxy result reporting is working
 *
 * Revision 1.15  1996/12/17 22:13:22  jim
 * implemented ProxyDataMsg use
 *
 * Revision 1.14  1996/12/17 08:55:25  jim
 * added node argument to sendProxyAtoms
 *
 * Revision 1.13  1996/12/16 22:52:43  jim
 * added placement new and explicit destructor calls to ProxyAtomsMsg
 *
 * Revision 1.12  1996/12/14 00:02:42  jim
 * debugging ProxyAtomsMsg path to make compute creation work
 *
 * Revision 1.11  1996/12/11 22:31:41  jim
 * added integration methods for Sequencer
 *
 * Revision 1.10  1996/12/05 23:45:09  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/12/05 01:44:16  ari
 * started toward proxy management
 *
 * Revision 1.8  1996/12/01 02:35:47  jim
 * added patchID reporting
 *
 * Revision 1.7  1996/11/30 00:35:51  jim
 * implemented boxClosed(), useSequencer(), runSequencer()
 *
 * Revision 1.6  1996/11/22 01:44:53  jim
 * added calls to service AtomMap
 *
 * Revision 1.5  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/04 21:14:47  jim
 * Moved functionality to Patch
 *
 * Revision 1.3  1996/09/03 22:54:09  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.7  1996/07/15 21:18:15  gursoy
 * *** empty log message ***
 *
 * Revision 1.6  1996/07/10 16:25:29  gursoy
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/10 16:19:52  gursoy
 * *** empty log message ***
 *
 * Revision 1.4  1996/07/02 15:10:27  gursoy
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/11 22:36:23  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/11 20:07:22  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/05/30 21:31:36  gursoy
 * Initial revision
 *
 ***************************************************************************/
