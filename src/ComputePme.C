/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#include "Namd.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputePme.h"
#include "ComputePmeMsgs.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "SimParameters.h"

#include "PmeCoulomb.h"

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif

class ComputePmeMaster {
private:
  friend class ComputePme;
  ComputePme *host;
  ComputePmeMaster(ComputePme *);
  ~ComputePmeMaster();
  void recvData(ComputePmeDataMsg *);
  ResizeArray<int> homeNode;
  ResizeArray<int> endForNode;
  int numWorkingPes;
  int numLocalAtoms;
  PmeParticle *localData;
  SubmitReduction *reduction;
  int runcount;
  PmeCoulomb *myPme;
};

ComputePme::ComputePme(ComputeID c, ComputeMgr *m) :
  ComputeHomePatches(c), comm(m)
{
  DebugM(4,"ComputePme created.\n");
  useAvgPositions = 1;

  int numWorkingPes = CkNumPes();
  {
    int npatches=(PatchMap::Object())->numPatches();
    if ( numWorkingPes > npatches ) numWorkingPes = npatches;
  }

  masterNode = numWorkingPes - 1;
  if ( CkMyPe() == masterNode ) {
    master = new ComputePmeMaster(this);
    master->numWorkingPes = numWorkingPes;
  }
  else master = 0;
}

ComputePme::~ComputePme()
{
  delete master;
}

void ComputePme::doWork()
{
  DebugM(4,"Entering ComputePme::doWork().\n");

  PmeParticle *localData;

  ResizeArrayIter<PatchElem> ap(patchList);

  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics )
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      Position *x = (*ap).positionBox->open();
      AtomProperties *a = (*ap).atomBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).atomBox->close(&a);
      (*ap).forceBox->close(&r);
    }
    if ( master ) {
      master->reduction->submit();
    }
    return;
  }

  // allocate storage
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  Lattice lattice = patchList[0].p->flags.lattice;

  localData = new PmeParticle[numLocalAtoms];  // given to message

  // get positions and charges
  PmeParticle * data_ptr = localData;
  const BigReal coloumb_sqrt = sqrt( COLOUMB * ComputeNonbondedUtil::scaling
				* ComputeNonbondedUtil::dielectric_1 );

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Position *x = (*ap).positionBox->open();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    AtomProperties *a = (*ap).atomBox->open();
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      Vector tmp = lattice.delta(x[i]);
      data_ptr->x = tmp.x;
      data_ptr->y = tmp.y;
      data_ptr->z = tmp.z;
      data_ptr->cg = coloumb_sqrt * a[i].charge;
      ++data_ptr;
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
    (*ap).atomBox->close(&a);
  }

  // send data to master
  ComputePmeDataMsg *msg = new ComputePmeDataMsg;
  msg->node = CkMyPe();
  msg->numParticles = numLocalAtoms;
  msg->particles = localData;
  comm->sendComputePmeData(msg);
}

void ComputePme::recvData(ComputePmeDataMsg *msg)
{
  if ( master ) {
    master->recvData(msg);
  }
  else NAMD_die("ComputePme::master is NULL!");
}

ComputePmeMaster::ComputePmeMaster(ComputePme *h) :
  host(h), numLocalAtoms(0), runcount(0)
{
  DebugM(4,"ComputePmeMaster created.\n");

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  Molecule * molecule = Node::Object()->molecule;
  localData = new PmeParticle[molecule->numAtoms];
  SimParameters * simParams = Node::Object()->simParameters;
  PmeGrid grid;
  grid.K1 = simParams->PMEGridSizeX;
  grid.K2 = simParams->PMEGridSizeY;
  grid.K3 = simParams->PMEGridSizeZ;
  grid.order = simParams->PMEInterpOrder;
  myPme = new PmeCoulomb(grid, molecule->numAtoms);
}

ComputePmeMaster::~ComputePmeMaster()
{
  delete reduction;
  delete [] localData;
  delete myPme;
}

void ComputePmeMaster::recvData(ComputePmeDataMsg *msg)
{ 
  DebugM(4,"ComputePmeMaster::recvData() " << msg->numParticles
	<< " particles from node " << msg->node << "\n");

  {
    homeNode.add(msg->node);
    PmeParticle *data_ptr = localData + numLocalAtoms;
    for ( int j = 0; j < msg->numParticles; ++j, ++data_ptr ) {
      *data_ptr = msg->particles[j];
    }
    numLocalAtoms += msg->numParticles;
    endForNode.add(numLocalAtoms);
    delete msg;
  }

  if ( homeNode.size() < numWorkingPes ) return;  // messages outstanding

  DebugM(4,"ComputePmeMaster::recvData() running serial code.\n");

  // single processor version

  Lattice lattice = host->getFlags()->lattice;
  SimParameters * simParams = Node::Object()->simParameters;
  int i;
  Vector *localResults;
  double recip_vir[6];
  BigReal ewaldcof = ComputeNonbondedUtil::ewaldcof;
  localResults = new Vector[numLocalAtoms];

  // perform calculations
  BigReal electEnergy = 0;

  // calculate self energy
  PmeParticle *data_ptr = localData;
  for(i=0; i<numLocalAtoms; ++i)
  {
    electEnergy += data_ptr->cg * data_ptr->cg;
    ++data_ptr;
  }
  electEnergy *= -1. * ewaldcof / SQRT_PI;

  DebugM(4,"Ewald self energy: " << electEnergy << "\n");

  double recipEnergy;
  DebugM(4,"Calling compute_recip.\n");
  double pme_start_time = 0;
  if ( runcount == 1 ) pme_start_time = CmiTimer();
  // Last argument should be nonzero if this is the last call
  // Compute it using mytime, tsteps, etc.
  // I'll just let the runtime environment take care of it for now.
  recipEnergy = myPme->compute_recip(localData, lattice, ewaldcof, recip_vir, localResults);
  if ( runcount == 1 ) {
    iout << iINFO << "PME reciprocal sum CPU time per evaluation: "
         << (CmiTimer() - pme_start_time) << "\n" << endi;
  }
  electEnergy += recipEnergy;
  DebugM(4,"Returned from PmeRecipCoulomb->calc_recip.\n");
  // No need to reverse sign from new code

  // send out reductions
  DebugM(4,"Timestep : " << host->getFlags()->step << "\n");
  DebugM(4,"Reciprocal sum energy: " << electEnergy << "\n");
  DebugM(4,"Reciprocal sum virial: " << recip_vir[0] << " " <<
	recip_vir[1] << " " << recip_vir[2] << " " << recip_vir[3] << " " <<
	recip_vir[4] << " " << recip_vir[5] << "\n");
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += electEnergy;
  reduction->item(REDUCTION_VIRIAL_SLOW_X) += (BigReal)(recip_vir[0]);
  reduction->item(REDUCTION_VIRIAL_SLOW_Y) += (BigReal)(recip_vir[3]);
  reduction->item(REDUCTION_VIRIAL_SLOW_Z) += (BigReal)(recip_vir[5]);
  reduction->submit();

  Vector *results_ptr;
  results_ptr = localResults;

  numLocalAtoms = 0;
  for ( i = 0; i < homeNode.size(); ++i ) {
    ComputePmeResultsMsg *msg = new ComputePmeResultsMsg;
    msg->node = homeNode[i];
    msg->numParticles = endForNode[i] - numLocalAtoms;
    msg->forces = new Vector[msg->numParticles];
    for ( int j = 0; j < msg->numParticles; ++j, ++results_ptr ) {
      msg->forces[j] = *results_ptr;
    }
    numLocalAtoms = endForNode[i];
    host->comm->sendComputePmeResults(msg,homeNode[i]);
  }
  delete [] localResults;

  // reset
  runcount += 1;
  numLocalAtoms = 0;
  homeNode.resize(0);
  endForNode.resize(0);

}

void ComputePme::recvResults(ComputePmeResultsMsg *msg)
{
  if ( CkMyPe() != msg->node ) {
    NAMD_die("ComputePme results sent to wrong node!\n");
    return;
  }
  if ( numLocalAtoms != msg->numParticles ) {
    NAMD_die("ComputePme sent wrong number of results!\n");
    return;
  }

  Vector *results_ptr = msg->forces;
  ResizeArrayIter<PatchElem> ap(patchList);

  // add in forces
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      f[i].x += results_ptr->x;
      f[i].y += results_ptr->y;
      f[i].z += results_ptr->z;
      ++results_ptr;
    }
  
    (*ap).forceBox->close(&r);
  }
   
  delete msg;

}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputePme.C,v $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.6 $	$Date: 1999/09/03 20:46:10 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputePme.C,v $
 * Revision 1.6  1999/09/03 20:46:10  jim
 * Support for non-orthogonal periodic boundary conditions.
 *
 * Revision 1.5  1999/08/20 19:11:11  jim
 * Added MOLLY - mollified impluse method.
 *
 * Revision 1.4  1999/06/17 17:05:43  jim
 * Renamed seq to step in most places.  Now has meaning only to user.
 *
 * Revision 1.3  1999/06/17 15:46:12  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.2  1999/06/09 15:02:07  jim
 * Added nonbondedScaling parameter.
 *
 * Revision 1.1  1999/06/08 14:52:06  jim
 * Incorporated Justin's faster PME code along side DPME.
 *
 * Revision 1.11  1999/05/11 23:56:17  brunner
 * Changes for new charm version
 *
 * Revision 1.10  1999/03/22 20:55:43  jim
 * DPME now compiles using C compiler.
 *
 * Revision 1.9  1999/03/10 00:52:26  jim
 * Adding timing output for PME reciprocal sum.
 *
 * Revision 1.8  1999/02/17 04:09:54  jim
 * Fixes to make optional force modules work with more nodes than patches.
 *
 * Revision 1.7  1999/02/12 21:55:15  jim
 * Fixed sign of virial from DPME.
 *
 * Revision 1.6  1999/01/06 00:56:20  jim
 * All compute objects except DPMTA now return diagonal of virial tensor.
 *
 * Revision 1.5  1998/09/14 21:45:25  jim
 * Turned off DPME and load balancer verbose output.
 *
 * Revision 1.4  1998/06/18 14:48:00  jim
 * Split virial into NORMAL, NBOND, and SLOW parts to match force classes.
 *
 * Revision 1.3  1998/04/15 22:13:49  jim
 * Make depends returns same results regardless of DPME, DPMTA, TCL or MDCOMM.
 *
 * Revision 1.2  1998/04/10 04:15:57  jim
 * Finished incorporating DPME.
 *
 * Revision 1.1  1998/04/07 00:52:37  jim
 * Added DPME interface.
 *
 *
 ***************************************************************************/
