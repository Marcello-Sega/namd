/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifdef DPME

#include "Namd.h"
#include "Node.h"
#include "PatchMap.h"
#include "AtomMap.h"
#include "ComputeDPME.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "Communicate.h"
#include "dpme2.h"
//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

ComputeDPME::ComputeDPME(ComputeID c) : ComputeHomePatches(c)
{
  DebugM(4,"ComputeDPME created.\n");
  reduction->Register(REDUCTION_ELECT_ENERGY);
  reduction->Register(REDUCTION_VIRIAL);
}

ComputeDPME::~ComputeDPME()
{
  reduction->unRegister(REDUCTION_ELECT_ENERGY);
  reduction->unRegister(REDUCTION_VIRIAL);
}


BigReal calc_fakeDPME(BigReal *data1, BigReal *results1, int n1,
                        BigReal *data2, BigReal *results2, int n2, int selfmode)
{
  const BigReal coloumb = COLOUMB * ComputeNonbondedUtil::dielectric_1;
  BigReal *dp1 = data1;
  BigReal *rp1 = results1;
  int j_begin = 0;
  register BigReal electEnergy = 0.;
  for(int i=0; i<n1; ++i)
  {
    register BigReal p_i_x = *(dp1++);
    register BigReal p_i_y = *(dp1++);
    register BigReal p_i_z = *(dp1++);
    register BigReal kq_i = coloumb * *(dp1++);
    register BigReal f_i_x = 0.;
    register BigReal f_i_y = 0.;
    register BigReal f_i_z = 0.;
    if ( selfmode )
    {
      ++j_begin; data2 += 4; results2 += 3;
    }
    register BigReal *dp2 = data2;
    register BigReal *rp2 = results2;
    register int n2c = n2;
    register int j;
    for( j = j_begin; j<n2c; ++j)
    {
      register BigReal p_ij_x = p_i_x - *(dp2++);
      register BigReal p_ij_y = p_i_y - *(dp2++);
      register BigReal p_ij_z = p_i_z - *(dp2++);

      register BigReal r_1;
      r_1 = 1./sqrt(p_ij_x * p_ij_x + p_ij_y * p_ij_y + p_ij_z * p_ij_z);
      register BigReal f = *(dp2++) * kq_i * r_1;
      electEnergy += f;
      f *= r_1*r_1;
      p_ij_x *= f;
      p_ij_y *= f;
      p_ij_z *= f;
      f_i_x += p_ij_x;
      f_i_y += p_ij_y;
      f_i_z += p_ij_z;
      *(rp2++) -= p_ij_x;
      *(rp2++) -= p_ij_y;
      *(rp2++) -= p_ij_z;
    }
    *(rp1++) += f_i_x;
    *(rp1++) += f_i_y;
    *(rp1++) += f_i_z;
  }

  return electEnergy;
}

// These are needed to fix up argument mismatches in DPME.

extern int cfftf(int *, double *, double *);

int cfftf(int *n, doublecomplex *c, double *wsave) {
  // Casting (doublecomplex*) to (double*) is probably OK.
  return cfftf(n, (double *)c, wsave);
}

extern int cfftb(int *, double *, double *);

int cfftb(int *n, doublecomplex *c, double *wsave) {
  // Casting (doublecomplex*) to (double*) is probably OK.
  return cfftb(n, (double *)c, wsave);
}

extern int cffti1(int *, double *, int *);

int cffti1(int *n, double *wa, double *ifac) {
  // Casting (double*) to (int*) is dangerous if sizes differ!!!
  return cffti1(n, wa, (int *)ifac);
}

extern int cfftf1(int *n, double *c, double *ch, double *wa, int *ifac);

int cfftf1(int *n, double *c, double *ch, double *wa, double *ifac) {
  // Casting (double*) to (int*) is dangerous if sizes differ!!!
  return cfftf1(n, c, ch, wa, (int *)ifac);
}

extern int cfftb1(int *n, double *c, double *ch, double *wa, int *ifac);

int cfftb1(int *n, double *c, double *ch, double *wa, double *ifac) {
  // Casting (double*) to (int*) is dangerous if sizes differ!!!
  return cfftb1(n, c, ch, wa, (int *)ifac);
}

void ComputeDPME::doWork()
{
  DebugM(4,"Entering ComputeDPME::doWork().\n");

  int i;
  int numLocalAtoms;
  Pme2Particle *localData;

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
    reduction->submit(patchList[0].p->flags.seq, REDUCTION_ELECT_ENERGY, 0.);
    reduction->submit(patchList[0].p->flags.seq, REDUCTION_VIRIAL, 0.0);
    return;
  }

  // allocate storage
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  Lattice lattice = patchList[0].p->flags.lattice;

  localData = new Pme2Particle[numLocalAtoms];  // freed at end of this method

  // get positions and charges
  Pme2Particle * data_ptr = localData;
  const BigReal coloumb_sqrt = sqrt( COLOUMB * ComputeNonbondedUtil::dielectric_1 );
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Position *x = (*ap).positionBox->open();
    AtomProperties *a = (*ap).atomBox->open();
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      Vector tmp = lattice.delta(x[i]);
      data_ptr->x = tmp.x;
      data_ptr->y = tmp.y;
      data_ptr->z = tmp.z;
      data_ptr->cg = coloumb_sqrt * a[i].charge;
      data_ptr->id = a[i].id;
      ++data_ptr;
    }

    (*ap).positionBox->close(&x);
    (*ap).atomBox->close(&a);
  }

/*
  // zero out forces
  PmeVector *results_ptr = localResults;
  for(int j=0; j<numLocalAtoms; ++j)
  {
    results_ptr->x = 0.;
    results_ptr->y = 0.;
    results_ptr->z = 0.;
    ++results_ptr;
  }

#define PEMOD(N) (((N)+CNumPes())%CNumPes())

  int numStages = CNumPes() / 2 + 2;
  int lastStage = numStages - 2;
  int sendDataPE = PEMOD(CMyPe()+1);
  int recvDataPE = PEMOD(CMyPe()-1);
  int sendResultsPE = PEMOD(CMyPe()-1);
  int recvResultsPE = PEMOD(CMyPe()+1);
  int numRemoteAtoms = numLocalAtoms;
  int oldNumRemoteAtoms;
  BigReal *remoteData = 0;
  BigReal *remoteResults = 0;
  register BigReal *remote_ptr;
  register BigReal *end_ptr;

  for ( int stage = 0; stage < numStages; ++stage )
  {
    // send remoteResults to sendResultsPE
    if ( stage > 1 )
    {
      DebugM(4,"send remoteResults to sendResultsPE " << sendResultsPE << "\n");
      MOStream *msg=CpvAccess(comm)->
		newOutputStream(sendResultsPE, FULLFORCETAG, BUFSIZE);
      msg->put(3*oldNumRemoteAtoms,remoteResults);
      delete [] remoteResults;
      msg->end();
      delete msg;
      sendResultsPE = PEMOD(sendResultsPE-1);
    }

    // send remoteData to sendDataPE
    if ( stage < lastStage )
    {
      DebugM(4,"send remoteData to sendDataPE " << sendDataPE << "\n");
      MOStream *msg=CpvAccess(comm)->
		newOutputStream(sendDataPE, FULLTAG, BUFSIZE);
      msg->put(numRemoteAtoms);
      msg->put(4*numRemoteAtoms,(stage?remoteData:localData));
      msg->end();
      delete msg;
    }

    // allocate new result storage
    if ( stage > 0 && stage <= lastStage )
    {
      DebugM(4,"allocate new result storage\n");
      remoteResults = new BigReal[3*numRemoteAtoms];
      remote_ptr = remoteResults;
      end_ptr = remoteResults + 3*numRemoteAtoms;
      for ( ; remote_ptr != end_ptr; ++remote_ptr ) *remote_ptr = 0.;
    }

    // do calculation
    if ( stage == 0 )
    {  // self interaction
      DebugM(4,"self interaction\n");
      electEnergy += calc_fakeDPME(
        localData,localResults,numLocalAtoms,
        localData,localResults,numLocalAtoms,1);
    }
    else if ( stage < lastStage ||
            ( stage == lastStage && ( CNumPes() % 2 ) ) )
    {  // full other interaction
      DebugM(4,"full other interaction\n");
      electEnergy += calc_fakeDPME(
        localData,localResults,numLocalAtoms,
        remoteData,remoteResults,numRemoteAtoms,0);
    }
    else if ( stage == lastStage )
    {  // half other interaction
      DebugM(4,"half other interaction\n");
      if ( CMyPe() < ( CNumPes() / 2 ) )
        electEnergy += calc_fakeDPME(
          localData,localResults,numLocalAtoms/2,
          remoteData,remoteResults,numRemoteAtoms,0);
      else
        electEnergy += calc_fakeDPME(
          localData,localResults,numLocalAtoms,
          remoteData + 4*(numRemoteAtoms/2),
          remoteResults + 3*(numRemoteAtoms/2),
          numRemoteAtoms - (numRemoteAtoms/2), 0);
    }

    delete [] remoteData;  remoteData = 0;
    oldNumRemoteAtoms = numRemoteAtoms;

    // receive newLocalResults from recvResultsPE
    if ( stage > 1 )
    {
      DebugM(4,"receive newLocalResults from recvResultsPE "
						<< recvResultsPE << "\n");
      MIStream *msg=CpvAccess(comm)->
		newInputStream(recvResultsPE, FULLFORCETAG);
      msg->get(3*numLocalAtoms,newLocalResults);
      delete msg;
      recvResultsPE = PEMOD(recvResultsPE+1);
      remote_ptr = newLocalResults;
      local_ptr = localResults;
      end_ptr = localResults + 3*numLocalAtoms;
      for ( ; local_ptr != end_ptr; ++local_ptr, ++remote_ptr )
	*local_ptr += *remote_ptr;
    }

    // receive remoteData from recvDataPE
    if ( stage < lastStage )
    {
      DebugM(4,"receive remoteData from recvDataPE "
						<< recvDataPE << "\n");
      MIStream *msg=CpvAccess(comm)->
		newInputStream(recvDataPE, FULLTAG);
      msg->get(numRemoteAtoms);
      remoteData = new BigReal[4*numRemoteAtoms];
      msg->get(4*numRemoteAtoms,remoteData);
      delete msg;
    }

  }

*/

  // single processor version

  SimParameters * simParams = Node::Object()->simParameters;

  AtomInfo atom_info;
  atom_info.numatoms = numLocalAtoms;
  atom_info.atompnt = 0;  // not used
  atom_info.freepnt = 0;  // not used
  atom_info.nlocal = numLocalAtoms;
  atom_info.nother = 0;

  BoxInfo box_info;
  box_info.box.x = lattice.a();
  box_info.box.y = lattice.b();
  box_info.box.z = lattice.c();
  box_info.box.origin = 0.;  // why only one number?
  box_info.prd.x = box_info.box.x;
  box_info.prd.y = box_info.box.y;
  box_info.prd.z = box_info.box.z;
  box_info.prd2.x = 0.5 * box_info.box.x;
  box_info.prd2.y = 0.5 * box_info.box.y;
  box_info.prd2.z = 0.5 * box_info.box.z;
  box_info.cutoff = simParams->cutoff;
  box_info.rs = simParams->cutoff;
  box_info.mc2.x = 2. * ( box_info.prd.x - box_info.cutoff );
  box_info.mc2.y = 2. * ( box_info.prd.y - box_info.cutoff );
  box_info.mc2.z = 2. * ( box_info.prd.z - box_info.cutoff );
  box_info.ewaldcof = ComputeNonbondedUtil::ewaldcof;
  box_info.dtol = simParams->PMETolerance;
  for (i = 0; i <= 8; i++) {
    box_info.recip[i] = 0.; /* assume orthogonal box */
  }
  box_info.recip[0] = 1.0/box_info.box.x;
  box_info.recip[4] = 1.0/box_info.box.y;
  box_info.recip[8] = 1.0/box_info.box.z;

  GridInfo grid_info;
  grid_info.order = simParams->PMEInterpOrder;
  grid_info.nfftgrd.x = simParams->PMEGridSizeX;
  grid_info.nfftgrd.y = simParams->PMEGridSizeY;
  grid_info.nfftgrd.z = simParams->PMEGridSizeZ;
  grid_info.nfft = 0;
  grid_info.volume = lattice.volume();

  PeInfo pe_info;  // hopefully this isn't used for anything
  pe_info.myproc.node = 0;
  pe_info.myproc.nprocs = 1;
  pe_info.myproc.ncube = 0;
  pe_info.inst_node[0] = 0;
  pe_info.igrid = 0;

  PmeVector *localResults;
  double recip_vir[6];
  int time_count = 0;
  int tsteps = 1;
  double mytime = 0.;

  // perform calculations
  BigReal electEnergy = 0;

  // calculate self energy
  data_ptr = localData;
  for(i=0; i<numLocalAtoms; ++i)
  {
    electEnergy += data_ptr->cg * data_ptr->cg;
    ++data_ptr;
  }
  electEnergy *= -1. * box_info.ewaldcof / SQRT_PI;

  DebugM(4,"Ewald self energy: " << electEnergy << "\n");

  DebugM(4,"Calling dpme_eval_recip().\n");

  electEnergy += dpme_eval_recip( atom_info, localData - 1, &localResults,
			recip_vir, grid_info, box_info, pe_info,
			time_count, tsteps, &mytime );

  DebugM(4,"Returned from dpme_eval_recip().\n");

  // send out reductions
  DebugM(4,"Timestep : " << patchList[0].p->flags.seq << "\n");
  DebugM(4,"Reciprocal sum energy: " << electEnergy << "\n");
  DebugM(4,"Reciprocal sum virial: " << recip_vir[0] << " " <<
	recip_vir[1] << " " << recip_vir[2] << " " << recip_vir[3] << " " <<
	recip_vir[4] << " " << recip_vir[5] << "\n");
  reduction->submit(patchList[0].p->flags.seq, REDUCTION_ELECT_ENERGY, electEnergy);
  reduction->submit(patchList[0].p->flags.seq, REDUCTION_VIRIAL, 
			(BigReal)(recip_vir[0] + recip_vir[3] + recip_vir[5]) );

  // add in forces
  PmeVector *results_ptr = localResults + 1;
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

  // free storage
  delete [] localData;		// allocated at beginning of this method

  DebugM(4,"Leaving ComputeDPME::doWork().\n");
}

#endif  // DPME
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeDPME.C,v $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1998/04/07 00:52:37 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDPME.C,v $
 * Revision 1.1  1998/04/07 00:52:37  jim
 * Added DPME interface.
 *
 *
 ***************************************************************************/
