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
#include "ComputeFullDirect.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "Communicate.h"
#include "Lattice.h"
//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

ComputeFullDirect::ComputeFullDirect(ComputeID c) : ComputeHomePatches(c)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  useAvgPositions = 1;
}

ComputeFullDirect::~ComputeFullDirect()
{
  delete reduction;
}

BigReal calc_fulldirect(BigReal *data1, BigReal *results1, int n1,
                        BigReal *data2, BigReal *results2, int n2,
			int selfmode, Lattice *lattice, Vector &virial)
{
  if ( lattice->a_p() || lattice->b_p() || lattice->c_p() ) {
    #define FULLDIRECT_PERIODIC
    #include "ComputeFullDirectBase.h"
  } else {
    #undef FULLDIRECT_PERIODIC
    #include "ComputeFullDirectBase.h"
  }
}

void ComputeFullDirect::doWork()
{
  int numLocalAtoms;
  BigReal *localData;
  BigReal *localResults;
  BigReal *newLocalResults;
  register BigReal *local_ptr;
  Lattice *lattice;

  int numWorkingPes = CkNumPes();
  {
    int npatches=(PatchMap::Object())->numPatches();
    if ( numWorkingPes > npatches ) numWorkingPes = npatches;
  }
  

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
    reduction->submit();
    return;
  }

  // allocate storage
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  localData = new BigReal[4*numLocalAtoms];	// freed at end of this method
  localResults = new BigReal[3*numLocalAtoms];	// freed at end of this method
  newLocalResults = new BigReal[3*numLocalAtoms];  // freed at end of this method

  lattice = &((*(ap.begin())).p->lattice);

  // get positions and charges
  local_ptr = localData;
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
      *(local_ptr++) = x[i].x;
      *(local_ptr++) = x[i].y;
      *(local_ptr++) = x[i].z;
      *(local_ptr++) = a[i].charge;
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
    (*ap).atomBox->close(&a);
  } 

  // zero out forces
  local_ptr = localResults;
  for(int j=0; j<numLocalAtoms; ++j)
  {
    *(local_ptr++) = 0.;
    *(local_ptr++) = 0.;
    *(local_ptr++) = 0.;
  }

  // perform calculations
  BigReal electEnergy = 0;
  Vector virial(0.,0.,0.);

#define PEMOD(N) (((N)+numWorkingPes)%numWorkingPes)

  int numStages = numWorkingPes / 2 + 2;
  int lastStage = numStages - 2;
  int sendDataPE = PEMOD(CkMyPe()+1);
  int recvDataPE = PEMOD(CkMyPe()-1);
  int sendResultsPE = PEMOD(CkMyPe()-1);
  int recvResultsPE = PEMOD(CkMyPe()+1);
  int numRemoteAtoms = numLocalAtoms;
  int oldNumRemoteAtoms = 0;
  BigReal *remoteData = 0;
  BigReal *remoteResults = 0;
  register BigReal *remote_ptr;
  register BigReal *end_ptr;

  MOStream *sendDataMsg=CpvAccess(comm)->
		newOutputStream(sendDataPE, FULLTAG, BUFSIZE);
  MIStream *recvDataMsg=CpvAccess(comm)->
		newInputStream(recvDataPE, FULLTAG);

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
      sendDataMsg->put(numRemoteAtoms);
      sendDataMsg->put(4*numRemoteAtoms,(stage?remoteData:localData));
      sendDataMsg->end();
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
      electEnergy += calc_fulldirect(
        localData,localResults,numLocalAtoms,
        localData,localResults,numLocalAtoms,1,lattice,virial);
    }
    else if ( stage < lastStage ||
            ( stage == lastStage && ( numWorkingPes % 2 ) ) )
    {  // full other interaction
      DebugM(4,"full other interaction\n");
      electEnergy += calc_fulldirect(
        localData,localResults,numLocalAtoms,
        remoteData,remoteResults,numRemoteAtoms,0,lattice,virial);
    }
    else if ( stage == lastStage )
    {  // half other interaction
      DebugM(4,"half other interaction\n");
      if ( CkMyPe() < ( numWorkingPes / 2 ) )
        electEnergy += calc_fulldirect(
          localData,localResults,numLocalAtoms/2,
          remoteData,remoteResults,numRemoteAtoms,0,lattice,virial);
      else
        electEnergy += calc_fulldirect(
          localData,localResults,numLocalAtoms,
          remoteData + 4*(numRemoteAtoms/2),
          remoteResults + 3*(numRemoteAtoms/2),
          numRemoteAtoms - (numRemoteAtoms/2), 0,lattice,virial);
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
      recvDataMsg->get(numRemoteAtoms);
      remoteData = new BigReal[4*numRemoteAtoms];
      recvDataMsg->get(4*numRemoteAtoms,remoteData);
    }

  }

  delete sendDataMsg;
  delete recvDataMsg;

  // send out reductions
  DebugM(4,"Full-electrostatics energy: " << electEnergy << "\n");
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += electEnergy;
  reduction->item(REDUCTION_VIRIAL_SLOW_X) += virial.x;
  reduction->item(REDUCTION_VIRIAL_SLOW_Y) += virial.y;
  reduction->item(REDUCTION_VIRIAL_SLOW_Z) += virial.z;
  reduction->submit();

  // add in forces
  local_ptr = localResults;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      f[i].x += *(local_ptr++);
      f[i].y += *(local_ptr++);
      f[i].z += *(local_ptr++);
    }

    (*ap).forceBox->close(&r);
  }

  // free storage
  delete [] localData;		// allocated at beginning of this method
  delete [] localResults;	// allocated at beginning of this method
  delete [] newLocalResults;	// allocated at beginning of this method
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1023 $	$Date: 1999/09/03 20:46:08 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeFullDirect.C,v $
 * Revision 1.1023  1999/09/03 20:46:08  jim
 * Support for non-orthogonal periodic boundary conditions.
 *
 * Revision 1.1022  1999/08/20 19:11:08  jim
 * Added MOLLY - mollified impluse method.
 *
 * Revision 1.1021  1999/06/17 15:46:05  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.1020  1999/05/11 23:56:20  brunner
 * Changes for new charm version
 *
 * Revision 1.1019  1999/04/13 08:17:58  jim
 * Eliminated memory leaks and corruption in MStream and ComputeFullDirect.
 *
 * Revision 1.1018  1999/02/17 04:09:55  jim
 * Fixes to make optional force modules work with more nodes than patches.
 *
 * Revision 1.1017  1999/01/06 00:56:21  jim
 * All compute objects except DPMTA now return diagonal of virial tensor.
 *
 * Revision 1.1016  1998/10/24 19:57:23  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1015  1998/06/18 14:48:02  jim
 * Split virial into NORMAL, NBOND, and SLOW parts to match force classes.
 *
 * Revision 1.1014  1998/03/30 21:01:16  jim
 * Added nearest-image support for periodic boundary conditions to full direct.
 *
 * Revision 1.1013  1997/12/17 10:28:07  jim
 * Full direct electrostatics now works in parallel.
 *
 * Revision 1.1012  1997/04/06 22:44:58  ari
 * Add priorities to messages.  Mods to help proxies without computes.
 * Added quick enhancement to end of list insertion of ResizeArray(s)
 *
 * Revision 1.1011  1997/03/20 23:53:36  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1010  1997/03/19 11:54:07  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
