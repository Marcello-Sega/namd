/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#include "Namd.h"
#include "Node.h"
#include "ComputeMDComm.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.top.h"
#include "mdcomm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

#ifdef MDCOMM
#include <rapp.h>
#include <rapp_app.h>

//////////////////////////////// receive simple force vectors 

//  Get a new force (int number of atoms, a list of indicies, a list of 
//    X coordinates, Y coords, and Z coords
//
int mdcomm_app_force_cmd(rapp_active_socket_t *sock)
{
  int tag;
  int n;
  int *indicies;
  float *fx, *fy, *fz;


  /* get the number of elements */
  rapp_recv(sock, &tag, &n, RAPP_INT);
  //  namdInfo << "Got " << n << " atoms from VMD" << sendmsg;

  /* make room for them */
  indicies = (int *)   malloc(n * sizeof(int));
  fx = (float *) malloc(3 * n * sizeof(float));
  fy = fx + n;
  fz = fy + n;

  /* get the indicies */
  rapp_recv(sock, &tag, indicies, RAPP_INT);
  /* get the x forces */
  rapp_recv(sock, &tag, fx, RAPP_FLOAT);
  /* get the y forces */
  rapp_recv(sock, &tag, fy, RAPP_FLOAT);
  /* get the z forces */
  rapp_recv(sock, &tag, fz, RAPP_FLOAT);


  // transfer the data to namd
  // note that this function will be called on node 0 only

  mdcomm_transfer_vmdForceData(n,indicies,fx);

  /* and print everything out
  fprintf(stderr, "Got new forces\n");
  for (i=0; i<n; i++) {
    fprintf(stderr, "%d \t: index %d\t %f\t %f\t %f\n",
	    i, indicies[i], fx[i], fy[i], fz[i]);
  }
  fflush(stderr);

  */

  // free indicies and fx
  free(indicies);
  free(fx);

  return 0;
}

void mdcomm_transfer_vmdForceData(int n,int *indicies, float *fx)
{

  ComputeMDComm::vmd_atoms.resize(n);
  ComputeMDComm::vmd_forces.resize(n);

  int i;
  AtomIDList::iterator a_i = ComputeMDComm::vmd_atoms.begin();
  ForceList::iterator f_i = ComputeMDComm::vmd_forces.begin();

  for ( i = 0; i < n; ++i, ++a_i, ++f_i ) {
    *a_i = indicies[i];
    f_i->x = fx[i];
    f_i->y = fx[n+i];
    f_i->z = fx[2*n+i];
  }

}

//////////////////////////////////////////// end of force vectors
#endif  // MDCOMM

AtomIDList ComputeMDComm::vmd_atoms;
ForceList ComputeMDComm::vmd_forces;

ComputeMDComm::ComputeMDComm(ComputeGlobal *h) : ComputeGlobalMaster(h) {
  DebugM(3,"Constructing ComputeMDComm\n");
  configMsg = 0;  resultsMsg = 0;
}

ComputeMDComm::~ComputeMDComm() {
  DebugM(3,"Destructing ComputeMDComm\n");
}


void ComputeMDComm::initialize() {
  DebugM(4,"Initializing master\n");

  configMsg = new (MsgIndex(ComputeGlobalConfigMsg)) ComputeGlobalConfigMsg;

  iout << iDEBUG << "MDComm - initialize()\n" << endi; 

/*  Done in Output.C
  if (rapp_app_set_handler(Node::Object()->output->vmdHandle, MDCOMM_FORCE_CMD,
                            (int (*)())mdcomm_app_force_cmd) == -1)
    rapp_perror("rapp_appd_set_handler for mdcomm_app_force_cmd");
*/

  // Send config to clients
  host->comm->sendComputeGlobalConfig(configMsg);
  configMsg = 0;
}


void ComputeMDComm::calculate() {
  DebugM(4,"Calculating forces on master\n");

  resultsMsg = new (MsgIndex(ComputeGlobalResultsMsg)) ComputeGlobalResultsMsg;
  resultsMsg->gforce.resize(gmass.size());

  iout << iDEBUG << "MDComm - calculate()\n" << endi; 

  // Copy forces from static array
  int numForces = vmd_atoms.size();
  resultsMsg->aid.resize(numForces);
  resultsMsg->f.resize(numForces);
  int i;
  AtomIDList::iterator aid_i = resultsMsg->aid.begin();
  AtomIDList::iterator aid_i2 = vmd_atoms.begin();
  ForceList::iterator f_i = resultsMsg->f.begin();
  ForceList::iterator f_i2 = vmd_forces.begin();
  for ( i = 0; i < numForces; ++i ) {
    aid_i[i] = aid_i2[i];
    f_i[i] = f_i2[i];
  }

  // Send results to clients
  DebugM(3,"Sending results (" << resultsMsg->aid.size() << " forces) on master\n");
  if ( resultsMsg->reconfig ) {
    DebugM(4,"Sending new configuration (" <<
			resultsMsg->newaid.size() << " atoms) on master\n");
  }
  host->comm->sendComputeGlobalResults(resultsMsg);
  resultsMsg = 0;
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1999/02/17 04:09:56 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMDComm.C,v $
 * Revision 1.2  1999/02/17 04:09:56  jim
 * Fixes to make optional force modules work with more nodes than patches.
 *
 * Revision 1.1  1998/04/30 04:53:23  jim
 * Added forces from MDComm and other improvements to ComputeGlobal.
 *
 *
 *
 *
 ***************************************************************************/
