/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "vmdsock.h"
#include "Node.h"
#include "IMDOutput.h"
#include "imd.h"
#include "SimParameters.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "ComputeIMD.h"
#include "ComputeMgr.h"
#include "UniqueSortedArray.h"
#include <errno.h>

struct vmdforce {
  int index;
  Vector force;
  int operator <(const vmdforce& v) {return index < v.index;}
  int operator ==(const vmdforce& v) {return index == v.index;}
  vmdforce& operator=(const vmdforce& v) {
    index=v.index;
    force=v.force; 
    return *this;
  }  
};

static UniqueSortedArray<vmdforce> vmdforces;

// Search for a free port in the range 1025-4096; return the successful port,
// or -1 if we failed.

static int find_free_port(void *sock, int defport) {
  if (vmdsock_bind(sock, defport)==0) return defport; // success
  for (int port=1025; port < 4096; port++) 
    if (vmdsock_bind(sock, port)==0) return port;
  return -1;
}
 
ComputeIMD::ComputeIMD(ComputeGlobal *h)
: ComputeGlobalMaster(h) {
  SimParameters *simparams = Node::Object()->simParameters;
  int rc;
  int port = simparams->IMDport;

  sock = vmdsock_create();
  port = find_free_port(sock, port);
  if (port < 0) {
    vmdsock_destroy(sock);
    NAMD_die("Unable to find free port\n");
  }
  iout << iINFO << "INTERACTIVE MD BIND    "<< port << "\n" << endi;
  rc = vmdsock_listen(sock); 
  // Wait for VMD to connect
  rc = vmdsock_selread(sock,3600);
  if (rc < 0) {
    iout << iDEBUG << "select: " << strerror(errno) << '\n' << endi;
    vmdsock_destroy(sock);
    NAMD_die("Connection error\n");
  }
  if (rc == 0) {
    vmdsock_destroy(sock);
    NAMD_die("No connection\n");
  }
  rc = vmdsock_accept(sock);
  if (rc < 0) {
    iout << iDEBUG << "Accept: " << strerror(errno) << '\n' << endi;
    NAMD_die("Connection error\n");
  } 
  Node::Object()->IMDinit(sock);

}

ComputeIMD::~ComputeIMD() {
  if (sock) 
    vmdsock_destroy(sock);
}

void ComputeIMD::initialize() {
  configMsg = new ComputeGlobalConfigMsg;
  host->comm->sendComputeGlobalConfig(configMsg);
  configMsg = 0;
}

void ComputeIMD::calculate() {
  // Assume for now that the only thing we get from VMD is a set of forces.
  // Later we'll want to look for and implement more sophisticated control
  // parameters. 
 
  resultsMsg = new ComputeGlobalResultsMsg;
  resultsMsg->gforce.resize(gmass.size());
  resultsMsg->gforce.setall(Vector(0,0,0));

  // Check/get new forces from VMD
  if (sock) {
    get_vmd_forces();
  }

  // Right now I don't check to see if any new forces were obtained.
  // An optimization would be cache the results message.  However, there
  // would still be copying since it looks like the messages get deleted
  // by the receiver.

    int num = vmdforces.size();
    resultsMsg->aid.resize(num);
    resultsMsg->f.resize(num);
    int i;
    AtomIDList::iterator aid_i = resultsMsg->aid.begin();
    ForceList::iterator f_i = resultsMsg->f.begin();
    UniqueSortedArray<vmdforce>::iterator v_i = vmdforces.begin();
    for ( i = 0; i < num; ++i, ++v_i) {
      aid_i[i] = v_i->index; 
      f_i[i] = v_i->force;
    }
  // Send results to clients
  host->comm->sendComputeGlobalResults(resultsMsg);
  resultsMsg = 0;
}

int ComputeIMD::get_vmd_forces() {
  IMDType type;
  int32 length;
  int32 *vmd_atoms;
  float *vmd_forces;
  int retval; 
  vmdforce *vtest, vnew;

  while (vmdsock_selread(sock,0) > 0)  {     // Drain the socket
    type = imd_recv_header(sock, &length);
    int i;
    switch (type) {
      case IMD_MDCOMM:
        // Expect the msglength to give number of indicies, and the data
        // message to consist of first the indicies, then the coordinates
        // in xyz1 xyz2... format.
        vmd_atoms = new int32[length];
        vmd_forces = new float[3*length];
        if (imd_recv_mdcomm(sock, length, vmd_atoms, vmd_forces)) {
          NAMD_die("Error reading MDComm forces\n");
        } 
        for (i=0; i<length; i++) {
          vnew.index = vmd_atoms[i];
          if ( (vtest=vmdforces.find(vnew)) != NULL) {
            // find was successful, so overwrite the old force values
            vtest->force.x = vmd_forces[3*i];
            vtest->force.y = vmd_forces[3*i+1];
            vtest->force.z = vmd_forces[3*i+2];
          }
          else {
            // Create a new entry in the table
            vnew.force.x = vmd_forces[3*i];
            vnew.force.y = vmd_forces[3*i+1];
            vnew.force.z = vmd_forces[3*i+2];
            vmdforces.add(vnew);
          }
        } 
        retval = 1;
        delete [] vmd_atoms;
        delete [] vmd_forces;
        break;
      case IMD_TRATE:
        iout << iINFO << "Setting transfer rate to " << length<<'\n'<<endi;	
        Node::Object()->imd->set_transrate(length);
        break;
      case IMD_DISCONNECT:
        iout<<iDEBUG<<"Detaching simulation from remote connection\n" << endi;
        vmdsock_destroy(sock);
        sock = 0;
        Node::Object()->imd->close();
        goto vmdEnd;
      case IMD_KILL:
        NAMD_quit(1);
        break;
      case IMD_ENERGIES:
        IMDEnergies junk;
        imd_recv_energies(sock, &junk);
        break;
      case IMD_FCOORDS:
        vmd_forces = new float[3*length];
        imd_recv_fcoords(sock, length, vmd_forces);
        delete [] vmd_forces;
        break;
      case IMD_IOERROR:
        vmdsock_destroy(sock);
        NAMD_die("IMD connection lost\n");
        break;    
      default: ;
    }
  }
vmdEnd:
  return retval;
}

