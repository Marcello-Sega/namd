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
#include "UniqueSortedArray.h"
#include "GlobalMaster.h"
#include "GlobalMasterIMD.h"
#include "Vector.h"
#include "InfoStream.h"

#include <errno.h>

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

struct vmdforce {
  int index;
  Vector force;
  int operator <(const vmdforce& v) {return index < v.index;}
  // XXX the following is an abuse of overloading!
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
 
GlobalMasterIMD::GlobalMasterIMD() {
  DebugM(3,"Constructing\n");
  SimParameters *simparams = Node::Object()->simParameters;
  int port = simparams->IMDport;
  IMDwait = simparams->IMDwait;
  IMDignore = simparams->IMDignore;
  coordtmp = NULL;
  coordtmpsize = 0;

  if ( vmdsock_init() ) {
    NAMD_die("Unable to initialize socket interface for IMD.\n");
  }
  sock = vmdsock_create();
  clientsock = NULL;
  int newport = find_free_port(sock, port);
  if (newport != port) {
    iout << iWARN << "Interactive MD failed to bind to port "
                  << port << ".\n" << endi;
  }
  if (newport < 0) {
    vmdsock_destroy(sock);
    NAMD_die("Interactive MD failed to find free port.\n");
  }
  vmdsock_listen(sock); 
  iout << iINFO << "Interactive MD listening on port "
                  << newport << ".\n" << endi;
  DebugM(2,"Done constructing ("<<requestedGroups().size()<<" initial groups)\n");

  Node::Object()->imd->use_imd(this);
}

GlobalMasterIMD::~GlobalMasterIMD() {
  if (sock) 
    vmdsock_destroy(sock);
  if (clientsock)
    vmdsock_destroy(clientsock);
  delete [] coordtmp;
}

static int my_imd_connect(void *s) {
  if (imd_handshake(s)) {
    iout << iWARN << "IMD handshake failed\n" << endi;
    return 0;
  }

  // Wait a second, then see if VMD has responded.
  double t = CmiWallTimer();
  while (CmiWallTimer()-t < 1.0);
  int32 length;
  if (vmdsock_selread(s,0) != 1 || imd_recv_header(s, &length) != IMD_GO) {
    iout << iWARN << "Incompatible Interactive MD, use VMD v1.4b2 or higher\n"
         << endi;
    return 0;
  }
  return 1;
}

void GlobalMasterIMD::calculate() {
  /* clear out the requested forces first! */
  if (!IMDignore) {
    modifyAppliedForces().resize(0);
    modifyForcedAtoms().resize(0);
    modifyGroupForces().resize(0);
  }

  if (!clientsock) {
    // check for incoming connection
    int rc;
    if (IMDwait) {
      iout << iINFO << "INTERACTIVE MD AWAITING CONNECTION\n" << endi;
      do { rc = vmdsock_selread(sock, 3600); } while (rc <= 0);
    } else {
      rc = vmdsock_selread(sock, 0);
    } 
    if (rc > 0) {
      clientsock = vmdsock_accept(sock);
      if (!clientsock) {
        iout << iWARN << "GlobalMasterIMD: Accept failed\n" << endi;
      } else {
        if (!my_imd_connect(clientsock)) {
          iout << iWARN << "GlobalMasterIMD: IMD handshake failed\n" << endi;
          vmdsock_destroy(clientsock);
          clientsock = NULL;
        } else {
          // success 
          Node::Object()->imd->set_transrate(1);
        }
      }
    }
  }

  // Assume for now that the only thing we get from VMD is a set of forces.
  // Later we'll want to look for and implement more sophisticated control
  // parameters. ie have a specified protocol

  // Check/get new forces from VMD
  if (clientsock) {
    get_vmd_forces();
  }

  // Right now I don't check to see if any new forces were obtained.
  // An optimization would be cache the results message.  However, there
  // would still be copying since it looks like the messages get deleted
  // by the receiver.

  /* set our arrays to be big enough to hold all of the forces */
  int num = vmdforces.size();

  DebugM(2,"Setting " << num << " forces.\n");
  
  if (!IMDignore) {
    modifyForcedAtoms().resize(num);
    modifyAppliedForces().resize(num);
  
    int i;
    UniqueSortedArray<vmdforce>::iterator v_i = vmdforces.begin();
    for ( i = 0; i < num; ++i, ++v_i) {
      modifyForcedAtoms().item(i) = v_i->index;
      modifyAppliedForces().item(i) = v_i->force;
    }
  }
}

void GlobalMasterIMD::get_vmd_forces() {
  IMDType type;
  int32 length;
  int32 *vmd_atoms;
  float *vmd_forces;
  int paused = 0;
  vmdforce *vtest, vnew;

  while (vmdsock_selread(clientsock,0) > 0 || paused) {  // Drain the socket
    type = imd_recv_header(clientsock, &length);
    int i;
    switch (type) {
      case IMD_MDCOMM:
        // Expect the msglength to give number of indicies, and the data
        // message to consist of first the indicies, then the coordinates
        // in xyz1 xyz2... format.
        vmd_atoms = new int32[length];
        vmd_forces = new float[3*length];
        if (imd_recv_mdcomm(clientsock, length, vmd_atoms, vmd_forces)) {
          NAMD_die("Error reading MDComm forces\n");
        } 
        if (IMDignore) {
          iout << iWARN << "Ignoring IMD forces due to IMDignore\n" << endi;
        } else {
          for (i=0; i<length; i++) {
            vnew.index = vmd_atoms[i];
            if ( (vtest=vmdforces.find(vnew)) != NULL) {
              // find was successful, so overwrite the old force values
              if (vmd_forces[3*i] != 0.0f || vmd_forces[3*i+1] != 0.0f
                  || vmd_forces[3*i+2] != 0.0f) {
                vtest->force.x = vmd_forces[3*i];
                vtest->force.y = vmd_forces[3*i+1];
                vtest->force.z = vmd_forces[3*i+2];
              } else {
                // or delete it from the list if the new force is ZERO
                vmdforces.del(vnew);
              }
            }
            else {
              // Create a new entry in the table if the new force isn't ZERO
              if (vmd_forces[3*i] != 0.0f || vmd_forces[3*i+1] != 0.0f
                  || vmd_forces[3*i+2] != 0.0f) {
                vnew.force.x = vmd_forces[3*i];
                vnew.force.y = vmd_forces[3*i+1];
                vnew.force.z = vmd_forces[3*i+2];
                vmdforces.add(vnew);
              }
            }
          } 
        }
        delete [] vmd_atoms;
        delete [] vmd_forces;
        break;
      case IMD_TRATE:
        iout << iINFO << "Setting transfer rate to " << length<<'\n'<<endi;	
        Node::Object()->imd->set_transrate(length);
        break;
      case IMD_PAUSE:
        if (IMDignore) {
          iout << iWARN << "Ignoring IMD pause due to IMDignore\n" << endi;
          break;
        }
        if ( paused ) {
          iout << iINFO << "Resuming IMD\n" << endi;
          IMDwait = Node::Object()->simParameters->IMDwait;
        }
        paused = ! paused;
        if ( paused ) {
          iout << iINFO << "Pausing IMD\n" << endi;
          IMDwait = 1;
        }
        break;
      case IMD_IOERROR:
        iout << iWARN << "IMD connection lost\n" << endi;
      case IMD_DISCONNECT:
        iout<<iDEBUG<<"Detaching simulation from remote connection\n" << endi;
        vmdsock_destroy(clientsock);
        clientsock = 0;
        goto vmdEnd;
      case IMD_KILL:
        if (IMDignore) {
          iout << iWARN << "Ignoring IMD kill due to IMDignore\n" << endi;
          break;
        }
        NAMD_quit(1);
        break;
      case IMD_ENERGIES:
        IMDEnergies junk;
        imd_recv_energies(clientsock, &junk);
        break;
      case IMD_FCOORDS:
        vmd_forces = new float[3*length];
        imd_recv_fcoords(clientsock, length, vmd_forces);
        delete [] vmd_forces;
        break;
      default: ;
    }
  }
vmdEnd: ;
}

void GlobalMasterIMD::send_energies(IMDEnergies *energies) {
  if (!clientsock || !vmdsock_selwrite(clientsock,0)) return;
  imd_send_energies(clientsock, energies);
}

void GlobalMasterIMD::send_fcoords(int N, FloatVector *coords) {
  if (!clientsock || !vmdsock_selwrite(clientsock,0)) return;
  if (sizeof(FloatVector) == 3*sizeof(float)) {
    imd_send_fcoords(clientsock, N, (float *)coords);
  } else {
    if (coordtmpsize < N) {
      delete [] coordtmp;
      coordtmp = new float[3*N];
      coordtmpsize = N;
    }
    for (int i=0; i<N; i++) {
      coordtmp[3*i] = coords[i].x; 
      coordtmp[3*i+1] = coords[i].y; 
      coordtmp[3*i+2] = coords[i].z; 
    } 
    imd_send_fcoords(clientsock, N, coordtmp);
  }
}
