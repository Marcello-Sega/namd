/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "IMDOutput.h"
#include "InfoStream.h"
#include "Molecule.h"
#include "Node.h"
#include "SimParameters.h"
#include "Vector.h"
#include "vmdsock.h"
#include "imd.h"

void IMDupdate(void *v) {
  IMDOutput* imd = (IMDOutput *)v;
  imd->update();
}

IMDOutput::IMDOutput() {
  sock = NULL;
  serversock = NULL;
  coordtmp = NULL;
  coordtmpsize = 0;
}

int IMDOutput::connect(void *s) {
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
  sock = s;
  transrate = 1;
  return 1;
}

void IMDOutput::disconnect() {
  sock = NULL;
}

IMDOutput::~IMDOutput() {
  delete [] coordtmp;
  if (serversock) vmdsock_destroy(serversock);
  if (sock) vmdsock_destroy(sock);
}

// used only with IMDignore option
void IMDOutput::manage_sockets() {
  if (!sock || !vmdsock_selwrite(sock,0)) {
    if ( ! serversock ) {
      SimParameters *simParams = Node::Object()->simParameters;
      if ( ! simParams->IMDignore ) return;
      if ( vmdsock_init() ) {
        iout << iWARN << "Unable to initialize socket interface for IMD.\n"
                      << endi;
        return;
      }
      serversock = vmdsock_create();
      int port = simParams->IMDport;
      if ( vmdsock_bind(serversock, port) ) {
        iout << iWARN << "Interactive MD failed to bind to port "
                        << port << ".\n" << endi;
        for (port=1025; port<4096 && vmdsock_bind(serversock,port); port++);
        if ( port < 4096 ) {
          iout << iWARN << "Interactive MD listening on port "
                        << port << ".\n" << endi;
        } else {
          vmdsock_destroy(serversock);
          serversock = NULL;
          iout << iWARN << "Interactive MD failed to find free port.\n" << endi;
        }
      }
      if ( serversock ) vmdsock_listen(serversock);
    }
    // check for incoming connection
    int rc;
    rc = vmdsock_selread(serversock, 0);
    if (rc > 0) {
      sock = vmdsock_accept(serversock);
      if (!sock) {
        iout << iWARN << "IMDOutput: Accept failed\n" << endi;
      } else {
        if (! connect(sock)) {
          iout << iWARN << "IMDOutput: IMD handshake failed\n" << endi;
          vmdsock_destroy(sock);
          sock = NULL;
        }
      }
    }
  }

  IMDType type;
  int32 length;
  int32 *vmd_atoms;
  float *vmd_forces;

  while (sock && vmdsock_selread(sock,0) > 0) {  // Drain the socket
    type = imd_recv_header(sock, &length);
    int i;
    switch (type) {
      case IMD_TRATE:
        iout << iINFO << "Setting transfer rate to " << length<<'\n'<<endi;
        set_transrate(length);
        break;
      case IMD_PAUSE:
        iout << iWARN << "IMD ignoring pause command.\n" << endi;
        break;
      case IMD_MDCOMM:
        // Expect the msglength to give number of indicies, and the data
        // message to consist of first the indicies, then the coordinates
        // in xyz1 xyz2... format.
        vmd_atoms = new int32[length];
        vmd_forces = new float[3*length];
        i = imd_recv_mdcomm(sock, length, vmd_atoms, vmd_forces);
        delete [] vmd_atoms;
        delete [] vmd_forces;
        if ( ! i ) {
          iout << iWARN << "IMD ignoring forces.\n" << endi;
          break;
        }
        iout<<iWARN<<"Error reading IMD forces.\n";
      case IMD_IOERROR:
        iout << iWARN << "IMD connection lost.\n" << endi;
      case IMD_DISCONNECT:
      case IMD_KILL:
        iout<<iINFO<<"Detaching simulation from remote connection.\n" << endi;
        vmdsock_destroy(sock);
        sock = 0;
        disconnect();
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
      default: ;
    }
  }
}

void IMDOutput::gather_energies(IMDEnergies *energies) { 
  if (Node::Object()->simParameters->IMDignore) manage_sockets();
  if (!sock || !vmdsock_selwrite(sock,0)) return;
  if (energies->tstep % transrate) return;
  imd_send_energies(sock, energies);
}

void IMDOutput::gather_coordinates(int timestep, int N, FloatVector *coords) {
  if (!sock || !vmdsock_selwrite(sock,0)) return;
  if (timestep % transrate) return;
  if (sizeof(FloatVector) == 3*sizeof(float)) {
    imd_send_fcoords(sock, N, (float *)coords);
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
    imd_send_fcoords(sock, N, coordtmp);
  }
}

void IMDOutput::update() { }

