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

IMDOutput::~IMDOutput() {}

void IMDOutput::gather_energies(IMDEnergies *energies) { 
  if (!sock || !vmdsock_selwrite(sock,0)) 
    return; 
  if (energies->tstep % transrate) return;
  imd_send_energies(sock, energies);
}

void IMDOutput::gather_coordinates(int timestep, int N, FloatVector *coords) {
  if (!sock || !vmdsock_selwrite(sock,0)) return;
  if (timestep % transrate) return;
  imd_send_fcoords(sock, N, (float *)coords);
}

void IMDOutput::update() { }

