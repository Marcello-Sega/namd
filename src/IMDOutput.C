/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "IMDOutput.h"
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

IMDOutput::IMDOutput(void *s) {
  sock = s;
  imd_handshake(sock);
  transrate = 1;
}

IMDOutput::~IMDOutput() { 
  //delete [] fcoords;
}

void IMDOutput::gather_energies(IMDEnergies *energies) { 
  if (!sock) return; 
  if (energies->tstep % transrate) return;
  imd_send_energies(sock, energies);
}

void IMDOutput::gather_coordinates(int timestep, int N, FloatVector *coords) {
  if (!sock) return;
  if (timestep % transrate) return;
  imd_send_fcoords(sock, N, (float *)coords);
}

void IMDOutput::update() { }

