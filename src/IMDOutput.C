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
  int numatoms = Node::Object()->molecule->numAtoms;
  imd_sendheader(sock, IMD_HANDSHAKE, numatoms, sizeof(int));
  int one = 1;
  imd_writen(sock, (char *)&one, sizeof(int));
  //fcoords = new float[3*numatoms]; 
  transrate = 1;
}

IMDOutput::~IMDOutput() { 
  //delete [] fcoords;
}

void IMDOutput::gather_energies(int timestep, IMDEnergies *energies) { 
  if (!sock) return; 
  if (timestep % transrate) return;
  size_t msize = sizeof(IMDEnergies);
  imd_sendheader(sock, IMD_ENERGIES, timestep, msize); 
  int rc = imd_writen(sock, (char *)energies, msize); 
  if (rc != msize) 
    iout << iWARN << "Error sending energies to VMD\n" << endi; 
}

void IMDOutput::gather_coordinates(int timestep, int N, FloatVector *coords) {
  if (!sock) return;
  if (timestep % transrate) return;
/*
  for (int i=0; i<N; i++) {
    fcoords[3*i] = coords[i].x;     // If only these coordinates were floats 
    fcoords[3*i+1] = coords[i].y;
    fcoords[3*i+2] = coords[i].z;
  }
*/
  size_t msize = 3*N*sizeof(float); 
  imd_sendheader(sock, IMD_FCOORDS, N, msize);
  //int rc = imd_writen(sock, (char *)fcoords, msize);
  int rc = imd_writen(sock, (char *)coords, msize);
}

void IMDOutput::update() { }

