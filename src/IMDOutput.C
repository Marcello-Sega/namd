
#include "IMDOutput.h"
#include "imd.h"
#include "Vector.h"
#include "vmdsock.h"

void IMDupdate(void *v) {
  IMDOutput* imd = (IMDOutput *)v;
  imd->update();
}

IMDOutput::IMDOutput(void *s) {
  sock = s;
}

IMDOutput::~IMDOutput() { }

void IMDOutput::gather_energies(int timestep, BigReal *energies, BigReal T, 
                       BigReal totalEnergy) {
  IMDEnergies nrgs;
  nrgs.Etot = totalEnergy;
  nrgs.Epot = totalEnergy - energies[7];
  nrgs.Ebond = energies[0];
  nrgs.Eangle = energies[1];
  nrgs.Edihe = energies[2];
  nrgs.Eimpr = energies[3];
  nrgs.Eelec = energies[4];
  nrgs.Evdw = energies[5];

  size_t msize = sizeof(IMDEnergies);
  imd_sendheader(sock, IMD_ENERGIES, timestep, msize); 
  int rc = imd_writen(sock, (char *)&nrgs, msize); 
  if (rc != msize) 
    iout << iWARN << "Error sending energies to VMD\n" << endi; 
}

void IMDOutput::gather_coordinates(int timestep, int N, Vector *coords) {
  float *fcoords = new float[3*N];
  for (int i=0; i<N; i++) {
    fcoords[3*i] = coords[i].x;     // If only these coordinates were pointers
    fcoords[3*i+1] = coords[i].y;
    fcoords[3*i+2] = coords[i].z;
  }
  
  size_t msize = 3*N*sizeof(float); 
  imd_sendheader(sock, IMD_FCOORDS, N, msize);
  int rc = imd_writen(sock, (char *)fcoords, msize);
  delete [] fcoords; 
  if (rc != msize) 
    iout << iWARN << "Error sending coordinates to VMD\n" << endi; 
}

void IMDOutput::update() { }

