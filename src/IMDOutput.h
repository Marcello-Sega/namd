
#ifndef IMD_OUTPUT_H__
#define IMD_OUTPUT_H__

#include "common.h"
#include "imd.h"
class Vector;

// IMDOutput
// An object that sits in Node and does all the outputting of coordinates,
// energies, etc.  Its socket connection gets initialized by ComputeIMD. 
// If there's anything that needs to be updated periodically by converse,
// it'll be here.

extern "C" {
  void IMDupdate(void*);
}

class IMDOutput {

public:
  IMDOutput(void *);   // Initialize with socket handle
  ~IMDOutput();

  void gather_energies(int timestep, IMDEnergies *energies); 
  void gather_coordinates(int timestep, int N, Vector *coords);

  void set_transrate(int newrate) {transrate = newrate; }
  void update();
  void close() {sock = 0;}

private:
  void *sock;

  float *fcoords;
  int curstep;   
  int transrate;
  int haveEnergies;
  int haveCoords;
};

#endif

