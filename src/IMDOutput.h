/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef IMD_OUTPUT_H__
#define IMD_OUTPUT_H__

#include "common.h"
#include "imd.h"
class FloatVector;

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

  void gather_energies(IMDEnergies *energies); 
  void gather_coordinates(int timestep, int N, FloatVector *coords);

  void set_transrate(int newrate) {transrate = newrate; }
  void update();
  void close() {sock = 0;}

private:
  void *sock;

  //float *fcoords;
  int curstep;   
  int transrate;
  int haveEnergies;
  int haveCoords;
};

#endif

