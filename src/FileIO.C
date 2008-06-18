/**
***  Copyright (c) 1995-2007 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
#include "FileIO.decl.h"
#include "FileIO.h"

FileIO::FileIO() { numAtomsToWriteAvg = 1; }
FileIO::FileIO(CkMigrateMessage* msg) { }

/* body here */
void FileIO::writeDCD(int base) 
{
  // calculate our offsets
  int xoff=numAtomsToWriteAvg*CkMyPe();
  int yoff=xoff*2;
  int zoff=xoff*3;
  xoff+=base;
  yoff+=base;
  zoff+=base;
  
  // force sending of coordinates for this proc.
  
  // smarter way is to have each patch do its own output.
  // why bother sending this info anywhere?
  // seeking around solves sequential ordering of output by
  // atom number
  
  // the only global issue is the clusters so we need to know which
  // parts of split clusters have to get their coordinates shifted to
  // match the unit box with the idea being that we don't break a
  // cluster over the box boundary.

}

void FileIO::doneDCD(CkReductionMsg *msg) {}
void FileIO::recvCoordinates() 
{
  // if array not allocated make it
  // copy coordinates into local array

  // when we have everything do box corrections

  // write atoms

  //contribute to reduction
}


#include "FileIO.def.h"
