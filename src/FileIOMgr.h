/**
***  Copyright (c) 1995-2007 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/** This is a group because it would be silly to have more writers
    than there are processors. */

/** Each FileIO object wraps a sequential Output object.
    This allows us to reuse the existing sequential code.

    FileIO handles the coordination of 
 */ 

#ifndef _FILEIO_H
#define _FILEIO_H

#include "charm++.h"


class FileIOMgr: public Group
{
 public:
  FileIO(int totalAtms, int numAtoms, int myAtoms)
    {
      totalAtoms=totalAtms;
      numAtomsToWriteAvg=numAtoms;
      myNumAtomsToWrite= myAtoms;
    }
  void writeDCD(int base);
  void doneDCD(CkReductionMsg *msg);
  void recvCoordinates();
 private:
  int totalAtoms;                  //! total number of atoms 
  int numAtomsToWriteAvg;          //! typical number per writer
  int myNumAtomsToWrite;           //! exact number for this writer  
};

#endif /* FILEIO_H*/