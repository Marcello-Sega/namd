/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GLOBALMASTERIMD_H
#define GLOBALMASTERIMD_H

class GlobalMasterIMD : public GlobalMaster {
 public: 
  /* initializes this according to the simulation parameters */
  GlobalMasterIMD();
  ~GlobalMasterIMD();

 protected:

  virtual void calculate();

  // Simple function for getting MDComm-style forces from VMD
  int get_vmd_forces();

  // flag for whether to proceed with simulation when there are no connections
  int IMDwait;

  // My socket handle
  void *sock;
  void *clientsock;
};

#endif

