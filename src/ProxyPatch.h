/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PROXYPATCH_H
#define PROXYPATCH_H

#include "Patch.h"

class ProxyDataMsg;
class ProxyAtomsMsg;
class ProxyAllMsg;

class ProxyPatch : public Patch
{
  public:

     ProxyPatch(PatchID pd);
     virtual ~ProxyPatch(void);

     void receiveAtoms(ProxyAtomsMsg*);
     void receiveData(ProxyDataMsg*);
     void receiveAll(ProxyAllMsg*);

  protected:

     virtual void boxClosed(int);

  private:

     void sendResults(void);
     ProxyDataMsg* msgBuffer;
     ProxyAllMsg* msgAllBuffer;

};


#endif

