//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

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
     virtual ~ProxyPatch(void) { };

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

