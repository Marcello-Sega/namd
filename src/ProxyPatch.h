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

class ProxyPatch : public Patch
{
  public:

     ProxyPatch(PatchID pd);
     virtual ~ProxyPatch(void) { };

     void receiveAtoms(ProxyAtomsMsg*);
     void receiveData(ProxyDataMsg*);

  protected:

     virtual void boxClosed(int);

  private:

     void sendResults(void);
     ProxyDataMsg* msgBuffer;

};


#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyPatch.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/12/17 22:13:22 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyPatch.h,v $
 * Revision 1.2  1996/12/17 22:13:22  jim
 * implemented ProxyDataMsg use
 *
 * Revision 1.1  1996/12/05 01:44:16  ari
 * Initial revision
 *
 *
 ***************************************************************************/
