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
#include "Templates/Queue.h"

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
     ShortQueue<ProxyDataMsg*> msgBuffer;

};


#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyPatch.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/12/05 01:44:16 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyPatch.h,v $
 * Revision 1.1  1996/12/05 01:44:16  ari
 * Initial revision
 *
 *
 ***************************************************************************/
