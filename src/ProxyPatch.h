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

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ProxyPatch.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:59:14 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ProxyPatch.h,v $
 * Revision 1.1000  1997/02/06 15:59:14  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:20  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:40  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:55  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.2  1996/12/17 22:13:22  jim
 * implemented ProxyDataMsg use
 *
 * Revision 1.1  1996/12/05 01:44:16  ari
 * Initial revision
 *
 *
 ***************************************************************************/
