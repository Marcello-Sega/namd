/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/HomePatch.C,v 1.14 1996/12/17 08:55:25 jim Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "HomePatch.h"
#include "AtomMap.h"
#include "main.h"
#include "ProxyMgr.top.h"
#include "ProxyMgr.h"

#define MIN_DEBUG_LEVEL 3
#define DEBUGM
#include "Debug.h"

HomePatch::HomePatch(PatchID pd, AtomIDList al, 
  PositionList pl, VelocityList vl)
  : Patch(pd,al,pl), pInit(pl), v(vl) {
    if (atomIDList.size() != v.size()) {
      CPrintf(
      "HomePatch::HomePatch(...) : Different numbers of Velocities and IDs!\n");
    }
    AtomMap::Object()->registerIDs(pd,al);
}



HomePatch::~HomePatch()
{
}


void HomePatch::boxClosed(int)
{
  if ( ! --boxesOpen )
  {
    DebugM(2,patchID << ": " << "Trying to awaken sequencer.\n");
    sequencer->awaken();
  }
  else
  {
    DebugM(2,patchID << ": " << boxesOpen << " boxes left to close.\n");
  }
}

void HomePatch::registerProxy(RegisterProxyMsg *msg) {
  proxy.add(ProxyListElem(msg->node,forceBox.checkOut()));
  ProxyAtomsMsg *nmsg = new (MsgIndex(ProxyAtomsMsg)) ProxyAtomsMsg;
  nmsg->patch = patchID;
  nmsg->atomIDList = atomIDList;
  ProxyMgr::Object()->sendProxyAtoms(nmsg,msg->node);
}

void HomePatch::unregisterProxy(UnregisterProxyMsg *msg) {
  int i = proxy.findIndex(ProxyListElem(msg->node));
  forceBox.checkIn(proxy[i].forceBox);
  proxy.del(i);
}

void HomePatch::addForceToMomentum(const BigReal timestep)
{
  const BigReal dt = TIMEFACTOR * timestep;
  for ( int i = 0; i < numAtoms; ++i )
  {
    v[i] += f[i] * ( dt / a[i].mass );
  }
}

void HomePatch::addVelocityToPosition(const BigReal timestep)
{
  const BigReal dt = TIMEFACTOR * timestep;
  for ( int i = 0; i < numAtoms; ++i )
  {
    p[i] += v[i] * dt;
  }
}







// direct local calculations

/*
void HomePatch::update_f_at_cycle_begin()
{
     // f = a[0] * f_long + f_short
}



void HomePatch::update_f_at_cycle_end()
{
     // f = f + b*f_long 
}



void HomePatch::update_f_at_step(int k)
{
    // f = a[k] * f_long + f_short
}



void HomePatch::advance_x()
{
     // integrate();
     // shake();

     // positions changed.
     // proxies must be updated
     // inform proxy communicator
}



void HomePatch::update_v()
{
     // v = (v+0.5*dt*f/m)/(1+0.5dt*b) 
}





// callback sequencer when the concurrent computations are done

void HomePatch::f_short_done()
{
    // f_short = f_bonded + f_elshort
    CthAwaken(sequencer); 
}



void HomePatch::f_long_done()
{
    // f_long = f_elfull - f_elshort
    CthAwaken(sequencer); 
}




// trigger concurrent computations

void HomePatch::compute_f_short()
{
     // inform bonded force objects;

     inform_bondedForce(); 

     // inform short_range electrostatic force object
     inform_elShortForce();
     CthSuspend(); 

}


void HomePatch::compute_f_long()
{
     // inform full-electrostatic object
     inform_elFullForce();

     // inform short_range electrostatic force object if not informed
     inform_bondedForce();
     CthSuspend(); 
}






// misc bookeeping methods that require romote information
// such as atom redistribution at the end of teh cycle

void HomePatch::prepare_for_next_cycle()
{
     // atom redistribution

     // initiate atom redistribution

     // figure out atoms movoing out

     // wait for confirmation from PatchManager
     CthSuspend();
}




// this function is invoked to inform the patch that atom redistribution
// phase is complete with new atoms moving in (if any)
void HomePatch::atom_redist_data(int n, int *newAtoms )
{
   // adjust internal data structure


   // call back the sequencer object so that the next cycle resumes
   CthAwaken(sequencer); 
}




void HomePatch::prepare_for_next_step(void)
{
   CthSuspend();
   CthAwaken(sequencer); 
}

*/


/* ************************************************************************ */
/* unpack the patch data that is sent to this patch                         */
/* "data" points the beginning of incoming data.                            */
/* return the adjusted pointer (which points to the next data block)        */
/*                                                                          */
/* a home patch receives forces back                                        */
/* ************************************************************************ */

/*
void HomePatch::updateData(char *&data)
{
}
*/





/* ************************************************************************ */
/* pack the patch data that is gonna be sent to the proxy of this home patch*/
/* at the beginning of a cycle (that is include atom list                   */
/* "data" points the beginning of area where patch data to be put.          */
/* a home patch sends atom list, x, and xbegin                              */
/* ************************************************************************ */

/*
void HomePatch::packInitData(char *data)
{
}
*/



/* ************************************************************************ */
/* pack the patch data that is gonna be sent to the proxy of this home patch*/
/* in the mid cycle (that is dont sent atom list           j                */
/* a home patch sends only x at mid cycle x                                 */
/* ************************************************************************ */
/*
void HomePatch::packData(char *&data)
{
}
*/



/*
void HomePatch::dispose(char *&data)
{
}
*/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: HomePatch.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.14 $	$Date: 1996/12/17 08:55:25 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HomePatch.C,v $
 * Revision 1.14  1996/12/17 08:55:25  jim
 * added node argument to sendProxyAtoms
 *
 * Revision 1.13  1996/12/16 22:52:43  jim
 * added placement new and explicit destructor calls to ProxyAtomsMsg
 *
 * Revision 1.12  1996/12/14 00:02:42  jim
 * debugging ProxyAtomsMsg path to make compute creation work
 *
 * Revision 1.11  1996/12/11 22:31:41  jim
 * added integration methods for Sequencer
 *
 * Revision 1.10  1996/12/05 23:45:09  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/12/05 01:44:16  ari
 * started toward proxy management
 *
 * Revision 1.8  1996/12/01 02:35:47  jim
 * added patchID reporting
 *
 * Revision 1.7  1996/11/30 00:35:51  jim
 * implemented boxClosed(), useSequencer(), runSequencer()
 *
 * Revision 1.6  1996/11/22 01:44:53  jim
 * added calls to service AtomMap
 *
 * Revision 1.5  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/10/04 21:14:47  jim
 * Moved functionality to Patch
 *
 * Revision 1.3  1996/09/03 22:54:09  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.7  1996/07/15 21:18:15  gursoy
 * *** empty log message ***
 *
 * Revision 1.6  1996/07/10 16:25:29  gursoy
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/10 16:19:52  gursoy
 * *** empty log message ***
 *
 * Revision 1.4  1996/07/02 15:10:27  gursoy
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/11 22:36:23  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/11 20:07:22  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/05/30 21:31:36  gursoy
 * Initial revision
 *
 ***************************************************************************/
