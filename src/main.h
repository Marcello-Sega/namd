//-*-c++-*-
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

#ifndef MAIN_H
#define MAIN_H

#include "charm++.h"

#include "BOCgroup.h"
#include "NamdTypes.h"

#include "main.decl.h"

// Message to send our per processor BOC's list of groupIDs of 
// all other BOC's
class GroupInitMsg : public CMessage_GroupInitMsg
{
public:
  BOCgroup group;
};

class SlaveInitMsg : public GroupInitMsg
{
public:
  CkChareID master;
};

class Compute;

// For Compute objects to enqueue themselves when ready to compute
class LocalWorkMsg : public CMessage_LocalWorkMsg
{
public:
  Compute *compute;
};

#endif /* MAIN_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: main.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1004 $	$Date: 1999/05/11 23:56:57 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.h,v $
 * Revision 1.1004  1999/05/11 23:56:57  brunner
 * Changes for new charm version
 *
 * Revision 1.1003  1998/03/03 23:05:32  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1002  1997/12/26 23:11:08  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1001  1997/02/26 16:53:22  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1000  1997/02/06 15:59:37  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:39  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:46  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:37:14  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.20  1997/01/13 19:18:37  jim
 * added init message for slave/master pattern
 *
 * Revision 1.19  1996/12/17 17:07:41  jim
 * moved messages from main to ProxyMgr
 *
 * Revision 1.18  1996/12/17 08:54:40  jim
 * fixed several bugs, including not saving patch
 *
 * Revision 1.17  1996/12/16 22:52:43  jim
 * added placement new and explicit destructor calls to ProxyAtomsMsg
 *
 * Revision 1.16  1996/12/16 22:19:26  jim
 * added placement new and destructor to messages
 *
 * Revision 1.15  1996/12/14 00:02:42  jim
 * debugging ProxyAtomsMsg path to make compute creation work
 *
 * Revision 1.14  1996/12/05 20:26:27  jim
 * missing ;
 *
 * Revision 1.13  1996/12/05 20:25:03  jim
 * forgot to fill in length of message
 *
 * Revision 1.12  1996/12/05 17:56:58  jim
 * added pack and unpack functions for proxy messages
 *
 * Revision 1.11  1996/12/05 17:15:23  jim
 * added proxy messages
 *
 * Revision 1.10  1996/12/05 17:00:05  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/12/05 01:47:40  ari
 * added messages for proxy management
 *
 * Revision 1.8  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/09/15 20:43:11  jim
 * fixed missing semicolon
 *
 * Revision 1.6  1996/09/10 04:45:06  ari
 * added LocalWorkMsg
 *
 * Revision 1.5  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/22 16:15:07  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/19 22:05:31  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 20:52:30  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/02 19:20:13  gursoy
 * Initial revision
 *
 ***************************************************************************/
