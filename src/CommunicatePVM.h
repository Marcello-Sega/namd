/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: CommunicatePVM.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/12/06 19:52:20 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * CommunicatePVM - PVM version of the Communicate object.  Allows
 * user to establish id's for available nodes, establish connections, and
 * send/receive data.
 *
 * When created, this object checks to see if other nodes have been
 * spawned yet; if no, this spawns them and confirms everthing is running
 * well.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: CommunicatePVM.h,v $
 * Revision 1.1  1996/12/06 19:52:20  ari
 * Initial revision
 *
 * Revision 1.11  1995/10/06 20:35:34  nelson
 * Added PVMe definition for include of pvm3.h
 *
 * Revision 1.10  95/04/10  11:27:51  11:27:51  nelson (Mark T. Nelson)
 * Changed handling of previously static variables
 * 
 * Revision 1.9  95/03/20  11:54:03  11:54:03  nelson (Mark T. Nelson)
 * Removed numActive flag
 * 
 * Revision 1.8  95/03/18  02:42:46  02:42:46  nelson (Mark T. Nelson)
 * Reworked extensively to improve performance
 * 
 * Revision 1.7  95/03/08  14:46:19  14:46:19  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.6  94/12/14  16:13:14  16:13:14  nelson (Mark T. Nelson)
 * Added get_tids for FMA interface
 * 
 * Revision 1.5  94/07/26  16:52:28  16:52:28  billh (Bill Humphrey)
 * (Hopefully) fixed problem with broadcast ... now just loops through nodes
 * and sends copy of message to each node, instead of having specific
 * broadcast routine in child class.
 * 
 * Revision 1.4  94/07/03  00:15:02  00:15:02  billh (Bill Humphrey)
 * New format ... user creates a Message object, and given that to a
 * Communicate object; these can be kept in a list to be all sent at once
 * or sent as soon as requested.
 * 
 ***************************************************************************/

#ifndef COMMUNICATEPVM_H
#define COMMUNICATEPVM_H

#ifndef PVMe
#include <pvm3.h>
#else
#include "/usr/include/pvm3.h"		//  PVMe include file is in /usr/include
#endif  //  PVMe

#include "Communicate.h"
#include "Message.h"


class CommunicatePVM : public Communicate {

private:
  char *execName; // name of executable (without path)
  int maxHosts;	 // total number of hosts available for running jobs
  int myHost;	 // my node number (not tid, but index into tid list)
  int *tids;	 // task ID's for the tasks used by this application

  // take data from the given Message, and pack it into the current send buf
  void pack_message(Message *, int);

  Message *unpack_message(int &tag, int &node);

protected:
  //
  // protected virtual functions
  //
  
  virtual int do_send_queue(int);
  virtual int do_send_msg(Message *, int, int, int delmsg=TRUE);

  // check for a message from the given node and tag.  Return a NEW
  // Message object if a message arrives, or NULL if no message available.
  // If another message is found while searching for the given message, add
  // it to the received queue, and keep looking.
  // If node = (-1), checks for messages from any node.
  // If tag = (-1), checks for messages with any tag.
  virtual Message *do_receive(int &, int &);

public:
  // constructor and destructor
  // arguments: command line options, debugging flag 
  CommunicatePVM(int, char *[], int = 0);
  virtual ~CommunicatePVM(void);

  //
  // public virtual routines
  //
  
  // add a node to the list of communications.  Return node number if OK,
  //  (-1) if problem.
  virtual int add_node(void *id);

  virtual int *get_tids();
  
};

#endif
