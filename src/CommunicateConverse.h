//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 * CommunicateConverse - Converse version of the Communicate object.  Allows
 * user to establish id's for available nodes, establish connections, and
 * send/receive data.
 * 
 ***************************************************************************/

#ifndef COMMUNICATECONVERSE_H
#define COMMUNICATECONVERSE_H

extern "C" {
#include "converse.h"
extern void CmiGrabBuffer(char **msg);
extern void ConverseInit(char *argv[]);
extern void ConverseExit(void);
}

#include "Communicate.h"
#include "Message.h"


class CommunicateConverse : public Communicate {

private:
  int CsmHandlerIndex;

  int size_message(Message *);
  char *pack_message(Message *, char *);
  Message *unpack_message(char *, int &tag, int &node);

protected:
  
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
  CommunicateConverse(int, char *[], int = 0);
  virtual ~CommunicateConverse(void);

  // public virtual routines
  
  // add a node to the list of communications.  Return node number if OK,
  //  (-1) if problem.
  virtual int add_node(void *id);

  virtual int *get_tids();
  
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: CommunicateConverse.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/04 22:37:06 $
 *
 ***************************************************************************/

