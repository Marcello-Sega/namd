//-*-c++-*-
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
 *	$RCSfile: Communicate.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.778 $	$Date: 1997/01/28 00:29:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * Communicate - base communications object.  Allows user to establish id's
 * for available nodes, establish connections, and send/receive data.
 * Each Communicate object contains a list of current sending and received
 * messages.  To send a message, a Message object must be created and loaded
 * with data, which is then sent via a call to send(msg, node, tag); to receive
 * a message, one calls Message *receive(node, tag), which returns a message
 * if available, or NULL otherwise.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Communicate.h,v $
 * Revision 1.778  1997/01/28 00:29:54  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:44:53  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:35:27  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.4  1996/12/06 19:52:20  ari
 * *** empty log message ***
 *
 * Revision 1.14  1996/01/28 21:50:08  jean
 * Attempting to make stable RCS without Mark Nelson's
 * fma/pairlist decoupling
 *
 * Revision 1.15  1995/12/04 20:46:19  brunner
 *  More communication stats - message sizes
 *
 * Revision 1.14  1995/11/22 12:07:06  brunner
 * Added print_comm_stats(), and variables to store message traffic
 * information.  Activated by -DCHEAP_PERFORMANCE
 *
 * Revision 1.13  95/10/06  17:53:51  17:53:51  hazen (Brett Hazen)
 * Memory Allocation error-checking added
 * 
 * Revision 1.12  1995/03/20  11:44:27  nelson
 * Made minor modification to send_method for 1 node case
 *
 * Revision 1.11  95/03/20  09:26:27  09:26:27  nelson (Mark T. Nelson)
 * Added single node checks for send_all
 * 
 * Revision 1.10  95/03/18  02:41:50  02:41:50  nelson (Mark T. Nelson)
 * Reworked extensively to improve performance
 * 
 * Revision 1.9  95/03/08  14:37:41  14:37:41  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.8  94/12/14  16:12:38  16:12:38  nelson (Mark T. Nelson)
 * Added get_tids virtual function
 * 
 * Revision 1.7  94/10/24  09:25:07  09:25:07  nelson (Mark T. Nelson)
 * Added MessageManagers to make looking for received messages much
 * more efficient
 * 
 * Revision 1.6  94/07/26  16:52:24  16:52:24  billh (Bill Humphrey)
 * (Hopefully) fixed problem with broadcast ... now just loops through nodes
 * and sends copy of message to each node, instead of having specific
 * broadcast routine in child class.
 * 
 * Revision 1.5  94/07/03  01:22:53  01:22:53  billh (Bill Humphrey)
 * Made Communicate enum's and MsgItem struct public members of
 * Communicate, instead of global.
 * 
 * Revision 1.4  94/07/03  00:14:14  00:14:14  billh (Bill Humphrey)
 * New format ... user creates a Message object, and given that to a
 * Communicate object; these can be kept in a list to be all sent at once
 * or sent as soon as requested.
 * 
 ***************************************************************************/

#ifndef COMMUNICATE_OBJ_H
#define COMMUNICATE_OBJ_H

#include "Message.h"
#include "IntTree.h"

class MessageManager;
class MessageQueue;

class Communicate 
{

public:
  // send method
  enum SendMethod { NOW, WAIT };

  // error codes
  enum CommError { NOERROR, ERROR, NONODES, NOSEND, NORECEIVE };


private:
  // print debugging messages flag
  int Debug;

  // sending method ... either send all messages as soon as they arrive,
  // or keep them until a signal is given to send all at once
  SendMethod sendMethod;

protected:
#ifdef CHEAP_PERFORMANCE
  unsigned long num_local_sends;
  unsigned long num_local_recvs;
  
  unsigned long num_remote_sends;
  unsigned long num_remote_recvs;
  unsigned long remote_bytes_sent;
  unsigned long remote_bytes_rcvd;
#endif /* CHEAP_PERFORMANCE */

  //  There is one MessageQueue object for each node.  These are used to store
  //  messages that will be combined into one large message
  static MessageQueue *sendQueues;

  //  MessageManagers that are used to store received messages in such
  //  a way that they can be queried very efficiently
  MessageManager *localMsgs;
  MessageManager *remoteMsgs;

  //  This object is used to track the nodes that have messages
  //  waiting for them during the WAIT mode
  IntTree *sendNodes;

  // number of nodes available (numbers 0 ... # nodes - 1)
  int TotalNodes;

  // which node am I?
  int thisNode;
  
  // error status.  Value depends on inherited object.
  CommError ErrorStatus;

  // return whether the given node is available
  int node_available(int n) 
  {
    return(n >= 0 && n < TotalNodes);
  }

  //  Delete message, if appropriate
  void del_msg(Message *, int);

  //
  // protected virtual functions
  //
  
  //  Send all messages that are queued up to go to a specified node
  virtual int do_send_queue(int) = 0;

  //  Send an indivdual message right now
  virtual int do_send_msg(Message *, int, int, int delmsg=TRUE) = 0;

  // check for a message from the given node and tag.  Return a NEW
  // MsgList object if a message arrives, or NULL if no message available.
  // If another message is found while searching for the given message, add
  // it to the received queue, and keep looking.
  // If node = (-1), checks for messages from any node.
  // If tag = (-1), checks for messages with any tag.
  virtual Message *do_receive(int &node, int &tag) = 0;

public:
  // constructor and destructor
  Communicate(void);
  virtual ~Communicate(void);

#ifdef CHEAP_PERFORMANCE
  void print_comm_stats(void);
#endif /* CHEAP_PERFORMANCE */

  // return info about connections in general
  int nodes(void) { return TotalNodes; }
  int this_node(void) { return thisNode; }
  CommError errorno(void) { return ErrorStatus; }

  // enable/disable debugging
  void debug(int newsetting) { Debug = newsetting; }
  int debug(void) { return Debug; }

  // query/set send method
  void send_method(SendMethod sm) 
  {
    sendMethod = sm;

    if (TotalNodes>1)
    {
    	if( sm == NOW )
		send_all();

    	if ( (sm == WAIT) && (sendNodes == NULL) )
	{
		sendNodes = new IntTree();
		if ( sendNodes == NULL )
		{
		  NAMD_die("new failed in Communicate::send_method");
		}
	}
   }
  }

  SendMethod send_method(void) { return sendMethod; }

  //
  //    routines to send/receive data
  //

  // send data to another node.  Returns success (T or F).
  // last argument specifieswhether to delete the Message after sending
  // (if message is for another node).
  int send(Message *, int node, int tag, int delmsg = TRUE);

  // send ALL messages that are currently in the send queues.  
  // Return number of messages sent out.
  int send_all(void);


  // sends the given message, and sends it NOW, regardless of sendMethod
  // last argument specifieswhether to delete the Message after sending
  // (if message is for another node).
  int send_now(Message *msg, int node, int tag, int delmsg = TRUE) 
  {
    int retval;
    SendMethod sm=sendMethod;

    sendMethod = NOW;
    retval = send(msg, node, tag, delmsg);
    sendMethod = sm;

    return retval;
  }
  
  // receive data from another node.  Returns newly created Message object
  // with received message, or NULL if no message is available.
  // If node is < 0, this will receive the next message with the given tag
  // from any node.  If tag < 0, this will receive the next message with
  // any tag from the given node.  If both are < 0, this will receive the
  // next message, period.  node and tag are passed by reference; if either
  // is < 0, and a message is received, they are changed to their actual
  // values.
  Message *receive(int& node, int& tag);

  // broadcast the current message to other nodes.
  // Return number of nodes actually sent to.
  // The first version sends to all nodes including this node.
  // The second version sends to all nodes except this node.
  // The first argument is the Message; the last argument is the tag.
  int broadcast_all(Message *, int);
  int broadcast_others(Message *, int, int delmsg=TRUE);

  //
  // public virtual routines
  //
  
  // add a node to the list of communications.  Return node number if OK,
  //  (-1) if problem.
  virtual int add_node(void *id) = 0;

  //  Get the PVM tids for FMA.  This really shouldn't be a vritual
  //  function, but it has to be in the current setup.  Oh Well.
  virtual int *get_tids() = 0;
};

#endif
