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
 *	$RCSfile: Communicate.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:26 $
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
 * $Log: Communicate.C,v $
 * Revision 1.777  1997/01/17 19:35:26  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/12/06 19:52:20  ari
 * Initial revision
 *
 * Revision 1.10  1996/01/28 21:50:08  jean
 * Attempting to make stable RCS without Mark Nelson's
 * fma/pairlist decoupling
 *
 * Revision 1.11  1995/12/04 20:45:46  brunner
 * More communication stats - message sizes
 *
 * Revision 1.10  1995/11/22 12:10:05  brunner
 * Added print_comm_stats(), and message traffic accounting.
 *
 * Revision 1.9  95/10/06  17:54:12  17:54:12  hazen (Brett Hazen)
 * Memory Allocation error-checking added
 * 
 * Revision 1.8  1995/03/18  13:04:49  nelson
 * Fixed memory leak in do_send_queue
 *
 * Revision 1.7  95/03/18  02:42:09  02:42:09  nelson (Mark T. Nelson)
 * Reworked extensively to improve performance
 * 
 * Revision 1.6  95/03/08  14:36:32  14:36:32  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.5  94/10/24  09:25:27  09:25:27  nelson (Mark T. Nelson)
 * Added MessageManagers to make looking for received messages much
 * more efficient
 * 
 * Revision 1.4  94/07/26  16:51:34  16:51:34  billh (Bill Humphrey)
 * (Hopefully) fixed problem with broadcast ... now just loops through nodes
 * and sends copy of message to each node, instead of having specific
 * broadcast routine in child class.
 * 
 * Revision 1.3  94/07/07  06:13:14  06:13:14  dalke (Andrew Dalke)
 * Added instances for the static members "sendList" and "recList"
 * 
 * Revision 1.2  1994/07/03  01:22:24  billh
 * Made Communicate enum's and MsgItem struct public members of
 * Communicate, instead of global.
 *
 * Revision 1.1  94/07/03  00:18:51  00:18:51  billh (Bill Humphrey)
 * Initial revision
 * 
 ***************************************************************************/
static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Communicate.C,v 1.777 1997/01/17 19:35:26 ari Exp $";

#include "Communicate.h"
#include "Message.h"
#include "common.h"
#include "MessageManager.h"

#ifdef CHEAP_PERFORMANCE
#include "Inform.h"
#endif /* CHEAP_PERFORMANCE */

// static members
MessageQueue *Communicate::sendQueues = NULL;

// constructor
Communicate::Communicate(void) 
{

#ifdef CHEAP_PERFORMANCE
  num_local_sends=0;
  num_local_recvs=0;

  num_remote_sends=0;
  num_remote_recvs=0;
  remote_bytes_sent=0;
  remote_bytes_rcvd=0;
#endif

  sendQueues = NULL;
  
  TotalNodes = 0;
  thisNode = 0;
  ErrorStatus = NOERROR;
  Debug = 0;
  sendMethod = NOW;

  sendNodes = NULL;

  //  Allocate MessageManagers for local and remote messages
  localMsgs  = new MessageManager();
  remoteMsgs = new MessageManager();
  if ( localMsgs == NULL || remoteMsgs == NULL )
  {
    NAMD_die("new failed in Communicate::Communicate");
  }
}

// destructor
Communicate::~Communicate(void) 
{
    //  Delete the array of send queues
    if (sendQueues != NULL)
    {
	delete [] sendQueues;
    }

    //  Delete the MessageManager objects
    delete localMsgs;
    delete remoteMsgs;

    if (sendNodes != NULL)
	delete sendNodes;
}

#ifdef CHEAP_PERFORMANCE
void Communicate::print_comm_stats(void)
{
  namdInfo << "COMM: " << num_local_sends << " LOCAL SENDS\n";
  namdInfo << "COMM: " << num_local_recvs << " LOCAL RECEIVES\n";

  namdInfo << "COMM: " << num_remote_sends << " REMOTE SENDS\n";
  namdInfo << "COMM: " << num_remote_recvs << " REMOTE RECEIVES\n";
  namdInfo << "COMM: " << remote_bytes_sent << " REMOTE BYTES SENT\n";
  namdInfo << "COMM: " << remote_bytes_rcvd << " REMOTE BYTES RECEIVED\n";
  namdInfo << "COMM: AVG SEND SIZE " 
	   << remote_bytes_sent/(float)num_remote_sends << " BYTES\n";
  namdInfo << "COMM: AVG RECEIVE SIZE " 
	   << remote_bytes_rcvd/(float)num_remote_recvs << " BYTES\n"
	   << sendmsg;
}
#endif /* CHEAP_PERFORMANCE */

// send data to another node.  Returns success (T or F).
// last argument specfies whether to delete the Message object after the
// msg is sent to another node.
int Communicate::send(Message *msg, int node, int tag, int delmsg) 
{
  int retval = TRUE;

  if(node < 0 || node >= TotalNodes || tag < 0)
    return FALSE;

  // If this message is addressed to this node, put in receive queue, otherwise
  // put it in the send queue
  if(node == thisNode) 
  {
	localMsgs->add_message(msg, tag, node);
#ifdef CHEAP_PERFORMANCE
	num_local_sends++;
#endif
  } 
  else 
  {
    //  Check to see if we should send it now or later
    if(sendMethod == NOW) 
    {
      retval = do_send_msg(msg, node, tag, delmsg);

    } 
    else 
    {
      sendQueues[node].add_msg(msg, tag, node);
      sendNodes->add_value(node);
    }
  }
  
  return TRUE;
}


// receive data from another node.  Returns newly created Message object
// with received message, or NULL if no message is available.
// If node is < 0, this will receive the next message with the given tag
// from any node.  If tag < 0, this will receive the next message with
// any tag from the given node.  If both are < 0, this will receive the
// next message, period.  node and tag are passed by reference; if either
// is < 0, and a message is received, they are changed to their actual
// values.
// Messages are searched for in this order (if node == -1) :
//	1. In receive queue, from another node
//	2. Pending in network
//	3. In receive queue, from this node
Message* Communicate::receive(int& node, int& tag) 
{
  Message *msg = NULL;

  if ((TotalNodes > 1) && ((node == -1) || ((node!=-1) && (node!=thisNode))))
  {

	//  we have multiple processors and are looking either for any
	//  node or a remote node.  So first, check the remote messages
	//  that have already been received
	msg = remoteMsgs->find_message(tag, node);

	//  If we didn't find a match, check for new messages from the
	//  network
	if (msg == NULL)
	{
      	    msg = do_receive(node, tag);
	}
  }

  //  If we haven't already found a match, and we're not looking for
  //  something specifically from another node, then check the local
  //  messages
  if ( (msg==NULL) && !( (node!=-1) && (node != thisNode) ) )
  {
	msg = localMsgs->find_message(tag, node);
#ifdef CHEAP_PERFORMANCE
	if (msg != NULL )
	  num_local_recvs++;
#endif /* CHEAP_PERFORMANCE */
  }

  if (msg != NULL)
	msg->reset();

  return(msg);
}


//  Send all messages.  Go through all of the sendQueues and empty out
//  any messages that might be in the queues
int Communicate::send_all(void) 
{
  int sent = 0;		//  Total number of messages sent
  int nodenum;		//  Current node being checked
  int i;		//  Loop counter
  int *nodeList;	//  List of nodes to send to

  if ( (sendNodes == NULL) || (sendNodes->size() == 0) )
  {
	return(0);
  }

  nodeList = sendNodes->make_array();

  for (i=0; i<sendNodes->size(); i++)
  {
	nodenum = nodeList[i];

	sent += do_send_queue(nodenum);
  }

  delete [] nodeList;
  delete sendNodes;
  sendNodes = NULL;

  if (sendMethod == WAIT)
  {
    sendNodes = new IntTree;
    if ( sendNodes == NULL )
    {
      NAMD_die("new failed in Communicate::Communicate");
    }
  }

  // report total number of messages sent
  return sent;
}


// broadcast the given message to ALL nodes, including this node.
// return number of nodes sent to.
// arguments are the Message, and the tag for the message.
int Communicate::broadcast_all(Message *msg, int tag) 
{

  // send message to all other nodes
  for(int i=0; i < TotalNodes; i++) 
  {
    if(i != thisNode) {
      do_send_msg(msg, i, tag, FALSE);
    }
  }
    
  // send message to this node; since we do this, don't need to delete msg
  localMsgs->add_message(msg, tag, thisNode);
#ifdef CHEAP_PERFORMANCE
      num_local_sends++;
#endif /* CHEAP_PERFORMANCE */
 
  return TotalNodes;
}


// broadcast the given message to all OTHER nodes, but not this node.
// return number of nodes sent to.
// arguments are the Message, and the tag for the message.
int Communicate::broadcast_others(Message *msg, int tag, int delmsg) 
{
  // send message to all other nodes
  for(int i=0; i < TotalNodes; i++) 
  {
    if(i != thisNode) {
      do_send_msg(msg, i, tag, FALSE);
    }
  }
    
  // delete message
  if (delmsg)
  	delete msg;
  
  return TotalNodes - 1;
}
