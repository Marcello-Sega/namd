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
 *	$RCSfile: MessageManager.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.778 $	$Date: 1997/01/28 00:30:49 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *	MessageManager is a completely inlined class that provides an
 * efficient way to store messages.  The purpose is to be able to return
 * messages with a specific tag in constant time, since that is what is
 * done during the processing loops of namd.
 *	To accomplish this, the class maintains a list of MessageQueues
 * where each queue contains messages with a specifc tag.  The queues are
 * arranged in tag order so that a binary search can be used to find the
 * correct queue.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: MessageManager.h,v $
 * Revision 1.778  1997/01/28 00:30:49  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:23  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:24  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/12/06 19:54:12  ari
 * Initial revision
 *
 * Revision 1.5  1995/10/09 04:51:08  hazen
 * Updated memory allocation to use C++ new/delete
 *
 * Revision 1.4  1995/03/18  13:06:18  nelson
 * Changed scheme for storing to rely on tag values
 *
 * Revision 1.3  95/03/18  02:39:47  02:39:47  nelson (Mark T. Nelson)
 * Changed so that Message is stored rather than MsgList
 * 
 * Revision 1.2  95/03/08  14:35:52  14:35:52  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.1  94/10/24  09:17:31  09:17:31  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef MESSAGEMANAGER_H

#define MESSAGEMANAGER_H

#include "MessageQueue.h"
#include "Communicate.h"
#include "Message.h"
#include "string.h"

#define MINTAGVALUE	100

class MessageManager
{
private:
   int numMessages;		//  Number of messages stored
   int numTags;			//  Number of tag queues
   MessageQueue **queues;	//  Array of MessageQueue pointers
   
public:
   MessageManager();
   
   //  Destructor, remove the queues then the arrays
   ~MessageManager()
   {
      int i;
      
      for (i=0; i<numTags; i++)
      {
	 delete queues[i];
      }
      
      delete [] queues;
   }
   
   int num() {return (numMessages);}	//  Return the number of messages

   //  Add the specified message to the correct queue
   void add_message(Message *new_msg, int tag, int node)
   {
      int queueNum;	//  queue number to add message to

      //  Find the right queue
      queueNum = tag - MINTAGVALUE;

      if (queueNum > numTags)
      {
	  NAMD_die("Maximum tag size exceeded in MessageManager. "
		   "Did someone forget to update MAXTAGVALUE?");
      }

      //  add the message
      queues[queueNum]->add_msg(new_msg, tag, node);
	 
      numMessages++;
   }
   
   //  Find a message given a tag and node value
   Message *find_message(int &tag, int &node)
   {
      int queue_num;		//  queue number to look in
      Message *ret_val;		//  Message to return
      int i;			//  loop counter

      if (tag != -1)
      {
	 //  Specific tag given, find out the queue number for this tag
	 queue_num = tag - MINTAGVALUE;

         if (queue_num > numTags)
         {
	    NAMD_die("Maximum tag size exceeded in MessageManager.  Did someone forget to update MAXTAGVALUE?");
         }
	 else
	 {
	    if (node == -1)
	    {
	       //  No node value given, so just take the first message
	       //  off the queue
	       ret_val = queues[queue_num]->get_head(tag, node);
	    }
	    else
	    {
	       //  A specific node was given, so search in the queue
	       //  for this node
	       ret_val = queues[queue_num]->get_msg_by_node(tag, node);
	    }
	 }
      }
      else
      {
	 //  No tag was given, so we'll have to search all the queues
	 //  until we find a message.  This is not the most efficient means,
	 //  but at the moment, we never do this anyway.
	 if (node == -1)
	 {
	    //  No node value given.  So just go down the queues looking
	    //  for the first message we can find
	    ret_val = NULL;
	    
	    for (i=0; ( (i<numTags) && (ret_val == NULL) ); i++)
	    {
	       ret_val = queues[i]->get_head(tag, node);
	    }
	 }
	 else
	 {
	    //  A specific node was given.  So search each queue for
	    //  a message from this node and return the first one
	    //  we find
	    ret_val = NULL;

	    for (i=0; ( (i<numTags) && (ret_val == NULL) ); i++)
	    {
	       ret_val = queues[i]->get_msg_by_node(tag, node);
	    }
	 }
      }

      //  If we are really returning a message, decrement the message count
      if (ret_val != NULL)
	numMessages--;

      return(ret_val);
   }
};

#endif
