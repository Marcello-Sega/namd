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
 *	$RCSfile: MessageQueue.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:29 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *	This is a mainly inlined class that provides a FIFO queue for storing
 * received messages.  It will add a new message and retrieve the head of
 * the queue in constant time.  It can also search for a message by node
 * or tag in O(number of messages) time.
 *	This class is currently implemented as a doubly linked list.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: MessageQueue.h,v $
 * Revision 1.1001  1997/03/19 11:54:29  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:58:42  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:51  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:24  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:26  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/12/06 19:54:12  ari
 * Initial revision
 *
 * Revision 1.4  1995/10/09 04:53:37  hazen
 * Updated memory allocation to use C++ new/delete
 *
 * Revision 1.3  1995/03/18  02:40:04  nelson
 * Changed so that Message is stored rather than MsgList
 *
 * Revision 1.2  95/03/08  14:38:20  14:38:20  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.1  94/10/24  09:16:45  09:16:45  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef MESSAGEQUEUE_H

#define MESSAGEQUEUE_H

#include "common.h"
#include "Message.h"
#include "Communicate.h"

class MessageQueue
{
typedef struct msglist
{
	Message *msg;
	int node;
	int tag;
	struct msglist *next;
	struct msglist *prev;
} MsgList;

private:
   int numMsgs;		//  Number of messages currently stored
   MsgList *head;	//  Head of the list
   MsgList *tail;	//  Tail of the list
   
public:
   int num() {return numMsgs;}		//  Return  the number of messages
   
   //  This is the constructor, just NULL everything out
   MessageQueue()
   {
      numMsgs=0;
      head=NULL;
      tail=NULL;
   }
   
   //  This is the destructor.  Walk down the list and free each Message
   ~MessageQueue()
   {
      MsgList *ptr, *next;
      
      ptr=head;
      
      while (ptr!=NULL)
      {
	 delete ptr->msg;
	 next = ptr->next;
	 delete ptr;
	 ptr=next;
      }
   }    
   
   //  Add a message to the end of the queue
   void add_msg(Message *new_msg, int tag, int node)
   {
      MsgList *newlink = new MsgList;

      if (newlink == NULL)
      {
	NAMD_die("Memory allocation failed in MessageQueue::add_msg");
      }

      newlink->msg = new_msg;
      newlink->node = node;
      newlink->tag = tag;

      if (head == NULL)
      {
	 //  The list was empty
	 newlink->next=NULL;
	 newlink->prev=NULL;
	 head=newlink;
	 tail=newlink;
      }
      else
      {
	 //  The list wasn't empty, add to the end
	 newlink->next=NULL;
	 newlink->prev=tail;
	 tail->next=newlink;
	 tail=newlink;
      }
      
      numMsgs++;
   }
   
   //  This routine returns a pointer to the message at the head
   //  of the queue and removes this message from the queue.
   Message *get_head(int &tag, int &node)
   {
      MsgList *ptr;	//  Message to return
      Message *retptr;  //  Pointer to returned Message
      
      //  If the list is empty, just return null
      if (head == NULL)
      {
	 retptr=NULL;
      }
      else
      {
         ptr=head;
      
         //  Reassign the head to the next node
         head=head->next;
      
         if (head != NULL)
         {
	    //  Reset the previous pointer of the head node
            head->prev=NULL;
         }
         else
         {
	    //  The list is empty, set the tail pointer to NULL
	    tail = NULL;
         }
      
	 retptr = ptr->msg;
	 node = ptr->node;
	 tag = ptr->tag;

	 delete ptr;

         numMsgs--;
      }
      
      return(retptr);
   }
   
   Message *get_msg_by_tag(int&, int&);	//  Search the queue for a message
					//  with a specific tag
   Message *get_msg_by_node(int&, int&);//  Search the queue for a message
					//  from a specfic node
};

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:29 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: MessageQueue.h,v $
 * Revision 1.1001  1997/03/19 11:54:29  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
