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
 *	$RCSfile: MessageQueue.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:28 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *	These are the non-inlined functions of MessageQueue.  They are not
 * inlined since they have multiple exit points, which is currently not
 * supported on inline functions on the HPs.  This is OK since these
 * functions should not be used very often, and if they are called, they
 * must do a linear search, so the function overhead won't be significant.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: MessageQueue.C,v $
 * Revision 1.1001  1997/03/19 11:54:28  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:58:41  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:50  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:36:24  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/12/06 19:54:12  ari
 * Initial revision
 *
 * Revision 1.3  1995/03/18 02:40:03  nelson
 * Changed so that Message is stored rather than MsgList
 *
 * Revision 1.2  95/03/08  14:36:59  14:36:59  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.1  94/10/24  09:16:45  09:16:45  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

#include "MessageQueue.h"

/************************************************************************/
/*									*/
/*			FUNCTION get_msg_by_node			*/
/*									*/
/*  INPUT:								*/
/*	node - node value to search for					*/
/*									*/
/*	This function searches for the first message on the queue 	*/
/*  with the specifed node value.  If such a message is found, a 	*/
/*  pointer to it is returned and it is removed from the list.		*/
/*	If no message with the given node is found, NULL is returned.	*/
/*									*/
/************************************************************************/

Message *MessageQueue::get_msg_by_node(int &tag, int &node)
{
      MsgList *ptr;		//  Current position in list
      Message *retptr;		//  Message pointer to return
      
      ptr=head;
      
      //  Loop until we find the end of the list or a match
      while ( (ptr != NULL) && (ptr->node != node) )
      {
	 ptr = ptr->next;
      }
      
      if (ptr == NULL)
      {
	 //  Didn't find a match, so return NULL
	 return (NULL);
      }
      
      //  Did find a match, now take it out of the list
      if (ptr->prev != NULL)
      {
	 //  This element wasn't the head, so link the previous node
	 //  to the following node
	 ptr->prev->next = ptr->next;

      }
      else
      {
	 //  This element was the head of the list, so reassign the head
	 head = ptr->next;

	 if (head != NULL)
	 {
	    head->prev=NULL;
	 }
      }

      //  If this wasn't the tail of the list, relink the previous pointer
      //  of the next node to the correct value
      if (ptr->next != NULL)
      {
	 ptr->next->prev = ptr->prev;
      }
      
      //  If this was the tail of the list, then reset the tail ptr
      if (tail == ptr)
      {
	 tail = ptr->prev;
      }
      
      node = ptr->node;
      tag = ptr->tag;
      retptr = ptr->msg;

      delete ptr;
      
      numMsgs--;
      
      return(retptr);
}
/*			END OF FUNCTION get_msg_by_node			*/
   
/************************************************************************/
/*									*/
/*			FUNCTION get_msg_by_tag				*/
/*									*/
/*   INPUTS:								*/
/*	tag - tag to search for						*/
/*									*/
/*	This function is equivalent to get_msg_by_node, except that	*/
/*   it searches for a message with a specific tag.  In the current	*/
/*   implementation where the only class the uses MessageQueue's is     */
/*   the MessageManager class and it maintains a different queue	*/
/*   for each tag, this function is currently never used.  But hey,	*/
/*   it was easy to add, and it will allow this class to be used	*/
/*   more generically later.						*/
/*									*/
/************************************************************************/

Message *MessageQueue::get_msg_by_tag(int &tag, int &node)
{
      MsgList *ptr;	//  Current position in list
      Message *retptr;
      
      ptr=head;
      
      //  Loop until we find the end of the list or a match
      while ( (ptr != NULL) && (ptr->tag != tag) )
      {
	 ptr = ptr->next;
      }
      
      if (ptr == NULL)
      {
	 //  Didn't find a match, so return NULL
	 return (NULL);
      }
      
      //  Did find a match, now take it out of the list
      if (ptr->prev != NULL)
      {
	 //  This element wasn't the head, so link the previous node
	 //  to the following node
	 ptr->prev->next = ptr->next;
      }
      else
      {
	 //  This element was the head of the list, so reassign the head
	 head = ptr->next;
	 if (head != NULL)
	 {
	    head->prev=NULL;
	 }
      }

      //  If this wasn't the tail, relink the previous pointer of the next
      //  node
      if (ptr->next != NULL)
      {
	 ptr->next->prev = ptr->prev;
      }
      
      //  If this was the tail, then reset the tail pointer
      if (tail == ptr)
      {
	 tail = ptr->prev;
      }
      
      node = ptr->node;
      tag = ptr->tag;
      retptr = ptr->msg;

      delete ptr;
      
      numMsgs--;
      
      return(retptr);
}
/*			END OF FUNCTION get_msg_by_tag			*/


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:28 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: MessageQueue.C,v $
 * Revision 1.1001  1997/03/19 11:54:28  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
