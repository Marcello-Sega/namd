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
 *	LintList implements a linked list of non-negative integers.  It is a 
 * completely inlined class.
 *
 ***************************************************************************/

#ifndef LINTLIST_H

#define LINTLIST_H

#include "common.h"
#include <values.h>

//  This is the structure used to actually store the list
typedef struct lint_node
{
   int value;
   struct lint_node *next;
} LintNode;

//  This is the value returned when the list is empty
#define LIST_EMPTY	( -1 * ( MAXINT - 1 ) )

class LintList
{
private:
   LintNode *llhead;		//  head of the list
   LintNode *tail;		//  tail of the list
   LintNode *cur;		//  cur value in traversal of list
   
public:
   LintList() {llhead=NULL; tail=NULL; cur=NULL;}
   
   ~LintList()
   {
      LintNode *ptr=llhead;	//  Current node
      LintNode *next;		//  Next node
      
      //  Walk down and destroy each node one by one
      while (ptr != NULL)
      {
	 next=ptr->next;
	 delete ptr;
	 ptr=next;
      }
   }
   
   //  This function returns the first integer in the list and
   //  resets the current position in the list to the head.  If
   //  the list is empty, then LIST_EMPTY is returned
   int head()
   {
      if (llhead == NULL)
      	return(LIST_EMPTY);
      else
      {
	 cur=llhead;
	 return(cur->value);
      }
   }
   
   //  This function returns the next integer in the list from
   //  the current traversal position.  If the list is empty or
   //  the traversal has reached the end of the list, then 
   //  LIST_EMPTY is returned
   int next()
   {
      if ( (cur == NULL) || (cur->next == NULL) )
      {
	 return(LIST_EMPTY);
      }
      else
      {
	 cur = cur->next;
	 return(cur->value);
      }
   }
   
   //  Add a new value to the tail of the list
   void add(int new_value)
   {
      LintNode *new_node;	//  New node to add to list
      
      new_node = new LintNode;
      
      if (new_node == NULL)
      {
	 NAMD_die("memory allocation failed in LintList::add");
      }
      
      new_node->value = new_value;
      new_node->next = NULL;
      
      if (tail == NULL)
      {
	 //  The list was empty
	 llhead = new_node;
	 tail = new_node;
      }
      else
      {
	 //  add it to the tail
	 tail->next = new_node;
	 tail = new_node;
      }
   }

   //  Return the number of elements in the list in O(n) time
   int num()
   {
	int nele=0;
	LintNode *ptr = llhead;

	while (ptr != NULL)
	{
		nele++;
		ptr=ptr->next;
	}

	return(nele);
   }

   //  Return the number of occurrences of a in the list in O(n) time
   int num(int match)
   {
	int nele=0;
	LintNode *ptr = llhead;

	while (ptr != NULL)
	{
		if ( ptr->value == match ) nele++;
		ptr=ptr->next;
	}

	return(nele);
   }

   //  Delete the first occurrence of a in the list in O(n) time
   void del(int match)
   {
	LintNode **ptr = &llhead;

	for(; (*ptr) != NULL && (*ptr)->value != match; ptr=&((*ptr)->next));

	if ( (*ptr) != NULL )
	{
		LintNode *ptr2 = (*ptr)->next;
		int newtail = ( *ptr == tail );
		delete *ptr;
		(*ptr) = ptr2;
		if ( newtail )
		{
			LintNode *ptr3 = llhead;
			if ( ptr3 ) while ( ptr3->next ) ptr3 = ptr3->next;
			tail = ptr3;
		}
	}

   }

   //  Ditto for all elements of another LintList
   void del(LintList &ll)
   {
	int match;
	match = ll.head();
	while ( match != LIST_EMPTY )
	{
		del(match);
		match = ll.next();
	}
   }

   //  Add a new value only if it doesn't already exist in O(n) time
   void merge(int match)
   {
	if ( ! num(match) ) add(match);
   }

   //  Ditto for all elements of another LintList
   void merge(LintList &ll)
   {
	int match;
	match = ll.head();
	while ( match != LIST_EMPTY )
	{
		merge(match);
		match = ll.next();
	}
   }

};

#endif   

