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
 *	$RCSfile: IntTree.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:36:17 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *	IntTree maintains a binary tree of integers.  It is currently only
 * used for storing the bonded coordinates to be returned by a patch.  It
 * only stores DISTINCT values.  If a duplicate value is passsed to add_value,
 * the tree remains unchanged.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: IntTree.h,v $
 * Revision 1.777  1997/01/17 19:36:17  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/12/06 19:52:20  ari
 * Initial revision
 *
 * Revision 1.2  1995/03/08 14:46:34  nelson
 * Added copyright
 *
 * Revision 1.1  94/09/04  20:24:46  20:24:46  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef INTTREE_H

#define INTTREE_H

#include "common.h"

typedef struct int_node
{
	int value;
	struct int_node *left;
	struct int_node *right;
} IntNode;


class IntTree
{
private:

	int numNodes;			//  Number of nodes in the tree
	IntNode *tree;			//  Tree itself

	IntNode *insert_value(int, IntNode *);
					//  Recursive function to add a
					//  value to the tree
	int populate_array(IntNode *, int *, int);
					//  Traverse the tree in order
					//  placing each element into an array
	void destroy_tree(IntNode *); 	//  Recursive function used by
					//  the destructor to tear down the
					//  tree

public:
	IntTree();			//  Constructor
	~IntTree();			//  Destructor
	int size() {return numNodes;}   //  Return the size of the tree

	//  Insert a value into the tree
	void add_value(int intval)
	{ tree = insert_value(intval, tree); }

	int *make_array();		//  Convert the tree into an array
};

#endif
