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
 *	$RCSfile: IntTree.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/12/06 19:52:20 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: IntTree.C,v $
 * Revision 1.1  1996/12/06 19:52:20  ari
 * Initial revision
 *
 * Revision 1.3  1995/10/09 03:45:52  hazen
 * Updated memory allocation to use C++ new/delete
 *
 * Revision 1.2  1995/03/08  14:46:27  nelson
 * Added copyright
 *
 * Revision 1.1  94/09/04  20:24:46  20:24:46  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/IntTree.C,v 1.1 1996/12/06 19:52:20 ari Exp $";

#include "IntTree.h"

/************************************************************************/
/*									*/
/*			FUNCTION IntTree				*/
/*									*/
/*	This is the constructor for the IntTree class			*/
/*									*/
/************************************************************************/

IntTree::IntTree()
{
	numNodes = 0;
	tree = NULL;
}
/*			END OF FUNCTION IntTree				*/

/************************************************************************/
/*									*/
/*			FUNCTION ~IntTree				*/
/*									*/
/*	This is the destructor for the IntTree class			*/
/*									*/
/************************************************************************/

IntTree::~IntTree()
{
	if (tree != NULL)
	{
		destroy_tree(tree);
	}
}
/*			END OF FUNCTION ~IntTree			*/

/************************************************************************/
/*									*/
/*			FUNCTION destroy_tree				*/
/*									*/
/*   INPUTS:								*/
/*	ptr - Pointer to the head of the tree to destroy		*/
/*									*/
/*	This is a recusive function to free the binary tree.  		*/
/*									*/
/************************************************************************/

void IntTree::destroy_tree(IntNode *ptr)

{
	if (ptr->left != NULL)
	{
		destroy_tree(ptr->left);
	}

	if (ptr->right != NULL)
	{
		destroy_tree(ptr->right);
	}

	CmiFree(ptr);
}
/*			END OF FUNCTION destroy_tree			*/

/************************************************************************/
/*									*/
/*			FUNCTION insert_value				*/
/*									*/
/*   INPUTS:								*/
/*	intval - value to add						*/
/*	itree - tree to add value to					*/
/*									*/
/*	This is a recursive function that adds a value to the given     */
/*   tree.  As mentioned before, only distinct values are stored, so    */
/*   if intval already exists in the tree, the tree will remain 	*/
/*   unchanged.								*/
/*									*/
/************************************************************************/

IntNode *IntTree::insert_value(int intval, IntNode *itree)

{
	IntNode *newNode;		//  New node

	if (itree == NULL)
	{
		//  This tree is currently empty, so make a new node
		newNode = new IntNode;

		if (newNode == NULL)
		{
			NAMD_die("memory allocation failed in IntTree::insert_value");
		}

		newNode->value = intval;
		newNode->left = NULL;
		newNode->right = NULL;

		numNodes++;

		return(newNode);
	}

	if (itree->value == intval)
	{
		//  Duplicate value, do nothing
		return(itree);
	}

	if (itree->value > intval)
	{
		//  Add to left child
		itree->left = insert_value(intval, itree->left);
	}
	else
	{
		//  Add to right child
		itree->right = insert_value(intval, itree->right);
	}

	return(itree);
}
/*			END OF FUNCTION insert_value			*/

/************************************************************************/
/*									*/
/*			FUNCTION make_array				*/
/*									*/
/*	This function converts all the elements stored in the tree into */
/*  and array of integers.  If the tree is currently empty, NULL is     */
/*  returned.								*/
/*									*/
/************************************************************************/

int *IntTree::make_array()

{
	int *new_array;		//  Array to be created

	if (numNodes == 0)
	{
		return(NULL);
	}

	//  Allocate the array
	new_array = new int[numNodes];

	if (new_array == NULL)
	{
		NAMD_die("memory allocation failed in IntTree::make_array");
	}

	//  Populate the array
	populate_array(tree, new_array, 0);

	return(new_array);
}
/*			END OF FUNCTION make_array			*/

/************************************************************************/
/*									*/
/*			FUNCTION populate_array				*/
/*									*/
/*   INPUTS:								*/
/*	itree - Tree to populate array from				*/
/*	iarray - Array of integers to populate				*/
/*	cur_index - current index into array				*/
/*									*/
/*	This is a recursive function to populate the given array of     */
/*   integers with the values from the tree.				*/
/*									*/
/************************************************************************/

int IntTree::populate_array(IntNode *itree, int *iarray, int cur_index)

{
	//  Get the left subtree, if there is one
	if (itree->left != NULL)
	{
		cur_index = populate_array(itree->left, iarray, cur_index);
	}

	//  Assign the current node
	iarray[cur_index] = itree->value;
	cur_index++;

	//  Get the right subtree, if there is one
	if (itree->right != NULL)
	{
		cur_index = populate_array(itree->right, iarray, cur_index);
	}

	return(cur_index);
}
/*			END OF FUNCTION populate_array			*/
