
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *	VoidTree maintains a binary tree, with two integers as keys. It
 * only stores DISTINCT values.  If a duplicate value is passsed to
 * add_value, the tree remains unchanged.  For each value, there is an
 * associated data item, stored as a void pointer.
 *
 ***************************************************************************/

#include "VoidTree.h"

typedef void *voidptr;


/************************************************************************/
/*									*/
/*			FUNCTION VoidTree				*/
/*									*/
/*	This is the constructor for the VoidTree class			*/
/*									*/
/************************************************************************/

VoidTree::VoidTree(void)
{
	numNodes = 0;
	tree = NULL;
}
/*			END OF FUNCTION VoidTree			*/

/************************************************************************/
/*									*/
/*			FUNCTION ~VoidTree				*/
/*									*/
/*	This is the destructor for the VoidTree class			*/
/*									*/
/************************************************************************/

VoidTree::~VoidTree()
{
	clear();
}
/*			END OF FUNCTION ~VoidTree			*/

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

void VoidTree::destroy_tree(VoidNode *ptr)

{
        if(!ptr)
          return;

	destroy_tree(ptr->left);
	destroy_tree(ptr->right);

	delete ptr;
}
/*			END OF FUNCTION destroy_tree			*/

/************************************************************************/
/*									*/
/*			FUNCTION insert_value				*/
/*									*/
/*   INPUTS:								*/
/*	intval1, intval2 - indices of value to add			*/
/*      vitem - data to store (may be NULL)                             */
/*	itree - tree to add value to					*/
/*									*/
/*	This is a recursive function that adds a value to the given     */
/*   tree.  As mentioned before, only distinct values are stored, so    */
/*   if intval already exists in the tree, the tree will remain 	*/
/*   unchanged.								*/
/*									*/
/************************************************************************/

VoidNode *VoidTree::insert_value(int intval1, int intval2,
				 void *vitem, VoidNode *itree)

{
	VoidNode *newNode;		//  New node

	if (itree == NULL)
	{
		//  This tree is currently empty, so make a new node
		newNode = new VoidNode;
		newNode->value1 = intval1;
		newNode->value2 = intval2;
		newNode->item = vitem;
		newNode->left = NULL;
		newNode->right = NULL;

		numNodes++;

		return(newNode);
	}

	if (itree->value1 == intval1 && itree->value2 == intval2)
	{
		//  Duplicate value, only change is if a new item is given
		if(!itree->item)
			itree->item = vitem;

		// but other than that, no change
		return(itree);
	}

	if (itree->value1 > intval1 || 
	    (itree->value1 == intval1 && itree->value2 > intval2))
	{
		//  Add to left child
		itree->left=insert_value(intval1,intval2,vitem,itree->left);
	}
	else
	{
		//  Add to right child
		itree->right=insert_value(intval1,intval2,vitem,itree->right);
	}

	return(itree);
}
/*			END OF FUNCTION insert_value			*/

/************************************************************************/
/*									*/
/*			FUNCTION locate_value				*/
/*									*/
/*   INPUTS:								*/
/*	intval1, intval2 - values to find				*/
/*	itree - tree to find value in					*/
/*									*/
/*   OUTPUT:								*/
/*	The item stored with the given value, or NULL if not found.	*/
/*									*/
/*	This is a recursive function that finds a value in the given    */
/*   tree.  As mentioned before, only distinct values are stored, so    */
/*   if intval already exists in the tree, the tree will remain 	*/
/*   unchanged.								*/
/*									*/
/************************************************************************/

void *VoidTree::locate_value(int intval1, int intval2, VoidNode *itree)
{
	if(!itree)
		return NULL;

	if(itree->value1 == intval1 && itree->value2 == intval2)
		return itree->item;
	else if(itree->value1 > intval1 ||
		(itree->value1 == intval1 && itree->value2 > intval2))
		return locate_value(intval1, intval2, itree->left);
	else
		return locate_value(intval1, intval2, itree->right);
}
/*			END OF FUNCTION locate_value			*/


/************************************************************************/
/*                                                                      */
/*                      FUNCTION make_data_array                        */
/*                                                                      */
/*      This function converts all the elements stored in the tree into */
/*  and array of void *'s.  If the tree is currently empty, NULL is     */
/*  returned.                                                           */
/*                                                                      */
/************************************************************************/
 
void **VoidTree::make_data_array(void **vlist)
 
{
        void **new_array;         //  Array to be created
 
        if (numNodes == 0)
        {
                return(vlist);
        }
 
        //  Allocate the array
        new_array = (vlist ? vlist : new voidptr[numNodes]);
 
        //  Populate the array
        populate_data_array(tree, new_array, 0);
 
        return(new_array);
}
/*                      END OF FUNCTION make_data_array                 */


/************************************************************************/
/*									*/
/*			FUNCTION populate_data_array			*/
/*									*/
/*   INPUTS:								*/
/*	itree - Tree to populate array from				*/
/*	iarray - Array of void *'s to populate				*/
/*	cur_index - current index into array				*/
/*									*/
/*	This is a recursive function to populate the given array of     */
/*   void *'s with the values from the tree.				*/
/*									*/
/************************************************************************/

int VoidTree::populate_data_array(VoidNode *itree, void **iarray, int cur_index)

{
	//  Get the left subtree, if there is one
	if (itree->left != NULL)
	{
		cur_index=populate_data_array(itree->left,iarray,cur_index);
	}

	//  Assign the current node
	iarray[cur_index] = itree->item;
	cur_index++;

	//  Get the right subtree, if there is one
	if (itree->right != NULL)
	{
		cur_index=populate_data_array(itree->right,iarray,cur_index);
	}

	return(cur_index);
}
/*			END OF FUNCTION populate_data_array		*/

