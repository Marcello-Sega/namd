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
 *	VoidTree maintains a binary tree, with two integers as keys. It
 * only stores DISTINCT values.  If a duplicate value is passsed to
 * add_value, the tree remains unchanged.  For each value, there is an
 * associated data item, stored as a void pointer.
 *
 ***************************************************************************/
#ifndef VOIDTREE_H
#define VOIDTREE_H

#include "common.h"


// structure used for each element in the tree
typedef struct void_node {
	int value1, value2;
	void *item;
	struct void_node *left;
	struct void_node *right;
} VoidNode;


////////////////////////// VoidTree class definition
class VoidTree {

private:
  int numNodes;				//  Number of nodes in the tree
  VoidNode *tree;			//  Tree itself

  // Recursive function to add a value to the tree
  VoidNode *insert_value(int, int, void *, VoidNode *);

  // Recursive function to find a value in the tree, and return
  // its item data
  void *locate_value(int, int, VoidNode *);

  //  Traverse the tree in order, placing each data item into an array
  int populate_data_array(VoidNode *, void **, int);

  // Recursive function used by the destructor to tear down the tree
  void destroy_tree(VoidNode *);

public:
  // constructor/destructor
  VoidTree(void);
  ~VoidTree(void);

  // return the size of the tree
  int size(void) { return numNodes; }

  // clear the currently stored tree
  void clear(void) { destroy_tree(tree); tree = NULL; numNodes = 0; }

  // Insert a value into the tree, only giving one value (the other will be 0)
  void add_data(int intval, void *item = NULL) {
    tree = insert_value(intval, 0, item, tree);
  }

  // Insert a value into the tree
  void add_data(int intval1, int intval2, void *item = NULL) {
    tree = insert_value(intval1, intval2, item, tree);
  }

  // return the item stored for the given integer value, or NULL
  // if not found or not stored
  void *get_data(int intval) { return locate_value(intval, 0, tree); }
  void *get_data(int iv1, int iv2) { return locate_value(iv1, iv2, tree); }

  //  Convert the tree into an array of the data items.  Use the given
  //  array if possible, or create a new one if necessary.
  void **make_data_array(void ** = (void **)NULL);

};

#endif

