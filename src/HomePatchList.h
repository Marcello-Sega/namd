//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: Helper class for PatchMgr to manager HomePatch(es)
 *              It is object contained in the container HomePatchList
 *
 ***************************************************************************/

#ifndef HOMEPATCHLIST_H
#define HOMEPATCHLIST_H
#include "NamdTypes.h"

class HomePatch;

class HomePatchElem {
public:
  PatchID   pid;
  HomePatch *patch;

  int operator<(HomePatchElem e) { return (pid < e.pid); }
  int operator==(HomePatchElem e) { return (pid == e.pid); }

  HomePatchElem(PatchID id=-1, HomePatch *p=NULL) : pid(id), patch(p) {};
  ~HomePatchElem() { };
  HomePatchElem& operator=(const HomePatchElem &e) { 
    pid = e.pid; 
    patch = e.patch;  // Do not delete patch!  This op only used to shuffle
		      // we delete the patch here only when the HomePatch is 
		      // moved off!
    return(*this);
  };
};

typedef SortedArray<HomePatchElem> HomePatchList;
typedef ResizeArrayIter<HomePatchElem> HomePatchListIter;

#endif /* HOMEPATCHLIST_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: HomePatchList.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1998/01/15 05:40:52 $
 *
 ***************************************************************************
 * $Log: HomePatchList.h,v $
 * Revision 1.1003  1998/01/15 05:40:52  jim
 * Added int return type to comparison operators.
 *
 * Revision 1.1002  1997/03/06 22:06:03  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 * Revision 1.1001  1997/02/26 16:53:10  ari
 * Cleaning and debuging for memory leaks.
 * Adding comments.
 * Removed some dead code due to use of Quiescense detection.
 *
 * Revision 1.1000  1997/02/06 15:58:27  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.2  1997/02/06 15:53:13  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.1.2.1  1997/02/05 22:20:02  ari
 * Additional File for migration -
 * this is along the lines of a clean-up - HomePatchList
 * is shared among different objects, and this deserves
 * a file of its own.
 *
 *
 ***************************************************************************/
