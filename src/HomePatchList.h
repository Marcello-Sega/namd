//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef HOMEPATCHLIST_H
#define HOMEPATCHLIST_H
#include "NamdTypes.h"

class HomePatch;

class HomePatchElem {
public:
  PatchID   pid;
  HomePatch *p;

  operator<(HomePatchElem e) { return (pid < e.pid); }
  operator==(HomePatchElem e) { return (pid == e.pid); }

  HomePatchElem(PatchID id=-1, HomePatch *patch=NULL) : pid(id), p(patch) {};
  ~HomePatchElem() { };
  HomePatchElem& operator=(const HomePatchElem &e) { 
    pid = e.pid; p = e.p;  // Do not delete p!  This op only used to shuffle
			   // we delete the p here only when the HomePatch is 
			   // moved off!
    return(*this);
  };
};

typedef SortedArray<HomePatchElem> HomePatchList;

#endif /* HOMEPATCHLIST_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: HomePatchList.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1000 $	$Date: 1997/02/06 15:58:27 $
 *
 ***************************************************************************
 * $Log: HomePatchList.h,v $
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
