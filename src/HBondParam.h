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
 *	$RCSfile: HBondParam.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:17 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *	
 * HBondParam stores parameters needed for hydrogen bond calculations, and
 * provides methods for accessing this data.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HBondParam.h,v $
 * Revision 1.1001  1997/03/19 11:54:17  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:58:25  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:34  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.2  1997/01/27 22:45:13  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777.2.1  1997/01/24 02:29:50  jim
 * Fixed bug where only first parameter file was read!
 * Added files for hydrogen bond parameter reading.
 *
 * Revision 1.1  1996/04/18 18:44:29  billh
 * Initial revision
 *
 ***************************************************************************/

#ifndef HBONDPARAM_H
#define HBONDPARAM_H

#include "HBondPairData.h"
#include "VoidTree.h"
#include "structures.h"
#include "common.h"

// forward references
class HBondParam;
class Parameters;
class Message;


////////////////////////// HBondParam class definition
class HBondParam {

public:
  //
  // constructor and destructor.
  // ctor args: none
  //
  HBondParam(void);
  ~HBondParam(void);

  //
  // query methods
  //

  // return number of hbond pair data items we have
  int num(void) { return numParams; }

  // given a pair data structure, extract the stored values
  Real Emin(HBondPairData *hbpd) { return hbpd->emin(); }
  Real Rmin(HBondPairData *hbpd) { return hbpd->rmin(); }
  Real A(HBondPairData *hbpd) { return hbpd->A(); }
  Real B(HBondPairData *hbpd) { return hbpd->B(); }

  //
  // action methods
  //

  // clear the cached atom type index binary tree
  void clear_indices(void);

  // delete all the stored h-bond pair data (also clears the cache)
  void delete_hbond_data(void);

  // Stores the provided hbond pair data in the linked list (if it has
  // not already been added).  If a duplicate is given, the new
  // values are stored.  Returns TRUE if added element was unique, FALSE
  // if the pair already had been added.
  Bool add_hbond_pair(char *, char *, Real, Real);

  // Search through the list of hydrogen bond pairs and find the data
  // for the given atoms.  The user provides the indices of the donor
  // and acceptor atom types, and access to the proper Parameters
  // object.  The first time through, the strings for the atom types
  // are compared to the hbond pair names, and if a match is found,
  // the indices are stored in a binary tree for fast retrieval in the
  // future.  If no match is found, this returns NULL.  Otherwise,
  // a pointer to the HBondPairData structure with the relevant data
  // is returned.
  HBondPairData *get_hbond_data(Index, Index, Parameters *);

  // Put necessary info in the given Message, in order to send it to
  // another node.  Return number of pair list items put in message, or
  // -1 if there is an error.
  int create_message(Message *);

  // retrieve data from the given message and store it; this clears out
  // any previously stored data.  Return number of items read from message, or
  // -1 if there is an error.
  int receive_message(Message *);

private:
  // linked list storing the pair parameter data
  HBondPairData *pairData, *pairDataTail;

  // number of pair data items
  int numParams;

  // tree with hbond pair information
  VoidTree pairDataTree;

  //
  // private methods
  //

  // recursively check the given normal string against the wildcard
  // string.  The wildcard matching follows these rules:
  //      * - match any string (can be of zero length)
  //      % - match a single char
  //      # - match any set of digits (0 ... 9) (can be of zero length)
  //      + - match a single digit (0 ... 9)
  //      all other characters - match exactly (case insensitive)
  // st is the string to check, and wc is the wildcard string; stbeg is
  // a pointer to the start of the string being checked (st and wc change
  // as we recursively call this routine).
  // Returns TRUE if match found, FALSE otherwise.
  Bool compare_string_wildcard(char *stbeg, char *st, char *wc);

  // scan the list of stored pairs if one with the given types have
  // been already added.  If so, return the pair data structure, otherwise
  // return NULL.
  // The last argument indicates whether to use wildcard searching.
  HBondPairData *search_pair_list(char *, char *, Bool = FALSE);

  // put an new pair structure at the end of our linked list
  void add_to_pair_list(HBondPairData *);

  // find the pair data object for  the given type indices from our
  // binary tree cache.  Return pointer, or NULL if not found.
  HBondPairData *search_binary_tree(Index, Index);

  // add data to the binary tree used to cache previously-found pairs
  void add_to_binary_tree(Index, Index, HBondPairData *);
};


#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:17 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HBondParam.h,v $
 * Revision 1.1001  1997/03/19 11:54:17  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
