/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   HBondParam stores parameters needed for hydrogen bond calculations, and
   provides methods for accessing this data.
*/

#ifndef HBONDPARAM_H
#define HBONDPARAM_H

#include "HBondPairData.h"
#include "VoidTree.h"
#include "structures.h"
#include "common.h"

// forward references
class HBondParam;
class Parameters;
class MIStream;
class MOStream;


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
  int create_message(MOStream *);

  // retrieve data from the given message and store it; this clears out
  // any previously stored data.  Return number of items read from message, or
  // -1 if there is an error.
  int receive_message(MIStream *);

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

