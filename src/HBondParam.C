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
 *	$RCSfile: HBondParam.C,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/10/01 16:46:49 $
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
 * $Log: HBondParam.C,v $
 * Revision 1.1002  1997/10/01 16:46:49  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.1001  1997/03/19 11:54:16  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:58:24  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 01:04:40  ari
 * uplevel
 *
 * Revision 1.2  1997/01/28 00:51:48  ari
 * uplevel
 *
 * Revision 1.1.2.1  1997/01/24 02:29:49  jim
 * Fixed bug where only first parameter file was read!
 * Added files for hydrogen bond parameter reading.
 *
 * Revision 1.2  1996/04/27 23:42:13  billh
 * commented out debugging messages.
 *
 * Revision 1.1  1996/04/18 18:44:29  billh
 * Initial revision
 *
 ***************************************************************************/

#include "HBondParam.h"
#include "Parameters.h"
#include "Inform.h"
#include "MStream.h"


////////////////////////////  constructor
HBondParam::HBondParam(void) {

  // initialize linked list
  pairData = NULL;
  pairDataTail = NULL;
  numParams = 0;
}


////////////////////////////  destructor
HBondParam::~HBondParam(void) {

  // delete all the previously stored pair data, if any
  delete_hbond_data();
}


/***************************************************************************
 *  private methods
 ***************************************************************************/

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
Bool HBondParam::compare_string_wildcard(char *stbeg, char *st, char *wc) {

  // make sure none of the strings are null ...
  if(st == NULL || wc == NULL)
    return FALSE;

  // if we've reached the end of the wildcard string and the end of the
  // comparison string, we have a success
  if(*st == '\0' && *wc == '\0') {
    return TRUE;
  }

/*
  namdWarn << "HBP: Comparing ";
  if(*st == '\0')
    namdWarn << "<nil>";
  else
    namdWarn << *st;
  namdWarn << " to ";
  if(*wc == '\0')
    namdWarn << "<nil>";
  else
    namdWarn << *wc;
  namdWarn << sendmsg;
*/

  // now check for wildcard characters ...
  if (*wc == '+') {             // look for a single digit
    if(isdigit(*st))
      return compare_string_wildcard(st, st+1, wc+1);
    else if(st == stbeg)
      return FALSE;
    else
      return compare_string_wildcard(stbeg, st-1, wc);

  } else if (*wc == '#') {      // look for any number of digits
    if(isdigit(*st))
      return compare_string_wildcard(stbeg, st+1, wc);
    else
      return compare_string_wildcard(stbeg, st, wc+1);

  } else if (*wc == '%') {      // look for a single character
    if(*st != '\0')
      return compare_string_wildcard(st, st+1, wc+1);
    else if(st == stbeg)
      return FALSE;
    else
      return compare_string_wildcard(stbeg, st-1, wc);

  } else if (*wc == '*') {      // look for any number of characters
    if(*st != '\0')
      return compare_string_wildcard(stbeg, st+1, wc);
    else
      return compare_string_wildcard(stbeg, st, wc+1);

  } else {
    // now, since the wc character isn't one of the special ones, look
    // for a simple match against the regular string.  But we must be careful
    // about the end of the regular string.
    if(*st == '\0' && st != stbeg)
      return compare_string_wildcard(stbeg, st-1, wc);
    else if((toupper(*st) == toupper(*wc)))
      return compare_string_wildcard(st, st+1, wc+1);
    else
      return FALSE;
  }
}


// scan the list of stored pairs if one with the given types have
// been already added.  If so, return the pair data structure, otherwise
// return NULL.
// The last argument indicates whether to use wildcard searching.
HBondPairData *HBondParam::search_pair_list(char *d1, char *a1, Bool useWC) {
  Bool found = FALSE;
  HBondPairData *ptr = pairData;

  while(ptr) {

    if(!useWC) {
      found = (! strncasecmp(ptr->name(0), d1, MAX_ATOMTYPE_CHARS) &&
	       ! strncasecmp(ptr->name(1), a1, MAX_ATOMTYPE_CHARS));
    } else {
      found = (compare_string_wildcard(d1, d1, ptr->name(0)) &&
	       compare_string_wildcard(a1, a1, ptr->name(1)));
/*
      namdWarn << "HBP: Wildcard check for " << d1 << " == ";
      namdWarn << ptr->name(0) << " and " << a1 << " == ";
      namdWarn << ptr->name(1) << " ... ";
      namdWarn << found << sendmsg;
*/
    }

    if(found)
      return ptr;

    // not found; look at next item in list
    ptr = ptr->next;
  }

  //  if we're here, nothing was found.
  return NULL;
}


// put an new pair structure at the end of our linked list
void HBondParam::add_to_pair_list(HBondPairData *data) {

  // check to see if the list has been previously started
  if(pairDataTail)
    pairDataTail->next = data;
  else
    pairData = data;

  // put data at end of list
  pairDataTail = data;
  pairDataTail->next = NULL;

  numParams++;
}
  
  
// find the pair data object for  the given type indices from our
// binary tree cache.  Return pointer, or NULL if not found.
HBondPairData *HBondParam::search_binary_tree(Index d1, Index a1) {

  // scan the binary tree, return the object found (if any)
  return (HBondPairData *)pairDataTree.get_data((int)d1, (int)a1);
}


// add data to the binary tree used to cache previously-found pairs
void HBondParam::add_to_binary_tree(Index d1, Index a1, HBondPairData *hbp) {

  // put value in our binary tree.  This should only be called when the
  // data is new, but it will work even if this data already has been added
  // to the tree.
  pairDataTree.add_data((int)d1, (int)a1, (void *)hbp);
}


/***************************************************************************
 *  action methods
 ***************************************************************************/


// clear the cached atom type index binary tree
void HBondParam::clear_indices(void) {
  pairDataTree.clear();
}


// delete all the stored h-bond pair data (also clears the cache)
void HBondParam::delete_hbond_data(void) {
  HBondPairData *ptr = pairData;

  while(ptr) {
    HBondPairData *dptr = ptr;
    ptr = ptr->next;
    delete dptr;
  }

  pairData = pairDataTail = NULL;
  numParams = 0;

  clear_indices();
}


// Stores the provided hbond pair data in the linked list (if it has
// not already been added).  If a duplicate is given, the new
// values are stored.  Returns TRUE if added element was unique, FALSE
// if the pair already had been added.
Bool HBondParam::add_hbond_pair(char *d1, char *a1, Real emin, Real rmin) {
  Bool isUnique = TRUE;

  // scan the list of pairs to see if we already have this one ...
  HBondPairData *data = search_pair_list(d1, a1);

  if (data == NULL) {
    // if the data has not yet been added, do so now
    add_to_pair_list(new HBondPairData(d1, a1, emin, rmin));

  } else {
    // data for these two donor, acceptor names have been previously added.
    // replace the data and indicate we've done this before.
    data->set_emin(emin);
    data->set_rmin(rmin);
    isUnique = FALSE;
  }

  // return whether we had to replace the dataa
  return isUnique;
}


// Search through the list of hydrogen bond pairs and find the data
// for the given atoms.  The user provides the indices of the donor
// and acceptor atom types, and access to the proper Parameters
// object.  The first time through, the strings for the atom types
// are compared to the hbond pair names, and if a match is found,
// the indices are stored in a binary tree for fast retrieval in the
// future.  If no match is found, this returns NULL.  Otherwise,
// a pointer to the HBondPairData structure with the relevant data
// is returned.
HBondPairData *HBondParam::get_hbond_data(Index d1, Index a1, Parameters *p) {
  HBondPairData *data;

  // first,  check to see if we have this item in our cache
  if((data = search_binary_tree(d1, a1)) != NULL)
    return data;

  // we don't have it yet; so try to find a pair of names from our list
  // of donor-acceptor pairs that will match type names.  d1 and a1 are
  // indices of type names as stored in the Parameters object, so we must
  // first get the names from there.
  char *d1name = p->atom_type_name(d1);
  char *a1name = p->atom_type_name(a1);

  // check for a match against these names
  if((data = search_pair_list(d1name, a1name, TRUE)) != NULL) {
    // we found a match; save it and return data
    add_to_binary_tree(d1, a1, data);
    return data;
  }

  // no match at all; return null pointer
  return NULL;
}


// Put necessary info in the given Message, in order to send it to
// another node.  Return number of pair list items put in message, or
// -1 if there is an error.
int HBondParam::create_message(MOStream *msg) {

  if (msg == NULL)
    return (-1);

  // store the number of pair list items
  msg->put(num());

  // put name and parameter data in as separate arrays
  if(num() > 0) {

    // allocate space for data 
    char *names = new char[2 * num() * (MAX_ATOMTYPE_CHARS + 1)];
    Real *rvals = new Real[2 * num()];
    if (names == NULL || rvals == NULL)
      return (-1);

    // put in data
    HBondPairData *ptr = pairData;
    char *nmptr = names;
    Real *rvptr = rvals;
    while(ptr) {
      strcpy(nmptr, ptr->name(0));
      nmptr += (MAX_ATOMTYPE_CHARS + 1);
      strcpy(nmptr, ptr->name(1));
      nmptr += (MAX_ATOMTYPE_CHARS + 1);
      *(rvptr++) = ptr->emin();
      *(rvptr++) = ptr->rmin();
      ptr = ptr->next;
    }

    msg->put(num() * 2 * (MAX_ATOMTYPE_CHARS + 1));
    msg->put(num() * 2 * (MAX_ATOMTYPE_CHARS + 1), names);
    delete[] names;
    msg->put(num() * 2);
    msg->put(num() * 2, rvals);
    delete[] rvals;
  }

  // if here, everything worked
  return num();
}


// retrieve data from the given message and store it; this clears out
// any previously stored data.  Return number of items read from message, or
// -1 if there is an error.
int HBondParam::receive_message(MIStream *msg) {
  int np;

  // get number of items from message
  msg->get(np);

  if (np < 0)
    return (-1);

  // clear out our current values
  delete_hbond_data();

  // get the data itself now
  if(np > 0) {
    
    // first get the pointer to the list of names, and then the pointer
    // to the list of parameters
    int numNames;
    msg->get(numNames);
    char *names = new char[numNames];
    msg->get(numNames, names);
    int numRvals;
    msg->get(numRvals);
    Real *rvals = new Real[numRvals];
    msg->get(numRvals, rvals);

    // now go through the lists and create new parameter items
    int i;
    char *nmptr = names;
    Real *rvptr = rvals;
    for(i=0; i < np; i++) {
      add_hbond_pair(nmptr, nmptr + MAX_ATOMTYPE_CHARS + 1,
		     *rvptr, *(rvptr + 1));

/*
      if(i % 4 == 0) {
	namdWarn << "HBP: Recvd hbond pair " << i << " = " << nmptr << " and ";
	namdWarn << nmptr + MAX_ATOMTYPE_CHARS + 1 << ": emin=" << *rvptr;
	namdWarn << ", rmin=" << *(rvptr + 1) << sendmsg;
      }
*/

      nmptr += (2 * (MAX_ATOMTYPE_CHARS + 1));
      rvptr += 2;
    }
  }

  // done; return number of unique items actually retrieved
  return num();
}


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/10/01 16:46:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: HBondParam.C,v $
 * Revision 1.1002  1997/10/01 16:46:49  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.1001  1997/03/19 11:54:16  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
