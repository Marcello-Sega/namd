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
 *	$RCSfile: IntList.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/06 20:38:38 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Generic list of integers implemented as an extendable array.  This
 * class is completely inlined.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: IntList.h,v $
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.14  1996/04/17 23:40:14  billh
 * added reset routine.
 *
 * Revision 1.13  1996/03/29 23:00:29  jean
 * fixed excise to shrink buffer if possible
 * if not all the atoms in the input list are removed, then
 * there is a possibility that the buffer will shrink too much.
 * this is caught and causes namd to exit with an error message.
 *
 * Revision 1.12  1996/03/21 23:01:50  jean
 * modified excise to take a list of indices to be deleted
 * still need to shrink list if possible
 *
 * Revision 1.11  1996/02/23 23:29:29  jean
 * minor twiddles on variable typing
 *
 * Revision 1.10  1996/02/22 23:04:27  jean
 * added NAMD_bsearch and IntList::excise
 * excise is for moving hydrogens to same patch as their mother
 *
 * Revision 1.9  1995/10/27 21:27:20  jim
 * Fixed compile-time error on SGI's introduced in previous revision.
 *
 * Revision 1.8  95/10/27  12:10:33  12:10:33  jim (Jim Phillips)
 * Added several new functions to merge, clear, and compare lists.
 * 
 * Revision 1.7  95/10/09  03:53:18  03:53:18  hazen (Brett Hazen)
 * Memory allocation error-checking added
 * 
 * Revision 1.6  1995/09/21  17:43:33  billh
 * Removed IntList:: from definition of operator[] function.
 *
 * Revision 1.5  95/03/08  14:46:21  14:46:21  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.4  95/02/05  16:13:53  16:13:53  nelson (Mark T. Nelson)
 * Changed binary search from hand written code to use lirbary function bsearch
 * 
 * Revision 1.3  95/01/19  13:32:07  13:32:07  brunner (Robert Brunner)
 * Fixed the delete [] intarray to check for a null pointer first.  Also
 * if no elements are in the array, NOTFOUND is returned for a find().
 * 
 * Revision 1.2  94/10/12  09:59:04  09:59:04  nelson (Mark T. Nelson)
 * Reworked, added sort() and find(), changed the number of
 * elements allocated at a time, made completely inlined,
 * removed merge, and made the reallocation use memcpy rather
 * than loop
 * 
 * Revision 1.1  94/07/05  12:52:16  12:52:16  dalke (Andrew Dalke)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef INTLIST_H
#define INTLIST_H

//  This adds the following names to the global namespace
//    the class IntList

#include "common.h"
#include "Inform.h"
#include <stddef.h>  // for NULL
#include <stdlib.h>  // for qsort and bsearch

#define INTLIST_NOTFOUND	-100	/*  Value for unsuccessful searches	*/
#define SIZE_INCREMENT		5	/*  # of elements to add at each alloc  */

class IntList 
{
  private:
    int count;       // number of elements in the current list
    int maxcount;    // max number of elements allowed before realloc
    int *intarray;   // the list of all integers
    
    IntList( const IntList&) ; 		// prevent copying
    IntList& operator=(const IntList&); // prevent copying
    
  public:
    IntList() 			// init everything
    {       
       count = 0;
       maxcount = 0;
       intarray = NULL;
    } 

    ~IntList()    		// dealloc everything
    {   
      if (intarray) 
        delete [] intarray;
    }

    //  Sort the array in ascending order
    void sort()
    {
	qsort( (void *) intarray, count, sizeof(int), NAMD_compare_ints);
    }

    //  Return the size of the array
    int num() const
    {
      return count;
    }

    //  Add a new element to the array.  If the array has out grown its size,
    //  then allocate more storage.
    void add(int newint) 
    {          	
      if (count==maxcount)  
      {                   	
         int *newarr = new int[maxcount + SIZE_INCREMENT];
	 if ( newarr == NULL )
	 {
	   NAMD_die("memory allocation failed in IntList::add");
	 }
	   

	 memcpy(newarr, intarray, count*sizeof(int));

	 if (intarray != NULL)
	   delete [] intarray;

         intarray = newarr;
         maxcount += SIZE_INCREMENT;
      } 

      intarray[count++] = newint;  // add the new element to the array
    }

    //  Use a binary search to find the index of an element in the
    //  array.  This function assumes that sort() has already been
    //  used to sort the array into ascending order.  If the target
    //  isn't found, then INTLIST_NOTFOUND is returned
    int find(int target)
    {
	int *answer;
	int ret_val;

	if (intarray==NULL)
	{
                ret_val = INTLIST_NOTFOUND;
	}
	else
	{
		answer = (int *) NAMD_bsearch((int *) &target,
			 (int *) intarray, count, sizeof(int), NULL);

		if (answer == NULL)
			ret_val = INTLIST_NOTFOUND;
		else
			ret_val = answer-intarray;
	}

	return(ret_val);
    }

    // get an element using []'s
    const int& operator[](int getint) const
    { 
        if (getint<0 || getint>=count)
		NAMD_die("Subscript out of range in a reference to an IntList\n");

	return intarray[getint];
    }

   //  Add a new value only if it doesn't already exist.
   //  Assumes the list is already sorted, and sorts it when done.
   void merge(int match)
   {
	if ( find(match) == INTLIST_NOTFOUND ) add(match);
	sort();
   }

   //  Ditto for all elements of another IntList (assumed unique).
   //  Assumes the list is already sorted, and sorts it when done.
   void merge(IntList &il)
   {
	int cstart = count;
	int *answer;
	int i;

	if ( cstart )
	    for ( i = 0; i < il.count; ++i )
	    {
		int target = il[i];
		answer = (int *) NAMD_bsearch((int *) &target,
			 (int *) intarray, cstart, sizeof(int), NULL);
		if ( ! answer ) add(il[i]);
	    }
	else
	    for ( i = 0; i < il.count; ++i )
		add(il[i]);

	sort();
   }

   //  Looks for n common atoms between two lists.
   //  Assumes the lists are unique and sorted.
   int hasany(int n, IntList &il)
   {
	int i,j,k;
	int rval = 1;
	i = -1;
	for ( k = 0; k < n && rval; ++k )
	{
		j = il.num();
		while ( j == il.num() && ++i < num() )
			for ( j = 0; j < il.num() && intarray[i] != il[j]; ++j );
		if ( j == il.num() ) rval = 0;
	}
	return rval;
   }

   //  Looks for exactly n common atoms between two lists.
   //  Assumes the lists are unique and sorted.
   int hasexactly(int n, const IntList &il)
   {
	int i,j;
	int c = 0;
	j = 0;
	for ( i = 0; j < il.num() && i < num(); ++i )
	{
		while ( j < il.num() && il[j] < intarray[i] ) ++j;
		if ( j < il.num() && il[j] == intarray[i] ) ++c;
	}
	return ( c == n );
   }

   //  excise a specified number of elements from the list
   int excise(int *indices,int how_many)
   {
	if (how_many <= 0 || how_many > count)
		return 0;
	int target_count = count - how_many;
	int target_chunks = target_count / SIZE_INCREMENT;
	if (target_count % SIZE_INCREMENT)
	    target_chunks++;
	target_count = target_chunks * SIZE_INCREMENT;
	int *newarray = new int[maxcount];
	int k = 0; // index into new array
	int i, j, skip;
	for (i = 0; i < count; i++)
	{
	    skip = 0;
	    for (j= 0; j < how_many; j++)
	    {
		if (indices[j] == i)
		{ // remove this one
		    skip = 1;
		}
	    }
	    if (!skip)
	    {
	        if (k > target_count)
	        {
		    char errmsg[128];
		    sprintf(errmsg,"Buffer overrun in IntList::excise; count = %d, target_count = %d, how_many = %d, k = %d\n",count,target_count,how_many,k);
	  	    NAMD_die(errmsg);
	        }
		newarray[k++] = intarray[i];
	    }
	}
	delete [] intarray;
	intarray = newarray;
	// return number of elements removed
	i = count - k;
	// set new count
	count = k;
	return i;
   }

   //  Removes all elements from the list.
   void clear()
   {
	if (intarray) 
		delete [] intarray;
	intarray = NULL;
	count = 0;
	maxcount = 0;
   }

   //  Just reset the array, do not deallocate anything
   void reset(void)
   {
	count = 0;
   }

   //  Writes out the elements of the list.
   friend Inform& operator<<(Inform& o, const IntList& il)
   {
	for ( int i=0; i<il.num(); ++i)
		o << " " << il[i];
	return o;
   }

};

#endif // INTLIST_H
