/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Generic list of integers implemented as an extendable array.  This
   class is completely inlined.
*/

#ifndef INTLIST_H
#define INTLIST_H

//  This adds the following names to the global namespace
//    the class IntList

#include "common.h"
#include "InfoStream.h"
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
   friend infostream& operator<<(infostream& o, const IntList& il)
   {
	for ( int i=0; i<il.num(); ++i)
		o << " " << il[i];
	return o;
   }

};

#endif // INTLIST_H

