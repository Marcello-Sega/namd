/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef USRTARRAY_H
#define USRTARRAY_H

#include "SortedArray.h"

template <class Elem> class UniqueSortedArray : public SortedArray<Elem> {

  public:

    UniqueSortedArray(int s=0) : SortedArray<Elem>(s) { }

    UniqueSortedArray(const UniqueSortedArray<Elem> &ua) : 
      SortedArray<Elem>(ua) { }

    UniqueSortedArray(const SortedArray<Elem> &sa) : SortedArray<Elem>(sa) { 
      this->uniq(); 
    }

    UniqueSortedArray(const ResizeArray<Elem> &ra) : SortedArray<Elem>(ra) {
      this->uniq();
    }
  
    UniqueSortedArray<Elem>& operator =(const UniqueSortedArray<Elem> & ua) {
      SortedArray<Elem>::operator=(ua);
      return(*this);
    }
  
    UniqueSortedArray<Elem>& operator =(const SortedArray<Elem> &sa) {
      SortedArray<Elem>::operator=(sa);
      this->uniq();
      return(*this);
    }

    UniqueSortedArray<Elem>& operator =(const ResizeArray<Elem> &ra) {
      SortedArray<Elem>::operator=(ra);
      this->uniq();
      return(*this);
    }
  
    int add(const Elem& elem) { return(insert(elem)); }

    int insert(const Elem& elem);

};

template <class Elem>
int 
UniqueSortedArray<Elem>::insert(const Elem& elem) {
  this->found = bsearch(elem);
  if (this->found == -1) {
    return ResizeArray<Elem>::insert(elem, 0);
  }
  if (this->found < this->size() && (*(this->rep))[this->found] == elem) {
    return -2;
  }
  if (this->found == (this->size()-1) && (*(this->rep))[this->found] < elem) {
    return ResizeArray<Elem>::insert(elem, this->size());
  } else {
    return ResizeArray<Elem>::insert(elem, this->found);
  }
}
#endif
