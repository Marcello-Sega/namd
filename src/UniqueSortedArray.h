#ifndef USRTARRAY_H
#define USRTARRAY_H

#include "SortedArray.h"

template <class Elem> class UniqueSortedArray : public SortedArray<Elem> {

  public:

    UniqueSortedArray(int s=0) : SortedArray<Elem>(s) { }

    UniqueSortedArray(const UniqueSortedArray<Elem> &ua) : 
      SortedArray<Elem>(ua) { }

    UniqueSortedArray(const SortedArray<Elem> &sa) : SortedArray<Elem>(sa) { 
      uniq(); 
    }

    UniqueSortedArray(const ResizeArray<Elem> &ra) : SortedArray<Elem>(ra) {
      uniq();
    }
  
    UniqueSortedArray<Elem>& operator =(const UniqueSortedArray<Elem> & ua) {
      SortedArray<Elem>::operator=(ua);
      return(*this);
    }
  
    UniqueSortedArray<Elem>& operator =(const SortedArray<Elem> &sa) {
      SortedArray<Elem>::operator=(sa);
      uniq();
      return(*this);
    }

    UniqueSortedArray<Elem>& operator =(const ResizeArray<Elem> &ra) {
      SortedArray<Elem>::operator=(ra);
      uniq();
      return(*this);
    }
  
    int add(const Elem& elem) { return(insert(elem)); }

    int insert(const Elem& elem);

};

template <class Elem>
int 
UniqueSortedArray<Elem>::insert(const Elem& elem) {
  found = bsearch(elem);
  if (found == -1) {
    return ResizeArray<Elem>::insert(elem, 0);
  }
  if (found < size() && (*rep)[found] == elem) {
    return -2;
  }
  if (found == (size()-1) && (*rep)[found] < elem) {
    return ResizeArray<Elem>::insert(elem, size());
  } else {
    return ResizeArray<Elem>::insert(elem, found);
  }
}
#endif
