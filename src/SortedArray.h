#ifndef SRTARRAY_H
#define SRTARRAY_H

#include "SortableResizeArray.h"

template <class Elem> class SortedArray: public SortableResizeArray<Elem> {

  protected:

    int found;
    int isSorted;

  public:

    SortedArray(void) : SortableResizeArray<Elem>() { 
      found = -1; isSorted = 1;
    }

    SortedArray(int s) : SortableResizeArray<Elem>(s) { 
      found = -1; isSorted = 1;
    }

    SortedArray(const SortedArray<Elem> &sa) : SortableResizeArray<Elem>(sa) { 
      found = -1; if(!(isSorted = sa.isSorted)) sort();
      isSorted = 1;
    }

    SortedArray(const SortableResizeArray<Elem> &ra) : 
        SortableResizeArray<Elem>(ra) {
      found = -1; sort(); isSorted = 1;
    }

    SortedArray<Elem>& operator =(const SortedArray<Elem> & sa) {
      SortableResizeArray<Elem>::operator=(sa);
      found = -1; isSorted = sa.isSorted;
      return(*this);
    }

    SortedArray<Elem>& operator =(const SortableResizeArray<Elem> &ra) {
      SortableResizeArray<Elem>::operator=(ra);
      found = -1; sort(); isSorted = 1;
      return(*this);
    }

    int load(const Elem& elem) {
      isSorted = 0;
      return(ResizeArray<Elem>::add(elem));
    }

    int add(const Elem& elem) {
      return(insert(elem));
    }

    int del(const Elem & elem) {
      found = bsearch(elem);
      if (size() != 0 && (*rep)[found] == elem) {
        return(SortableResizeArray<Elem>::del(found,1));
      } else {
        found = -1;
        return(-1);
      }
    }

    void sort(void) { SortableResizeArray<Elem>::sort(); isSorted = 1; }

    int bsearch(const Elem& elem) { 
      if (!isSorted) sort();
      return (SortableResizeArray<Elem>::bsearch(elem));
    }

    int insert(const Elem& elem);

    int index(const Elem& elem) { return(found = bsearch(elem)); }

    Elem *find(const Elem& elem);

    Elem * find(void) {
      if (found < 0 || ++found >= size()) {
        return((Elem *)NULL);
      }
      else if ((*rep)[found] == (*rep)[found-1]) {
         return (&(*rep)[found]);
      } else {
        return ((Elem *)NULL);
      }
    }
};

template <class Elem>
int SortedArray<Elem>::insert(const Elem& elem) {
  found = bsearch(elem);
  if (found == -1) {
    return (ResizeArray<Elem>::insert(elem, 0));
  }
  if (found == (size()-1) && (*rep)[found] < elem) {
    return (ResizeArray<Elem>::insert(elem, size()));
  } else {
    return (ResizeArray<Elem>::insert(elem, found));
  }
}

template <class Elem>
Elem * SortedArray<Elem>::find(const Elem& elem) {
  if ( (found = bsearch(elem)) < 0) 
    return ((Elem *)NULL);
  if ((*rep)[found] == elem) {
    return (&(*rep)[found]);
  } else {
    found = -1;
    return ((Elem *)NULL);
  }
}

#endif
