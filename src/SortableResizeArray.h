#ifndef SRTABARRAY_H
#define SRTABARRAY_H

#include <string.h>
#include "ResizeArray.h"

template <class Elem> class SortableResizeArray : public ResizeArray<Elem> {

  private:

    inline void swap(int offset, int i, int j) {
      register Elem *r = (rep->array+offset);
      memcpy((char *)&obj, (char *)&r[i], sizeof(Elem));
      memcpy((char *)&(r[i]), (char *)&(r[j]), sizeof(Elem));
      memcpy((char *)&r[j], (char *)&obj, sizeof(Elem));
    }


    inline void siftup(int offset, int i, int size) {
      register int j;
      register Elem *r = (rep->array+offset);
    
      while ((j = 2*i+1) < size) {
        if (j+1 < size) {
          if (r[j] < r[j+1])
            j = j+1;
          }
          if (r[i] < r[j]) {
            memcpy((char *)&obj, (char *)&r[i], sizeof(Elem));
            memcpy((char *)&(r[i]), (char *)&(r[j]), sizeof(Elem));
            memcpy((char *)&r[j], (char *)&obj, sizeof(Elem));
            i = j;
          } else {
            break;
          }
      }
    }

    Elem obj;

  public:

    SortableResizeArray(void) : ResizeArray<Elem>() { init(); }

    SortableResizeArray(int size) : ResizeArray<Elem>(size) { init(); }

    SortableResizeArray(const ResizeArray<Elem> &ra) : 
      ResizeArray<Elem>(ra) { init(); }

    SortableResizeArray(const SortableResizeArray<Elem> &ra) :
      ResizeArray<Elem>(ra) { init(); }

    SortableResizeArray(ResizeArray<Elem>* const ra) : 
      ResizeArray<Elem>(ra) { init(); }

    SortableResizeArray(SortableResizeArray<Elem>* const ra) : 
      ResizeArray<Elem>(ra) { init(); }

    SortableResizeArray(Elem* * const r, int numElem, int maxElem) : 
      ResizeArray<Elem>(r,numElem,maxElem) { init(); }

    SortableResizeArray<Elem>& operator =(const SortableResizeArray<Elem>& sa) {
      ResizeArray<Elem>::operator=(sa);
      return(*this);
    }

    ~SortableResizeArray(void) { }
  
    void init(void) { }
  
    void sort(void) { sort(0, rep->size()-1); }

    // Heap Sort - worst case is O(n log n)
    //      bot = bottom element of sort range
    //      top = top element of sort range
    void sort(int bot, int top) {
      int index, size;
      if (top > rep->size()) top = rep->size();
      size = top - bot + 1;
    
      // Make all sub-heaps
      for ( index = size/2-1; index > 0; index-- )
        siftup(bot, index, size);
    
      // Take top element of overall heap, and put on top
      for ( index = size; index > 1; index-- ) {
        siftup(bot, 0, index);
        swap(bot, 0, index-1);
      }
    }

    void uniq(void);

    // Search returns index of position where elem should be inserted.
    // This is equal to the first position equal to elem
    // or the first item greater than elem
    // if elem is larger than any item, it returns
    // the index just beyond the end of the list
    // We stick with the < operator only
    int bsearch(const Elem& elem) const {
      int test;
      int bot = -1;
      int top = size();
      if (size() == 0) return (-1);
      while (top - bot > 1) {
        if ( rep->array[test = (bot+top)/2] < elem )
          bot = test;
        else
          top = test;
      }
      return(top);
    }

};

template <class Elem>
void SortableResizeArray<Elem>::uniq(void) {
  if (size()) {
    int oldIndex=0;
    int newIndex=0;
    while (++oldIndex < size()) {
      if ( ! ( rep->array[oldIndex] == rep->array[newIndex] ) ) {
        if (++newIndex != oldIndex)
          memcpy((void *)&(rep->array[newIndex]),
                 (void *)&(rep->array[oldIndex]),
                 sizeof(Elem));
      } else {
        rep->array[oldIndex].~Elem();
      }
    }
    rep->arraySize = ++newIndex;
  }
}

#endif
