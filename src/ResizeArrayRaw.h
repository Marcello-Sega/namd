//-*-c++-*-
/************************************************************************/
/*                                                                      */
/*              (C) Copyright 1996 The Board of Trustees of the         */
/*                          University of Illinois                      */
/*                           All Rights Reserved                        */
/*                                                                      */
/************************************************************************/

/************************************************************************
 * DESCRIPTION: ResizeArrayRaw template					*
 * Object Requirements: new                                             *
 *                      Elem(Elem &)                                    *
 *			~Elem()						*
 *                      Elem & operator= (Elem &)                       *
 ************************************************************************/


#ifndef RESIZEARRAYRAW_H
#define RESIZEARRAYRAW_H

#include <new.h>
#include <stream.h>
#include <string.h>
#include "ckdefs.h"

// Undefined if src,dest overlap
extern "C" void *memcpy(void *, const void *, size_t);
// Use memmove on possibly overlapping memory moves.
extern "C" void *memmove(void *, const void *, size_t);

#define GrowthFactor 1.5
#define MinSize 8

// Need this juju to use templated friend below
template <class Type> class ResizeArray;
template <class Type> class SortableResizeArray;
template <class Type> class ResizeArrayIter;

// Class assumes that one can bit move objects
// around on array.  This will be true
// as long as object has no pointers to itself.
template <class Elem> class ResizeArrayRaw {

  private:
    Elem *array;
    void *varray;

    int arraySize;
    int allocSize;

    int refCount;

    float growthFactor;
    int minSize;

    // No constructor run on new elements
    // arraySize is not adjusted, only allocSize
    void resizeRaw(int size) {
      if (size <= allocSize) return;
  
      if (size < (int)(allocSize*growthFactor)) 
        size = (int)(allocSize*growthFactor);
      if ( (size-allocSize) < minSize) 
        size = allocSize+minSize;
  
      void *tmp = (void *) new char[size*sizeof(Elem)];
      memcpy(tmp, varray, sizeof(Elem)*arraySize);
  
      if (allocSize) 
        delete[] varray;
      array = (Elem *)(varray = tmp);
      allocSize = size;
    }

    // Empty function for now.
    // eventually, this should get smaller storage and free
    void reduce(void) {}; 

  public:
    friend class ResizeArray<Elem>;
    friend class SortableResizeArray<Elem>;
    friend class ResizeArrayIter<Elem>;

    inline int size(void) const { return arraySize; }
    inline Elem &operator[](int index) const { return array[index]; }

    // Default constructor 
    ResizeArrayRaw(void) : 
      array((Elem *)0), varray((void *)0), arraySize(0), allocSize(0) { 
      growthFactor = GrowthFactor;
      minSize = MinSize;
    }

    // Copy constructor - true copy on construction.
    ResizeArrayRaw(const ResizeArrayRaw<Elem> &rar ) {
      growthFactor = rar.growthFactor;
      minSize = rar.minSize;
      // We want rar.size() slots, but no constructor run on the elements
      resizeRaw(rar.size());
      memcpy(varray, rar.varray, sizeof(Elem)*rar.size());
      arraySize = rar.size();
    }
  
    // Encap a pre-existing array
    ResizeArrayRaw( Elem * * const array, int arraySize, int allocSize) {
      if (allocSize < arraySize) allocSize = arraySize;
      this->allocSize = allocSize;
      this->arraySize = arraySize;
      varray = *array;
      this->array = (Elem *)*array;
      *array = 0;
    }
  
    ~ResizeArrayRaw(void);
  
    // minSize = minimum growth size - (also initial size of array)
    // growthFactor = mulplicative factor by which to grow array.
    void setResizeParams(int min, float growth) {
      minSize = min;
      growthFactor = growth;
    }
  
  
    // True copy made on assignment.
    ResizeArrayRaw<Elem> & operator=(const ResizeArrayRaw<Elem> &rar ) {
      growthFactor = rar.growthFactor;
      minSize = rar.minSize;
  
      // Clean up this array
      resize(0);
      resizeRaw(rar.size());
  
      memcpy(varray, rar.varray, sizeof(Elem)*rar.size());
      arraySize = rar.size();
      return *this;
    }
  
    // Properly constructs default object on new elements
    // Properly destructs on removed elements
    // arraySize is properly updated
    void resize(int size) {
      int i;
  
      if (size < arraySize) {
        for (i=size; i<arraySize; i++) {
          array[i].~Elem();
        }
      } else if (size > arraySize) {
        resizeRaw(size);
        for (i=arraySize; i<size; i++) {
          new ((void *)&array[i]) Elem;
        }
      }
      arraySize = size;
    }
  
    inline int del(int index, int number) {
      int i;
  
      // Fix up number to delete if deletes off end of array
      if (index >= arraySize)
        number=0; // for inline sake, don't have multiple returns
      else if (index+number-1 > arraySize) {
        number = index-arraySize+1;
      }
  
      // Destruct objects to be deleted
      for (i=index; i < index+number; i++) {
        array[i].~Elem();
      }
  
      // Shift down
      memmove((void *)((char *)varray+index*sizeof(Elem)),
         (void *)((char *)varray+(index+number)*sizeof(Elem)),
         (arraySize-number-index)*sizeof(Elem));
      
      // fixup size of array
      arraySize -= number;
      return(number);
    }
  
    
    // Insert element in array
    // If index is over the end of array, default
    // constructor should be run over blank elements to pad.
    inline void ins(const Elem &e, int index) {
      // Size array depending if index is in current array or reaches beyond.
      if (index < arraySize) {
        resizeRaw(arraySize+1);
        // Shift up
        memmove((void *)((char *)varray+(index+1)*sizeof(Elem)),
          (void *)((char *)varray+index*sizeof(Elem)),
          (arraySize-index)*sizeof(Elem));
      } else {
        resizeRaw(index+1);
      }
      
      // Write in new element via assignment - allows any refcounting
      // etc. to take place correctly!
      new((void *)&array[index]) Elem;
      array[index] = e;
    
      // Take care of fill and setting correct arraySize 
      if (index > arraySize) {
        for (Elem *tmp = array+arraySize; tmp < array+index; tmp++) {
          new ((void *)tmp) Elem;
        }
        arraySize = index+1;
      } else
        arraySize++;
    }
};	// end template definition

// If ~Elem() is done properly, on light classes, the loop should
// optimize to nothing.
template <class Elem>
ResizeArrayRaw<Elem>::~ResizeArrayRaw() {
  for (int i=0; i < size(); i++) {
    array[i].~Elem();
  }
  delete[] varray;
}
#endif
