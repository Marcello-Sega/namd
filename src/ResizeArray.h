/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   ResizeArray template
   Uses simple contingous array allocation in a hidden manner
   so that array object can have items added without limit
   Suffers from memory fragmentation during resizing
   Fast access, safe and efficient passing of encapsulated array thru
   function arguments.
*/

#ifndef RESIZEARRAY_H
#define RESIZEARRAY_H

#include "ResizeArrayRaw.h"
#include "charm++.h"

// Need this juju to use templated friend below
template <class Type> class ResizeArrayIter;
template <class Type> class ResizeArrayRaw;

template <class Elem> class ResizeArray {
  friend class ResizeArrayIter<Elem>;
  private:
    Elem *secret;
    int   secretArraySize;
    int   secretAllocSize;

  protected:
    ResizeArrayRaw<Elem> *rep;

  public:
    // STL style iterators
    typedef Elem* iterator;
    iterator begin(void) { return rep->array; }
    iterator end(void) { return rep->array + rep->arraySize; }

    // Various Constructors
    ResizeArray(void);

    // Constructor make ResizeArray of predefined size
    ResizeArray(int s) {
      rep = new ResizeArrayRaw<Elem>();
      rep->resize(s);
      rep->refCount = 1;
      secret = NULL;
    }

    // Contructor makes ResizeArray which points to same ResizeArrayRaw
    ResizeArray(const ResizeArray<Elem> &ra) {
      rep = ra.rep;
      rep->refCount++;
      secret = NULL;
    }

    // Constructor makes a copy of ResizeArrayRaw
    ResizeArray(ResizeArray<Elem>* const ra) {
      rep = new ResizeArrayRaw<Elem>(*(ra->rep));
      rep->refCount = 1;
      secret = NULL;
    }

    // Constructor to take-in pre-existing array
    ResizeArray(Elem * * const array, int arraySize, int allocSize=0) {
      rep = new ResizeArrayRaw<Elem>(array, arraySize, allocSize);
      rep->refCount = 1;
      secret = NULL;
    }

    virtual ~ResizeArray(void) {
      if (!--rep->refCount) delete rep;
    }

    // We copy reference to ResizeArrayRaw
    ResizeArray<Elem> & operator= (const ResizeArray<Elem> &ra) {
      secret = NULL;
      if (rep != NULL && !(--rep->refCount) )
        delete rep;
      rep = ra.rep;
      rep->refCount++;
      return (*this);
    }

    // If array is expanded - new elements are default constructed
    // if array is reduced, removed elements have ~Elem() run
    void resize(int i) { rep->resize(i); }

    // Set all elements to a given value (like 0).
    void setall(const Elem &elem) {
      iterator i = begin();
      iterator e = end();
      for ( ; i != e; ++i ) *i = elem;
    }
  
    // Add element to end of array
    int add (const Elem &elem) {
      int end=rep->size();
      rep->ins(elem, end);
      return(end);
    }
  
    // delete num elements from current index
    int del(int index, int num=1) {
      return(rep->del(index,num));
    }

    // delete elements that == e
    virtual int del(const Elem &e) {
      for (int i=0; i < rep->size(); i++) {
        if (rep->array[i] == e) {
          return(del(i,1));
        }
      }
      return(0);
    }

    // insert element at index
    int insert (const Elem& elem, int index) {
      rep->ins(elem,index);
      return (index);
    }

    // array member access (can be lvalue) that grows array.
    inline Elem & item(int i) {
      i = ( i < 0 ? 0 : i );
      if ((i+1) > size())
          resize(i+1);
      return rep->array[i];
    }

    // array member access (can be lvalue) no checks.
    inline Elem & operator[](int index) const { return rep->array[index]; }

    // returns size of ResizeArray
    inline int size(void) const { return rep->size(); }

    // Find and return pointer to element that == e
    Elem * find(const Elem &e) {
      for (int i=0; i < rep->size(); i++) {
        if (rep->array[i] == e) return ((rep->array)+i);
      }
      return NULL;
    }
  
    // Find and return index to element that == e
    int findIndex(const Elem &e) {
      for (int i=0; i < rep->size(); i++) {
        if (rep->array[i] == e) return (i);
      }
      return -1;
    }

    // reduce storage size
    void reduce(void) { rep->reduce(); }

    // Unencap the internal array for use
    // We no-longer control encapsulated array
    Elem * unencap(void) {
      secret = rep->array;
      secretArraySize = rep->arraySize;
      secretAllocSize = rep->allocSize;
      rep->array = NULL;
      rep->varray = NULL;
      rep->arraySize = 0;
      rep->allocSize = 0;
      return secret;
    }

    // Unencap with all the info
    Elem * unencap(int *numElem, int *numAlloc) {
      *numElem = rep->arraySize;
      *numAlloc = rep->allocSize;
      return(unencap());
    }

    // Encap the array back into the ResizeArrayRaw
    void encap(Elem * *const array, int arraySize, int allocSize=0) {
      if (allocSize < arraySize) allocSize = arraySize;
    
      if (secret == *array) {
        rep->arraySize = secretArraySize;
        rep->allocSize = secretAllocSize;
        rep->array = *array;
        rep->varray = (void *)*array;
        secret = NULL;
        *array = NULL;
      } else {
        // Out with the old (but with a proper burial)
        ResizeArrayRaw<Elem> *tmprep
          = new ResizeArrayRaw<Elem>(&secret, secretArraySize, secretAllocSize);
        delete tmprep;
    
        // In with the new
        rep = new ResizeArrayRaw<Elem>(array, arraySize, allocSize);
      }
    }
};

template <class Elem>
ResizeArray<Elem>::ResizeArray(void) {
  rep = new ResizeArrayRaw<Elem>();
  rep->resize(0);
  rep->refCount = 1;
  secret = NULL;
}

#endif
