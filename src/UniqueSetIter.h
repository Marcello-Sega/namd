/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: UniqueSet template - (hashtable)
 *
 ***************************************************************************/

#ifndef USITER_H
#define USITER_H

#include "iostream.h"
#include "UniqueSetRaw.h"
#include "UniqueSet.h"


template <class T> class UniqueSetIter {

  private:

    UniqueSet<T> *us;
    EntryGlob<T> *eg;
    int index;
  
  public:
  
    UniqueSetIter(void) { us = NULL; eg = NULL; index = 0; }

    UniqueSetIter(UniqueSet<T>& us) { 
       this->us = &us; eg = us.rep->globHead; index = 0;
    }

    UniqueSetIter(const UniqueSetIter<T>& iter) {
       us = iter.us; eg = iter.eg; index = iter.index;
    }

    UniqueSetIter<T>& operator= (const UniqueSetIter<T>& iter) {
       us = iter.us; eg = iter.eg; index = iter.index; return (*this);
    }

    ~UniqueSetIter(void) {}
  
    T *operator->(void) { 
      gotoUsed();
      if (eg)
        return (T *)&(eg->glob[index].obj);
      else { 
        index = 0;
        return(NULL);
      }
    }

    UniqueSetIter<T> begin(void) const {
        UniqueSetIter<T> iter;
        iter.us = us;
        iter.index = 0;
        iter.eg = us->rep->globHead;
        iter.gotoUsed();
        return(iter);
    }

    UniqueSetIter<T> end(void) const {
        UniqueSetIter<T> iter;
        iter.us = us;
        iter.index = 0;
        iter.eg = NULL;
        return(iter);
    }
        
    int operator!= (const UniqueSetIter<T> &iter) const {
        return (iter.index != index || iter.eg != eg);
    }

    int operator== (const UniqueSetIter<T> &iter) const {
        return (!operator!=(iter));
    }

    UniqueSetIter<T> operator++(void) {
      index++;
      gotoUsed();
      return (*this);
    }

    UniqueSetIter<T> operator++(int) {
       UniqueSetIter<T> tmp(*this);
       index++;
       gotoUsed();
       return (tmp);
    }
  
    T& operator* (void) {
       gotoUsed();
       return *((T *)&(eg->glob[index].obj));
    }
  
    void status(void) {
      cout << "Index is " << index << " addr is " << eg << endl;
    }
  
    void gotoUsed(void) {
      while (eg) {
        for(;index < us->rep->globSize; index++) {
	  if (eg->glob[index].isUsed()) break;
        }
        if (index < us->rep->globSize) break;
        index = 0;
        eg = eg->next();
      }
      if (!eg) index = 0;
    }
};

#endif
