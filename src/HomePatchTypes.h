/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef HOMEPATCHTYPES_H
#define HOMEPATCHTYPES_H

#include "PatchTypes.h"
#include "Box.h"

class ProxyListElem {
public:
  ProxyListElem() {};
  ProxyListElem(NodeID n) : node(n), forceBox(0) {};
  ProxyListElem(NodeID n, Box<Patch,Results> *f ) : node(n), forceBox(f) {};
  ~ProxyListElem() {};

  int operator==(const ProxyListElem &p) {
    return (node == p.node);
  }

  NodeID node;
  Box<Patch,Results> *forceBox;
};

typedef ResizeArray<ProxyListElem> ProxyList;
typedef ResizeArrayIter<ProxyListElem> ProxyListIter;


#endif

