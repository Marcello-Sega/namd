//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

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

