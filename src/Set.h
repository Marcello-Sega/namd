/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef IRSET_DEFS_H
#define IRSET_DEFS_H
class InfoRecord;


class listNode {
public:
listNode *next;
InfoRecord *const info;
listNode(InfoRecord *i) : info(i) {;}
};

class Iterator{
public:
  int id; // for debugging
  listNode* next;
};

class IRSet {

private:
 listNode *head;

public:
 IRSet();
 ~IRSet();
 void insert(InfoRecord *);
 int find(InfoRecord *) ;
 int remove(InfoRecord *);
 void myRemove(listNode **n, InfoRecord *r);
 InfoRecord *iterator(Iterator *);
 InfoRecord *next(Iterator *);
 int numElements();
 void print();
};

#endif
