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

#ifndef MAIN_H
#define MAIN_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include <new.h>

#include "BOCgroup.h"
#include "NamdTypes.h"

class GroupInitMsg : public comm_object
{
public:
  BOCgroup group;
};

class EmptyMsg : public comm_object {
  int dummy;
};

class InitMsg : public EmptyMsg { 
};

class RunMsg : public EmptyMsg { 
};

class ReadyMsg : public EmptyMsg { 
};

class DoneMsg : public EmptyMsg { 
};

class RegisterProxyMsg : public comm_object {
public:
  NodeID node;
  PatchID patch;
};

class UnregisterProxyMsg : public comm_object {
public:
  NodeID node;
  PatchID patch;
};

class ProxyAtomsMsg : public comm_object {
public:
  PatchID patch;
  AtomIDList atomIDList;
  void * pack (int *length)
  {
    int size = atomIDList.size();
    *length = 2 * sizeof(int) + size * sizeof(AtomID);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = patch;
    *((int*)(buffer+sizeof(int))) = size;
    AtomID *data = (AtomID*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = atomIDList[i];
    this->~ProxyAtomsMsg();
    return buffer;
  }
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }
  void unpack (void *in)
  {
    new((void*)this) ProxyAtomsMsg;
    char *buffer = (char*)in;
    int patch = *((int*)buffer);
    int size = *((int*)(buffer+sizeof(int)));
    atomIDList.resize(size);
    AtomID *data = (AtomID*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      atomIDList[i] = data[i];
  }
};

class ProxyDataMsg : public comm_object {
public:
  PatchID patch;
  PositionList positionList;
  void * pack (int *length)
  {
    int size = positionList.size();
    *length = 2 * sizeof(int) + size * sizeof(Position);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = patch;
    *((int*)(buffer+sizeof(int))) = size;
    Position *data = (Position*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = positionList[i];
    this->~ProxyDataMsg();
    return buffer;
  }
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t s) { return comm_object::operator new(s); }
  void * operator new(size_t, void *ptr) { return ptr; }
  void unpack (void *in)
  {
    new((void*)this) ProxyDataMsg;
    char *buffer = (char*)in;
    int patch = *((int*)buffer);
    int size = *((int*)(buffer+sizeof(int)));
    positionList.resize(size);
    Position *data = (Position*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      positionList[i] = data[i];
  }
};

class ProxyResultMsg : public comm_object {
public:
  PatchID patch;
  ForceList forceList;
  void * pack (int *length)
  {
    int size = forceList.size();
    *length = 2 * sizeof(int) + size * sizeof(Force);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = patch;
    *((int*)(buffer+sizeof(int))) = size;
    Force *data = (Force*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = forceList[i];
    this->~ProxyResultMsg();
    return buffer;
  }
  void * operator new(size_t s, int i) {return comm_object::operator new(s,i);}
  void * operator new(size_t size) { return comm_object::operator new(size); }
  void * operator new(size_t, void *ptr) { return ptr; }
  void unpack (void *in)
  {
    new((void*)this) ProxyResultMsg;
    char *buffer = (char*)in;
    int patch = *((int*)buffer);
    int size = *((int*)(buffer+sizeof(int)));
    forceList.resize(size);
    Force *data = (Force*)(buffer+2*sizeof(int));
    for ( int i = 0; i < size; ++i )
      forceList[i] = data[i];
  }
};

class Compute;

class LocalWorkMsg : public comm_object
{
public:
  Compute *compute;
};

#endif /* MAIN_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: main.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.18 $	$Date: 1996/12/17 08:54:40 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.h,v $
 * Revision 1.18  1996/12/17 08:54:40  jim
 * fixed several bugs, including not saving patch
 *
 * Revision 1.17  1996/12/16 22:52:43  jim
 * added placement new and explicit destructor calls to ProxyAtomsMsg
 *
 * Revision 1.16  1996/12/16 22:19:26  jim
 * added placement new and destructor to messages
 *
 * Revision 1.15  1996/12/14 00:02:42  jim
 * debugging ProxyAtomsMsg path to make compute creation work
 *
 * Revision 1.14  1996/12/05 20:26:27  jim
 * missing ;
 *
 * Revision 1.13  1996/12/05 20:25:03  jim
 * forgot to fill in length of message
 *
 * Revision 1.12  1996/12/05 17:56:58  jim
 * added pack and unpack functions for proxy messages
 *
 * Revision 1.11  1996/12/05 17:15:23  jim
 * added proxy messages
 *
 * Revision 1.10  1996/12/05 17:00:05  ari
 * *** empty log message ***
 *
 * Revision 1.9  1996/12/05 01:47:40  ari
 * added messages for proxy management
 *
 * Revision 1.8  1996/11/22 00:18:51  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/09/15 20:43:11  jim
 * fixed missing semicolon
 *
 * Revision 1.6  1996/09/10 04:45:06  ari
 * added LocalWorkMsg
 *
 * Revision 1.5  1996/08/29 00:50:42  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/22 16:15:07  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/19 22:05:31  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 20:52:30  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/02 19:20:13  gursoy
 * Initial revision
 *
 ***************************************************************************/
