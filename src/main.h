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
  AtomIDList *atomIDList;
  ProxyAtomsMsg(PatchID pid, AtomIDList a) : patch(pid)
  {
    atomIDList = new AtomIDList(a);
  }
  ~ProxyAtomsMsg()
  {
    delete atomIDList;
  }
  void * pack (int *length)
  {
    int size = atomIDList->size();
    *length = sizeof(int) + size * sizeof(AtomID);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = size;
    AtomID *data = (AtomID*)(buffer+sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = (*atomIDList)[i];
    atomIDList = NULL;
    return buffer;
  }
  void unpack (void *in)
  {
    char *buffer = (char*)in;
    int size = *((int*)buffer);
    atomIDList = new AtomIDList; atomIDList->resize(size);
    AtomID *data = (AtomID*)(buffer+sizeof(int));
    for ( int i = 0; i < size; ++i )
      (*atomIDList)[i] = buffer[i];
  }
};

class ProxyDataMsg : public comm_object {
public:
  PatchID patch;
  PositionList positionList;
  void * pack (int *length)
  {
    int size = positionList.size();
    *length = sizeof(int) + size * sizeof(Position);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = size;
    Position *data = (Position*)(buffer+sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = positionList[i];
    return buffer;
  }
  void unpack (void *in)
  {
    char *buffer = (char*)in;
    int size = *((int*)buffer);
    positionList.resize(size);
    Position *data = (Position*)(buffer+sizeof(int));
    for ( int i = 0; i < size; ++i )
      positionList[i] = buffer[i];
  }
};

class ProxyResultMsg : public comm_object {
public:
  PatchID patch;
  ForceList forceList;
  void * pack (int *length)
  {
    int size = forceList.size();
    *length = sizeof(int) + size * sizeof(Force);
    char *buffer = (char*)new_packbuffer(this,*length);
    *((int*)buffer) = size;
    Force *data = (Force*)(buffer+sizeof(int));
    for ( int i = 0; i < size; ++i )
      data[i] = forceList[i];
    return buffer;
  }
  void unpack (void *in)
  {
    char *buffer = (char*)in;
    int size = *((int*)buffer);
    forceList.resize(size);
    Force *data = (Force*)(buffer+sizeof(int));
    for ( int i = 0; i < size; ++i )
      forceList[i] = buffer[i];
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
 *	$Revision: 1.15 $	$Date: 1996/12/14 00:02:42 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.h,v $
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
