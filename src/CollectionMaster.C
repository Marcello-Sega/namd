#include "CollectionMaster.top.h"
#include "CollectionMaster.h"
#include "InfoStream.h"
#include "Node.h"
#include "Output.h"
#include "ProcessorPrivate.h"

// #define DEBUGM
#include "Debug.h"

CollectionMaster::CollectionMaster(InitMsg *msg)
{
  delete msg;
  if (CpvAccess(CollectionMaster_instance) == 0) {
    CpvAccess(CollectionMaster_instance) = this;
  } else {
    DebugM(1, "CollectionMaster::CollectionMaster() - another instance of CollectionMaster exists!\n");
  }
}


CollectionMaster::~CollectionMaster(void)
{
}


void CollectionMaster::receivePositions(CollectVectorMsg *msg)
{
  CollectVectorInstance *c;
  if ( c = positions.submitData(msg->seq,msg->aid,msg->data) )
  {
    disposePositions(c);
  }
  delete msg;
}

void CollectionMaster::enqueuePositions(int seq)
{
  CollectVectorInstance *c;
  if ( c = positions.enqueue(seq) )
  {
    disposePositions(c);
  }
}

void CollectionMaster::disposePositions(CollectVectorInstance *c)
{
    DebugM(3,"Collected positions at " << c->seq << endl);
    int seq = c->seq;
    int size = c->data.size();
    Vector *data = new Vector[size];
    for ( int i = 0; i < size; ++i )
    {
      data[c->data[i].aid] = c->data[i].data;
    }
    delete c;
    Node::Object()->output->coordinate(seq,size,data);
    delete [] data;
}


void CollectionMaster::receiveVelocities(CollectVectorMsg *msg)
{
  CollectVectorInstance *c;
  if ( c = velocities.submitData(msg->seq,msg->aid,msg->data) )
  {
    disposeVelocities(c);
  }
  delete msg;
}

void CollectionMaster::enqueueVelocities(int seq)
{
  CollectVectorInstance *c;
  if ( c = velocities.enqueue(seq) )
  {
    disposeVelocities(c);
  }
}

void CollectionMaster::disposeVelocities(CollectVectorInstance *c)
{
    DebugM(3,"Collected velocities at " << c->seq << endl);
    int seq = c->seq;
    int size = c->data.size();
    Vector *data = new Vector[size];
    for ( int i = 0; i < size; ++i )
    {
      data[c->data[i].aid] = c->data[i].data;
    }
    delete c;
    Node::Object()->output->velocity(seq,size,data);
    delete data;
}


void CollectionMaster::receiveForces(CollectVectorMsg *msg)
{
  CollectVectorInstance *c;
  if ( c = forces.submitData(msg->seq,msg->aid,msg->data) )
  {
    disposeForces(c);
  }
  delete msg;
}

void CollectionMaster::enqueueForces(int seq)
{
  CollectVectorInstance *c;
  if ( c = forces.enqueue(seq) )
  {
    disposeForces(c);
  }
}

void CollectionMaster::disposeForces(CollectVectorInstance *c)
{
    DebugM(3,"Collected forces at " << c->seq << endl);
    int size = c->data.size();
    Vector *data = new Vector[size];
    for ( int i = 0; i < size; ++i )
    {
      data[c->data[i].aid] = c->data[i].data;
    }
    delete c;
    Node::Object()->output->all_force(c->seq,size,data);
    delete data;
}


void * CollectVectorMsg::pack(int *length)
{
  *length = sizeof(int) + sizeof(int) +
		aid.size() * sizeof(AtomID) +
		data.size() * sizeof(Vector);
  char *buffer = (char*)new_packbuffer(this,*length);
  char *b = buffer;
  memcpy(b, &seq, sizeof(int)); b += sizeof(int);
  int size = aid.size(); memcpy(b, &size, sizeof(int)); b += sizeof(int);
  memcpy(b, aid.begin(), size*sizeof(AtomID)); b += size*sizeof(AtomID);
  memcpy(b, data.begin(), size*sizeof(Vector)); b += size*sizeof(Vector);
  this->~CollectVectorMsg();
  return buffer;
}


void CollectVectorMsg::unpack(void *in)
{
  new((void*)this) CollectVectorMsg;
  char *b = (char*)in;
  memcpy(&seq, b, sizeof(int)); b += sizeof(int);
  int size; memcpy(&size, b, sizeof(int)); b += sizeof(int);
  aid.resize(size);
  memcpy(aid.begin(), b, size*sizeof(AtomID)); b += size*sizeof(AtomID);
  data.resize(size);
  memcpy(data.begin(), b, size*sizeof(Vector)); b += size*sizeof(Vector);
}


#include "CollectionMaster.bot.h"


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1012 $	$Date: 1997/12/19 23:42:38 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: CollectionMaster.C,v $
 * Revision 1.1012  1997/12/19 23:42:38  jim
 * Replaced assignments with memcpys and reordered memcpys for efficiency.
 *
 * Revision 1.1011  1997/12/10 17:53:33  milind
 * Removed the dcd file already exists error. Now, if a dcd file already exists,
 * it is moved to a .bak before writing new dcd file.
 *
 * Revision 1.1010  1997/12/02 22:04:00  milind
 * Fixed a silly bug that accessed memory after it was freed.
 *
 * Revision 1.1009  1997/11/07 20:17:33  milind
 * Made NAMD to run on shared memory machines.
 *
 * Revision 1.1008  1997/04/04 23:34:13  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1007  1997/03/19 11:53:56  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 *
 ***************************************************************************/
