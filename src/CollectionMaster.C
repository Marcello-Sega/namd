#include "CollectionMaster.top.h"
#include "CollectionMaster.h"
#include "InfoStream.h"
#include "Node.h"
#include "Output.h"

// #define DEBUGM
#include "Debug.h"

CollectionMaster::CollectionMaster(InitMsg *msg)
{
  delete msg;
}


CollectionMaster::~CollectionMaster(void)
{
}


void CollectionMaster::receivePositions(CollectVectorMsg *msg)
{
  CollectVectorInstance *c;
  if ( c = positions.submitData(msg->seq,msg->aid,msg->data) )
  {
    DebugM(1,"Collected positions at " << c->seq << endl);
    int size = c->data.size();
    Vector *data = new Vector[size];
    for ( int i = 0; i < size; ++i )
    {
      data[c->data[i].aid] = c->data[i].data;
    }
    delete c;
    Node::Object()->output->coordinate(msg->seq,size,data);
    delete data;
  }
  delete msg;
}


void CollectionMaster::receiveVelocities(CollectVectorMsg *msg)
{
  CollectVectorInstance *c;
  if ( c = velocities.submitData(msg->seq,msg->aid,msg->data) )
  {
    DebugM(3,"Collected velocities at " << c->seq << endl);
    int size = c->data.size();
    Vector *data = new Vector[size];
    for ( int i = 0; i < size; ++i )
    {
      data[c->data[i].aid] = c->data[i].data;
    }
    delete c;
    Node::Object()->output->velocity(msg->seq,size,data);
    delete data;
  }
  delete msg;
}


void CollectionMaster::receiveForces(CollectVectorMsg *msg)
{
  CollectVectorInstance *c;
  if ( c = forces.submitData(msg->seq,msg->aid,msg->data) )
  {
    DebugM(3,"Collected forces at " << c->seq << endl);
    int size = c->data.size();
    Vector *data = new Vector[size];
    for ( int i = 0; i < size; ++i )
    {
      data[c->data[i].aid] = c->data[i].data;
    }
    delete c;
    Node::Object()->output->all_force(msg->seq,size,data);
    delete data;
  }
  delete msg;
}


void * CollectVectorMsg::pack(int *length)
{
  *length = sizeof(int) + sizeof(int) +
		aid.size() * sizeof(AtomID) +
		data.size() * sizeof(Vector);
  char *buffer = (char*)new_packbuffer(this,*length);
  char *b = buffer;
  *((int*)b) = seq; b += sizeof(int);
  *((int*)b) = aid.size(); b += sizeof(int);
  for ( int i = 0; i < aid.size(); i++ )
  {
    *((AtomID*)b) = aid[i]; b += sizeof(AtomID);
    *((Vector*)b) = data[i]; b += sizeof(Vector);
  }
  this->~CollectVectorMsg();
  return buffer;
}


void CollectVectorMsg::unpack(void *in)
{
  new((void*)this) CollectVectorMsg;
  char *b = (char*)in;
  seq = *((int*)b); b += sizeof(int);
  int size = *((int*)b); b += sizeof(int);
  aid.resize(size);
  data.resize(size);
  for ( int i = 0; i < size; i++ )
  {
    aid[i] = *((AtomID*)b); b += sizeof(AtomID);
    data[i] = *((Vector*)b); b += sizeof(Vector);
  }
}


#include "CollectionMaster.bot.h"
