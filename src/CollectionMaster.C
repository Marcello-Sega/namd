#include "CollectionMaster.top.h"
#include "CollectionMaster.h"


CollectionMaster::CollectionMaster(void)
{
}


CollectionMaster::~CollectionMaster(void)
{
}


void CollectionMaster::receivePositions(CollectVectorMsg *msg)
{
}


void CollectionMaster::receiveVelocities(CollectVectorMsg *msg)
{
}


void CollectionMaster::receiveForces(CollectVectorMsg *msg)
{
}


void * CollectVectorMsg::pack(int *length)
{
  return 0;
}


void CollectVectorMsg::unpack(void *in)
{
}


#include "CollectionMaster.bot.h"
