#ifndef COLLECTIONMASTER_H
#define COLLECTIONMASTER_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"
#include "NamdTypes.h"

class CollectVectorMsg;

class CollectionMaster : public chare_object
{
public:

  CollectionMaster(void);
  ~CollectionMaster(void);

  void receivePositions(CollectVectorMsg *msg);
  void receiveVelocities(CollectVectorMsg *msg);
  void receiveForces(CollectVectorMsg *msg);

private:



};


class CollectVectorMsg : public comm_object
{
public:

  void * pack(int *length);
  void unpack(void *in);

};



#endif
