#ifndef COLLECTIONMGR_H
#define COLLECTIONMGR_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"


class CollectionMgr : public BOCclass
{
public:

  CollectionMgr(SlaveInitMsg *msg);
  ~CollectionMgr(void);

  void submitPositions(AtomIDList &i, PositionList &d);
  void submitVelocities(AtomIDList &i, VelocityList &d);
  void submitForces(AtomIDList &i, ForceList &d);

private:

};


#endif
