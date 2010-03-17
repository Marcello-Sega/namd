/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/NamdHybridLB.h,v $
 * $Author: gzheng $
 * $Date: 2010/03/17 07:48:26 $
 * $Revision: 1.4 $
 *****************************************************************************/

#ifndef _NAMDHYBRIDLB_H_
#define _NAMDHYBRIDLB_H_

#include <HybridBaseLB.h>
#include "NamdHybridLB.decl.h"

#include "Node.h"
#include "PatchMap.h"
#include "SimParameters.h"
#include "RefineOnly.h"
#include "Alg7.h"
#include "AlgRecBisection.h"
#include "InfoStream.h"
#include "NamdCentLB.h"
#include "NamdDummyLB.h"
#include "TorusLB.h"
#include "RefineTorusLB.h"

void CreateNamdHybridLB();

class NamdHybridLB : public HybridBaseLB {

public:
  NamdHybridLB();
  NamdHybridLB(CkMigrateMessage *m):HybridBaseLB(m) {}
  void UpdateComputeMap(CLBMigrateMsg *msg);
  void CollectInfo(Location *loc, int n, int fromlevel);

private:
  CProxy_NamdHybridLB thisProxy;
  int updateCount;
  bool collectFlag;
  bool updateFlag;
  int parent_backup;
  Location *loc_backup;
  int n_backup;
  int fromlevel_backup;

  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;

  CmiBool QueryBalanceNow(int step);
  CmiBool QueryDumpData();
//  LBVectorMigrateMsg* VectorStrategy(LDStats* stats,int count);
  LBMigrateMsg* Strategy(LDStats* stats,int count);
  LBMigrateMsg* GrpLevelStrategy(LDStats* stats,int count);
  
  int buildData(CentralLB::LDStats* stats, int count);
  int requiredProxies(PatchID id, int neighborNodes[]);
  void dumpDataASCII(char *file, int numProcessors, int numPatches,
                int numComputes);

  // centralized load balancer for load balancing all the children processors
  NamdCentLB *centralLB;
  NamdDummyLB *dummyLB;
};

#endif /* _NAMDHYBRIDLB_H_ */
