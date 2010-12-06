/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/NamdHybridLB.h,v $
 * $Author: bhatele $
 * $Date: 2010/12/06 21:51:22 $
 * $Revision: 1.8 $
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
  //void CollectInfo(Location *loc, int n, int fromlevel);

private:
  CProxy_NamdHybridLB thisProxy;
  int updateCount;
  bool collectFlag;
  bool updateFlag;
  int parent_backup;
  Location *loc_backup;
  int n_backup;
  int fromlevel_backup;

  int *from_procs;
  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;

  CmiBool QueryBalanceNow(int step);
  CmiBool QueryDumpData();
  // LBVectorMigrateMsg* VectorStrategy(LDStats* stats);

#if CHARM_VERSION > 60301
  CLBMigrateMsg* Strategy(LDStats* stats);
  LBMigrateMsg* GrpLevelStrategy(LDStats* stats);
#else
  CLBMigrateMsg* Strategy(LDStats* stats, int n_pes);
  LBMigrateMsg* GrpLevelStrategy(LDStats* stats, int n_pes);
#endif
  
  int buildData(LDStats* stats);
  int requiredProxies(PatchID id, int neighborNodes[]);
  void dumpDataASCII(char *file, int numProcessors, int numPatches,
                int numComputes);

  // centralized load balancer for load balancing all the children processors
  NamdCentLB *centralLB;
  NamdDummyLB *dummyLB;
};

#endif /* _NAMDHYBRIDLB_H_ */
