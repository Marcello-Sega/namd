#ifndef _NAMDNBORLB_H_
#define _NAMDNBORLB_H_

#include "NeighborLB.h"
#include "NamdNborLB.decl.h"

#include "Node.h"
#include "PatchMap.h"
#include "SimParameters.h"
#include "AlgNbor.h"
#include "InfoStream.h"

void CreateNamdNborLB();

class NamdNborLB : public NeighborLB {

public:
  NamdNborLB();
private:
  CmiBool QueryBalanceNow(int step);
  NLBMigrateMsg* Strategy(NborBaseLB::LDStats* stats, int count);
  int buildData(NborBaseLB::LDStats* stats, int count);
  int requiredProxies(PatchID id, int neighborNodes[]);

  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;
};

#endif /* _NAMDCENTLB_H_ */
