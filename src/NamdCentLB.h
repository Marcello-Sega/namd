#ifndef _NAMDCENTLB_H_
#define _NAMDCENTLB_H_

#include <CentralLB.h>
#include "NamdCentLB.decl.h"

#include "Node.h"
#include "PatchMap.h"
#include "SimParameters.h"
#include "RefineOnly.h"
#include "Alg7.h"
#include "InfoStream.h"

void CreateNamdCentLB();

class NamdCentLB : public CentralLB {

public:
  NamdCentLB();
private:
  CmiBool QueryBalanceNow(int step);
  CLBMigrateMsg* Strategy(CentralLB::LDStats* stats, int count);
  int buildData(CentralLB::LDStats* stats, int count);
  int requiredProxies(PatchID id, int neighborNodes[]);

  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;
};

#endif /* _NAMDCENTLB_H_ */
