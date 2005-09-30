#ifndef _NAMDCENTLB_H_
#define _NAMDCENTLB_H_

#include <CentralLB.h>
#include "NamdCentLB.decl.h"

#include "Node.h"
#include "PatchMap.h"
#include "SimParameters.h"
#include "RefineOnly.h"
#include "Alg7.h"
#include "AlgRecBisection.h"
#include "InfoStream.h"

void CreateNamdCentLB();

class NamdCentLB : public CentralLB {

public:
  NamdCentLB();
private:
  CmiBool QueryBalanceNow(int step);
  CmiBool QueryDumpData();
  CLBMigrateMsg* Strategy(CentralLB::LDStats* stats, int count);
  int buildData(CentralLB::LDStats* stats, int count);
  int requiredProxies(PatchID id, int neighborNodes[]);
#if CMK_VERSION_BLUEGENE
  int requiredProxiesOnProcGrid(PatchID id, int neighborNodes[]);
#endif
  void dumpData(char *file, int numProcessors, int numPatches,
		int numComputes);
  void loadData(char *file, int &numProcessors, int &numPatches,
		int &numComputes);
  void dumpDataASCII(char *file, int numProcessors, int numPatches,
		int numComputes);
  void loadDataASCII(char *file, int &numProcessors, int &numPatches,
		int &numComputes);

  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;
};

#endif /* _NAMDCENTLB_H_ */
