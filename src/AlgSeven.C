#include <iostream.h>
#include "Alg7.h"

Alg7::Alg7(computeInfo *computeArray, patchInfo *patchArray, 
	   processorInfo *processorArray, int nComps, 
	   int nPatches, int nPes) :
  Rebalancer(computeArray, patchArray, 
	     processorArray, nComps, 
	     nPatches, nPes)
{
strategyName = "Alg7";
strategy();
}

void Alg7::strategy()
{
   double bestSize0, bestSize1, bestSize2;
   computeInfo *c;
   processorInfo *bestP0, *bestP1, *bestP2;
   bestP0 = bestP1 = bestP2 = (processorInfo *)0;
   makeHeaps();
   computeAverage();
   printLoads();
	      
  for (int i=0; i<numComputes; i++){
     c = (computeInfo *) computesHeap->deleteMax();
    heapIterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      pes->iterator((heapIterator *) &nextProcessor);
    bestSize0 = bestSize1 = bestSize2 = 0;
    bestP0 = p;
    while (p){
       int n=0;
       n = numAvailable(c,p);
       switch(n){
       case 0:
	 if (!bestP0){
	   if (c->load + p->load < overLoad*averageLoad)
	       bestP0 = p;
	 }
	 else
	   if ( ( c->load + p->load < overLoad*averageLoad) && 
		(p->load<bestP0->load))
	     bestP0 = p;
     
       break;
       case 1:
	 if (!bestP1){
	   if (c->load + p->load < overLoad*averageLoad)
	       bestP1 = p;
	 }
	 else
	   if (( c->load + p->load < overLoad*averageLoad) &&
		    (p->load < bestP1->load))
      	 bestP1 = p;
       break;
       case 2: 
	 if (!bestP2){
	   if (c->load + p->load < overLoad*averageLoad)
	       bestP2 = p;
	 }
	 else
	   if (( c->load + p->load < overLoad*averageLoad) &&
		    (p->load < bestP2->load))
	 bestP2 = p;
       break;
       default: cout << "Error. Illegal number of proxies.\n";    
       }
     p = (processorInfo *) pes->next(&nextProcessor);
    }
    
    if ((bestP2) && (bestP2->load < 1.2*bestP0->load)){
      assign(c, bestP2);
    }
    else if ((bestP1) && (bestP1->load < 1.2*bestP0->load)){
      assign(c, bestP1);
    }
    else if (bestP0){
     assign(c, bestP0);
    }
    else { 
      cout << "No receiver found" << "\n";
      break;
    }

  }
  printLoads();
  overLoad = 1.01;
  cout << "Starting overLoad = " << overLoad << endl;
  for (; !refine(); overLoad += .01);
  cout << "Ending overLoad = " << overLoad << endl;
  printLoads();
}










