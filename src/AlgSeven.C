#include <iostream.h>
#include <InfoStream.h>
#include "Alg7.h"

#define TINYLOAD 0.0005

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
   int numAssigned;
   processorInfo *bestP, *bestP0, *bestP1, *bestP2;
   bestP0 = bestP1 = bestP2 = (processorInfo *)0;
   //   iout << iINFO  << "calling makeHeaps. \n";
   makeHeaps();
   computeAverage();
   //   iout << iINFO
   //	<< "Before assignment\n" << endi;
   //   printLoads();
	      
   numAssigned = 0;
   int numAssigned1 = 0;
   int numAssigned2 = 0;
   int numAssignedP2 = 0;
   int numAssignedP1 = 0;
   int numAssignedP0 = 0;
   int numAssignedP4 = 0;

   //   for (int i=0; i<numPatches; i++)
   //     { cout << "(" << patches[i].Id << "," << patches[i].processor ;}
   int i;
#if 0
   for (i=0; i<numComputes; i++){
     /* if compute i is a self-interaction, assign it to the home processor.
	Also, if compute i represents load below a threshold ("zero load")
	assign it to the old processor it was on */
    
     computeInfo * c = (computeInfo *) &(computes[i]);
     if (c->patch1 == c->patch2)
       {
	 //	 assign(c, patches[c->patch1].processor);
	 assign(c, c->oldProcessor);

	 /*
	   cout << "assigning load " << c->load << "," << "to processor:" <<
	   patches[c->patch1].processor << "because of patch:" << c->patch1 
	      << "patchID" << patches[c->patch1].Id 
	      << endl; 
	      */

	 numAssigned++;
	 numAssigned1++;
       }
     else if (c->load < TINYLOAD) {
       assign(c, c->oldProcessor);
       numAssigned++;
       numAssigned2++;
     }
   }
#endif
   //   iout << iINFO  << numAssigned <<  "done initial assignments.\n";

   //   printLoads();
   overLoad = 1.7;
  for (i=0; i<numComputes; i++){
     c = (computeInfo *) computesHeap->deleteMax();
     if (c->processor != -1) continue; // skip to the next compute;
     heapIterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      pes->iterator((heapIterator *) &nextProcessor);
    bestSize0 = bestSize1 = bestSize2 = 0;
    bestP0 = bestP1 = bestP2 = (processorInfo *)0;
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
       default:
	 iout << iINFO  << "Error. Illegal number of proxies.\n" << endi;    
       }
     p = (processorInfo *) pes->next(&nextProcessor);
    }

    if (numAssigned >= 0)
    {
      if (bestP2)
      {
//      if ((bestP0==NULL) || (bestP2->load < 1.2*bestP0->load)) {
	assign(c, bestP2);
	numAssigned++;
	numAssignedP2++;
//      }
      }
      else if (bestP1)
      {
//      if ((bestP0==NULL) || (bestP1->load < 1.2*bestP0->load)){
	assign(c, bestP1);
	numAssigned++;
	numAssignedP1++;
//      }
      }
      else if (bestP0){
        assign(c, bestP0);
        numAssigned++;
        numAssignedP0++;
      }
      else { 
        iout << iINFO  << "Alg 7 No receiver found 1" << "\n" <<endi;
        break;
      }
    }
    else 
    {
      // At start, load is most important, rather than communications
      int *numAssignedptr = &numAssignedP2;
      bestP = bestP2;
      if (!bestP || (bestP1 && (bestP1->load < 0.8 * bestP->load)) )
      {
	bestP=bestP1;
	numAssignedptr = &numAssignedP1;
      }
      if (!bestP || (bestP0 && (bestP0->load < 0.75 * bestP->load)))
      {
	bestP=bestP0;
	numAssignedptr = &numAssignedP0;
      }
      if (!bestP)
      {
	iout << iINFO  << "Alg7 No receiver found 2" << "\n" <<endi;
	break;
      }
      assign(c,bestP);
      (*numAssignedptr)++;
      numAssigned++;
      numAssignedP4++;
      
    }
  }

#ifdef DEBUG
  iout << iINFO
       << "numAssigned = " << numAssigned
       << "\nnumAssigned1 = " << numAssigned1
       << "\nnumAssigned2 = " << numAssigned2
       << "\nnumAssignedP2 = " << numAssignedP2
       << "\nnumAssignedP1 = " << numAssignedP1
       << "\nnumAssignedP0 = " << numAssignedP0
       << "\nnumAssignedP4 = " << numAssignedP4
       << "\n" << endi;
#endif

//

//  p = (processorInfo *) pes->iterator((heapIterator *) &nextProcessor);
//  while(p)
 // {
//    iout << iINFO 
//	 << "ID = " << p->Id
//	 << "  NComputes = " << p->computeSet->numElements()
//	 << "\n" << endi;
//     p = (processorInfo *) pes->next(&nextProcessor);
//  }
    
  // printLoads();
  overLoad = 1.02;
  //  iout << iINFO  << "num assigned: " << numAssigned << endi;
  //  iout << iINFO  << "Starting overLoad = " << overLoad << endi;
  for (; !refine(); overLoad += .01);
  //  iout << iINFO  << "Ending overLoad = " << overLoad << endi;
  //  iout << iINFO
  //   << "After assignment\n" << endi;
  // printLoads();
}















