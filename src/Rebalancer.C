#include <iostream.h>
#include "InfoStream.h"
#include "Rebalancer.h"

Rebalancer::Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
		       processorInfo *processorArray,
		       int nComps, int nPatches, int nPes)
{
  bytesPerAtom = 32;
  strategyName = "dummy";
  overLoad = 1.1;
  computes = computeArray;
  patches =  patchArray;
  processors =  processorArray;
  numComputes = nComps;
  numPatches = nPatches;
  P = nPes;
  int i;
  for (i=0; i<P; i++){
    processors[i].load = processors[i].backgroundLoad;
    processors[i].computeLoad = 0;
    processors[i].patchSet = new Set();
    processors[i].computeSet = new Set();
    processors[i].computesWithBoth = new maxHeap(numComputes);
    processors[i].computesWithOne = new maxHeap(numComputes);
  }
  for (i=0; i<nPatches; i++){
    if (!patches[i].proxiesOn->find(&(processors[patches[i].processor]))) {
      patches[i].proxiesOn->insert(&(processors[patches[i].processor]));
      processors[patches[i].processor].proxies->insert(&(patches[i]));
    }
    processors[ patches[i].processor].patchSet->insert(&patches[i]);
  }		          

  for (i=0; i<numComputes; i++){
    computeArray[i].processor = -1;
  }

  for (i=0; i <P; i++) {
    //    cout << "\n proxies on proc. " << i << " are for patches:";
    //    processorArray[i].proxies->print();
  }
  cout <<"\n";

  //strategy();

}

void Rebalancer::strategy(){
  cout << "Strategy not implemented for the base class.\n";
}

void Rebalancer::makeHeaps()
{
  int i, j;

 pes = new minHeap(P+2);
 for (i=0; i<P; i++)
   pes->insert((InfoRecord *) &(processors[i]));

 computesHeap = new maxHeap(numComputes+2);
 for (i=0; i<numComputes; i++)
   computesHeap->insert( (InfoRecord *) &(computes[i]));

 for (i=0; i<P; i++) {
   processors[i].computesWithBoth = new maxHeap(numComputes);
   processors[i].computesWithOne = new maxHeap(numComputes);
    for (j=0; j<numComputes; j++) {
      int count = 0;
      if (patches[computes[j].patch1].processor = i) count ++;
      if (patches[computes[j].patch2].processor = i) count ++;
      if (count ==2) processors[i].computesWithBoth->
		       insert( (InfoRecord *) &(computes[j]));

      if (count ==1) processors[i].computesWithOne->
		       insert( (InfoRecord *) &(computes[j]));
    }
 }
}

void Rebalancer::assign(computeInfo *c, processorInfo *p)
{
  // cout << "assigning " << c->Id << "with work = "
  //   << c->load << "to processor " << p->processorNum << "\n";
  c->processor = p->Id;
  p->computeSet->insert((InfoRecord *) c);
  p->computeLoad += c->load;
  p->load = p->computeLoad + p->backgroundLoad;
  patchInfo * patch1 =   (patchInfo *) &(patches[c->patch1]);
  patchInfo * patch2 =   (patchInfo *) &(patches[c->patch2]);

  if (!p->proxies->find(patch1))   p->proxies->insert(patch1); 
  if (!patch1->proxiesOn->find(p))
    patch1->proxiesOn->insert(p);

  if (!p->proxies->find(patch2))   p->proxies->insert(patch2); 
  if (!patch2->proxiesOn->find(p)) patch2->proxiesOn->insert(p);
  


}

void  Rebalancer::deAssign(computeInfo *c, processorInfo *p)
{
  c->processor = -1;
  p->computeSet->remove(c);
  p->computeLoad -= c->load;
  p->load = p->computeLoad + p->backgroundLoad;
  // we need to maintain counts of usage if we want to remove patches from
  // the proxies Sets
}

int Rebalancer::refine()
{
  int finish = 1;
  maxHeap *heavyProcessors = new maxHeap(P);

  Set *lightProcessors = new Set();
  int i;
  for (i=0; i<P; i++){
    
    cout << "\n Computes on processor " << i << " ";
    processors[i].computeSet->print();
    cout << "\n";
    if (processors[i].load > overLoad*averageLoad)
      heavyProcessors->insert((InfoRecord *) &(processors[i]));
    else
      if (processors[i].load < averageLoad)
	lightProcessors->insert((InfoRecord *) &(processors[i]));
  }
  int done = 0;
  while (!done){
    double bestSize0, bestSize1, bestSize2;
    computeInfo *bestCompute0, *bestCompute1, *bestCompute2;
    processorInfo *bestP;
    
    processorInfo *donor = (processorInfo *) heavyProcessors->deleteMax();
    if (!donor) break;
    //find the best pair (c,receiver)
    Iterator nextProcessor;
    processorInfo *p = (processorInfo *) 
      lightProcessors->iterator((Iterator *) &nextProcessor);
    bestSize0 = bestSize1 = bestSize2 = 0;
    bestCompute0 = bestCompute1 = bestCompute2 = 0;
    cout << "Finding receiver for processor " << donor->Id << "\n";
    while (p){
      Iterator nextCompute;
      nextCompute.id = 0;
      computeInfo *c = (computeInfo *) donor->computeSet->iterator((Iterator *)&nextCompute);
      cout << "Considering Procsessor : " << p->Id << "\n";
      while (c){
	int n=0;
	n = numAvailable(c,p);
	cout << "Considering Compute : " << c->Id << " with load " << c->load << "\n";
	switch(n){
	case 0: if (( c->load + p->load < overLoad*averageLoad) &&
		    (c->load > bestSize0)) {
	  bestSize0 = c->load;
	  bestCompute0 = c;
	  bestP = p;
	}
	break;
	case 1: if (( c->load + p->load < overLoad*averageLoad) &&
		    (c->load > bestSize1)){
	  bestSize1 = c->load;
	  bestCompute1 = c;
	  bestP = p;
	}
	break;
	case 2: if (( c->load + p->load < overLoad*averageLoad) &&
		    (c->load > bestSize2)){
	  bestSize2 = c->load;
	  bestCompute2 = c;
	  bestP = p;
	}
	break;
	default: cout << "Error. Illegal number of proxies.\n";    
	}
	nextCompute.id++;
	c = (computeInfo *) donor->computeSet->next((Iterator *)&nextCompute);
      }
      p = (processorInfo *) 
	lightProcessors->next((Iterator *) &nextProcessor);
    }
    //we have narrowed the choice to 3 candidates.
    if (bestCompute2){
      deAssign(bestCompute2, donor);      
      assign(bestCompute2, bestP);
    }
    else if (bestCompute1){
      deAssign(bestCompute1, donor);
      assign(bestCompute1, bestP);
    }
    else if (bestCompute0){
      deAssign(bestCompute0, donor);
      assign(bestCompute0, bestP);
    }
    else { 
      cout << "No receiver found" << "\n";
      finish = 0;
      break;
    }
    if (bestP->load > averageLoad)
      lightProcessors->remove(bestP);
    
    if (donor->load > overLoad*averageLoad)
      heavyProcessors->insert((InfoRecord *) donor);
    else if (donor->load < averageLoad)
      lightProcessors->insert((InfoRecord *) donor);
    
  }  
  return finish;
  
}


void Rebalancer::printResults()
{
  cout << "ready to print result \n";
}


void Rebalancer::printLoads()
{
int i, total = 0, numBytes = 0;
double max;

  for (i=0; i<P; i++){
    CPrintf("load on %d is : %f [%f,%f]\n",i,processors[i].load,
	    processors[i].backgroundLoad,processors[i].computeLoad);
    cout << "# Messages received: " << 
      processors[i].proxies->numElements() - processors[i].patchSet->numElements();
    Iterator p;
    int count = 0;
    
    patchInfo *patch = (patchInfo *) processors[i].patchSet->iterator(&p);
    while (patch){
      int myProxies;
      myProxies = patch->proxiesOn->numElements()-1;
      numBytes += myProxies *patch->numAtoms*bytesPerAtom;
      count += myProxies;
      patch = (patchInfo *)processors[i].patchSet->next(&p);
   
    }
    cout << " # Messages sent: " << count << "\n";
    total += count;
  }
  computeAverage();
  max = computeMax();
  CPrintf("Summary: (%s: %d,%d,%d),avg=%lf max=%lf messages = %d [%d bytes]\n",
	  strategyName,P,numPatches,numComputes,
	  averageLoad,max,total,numBytes);
}

void Rebalancer::computeAverage()
{
  int i;
  double total = 0;
  for (i=0; i<numComputes; i++){
    total += computes[i].load;
  }
  for (i=0; i<P; i++){
    total += processors[i].backgroundLoad;
  }
  averageLoad = total/P;
}

double Rebalancer::computeMax()
{
  int i;
  double max = 0;
  for (i=0; i<P; i++){
    if (processors[i].load > max) 
    max = processors[i].load;
  }
  return max;
}


int Rebalancer::numAvailable(computeInfo *c, processorInfo *p)
{
  //return the number of proxy/home patches available on p for c (0,1,2)
  int p1, p2;
  p1 = c->patch1;
  p2 = c->patch2;
  int count = 0;
  if (isAvailableOn((patchInfo *)&(patches[p1]), p))
      count++;
  if (isAvailableOn((patchInfo *)&(patches[p2]), p))
      count++;
  return count;   
  
}

int Rebalancer::isAvailableOn(patchInfo *patch, processorInfo *p)
{
  return  p->proxies->find(patch);
}


