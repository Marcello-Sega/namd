/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include <iostream.h>
#include <iomanip.h>
#include "Rebalancer.h"

Rebalancer::Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
      processorInfo *processorArray, int nComps, int nPatches, int nPes)
{
   bytesPerAtom = 32;
   strategyName = "dummy";
   computes = computeArray;
   patches =  patchArray;
   processors =  processorArray;
   numComputes = nComps;
   numPatches = nPatches;
   P = nPes;
   int i;
   for (i=0; i<P; i++)
   {
      // For testing only...
      // processors[i].backgroundLoad = 0;
      // End of test section
      processors[i].load = processors[i].backgroundLoad;
      processors[i].computeLoad = 0;
      processors[i].patchSet = new Set();
      processors[i].computeSet = new Set();
      processors[i].computesWithBoth = NULL; // new maxHeap(numComputes);
      processors[i].computesWithOne = NULL; // new maxHeap(numComputes);
   }

   InitProxyUsage();

   for (i=0; i<nPatches; i++)
   {
      if (!patches[i].proxiesOn->find(&(processors[patches[i].processor]))) 
      {
         patches[i].proxiesOn->insert(&(processors[patches[i].processor]));
         processors[patches[i].processor].proxies->insert(&(patches[i]));
      }
      processors[patches[i].processor].patchSet->insert(&patches[i]);
   }		          

   for (i=0; i<numComputes; i++)
      computeArray[i].processor = -1;

   for (i=0; i < numComputes; i++)
      processors[computes[i].oldProcessor].computeLoad += computes[i].load;

   // Added 4-29-98: Temporarily adds the compute load to the background
   // load so that the correct value for the total load can be displayed.
   float *temploads = new float[P];
   for(i=0; i<P; i++)
   {
      temploads[i] = processors[i].load;
      processors[i].load += processors[i].computeLoad;
   }

   // iout << iINFO << "Initial load" << "\n";
   // printLoads();

   for(i=0;i<P; i++)
   {
      processors[i].load = temploads[i];
      processors[i].computeLoad = 0;
   }
   
   delete [] temploads;

   // int count1=0, count2=0;
   // for (i=0; i<nPatches; i++)
   // {
   //    if (patches[i].proxiesOn->numElements() <= 1)
   //    count1++;
   //    else count2++;
   // }		          
   // iout << iINFO << "Count1 = " << count1
   //      << "Count2 = " << count2
   //      << "\n" << endl;
   // 
   // for (i=0; i <P; i++) 
   // {
   //    iout << iINFO << "\n proxies on proc. " << i << " are for patches:";
   //    processorArray[i].proxies->print();
   // }
   // 
   // iout << iINFO <<"\n" << endi;
   // strategy();
}

Rebalancer::~Rebalancer()
{
   for(int i=0; i<P; i++)
      delete [] processors[i].proxyUsage;
}

// Added 4-29-98: array proxyUsage on each processor keeps track of 
// how many computes are accessing each proxy on the processor.  If
// no computes are accessing it, the proxy can be removed in DeAssign
void Rebalancer::InitProxyUsage()
{
   for(int i=0; i<P; i++)
   {
      processors[i].proxyUsage = new int[numPatches];
      for(int j=0; j<numPatches; j++)
      {
         processors[i].proxyUsage[j] = 0;
      }

      Iterator nextCompute;
      nextCompute.id = 0;

      computeInfo *c = (computeInfo *)
         processors[i].computeSet->iterator((Iterator *)&nextCompute);

      while(c)
      {
         /* int n1 = */ processors[i].proxyUsage[c->patch1]++;
         /* int n2 = */ processors[i].proxyUsage[c->patch2]++;

         // iout << iINFO  
         // << "Assigning compute " << c->Id << " with work = " << c->load 
         // << " to processor " << processors[i].Id << "\n"
         // << "\tproxyUsage[" << c->patch1 << "]: " << n1 << " --> " << n1+1 << "\n";
         // << "\tproxyUsage[" << c->patch2 << "]: " << n2 << " --> " << n2+1 << "\n";
         // << endl;

         nextCompute.id++;
         c = (computeInfo *) processors[i].computeSet->next((Iterator *)&nextCompute);
      }
   }
}


void Rebalancer::strategy()
{
   iout << iINFO << "Strategy not implemented for the base class.\n" << "\n";
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

   for (i=0; i<P; i++) 
   {
      processors[i].computesWithBoth = NULL; // new maxHeap(numComputes+2);
      processors[i].computesWithOne = NULL; // new maxHeap(numComputes+2);
      for (j=0; j<numComputes; j++) 
      {
         int count = 0;
         if ( (patches[computes[j].patch1].processor = i) ) count ++;
         if ( (patches[computes[j].patch2].processor = i) ) count ++;

         // if (count ==2) processors[i].computesWithBoth->
         //    insert( (InfoRecord *) &(computes[j]));
         // if (count ==1) processors[i].computesWithOne->
         //    insert( (InfoRecord *) &(computes[j]));
      }
   }

   // for (int ii=0; ii<numPatches; ii++)
   // { iout << iINFO << "(3:" << patches[ii].Id << "," << patches[ii].processor <<"]" ;}
   //
}

void Rebalancer::assign(computeInfo *c, int processor)
{
   assign(c, &(processors[processor]));
}

void Rebalancer::assign(computeInfo *c, processorInfo *p)
{
   c->processor = p->Id;
   p->computeSet->insert((InfoRecord *) c);
   p->computeLoad += c->load;
   p->load = p->computeLoad + p->backgroundLoad;
   patchInfo* patch1 = (patchInfo *) &(patches[c->patch1]);
   patchInfo* patch2 = (patchInfo *) &(patches[c->patch2]);

   if (!p->proxies->find(patch1))   p->proxies->insert(patch1); 
   if (!patch1->proxiesOn->find(p)) patch1->proxiesOn->insert(p);

   if (!p->proxies->find(patch2))   p->proxies->insert(patch2); 
   if (!patch2->proxiesOn->find(p)) patch2->proxiesOn->insert(p);
   
   // 4-29-98: Added the following code to keep track of how many proxies
   // on each processor are being used by a compute on that processor
   /* int n1 = */ p->proxyUsage[c->patch1]++;
   /* int n2 = */ p->proxyUsage[c->patch2]++;

   // iout << iINFO  
   // << "Assigning compute " << c->Id << " with work = " << c->load 
   // << " to processor " << p->Id << "\n"
   // << "\tproxyUsage[" << c->patch1 << "]: " << n1 << " --> " << n1+1 << "\n"
   // << "\tproxyUsage[" << c->patch2 << "]: " << n2 << " --> " << n2+1 << "\n"
   // << endl;
}

void  Rebalancer::deAssign(computeInfo *c, processorInfo *p)
{
   c->processor = -1;
   p->computeSet->remove(c);
   p->computeLoad -= c->load;
   p->load = p->computeLoad + p->backgroundLoad;

   // 4-29-98: Added the following code to keep track of how many proxies 
   // on each processor are being used by a compute on that processor.
   // If no computes are using the proxy, it should be removed if it is not
   // on the processor that its patch is on.
   /* int n1 = */ p->proxyUsage[c->patch1]--;
   /* int n2 = */ p->proxyUsage[c->patch2]--;

   // iout << iINFO
   // << "De-assigning compute " << c->Id << " from processor " << p->Id << "\n"
   // << "\tproxyUsage[" << c->patch1 << "]: " << n1 << " --> " << n1-1 << "\n"
   // << "\tproxyUsage[" << c->patch2 << "]: " << n2 << " --> " << n2-1 << "\n"
   // << endl;

   if(p->proxyUsage[c->patch1] <= 0 && p->Id != patches[c->patch1].processor)
   {
      // iout << iINFO 
      // << "REMOVING PROXY " << c->patch1 << " FROM PROCESSOR " << p->Id 
      // << endl << endl;

      patchInfo* patch1 = (patchInfo *) &(patches[c->patch1]);
      p->proxies->remove(patch1);
      patch1->proxiesOn->remove(p);
   }
   if(p->proxyUsage[c->patch2] <= 0 && p->Id != patches[c->patch2].processor)
   {
      // iout << iINFO
      // << "REMOVING PROXY " << c->patch1 << " FROM PROCESSOR " << p->Id 
      // << endl << endl;

      patchInfo* patch2 = (patchInfo *) &(patches[c->patch2]);
      p->proxies->remove(patch2);
      patch2->proxiesOn->remove(p);
   }
}

int Rebalancer::oldrefine()
{
   int finish = 1;
   maxHeap *heavyProcessors = new maxHeap(P);

   Set *lightProcessors = new Set();
   int i;
   int overloaded = 0;
   int underloaded = 0;
   double thresholdLoad = overLoad * averageLoad;
   for (i=0; i<P; i++) {
      // iout << iINFO << "\n Computes on processor " << i << " ";
      // processors[i].computeSet->print();
      // iout << iINFO << "\n" << endi;
     if (processors[i].load >= thresholdLoad ) {
       heavyProcessors->insert((InfoRecord *) &(processors[i]));
       overloaded++;
     } else {
       lightProcessors->insert((InfoRecord *) &(processors[i]));
       underloaded++;
     }
   }
   iout << iINFO << overloaded << " overloaded and " 
	<< underloaded << " underloaded processors\n" << endi;

   int done = 0;
   while (!done)
   {
      computeInfo *bestCompute0, *bestCompute1, *bestCompute2;
      processorInfo *bestP0,*bestP1,*bestP2;
    
      processorInfo *donor = (processorInfo *) heavyProcessors->deleteMax();
      if (!donor) break;

      //find the best pair (c,receiver)
      //      iout << iINFO << "Finding receiver for processor " << donor->Id 
      //	   << "\n" << endi;
      selectComputeCandidates(lightProcessors, donor, thresholdLoad,
			      &bestCompute2, &bestP2,
			      &bestCompute1, &bestP1,
			      &bestCompute0, &bestP0);

      //we have narrowed the choice to 3 candidates.
      processorInfo* bestP;

      if (bestCompute2) {
         deAssign(bestCompute2, donor);      
         assign(bestCompute2, bestP2);
         bestP = bestP2;
      } else if (bestCompute1) {
         deAssign(bestCompute1, donor);
         assign(bestCompute1, bestP1);
         bestP = bestP1;
      } else if (bestCompute0) {
         deAssign(bestCompute0, donor);
         assign(bestCompute0, bestP0);
         bestP = bestP0;
      } else {
         // iout << iINFO << "Refine: No receiver found" << "\n" << endl;
         finish = 0;
         break;
      }

      if (bestP->load > thresholdLoad) {
         lightProcessors->remove(bestP);
	 heavyProcessors->insert((InfoRecord*) bestP);
      }
    
      if (donor->load > thresholdLoad)
         heavyProcessors->insert((InfoRecord *) donor);
      else lightProcessors->insert((InfoRecord *) donor);
   }
  
#if 0
   // After refining, compute min, max and avg processor load
   double total = processors[0].load;
   double min = processors[0].load;
   int min_proc = 0;
   double max = processors[0].load;
   int max_proc = 0;
   for (i=1; i<P; i++) {
     total += processors[i].load;
     if (processors[i].load < min) {
       min = processors[i].load;
       min_proc = i;
     }
     if (processors[i].load > max) {
       max = processors[i].load;
       max_proc = i;
     }
   }
   iout << iINFO << "Refinement at overLoad=" << overLoad << "\n";
   iout << iINFO << "  min = " << min << " processor " << min_proc << "\n";
   iout << iINFO << "  max = " << max << " processor " << max_proc << "\n";
   iout << iINFO << "  total = " << total << " average = " << total/P << "\n"
	<< endi;
   
   if (!finish) {
     iout << iINFO << "Refine: No solution found for overLoad = " 
	  << overLoad << "\n" << endi;
   }
#endif

   return finish;
}

void 
Rebalancer::selectComputeCandidates(Set* lightProcessors,
				    processorInfo* donor,
				    double thresholdLoad, 
				    computeInfo** bestCompute2,
				    processorInfo** bestP2,
				    computeInfo** bestCompute1,
				    processorInfo** bestP1,
				    computeInfo** bestCompute0,
				    processorInfo** bestP0)
{
  Iterator nextProcessor;
  processorInfo *p = (processorInfo *)lightProcessors->
    iterator((Iterator *) &nextProcessor);
  double bestSize0=0;
  double bestSize1=0;
  double bestSize2=0;
  *bestP0 = *bestP1 = *bestP2 = 0;
  *bestCompute0 = *bestCompute1 = *bestCompute2 = 0;


  //  iout << iINFO << "Starting with processor " << (int)p << endi;
  while (p) {
    Iterator nextCompute;
    nextCompute.id = 0;
    computeInfo *c = (computeInfo *) donor->computeSet->
      iterator((Iterator *)&nextCompute);
    //    iout << iINFO << "Considering Procsessor : " << p->Id << "\n" << endi;
    while (c) {
      if ( c->load + p->load < donor->load - p->load)  {
	int n= numAvailable(c,p);
	//	iout << iINFO << "Considering Compute : " << c->Id << " with load " 
	//	     << c->load << "\n" << endi;
	switch(n) {
	case 0: 
	  if(c->load > bestSize0 
	     && (!(*bestP0) || p->load<(*bestP0)->load)) {
	    bestSize0 = c->load;
	    *bestCompute0 = c;
	    *bestP0 = p;
	  }
	  break;
	case 1: 
	  if(c->load > bestSize1 
	     && (!(*bestP1) || p->load< (*bestP1)->load)) {
	    bestSize1 = c->load;
	    *bestCompute1 = c;
	    *bestP1 = p;
	  }
	  break;
	case 2: 
	  if(c->load > bestSize2 
	     && (!(*bestP2) || p->load<(*bestP2)->load)) {
	    bestSize2 = c->load;
	    *bestCompute2 = c;
	    *bestP2 = p;
	  }
	  break;
	default:
	  iout << iINFO <<  "Error. Illegal number of proxies.\n" << "\n";    
	}
      }
      nextCompute.id++;
      c = (computeInfo *) donor->computeSet->next((Iterator *)&nextCompute);
    }
    p = (processorInfo *) lightProcessors->next((Iterator *)&nextProcessor);
  }
}

int Rebalancer::refine()
{
   int finish = 1;
   maxHeap *heavyProcessors = new maxHeap(P);

   Set *lightProcessors = new Set();
   int i;
   double thresholdLoad = overLoad * averageLoad;
   for (i=0; i<P; i++)
   {
      // iout << iINFO << "\n Computes on processor " << i << " ";
      // processors[i].computeSet->print();
      // iout << iINFO << "\n" << endi;
      if (processors[i].load > thresholdLoad)
         heavyProcessors->insert((InfoRecord *) &(processors[i]));
      else lightProcessors->insert((InfoRecord *) &(processors[i]));
   }
   int done = 0;
   while (!done)
   {
      double bestSize0, bestSize1, bestSize2;
      computeInfo *bestCompute0, *bestCompute1, *bestCompute2;
      processorInfo *bestP,*bestP0,*bestP1,*bestP2;
    
      processorInfo *donor = (processorInfo *) heavyProcessors->deleteMax();
      /* Keep selecting new donors, until we find one with some compute to
       * migrate
       */
      computeInfo* c=0;
      while (donor && !c) {
        Iterator nextCompute;
        nextCompute.id = 0;
        c = (computeInfo *) donor->
            computeSet->iterator((Iterator *)&nextCompute);
        if (!c) {
          iout << iINFO << "Ignoring donor " << donor->Id
               << " because no computes\n" << endi;
	  donor = (processorInfo*)heavyProcessors->deleteMax();
        }
      };
  
      if (!donor) break;  // No donors found at all! Give up 

      //find the best pair (c,receiver)
      Iterator nextProcessor;
      processorInfo *p = (processorInfo *) 
      lightProcessors->iterator((Iterator *) &nextProcessor);
      bestSize0 = bestSize1 = bestSize2 = 0;
      bestP0 = bestP1 = bestP2 = 0;
      bestCompute0 = bestCompute1 = bestCompute2 = 0;

      // iout << iINFO << "Finding receiver for processor " << donor->Id << "\n" << endi;
      while (p)
      {
         Iterator nextCompute;
         nextCompute.id = 0;
         computeInfo *c = (computeInfo *) 
            donor->computeSet->iterator((Iterator *)&nextCompute);
         // iout << iINFO << "Considering Procsessor : " << p->Id << "\n" << endi;
         while (c)
         {
            if ( c->load + p->load < thresholdLoad) 
            {
               int n= numAvailable(c,p);
               // iout << iINFO << "Considering Compute : " << c->Id << " with load " 
               //      << c->load << "\n" << endi;
               switch(n)
               {
                  case 0: 
                     if(c->load > bestSize0 && (!bestP0 || p->load<bestP0->load)) 
                     {
                        bestSize0 = c->load;
                        bestCompute0 = c;
                        bestP0 = p;
                     }
                     break;
                  case 1: 
                     if(c->load > bestSize1 && (!bestP1 || p->load<bestP1->load))
                     {
                        bestSize1 = c->load;
                        bestCompute1 = c;
                        bestP1 = p;
                     }
                     break;
                  case 2: 
                     if(c->load > bestSize2 && (!bestP2 || p->load<bestP2->load))
                     {
                        bestSize2 = c->load;
                        bestCompute2 = c;
                        bestP2 = p;
                     }
                     break;
                  default:
                     iout << iINFO <<  "Error. Illegal number of proxies.\n" << "\n";    
               }
            }
            nextCompute.id++;
            c = (computeInfo *) donor->computeSet->next((Iterator *)&nextCompute);
         }
         p = (processorInfo *) 
         lightProcessors->next((Iterator *) &nextProcessor);
      }

      //we have narrowed the choice to 3 candidates.
      if (bestCompute2)
      {
         deAssign(bestCompute2, donor);      
         assign(bestCompute2, bestP2);
         bestP = bestP2;
      }
      else if (bestCompute1)
      {
         deAssign(bestCompute1, donor);
         assign(bestCompute1, bestP1);
         bestP = bestP1;
      }
      else if (bestCompute0)
      {
         deAssign(bestCompute0, donor);
         assign(bestCompute0, bestP0);
         bestP = bestP0;
      }
      else 
      {
         // iout << iINFO << "Refine: No receiver found" << "\n" << endl;
         finish = 0;
         break;
      }

      if (bestP->load > averageLoad)
         lightProcessors->remove(bestP);
    
      if (donor->load > thresholdLoad)
         heavyProcessors->insert((InfoRecord *) donor);
      else lightProcessors->insert((InfoRecord *) donor);
   }  
#if 1
   // After refining, compute min, max and avg processor load
   double total = processors[0].load;
   double min = processors[0].load;
   int min_proc = 0;
   double max = processors[0].load;
   int max_proc = 0;
   for (i=1; i<P; i++) {
     total += processors[i].load;
     if (processors[i].load < min) {
       min = processors[i].load;
       min_proc = i;
     }
     if (processors[i].load > max) {
       max = processors[i].load;
       max_proc = i;
     }
   }
   iout << iINFO << "Refinement at overLoad=" << overLoad << "\n";
   iout << iINFO << "  min = " << min << " processor " << min_proc << "\n";
   iout << iINFO << "  max = " << max << " processor " << max_proc << "\n";
   iout << iINFO << "  total = " << total << " average = " << total/P << "\n"
	<< endi;
   
   if (!finish) {
     iout << iINFO << "Refine: No solution found for overLoad = " 
	  << overLoad << "\n" << endi;
   }
#endif

   delete heavyProcessors;
   delete lightProcessors;

   return finish;
}


void Rebalancer::printResults()
{
  iout << iINFO << "ready to print result \n" << "\n";
}


void Rebalancer::printLoads()
{
#if 0  // Something evil in these print statements.  -JCP

   int i, total = 0, numBytes = 0;
   double max;

   iout << iINFO << "\n" << iINFO;
   for(i=0; i<3; i++) iout << "     TOTAL  BACKGRD COMPUTE | ";
   iout << "\n" << iINFO;
   for(i=0; i<3; i++) iout << " P#  LOAD    LOAD    LOAD   | ";
   iout << "\n" << iINFO;
   for(i=0; i<3; i++) iout << "--- ------- ------- ------- | ";
   iout << "\n" << endi;

   iout.setf(ios::right | ios::fixed);
   iout.precision(3);
   for (i=0; i<P; i++)
   {
      if (i == 0 ) iout << iINFO;
      if (i != 0 && i % 3 == 0) iout << "\n" << endi << iINFO;
      iout << setw(3) << i << " "
           << setw(7) << processors[i].load << " "
           << setw(7) << processors[i].backgroundLoad << " "
           << setw(7) << processors[i].computeLoad << " | ";

      // iout << iINFO << "# Messages received: "
      //	     << processors[i].proxies->numElements() - 
      //         processors[i].patchSet->numElements() 
      //	     << endl << endi;
      // iout << iINFO << "load on "<< i << " is :" << processors[i].load 
      //      << "[ " << processors[i].backgroundLoad << "," 
      //	     << processors[i].computeLoad << "]. ";
      // iout << iINFO << "# Messages received: " 
      //      << processors[i].proxies->numElements() - 
      //         processors[i].patchSet->numElements();

      Iterator p;
      int count = 0;
    
      patchInfo *patch = (patchInfo *) processors[i].patchSet->iterator(&p);
      while (patch)
      {
         int myProxies;
         myProxies = patch->proxiesOn->numElements()-1;
         numBytes += myProxies *patch->numAtoms*bytesPerAtom;
         count += myProxies;
         patch = (patchInfo *)processors[i].patchSet->next(&p);
      }
      total += count;

      // iout << iINFO << " # Messages sent: " << count << "\n" << endi;
   }

   iout << "\n" << endi;

   computeAverage();
   max = computeMax();

   iout << iINFO << "\n" << endi;
   iout << iINFO << "------------------------------------------------------------\n" << endi; 
   iout << iINFO << "          LOAD SUMMARY FOR STRATEGY \"" << strategyName << "\"\n\n" << endi;
   iout << iINFO << "Processors = " << setw(5) << P << "\t"
        << "  Overload = " << setw(7) << overLoad << "\n";
   iout << iINFO << "Patches    = " << setw(5) << numPatches << "\t"
        << "  Avg load = " << setw(7) << averageLoad << "\n";
   iout << iINFO << "Computes   = " << setw(5) << numComputes << "\t"
        << "  Max load = " << setw(7) << max << "\n";
   iout << iINFO << "# messages = " << setw(5) << total << "\t"
        << "  Msg size = " << numBytes << " bytes" << "\n" << "\n";
  iout << iINFO <<"============================================================\n"
       << "\n" << endi;
   iout.unsetf(ios::right);

#endif

}

double Rebalancer::computeAverage()
{
   int i;
   double total = 0;
   for (i=0; i<numComputes; i++)
      total += computes[i].load;

   for (i=0; i<P; i++)
      total += processors[i].backgroundLoad;
  
   averageLoad = total/P;
   return averageLoad;
}

double Rebalancer::computeMax()
{
   int i;
   double max = processors[0].load;
   for (i=1; i<P; i++)
   {
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
