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
   pes = NULL;
   computesHeap = NULL;
   int i;
   for (i=0; i<P; i++)
   {
      // For testing only...
      // processors[i].backgroundLoad = 0;
      // End of test section
      processors[i].load = processors[i].backgroundLoad;
      processors[i].computeLoad = 0;
   }

   for (i=0; i<nPatches; i++)
   {
      if (!patches[i].proxiesOn.find(&(processors[patches[i].processor]))) 
      {
         patches[i].proxiesOn.insert(&(processors[patches[i].processor]));
         processors[patches[i].processor].proxies.insert(&(patches[i]));
      }
      processors[patches[i].processor].patchSet.insert(&patches[i]);
   }		          

   InitProxyUsage();

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
   printLoads();

   for(i=0;i<P; i++)
   {
      processors[i].load = temploads[i];
      processors[i].computeLoad = 0;
   }
   
   delete [] temploads;

   // int count1=0, count2=0;
   // for (i=0; i<nPatches; i++)
   // {
   //    if (patches[i].proxiesOn.numElements() <= 1)
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

   // for (i=0; i<nPatches; i++)
   // {
   //    iout << "patch " << i << " on processor " << patches[i].processor << "\n" << endi;
   // }
}

Rebalancer::~Rebalancer()
{
   for(int i=0; i<P; i++)
      delete [] processors[i].proxyUsage;
   delete pes;
   delete computesHeap;
}

// Added 4-29-98: array proxyUsage on each processor keeps track of 
// how many computes are accessing each proxy on the processor.  If
// no computes are accessing it, the proxy can be removed in DeAssign
void Rebalancer::InitProxyUsage()
{
   int i;

   for(i=0; i<P; i++) {
      processors[i].proxyUsage = new int[numPatches];
      for(int j=0; j<numPatches; j++)
      {
         processors[i].proxyUsage[j] = 0;
      }

      Iterator nextCompute;
      nextCompute.id = 0;

      computeInfo *c = (computeInfo *)
         processors[i].computeSet.iterator((Iterator *)&nextCompute);

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
         c = (computeInfo *) processors[i].computeSet.next((Iterator *)&nextCompute);
      }
   }

  for (i=0; i<numPatches; i++)
  {
      Iterator nextProc;
      processorInfo *p = (processorInfo *)patches[i].proxiesOn.iterator((Iterator *)&nextProc);
      while (p) {
          p->proxyUsage[i] += 1;
          p = (processorInfo *)patches[i].proxiesOn.next((Iterator*)&nextProc);
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

   delete pes;
   pes = new minHeap(P+2);
   for (i=0; i<P; i++)
      pes->insert((InfoRecord *) &(processors[i]));

   delete computesHeap;
   computesHeap = new maxHeap(numComputes+2);
   for (i=0; i<numComputes; i++)
      computesHeap->insert( (InfoRecord *) &(computes[i]));

   for (i=0; i<P; i++) 
   {
      for (j=0; j<numComputes; j++) 
      {
         int count = 0;
         if ( (patches[computes[j].patch1].processor == i) ) count ++;
         if ( (patches[computes[j].patch2].processor == i) ) count ++;
      }
   }
}

void Rebalancer::assign(computeInfo *c, int processor)
{
   assign(c, &(processors[processor]));
}

void Rebalancer::assign(computeInfo *c, processorInfo *p)
{
   c->processor = p->Id;
   p->computeSet.insert((InfoRecord *) c);
   p->computeLoad += c->load;
   p->load = p->computeLoad + p->backgroundLoad;
   patchInfo* patch1 = (patchInfo *) &(patches[c->patch1]);
   patchInfo* patch2 = (patchInfo *) &(patches[c->patch2]);

   if (!p->proxies.find(patch1))   p->proxies.insert(patch1); 
   if (!patch1->proxiesOn.find(p)) patch1->proxiesOn.insert(p);

   if (!p->proxies.find(patch2))   p->proxies.insert(patch2); 
   if (!patch2->proxiesOn.find(p)) patch2->proxiesOn.insert(p);
   
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
   p->computeSet.remove(c);
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
      p->proxies.remove(patch1);
      patch1->proxiesOn.remove(p);
   }
   if(p->proxyUsage[c->patch2] <= 0 && p->Id != patches[c->patch2].processor)
   {
      // iout << iINFO
      // << "REMOVING PROXY " << c->patch1 << " FROM PROCESSOR " << p->Id 
      // << endl << endl;

      patchInfo* patch2 = (patchInfo *) &(patches[c->patch2]);
      p->proxies.remove(patch2);
      patch2->proxiesOn.remove(p);
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
      processorInfo* goodP[3][3];  // goodP[# of real patches][# of proxies]
      computeInfo* goodCompute[3][3];
      double goodSize[3][3];
      processorInfo* bestP;
    
      processorInfo *donor = (processorInfo *) heavyProcessors->deleteMax();
      /* Keep selecting new donors, until we find one with some compute to
       * migrate
       */
/*
      computeInfo* c=0;
      while (donor && !c) {
        Iterator nextCompute;
        nextCompute.id = 0;
        c = (computeInfo *) donor->
            computeSet.iterator((Iterator *)&nextCompute);
        if (!c) {
          iout << iINFO << "Ignoring donor " << donor->Id
               << " because no computes\n" << endi;
	  donor = (processorInfo*)heavyProcessors->deleteMax();
        }
      };
*/
      while (donor) {
	if (donor->computeSet.numElements()) break;
        iout << iINFO << "Ignoring donor " << donor->Id
              << " because no computes\n" << endi;
	 donor = (processorInfo*)heavyProcessors->deleteMax();
      }
  
      if (!donor) break;  // No donors found at all! Give up 

      //find the best pair (c,receiver)
      Iterator nextProcessor;
      processorInfo *p = (processorInfo *) 
      lightProcessors->iterator((Iterator *) &nextProcessor);

      int i,j;
      for(i=0; i < 3; i++)
	for(j=0; j<3; j++) {
	  goodP[i][j] = 0;
	  goodCompute[i][j] = 0;
	  goodSize[i][j] = 0.;
	}

      // iout << iINFO << "Finding receiver for processor " << donor->Id << "\n" << endi;
      while (p)
      {
         Iterator nextCompute;
         nextCompute.id = 0;
         computeInfo *c = (computeInfo *) 
            donor->computeSet.iterator((Iterator *)&nextCompute);
         // iout << iINFO << "Considering Procsessor : " << p->Id << "\n" << endi;
         while (c)
         {
            if ( c->load + p->load < thresholdLoad) 
            {
               int nPatches = numPatchesAvail(c,p);
	       int nProxies = numProxiesAvail(c,p);
	       
	       if (nPatches < 0 || nPatches > 2)
		 iout << iERROR << "Too many patches: " << nPatches 
		      << "\n" << endi;
	       if (nProxies < 0 || nProxies > 2)
		 iout << iERROR << "Too many proxies: " << nProxies 
		      << "\n" << endi;
	       if (nProxies + nPatches > 2)
		 iout << iERROR << "Too many patches (" << nPatches
		      << ") + proxies (" << nProxies << ")\n" << endi;

	       if ((c->load > goodSize[nPatches][nProxies]) 
		   && (!goodP[nPatches][nProxies] 
		       || p->load < goodP[nPatches][nProxies]->load) ) {
		 goodSize[nPatches][nProxies] = c->load;
		 goodCompute[nPatches][nProxies] = c;
		 goodP[nPatches][nProxies] = p;
	       }
	    }
            nextCompute.id++;
            c = (computeInfo *) donor->computeSet.
	      next((Iterator *)&nextCompute);
         }
         p = (processorInfo *) 
	   lightProcessors->next((Iterator *) &nextProcessor);
      }

      //we have narrowed the choice to 6 candidates.
      if (goodCompute[2][0]) {
         deAssign(goodCompute[2][0], donor);      
         assign(goodCompute[2][0], goodP[2][0]);
         bestP = goodP[2][0];
      } else if (goodCompute[1][1]) {
         deAssign(goodCompute[1][1], donor);      
         assign(goodCompute[1][1], goodP[1][1]);
         bestP = goodP[1][1];
      } else if (goodCompute[0][2]) {
         deAssign(goodCompute[0][2], donor);      
         assign(goodCompute[0][2], goodP[0][2]);
         bestP = goodP[0][2];
      } else if (goodCompute[1][0]) {
         deAssign(goodCompute[1][0], donor);      
         assign(goodCompute[1][0], goodP[1][0]);
         bestP = goodP[1][0];
      } else if (goodCompute[0][1]) {
         deAssign(goodCompute[0][1], donor);      
         assign(goodCompute[0][1], goodP[0][1]);
         bestP = goodP[0][1];
      } else if (goodCompute[0][0]) {
         deAssign(goodCompute[0][0], donor);      
         assign(goodCompute[0][0], goodP[0][0]);
         bestP = goodP[0][0];
      } else {
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

// this binary search refinement procedure assume you already assigned computes
// to their processors before calling this!!
void Rebalancer::multirefine()
{
  // The New refinement procedure.  This is identical to the code in
  // RefineOnly.C, and probably should be merged with that code to form
  // a binary-search function

  double avg = computeAverage();
  double max = computeMax();

  const double overloadStep = 0.01;
  const double overloadStart = 1.02;
  double dCurOverload = max / avg;

  int minOverload = 0;
  int maxOverload = (int)((dCurOverload - overloadStart)/overloadStep + 1);
  double dMinOverload = minOverload * overloadStep + overloadStart;
  double dMaxOverload = maxOverload * overloadStep + overloadStart;

  iout << iINFO
       << "Balancing from " << minOverload << " = " << dMinOverload 
       << " to " << maxOverload << "=" << dMaxOverload 
       << " dCurOverload=" << dCurOverload << " max=" << max << " avg=" << avg
       << "\n" << endi;

  int curOverload;
  int refineDone = 0;

  overLoad = dMinOverload;
  if (refine())
    refineDone = 1;
  else {
    overLoad = dMaxOverload;
    if (!refine()) {
      iout << iINFO << "ERROR: Could not refine at max overload\n" << endi;
      refineDone = 1;
    }
  }

  // Scan up, until we find a refine that works
  while (!refineDone) {
    if (maxOverload - minOverload <= 1)
      refineDone = 1;
    else {
      curOverload = (maxOverload + minOverload ) / 2;

      overLoad = curOverload * overloadStep + overloadStart;
      iout << iINFO << "Testing curOverload " << curOverload 
	   << "=" << overLoad << " [min,max]=" 
	   << minOverload << ", " << maxOverload
	   << "\n" << endi;
      if (refine())
	maxOverload = curOverload;
      else
	minOverload = curOverload;
    }
  }

}

void Rebalancer::printResults()
{
  iout << iINFO << "ready to print result \n" << "\n";
}


void Rebalancer::printLoads()
{
#if 1  // Something evil in these print statements.  -JCP

   int i, total = 0, numBytes = 0;
   double max;

#if 0
   iout << iINFO << "\n" << iINFO;
   for(i=0; i<3; i++) iout << "     TOTAL  BACKGRD COMPUTE | ";
   iout << "\n" << iINFO;
   for(i=0; i<3; i++) iout << " P#  LOAD    LOAD    LOAD   | ";
   iout << "\n" << iINFO;
   for(i=0; i<3; i++) iout << "--- ------- ------- ------- | ";
   iout << "\n" << endi;
#endif

   iout.setf(ios::right | ios::fixed);
   iout.precision(3);
   for (i=0; i<P; i++)
   {
#if 0
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
#endif

      Iterator p;
      int count = 0;
    
      patchInfo *patch = (patchInfo *) processors[i].patchSet.iterator(&p);
      while (patch)
      {
         int myProxies;
         myProxies = patch->proxiesOn.numElements()-1;
         numBytes += myProxies *patch->numAtoms*bytesPerAtom;
         count += myProxies;
         patch = (patchInfo *)processors[i].patchSet.next(&p);
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

void Rebalancer::printSummary()
{
   int i;
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
   iout << iINFO << "SUMMARY" << "\n";
   iout << iINFO << "  min = " << min << " processor " << min_proc << "\n";
   iout << iINFO << "  max = " << max << " processor " << max_proc << "\n";
   iout << iINFO << "  total = " << total << " average = " << total/P << "\n"
	<< endi;
   
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
   // self computes get one patch for free here
   if (p1 == p2 || isAvailableOn((patchInfo *)&(patches[p2]), p))
      count++;
   return count;   
}

int Rebalancer::numProxiesAvail(computeInfo *c, processorInfo *p)
{
   //return the number of proxy/home patches available on p for c (0,1,2)
   int p1, p2;
   p1 = c->patch1;
   p2 = c->patch2;
   int count = 0;
   if (isAvailableOn((patchInfo *)&(patches[p1]), p) 
       && patches[p1].processor != p->Id ) {
      count++;
      //iout << iINFO << "Patch " << patches[p1].Id << " has a proxy on " 
      //     << p->Id << "\n" << endi;
   }
   if (p1 != p2  // self computes get one patch for free so don't allow 2
       && isAvailableOn((patchInfo *)&(patches[p2]), p)
       && patches[p2].processor != p->Id ) {
      count++;
      //iout << iINFO << "Patch " << patches[p2].Id << " has a proxy on " 
      //   << p->Id << "\n" << endi;
   }

   //iout << iINFO << "Returning " << count << " proxies\n" << endi;
   return count;   
}

int Rebalancer::numPatchesAvail(computeInfo *c, processorInfo *p)
{
   //return the number of proxy/home patches available on p for c (0,1,2)
   const int p1 = c->patch1;
   const int p2 = c->patch2;

   int count = 0;

   if (patches[p1].processor == p->Id) {
     count++;
     //iout << iINFO << "Patch " << patches[p1].Id << " is on " 
     //  << patches[p1].processor << "\n" << endi;
   }
   // self computes get one patch for free here
   if (p1 == p2 || patches[p2].processor == p->Id) {
     count++;
     //iout << iINFO << "Patch " << patches[p2].Id << " is on " 
     //  << patches[p2].processor << "\n" << endi;
   }
     
   //iout << iINFO << "Returning " << count << " patches\n" << endi;
   return count;   
}

int Rebalancer::isAvailableOn(patchInfo *patch, processorInfo *p)
{
   return  p->proxies.find(patch);
}
