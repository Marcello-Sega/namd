/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "Rebalancer.h"


Rebalancer::Rebalancer(computeInfo *computeArray, patchInfo *patchArray,
      processorInfo *processorArray, int nComps, int nPatches, int nPes)
{
   bytesPerAtom = 32;
   strategyName = "None";
   computes = computeArray;
   patches =  patchArray;
   processors =  processorArray;
   numComputes = nComps;
   numPatches = nPatches;
   P = nPes;
   pes = NULL;
   computePairHeap = NULL;
   computeSelfHeap = NULL;
   computeBgPairHeap = NULL;
   computeBgSelfHeap = NULL;
   overLoad = 0.;
   numPesAvailable = 0;
   int i;
   for (i=0; i<P; i++)
   {
      // For testing only...
      // processors[i].backgroundLoad = 0;
      // End of test section
      processors[i].load = processors[i].backgroundLoad;
      processors[i].computeLoad = 0;
      if (processors[i].available) {
        numPesAvailable += 1;
      }
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
   delete computePairHeap;
   delete computeSelfHeap;
   delete computeBgPairHeap;
   delete computeBgSelfHeap;
}

// Added 4-29-98: array proxyUsage on each processor keeps track of 
// how many computes are accessing each proxy on the processor.  If
// no computes are accessing it, the proxy can be removed in DeAssign
void Rebalancer::InitProxyUsage()
{
   int i;
   numProxies = 0;

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
      numProxies += ( patches[i].proxiesOn.numElements() - 1 );
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

   delete computePairHeap;
   delete computeSelfHeap;
   delete computeBgPairHeap;
   delete computeBgSelfHeap;

   double bgLoadLimit = 0.5 * averageLoad;
   /*
   iout << iINFO << "Background load limit = " << bgLoadLimit << "\n";
   for (i=0; i<P; i++)
     if ( processors[i].backgroundLoad > bgLoadLimit )
       iout << iINFO << "Processor " << i << " background load = "
            << processors[i].backgroundLoad << "\n";
   iout << endi;
   */

   int numSelfComputes, numPairComputes, numBgSelfComputes, numBgPairComputes;

   while ( 1 ) {
    numSelfComputes = 0;
    numPairComputes = 0;
    numBgSelfComputes = 0;
    numBgPairComputes = 0;
    for (i=0; i<numComputes; i++) {
     int pa1 = computes[i].patch1;
     int pa2 = computes[i].patch2;
     if ( pa1 == pa2 ) {
        if ( processors[patches[pa1].processor].backgroundLoad > bgLoadLimit) {
          ++numBgSelfComputes;
        } else {
          ++numSelfComputes;
        }
     } else {
        if ( processors[patches[pa1].processor].backgroundLoad > bgLoadLimit
          || processors[patches[pa2].processor].backgroundLoad > bgLoadLimit) {
          ++numBgPairComputes;
        } else {
          ++numPairComputes;
        }
     }
    }

    int numBgComputes = numBgPairComputes + numBgSelfComputes;

    if ( numBgComputes ) {
        iout << iINFO << numBgComputes << " of " << numComputes
        << " computes have background load > " << bgLoadLimit << "\n" << endi;
    }

    if ( numBgComputes < 0.3 * numComputes ) break;
    else bgLoadLimit += 0.1 * averageLoad;
   }

   computePairHeap = new maxHeap(numPairComputes+2);
   computeSelfHeap = new maxHeap(numSelfComputes+2);
   computeBgPairHeap = new maxHeap(numBgPairComputes+2);
   computeBgSelfHeap = new maxHeap(numBgSelfComputes+2);

   for (i=0; i<numComputes; i++) {
     int pa1 = computes[i].patch1;
     int pa2 = computes[i].patch2;
     if ( pa1 == pa2 ) {
        if ( processors[patches[pa1].processor].backgroundLoad > bgLoadLimit) {
          computeBgSelfHeap->insert( (InfoRecord *) &(computes[i]));
        } else {
          computeSelfHeap->insert( (InfoRecord *) &(computes[i]));
        }
     } else {
        if ( processors[patches[pa1].processor].backgroundLoad > bgLoadLimit
          || processors[patches[pa2].processor].backgroundLoad > bgLoadLimit) {
          computeBgPairHeap->insert( (InfoRecord *) &(computes[i]));
        } else {
          computePairHeap->insert( (InfoRecord *) &(computes[i]));
        }
     }
   }

/*
   delete computePairHeap;
   delete computeSelfHeap;

   int numSelfComputes = 0;
   for (i=0; i<numComputes; i++)
      if ( computes[i].patch1 == computes[i].patch2 ) ++numSelfComputes;

   computeSelfHeap = new maxHeap(numSelfComputes+2);
   computePairHeap = new maxHeap(numComputes-numSelfComputes+2);

   for (i=0; i<numComputes; i++)
      if ( computes[i].patch1 == computes[i].patch2 )
         computeSelfHeap->insert( (InfoRecord *) &(computes[i]));
      else
         computePairHeap->insert( (InfoRecord *) &(computes[i]));
*/
}

void Rebalancer::assign(computeInfo *c, int processor)
{
   assign(c, &(processors[processor]));
}

void Rebalancer::assign(computeInfo *c, processorInfo *p)
{
#ifdef LDB_VERBOSE
   int nPatches, nProxies, badForComm;
   numAvailable(c,p,&nPatches,&nProxies,&badForComm);
#endif

   c->processor = p->Id;
   p->computeSet.insert((InfoRecord *) c);
   p->computeLoad += c->load;
   p->load = p->computeLoad + p->backgroundLoad;
   patchInfo* patch1 = (patchInfo *) &(patches[c->patch1]);
   patchInfo* patch2 = (patchInfo *) &(patches[c->patch2]);

   if (!p->proxies.find(patch1))   p->proxies.insert(patch1); 
   if (!patch1->proxiesOn.find(p)) {patch1->proxiesOn.insert(p); numProxies++;}

   if (!p->proxies.find(patch2))   p->proxies.insert(patch2); 
   if (!patch2->proxiesOn.find(p)) {patch2->proxiesOn.insert(p); numProxies++;}
   
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

#ifdef LDB_VERBOSE
   iout << "Assign " << c->Id << " patches " << c->patch1 << " " << c->patch2
        << " load " << c->load << " to " << p->Id << " new load "
        << p->load << " background " << p->backgroundLoad
        << " nPatches " << nPatches << " nProxies " << nProxies;
   if ( nPatches + nProxies < 2 ) iout << " addProxy";
   if ( badForComm ) iout << " badForComm";
   iout << "\n" << endi;
#endif
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
      numProxies--;
   }
   if(p->proxyUsage[c->patch2] <= 0 && p->Id != patches[c->patch2].processor)
   {
      // iout << iINFO
      // << "REMOVING PROXY " << c->patch1 << " FROM PROCESSOR " << p->Id 
      // << endl << endl;

      patchInfo* patch2 = (patchInfo *) &(patches[c->patch2]);
      p->proxies.remove(patch2);
      patch2->proxiesOn.remove(p);
      numProxies--;
   }
}

void Rebalancer::refine_togrid(pcgrid &grid, double thresholdLoad,
			processorInfo *p, computeInfo *c) {

  if(p->available == CmiFalse) return;

  if ( c->load + p->load < thresholdLoad) {
    int nPatches, nProxies, badForComm;
    numAvailable(c,p,&nPatches,&nProxies,&badForComm);

    // if ( badForComm ) return;

    pcpair *pair = &grid[nPatches][nProxies][badForComm];

    if (! pair->c) {
      pair->c = c;
      pair->p = p;
    } else {
      double newval = p->load - c->load;
      if ( c->load + p->load < averageLoad ) {
         newval -= averageLoad;
      }
      double oldval = pair->p->load - pair->c->load;
      if ( pair->c->load + pair->p->load < averageLoad ) {
         oldval -= averageLoad;
      }
      if (newval < oldval) {
	pair->c = c;
	pair->p = p;
      }
    }
  }
}

int Rebalancer::refine()
{
   int finish = 1;
   int no_new_proxies = 0;  // set to true if new proxies are futile
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
      // processorInfo *donor = (processorInfo *) heavyProcessors->deleteMax();
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

      processorInfo *donor;
      while (donor = (processorInfo*)heavyProcessors->deleteMax()) {
	if (donor->computeSet.numElements()) break;
        if ( ! no_new_proxies ) {
          /*
          iout << iINFO << "Most-loaded processor " << donor->Id
               << " (" << donor->patchSet.numElements() << " patches, "
               << donor->proxies.numElements() << " proxies)"
               << " has no migratable work.\n" << endi;
          */
          no_new_proxies = 1;  // New proxies would not improve load balance.
        }
      }
  
      if (!donor) break;  // No donors found at all! Give up 

      pcgrid grid;
#define REASSIGN(GRID) if (GRID.c) { \
	   deAssign(GRID.c, donor); \
           assign(GRID.c, GRID.p); \
           bestP = GRID.p; \
        }

      // try for at least one proxy
      {
        Iterator nextCompute;
        nextCompute.id = 0;
        computeInfo *c = (computeInfo *) 
           donor->computeSet.iterator((Iterator *)&nextCompute);
        while (c)
        {
          Iterator nextProc;
          processorInfo *p;

	  p = &processors[patches[c->patch1].processor];
	  refine_togrid(grid, thresholdLoad, p, c);

	  if (c->patch1 != c->patch2)
	  {
	  p = &processors[patches[c->patch2].processor];
	  refine_togrid(grid, thresholdLoad, p, c);
	  }

	  p = (processorInfo *)patches[c->patch1].
				proxiesOn.iterator((Iterator *)&nextProc);
          while (p) {
	    refine_togrid(grid, thresholdLoad, p, c);
            p = (processorInfo *)patches[c->patch1].
				proxiesOn.next((Iterator*)&nextProc);
          }

	  if (c->patch1 != c->patch2) 
	  {
	  p = (processorInfo *)patches[c->patch2].
				proxiesOn.iterator((Iterator *)&nextProc);
          while (p) {
	    refine_togrid(grid, thresholdLoad, p, c);
            p = (processorInfo *)patches[c->patch2].
				proxiesOn.next((Iterator*)&nextProc);
          }
	  }

          nextCompute.id++;
          c = (computeInfo *) donor->computeSet.
	    next((Iterator *)&nextCompute);
        }
        processorInfo* bestP = 0;
        // prefer proxies to home patches
	REASSIGN(grid[0][2][0])
	else REASSIGN(grid[1][1][0])
	else REASSIGN(grid[2][0][0])
        else if ( no_new_proxies ) { finish = 0; break; }
	else REASSIGN(grid[0][1][0])
	else REASSIGN(grid[1][0][0])
	else REASSIGN(grid[0][0][0])
	// else REASSIGN(grid[0][1][1])
	// else REASSIGN(grid[1][0][1])
	// else REASSIGN(grid[0][0][1])
        if (bestP) {
	  if (bestP->load > averageLoad) lightProcessors->remove(bestP);
	  if (donor->load > thresholdLoad)
		heavyProcessors->insert((InfoRecord *) donor);
	  else lightProcessors->insert((InfoRecord *) donor);
	  continue;
        }
      }

      if ( no_new_proxies ) iout << iINFO
         << "ERROR: Rebalancer::refine() algorithm is broken.\n" << endi;

      // no luck, do it the long way

      //find the best pair (c,receiver)
      Iterator nextProcessor;
      processorInfo *p = (processorInfo *) 
      lightProcessors->iterator((Iterator *) &nextProcessor);

      while (p)
      {
         Iterator nextCompute;
         nextCompute.id = 0;
         computeInfo *c = (computeInfo *) 
            donor->computeSet.iterator((Iterator *)&nextCompute);
         while (c)
         {
#if CMK_VERSION_BLUEGENE
	   BGLTorusManager *tmgr = BGLTorusManager::getObject();
	   if(tmgr->isNeighborOfBoth(p->Id, patches[c->patch1].processor, 
				     patches[c->patch2].processor, 2))
#endif
	     {	     
	       refine_togrid(grid, thresholdLoad, p, c);
	     }
	   nextCompute.id++;
	   c = (computeInfo *) donor->computeSet.
	     next((Iterator *)&nextCompute);
         }
         p = (processorInfo *) 
	   lightProcessors->next((Iterator *) &nextProcessor);
      }

      //we have narrowed the choice to 6 candidates.
      // prefer proxies to home patches
      {
        processorInfo* bestP = 0;
	REASSIGN(grid[0][2][0])
	else REASSIGN(grid[1][1][0])
	else REASSIGN(grid[2][0][0])
	else REASSIGN(grid[0][1][0])
	else REASSIGN(grid[1][0][0])
	else REASSIGN(grid[0][0][0])
	// else REASSIGN(grid[0][1][1])
	// else REASSIGN(grid[1][0][1])
	// else REASSIGN(grid[0][0][1])
	else { finish = 0; break; }
	if (bestP->load > averageLoad) lightProcessors->remove(bestP);
	if (donor->load > thresholdLoad)
		heavyProcessors->insert((InfoRecord *) donor);
	else lightProcessors->insert((InfoRecord *) donor);
      }

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

   delete heavyProcessors;
   delete lightProcessors;

   return finish;
}

// this binary search refinement procedure assume you already assigned computes
// to their processors before calling this!!
void Rebalancer::multirefine(double overload_start)
{
  // The New refinement procedure.  This is identical to the code in
  // RefineOnly.C, and probably should be merged with that code to form
  // a binary-search function

  double avg = computeAverage();
  double max = computeMax();

  int numOverloaded = 0;
  for (int ip=0; ip<P; ip++) {
    if ( processors[ip].backgroundLoad > averageLoad ) ++numOverloaded;
  }
  if ( numOverloaded ) {
    iout << iWARN << numOverloaded
      << " processors are overloaded due to high background load.\n" << endi;
  }

  const double overloadStep = 0.01;
  const double overloadStart = overload_start;       //1.05;
  double dCurOverload = max / avg;
  
  int minOverload = 0;   //Min overload should be 1.05 ?
  int maxOverload = (int)((dCurOverload - overloadStart)/overloadStep + 1);
  double dMinOverload = minOverload * overloadStep + overloadStart;
  double dMaxOverload = maxOverload * overloadStep + overloadStart;

#if 0
  iout << iINFO
       << "Balancing from " << minOverload << " = " << dMinOverload 
       << " to " << maxOverload << "=" << dMaxOverload 
       << " dCurOverload=" << dCurOverload << " max=" << max << " avg=" << avg
       << "\n" << endi;
#endif

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
      /*
      iout << iINFO << "Testing curOverload " << curOverload 
	   << "=" << overLoad << " [min,max]=" 
	   << minOverload << ", " << maxOverload
	   << "\n" << endi;
      */
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

   int i, total = 0, numBytes = 0;
   double max;

#ifdef LDB_VERBOSE
   iout << iINFO << "Proc  Total  Background  Compute\n" << endi;
#endif

// Something evil in these print statements.  -JCP
#if 0
   iout << iINFO << "\n" << iINFO;
   for(i=0; i<3; i++) iout << "     TOTAL  BACKGRD COMPUTE | ";
   iout << "\n" << iINFO;
   for(i=0; i<3; i++) iout << " P#  LOAD    LOAD    LOAD   | ";
   iout << "\n" << iINFO;
   for(i=0; i<3; i++) iout << "--- ------- ------- ------- | ";
   iout << "\n" << endi;
#endif

#if 0
   iout.setf(ios::right | ios::fixed);
   iout.precision(3);
#endif
   int maxproxies = 0;
   int maxpatchproxies = 0;
   for (i=0; i<P; i++)
   {
#ifdef LDB_VERBOSE
   iout << iINFO << i << " "
           << processors[i].load << " "
           << processors[i].backgroundLoad << " "
           << processors[i].computeLoad << "\n" << endi;
#endif
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

      int nproxies = processors[i].proxies.numElements() - 
			processors[i].patchSet.numElements();
      if ( nproxies > maxproxies ) maxproxies = nproxies;

      Iterator p;
      int count = 0;
    
      patchInfo *patch = (patchInfo *) processors[i].patchSet.iterator(&p);
      while (patch)
      {
         int myProxies;
         myProxies = patch->proxiesOn.numElements()-1;
         if ( myProxies > maxpatchproxies ) maxpatchproxies = myProxies;
         numBytes += myProxies *patch->numAtoms*bytesPerAtom;
         count += myProxies;
         patch = (patchInfo *)processors[i].patchSet.next(&p);
      }
      total += count;

      // iout << iINFO << " # Messages sent: " << count << "\n" << endi;
   }

   computeAverage();
   max = computeMax();

#if 0
   // iout << iINFO << "------------------------------------------------------------\n" << endi; 
   iout << iINFO << "          LOAD SUMMARY FOR STRATEGY \"" << strategyName << "\"\n" << endi;
   // iout << iINFO << "Processors = " << setw(5) << P << "\t"
   //      << "  Overload = " ; setw(7); iout << overLoad << "\n";
   iout << iINFO << "Patches    = " << setw(5) << numPatches << "\t"
        << "  Avg load = " ; setw(7); iout << averageLoad << "\n";
   iout << iINFO << "Computes   = " << setw(5) << numComputes << "\t"
        << "  Max load = " ; setw(7); iout << max << "\n";
   iout << iINFO << "Messages   = " << setw(5) << total << "\t"
        << "  Max msgs = " << maxproxies << ", " << maxpatchproxies << "\n";
   // iout << iINFO << "------------------------------------------------------------\n" << endi; 
   iout << endi;
   iout.unsetf(ios::right);
#else
   iout << "LDB:  LOAD: AVG " << averageLoad << " MAX " << max
     << "  MSGS: TOTAL " << total
     << " MAXC " << maxproxies << " MAXP " << maxpatchproxies
     << "  " << strategyName << "\n" << endi;

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
   double total = 0.;
   for (i=0; i<numComputes; i++)
      total += computes[i].load;

   for (i=0; i<P; i++) {
      if (processors[i].available) {
        total += processors[i].backgroundLoad;
      }
   }
  
   if (numPesAvailable == 0) {
     CmiPrintf("Warning: no processors available for load balancing!\n");
     averageLoad = 0.0;
   }
   else 
     averageLoad = total/numPesAvailable;
   return averageLoad;
}

void Rebalancer::adjustBackgroundLoadAndComputeAverage()
{
  // useful for AlgSeven when some loads start out as zero

   if (numPesAvailable == 0) {
     computeAverage();  // because otherwise someone will forget
     return;
   }
  
   int i;
   double bgtotal = 0.;
   for (i=0; i<P; i++) {
      if (processors[i].available) {
        bgtotal += processors[i].backgroundLoad;
      }
   }
   double bgavg = bgtotal / numPesAvailable;

   int nadjusted = 0;
   for (i=0; i<P; i++) {
      if (processors[i].available) {
        double bgload = processors[i].backgroundLoad;
        if ( bgload < bgavg ) {
          processors[i].backgroundLoad = bgavg;
          ++nadjusted;
        }
      }
   }
   iout << iINFO << "Adjusted background load on " << nadjusted
        << " nodes.\n" << endi;

   computeAverage();  // because otherwise someone will forget
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

int Rebalancer::isAvailableOn(patchInfo *patch, processorInfo *p)
{
   return  p->proxies.find(patch);
}

void Rebalancer::numAvailable(computeInfo *c, processorInfo *p,
           int *nPatches, int *nProxies, int *isBadForCommunication)
{
   //return the number of proxy/home patches available on p for c (0,1,2)

   int patch_count = 0;
   int proxy_count = 0;

   patchInfo &pa1 = patches[c->patch1];
   patchInfo &pa2 = patches[c->patch2];

   int pa1_avail = 1;
   int pa2_avail = 1;

   if (pa1.processor == p->Id) {
     patch_count++;
   } else if ( p->proxies.find(&pa1) ) {
     proxy_count++;
   } else {
     pa1_avail = 0;
   }

   // self computes get one patch for free here
   if (c->patch1 == c->patch2 || pa2.processor == p->Id) {
     patch_count++;
   } else if ( p->proxies.find(&pa2) ) {
     proxy_count++;
   } else {
     pa2_avail = 0;
   }

   *nPatches = patch_count;
   *nProxies = proxy_count;

   int bad = 0;

   if ( patch_count + proxy_count < 2 ) {

     double bgLoadLimit = 0.8 * averageLoad;

     if ( p->backgroundLoad > bgLoadLimit ) bad = 1;
     else {

      int proxiesPerPeLimit = numProxies / numPesAvailable + 3;
      if ( proxiesPerPeLimit < 6 ) proxiesPerPeLimit = 6;

      if ( p->proxies.numElements() > proxiesPerPeLimit ) bad = 1;

      int proxiesPerPatchLimit = numProxies / numPatches + 3;
      if ( proxiesPerPatchLimit < 6 ) proxiesPerPatchLimit = 6;

      if ( ! bad && ! pa1_avail ) {
        if ( processors[pa1.processor].backgroundLoad > bgLoadLimit) bad = 1;
        else if ( pa1.proxiesOn.numElements() > proxiesPerPatchLimit ) bad = 1;
      }
      if ( ! bad && ! pa2_avail ) {
        if ( processors[pa2.processor].backgroundLoad > bgLoadLimit) bad = 1;
        else if ( pa2.proxiesOn.numElements() > proxiesPerPatchLimit ) bad = 1;
      }

     }
   }

   *isBadForCommunication = bad;
}


