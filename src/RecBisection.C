/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <math.h>
#include "RecBisection.h"
#include "PatchMap.inl"
#include "Patch.h"

/* ********************************************************************* */
/* Constructor for the RecBisection Class                                */
/* ********************************************************************* */
RecBisection::RecBisection(int numpartitions, PatchMap *thePatchMap)
{
    patchMap    = thePatchMap;
    npartition  = numpartitions;
    numPatches  = patchMap->numPatches();
    partitions  = new Partition[npartition];
    patchload   = new PatchLoad[numPatches];
    currentp    = 0;

    if ( partitions == NULL ||
	 patchload == NULL )
    {
      NAMD_die("memory allocation failed in RecBisection::RecBisection");
    }


    // the cost coeffiencients that is used to compute the load introduced
    // to the processor by a patch

    //    c_local0    = 0.1;
    //    c_local1    = 0.015;
    //    c_edge0     = 0.001;
    //    c_edge1     = 0.001;
    //    c_icompute0 = 0.001;
    //    c_icompute1 = 0.000035;

    c_local0    = 0.;
    c_local1    = 1.0;
    c_edge0     = 0.;
    c_edge1     = 0.;
    c_icompute0 = 0.;
    c_icompute1 = 0.;
}


/* ********************************************************************* */
/* Destructor for the RecBisection Class                                */
/* ********************************************************************* */
RecBisection::~RecBisection()
{
    delete [] partitions; 
    delete [] patchload;
}



/* *********************************************************************** */
/* This is a recursive function to partition a 3-D mesh into n cubical     */
/* subpartitions with approximately equal loads.                           */
/*                                                                         */
/* Input  n : total number of subpartitions that the partition "p" has to  */
/*            be divided.  Its value is any integer > 0                    */
/*        p : the current partition                                        */
/*                                                                         */
/* The partition p is divided into two partitions (say p1 and p2) such     */
/* that a) p1 and p2 equally loaded b) p1 and p2 further to be partitioned */
/* ino n1 and n2 subpartitions.                                            */
/* If n is even then n1=n2=n/2. Otherwise n1 is one more than n2.          */
/*                                                                         */
/* Since the subpartitions are rectangualr prisms (not an artibrary shape),*/
/* it is not always possible to find a bisection point where the load      */
/* is equally divided.                                                     */
/* The following strategy is used to get a good partitioning:              */
/* We divide the initial partition along  x,y, and z directions tentatively*/
/* and chose the one that gives the best division                          */
/* *********************************************************************** */

void RecBisection::rec_divide(int n, const Partition &p)
{
    int       i=0,j=0,k=0;              // general purpose index vars
    int       posi[3],posj[3],posk[3];  // division points along x,y,z direction
    int       mindir;                   // the best direction
    int       n1, n2;                   // number of subpartitions in p1 and p2
    int       p1_empty[3], p2_empty[3]; // true if a subpartition is empty
    float     load1;                   // desired load of the first subpartition
    float     prevload,currentload;    
    float     loadarray[3];            // actual loads of p1 (for each x,y,z
                                       // division)
    float     diff;                    // load diffenrence (actual - desired)
    float     mindiff;                 // minimum difference (among directions)
    Partition p1;                      // first subpartition
    Partition p2;                      // second subpartition


    if (n==1)
    {
       // no further subdivision
       // record teh partition p as a final partition
       partitions[currentp++] = p;
       return;
    }

    // calculate division ratio: 1/2 iff n is even, otherwise
    // first partition has more load

    n2 = n/2;
    n1 = n-n2;

    load1  = ( (float) n1/(float) n ) * p.load;

    for(i=XDIR; i<=ZDIR; i++) {p1_empty[i] = p2_empty[i] = 0;}

    p1 = p;
    p2 = p;

    // now try dividing along the x,y,z directions 
    
    // try x-axis 
    currentload = 0.0;
    i = p.origin.x;
    while(currentload < load1 && i<= p.corner.x) 
      {
        prevload = currentload;
        for(j=p.origin.y; j<=p.corner.y; j++)
          for(k=p.origin.z; k<=p.corner.z; k++)
            currentload += patchload[patchMap->pid(i,j,k)].total;
        if (currentload > load1)
           if ( prev_better(prevload,currentload,load1) ) 
               {currentload=prevload; break;}
        i++; 
      }
    posi[XDIR] = i; posj[XDIR] = j; posk[XDIR] = k;
    loadarray[XDIR] = currentload; 
    if (i == p.origin.x) p1_empty[XDIR] = 1;
    if (i > p.corner.x)  p2_empty[XDIR] = 1;


    // z axis
    currentload = 0.0;
    k = p.origin.z;
    while(currentload < load1 && k <= p.corner.z)
      {
        prevload = currentload;
        for(i=p.origin.x; i<=p.corner.x; i++)
          for(j=p.origin.y; j<=p.corner.y; j++)
            currentload += patchload[patchMap->pid(i,j,k)].total;
        if (currentload > load1)
           if ( prev_better(prevload,currentload,load1) )
               {currentload=prevload; break;}
        k++;
      }
    posi[ZDIR] = i; posj[ZDIR] = j; posk[ZDIR] = k;
    loadarray[ZDIR] = currentload; 
    if (k == p.origin.z) p1_empty[ZDIR] = 1;
    if (k > p.corner.z)  p2_empty[ZDIR] = 1;
    

    // y axis
    currentload = 0.0;
    j = p.origin.y;
    while(currentload < load1 && j <= p.corner.y) 
      {
        prevload = currentload;
        for(i=p.origin.x; i<=p.corner.x; i++)
          for(k=p.origin.z; k<=p.corner.z; k++)
            currentload += patchload[patchMap->pid(i,j,k)].total;
        if (currentload > load1)
           if ( prev_better(prevload,currentload,load1) )
               {currentload=prevload; break;}
        j++; 
      }
    posi[YDIR] = i; posj[YDIR] = j; posk[YDIR] = k;
    loadarray[YDIR] = currentload; 
    if (j == p.origin.y) p1_empty[YDIR] = 1;
    if (j > p.corner.y)  p2_empty[YDIR] = 1;

    // determine the best division direction
    mindiff = load1;
    mindir   = -1;
    for(i=XDIR; i<=ZDIR; i++) { 
       diff =  load1 - loadarray[i];
       if (diff < 0.0) diff = -diff;
       if (mindiff >= diff) {mindiff = diff; mindir = i;}
    }

    // always divide along x if possible
    if ( p.origin.x != p.corner.x ) mindir = XDIR;

    // divide along mindir
    switch (mindir) {
      case XDIR: p1.corner.x = posi[XDIR] - 1;
                 p2.origin.x = posi[XDIR];
                 break;
      case YDIR: p1.corner.y = posj[YDIR] - 1;
                 p2.origin.y = posj[YDIR];
                 break;
      case ZDIR: p1.corner.z = posk[ZDIR] - 1;
                 p2.origin.z = posk[ZDIR];
                 break;
      default:   NAMD_bug("RecBisection failing horribly!");
    }
    p1.load = loadarray[mindir];
    p2.load = p.load - p1.load;
    if (!p1_empty[mindir]) rec_divide(n1,p1); 
    if (!p2_empty[mindir]) rec_divide(n2,p2);
}


/* ************************************************************************ */
/* Compute the initial overhead/load of each patch to the processor. The    */
/* load of patches are computed as follows: Each patch has a fixed cost and */
/* a variable cost depending of number atoms it has; c_local0 and c_local1  */
/* respectively. Secondly, due to interaction with each patch, the patch    */
/* incurs some load determined by c_edge0 and c_edge1 cost parameters.      */
/* Finally, each patch has load due to computing forces between its certain */
/* neighbours, with cost coefficients c_icompute0 and c_icompute1           */
/* ************************************************************************ */

void RecBisection::compute_patch_load() 
{  
   int   i,nix,neighbour;
   int   numAtoms, numNeighAtoms;
   float total_icompute;

   for(i=0; i<numPatches; i++) { 

     numAtoms      = patchMap->patch(i)->getNumAtoms();

     patchload[i].total = 0.0;
     patchload[i].edge  = 0.0;

     patchload[i].local = c_local0 + c_local1 * numAtoms;
     patchload[i].local += c_icompute0 + c_icompute1*numAtoms*numAtoms;


     total_icompute = 0.0;

     PatchID neighbors[PatchMap::MaxOneAway + PatchMap::MaxTwoAway];
     int nNeighbors = patchMap->oneOrTwoAwayNeighbors(i,neighbors);
     
     for(nix=0; nix<nNeighbors; nix++) {
       neighbour = neighbors[nix];
       numNeighAtoms = patchMap->patch(neighbour)->getNumAtoms();
       patchload[i].icompute[nix] = 
	 c_icompute0 + c_icompute1*numNeighAtoms*numAtoms;
       total_icompute += patchload[i].icompute[nix]; 
       patchload[i].edge += c_edge0 + c_edge1 * numNeighAtoms; 
     }
     patchload[i].total+=patchload[i].local+total_icompute+patchload[i].edge; 
   }
}




/* *********************************************************************  */
/* Partitions a 3D space. First a recursive algorithm divides the initial */
/* space into rectangular prisms with approximately equal loads.          */
/* Then these rectangular prisms ar emodified to firther increase the     */
/* load balance and redice communication cost                             */
/* *********************************************************************  */

int RecBisection::partition(int *dest_arr)
{
    int i;
  
    top_partition.origin.x = 0; 
    top_partition.origin.y = 0; 
    top_partition.origin.z = 0; 
    top_partition.corner.x  = patchMap->gridsize_a()-1;
    top_partition.corner.y  = patchMap->gridsize_b()-1;
    top_partition.corner.z  = patchMap->gridsize_c()-1;
    top_partition.load      = 0.0;

    // calculate estimated computational load due to each patch
    compute_patch_load();


    for(i=0; i<numPatches; i++) top_partition.load += patchload[i].total;

    // divide into rectangular prisms with load as equal as possible
    rec_divide(npartition,top_partition);

    if (currentp != npartition) 
          return 0;
    else  {
      if (dest_arr==NULL)
	assignNodes();
      else
	assign_nodes_arr(dest_arr);
    }

    return 1;
}


/* ********************************************************************* */
/* partitioning done, update PatchDistib data structure to assign the    */
/* patches to nodes. (patches in partition i is assigned to node i)      */
/* ********************************************************************* */

void RecBisection::assignNodes()
{
    int i,j,k,pix;
    Partition p;

    for(pix=0; pix<npartition; pix++)
    {
       p = partitions[pix];
       for(i=p.origin.x; i<=p.corner.x; i++)
          for(j=p.origin.y; j<=p.corner.y; j++)
             for(k=p.origin.z; k<=p.corner.z; k++)
                patchMap->assignNode(patchMap->pid(i,j,k),pix);  
    }
}

/* ********************************************************************* */
/* partitioning done, save results to an array, rather than updating     */
/* the PatchDistib data structure.  Otherwise, thist is identical to     */
/* assignNodes()                                                         */
/* ********************************************************************* */

void RecBisection::assign_nodes_arr(int *dest_arr)
{
    int i,j,k,pix;
    Partition p;

    for(pix=0; pix<npartition; pix++)
    {
       p = partitions[pix];
       for(i=p.origin.x; i<=p.corner.x; i++)
          for(j=p.origin.y; j<=p.corner.y; j++)
             for(k=p.origin.z; k<=p.corner.z; k++)  {
                dest_arr[patchMap->pid(i,j,k)] = pix;  
	      }
    }
}


/* ************************************************************************ */
/* to be implemented:                                                       */
/* this functions revisits edge directions to refine the load distribution  */
/* ************************************************************************ */
void RecBisection::refine_edges()
{
}


/* ************************************************************************ */
/* to be implemented:                                                       */
/* this function refines the boundaries of subpartitions to redice the      */
/* communication across processors                                          */
/* ************************************************************************ */
void RecBisection::refine_boundaries()
{
}


/* ************************************************************************ */
/* to be implemented:                                                       */
/* refine boundries invokes this function for eeach surface                 */
/* ************************************************************************ */
void RecBisection::refine_surface()
{
}

/* ************************************************************************ */
/* return true if the difference between previous load (prev1) and desired  */
/* load (load1) is less than  teh difference between current an desired     */
/* ************************************************************************ */

int RecBisection::prev_better(float prev, float current, float load1)
{
   float diff1,diff2;

   diff1 = load1 - prev;
   diff2 = current - load1;

   if (diff1 < 0.0) diff1 = -diff1;
   if (diff2 < 0.0) diff2 = -diff2;

   return (diff1 <= diff2);
}

