/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/


#ifndef RECBISECTION_H

#define RECBISECTION_H

class PatchMap;

/* *********************************************************************** */
/* This class performs a recursive coordinate bisection partitioning       */
/* together with some communication and computation refinements            */
/* *********************************************************************** */

#define MAXNEIGHBOUR 26
class RecBisection 
{
    private:

      typedef struct {                          // a rectangular prism
         float  load;                           // is represented with
         struct { int x,y,z; } origin;          // origin and corner coordinates
         struct { int x,y,z; } corner;
      } Partition;


      typedef struct {                          // the cost of a patch
         float total;                           // is represented here
         float local;
         float edge;
         float icompute[MAXNEIGHBOUR];
      } PatchLoad;

      enum {XDIR=0,YDIR,ZDIR};

      // cost parameters 
      float c_local0;       // fixed cost per patch
      float c_local1;       // cost per atom in the patch
      float c_edge0;        // fixed cost for having a neighbor patch
      float c_edge1;        // cost per atom of the neighboring patch
      float c_icompute0;    // fixed cost per calc. forces for a neighbor
      float c_icompute1;    // cost per atom of the neihgbor that I calc force.

 

      int          numPatches;
      int          npartition; 
      int          currentp;
      Partition    *partitions;
      PatchLoad    *patchload;
      PatchMap *patchMap;     
      Partition    top_partition;

      void compute_patch_load();              // determine cost of each patch
      void rec_divide(int, const Partition&); // recursive  partitioning 
      void assignNodes();                     // assign partitions to nodes
      void assign_nodes_arr(int *);           // assign partitions to array
      void refine_edges();                    
      void refine_boundaries();
      void refine_surface();
      int  prev_better(float,float,float);   

    public:

      RecBisection(int, PatchMap *);
      ~RecBisection();
      int partition(int *);                    // perform partitioning.
					       // if parameter=NULL, store
					       // results in patchDistrib,
					       // otherwise, store in array
};

#endif

