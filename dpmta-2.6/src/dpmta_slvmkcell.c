/* 
*  dpmta_slvmkcell.c - routines to identify and allocate cell table 
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  these routines are a repackaging of the dpmta_rcvpart.c module,
*  version 2.4, and dpmta_rcvilist.c version 2.3
*/

static char rcsid[] = "$Id: dpmta_slvmkcell.c,v 1.1 1997/09/05 19:42:03 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvmkcell.c,v $
 * Revision 1.1  1997/09/05 19:42:03  jim
 * Original distribution.
 *
 * Revision 2.11  1997/05/12  18:06:09  wrankin
 * added routines to clean up dynamic memory allocation when
 * PMTAinit() is called
 *
 * Revision 2.10  1997/05/07  18:59:36  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.9  1997/02/26  16:54:32  wrankin
 * added support for serial non-pvm compilation of DPMTA
 *
 * Revision 2.8  1996/11/01  02:26:46  wrankin
 * modifications to support cray t3d compilation
 * version update for 2.5 release
 *
 * Revision 2.7  1996/09/24  18:43:01  wrankin
 * many changes for support of version 2.5
 *  - non cubic cells
 *  - non power-of-2 processors
 *  - fixed the resize code
 *  - new virial interface
 *  - changes to macroscopic code (still not working)
 *  - code cleanup
 *  - new test program that reads PDB file
 *  - update code for T3D support
 *
 * Revision 2.6  1996/09/10  18:07:06  wrankin
 * fixed macroscopic processing for non-cubic cells.
 * fixed passing of local expansions for non-2^n processors
 *
 * Revision 2.5  1996/08/20  17:12:51  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.4  1996/08/09  15:31:00  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.3  1995/10/01  21:46:18  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.2  1995/06/27  14:20:29  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.1  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 1.7  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 1.6  1995/01/13  16:46:45  wrankin
 * particle and force arrays are no longer allocated during
 *   cell table creation.
 * removed some extranious procedures
 *
 * Revision 1.5  1994/11/28  04:27:53  wrankin
 * added hooks to handle multiple slave calls to PMTAforce
 * added static allocation of particle data structures
 *
 * Revision 1.4  1994/11/24  13:54:24  wrankin
 * updates to use static particle structure allocation
 * in preparation for having multiple slave entry points
 *
 * Revision 1.3  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 1.2  1994/10/19  00:08:16  wrankin
 * allocation of cells is now done through bill elliotts routines
 * from dpmta_slvmultipole.c
 *
 * Revision 1.1  1994/10/14  05:10:51  wrankin
 * Initial revision
 *
 *
*/

/* include files */
#include <stdio.h>
#include "dpmta_cell.h"
#include "dpmta_slave.h"

/* a little prototyping */
int Cell2Cell( int, int, int, int *, int * );


/****************************************************************
*
*  this funtion will allocate and initialize the cell table for
*  a specific slave processor.  the individual cell and particle
*  structures will be allocated for the cells owned by that process.
*
*/

Alloc_Cell_Table()
{

   int  num_cells;         /* total number of cells */
   int  i, j, k;           /* loop counters */
   CellPtr celltemp;       /* temporary cell pointer */
   MdataPtr mdatatemp;     /* temporary mdata pointer */

   /*
   *  initialize the cell indices
   *  these indices identify which cells (for each level) that
   *  this processor owns.
   */

   cell_identify();

   /* allocate array of cell pointers */

   Dpmta_CellTbl = (CellPtrPtrPtr)malloc(Dpmta_NumLevels*sizeof(CellPtrPtr));
   if ( Dpmta_CellTbl == NULL ) {
      fprintf(stderr,"Alloc_Cell_Table(): malloc failed [1]\n");
      exit(-1);
   }
   num_cells = Dpmta_LevelLocate[Dpmta_NumLevels];
   Dpmta_CellTbl[0] = (CellPtrPtr)malloc(num_cells*sizeof(CellPtr));
   if ( Dpmta_CellTbl[0] == NULL ) {
      fprintf(stderr,"Alloc_Cell_Table(): malloc failed [2]\n");
      exit(-1);
   }
   for ( i=1; i<Dpmta_NumLevels; i++ )
      Dpmta_CellTbl[i] = &(Dpmta_CellTbl[0][Dpmta_LevelLocate[i]]);
   for ( i=0; i<num_cells; i++ )
      Dpmta_CellTbl[0][i] = (CellPtr) NULL;

   /*
   *  allocate the cell table entries owned by the process
   */

   /* cycle through each level */
   for ( i=0; i<Dpmta_NumLevels; i++ ) {
      if ( (Dpmta_Scell[i]) != -1 ) {

         num_cells = Dpmta_Ecell[i] - Dpmta_Scell[i] + 1;

         /*
         *  all allocation is pretty much done.
         *  cycle through and setup pointers and initial values
         */

         for ( j=0; j<num_cells; j++ ) {
            k = j + Dpmta_Scell[i];

	    celltemp = (CellPtr)malloc(sizeof(Cell));
	    if ( celltemp == NULL ) {
	       fprintf(stderr,"Alloc_Cell_Table(): malloc() failed [3]\n");
	       exit(-1);
	    }

	    mdatatemp = (MdataPtr)malloc(sizeof(Mdata));
	    if ( mdatatemp == NULL ) {
	       fprintf(stderr,"Alloc_Cell_Table(): malloc() failed [4]\n");
	       exit(-1);
	    }

	    Dpmta_CellTbl[i][k] = celltemp;
            Dpmta_CellTbl[i][k]->mdata = mdatatemp;

            if (Dpmta_FFT)
               CallocF(&(Dpmta_CellTbl[i][k]->m), Dpmta_Mp, Dpmta_FftBlock);
            else
               Calloc(&(Dpmta_CellTbl[i][k]->m), Dpmta_Mp);
            Calloc(&(Dpmta_CellTbl[i][k]->mdata->l), Dpmta_Mp);

            Dpmta_CellTbl[i][k]->pid = Dpmta_Pid;
            Dpmta_CellTbl[i][k]->id = Dpmta_LevelLocate[i] + k;
            Dpmta_CellTbl[i][k]->mvalid = FALSE;
            Dpmta_CellTbl[i][k]->plist = NULL;
            Dpmta_CellTbl[i][k]->psize = 0;
            Dpmta_CellTbl[i][k]->n = 0;
            Dpmta_CellTbl[i][k]->mdata->lvalid = FALSE;
            Dpmta_CellTbl[i][k]->mdata->flist = NULL;
            Dpmta_CellTbl[i][k]->mdata->part_id = NULL;
            Dpmta_CellTbl[i][k]->mdata->proc_id = NULL;

#ifdef COMP_LJ
            LJalloc(&(Dpmta_CellTbl[i][k]->m_lj), Dpmta_Mp_LJ);
            LJalloc(&(Dpmta_CellTbl[i][k]->mdata->l_lj), Dpmta_Mp_LJ);
            Dpmta_CellTbl[i][k]->mdata->f_lj = NULL;
#endif

         } /* for j */

         /*
         *  now, if we are not at the top level, we need to make sure
         *  that a parent has been allocated for each cell on the level
         *  that we just allocated.  we must have a local parent of our
         *  cells to hold the local mpe expansion.
         *
         *  the below implementation allocates only the parent cell
         *  one level above topmost cell the processor owns (if this
         *  cell exists - it doesn't for pid 0)
         */

         if ( i > 0 ) {

 	    /* cycle through all children and allocate parents */
 	    for ( j=Dpmta_Scell[i]; j<=Dpmta_Ecell[i]; j++ ) {
               k = getparent(j);
	       if ( Dpmta_CellTbl[i-1][k] == (CellPtr)NULL ) {

                  /* allocate cell */
                  celltemp = (CellPtr)malloc(sizeof(Cell));
		  if ( celltemp == (CellPtr)NULL ) {
		     fprintf(stderr,
                        "Alloc_Cell_Table(): malloc() failed [5]\n");
		     exit(-1);
                  }

                  mdatatemp = (MdataPtr)malloc(sizeof(Mdata));
                  if ( mdatatemp == NULL ) {
                     fprintf(stderr,
                        "Alloc_Cell_Table(): malloc() failed [6]\n");
                     exit(-1);
                  }

                  /*
                  *  all allocation is pretty much done.
                  *  cycle through and setup pointers and initial values
                  */

                  Dpmta_CellTbl[i-1][k] = celltemp;
                  Dpmta_CellTbl[i-1][k]->mdata = mdatatemp;

                  if (Dpmta_FFT)
                     CallocF(&(Dpmta_CellTbl[i-1][k]->m),
			     Dpmta_Mp, Dpmta_FftBlock);
		  else
		     Calloc(&(Dpmta_CellTbl[i-1][k]->m), Dpmta_Mp);

                  Calloc(&(Dpmta_CellTbl[i-1][k]->mdata->l), Dpmta_Mp);

		  Dpmta_CellTbl[i-1][k]->pid = -1;
                  Dpmta_CellTbl[i-1][k]->id = Dpmta_LevelLocate[i-1] + k;
		  Dpmta_CellTbl[i-1][k]->mvalid = FALSE;
		  Dpmta_CellTbl[i-1][k]->plist = NULL;
		  Dpmta_CellTbl[i-1][k]->psize = 0;
		  Dpmta_CellTbl[i-1][k]->n = 0;
		  Dpmta_CellTbl[i-1][k]->mdata->lvalid = FALSE;
		  Dpmta_CellTbl[i-1][k]->mdata->flist = NULL;
		  Dpmta_CellTbl[i-1][k]->mdata->part_id = NULL;
		  Dpmta_CellTbl[i-1][k]->mdata->proc_id = NULL;

#ifdef COMP_LJ
		  LJalloc(&(Dpmta_CellTbl[i-1][k]->m_lj), Dpmta_Mp_LJ);
		  LJalloc(&(Dpmta_CellTbl[i-1][k]->mdata->l_lj), Dpmta_Mp_LJ);
		  Dpmta_CellTbl[i-1][k]->mdata->f_lj = NULL;
#endif
	       } /* if Dpmta_CellTbl */
	    } /* for j */
         } /* if i */
      } /* if Dpmta_Scell */
   } /* for i */

   /*
   *  get size of the multipole and local expansions.  these are
   *  used later for sending and receiving MPE and local expansions.
   *  also allocate temporary buffers used to send cell MPEs
   *  between parent and child pids.
   *
   *  actually this should probably be moved to MultipoleSetup()
   */

   if (Dpmta_FFT) {
      Dpmta_MpeSize = CsizeF(Dpmta_Mp);
      CallocF(&(Dpmta_Temp_Mpe), Dpmta_Mp, Dpmta_FftBlock);
   }
   else {
      Dpmta_MpeSize = Csize(Dpmta_Mp);
      Calloc(&(Dpmta_Temp_Mpe), Dpmta_Mp);
   }

   Dpmta_LclSize = Csize(Dpmta_Mp);


#ifdef COMP_LJ
   Dpmta_MpeSize_LJ = LJsize(Dpmta_Mp_LJ);
   Dpmta_LclSize_LJ = LJsize(Dpmta_Mp_LJ);
   LJalloc(&(Dpmta_Temp_Mpe_LJ), Dpmta_Mp_LJ);
#endif

} /* Alloc_Cell_TAble */



/****************************************************************
*
*  the following procedure will compute an array of indices
*  into the cell table which identify, for a given level, 
*  which cells that processor owns.
*
*  the algorithm cycles through each level.  for a given level, 
*  if there are more processors than there are cells at that level,
*  then this processor may own at most a single cell.
*
*  otherwise, the process owns one or more cells at that level,
*  and we can compute the range of cells it owns.  note that
*  this makes the assumption that the cells owned by a process are
*  all adjacent within that level.
*
*  in addition, the procedure computes an array of the number of
*  multipole expansions that are received from other pids as a result
*  of the upward M2M pass of the algorithm.
*/

cell_identify()
{

   int i, j, k, l, m, n;  /* loop counters, what else? */
   int num_cells;         /* number of cells in that level */

   for (i=0; i<Dpmta_NumLevels; i++) {
      num_cells = 0x1 << (3*i);

      Dpmta_Scell[i] =
         num_cells - (num_cells*(Dpmta_Nproc-Dpmta_Pid))/Dpmta_Nproc;
      Dpmta_Ecell[i] =
         num_cells - (num_cells*(Dpmta_Nproc-Dpmta_Pid-1))/Dpmta_Nproc - 1;

      /*
      *  check to see if we don't have any cells at this level,
      *  in which case the pid computed from the Scell will not
      *  match the pid of this process
      *
      *  if this is the case, we should only have one cell per pid.
      *  and this is an assumption that we make.
      */

      if ( Dpmta_Nproc > num_cells )  {
         j = getslvpid(i,Dpmta_Scell[i]);
         if ( Dpmta_Pid != j )
            Dpmta_Scell[i] = -1;
         Dpmta_Ecell[i] = Dpmta_Scell[i];
      } /* if Dpmta_Nproc */

   } /* for i */

   /*
   *  compute multipole receive index -
   *
   *  cycle through all the levels (except the bottom).
   *  for each of the children of each cell, count the number of
   *  other pids which own those children.  this will be the number
   *  of processes that will send multipole expansions over the course
   *  of the upward pass.
   *
   *  note that like the start and end cell arrays computed above, this
   *  array will need to be recomputed anytime cells are reassigned to
   *  new processors.
   *
   *  if we are not doing PBC or Macroscopic Expansions, then
   *  there is no reason to propagate the MPEs up past one
   *  level above our max.
   */

   for (i=0; i<Dpmta_NumLevels-1; i++) {
      Dpmta_RMcell[i] = 0;
   }

   for (i=(Dpmta_DownPassStart-1); i<Dpmta_NumLevels-1; i++) {
      if ( Dpmta_Scell[i] != -1 ) {
         for (j=Dpmta_Scell[i]; j<=Dpmta_Ecell[i]; j++) {
            k = getfirstchild(j);
            m = -1;
            for (l=k; l<(k+8); l++) {
               n = getslvpid(i+1,l);
               if ( m != n ) {
                  m = n;
                  if ( n != Dpmta_Pid ) {
	             Dpmta_RMcell[i]++;
                  } /* if n */
               } /* if m */
            } /* for l */
         } /* for j */
      } /* if Dpmta_Scell */
   } /* for i */


   /*
   *  compute local exp receive index -
   *
   *  cycle through all the levels (except the top).  for each of
   *  the unique parents, count the number of parents which are
   *  owned by other pids.  this will be the number of processes that
   *  will send local expansions over the course of the downward
   *  pass.
   *
   *  note that like the start and end cell arrays computed above, this
   *  array will need to be recomputed anytime cells are reassigned to
   *  new processors.
   */

   for (i=0; i<Dpmta_NumLevels; i++) {
      Dpmta_RLcell[i] = 0;
   }

   for (i=Dpmta_DownPassStart; i<Dpmta_NumLevels; i++) {
      if ( Dpmta_Scell[i] != -1 ) {
         l = -1;
         for (j=Dpmta_Scell[i]; j<=Dpmta_Ecell[i]; j++) {
            k = getparent(j);
            if ( l != k ) {
	       l = k;
               m = getslvpid(i-1,k);
               if ( m != Dpmta_Pid ) {
                  Dpmta_RLcell[i]++;
               } /* if m */
            } /* if l */
         } /* for j */
      } /* if Dpmta_Scell */
   } /* for i */

   /*
    * there is no need ot pass locals at the upper levels
    * if we are not doing PBC or Macroscopic Expansions, so
    * filter these out.
    */

#ifdef MACROSCOPIC
   if ( Dpmta_K == 0 ) {
      Dpmta_RLcell[1] = 0;
   }
#else
   Dpmta_RLcell[1] = 0;
#endif

   if ( Dpmta_PBC == 0 ) {
      Dpmta_RLcell[2] = 0;
   }

} /* cell_identify */


/****************************************************************
*
*  cell_center - set the center of the cell table
*
*/

cell_center( int level, int cell )
{

   int i;                  /* loop counter */
   int mul;                /* number of cells per edge */
   int index;              /* cell id */
#ifndef PIPED
   int x,y,z;              /* cell coordinates */
   Vector cube;            /* size of the smallest cube */
#else
   int v1, v2, v3;
   double v1len, v2len, v3len, shiftc;
   double imul;
   Vector scaledv1, scaledv2, scaledv3;
   Vector tmpc;
#endif

#ifndef PIPED
   x = y = z = 0;

   mul = 1 << level;
   index = cell;
   for ( i=0; i<level; i++ ) {
      x |= (index & 01) << i;
      index >>= 1;
      y |= (index & 01) << i;
      index >>= 1;
      z |= (index & 01) << i;
      index >>= 1;
   }


   /* cube the length of a cube edge at that level. */

   cube.x = (1.0/((double) mul))* (Dpmta_CellLength.x / Dpmta_MaxCellLen);
   cube.y = (1.0/((double) mul))* (Dpmta_CellLength.y / Dpmta_MaxCellLen);
   cube.z = (1.0/((double) mul))* (Dpmta_CellLength.z / Dpmta_MaxCellLen);

   Dpmta_CellTbl[level][cell]->p.x = cube.x*((double)x + 0.5);
   Dpmta_CellTbl[level][cell]->p.y = cube.y*((double)y + 0.5);
   Dpmta_CellTbl[level][cell]->p.z = cube.z*((double)z + 0.5);

#else
   v1 = v2 = v3 = 0;

   mul = 1 << level;
   index = cell;
   for ( i=0; i<level; i++ ) {
      v1 |= (index & 01) << i;
      index >>= 1;
      v2 |= (index & 01) << i;
      index >>= 1;
      v3 |= (index & 01) << i;
      index >>= 1;
   }

   scaledv1.x = Dpmta_PVector1.x / Dpmta_MaxCellLen;
   scaledv1.y = Dpmta_PVector1.y / Dpmta_MaxCellLen;
   scaledv1.z = Dpmta_PVector1.z / Dpmta_MaxCellLen;
   scaledv2.x = Dpmta_PVector2.x / Dpmta_MaxCellLen;
   scaledv2.y = Dpmta_PVector2.y / Dpmta_MaxCellLen;
   scaledv2.z = Dpmta_PVector2.z / Dpmta_MaxCellLen;
   scaledv3.x = Dpmta_PVector3.x / Dpmta_MaxCellLen;
   scaledv3.y = Dpmta_PVector3.y / Dpmta_MaxCellLen;
   scaledv3.z = Dpmta_PVector3.z / Dpmta_MaxCellLen;

   tmpc.x = (scaledv1.x * ((double)v1 + 0.5)) +
            (scaledv2.x * ((double)v2 + 0.5)) +
            (scaledv3.x * ((double)v3 + 0.5));
   tmpc.y = (scaledv1.y * ((double)v1 + 0.5)) +
            (scaledv2.y * ((double)v2 + 0.5)) +
            (scaledv3.y * ((double)v3 + 0.5));
   tmpc.z = (scaledv1.z * ((double)v1 + 0.5)) +
            (scaledv2.z * ((double)v2 + 0.5)) +
            (scaledv3.z * ((double)v3 + 0.5));
   imul = 1.0 / (double)mul;

   Dpmta_CellTbl[level][cell]->p.x = tmpc.x * imul;
   Dpmta_CellTbl[level][cell]->p.y = tmpc.y * imul;
   Dpmta_CellTbl[level][cell]->p.z = tmpc.z * imul;
#endif
} /* cell_center */



/****************************************************************
*
*  Alloc_Ilist_Cells() - allocate the cell table entries for the
*     cells pointed to in the interaction lists.
*
*  this routine will traverse the set of 'private' cells in the
*  cell table and for each interaction list, check to make sure
*  that the complete set of cell table entries are allocated
*  for the remote cells that will be accessed by this processor.
*
*  note that we do not look at the double direct interaction list
*  since all the cells in this table are defined as being located
*  local to the current processor, and are already allocated.
*/

Alloc_Ilist_Cells()
{
   int i,j,k;               /* loop counters */
   int temp1, temp2;        /* temp cell ids */
   int sep;                 /* cell separation index */
   int pcell, plevel;       /* parent cell id */
   int ovfl;                /* overflow cell boundary flag */
   Complex *mpe_block;      /* temporary mpe block pointer */
   ParticlePtr ptemp;       /* temporary particle list pointer */


   /*  cycle through all the levels */
   for ( i=Dpmta_DownPassStart; i<Dpmta_NumLevels; i++ ) {

      /*
      *  see if we own any cells at this level. if we do
      *  then go through their interaction lists and allocate any
      *  needed cells
      */

      if ( Dpmta_Scell[i] != -1 ) {
         for ( j=Dpmta_Scell[i]; j<=Dpmta_Ecell[i]; j++ ) {


	    /*
	    *  determine relative cell position
	    */
	    temp1 = j & 0x07;

            /*
            *  check direct interaction lists
            */
            for ( k=0; k<Dpmta_Intlist[temp1].dcnt; k++ ) {
	       sep = Dpmta_Intlist[temp1].dlist[k];
               if ( Cell2Cell(i,j,sep,&temp2,&ovfl) ) {
                  if ( Dpmta_CellTbl[i][temp2] == (CellPtr)NULL ) {
                     Dpmta_CellTbl[i][temp2] = (CellPtr)malloc(sizeof(Cell));
                     if ( Dpmta_CellTbl[i][temp2] == (CellPtr)NULL ) {
                        fprintf(stderr,
				"Alloc_Ilist_Cells(): malloc() failed [1]\n");
                        exit(-1);
                     }
                     Dpmta_CellTbl[i][temp2]->pid = -1;
                     Dpmta_CellTbl[i][temp2]->id = Dpmta_LevelLocate[i] + temp2;
                     Dpmta_CellTbl[i][temp2]->plist = (ParticlePtr)NULL;
                     Dpmta_CellTbl[i][temp2]->psize = 0;
                     Dpmta_CellTbl[i][temp2]->n = 0;
                     Dpmta_CellTbl[i][temp2]->m = (Mtype)NULL;
                     Dpmta_CellTbl[i][temp2]->mvalid = FALSE;
                     Dpmta_CellTbl[i][temp2]->mdata = (MdataPtr)NULL;
#ifdef COMP_LJ
                     Dpmta_CellTbl[i][temp2]->m_lj = (MtypeLJ)NULL;
#endif

                 } /* if Dpmta_CellTbl */
	       } /* if Cell2Cell() */
            } /* for k */

            /*
            *  check multipole sibling interaction lists
            */
            for ( k=0; k<Dpmta_Intlist[temp1].scnt; k++ ) {
	       sep = Dpmta_Intlist[temp1].slist[k];
               if ( Cell2Cell(i,j,sep,&temp2,&ovfl) ) {
                  if ( Dpmta_CellTbl[i][temp2] == (CellPtr)NULL ) {
                     Dpmta_CellTbl[i][temp2] = (CellPtr)malloc(sizeof(Cell));
                     if ( Dpmta_CellTbl[i][temp2] == (CellPtr)NULL ) {
                        fprintf(stderr,
				"Alloc_Ilist_Cells(): malloc() failed [2]\n");
                        exit(-1);
                     }

                     Dpmta_CellTbl[i][temp2]->pid = -1;
                     Dpmta_CellTbl[i][temp2]->id = Dpmta_LevelLocate[i] + temp2;
                     Dpmta_CellTbl[i][temp2]->plist = (ParticlePtr)NULL;
                     Dpmta_CellTbl[i][temp2]->psize = 0;
                     Dpmta_CellTbl[i][temp2]->n = 0;
                     Dpmta_CellTbl[i][temp2]->mdata = (MdataPtr)NULL;
                     Dpmta_CellTbl[i][temp2]->mvalid = FALSE;
                     Dpmta_CellTbl[i][temp2]->m = (Mtype)NULL;
#ifdef COMP_LJ
                     Dpmta_CellTbl[i][temp2]->m_lj = (MtypeLJ)NULL;
#endif
   	          } /* if Dpmta_CellTbl[][] */

                  if ( Dpmta_CellTbl[i][temp2]->m == (Mtype)NULL ) {
                     if (Dpmta_FFT)
                        CallocF(&(Dpmta_CellTbl[i][temp2]->m), Dpmta_Mp,
				Dpmta_FftBlock);
                     else
                        Calloc(&(Dpmta_CellTbl[i][temp2]->m), Dpmta_Mp);
                  } /* if Dpmta_CellTbl[][]->m */

#ifdef COMP_LJ
                  if ( Dpmta_CellTbl[i][temp2]->m_lj == (MtypeLJ)NULL ) {
                     LJalloc(&(Dpmta_CellTbl[i][temp2]->m_lj), Dpmta_Mp_LJ);
                  } /* if Dpmta_CellTbl[][]->m_lj */
#endif

	       } /* if Cell2Cell */
	    } /* for k */

            /*
            *  check multipole parental interaction lists
	    *
	    *  NOTE: ONLY IF (i>0) !!!!!!
            */
            for ( k=0; k<Dpmta_Intlist[temp1].pcnt; k++ ) {
	       sep = Dpmta_Intlist[temp1].plist[k];
	       pcell = getparent(j);
	       plevel = i-1;
               if ( Cell2Cell(plevel,pcell,sep,&temp2,&ovfl) ) {
                  if ( Dpmta_CellTbl[plevel][temp2] == (CellPtr)NULL ) {
                     Dpmta_CellTbl[plevel][temp2] =
		       (CellPtr)malloc(sizeof(Cell));
                     if ( Dpmta_CellTbl[plevel][temp2] == (CellPtr)NULL ) {
                        fprintf(stderr,
				"Alloc_Ilist_Cells(): malloc() failed [3]\n");
                        exit(-1);
                     }

                     Dpmta_CellTbl[plevel][temp2]->pid = -1;
                     Dpmta_CellTbl[plevel][temp2]->id =
		       Dpmta_LevelLocate[plevel] + temp2;
                     Dpmta_CellTbl[plevel][temp2]->plist = (ParticlePtr)NULL;
                     Dpmta_CellTbl[plevel][temp2]->psize = 0;
                     Dpmta_CellTbl[plevel][temp2]->n = 0;
                     Dpmta_CellTbl[plevel][temp2]->mdata = (MdataPtr)NULL;
                     Dpmta_CellTbl[plevel][temp2]->mvalid = FALSE;
                     Dpmta_CellTbl[plevel][temp2]->m = (Mtype)NULL;

#ifdef COMP_LJ
                     Dpmta_CellTbl[plevel][temp2]->m_lj = (MtypeLJ)NULL;
#endif

   	          } /* if Dpmta_CellTbl[][] */

                  if ( Dpmta_CellTbl[plevel][temp2]->m == (Mtype)NULL ) {
                     if (Dpmta_FFT)
                        CallocF(&(Dpmta_CellTbl[plevel][temp2]->m), Dpmta_Mp,
				Dpmta_FftBlock);
                     else
                        Calloc(&(Dpmta_CellTbl[plevel][temp2]->m), Dpmta_Mp);
                  } /* if Dpmta_CellTbl[][]->m */

#ifdef COMP_LJ
                  if ( Dpmta_CellTbl[plevel][temp2]->m_lj == (MtypeLJ)NULL ) {
                     LJalloc(&(Dpmta_CellTbl[plevel][temp2]->m_lj),
			     Dpmta_Mp_LJ);
                  } /* if Dpmta_CellTbl[][]->m_lj */
#endif

	       } /* if Cell2Cell */
	    } /* for k */

	 } /* for j */
      } /* if Dpmta_Scell[] */
   } /* for i */
 
} /* Alloc_Ilist_Cells */


/****************************************************************
*
*  Make_Cell_Centers() - cycle through all cells and set centers
*    of the cells.
*/

void Make_Cell_Centers()
{

   int i,j,k;

   for (i=0; i<Dpmta_NumLevels; i++) {
      for (j=0; j<Dpmta_Power8[i]; j++) {
         if ( Dpmta_CellTbl[i][j] != NULL ) {
 	    cell_center(i,j);
	 } /* if */
      } /* for j */
   }/* for i */

} /* Make_Cell_Centers */



/****************************************************************
*
*  Delete_Cell_Table() - free all dynamic data structures
*
*  this funtion will delete the cell table and all associated
*  structures.  it traverses the cell table, free()-ing each data
*  structure as it finds them.
*
*  it's slow and plodding, but at this stage of the endgame, who cares?
*
*/

Delete_Cell_Table()
{

   int i, j, k;    /* loop counters */
   int num_cells;  /* size of cell table */

   num_cells = Dpmta_LevelLocate[Dpmta_NumLevels];

   for ( i=0; i<num_cells; i++ ) {
      if ( Dpmta_CellTbl[0][i] != (CellPtr) NULL ) {

	 /* check if we have an mdata structure for this cell */

	 if ( Dpmta_CellTbl[0][i]->mdata != (MdataPtr) NULL ) {

	    if ( Dpmta_CellTbl[0][i]->mdata->flist != NULL ) {
	       free( Dpmta_CellTbl[0][i]->mdata->flist );
	    }

#ifdef COMP_LJ
	    if ( Dpmta_CellTbl[0][i]->mdata->f_lj != NULL ) {
	       free( Dpmta_CellTbl[0][i]->mdata->f_lj );
	    }
#endif	    
            if ( Dpmta_CellTbl[0][i]->mdata->part_id != NULL ) {
	       free( Dpmta_CellTbl[0][i]->mdata->part_id );
	    }
	    
            if ( Dpmta_CellTbl[0][i]->mdata->proc_id != NULL ) {
	       free( Dpmta_CellTbl[0][i]->mdata->proc_id );
	    }
	    
            if ( Dpmta_CellTbl[0][i]->mdata->l != NULL ) {
	       Cfree( Dpmta_CellTbl[0][i]->mdata->l, Dpmta_Mp );
	    }

#ifdef COMP_LJ
            if ( Dpmta_CellTbl[0][i]->mdata->l_lj != NULL ) {
	       LJfree( Dpmta_CellTbl[0][i]->mdata->l_lj, Dpmta_Mp_LJ );
	    }
#endif

	    /* free up the mdata structure */
	    free( Dpmta_CellTbl[0][i]->mdata );

	 } /* if mdata */

	 /* check if we have a particle list and free it */
	 if ( Dpmta_CellTbl[0][i]->plist != NULL ) {
	    free( Dpmta_CellTbl[0][i]->plist );
	 }

	 /* free up multipoles */
	 if ( Dpmta_CellTbl[0][i]->m != NULL ) {
	    if (Dpmta_FFT) {
	       CfreeF( Dpmta_CellTbl[0][i]->m, Dpmta_Mp, Dpmta_FftBlock );
	    }
	    else {
       	       Cfree( Dpmta_CellTbl[0][i]->m, Dpmta_Mp );
	    }
	 }
	    
#ifdef COMP_LJ
	 if ( Dpmta_CellTbl[0][i]->m_lj != NULL ) {
	    LJfree( Dpmta_CellTbl[0][i]->m_lj, Dpmta_Mp_LJ );
	 }
#endif	    

	 /* finally, free the cell */
	 free( Dpmta_CellTbl[0][i] );

      } /* if cell */
   } /* for i */

   /* finally, free the cell table array */
   free( Dpmta_CellTbl[0] );
   free( Dpmta_CellTbl );

   /*
    *  here are some other data arrays that are allocated
    *  somewhere in this module and thus need to be free
    */

   if (Dpmta_FFT) {
      CfreeF(Dpmta_Temp_Mpe, Dpmta_Mp, Dpmta_FftBlock);
   }
   else {
      Cfree(Dpmta_Temp_Mpe, Dpmta_Mp);
   }
   
#ifdef COMP_LJ
   LJfree(Dpmta_Temp_Mpe_LJ, Dpmta_Mp_LJ);
#endif



} /* Delete_Cell_Table */

