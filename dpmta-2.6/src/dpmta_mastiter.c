/****************************************************************
*
*  dpmta_mastiter.c - master iteration routines
*
*  w.t.rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*  these routines handle the distribution and reception of data to
*  and from the slave processes during a single PMTA iteration.
*
*/

static char rcsid[]="$Id: dpmta_mastiter.c,v 1.4 1997/09/29 23:58:38 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_mastiter.c,v $
 * Revision 1.4  1997/09/29 23:58:38  jim
 * Incorporated changes from version 2.6.1 of DPMTA.
 *   - fixes for bad handling of empty/invalid multipoles when
 *     using large processor sets.
 *   - moved functions that provide data mapping to processors.  master
 *     and slave routines now call the same function in dpmta_distmisc.c
 * Also, switched pvmc.h back to pvm3.h.
 *
 * Revision 1.3  1997/09/28 10:24:08  milind
 * Fixed Makefiles to include arch params from a single file.
 *
 * Revision 1.2  1997/09/12 22:56:31  jim
 * Modifications to work with converse pvm.
 *
 * Revision 1.1  1997/09/05 19:41:54  jim
 * Original distribution.
 *
 * Revision 2.22  1997/05/12  18:06:03  wrankin
 * added routines to clean up dynamic memory allocation when
 * PMTAinit() is called
 *
 * Revision 2.21  1997/05/08  18:38:45  chumphre
 * ||-Piped bug fix in sorting code
 *
 * Revision 2.20  1997/05/07  21:27:49  chumphre
 * cubectr changed to cellctr and other minor changes related to the interface
 *
 * Revision 2.19  1997/05/07  18:59:21  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.18  1997/04/28  16:09:25  wrankin
 * implemented multi-buffer messaging for distribution particle information
 *   to DPMTA slave processes.  should be much more efficient.
 *
 * Revision 2.17  1997/02/26  20:43:26  wrankin
 * updated timing measurements and placed routines in single module
 *
 * Revision 2.16  1997/02/26  18:52:32  wrankin
 * cleaned up handling of cell resizing
 *
 * Revision 2.15  1997/02/24  19:12:20  wrankin
 * fixes to timing code
 *
 * Revision 2.14  1997/01/28  17:10:39  wrankin
 * added new timing codes for macroscopic expansion calculations
 *
 * Revision 2.13  1996/11/18  19:29:27  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.12  1996/11/14  17:50:15  wrankin
 * performance enhancements to multipole M2L routine.
 * additions to make slaves exit gracefully.
 *
 * Revision 2.11  1996/10/28  23:01:55  wrankin
 * additions to test routines to provide processing paramenters as
 *   part of output file.
 * dpmta_direct will now perform macroscopic computations, reading
 *   in processing parameters from particle position file.
 *
 * Revision 2.10  1996/10/18  17:04:43  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.9  1996/09/24  18:41:55  wrankin
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
 * Revision 2.8  1996/08/20  17:12:33  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.7  1996/08/09  15:30:43  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.6  1996/02/29  21:13:29  wrankin
 * New relaease: 2.4 (.1)
 *    - simplified calling structure for initialization
 *    - macroscopic periodic code
 *    - fixed PBC calculation - all particles are now stored in the
 *      cell as positions relative to the cell center.
 *    - virial preasure tensor computed
 *    - fix to allow particles on the outer cube boundary to be
 *      included. (UNC fix)
 *    - fix to order reception of particle data during the distributed
 *      calling sequence (UIUC fix)
 *    - removed M2L code that didn't use transfer matrices
 *    - early hooks in to perform interaction list sorting.
 *    - fixed LJ scaling factor for 1/r^12 potential.
 *    - cleaned up the LJ interface.
 *    - and of course, my continued efforts to ANSI-fy this beast.
 *
 * Revision 2.5  1995/12/08  23:00:16  wrankin
 * preliminary release of DPMTA 2.3
 *   - added working Periodic Boundary Conditions (PDC) flag
 *   - added code for Virial Computation (not working yet)
 *   - general cleanup of some modules
 *
 * Revision 2.4  1995/11/29  22:29:09  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.3  1995/10/01  21:45:53  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.2  1995/09/18  19:02:30  wrankin
 * updates to move the scaling of the simulation box to the slave
 * process.
 *
 * Revision 2.1  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 1.6  1995/04/09  22:13:02  wrankin
 * addition of some static global arrays for support of T3D port
 *
 * Revision 1.5  1994/12/31  21:36:42  wrankin
 * Send_Slave_Particles() now check for num_parts==0 before
 * freeing cell_index array
 *
 * Revision 1.4  1994/12/06  19:07:17  wrankin
 * added code to handle offset cell centers
 * moved there the cell_id is placed in slave message for consistency
 *
 * Revision 1.3  1994/12/01  21:16:52  wrankin
 * added check to handle zero particle call to PMTAforce()
 *
 * Revision 1.2  1994/11/30  18:03:33  wrankin
 * added processing to handle distributed calling sequence
 * timing information (if compiled in) outputs to stderr.
 *
 * Revision 1.1  1994/10/26  02:39:56  wrankin
 * Initial revision
 *
 *
*/

#include <stdio.h>
#include "pvm3.h"
#include "dpmta.h"
#include "dpmta_pvm.h"

/*
*  static indexing arrays
*/

static int Power8[] = { 1, 8, 64, 512, 4096, 32768, 262144, 2097152,
                        16777216, 134217728, 1073741824 };

static int LevelLocate[] = { 0, 1, 9, 73, 585, 4681, 37449, 299593,
                        2396745, 19173961, 153391689 };

/*
*  static communication buffers for sending and receiving
*  particle data.  Make them global so that we have them all
*  in one place in case we need to free() them later.
*/

static int         *SendPartCnt = NULL;
static int         SendPartCntSz = 0;
static int         *SendCellId = NULL;
static int         SendCellIdSz = 0;

static PmtaPartInfoPtr RecvBuf = NULL;
static int         RecvBufSz = 0;

static int         *SendMsgBuf = NULL;
static int         Msg_Term = -1;


#if defined VIRIAL || defined OLDVIRIAL
static double      VirPot;
static PmtaVector  VirForce;
#ifdef COMP_LJ
static double      VirPot_LJ;
static PmtaVector  VirForce_LJ;
#endif
#endif


/****************************************************************
*
*  Send_Slave_Particles() - sort and distribute particle lists
*     to the slave processes.
*
*  this process will take a list of particles, scale the particle
*  positions to fit within the unit cube used by the slave process,
*  and then send the particles to the appropriate slave processes.
*
*  in order to sort the particle data, the routine simply traverses
*  the entire list once for each slave process and packs those
*  particles belonging to that slave.  not the most efficient 
*  distribution scheme, but it is simple and it works.
*
*/

Send_Slave_Particles(
   int num_parts,               /* number of particles */
   PmtaParticlePtr part_table,  /* particle table */
   int num_proc,                /* number of slave processes */
   int *tids,                   /* array of slave process ids */
   int num_levels,              /* number of levels in oct-tree */
   int resize,                  /* resize flag */
#ifndef PIPED
   PmtaVector *cubelen,         /* length of cell cube sides */
#else
   PmtaVector *PVector1,
   PmtaVector *PVector2,
   PmtaVector *PVector3,
#endif
   PmtaVector *cubecen,         /* center coordinate of bounding cube */
   int procnum )                /* procnum of sending process */

{
   int    i,j,k;                /* loop counters */
   int    num_cells;            /* number of cells */
   int    cells_per_proc;       /* number f cells per processor */
   int    scell, ecell;         /* starting and ending cell indexes */
   int    cell_mask;            /* cell id mask */
   int    cell_x, cell_y, cell_z;
   int    cell_id;              /* cell identifier */
   int    proc_id;              /* processor id */
   int    wrap_mask;            /* mask to wrap cells around boundary */
   int    msg_type;             /* message type flag */
   int    orig_send_msg_buf;	/* to restore after sending messages */
   double cell_edge;            /* number of cells per edge */
   double pscale, poff;         /* offset and scaling for particle position */
#ifdef PIPED
   double v2xv3dotp, v3xv1dotp, v1xv2dotp; /* Dot products */
   double shiftc;               /* Proportional shift of center */
   PmtaVector v1xv2, v2xv3, v3xv1; /* Cross Products */
   PmtaVector tmpc;             /* Temporary Parallelepiped center */
#endif


   /*
    *  check to see if we have allocated space for the PVM
    *  message buffer id's.  also, create the message buffers.
    */

   if ( SendMsgBuf == (int *)NULL ) {
      SendMsgBuf = (int *)malloc(num_proc*sizeof(int));
      if ( SendMsgBuf == (int *)NULL ) {
	 fprintf(stderr,"ERROR: malloc() failed\n");
	 exit(-1);
      }
   }

   for (i=0; i<num_proc; i++ ) {
      SendMsgBuf[i] = pvm_mkbuf(DATA_NORMAL_PVM);
      if ( SendMsgBuf[i] < 0 ) {
	fprintf(stderr,"ERROR: pvm_mkbuk() failed\n");
	exit(-1);
      }
   }

   /*
   *  check to see if we have enough buffer space to hold the 
   *  particle and cell buffers and allocate extra if needed
   *
   */

   if ( num_parts > SendCellIdSz ) {
      SendCellId = (int *)realloc((void *)SendCellId,num_parts*sizeof(int));
      SendCellIdSz = num_parts;
      if ( SendCellId == (int *)NULL ) {
	 fprintf(stderr,"ERROR: realloc() failed\n");
	 exit(-1);
      }
   }

   /*
   *  allocate cells index -  note that this should never have
   *  to be reallocated until we get into some funky reallocation
   *  scheme.  but is always pays to be prepared.
   */

   num_cells = Power8[num_levels-1];
   if ( num_cells > SendPartCntSz ) {
      SendPartCnt = (int *)realloc((void *)SendPartCnt,num_cells*sizeof(int));
      if ( SendPartCnt == (int *)NULL ) {
	 fprintf(stderr,"ERROR: realloc() failed\n");
	 exit(-1);
      }
      SendPartCntSz = num_cells;
   }

   /* initailize counters */
   for (i=0; i<num_cells; i++) {
      SendPartCnt[i]=0;
   }

   /* compute cell numbers and size */

   cell_edge = (double)(0x1 << (num_levels-1));
   wrap_mask = (0x1 << (num_levels-1)) - 1;


   /* cycle through all particles */
   for (i=0; i<num_parts; i++) {

#ifndef PIPED
      /* compute integer cell coordinates */
      cell_x = (int)(((part_table[i].p.x - cubecen->x)/cubelen->x + 
               0.5) * cell_edge);
      cell_y = (int)(((part_table[i].p.y - cubecen->y)/cubelen->y + 
               0.5) * cell_edge);
      cell_z = (int)(((part_table[i].p.z - cubecen->z)/cubelen->z + 
               0.5) * cell_edge);

      /*
       * when running PBCs, it is a good thing if particles
       * do not get wrapped, but rather if they wander outside
       * the box, they are kept in the outter-most cell nearest
       * them.  otherwise it would look (to the simulation) that
       * a buch of particles were disappearing off of one side of
       * simulation space and appearing "magically" on the other.
       * this could play hell with constant energy sims.
       *
       * if the user wants wrapping, then simply define the PBC_WRAP
       * flag.
       */

#ifdef PBC_WRAP
      cell_x &= wrap_mask;
      cell_y &= wrap_mask;
      cell_z &= wrap_mask;
#else
      if ( cell_x < 0 )
	 cell_x = 0;
      if ( cell_x > wrap_mask )
	 cell_x = wrap_mask;
      if ( cell_y < 0 )
	 cell_y = 0;
      if ( cell_y > wrap_mask )
	 cell_y = wrap_mask;
      if ( cell_z < 0 )
	 cell_z = 0;
      if ( cell_z > wrap_mask )
	 cell_z = wrap_mask;
#endif

      /*
      *  build cell id from absolute coordinates
      *  there may be a more efficient way to do this
      */

      cell_mask = 0x1;
      cell_id = 0;
      cell_y = cell_y << 1;
      cell_z = cell_z << 2;
      for (j=0; j<num_levels-1; j++) {
         cell_id |= cell_x & cell_mask;
         cell_mask = cell_mask << 1;
         cell_x = cell_x << 2;
         cell_id |= cell_y & cell_mask;
         cell_mask = cell_mask << 1;
         cell_y = cell_y << 2;
         cell_id |= cell_z & cell_mask;
         cell_mask = cell_mask << 1;
         cell_z = cell_z << 2;
      } /* for j */

#else
      tmpc.x = cubecen->x;
      tmpc.y = cubecen->y;
      tmpc.z = cubecen->z;
      cell_id = 0;
      shiftc = 2.0;

      v2xv3.x = PVector2->y * PVector3->z - PVector2->z * PVector3->y;
      v2xv3.y = PVector2->z * PVector3->x - PVector2->x * PVector3->z;
      v2xv3.z = PVector2->x * PVector3->y - PVector2->y * PVector3->x;
 
      v3xv1.x = PVector3->y * PVector1->z - PVector3->z * PVector1->y;
      v3xv1.y = PVector3->z * PVector1->x - PVector3->x * PVector1->z;
      v3xv1.z = PVector3->x * PVector1->y - PVector3->y * PVector1->x;

      v1xv2.x = PVector1->y * PVector2->z - PVector1->z * PVector2->y;
      v1xv2.y = PVector1->z * PVector2->x - PVector1->x * PVector2->z;
      v1xv2.z = PVector1->x * PVector2->y - PVector1->y * PVector2->x;

      for(j=0;j<num_levels-1;j++)
      {
        cell_id = cell_id << 3;

        v2xv3dotp = ((part_table[i].p.x - tmpc.x) * v2xv3.x) +
                 ((part_table[i].p.y - tmpc.y) * v2xv3.y) +
                 ((part_table[i].p.z - tmpc.z) * v2xv3.z);
        v3xv1dotp = ((part_table[i].p.x - tmpc.x) * v3xv1.x) +
                 ((part_table[i].p.y - tmpc.y) * v3xv1.y) +
                 ((part_table[i].p.z - tmpc.z) * v3xv1.z);
        v1xv2dotp = ((part_table[i].p.x - tmpc.x) * v1xv2.x) +
                 ((part_table[i].p.y - tmpc.y) * v1xv2.y) +
                 ((part_table[i].p.z - tmpc.z) * v1xv2.z);

        if(v2xv3dotp<0) cell_x = 0;
        else cell_x = 1;

        if(v3xv1dotp<0) cell_y = 0;
        else cell_y = 1;

        if(v1xv2dotp<0) cell_z = 0;
        else cell_z = 1;

        cell_y = cell_y << 1;
        cell_z = cell_z << 2;

        cell_id |= (cell_x | cell_y | cell_z);

        shiftc *= 2.0;

        if(cell_x)
        {
          tmpc.x += PVector1->x/shiftc;
          tmpc.y += PVector1->y/shiftc;
          tmpc.z += PVector1->z/shiftc;
        }
        else
        {
          tmpc.x -= PVector1->x/shiftc;
          tmpc.y -= PVector1->y/shiftc;
          tmpc.z -= PVector1->z/shiftc;
        }

        if(cell_y)
        {
          tmpc.x += PVector2->x/shiftc;
          tmpc.y += PVector2->y/shiftc;
          tmpc.z += PVector2->z/shiftc;
        }
        else
        {
          tmpc.x -= PVector2->x/shiftc;
          tmpc.y -= PVector2->y/shiftc;
          tmpc.z -= PVector2->z/shiftc;
        }

        if(cell_z)
        {
          tmpc.x += PVector3->x/shiftc;
          tmpc.y += PVector3->y/shiftc;
          tmpc.z += PVector3->z/shiftc;
        }
        else
        {
          tmpc.x -= PVector3->x/shiftc;
          tmpc.y -= PVector3->y/shiftc;
          tmpc.z -= PVector3->z/shiftc;
        }
      }
#endif

      /* store the cell index for the particle */
      SendCellId[i] = cell_id;

      /* increment the partical counter for that cell */
      SendPartCnt[cell_id] += 1;

   } /* for i */


   /*
   *  cycle through all slave processes, construct the appropriate 
   *  particle tables for each cell, and package the whole thing 
   *  up to ship it.
   *
   *  the message format:
   *    0)   int     - message type
   *    1)   int     - sending processor number
   *    2)   int     - cell resize flag
   *    2.1) real[3] - length of cell edge (if needed)
   *    2.2) real[3] - position of cell center (if needed)
   *    3)   int     - the number of cells
   *    3.1) int     - the starting cell id
   *    4)   int[]   - number of particles per cell
   *    5)   int     - cell id for first cell
   *    6)   int     - particle id
   *    7)   real[]  - array of particle data
   *    8+)  misc    - repeat 2-4 for each particle
   *
   */

   orig_send_msg_buf = pvm_getsbuf();

   msg_type = MSG_PART1;

   for (i=0; i<num_proc; i++ ) {
      pvm_setsbuf(SendMsgBuf[i]);
      pvm_pkint(&msg_type,1,1);
      pvm_pkint(&procnum,1,1);

#ifndef EMBEDDED
      pvm_pkint(&resize,1,1);
      if (resize) {
#ifndef PIPED
	 pvm_pkdouble(&(cubelen->x),3,1);
#else
         pvm_pkdouble(&(PVector1->x),3,1);
         pvm_pkdouble(&(PVector2->x),3,1);
         pvm_pkdouble(&(PVector3->x),3,1);
#endif
	 pvm_pkdouble(&(cubecen->x),3,1);
      }
#endif

      scell = getscell( num_proc, i, num_levels-1 );
      ecell = getecell( num_proc, i, num_levels-1 );
      cells_per_proc = ecell - scell + 1;

      pvm_pkint(&(cells_per_proc),1,1);
      pvm_pkint(&(scell),1,1);
      /* pvm_pkint(&(SendPartCnt[scell]),cells_per_proc,1); */
      for (j=scell; j<=ecell; j++ ) {
        pvm_pkint(&(SendPartCnt[j]),1,1);
      }
      
   } /* for i */

   /*
   *  cycle through each cell for that processor
   *  big assumption here that the particle data only
   *  consists of four (or six) doubles and no other data
   *  in the structure
   *
   *  note we also assume that we have more cells than
   *  processors in the bottom level of the tree.
   */

   for (i=0; i<num_parts; i++) {
      cell_id = SendCellId[i];
      proc_id = getslvpid(num_proc, num_levels-1, cell_id );

      pvm_setsbuf(SendMsgBuf[proc_id]);

      pvm_pkint(&cell_id,1,1);
      pvm_pkint(&i,1,1);

#ifdef COMP_LJ
      pvm_pkdouble(&(part_table[i].p.x),6,1);
#else
      pvm_pkdouble(&(part_table[i].p.x),4,1);
#endif
   } /* for i */

   /*
    * terminate the messages and send them
    */

   for (i=0; i<num_proc; i++) {
      pvm_setsbuf(SendMsgBuf[i]);
      pvm_pkint(&(Msg_Term),1,1);
      pvm_send(tids[i],MSG_PART1);
   } /* for i */

   /*
    * free the send buffers after sending.
    */

   pvm_setsbuf(orig_send_msg_buf);

   for (i=0; i<num_proc; i++ ) {
      pvm_freebuf(SendMsgBuf[i]);
   } /* for i */

} /* end procedure */



/****************************************************************
*
*  Recv_Slave_Results() - receives completed data from slaves
*
*  we may wish to modify the message format to pass the information
*  back by cells.
*
*  message format:
*     0) int - pid of slave process
*     1) int - particle id number
*     2) array of force data for cell
*     6+) repeat 2-4 for each cell
*     4) message terminator
*
*  after reception of messages, we need to scale the results by the
*  length of the simulation space.
*
*/

Recv_Slave_Results(
   int num_parts,                /* number of particles to generate */
   PmtaPartInfoPtr fresults,     /* array of force results */
   PmtaPartInfoPtr fresults_lj,  /* array of LJ force results */
   int nprocs )                  /* number of slave processes */

   {
   int i,j,k;                /* loop counters */
   int pid, cellid, partid;  /* message reception counters */
   int numcells, numparts;   /* message unpacking counters */
   double fscale, pscale;    /* potential and force scaling */

#if defined VIRIAL || defined OLDVIRIAL
   double ftmp[4];           /* virial values from slaves */
#endif

   /*
   *  cycle through and receive particle message from each
   *  processor, placing the results in the receive buffer
   */

#if defined VIRIAL || defined OLDVIRIAL
   VirPot = 0.0;
   VirForce.x = 0.0;
   VirForce.y = 0.0;
   VirForce.z = 0.0;
#ifdef COMP_LJ
   VirPot_LJ = 0.0;
   VirForce_LJ.x = 0.0;
   VirForce_LJ.y = 0.0;
   VirForce_LJ.z = 0.0;
#endif
#endif

   for (i=0; i<nprocs; i++) {
      pvm_recv(-1, MSG_RESLT);

      pvm_upkint(&pid,1,1);


#if defined VIRIAL || defined OLDVIRIAL
      pvm_upkdouble(ftmp,4,1);
      VirPot += ftmp[0];
      VirForce.x += ftmp[1];
      VirForce.y += ftmp[2];
      VirForce.z += ftmp[3];
#ifdef COMP_LJ
      pvm_upkdouble(ftmp,4,1);
      VirPot_LJ += ftmp[0];
      VirForce_LJ.x += ftmp[1];
      VirForce_LJ.y += ftmp[2];
      VirForce_LJ.z += ftmp[3];
#endif
#endif

      pvm_upkint(&partid,1,1);
      while ( partid != -1 ) {
         pvm_upkdouble(&(fresults[partid].f.x),4,1);
#ifdef COMP_LJ
         pvm_upkdouble(&(fresults_lj[partid].f.x),4,1);
#endif
         pvm_upkint(&partid,1,1);
      } /* while partid */
   } /* for i */

} /* Recv_Slave_Results */




#if defined VIRIAL || defined OLDVIRIAL
/****************************************************************
 *
 *  Return_Virial() - returns the virial from the last particle
 *    reception.  This simple routine os provided to prevent direct
 *    access to the global virial variables, which could cause potential
 *    linking problems with the calling application.
 *
 */

int Return_Virial( double *vp, PmtaVector *vf,
   double *vp_lj, PmtaVector *vf_lj  ) {

   *vp = VirPot;
   vf->x = VirForce.x;
   vf->y = VirForce.y;
   vf->z = VirForce.z;

#ifdef COMP_LJ
   *vp_lj = VirPot_LJ;
   vf_lj->x = VirForce_LJ.x;
   vf_lj->y = VirForce_LJ.y;
   vf_lj->z = VirForce_LJ.z;
#endif

} /* Return_Virial() */

#endif


/****************************************************************
*
* Init_Local_Buffers() - initializes local datastructures
*
*/

void Init_Local_Buffers()
{

   SendPartCnt = NULL;
   SendCellId = NULL;
   RecvBuf = NULL;
   SendMsgBuf = NULL;
   
   SendPartCntSz = 0;
   SendCellIdSz = 0;
   RecvBufSz = 0;

} /* Init_Local_Buffers */
   

/****************************************************************
*
* Delete_Local_Buffers() - free up locally allocated dynamic data
*   structures.
*
*/

void Delete_Local_Buffers()
{

   free( SendPartCnt );
   free( SendCellId );
   free( RecvBuf );
   free( SendMsgBuf );

} /* Delete_Local_Buffers */
