/*
*  dpmta_slvcomm -  routines to send and receive multipole information
*     between slaves, send and receive particle information between 
*     slaves and send and receive interaction and particle
*     information to and from the master process.
*
*  w. t. rankin
*  w. elliott
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_slvcomm.c,v 1.1 1997/09/05 19:41:59 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvcomm.c,v $
 * Revision 1.1  1997/09/05 19:41:59  jim
 * Original distribution.
 *
 * Revision 2.35  1997/05/12  18:03:09  wrankin
 * fixed bug in force initialization
 * added malloc() checks
 *
 * Revision 2.34  1997/05/07  18:59:33  chumphre
 * major changes to implement ||-piped cells
 *
 * Revision 2.33  1997/04/28  16:05:54  wrankin
 * cleaned up handling of cell length and center
 *
 * Revision 2.32  1997/02/26  18:52:37  wrankin
 * cleaned up handling of cell resizing
 *
 * Revision 2.31  1997/01/13  22:12:32  wrankin
 * general cleanup for certain m4 definitions
 *
 * Revision 2.30  1996/11/18  19:29:32  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.29  1996/11/14  17:50:28  wrankin
 * performance enhancements to multipole M2L routine.
 * additions to make slaves exit gracefully.
 *
 * Revision 2.28  1996/11/01  02:26:44  wrankin
 * modifications to support cray t3d compilation
 * version update for 2.5 release
 *
 * Revision 2.27  1996/10/29  19:34:14  wrankin
 * changes to makefile to support other platforms
 * fix for multi-master communication interface
 * new global dumpdata routine to support direct pbc/macro verification
 *
 * Revision 2.26  1996/10/18  17:04:55  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.25  1996/09/24  18:42:16  wrankin
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
 * Revision 2.24  1996/09/10  18:07:01  wrankin
 * fixed macroscopic processing for non-cubic cells.
 * fixed passing of local expansions for non-2^n processors
 *
 * Revision 2.23  1996/08/20  17:12:41  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.22  1996/08/09  15:30:51  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.21  1996/02/29  21:13:36  wrankin
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
 * Revision 2.20  1996/02/12  15:21:03  wrankin
 * Added Macroscopic Assemblies code.
 *
 * Revision 2.19  1995/12/08  23:00:25  wrankin
 * preliminary release of DPMTA 2.3
 *   - added working Periodic Boundary Conditions (PDC) flag
 *   - added code for Virial Computation (not working yet)
 *   - general cleanup of some modules
 *
 * Revision 2.18  1995/11/29  22:29:21  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.17  1995/10/01  21:46:05  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.16  1995/09/18  19:02:41  wrankin
 * updates to move the scaling of the simulation box to the slave
 * process.
 *
 * Revision 2.15  1995/06/27  14:20:22  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.14  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 2.13  1995/04/09  22:21:06  wrankin
 * renaming of all global variables to prevent conflicts for
 * port of dpmta library to the T3D.
 *
 * Revision 2.12  1995/03/06  07:57:44  wrankin
 * slave process now receives Theta parameter from master
 * slave process now computes inv. interaction list locally
 *
 * Revision 2.11  1995/01/13  16:17:40  wrankin
 * removed use of 'msgtype' variable since it is redundant
 * updated routines that use global temp mpe to use correct structure
 * CMsum(F) procedure now called to sum multipoles
 * particle reception now realloc()'s particle array if not big enough
 * particle sending does not attempt to pack zero particles
 * .`
 *
 * Revision 2.10  1994/12/10  15:15:09  wrankin
 * forget to move message position of particle id field in slave receive
 * to correspond to changes made in previous version
 *
 * Revision 2.9  1994/12/06  19:13:22  wrankin
 * moved position if particle id within particle message
 *
 * Revision 2.8  1994/11/30  18:06:14  wrankin
 * added send and receive routines for distributed calling sequence
 * general cleanup.
 *
 * Revision 2.7  1994/11/28  04:27:53  wrankin
 * added hooks to handle multiple slave calls to PMTAforce
 *
 * Revision 2.6  1994/11/24  13:54:24  wrankin
 * updates to use static particle structure allocation
 * in preparation for having multiple slave entry points
 *
 * Revision 2.5  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 2.4  1994/10/19  00:05:55  wrankin
 * added fft blocking factor to initialization message
 * removed T3DBUG workaround.
 *
 * Revision 2.3  1994/10/14  05:02:28  wrankin
 * slowly incorporating all communication procedures into this file
 *
 * Revision 2.2  1994/08/25  19:01:22  wrankin
 * added T3DBUG message forwarding
 *
 * Revision 2.1  1994/06/02  12:52:27  wrankin
 * No change.  Update version number for Release 2
 *
 * Revision 1.9  1994/06/02  03:35:48  wrankin
 * general cleanup of code
 *
 * Revision 1.8  1994/05/27  06:35:53  wrankin
 * cleaned up debugging output
 *
 * Revision 1.7  1994/05/02  17:18:50  wrankin
 * replaced all pvm_[u]pkbyte with explicit type calls for CRAY-T3D
 *
 * Revision 1.6  1994/05/02  15:47:16  wrankin
 * change PvmDataRaw communications flag to PvmDataDefault to support CRAY-T3D
 *
 * Revision 1.5  1994/04/03  17:56:32  wrankin
 * initialized DownPassStart
 *
 * Revision 1.4  1994/03/31  13:32:58  wrankin
 * added new global CubeLength, passed in from master program
 *
 * Revision 1.3  1994/03/26  16:20:12  wrankin
 * added FFT and MP processing parameters
 *
 * Revision 1.2  1994/03/18  17:49:35  wrankin
 * added debugging statements
 *
 * Revision 1.1  1994/03/14  15:33:42  wrankin
 * Initial revision
 *
 *
*/

/* include files */
#include <stdio.h>
#include "pvm3.h"
#include "dpmta_pvm.h"
#include "dpmta_cell.h"
#include "dpmta_slave.h"

double Vec_Mag(Vector *);
double Max_CellLength();

/*
 * some communication globals
 */

static int Msg_Term = -1;


/****************************************************************
*
*  Recv_Master_Info() - receive initialization message from master
*
*  message format:
*    1) int - slave pid number
*    2) int - number of slave processes
*    3) int[] - tids of all the slaves
*    4) int - number of levels in oct tree
*    5) int - fft flag
*    6) int - mp number
*    7) dbl - theta, MAC separation parameter
*    8) int - fft blocking factor
*    9) int - slave interface parameter
*   10) int[] - tids of calling slaves (if any)
*/

Recv_Master_Info()
{

   pvm_recv(-1,MSG_INIT1);

   pvm_upkint(&Dpmta_Pid,1,1);
   pvm_upkint(&Dpmta_Nproc,1,1);
   pvm_upkint(Dpmta_Tids,Dpmta_Nproc,1);
   pvm_upkint(&Dpmta_NumLevels,1,1);
   pvm_upkint(&Dpmta_FFT,1,1);
   pvm_upkint(&Dpmta_PBC,1,1);
   pvm_upkint(&Dpmta_Mp,1,1);
#ifdef COMP_LJ
   pvm_upkint(&Dpmta_Mp_LJ,1,1);
#endif
   pvm_upkdouble(&Dpmta_Theta,1,1);
#ifdef MACROSCOPIC
   pvm_upkint(&Dpmta_K,1,1);
#endif
   pvm_upkint(&Dpmta_FftBlock,1,1);
   pvm_upkint(&Dpmta_CallingNum,1,1);
   if ( Dpmta_CallingNum > 0 ) {
      pvm_upkint(Dpmta_CallingTids,Dpmta_CallingNum,1);
      }

   if (Dpmta_PBC != 0)
      Dpmta_DownPassStart = 1;
   else
      Dpmta_DownPassStart = 2;

#ifndef EMBEDDED
   /*
    * force an initial resizing
    */

#ifndef PIPED
   Dpmta_CellLength.x = 1.0;
   Dpmta_CellLength.y = 1.0;
   Dpmta_CellLength.z = 1.0;
#else
   Dpmta_PVector1.x = 1.0;
   Dpmta_PVector1.y = 0.0;
   Dpmta_PVector1.z = 0.0;
   Dpmta_PVector2.x = 0.0;
   Dpmta_PVector2.y = 1.0;
   Dpmta_PVector2.z = 0.0;
   Dpmta_PVector3.x = 0.0;
   Dpmta_PVector3.y = 0.0;
   Dpmta_PVector3.z = 1.0;
   Dpmta_PV1Mag = 1.0;
   Dpmta_PV2Mag = 1.0;
   Dpmta_PV3Mag = 1.0;
#endif
   Dpmta_CellCenter.x = 0.0;
   Dpmta_CellCenter.y = 0.0;
   Dpmta_CellCenter.z = 0.0;

   /* can't set this to zero, since we divide by it later */
   Dpmta_MaxCellLen = 1.0;
#endif

} /* Recv_Master_Info() */



/****************************************************************
*
*  this procedure will pack and send the resulting force values
*  to their appropriate destination.  it check the value of the
*  CallingNum flag and calls the appropriate commounication/sorting
*  routine.
*
*  right now we use two seperate routines for message reception,
*  but in the future there is no reason that this couldn't be reduced.
*
*/


Send_Results()
{

   if ( Dpmta_CallingNum <= 1 ) {
      Send_Master_Results();
      }
   else {
      Send_Slave_Results();
      }

} /* Send_Results */


/****************************************************************
*
*  Send_Master_Results() - send completed data to parent
*
*  message format:
*    int pid - pid of the slave process
*    int num_parts - number of particle returned from pid
*    struct f_a[num_parts] - structure containing data of results
*
*  we may wish to modify the message format to pass the information
*  back by cells.
*
*  message format:
*     0) int - pid of slave process
*     (for each cell)
*     1) int - particle id number
*     2) array of force data for cell
*     3+) repeat  for each cell
*     4) int - termination value (-1)
*/

Send_Master_Results()
{
   int i,j;
   int num_parts;
   int level;
   int cell_id;

   level = Dpmta_NumLevels - 1;

   pvm_initsend(DATA_NORMAL_PVM);

   pvm_pkint(&Dpmta_Pid,1,1);

#if defined VIRIAL || defined OLDVIRIAL
   pvm_pkdouble(&Dpmta_Vpot,1,1);
   pvm_pkdouble((double *)(&Dpmta_Vf),3,1);
#ifdef COMP_LJ
   pvm_pkdouble(&Dpmta_Vpot_LJ,1,1);
   pvm_pkdouble((double *)(&Dpmta_Vf_LJ),3,1);
#endif
#endif

   for ( i=Dpmta_Scell[level]; i<=Dpmta_Ecell[level]; i++ ) {
      num_parts = Dpmta_CellTbl[level][i]->n;
      for ( j=0; j<num_parts; j++ ) {
         pvm_pkint(&(Dpmta_CellTbl[level][i]->mdata->part_id[j]),1,1);
         pvm_pkdouble(&(Dpmta_CellTbl[level][i]->mdata->flist[j].f.x),4,1);
#ifdef COMP_LJ
         pvm_pkdouble(&(Dpmta_CellTbl[level][i]->mdata->f_lj[j].f.x),4,1);
#endif
      }
   }

   /* terminate the message */
   pvm_pkint(&(Msg_Term),1,1);


   pvm_send(Dpmta_MasterTid,MSG_RESLT);

}



/****************************************************************
*
*  Send_Slave_Results() - send complete results to original slave
*    source.
*
*  message format:
*     0) int - pid of slave process
*     (for all particles in cell)
*     1) int - particle id number
*     2) array of force data for cell
*     3+) repeat for each cell
*     4) message terminator
*/

Send_Slave_Results()
{
   int destid;
   int i,j;
   int num_parts;
   int proc_num;
   int part_id;
   int level;
   int cell_id;


   level = Dpmta_NumLevels - 1;

   for ( destid=0; destid < Dpmta_CallingNum; destid++ ) {

      pvm_initsend(DATA_NORMAL_PVM);

      pvm_pkint(&Dpmta_Pid,1,1);

#if defined VIRIAL || defined OLDVIRIAL
      pvm_pkdouble(&Dpmta_Vpot,1,1);
      pvm_pkdouble((double *)(&Dpmta_Vf),3,1);
#ifdef COMP_LJ
      pvm_pkdouble(&Dpmta_Vpot_LJ,1,1);
      pvm_pkdouble((double *)(&Dpmta_Vf_LJ),3,1);
#endif
#endif

      for ( i=Dpmta_Scell[level]; i<=Dpmta_Ecell[level]; i++ ) {
         for ( j=0; j<(Dpmta_CellTbl[level][i]->n); j++ ) {
            proc_num = Dpmta_CellTbl[level][i]->mdata->proc_id[j];
            if ( proc_num == destid ) {
               pvm_pkint(&(Dpmta_CellTbl[level][i]->mdata->part_id[j]),1,1);
               pvm_pkdouble(&(Dpmta_CellTbl[level][i]->mdata->flist[j].f.x),
                  4,1);
#ifdef COMP_LJ
               pvm_pkdouble(&(Dpmta_CellTbl[level][i]->mdata->f_lj[j].f.x),
                  4,1);
#endif
            } /* if procnum */
         }  /* for j */
      } /* for i */

      pvm_pkint(&(Msg_Term),1,1);

      pvm_send(Dpmta_CallingTids[destid],MSG_RESLT);

   } /* for destid */

} /* Send_Slave_Results */



/****************************************************************
*
*  Send MPE Inverse Interaction List
* 
*  Send contents of mpe inverse interaction list to all other
*  processes. 
*
*  message format:
*
*    1) int - number of cells in message
*    2) int - size (in doubles) of the multipole expansion
*    3) int - cell id of first cell
*    4) double[] - array of mpe data, size determined from (2)
*   5+) repeat (3-4) for each cell in (1)
*
*  note that since the Complex data type is twice the size of
*  a double, we are packing MpeSize*2 doubles in each message.
*/
 
Slave_Send_Multipole()
   {
   int    rproc;             /* remote slave counter */
   int    i, j;              /* loop counters */
   int    cell_id;           /* cell id index */
   int    num_cells;         /* number of cells in iilist */
   double *mpe_ptr;          /* pointer to mpe data block */

   for (j=1; j<Dpmta_Nproc; j++) {
      rproc = ( Dpmta_Pid + j ) % Dpmta_Nproc;

      pvm_initsend(DATA_INPLACE_PVM);

      num_cells = Dpmta_IIlist.mlen[rproc];
      pvm_pkint(&(Dpmta_IIlist.mlen[rproc]),1,1);

      for ( i=0; i<num_cells; i++ ) {
         cell_id = Dpmta_IIlist.mlist[rproc][i];
         pvm_pkint(&(Dpmta_IIlist.mlist[rproc][i]),1,1);
         pvm_pkdouble(&(Dpmta_CellTbl[0][cell_id]->m[0][0].x),
		      Dpmta_MpeSize*2,1);
#ifdef COMP_LJ
         pvm_pkdouble(&(Dpmta_CellTbl[0][cell_id]->m_lj[0][0][0].x),
		      Dpmta_MpeSize_LJ*2,1);
#endif
      } /* for i */

      pvm_send(Dpmta_Tids[rproc], MSG_MDIST);

   } /* for rproc */
} /* Slave_Send_Multipole */


/****************************************************************
*
*  Receive MPE Interaction List
*  Receive contents of mpe interaction list from all children
*
*  note that msize tells us the number of Complex data in the message
*  so that we have to unpack twice the number of doubles.
*/

Slave_Recv_Multipole()
   {
   int    sproc;             /* remote slave counter */
   int    i, j;              /* loop counters */
   int    cell_id;           /* cell id index */
   int    num_cells;         /* number of cells in iilist */
   double *mpe_ptr;          /* pointer to mpe data block */


   for ( i=1; i<Dpmta_Nproc; i++ ) {

      pvm_recv(-1, MSG_MDIST);

      pvm_upkint(&num_cells,1,1);

      for ( j=0; j<num_cells; j++ ) {
         pvm_upkint(&cell_id,1,1);

         mpe_ptr = &(Dpmta_CellTbl[0][cell_id]->m[0][0].x);
         pvm_upkdouble(mpe_ptr,Dpmta_MpeSize*2,1);

#ifdef COMP_LJ
         mpe_ptr = &(Dpmta_CellTbl[0][cell_id]->m_lj[0][0][0].x);
         pvm_upkdouble(mpe_ptr,Dpmta_MpeSize_LJ*2,1);
#endif

/* patch - need a better implementation later - wtrankin */
         Dpmta_CellTbl[0][cell_id]->mvalid = TRUE;

      } /* for j */
   } /* for i */
} /* Slave_Recv_Multipole */


/****************************************************************
*
*  Send_Mpe_to_Parent(level) - Send MPE to Pid
* 
*  sends the all mpe for this level to the pids that own them.
*
*  Requires following global variables:
*  MpeSize -- number of Complex in a multipole expansion (fcn of p and
*             whether or not FFT is used)
*   
*/

void Send_Mpe_to_Parent( int level )
{
   int i,j,k;
   int cellid;
   int tpid;
   double *mpe_ptr;

   /* we do not need to send MPEs too far up the tree */
   if ( level < Dpmta_DownPassStart ) {
      return;
   }

   k = -1;
   for (i=Dpmta_Scell[level]; i<=Dpmta_Ecell[level]; i++) {
      j = getparent(i);
      if ( j != k ) {
         k = j;
         tpid = getslvpid(level-1,j);
         if ( tpid != Dpmta_Pid ) {

#ifdef CRAY
            pvm_initsend(DATA_NORMAL_PVM);
#else
            pvm_initsend(DATA_INPLACE_PVM);
#endif
            /* compute old-fashioned cellid */
            cellid = j + Dpmta_LevelLocate[level-1];
            pvm_pkint(&cellid,1,1);

            mpe_ptr = &(Dpmta_CellTbl[0][cellid]->m[0][0].x);
            pvm_pkdouble(mpe_ptr,Dpmta_MpeSize*2,1);

#ifdef COMP_LJ
            mpe_ptr = &(Dpmta_CellTbl[0][cellid]->m_lj[0][0][0].x);
            pvm_pkdouble(mpe_ptr,Dpmta_MpeSize_LJ*2,1);
#endif
            pvm_send(Dpmta_Tids[tpid], MSG_M2P);

	 } /* if tpid */
      } /* if (j!=k) */
   } /* for i */

} /* Send_Mpe_to_Parent */


/****************************************************************
*
*  Receive MPE from Child
* 
*  pid with parental responsibilities receives partial shifted
*  mpe from another processor(s) and adds it to the current mpe.
*  Same setup as send MPE.
*/

void Rcv_Mpe_from_Child( int level )
{
   int i,j,k;
   int cellid_tmp;
   double *mpe_ptr, *mpe_temp;

   for (i=0; i<Dpmta_RMcell[level]; i++) {

      /*
      *  receive a message and place the mpe in temporary (global)
      *  storage.  then add the mpe data to the specified cell's mpe.
      */

      pvm_recv(-1, MSG_M2P);
      pvm_upkint(&cellid_tmp,1,1);

      mpe_temp = &(Dpmta_Temp_Mpe[0][0].x);
      pvm_upkdouble(mpe_temp,Dpmta_MpeSize*2,1);

#ifdef COMP_LJ
      mpe_temp = &(Dpmta_Temp_Mpe_LJ[0][0][0].x);
      pvm_upkdouble(mpe_temp,Dpmta_MpeSize_LJ*2,1);
#endif

      /*
      *  add child multipole into existing parent
      */

      if (Dpmta_FFT)
         CMsumF(Dpmta_Temp_Mpe,Dpmta_CellTbl[0][cellid_tmp]->m,Dpmta_Mp);
      else
         CMsum(Dpmta_Temp_Mpe,Dpmta_CellTbl[0][cellid_tmp]->m,Dpmta_Mp);

#ifdef COMP_LJ
      LJMsum(Dpmta_Temp_Mpe_LJ,Dpmta_CellTbl[0][cellid_tmp]->m_lj,Dpmta_Mp_LJ);
#endif

   } /* for i */

} /* Rcv_Mpe_from_Child */

    
/****************************************************************
*
*  Send_Lcl_to_Child - Send Local Expansion of Parent to all children
*  
*  Pid without parental responsibilities receives local expansion.
*  (Child will shift the expansion).
*
*  so, in any case, all we have to do is to send the local expansion
*  from a single cell down to each of this slaves children
*
*  note also that we do not use a broadcast for messages to multiple 
*  cells, which we shold probably do.
*/

void Send_Lcl_to_Child( int level )
{
   int i, j, k;
   int tpid;
   int ccell, cell_id;
   double *lcl_ptr;

   if (( Dpmta_Scell[level] != -1 ) && ( level < Dpmta_NumLevels-1 )) {
      for ( i=Dpmta_Scell[level]; i<=Dpmta_Ecell[level]; i++ ) {

         /*
         *  for each slave, send a message containing the local expansion
         */

         cell_id = i + Dpmta_LevelLocate[level];
         ccell = getfirstchild(i);

         tpid = -1;

         for ( j=ccell; j<ccell+8; j++ ) {

	    k = getslvpid(level+1,j);

            if ( k != Dpmta_Pid ) {
 	       if ( k != tpid ) {
                  tpid = k;

#ifdef CRAY
                  pvm_initsend(DATA_NORMAL_PVM);
#else
                  pvm_initsend(DATA_INPLACE_PVM);
#endif

                  pvm_pkint(&cell_id,1,1);

                  lcl_ptr = &(Dpmta_CellTbl[level][i]->mdata->l[0][0].x);
                  pvm_pkdouble(lcl_ptr,Dpmta_LclSize*2,1);

#ifdef COMP_LJ
                  lcl_ptr = &(Dpmta_CellTbl[level][i]->mdata->l_lj[0][0][0].x);
                  pvm_pkdouble(lcl_ptr,Dpmta_LclSize_LJ*2,1);
#endif

                  pvm_send(Dpmta_Tids[tpid],MSG_L2C);

	       } /* if k */
	    } /* if dpmta_pid */
         } /* for j */
      } /* for i */
   } /* in Dpmta_Scell */

} /* Send_Lcl_to_Child */
    

/****************************************************************
* 
*  Rcv_Lcl_from_Parent - Receive local expansion(s) from parent(s)
* 
*  Pid with parental responsibilities sends local expansion of parent
*  to pids holding children.  The number of local expansion that we
*  will receive has been precomputed and stored in Dpmta_RLcell[].
*
*/

void Rcv_Lcl_from_Parent( int level )
{
   int i,j,k;
   int cell_id;
   double *lcl_ptr;

   for ( i=0; i<Dpmta_RLcell[level]; i++ ) {

      pvm_recv(-1,MSG_L2C);
      pvm_upkint(&cell_id,1,1);

      lcl_ptr = &(Dpmta_CellTbl[0][cell_id]->mdata->l[0][0].x);
      pvm_upkdouble(lcl_ptr,Dpmta_LclSize*2,1);

#ifdef COMP_LJ
      lcl_ptr = &(Dpmta_CellTbl[0][cell_id]->mdata->l_lj[0][0][0].x);
      pvm_upkdouble(lcl_ptr,Dpmta_LclSize_LJ*2,1);
#endif

      Dpmta_CellTbl[0][cell_id]->mdata->lvalid = TRUE;

   } /* for i */
} /* Rcv_Lcl_from_Parent */


/****************************************************************
*
*  this procedure will receive the particle information from the
*  other DPMTA slave processes and load the values into it's cell
*  table.  these particles will be used to compute the single direct
*  particle interactions.
*
*  the slave will (re)allocate additional space for the particle list
*  if it is needed and update the contents of the slaves cell table.
*
*  the message format:
*    1) int - number of cells
*    2) int - cell id for first cell
*    3) int - number of particles in cell
*    4) int[] - array of particle data for cell
*   5+) repeat 2-4 for each cell
*
*/

Slave_Recv_SDirect()
   {

   int i,j,k;                   /* loop counters */

   int num_cells;               /* number of cells in message */
   int cell_id;                 /* cell identifier */
   int num_parts;               /* message particle number */

   ParticlePtr part_tmp;        /* temporary pointer to list of particles */


   for (i=1; i<Dpmta_Nproc; i++) {
      pvm_recv(-1, MSG_PDIST);

      pvm_upkint(&num_cells,1,1);
      for (j=0; j<num_cells; j++ ) {
         pvm_upkint(&cell_id,1,1);

         /*
         *  check to see if we have allocated a cell here and bail out
         *  if we have not.  later, we should replace this with a call
         *  that will create a cell if one is needed.
         */

         if (Dpmta_CellTbl[0][cell_id] == NULL ) {
            fprintf(stderr,"Error: Slave_Recv_Part() - cell %d not alloc\n",
               cell_id);
            exit(-1);
            }

         /*
         *  here is where we need to check psize against the number 
         *  of particles in the list and realloc is it is not large
         *  enough.  since we are receiving all of out particles from
         *  a single source (the slave that owns this cell) we don't
         *  have to worry about preserving any existing data.
         */

         pvm_upkint(&num_parts,1,1);
         Dpmta_CellTbl[0][cell_id]->n = num_parts;

         if ( num_parts > Dpmta_CellTbl[0][cell_id]->psize ) {
            part_tmp = Dpmta_CellTbl[0][cell_id]->plist;
            part_tmp = (ParticlePtr)realloc(part_tmp,
                          num_parts*sizeof(Particle));
            Dpmta_CellTbl[0][cell_id]->plist = part_tmp;
            Dpmta_CellTbl[0][cell_id]->psize = num_parts;
            }

         if ( num_parts > 0 ) {
            part_tmp = &(Dpmta_CellTbl[0][cell_id]->plist[0]);
#ifdef COMP_LJ
            pvm_upkdouble((double *)part_tmp,num_parts*6,1);
#else
            pvm_upkdouble((double *)part_tmp,num_parts*4,1);
#endif
            }

	 } /* for j */
      } /* for i */
   } /* Slave_Recv_Sdirect */



/****************************************************************
*
*  this procedure will cycle through the other processors and
*  pack/send a message to each one that contains the particle 
*  lists needed by the target processor, as defined in the inverse 
*  interaction list.
*
*  message format:
*
*    1) int - number of cells
*    2) int - cell id of first cell
*    3) int - number of particles in first cell
*    4) particle[] - array of particle data, size determined from (2)
*   5+) repeat (2-4) for each cell in (1)
*
*/

Slave_Send_SDirect()
{

   int i, j, k;                 /* loop counters, what else? */
   int cell_id;                 /* temp var used to hold cell id */
   int num_cells;               /* temp ver to hold number of cells */
   int num_parts;               /* temp ver to hold particle number */
   int size_part;               /* computed length of the particle array */
   int proc_id;                 /* id of where to send message */

   ParticlePtr part_tmp;        /* temp pointer to particles */


   for (i=1; i<Dpmta_Nproc; i++) {
      proc_id = ( Dpmta_Pid + i ) % Dpmta_Nproc;

      pvm_initsend(DATA_INPLACE_PVM);

      num_cells = Dpmta_IIlist.dlen[proc_id];
      pvm_pkint(&(Dpmta_IIlist.dlen[proc_id]),1,1);

      /* cycle through each cell for that processor */
      for (j=0; j<num_cells; j++) {

         cell_id = Dpmta_IIlist.dlist[proc_id][j];
         pvm_pkint(&(Dpmta_IIlist.dlist[proc_id][j]),1,1);

         num_parts = Dpmta_CellTbl[0][cell_id]->n;         
         pvm_pkint(&(Dpmta_CellTbl[0][cell_id]->n),1,1);

         if ( num_parts > 0 ) {
#ifdef COMP_LJ
            pvm_pkdouble((double *)(Dpmta_CellTbl[0][cell_id]->plist),
			 num_parts*6,1);
#else
            pvm_pkdouble((double *)(Dpmta_CellTbl[0][cell_id]->plist),
			 num_parts*4,1);
#endif
         }

      } /* for j */

      pvm_send(Dpmta_Tids[proc_id],MSG_PDIST);

   } /* for i */

} /* Slave_Send_SDirect */



/****************************************************************
*
*  Recv_Particles() - receive particles from application
*
*  this procedure will receive the particle information from the
*  one or more slave processes and load the values into its cell
*   table.
*
*  all the particle data for all of the cells owned by
*  the slave processor is sent as a single packet (at least for
*  right now).
*
*  the message format:
*    0)   int     - message type
*    1)   int     - sending processor number
*    2)   int     - cell resize flag
*    2.1) real[3] - length of cell edge (if needed) or
*    2.1) real[9] - 3 Vectors defining ||-piped (if needed) 
*    2.2) real[3] - position of cell center (if needed)
*    3)   int     - the number of cells
*    3.1) int     - the starting cell id
*    4)   int[]   - number of particles per cell
*    5)   int     - cell id for first cell
*    6)   int     - particle id
*    7)   real[]  - array of particle data
*    8+)  misc    - repeat 2-4 for each particle
*
*  repeat for each sending slave
*
*  important fact: we store the process id of the sending process in the
*  upper bit positions of the particle ID.  the constant PID_SHIFT defines
*  the starting bit position of this field.  this limits the number of 
*  paticles that may be received from a single proc, as well as the largest
*  possible number of processors.
*
*/

Recv_Particles()
{

   int i,j,k;                   /* loop counters */

   int num_cells;               /* number of cells in message */
   int scell;                   /* starting cell identifier */
   int cell_id;                 /* cell identifier */
   int num_parts;               /* message particle number */
   int src_proc;                /* processor id of sending proc */
   int level;                   /* bottom level of tree */
   int msg_type;                /* flag for termination message */

   Vector cubelen;              /* length of cell edge */
   Vector cubecen;              /* position of cell center */

   CellPtr cell_tmp;            /* temp cell table pointer */
   ParticlePtr part_tmp;        /* temporary pointer to list of particles */
   PartInfoPtr flist_tmp;       /* temporary pointer to force vector */
   int *idlist_tmp;             /* temporary pointer to id vector */
   int *pidlist_tmp;            /* temporary pointer to pid vector */
   double *dbl_tmp;             /* temp pointer to particle data */


   level = Dpmta_NumLevels - 1;

   /*
   *  receive particle lists and check cell table
   *
   *  receives are ordered to ensure constant ordering of
   *  particles within cells, as per UIUC request.  we should
   *  probaly stagger both the sends and receives so as to minimize
   *  the waiting time.
   */

   for ( i=0; i<Dpmta_CallingNum; i++ ) {

      pvm_recv(Dpmta_CallingTids[i],MSG_PART1);

      /*
       *  determine if this is a termination message
       *  we should really call some sort of cleanup routine
       *  here to free up the data structures.
       */

      pvm_upkint(&msg_type,1,1);
      if ( msg_type == MSG_EXIT ) {
         pvm_exit();
         exit(0);
      }

      pvm_upkint(&src_proc,1,1);

#ifndef EMBEDDED
      /*
       * unpack the unit cell size and check if we need
       * to resize.  if we are an embedded application (T3E),
       * then we do not need to pass this information.  we
       * do not need to perform the resize call until after
       * all the processors have replied.
       */

      pvm_upkint(&(Dpmta_Resize),1,1);
      if ( Dpmta_Resize ) {
#ifndef PIPED
	 pvm_upkdouble(&(Dpmta_CellLength.x),3,1);
#else
         pvm_upkdouble(&(Dpmta_PVector1.x),3,1);
         pvm_upkdouble(&(Dpmta_PVector2.x),3,1);
         pvm_upkdouble(&(Dpmta_PVector3.x),3,1);
#endif
	 pvm_upkdouble(&(Dpmta_CellCenter.x),3,1);

#ifdef PIPED
         Dpmta_PV1Mag = Vec_Mag(&Dpmta_PVector1);
         Dpmta_PV2Mag = Vec_Mag(&Dpmta_PVector2);
         Dpmta_PV3Mag = Vec_Mag(&Dpmta_PVector3);
#endif
         Dpmta_MaxCellLen = Max_CellLength();
      }
#endif

      /*
       * unpack the cell size counters and determine if each
       * cell particle array is large enough to hold everything.
       */

      pvm_upkint(&num_cells,1,1);
      pvm_upkint(&scell,1,1);

      cell_id = scell;
      for ( j=0; j<num_cells; j++ ) {

	 /*
	  *  check to see if we have allocated a cell here and bail out
	  *  if we have not.  later, we should replace this with a call
	  *  that will create a cell if one is needed.
	  */

	 cell_tmp = Dpmta_CellTbl[level][cell_id];

	 if ( cell_tmp == NULL ) {
	    fprintf(stderr,"proc[%d] - cell[%d] not alloc\n",Dpmta_Pid,cell_id);
	    exit(-1);
	 }
	 if ( cell_tmp->mdata == NULL ) {
	    fprintf(stderr,"proc[%d] - cell[%d]->mdata not alloc\n",Dpmta_Pid,cell_id);
	    exit(-1);
	 }


	 pvm_upkint(&num_parts,1,1);

         /* keep track of old number of parts for indexing */
         num_parts += cell_tmp->n;

         if ( num_parts > cell_tmp->psize ) {

            part_tmp = cell_tmp->plist;
            part_tmp = (ParticlePtr)realloc(part_tmp, (num_parts)*sizeof(Particle));
	    if ( part_tmp == NULL ) {
	       fprintf(stderr,"ERROR: Recv_Particles() - alloc failed\n");
	       exit(-1);
	    }
            cell_tmp->plist = part_tmp;

            idlist_tmp = cell_tmp->mdata->part_id;
            idlist_tmp = (int *)realloc(idlist_tmp, (num_parts)*sizeof(int));
	    if ( idlist_tmp == NULL ) {
	       fprintf(stderr,"ERROR: Recv_Particles() - alloc failed\n");
	       exit(-1);
	    }
            cell_tmp->mdata->part_id = idlist_tmp;

            pidlist_tmp = cell_tmp->mdata->proc_id;
            pidlist_tmp = (int *)realloc(pidlist_tmp, (num_parts)*sizeof(int));
	    if ( pidlist_tmp == NULL ) {
	       fprintf(stderr,"ERROR: Recv_Particles() - alloc failed\n");
	       exit(-1);
	    }
            cell_tmp->mdata->proc_id = pidlist_tmp;

            flist_tmp = cell_tmp->mdata->flist;
            flist_tmp = (PartInfoPtr)realloc(flist_tmp,
					     (num_parts)*sizeof(PartInfo));
	    if ( flist_tmp == NULL ) {
	       fprintf(stderr,"ERROR: Recv_Particles() - alloc failed\n");
	       exit(-1);
	    }
            cell_tmp->mdata->flist = flist_tmp;

#ifdef COMP_LJ
            flist_tmp = cell_tmp->mdata->f_lj;
            flist_tmp = (PartInfoPtr)realloc(flist_tmp,
					     (num_parts)*sizeof(PartInfo));
	    if ( flist_tmp == NULL ) {
	       fprintf(stderr,"ERROR: Recv_Particles() - alloc failed\n");
	       exit(-1);
	    }
            cell_tmp->mdata->f_lj = flist_tmp;
#endif

            cell_tmp->psize = num_parts;
         } /* if num_parts */

	 /* access the next cell */
	 cell_id++;

      } /* for j */

      /*
       * now unpack eack particle one at a time
       * into the appropriate cell
       */

      pvm_upkint(&cell_id,1,1);
      
      while ( cell_id >= 0 ) {

	 cell_tmp = Dpmta_CellTbl[level][cell_id];

	 num_parts = cell_tmp->n;

	 dbl_tmp = &(cell_tmp->plist[num_parts].p.x);
	 idlist_tmp = &(cell_tmp->mdata->part_id[num_parts]);
	 pidlist_tmp = &(cell_tmp->mdata->proc_id[num_parts]);

	 /* unpack data into arrays */
	 pvm_upkint(idlist_tmp,1,1);
	 *pidlist_tmp = src_proc;
#ifdef COMP_LJ
	 pvm_upkdouble(dbl_tmp,6,1);
#else
	 pvm_upkdouble(dbl_tmp,4,1);
#endif

	 cell_tmp->n += 1;

	 pvm_upkint(&cell_id,1,1);

      } /* while cell_id */

   } /* for i */

   /*
    *  now that we have reached the end of the receive loop
    *  we need to initialize the various data arrays and
    *  set the cell size if we did a resize
    */

   cell_id = scell;

   for ( i=0; i<num_cells; i++ ) {

      cell_tmp = Dpmta_CellTbl[level][cell_id];

      num_parts = cell_tmp->n;

      /* zero out force array */
      flist_tmp = cell_tmp->mdata->flist;
      for ( j=0; j<num_parts; j++ ) {
	 flist_tmp->f.x = 0.0;
	 flist_tmp->f.y = 0.0;
	 flist_tmp->f.z = 0.0;
	 flist_tmp->v = 0.0;
	 flist_tmp++;
      } /* for j */

#ifdef COMP_LJ
      flist_tmp = cell_tmp->mdata->f_lj;
      for ( j=0; j<num_parts; j++ ) {
	 flist_tmp->f.x = 0.0;
	 flist_tmp->f.y = 0.0;
	 flist_tmp->f.z = 0.0;
	 flist_tmp->v = 0.0;
	 flist_tmp++;
      } /* for j */
#endif

      /* initialize next cell */
      cell_id++;

   } /* for i */

} /* Recv_Slave_Particles */
