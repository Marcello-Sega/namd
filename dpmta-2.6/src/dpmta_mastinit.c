/*
*  dpmta_mastinit - master initialization communications to slave processes
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*
*/

static char rcsid[] = "$Id: dpmta_mastinit.c,v 1.2 1997/09/12 22:56:30 jim Exp $";

/*
* revision history:
*
* $Log: dpmta_mastinit.c,v $
* Revision 1.2  1997/09/12 22:56:30  jim
* Modifications to work with converse pvm.
*
* Revision 1.1  1997/09/05 19:41:54  jim
* Original distribution.
*
 * Revision 2.8  1996/10/18  17:04:40  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.7  1996/09/24  18:41:49  wrankin
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
 * Revision 2.6  1996/08/09  15:30:40  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.5  1996/02/29  21:13:26  wrankin
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
 * Revision 2.4  1996/02/12  15:20:52  wrankin
 * Added Macroscopic Assemblies code.
 *
 * Revision 2.3  1995/11/29  22:29:06  wrankin
 * addition of periodic boundary code
 * beginning work on virial code.
 *
 * Revision 2.2  1995/10/01  21:45:50  wrankin
 * implemented LJ multipole calculations - compile time option only
 *
 * Revision 2.1  1995/05/16  04:02:16  wrankin
 * added code to allow resizing the sim cube on the fly
 * the calling process now scales the input to a unit cube
 * and scales the results back to the proper value
 * the slave assumes all data sent to it is in a unit cube
 * with the bottom corner at <0,0,0>
 *
 * Revision 1.5  1995/03/06  07:45:08  wrankin
 * master now sends Theta parameter to all slaves
 *
 * Revision 1.4  1994/11/30  18:02:49  wrankin
 * added processing to handle distributed calling sequence
 *
 * Revision 1.3  1994/11/28  04:27:53  wrankin
 * added hooks to handle multiple slave calls to PMTAforce
 *
 * Revision 1.2  1994/11/17  00:27:52  wrankin
 * added code to handled variable cell center
 *
 * Revision 1.1  1994/10/26  02:38:35  wrankin
 * Initial revision
 *
*
*/

/* include files */
#include <stdio.h>
#include "pvmc.h"
#include "dpmta.h"
#include "dpmta_pvm.h"


/****************************************************************
*
*  Send_Slave_Info() - send initialization message to slaves
*
*  message format:
*    1) int - slave pid number
*    2) int - number of slave processes
*    3) int[] - tids of all the slaves
*    4) int - number of levels in oct tree
*    5) int - fft flag
*    6) int - pbc flag
*    7) int - mp number
*    8) int - mp number for LJ (OPTIONAL)
*    9) dbl - theta, MAC separation parameter
*   10) int - number of macroscopic terms (OPTIONAL)
*   11) int - fft blocking factor
*   12) int - number of calling processes
*   13) int[] - tids of calling processes
*/


Send_Slave_Info(
   int nproc,
   int *tids,
   int num_levels,
   int fft,
   int pbc,
   int mp,
   int mp_lj,
   int macro_k,
   double theta,
   int fftblk,
   int callnum,
   int *calltids)
{

   int i;              /* loop counter */

   for (i=0; i<nproc; i++) {
      pvm_initsend(DATA_NORMAL_PVM);
      pvm_pkint(&i,1,1);
      pvm_pkint(&nproc,1,1);
      pvm_pkint(tids,nproc,1);
      pvm_pkint(&num_levels,1,1);
      pvm_pkint(&fft,1,1);
      pvm_pkint(&pbc,1,1);
      pvm_pkint(&mp,1,1);
#ifdef COMP_LJ
      pvm_pkint(&mp_lj,1,1);
#endif
      pvm_pkdouble(&theta,1,1);
#ifdef MACROSCOPIC
      pvm_pkint(&macro_k,1,1);
#endif
      pvm_pkint(&fftblk,1,1);
      pvm_pkint(&callnum,1,1);
      if ( callnum > 0 )
         pvm_pkint(calltids,callnum,1);

      pvm_send(tids[i],MSG_INIT1);
   } /* for i */
} /* Send_Slave_Info */
