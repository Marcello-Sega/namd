/* 
*  dpmta_slvcompute.c - routine to perform actuall DPMTA processing
*
*  w. t. rankin, w. elliott
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*  this routine implements the actaul DPMTA processing.
*  it was pulled out of the original dpmta_slave.c main() routine
*  for the purpose of t3d integration.
*
*/

static char rcsid[] = "$Id: dpmta_slvcompute.c,v 1.1 1997/09/05 19:41:59 jim Exp $";

/*
 * revision history:
 *
 * $Log: dpmta_slvcompute.c,v $
 * Revision 1.1  1997/09/05 19:41:59  jim
 * Original distribution.
 *
 * Revision 2.12  1997/03/04  19:18:07  wrankin
 * updates to timing codes
 *
 * Revision 2.11  1997/02/26  20:43:31  wrankin
 * updated timing measurements and placed routines in single module
 *
 * Revision 2.10  1997/02/26  16:54:28  wrankin
 * added support for serial non-pvm compilation of DPMTA
 *
 * Revision 2.9  1997/02/24  19:12:28  wrankin
 * fixes to timing code
 *
 * Revision 2.8  1997/01/28  17:10:45  wrankin
 * added new timing codes for macroscopic expansion calculations
 *
 * Revision 2.7  1996/10/18  17:04:59  wrankin
 * improved PVM message passing structures.
 * made current and cleaned up T3D library code
 * added additional test codes
 * implement new heirarchical make structure
 *
 * Revision 2.6  1996/09/24  18:42:38  wrankin
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
 * Revision 2.5  1996/08/20  17:12:43  wrankin
 * Enhancement to allow any number of processes to be utilized.
 * Formally the number of precesses was limited to a power of 2.
 *
 * Revision 2.4  1995/07/17  01:13:35  wrankin
 * updates to support SGI Power Challenge (IRIX 6.0.1)
 * cleanup of Makefile
 * initial work on T3D port (not yet supported in this release)
 *
 * Revision 2.3  1995/07/10  02:47:14  wrankin
 * multipole processing code modified to use precomputed mpe transfer
 *   matrices.
 *
 * in addition, all multippole calculation routines have been removed
 *   from this source distribution and will be placed in a separate
 *   multipole library to be supplied externally.
 *
 * Revision 2.2  1995/06/27  14:20:24  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.1  1995/06/13  04:26:10  wrankin
 * No updates.
 * Sources placed under CVS.
 *
 * Revision 1.1  1995/04/24  04:13:35  wrankin
 * Initial revision
 *
 *
*/

/* include files */

#include <stdio.h>
#include <unistd.h>
#include "dpmta_cell.h"          /* data type definitions */
#include "dpmta_slave.h"         /* global variable declarations */

#ifdef TIME
#include "dpmta_timer.h"
#endif


/****************************************************************
*
*  Slave_Compute() - perform all multipole computations
*
*  this routine performs all multipole comutations on the
*  particle data.
*
*/

Slave_Compute()
{
      /*
      *  based on inverse interaction lists, distribute particle
      *  cell information to all other processors
      */

#ifndef SERIAL
      Slave_Send_SDirect();
#endif

      /*
      *  perform the upward pass
      *
      *  each slave will process the multipole expansion for the cells
      *  that it owns.  if the higher level cell is owned by another slave,
      *  (based upon the slaved pid) then this slave will pass the multipole
      *  expansions to the appropriate parent.
      */

#ifdef TIME
      /* begin timing step 1 and 2 */
      times(&startbuf);
#endif

      Slave_Mpole_Exp();

#ifdef TIME
      times(&endbuf);
      times_arr[9] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif


      /*
      *  collect the particle information from the other processors
      *  and place the information in the local cell table.  note that
      *  we only receive nproc-1 messages.
      */

#ifndef SERIAL
      Slave_Recv_SDirect();
#endif

      /*
      *  based on the inverse interaction list, distribute mpe
      *  cell information to all other processors
      */

#ifndef SERIAL
      Slave_Send_Multipole();
#endif

#ifdef TIME
      gettimeofday(&runstruct,0);
      times_arr[2] = (double)runstruct.tv_sec +
	 ((double)runstruct.tv_usec/1000000);
      times_arr[2] -= times_arr[0];
#endif


      /*
      *  since we have all the local and remote particle information,
      *  do all the direct particle interaction calculations
      *  while we are waiting for the multipole data to 
      *  arrive.
      */

#ifdef TIME
      /* begin timing direct calculations */
      times(&startbuf);
#endif

	Slave_Direct_Calc();

#ifdef TIME
      /* end time for direct calculations */
      times(&endbuf);
      times_arr[10] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif


      /*
      *  collect the mpe data from the other processors
      */

#ifndef SERIAL
      Slave_Recv_Multipole();
#endif

#ifdef TIME
      gettimeofday(&runstruct,0);
      times_arr[3] = (double)runstruct.tv_sec +
	 ((double)runstruct.tv_usec/1000000);
      times_arr[3] -= times_arr[0];
#endif


      /*
      *  compute the mpe-particle interactions.  this is also known as
      *  the downward mpe pass.  after all the local expansions have
      *  been computed, then compute the particle forces based upon
      *  the local expansion.
      */

#ifdef TIME
      /* Time downward pass */
      times(&startbuf);
#endif

      Slave_MPE_Calc();

#ifdef TIME
      /* end time for direct calculations */
      times(&endbuf);
      times_arr[11] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif

#ifdef TIME
      /* Time downward pass */
      times(&startbuf);
#endif

      Slave_MPE_Force();

#ifdef TIME
      /* end time for direct calculations */
      times(&endbuf);
      times_arr[12] = (double)(endbuf.tms_utime - startbuf.tms_utime) / (double)CLK_TCK;
#endif

#ifdef TIME
      gettimeofday(&runstruct,0);
      times_arr[4] = (double)runstruct.tv_sec +
	 ((double)runstruct.tv_usec/1000000);
      times_arr[4] -= times_arr[0];
#endif


} /* Slave_Compute() */
