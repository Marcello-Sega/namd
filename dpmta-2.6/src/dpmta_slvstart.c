/*
*  dpmta_slvstart.c - procedures to perform any initialization needed
*     at the beginning of a DPMTA iteration
*
*  w. t. rankin
*
*  Copyright (c) 1995 Duke University
*  All rights reserved
*
*/

static char rcsid[] = "$Id: dpmta_slvstart.c,v 1.1 1997/09/05 19:42:08 jim Exp $";

/*
 * Revision history:
 *
 * $Log: dpmta_slvstart.c,v $
 * Revision 1.1  1997/09/05 19:42:08  jim
 * Original distribution.
 *
 * Revision 2.4  1996/11/18  19:29:40  wrankin
 * update to virial code to use multipole virial calculation
 * old virial code is turned on by defining OLDVIRIAL
 * direct test program now computes absolute virial tensor.
 *
 * Revision 2.3  1996/09/24  18:43:49  wrankin
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
 * Revision 2.2  1996/08/09  15:31:12  chumphre
 * Implemented noncubic (rectangular) boxes.
 *
 * Revision 2.1  1995/12/08  23:00:41  wrankin
 * preliminary release of DPMTA 2.3
 *   - added working Periodic Boundary Conditions (PDC) flag
 *   - added code for Virial Computation (not working yet)
 *   - general cleanup of some modules
 *
 *
*/

/* include files */

#include <stdio.h>
#include "dpmta_cell.h"
#include "dpmta_slave.h"


/****************************************************************
*
*  Slave_Start() - perform any needed initiializations
*
*/

Slave_Start()
{

#if defined VIRIAL || defined OLDVIRIAL
   Dpmta_Vpot = 0.0;
   Dpmta_Vf.x = 0.0;
   Dpmta_Vf.y = 0.0;
   Dpmta_Vf.z = 0.0;
#ifdef COMP_LJ
   Dpmta_Vpot_LJ = 0.0;
   Dpmta_Vf_LJ.x = 0.0;
   Dpmta_Vf_LJ.y = 0.0;
   Dpmta_Vf_LJ.z = 0.0;
#endif
#endif

} /* Slave_Start */

