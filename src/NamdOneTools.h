/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:	Extern's defined for NAMD 1.X functions used in odd places
 *              Mainly found in WorkDistrib.
 ***************************************************************************/

#include "common.h"
#include "NamdTypes.h"
#include "Vector.h"
#include "PDB.h"
#include "Molecule.h"

extern void velocities_from_PDB(char *filename, Vector *v, int totalAtoms);
extern void velocities_from_binfile(char *fname, Vector *vels, int n);

extern Vector gaussian_random_vector(void);

extern void random_velocities(BigReal Temp,
			Molecule *structure, Vector *v, int totalAtoms);
extern void remove_com_motion(Vector *vel, Molecule *structure, int n);
extern void read_binary_coors(char *fname, PDB *pdbobj);

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdOneTools.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1997/03/24 01:43:58 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdOneTools.h,v $
 * Revision 1.2  1997/03/24 01:43:58  jim
 * Added Langevin dynamics.
 *
 * Revision 1.1  1997/03/10 17:59:59  ari
 * Header for NamdOneTools.C - just externs
 *
 * Revision 1.1  1997/03/04 22:37:12  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 *
 ***************************************************************************/
