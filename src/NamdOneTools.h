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

extern void read_binary_coors(char *fname, PDB *pdbobj);
extern void vec_rotation_matrix(BigReal angle, Vector v, BigReal m[]);
extern Vector mat_multiply_vec(const Vector &v, BigReal m[]);


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdOneTools.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 1999/07/22 15:39:44 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdOneTools.h,v $
 * Revision 1.5  1999/07/22 15:39:44  jim
 * Eliminated last remnants of non-reentrant rand48 calls.
 *
 * Revision 1.4  1998/10/01 00:28:57  sergei
 * added vec_rotation_matrix (adopted from Mesa code) and mat_multiply_vec
 *
 * Revision 1.3  1998/01/05 20:34:09  sergei
 * added function BigReal gaussian_random_number(void);
 *
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
