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

