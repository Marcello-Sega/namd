/***************************************************************************/
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: NAMD 1.X functions copied without change for use primarily
 *		in WorkDistrib
 *
 ***************************************************************************/
#include "charm++.h"

#include "common.h"
#include "NamdTypes.h"
#include "Vector.h"
#include "PDB.h"
#include "Molecule.h"
#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"
extern "C" long int lrand48(void);

//#define DEBUGM
#include "Debug.h"

/************************************************************************/
/*   FUNCTION read_binary_coors						*/
/*   INPUTS:								*/
/*	fname - Filename to read coordinates from			*/
/*	pdbobj - PDB object to place coordinates into			*/
/*	This function reads initial coordinates from a binary 		*/
/*	restart file							*/
/************************************************************************/

void read_binary_coors(char *fname, PDB *pdbobj) {
  int32 n;			//  Number of atoms from file
  Vector *newcoords;	//  Array of vectors to hold coordinates from file
  FILE *fp;		//  File descriptor

  //  Open the file and die if the open fails
  if ( (fp = Fopen(fname, "r")) == NULL)
  {
      char errmsg[256];
      sprintf(errmsg, "Unable to open binary coordinate file %s", fname);
      NAMD_die(errmsg);
  }

  //  read the number of coordinates in this file
  fread(&n, sizeof(int32), 1, fp);

  //  Die if this doesn't match the number in the system
  if (n != pdbobj->num_atoms()) {
      NAMD_die("Number of coordinates in binary coordinate file incorrect");
  }

  //  Allocate an array to hold the new coordinates
  newcoords = new Vector[n];

  if (newcoords == NULL) {
    NAMD_die("Memory alloc of newcoords in read_binary_coors() failed");
  }

  //  Read the coordinate from the file
  if (fread(newcoords, sizeof(Vector), n, fp) != (size_t)n)
  {
    NAMD_die("Error reading binary coordinate file");
  }


  //  Set the coordinates in the PDB object to the new coordinates
  pdbobj->set_all_positions(newcoords);

  //  Clean up
  Fclose(fp);
  delete [] newcoords;
} // END OF FUNCTION read_binary_coors()


/* 
 * Generate a 3x3 transformation matrix for rotation about a given
 * vector v by an angle given in degrees
 */
void vec_rotation_matrix( BigReal angle, Vector v, BigReal m[] ) {
   /* This function contributed by Erich Boleyn (erich@uruk.org) */
   BigReal mag, s, c;
   BigReal xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

   s = sin(angle * PI/180.0);
   c = cos(angle * PI/180.0);

   mag = v.length();

   if (mag == 0.0) {
      /* generate an identity matrix and return */
      for ( int i = 0; i < 9; ++i ) m[i] = 0.0;
      m[0] = m[4] = m[8] = 1.0;
      return;
   }

   // normalize the vector 
   v /= mag;

   /*
    *     Arbitrary axis rotation matrix.
    *
    *  This is composed of 5 matrices, Rz, Ry, T, Ry', Rz', multiplied
    *  like so:  Rz * Ry * T * Ry' * Rz'.  T is the final rotation
    *  (which is about the X-axis), and the two composite transforms
    *  Ry' * Rz' and Rz * Ry are (respectively) the rotations necessary
    *  from the arbitrary axis to the X-axis then back.  They are
    *  all elementary rotations.
    *
    *  Rz' is a rotation about the Z-axis, to bring the axis vector
    *  into the x-z plane.  Then Ry' is applied, rotating about the
    *  Y-axis to bring the axis vector parallel with the X-axis.  The
    *  rotation about the X-axis is then performed.  Ry and Rz are
    *  simply the respective inverse transforms to bring the arbitrary
    *  axis back to it's original orientation.  The first transforms
    *  Rz' and Ry' are considered inverses, since the data from the
    *  arbitrary axis gives you info on how to get to it, not how
    *  to get away from it, and an inverse must be applied.
    *
    *  The basic calculation used is to recognize that the arbitrary
    *  axis vector (x, y, z), since it is of unit length, actually
    *  represents the sines and cosines of the angles to rotate the
    *  X-axis to the same orientation, with theta being the angle about
    *  Z and phi the angle about Y (in the order described above)
    *  as follows:
    *
    *  cos ( theta ) = x / sqrt ( 1 - z^2 )
    *  sin ( theta ) = y / sqrt ( 1 - z^2 )
    *
    *  cos ( phi ) = sqrt ( 1 - z^2 )
    *  sin ( phi ) = z
    *
    *  Note that cos ( phi ) can further be inserted to the above
    *  formulas:
    *
    *  cos ( theta ) = x / cos ( phi )
    *  sin ( theta ) = y / sin ( phi )
    *
    *  ...etc.  Because of those relations and the standard trigonometric
    *  relations, it is pssible to reduce the transforms down to what
    *  is used below.  It may be that any primary axis chosen will give the
    *  same results (modulo a sign convention) using thie method.
    *
    *  Particularly nice is to notice that all divisions that might
    *  have caused trouble when parallel to certain planes or
    *  axis go away with care paid to reducing the expressions.
    *  After checking, it does perform correctly under all cases, since
    *  in all the cases of division where the denominator would have
    *  been zero, the numerator would have been zero as well, giving
    *  the expected result.
    */

   // store matrix in (row, col) form, i.e., m(row,col) = m[row*3+col]

   xx = v.x * v.x;
   yy = v.y * v.y;
   zz = v.z * v.z;
   xy = v.x * v.y;
   yz = v.y * v.z;
   zx = v.z * v.x;
   xs = v.x * s;
   ys = v.y * s;
   zs = v.z * s;
   one_c = 1.0 - c;

   m[0] = (one_c * xx) + c;
   m[1] = (one_c * xy) - zs;
   m[2] = (one_c * zx) + ys;
   
   m[3] = (one_c * xy) + zs;
   m[4] = (one_c * yy) + c;
   m[5] = (one_c * yz) - xs;
   
   m[6] = (one_c * zx) - ys;
   m[7] = (one_c * yz) + xs;
   m[8] = (one_c * zz) + c;
}

// multiply vector v by a 3x3 matrix m stored so that m(row,col) = m[row*3+col]
Vector mat_multiply_vec(const Vector &v, BigReal m[]) {
      return Vector( m[0]*v.x + m[1]*v.y + m[2]*v.z,
		     m[3]*v.x + m[4]*v.y + m[5]*v.z,
		     m[6]*v.x + m[7]*v.y + m[8]*v.z );
}

