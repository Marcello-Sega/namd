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


/************************************************************************/
/*									*/
/*			FUNCTION velocities_from_PDB			*/
/*									*/
/*   INPUTS:								*/
/*	v - Array of vectors to populate				*/
/*	filename - name of the PDB filename to read in			*/
/*									*/
/*	This function reads in a set of initial velocities from a	*/
/*   PDB file.  It places the velocities into the array of Vectors      */
/*   passed to it.							*/
/*									*/
/************************************************************************/

void velocities_from_PDB(char *filename, Vector *v, int totalAtoms)
{
	PDB *v_pdb;		//  PDB info from velocity PDB
	int i;

	//  Read the PDB
	v_pdb = new PDB(filename);
	if ( v_pdb == NULL )
	{
	  NAMD_die("memory allocation failed in Node::velocities_from_PDB");
	}

	//  Make sure the number of velocities read in matches
	//  the number of atoms we have
	if (v_pdb->num_atoms() != totalAtoms)
	{
		char err_msg[129];

		sprintf(err_msg, "FOUND %d COORDINATES IN VELOCITY PDB!!",
		   v_pdb->num_atoms());

		NAMD_die(err_msg);
	}

	//  Get the entire list of atom info and loop through
	//  them assigning the velocity vector for each one
	v_pdb->get_all_positions(v);

	for (i=0; i<totalAtoms; i++)
	{
		v[i].x *= PDBVELINVFACTOR;
		v[i].y *= PDBVELINVFACTOR;
		v[i].z *= PDBVELINVFACTOR;
	}

	delete v_pdb;
}
/*		END OF FUNCTION velocities_from_PDB			*/


/************************************************************************/
/*			FUNCTION velocities_from_binfile		*/
/*   INPUTS:								*/
/*	fname - File name to write velocities to			*/
/*	n - Number of atoms in system					*/
/*	vels - Array of velocity vectors				*/
/*									*/
/*	This function writes out the velocities in binary format.	*/
/*      This is done to preserve accuracy between restarts of namd.	*/
/************************************************************************/
void velocities_from_binfile(char *fname, Vector *vels, int n)
{
	int32 filen;		//  Number of atoms read from file
	FILE *fp;		//  File descriptor

	//  Open the file and die if the open fails
	if ( (fp = Fopen(fname, "r")) == NULL)
	{
	  char errmsg[256];
  	  sprintf(errmsg, "Unable to open binary velocity file %s", fname);
	  NAMD_die(errmsg);
	}

	//  read the number of coordinates in this file
	fread(&filen, sizeof(int32), 1, fp);

	//  Die if this doesn't match the number in our system
	if (filen != n)
	{
	  NAMD_die("Number of coordinates in binary velocity file incorrect");
	}

	if (fread(vels, sizeof(Vector), n, fp) != (size_t)n)
	{
	  NAMD_die("Error reading binary velocity file");
	}


	Fclose(fp);
}
/*		END OF FUNCTION velocities_from_binfile		*/



/************************************************************************/
/*			FUNCTION gaussian_random_number			*/
/*	This function returns a normally distributed random number      */
/*   with zero mean and unit variance. See "Numerical Recipies in C"    */
/*   for details. To get a distribution with a variance s, multiply     */
/*   the returned value by s.                                           */
/************************************************************************/
BigReal gaussian_random_number(void) {
  static int iset = 0;
  static BigReal gset;
  BigReal fac, r, v1, v2;
  
  if (iset == 0) { // we do not have an extra result ready, so 
    r = 2.;                 // r >= 1.523e-8 ensures abs result < 6
    while (r >=1. || r < 1.523e-8) { // make sure we are within unit circle
      v1 = 2.0 * NAMD_random() - 1.0;
      v2 = 2.0 * NAMD_random() - 1.0;
      r = v1*v1 + v2*v2;
    }
    
    fac = sqrt(-2.0 * log(r)/r);
    // now make the Box-Muller transformation to get two normally 
    // distributed random numbers. Save one and return the other.
    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  }
  else { // use previously computed value
    iset = 0;
    return gset;
  }
}
/*			END OF FUNCTION gaussian_random_number		*/


Vector gaussian_random_vector(void)
{
	return Vector(  gaussian_random_number(),
			gaussian_random_number(),
			gaussian_random_number()  );
}


/*  This is very expensive so we'll use the above routine  */
/*
Vector gaussian_random_vector(void)
{
	int j;			//  Loop counter
	BigReal randnum;	//  Random number from -6.0 to 6.0
	Vector v;		//  Result

	//  This section generates a Gaussian random
	//  deviate of 0.0 mean and standard deviation RFD for
	//  each of the three spatial dimensions.
	//  The algorithm is a "sum of uniform deviates algorithm"
	//  which may be found in Abramowitz and Stegun,
	//  "Handbook of Mathematical Functions", pg 952.

	for (randnum = -6.0, j=0; j<12; j++) randnum += NAMD_random();

	v.x = randnum;

	for (randnum = -6.0, j=0; j<12; j++) randnum += NAMD_random();

	v.y = randnum;

	for (randnum = -6.0, j=0; j<12; j++) randnum += NAMD_random();

	v.z = randnum;

	return v;
}
*/

/************************************************************************/
/*			FUNCTION random_velocities			*/
/*   INPUTS:								*/
/*	v - array of vectors to populate				*/
/*	Temp - Temperature to acheive					*/
/*									*/
/*	This function assigns a random velocity distribution to a       */
/*   simulation to achieve a desired initial temperature.  The method   */
/*   used here was stolen from the program X-PLOR.			*/
/************************************************************************/

void random_velocities(BigReal Temp,
			Molecule *structure, Vector *v, int totalAtoms)
{
	int i, j;		//  Loop counter
	BigReal kbT;		//  Boltzman constant * Temp
	BigReal randnum;	//  Random number from -6.0 to 6.0
	BigReal kbToverM;	//  sqrt(Kb*Temp/Mass)

	kbT = Temp*BOLTZMAN;

	//  Loop through all the atoms and assign velocities in
	//  the x, y and z directions for each one
	for (i=0; i<totalAtoms; i++)
	{
		kbToverM = sqrt(kbT/structure->atommass(i));

		//  The following comment was stolen from X-PLOR where
		//  the following section of code was adapted from.

		//  This section generates a Gaussian random
		//  deviate of 0.0 mean and standard deviation RFD for
		//  each of the three spatial dimensions.
		//  The algorithm is a "sum of uniform deviates algorithm"
		//  which may be found in Abramowitz and Stegun,
		//  "Handbook of Mathematical Functions", pg 952.
		for (randnum=0.0, j=0; j<12; j++)
		{
			randnum += NAMD_random();
		}

		randnum -= 6.0;

		v[i].x = randnum*kbToverM;

		for (randnum=0.0, j=0; j<12; j++)
		{
			randnum += NAMD_random();
		}

		randnum -= 6.0;

		v[i].y = randnum*kbToverM;

		for (randnum=0.0, j=0; j<12; j++)
		{
			randnum += NAMD_random();
		}

		randnum -= 6.0;

		v[i].z = randnum*kbToverM;
	}
}
/*			END OF FUNCTION random_velocities		*/



/************************************************************************/
/*			FUNCTION remove_com_motion			*/
/*   INPUTS:								*/
/*	vel - Array of initial velocity vectors				*/
/*									*/
/*	This function removes the center of mass motion from a molecule.*/
/************************************************************************/

void remove_com_motion(Vector *vel, Molecule *structure, int n)
{
	Vector mv;		//  Sum of (mv)_i
	BigReal totalMass=0; 	//  Total mass of system
	int i;			//  Loop counter

	mv.x=0.0;
	mv.y=0.0;
	mv.z=0.0;

	//  Loop through and compute the net momentum
	for (i=0; i<n; i++)
	{
		mv.x += (structure->atommass(i))*vel[i].x;
		mv.y += (structure->atommass(i))*vel[i].y;
		mv.z += (structure->atommass(i))*vel[i].z;
		totalMass += structure->atommass(i);
	}

	mv.x = mv.x/totalMass;
	mv.y = mv.y/totalMass;
	mv.z = mv.z/totalMass;

	//  If any of the velocities really need to change, adjust them
	if ( (fabs(mv.x) > 0.0) || (fabs(mv.y) > 0.0) || (fabs(mv.y) > 0.0) )
	{
		iout << "ADJUSTING COM VELOCITY ("
			 << mv.x << ", " << mv.y << ", " << mv.z  
			 << ") TO REMOVE MOVEMENT\n" << endi;

		for (i=0; i<n; i++)
		{
			vel[i] -= mv;
		}
	}
}
/*			END OF FUNCTION remove_com_motion		*/


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


/***************************************************************************
* RCS INFORMATION:
*
*	$RCSfile: NamdOneTools.C,v $
*	$Author: jim $	$Locker:  $		$State: Exp $
*	$Revision: 1.18 $	$Date: 1999/06/04 14:09:08 $
*
***************************************************************************
* REVISION HISTORY:
*
* $Log: NamdOneTools.C,v $
* Revision 1.18  1999/06/04 14:09:08  jim
* Eliminated use of bzero().
*
* Revision 1.17  1999/05/11 23:56:37  brunner
* Changes for new charm version
*
* Revision 1.16  1998/10/24 19:57:42  jim
* Eliminated warnings generated by g++ -Wall.
*
* Revision 1.15  1998/10/01 00:28:56  sergei
* added vec_rotation_matrix (adopted from Mesa code) and mat_multiply_vec
*
* Revision 1.14  1998/09/01 23:10:34  brunner
* Fixed PDB velocity input conversion problem: PDBVELINVFACTOR
*
* Revision 1.13  1998/08/17 21:04:33  jim
* Added checks for short binary input files.
*
* Revision 1.12  1998/03/04 19:26:26  jim
* Switched to more efficient gaussian_random_vector implementation.
*
* Revision 1.11  1998/03/03 23:13:48  brunner
* Changing include files for new charm++ includes
*
* Revision 1.10  1998/01/05 20:34:08  sergei
* added function BigReal gaussian_random_number(void);
*
* Revision 1.9  1997/08/14 15:29:48  brunner
* More fixes for 32 bit ints in binary restart format
*
* Revision 1.8  1997/08/13 21:00:17  brunner
* Made binary files always use 32 bits for the number of atoms, so that it
* works on both 64 and 32-bit machines.  Also, I made Inform.C use CkPrintf,
* to fix the I/O buffering.
*
* Revision 1.7  1997/04/07 14:54:30  nealk
* Changed fclose() to Fclose() (found in common.[Ch]) to use with popen().
* Also corrected compilation warnings in Set.[Ch].
*
* Revision 1.6  1997/04/06 22:45:06  ari
* Add priorities to messages.  Mods to help proxies without computes.
* Added quick enhancement to end of list insertion of ResizeArray(s)
*
* Revision 1.5  1997/04/04 23:34:21  milind
* Got NAMD2 to run on Origin2000.
* Included definitions of class static variables in C files.
* Fixed alignment bugs by using memcpy instead of assignment in
* pack and unpack.
*
* Revision 1.4  1997/04/03 19:59:06  nealk
* 1) New Fopen() which handles .Z and .gz files.
* 2) localWaters and localNonWaters lists on each patch.
*
* Revision 1.3  1997/03/24 01:43:56  jim
* Added Langevin dynamics.
*
* Revision 1.2  1997/03/10 17:40:13  ari
* UniqueSet changes - some more commenting and cleanup
*
* Revision 1.1  1997/03/04 22:37:12  ari
* Clean up of code.  Debug statements removal, dead code removal.
* Minor fixes, output fixes.
* Commented some code from the top->down.  Mainly reworked Namd, Node, main.
*
*
***************************************************************************/

