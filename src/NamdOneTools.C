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
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

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
  int n;			//  Number of atoms from file
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
  fread(&n, sizeof(int), 1, fp);

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
  fread(newcoords, sizeof(Vector), n, fp);

  //  Set the coordinates in the PDB object to the new coordinates
  pdbobj->set_all_positions(newcoords);

  //  Clean up
  fclose(fp);
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
		v[i].x *= 0.05;
		v[i].y *= 0.05;
		v[i].z *= 0.05;
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
	int filen;		//  Number of atoms read from file
	FILE *fp;		//  File descriptor

	//  Open the file and die if the open fails
	if ( (fp = Fopen(fname, "r")) == NULL)
	{
	  char errmsg[256];
  	  sprintf(errmsg, "Unable to open binary velocity file %s", fname);
	  NAMD_die(errmsg);
	}

	//  read the number of coordinates in this file
	fread(&filen, sizeof(int), 1, fp);

	//  Die if this doesn't match the number in our system
	if (filen != n)
	{
	  NAMD_die("Number of coordinates in binary velocity file incorrect");
	}

	fread(vels, sizeof(Vector), n, fp);

	fclose(fp);
}
/*		END OF FUNCTION velocities_from_binfile		*/

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


/***************************************************************************
* RCS INFORMATION:
*
*	$RCSfile: NamdOneTools.C,v $
*	$Author: ari $	$Locker:  $		$State: Exp $
*	$Revision: 1.6 $	$Date: 1997/04/06 22:45:06 $
*
***************************************************************************
* REVISION HISTORY:
*
* $Log: NamdOneTools.C,v $
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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/NamdOneTools.C,v 1.6 1997/04/06 22:45:06 ari Exp $";
