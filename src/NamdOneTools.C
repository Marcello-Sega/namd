/***************************************************************************/
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "Vector.h"
#include "PDB.h"

//#define DEBUGM
#include "Debug.h"

/************************************************************************/
/*   FUNCTION read_binary_coors						*/
/*   INPUTS:								*/
/*	fname - Filename to read coordinates from			*/
/*	pdbobj - PDB object to place coordinates into			*/
/*	This function reads initial coordinates from a binary 		*/
/*	restart file							*
/************************************************************************/

void read_binary_coors(char *fname, PDB *pdbobj) {
  int n;			//  Number of atoms from file
  Vector *newcoords;	//  Array of vectors to hold coordinates from file
  FILE *fp;		//  File descriptor

  //  Open the file and die if the open fails
  if ( (fp = fopen(fname, "r")) == NULL)
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

/***************************************************************************
* RCS INFORMATION:
*
*	$RCSfile: NamdOneTools.C,v $
*	$Author: ari $	$Locker:  $		$State: Exp $
*	$Revision: 1.1 $	$Date: 1997/03/04 22:37:12 $
*
***************************************************************************
* DESCRIPTION:
*
***************************************************************************
* REVISION HISTORY:
*
* $Log: NamdOneTools.C,v $
* Revision 1.1  1997/03/04 22:37:12  ari
* Clean up of code.  Debug statements removal, dead code removal.
* Minor fixes, output fixes.
* Commented some code from the top->down.  Mainly reworked Namd, Node, main.
*
*
***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/NamdOneTools.C,v 1.1 1997/03/04 22:37:12 ari Exp $";
