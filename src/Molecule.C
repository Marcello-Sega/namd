/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Molecule.C,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1997/02/28 16:13:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *	The class Molecule is used to hold all of the structural information
 * for a simulation.  This information is read in from a .psf file and
 * cross checked with the Parameters object passed in.  All of the structural
 * information is then stored in arrays for use.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Molecule.C,v $
 * Revision 1.1002  1997/02/28 16:13:54  nealk
 * Turned off debugging code.
 *
 * Revision 1.1001  1997/02/10 08:14:37  jim
 * Fixed problem with exclusions where both modified and unmodified
 * versions of the same exclusion could be placed in the list, causing
 * one to be selected more or less randomly.  Also caused different
 * results on different numbers of processors.
 *
 * Revision 1.1000  1997/02/06 15:58:44  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.780  1997/02/06 15:53:15  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.2  1997/02/06 02:35:25  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778.2.1  1997/01/31 15:32:47  nealk
 * Corrected unititialized langForceVals which cause er-gre to core dump.
 *
 * Revision 1.778  1997/01/28 00:30:53  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:36:28  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.8  1996/12/06 06:54:41  jim
 * started from 1.3, added support for treating exclusions like bonds
 *
 * Revision 1.3  1996/11/11 19:54:09  nealk
 * Modified to use InfoStream instead of Inform.
 *
 * Revision 1.2  1996/08/16 04:39:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.29  1996/04/23 22:44:48  billh
 * Modified Atom structure to store a list of the hydrogen bond pairs
 * each atom participates in.
 *
 * Revision 1.28  1996/04/18 18:46:18  billh
 * Updated to read hydrogen bond information.
 *
 * Revision 1.27  1996/03/22 09:46:59  gursoy
 * corrected the commnents in the get_mother_atom
 *
 * Revision 1.26  96/03/05  11:08:02  11:08:02  gursoy (Attila Gursoy)
 * andrew fixed teh constrainf file related problem
 * 
 * Revision 1.25  96/02/16  17:06:34  17:06:34  gursoy (Attila Gursoy)
 * removed inline keywork from is_hydrogen and is_oxygen
 * 
 * Revision 1.24  96/02/16  16:39:28  16:39:28  gursoy (Attila Gursoy)
 * added atomType array and functions is_hydrogen, is_oxygen, get_mother_atom
 * and determine_atom_type
 * 
 * Revision 1.23  95/10/09  05:43:27  05:43:27  hazen (Brett Hazen)
 * Updated memory allocation to use C++ new/delete
 * 
 * Revision 1.22  1995/09/26  13:27:40  nelson
 * Added temperature coupling
 *
 * Revision 1.21  95/05/12  15:22:46  15:22:46  nelson (Mark T. Nelson)
 * Fixed problem where segment numbers like A000 caused namd to puke
 * on psf file, and fixed bug where separate conskfile caused
 * erroroneous error messages.
 * 
 * Revision 1.20  95/04/10  15:12:00  15:12:00  nelson (Mark T. Nelson)
 * Fixed bug where exclusions were freed even if never allocated
 * 
 * Revision 1.19  95/04/06  14:07:24  14:07:24  nelson (Mark T. Nelson)
 * Fixed problems where there are no bonds, angles, . . . etc.
 * 
 * Revision 1.18  95/03/08  14:44:58  14:44:58  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.17  95/02/23  11:37:35  11:37:35  nelson (Mark T. Nelson)
 * Fixed missing factor of 2 in random langevin forces
 * 
 * Revision 1.16  95/02/22  14:55:32  14:55:32  nelson (Mark T. Nelson)
 * Made changes for current working directory parameter
 * 
 * Revision 1.15  95/02/01  14:15:48  14:15:48  nelson (Mark T. Nelson)
 * Replaced references to namdDebug with DEBUG_MSG
 * 
 * Revision 1.14  95/01/26  14:17:28  14:17:28  nelson (Mark T. Nelson)
 * Made additions for Charm22 parameters
 * 
 * Revision 1.13  95/01/19  15:28:54  15:28:54  nelson (Mark T. Nelson)
 * Added langevin dynamics parameters
 * 
 * Revision 1.12  94/10/28  12:47:46  12:47:46  nelson (Mark T. Nelson)
 * Split out atom naming information into separate structure
 * from Atom
 * 
 * Revision 1.11  94/10/20  13:01:47  13:01:47  nelson (Mark T. Nelson)
 * Added code so that cosntraint PDB's will be shared if possible
 * 
 * Revision 1.10  94/10/19  21:59:11  21:59:11  nelson (Mark T. Nelson)
 * Added build_constraint_params for harmonic constrants
 * 
 * Revision 1.9  94/10/12  09:39:15  09:39:15  nelson (Mark T. Nelson)
 * Changed the way exclusions are handled so that bonded
 * exclusions are now pre-calculated and all exclusions
 * are stored in arrays by atom that can be searched via
 * binary searches.
 * 
 * Revision 1.8  94/10/04  10:35:44  10:35:44  nelson (Mark T. Nelson)
 * Corrected the routines that searched for 1-2, 1-3, and 1-4 interactions
 * 
 * Revision 1.7  94/09/24  20:11:05  20:11:05  nelson (Mark T. Nelson)
 * Added routines for pairlist generation that check for 1-2, 1-3, and 1-4
 * interactions.
 * 
 * Revision 1.6  94/09/13  14:32:36  14:32:36  gursoy (Attila Gursoy)
 * receive-message is moved to the Node object for charm++ integration
 * 
 * Revision 1.5  94/08/09  13:19:15  13:19:15  nelson (Mark T. Nelson)
 * Added code to deal with parameters that have no values in psf file
 * 
 * Revision 1.4  94/08/04  08:43:25  08:43:25  nelson (Mark T. Nelson)
 * Removed "RECEIVED MOLECULE" debug message
 * 
 * Revision 1.3  94/08/04  08:39:21  08:39:21  nelson (Mark T. Nelson)
 * Added send_Molecule and receive_Molecule functions
 * 
 * Revision 1.2  94/07/08  13:01:20  13:01:20  nelson (Mark T. Nelson)
 * I can't really remember what I did to it . . .
 * 
 * Revision 1.1  94/06/22  15:03:47  15:03:47  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Molecule.C,v 1.1002 1997/02/28 16:13:54 nealk Exp $";

#include "Molecule.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "strlib.h"
#include "InfoStream.h"
#include "IntList.h"
#include "Message.h"
// #include "Node.h"
#include "Parameters.h"
#include "PDB.h"
#include "SimParameters.h"

#define MIN_DEBUG_LEVEL 3
// #define DEBUGM
#include "Debug.h"

/************************************************************************/
/*									*/
/*			FUNCTION Molecule				*/
/*									*/
/*	This is the constructor for the Molecule class.  It simply sets */
/*  the counts for all the various parameters to 0 and sets the pointers*/
/*  to the arrays that will store these parameters to NULL, since they  */
/*  have not been allocated yet.					*/
/*									*/
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param, char *filename)
{
	this->simParams = simParams;
	/*  Initialize array pointers to NULL	*/
	atoms=NULL;
	atomNames=NULL;
	bonds=NULL;
	angles=NULL;
	dihedrals=NULL;
	impropers=NULL;
	donors=NULL;
	acceptors=NULL;
	exclusions=NULL;
	bondsByAtom=NULL;
	anglesByAtom=NULL;
	dihedralsByAtom=NULL;
	impropersByAtom=NULL;
	exclusionsByAtom=NULL;
	all_exclusions=NULL;
	onefour_exclusions=NULL;
	langevinParams=NULL;
	langForceVals=NULL;
	consIndexes=0;
	consParams=0;
	numMultipleDihedrals=0;
	numMultipleImpropers=0;

	/*  Initialize counts to 0 */
	numAtoms=0;
	numBonds=0;
	numAngles=0;
	numDihedrals=0;
	numImpropers=0;
	numDonors=0;
	numAcceptors=0;
	numExclusions=0;
	numConstraints=0;

	if (param != NULL && filename != NULL) {
	    read_psf_file(filename, param);
	}
}
/*			END OF FUNCTION Molecule			*/

/************************************************************************/
/*									*/
/*				FUNCTION Molecule			*/
/*									*/
/*	This is the destructor for the class Molecule.  It simply frees */
/*  the memory allocated for each of the arrays used to store the       */
/*  structure information.						*/
/*									*/
/************************************************************************/

Molecule::~Molecule()
{
	/*  Check to see if each array was ever allocated.  If it was   */
	/*  then free it						*/
	if (atoms != NULL)
		delete [] atoms;

	if (atomNames != NULL)
	{
	  for( int i=0; i<numAtoms; i++ )
	  {
	    delete [] atomNames[i].resname;
	    delete [] atomNames[i].atomname;
	    delete [] atomNames[i].atomtype;
	  }
	  delete [] atomNames;
	}

	if (bonds != NULL)
		delete [] bonds;

	if (angles != NULL)
		delete [] angles;

	if (dihedrals != NULL)
		delete [] dihedrals;

	if (impropers != NULL)
		delete [] impropers;

	if (donors != NULL)
		delete [] donors;

	if (acceptors != NULL)
		delete [] acceptors;

	if (exclusions != NULL)
		delete [] exclusions;
	
	if (bondsByAtom != NULL)
	   	delete [] bondsByAtom;
	
	if (anglesByAtom != NULL)
	   	delete [] anglesByAtom;
	
	if (dihedralsByAtom != NULL)
	   	delete [] dihedralsByAtom;
	
	if (impropersByAtom != NULL)
	   	delete [] impropersByAtom;
	
	if (exclusionsByAtom != NULL)
	   	delete [] exclusionsByAtom;
	
	if (all_exclusions != NULL)
	   	delete [] all_exclusions;
	
	if (onefour_exclusions != NULL)
	   	delete [] onefour_exclusions;
}
/*			END OF FUNCTION Molecule			*/

/************************************************************************/
/*									*/
/*				FUNCTION read_psf_file			*/
/*									*/
/*   INPUTS:								*/
/*	fname - Name of the .psf file to read				*/
/*	params - pointer to Parameters object to use to obtain          */
/*		 parameters for vdWs, bonds, etc.			*/
/*									*/
/*	This function reads a .psf file in.  This is where just about   */
/*   all of the structural information for this class comes from.  The  */
/*   .psf file contains descriptions of the atom, bonds, angles,        */
/*   dihedrals, impropers, and exclusions.  The parameter object is     */
/*   used to look up parameters for each of these entities.		*/
/*									*/
/************************************************************************/

void Molecule::read_psf_file(char *fname, Parameters *params)

{
	char err_msg[512];	//  Error message for NAMD_die
	char buffer[512];	//  Buffer for file reading
	int i;			//  Loop counter
	int NumTitle;		//  Number of Title lines in .psf file
	FILE *psf_file;		//  pointer to .psf file
	int ret_code;		//  ret_code from NAMD_read_line calls

	/* Try and open the .psf file 					*/
	if ( (psf_file = fopen(fname, "r")) == NULL)
	{
		sprintf(err_msg, "UNABLE TO OPEN .psf FILE %s", fname);
		NAMD_die(err_msg);
	}

	/*  Read till we have the first non-blank line of file		*/
	ret_code = NAMD_read_line(psf_file, buffer);

	while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
	{
		ret_code = NAMD_read_line(psf_file, buffer);
	}

	/*  Check to see if we dropped out of the loop because of a     */
	/*  read error.  This shouldn't happen unless the file is empty */
	if (ret_code!=0)
	{
		sprintf(err_msg, "EMPTY .psf FILE %s", fname);
		NAMD_die(err_msg);
	}

	/*  The first non-blank line should contain the word "psf".    */
	/*  If we can't find it, die.				       */
	if (!NAMD_find_word(buffer, "psf"))
	{
		sprintf(err_msg, "UNABLE TO FIND \"PSF\" STRING IN PSF FILE %s",
		   fname);
		NAMD_die(err_msg);
	}

	/*  Read until we find the next non-blank line			*/
	ret_code = NAMD_read_line(psf_file, buffer);

	while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
	{
		ret_code = NAMD_read_line(psf_file, buffer);
	}

	/*  Check to see if we dropped out of the loop because of a     */
	/*  read error.  This shouldn't happen unless there is nothing  */
	/*  but the PSF line in the file				*/
	if (ret_code!=0)
	{
		sprintf(err_msg, "MISSING EVERYTHING BUT PSF FROM %s", fname);
		NAMD_die(err_msg);
	}

	/*  This line should have the word "NTITLE" in it specifying    */
	/*  how many title lines there are				*/
	if (!NAMD_find_word(buffer, "NTITLE"))
	{
		sprintf(err_msg,"CAN NOT FIND \"NTITLE\" STRING IN PSF FILE %s",
		   fname);
		NAMD_die(err_msg);
	}

	sscanf(buffer, "%d", &NumTitle);

	/*  Now skip the next NTITLE non-blank lines and then read in the*/
	/*  line which should contain NATOM				*/
	i=0;

	while ( ((ret_code=NAMD_read_line(psf_file, buffer)) == 0) && 
		(i<NumTitle) )
	{
		if (!NAMD_blank_string(buffer))
			i++;
	}

	/*  Make sure we didn't exit because of a read error		*/
	if (ret_code!=0)
	{
		sprintf(err_msg, "FOUND EOF INSTEAD OF NATOM IN PSF FILE %s", 
		   fname);
		NAMD_die(err_msg);
	}

	while (NAMD_blank_string(buffer))
	{
		NAMD_read_line(psf_file, buffer);
	}

	/*  Check to make sure we have the line we want			*/
	if (!NAMD_find_word(buffer, "NATOM"))
	{
		sprintf(err_msg, "DIDN'T FIND \"NATOM\" IN PSF FILE %s",
		   fname);
		NAMD_die(err_msg);
	}

	/*  Read in the number of atoms, and then the atoms themselves	*/
	sscanf(buffer, "%d", &numAtoms);

	read_atoms(psf_file, params);

	/*  Read until we find the next non-blank line			*/
	ret_code = NAMD_read_line(psf_file, buffer);

	while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
	{
		ret_code = NAMD_read_line(psf_file, buffer);
	}

	/*  Check to make sure we didn't hit the EOF			*/
	if (ret_code != 0)
	{
		NAMD_die("EOF ENCOUNTERED LOOKING FOR NBONDS IN PSF");
	}

	/*  Look for the string "NBOND"					*/
	if (!NAMD_find_word(buffer, "NBOND"))
	{
		NAMD_die("DID NOT FIND NBOND AFTER ATOM LIST IN PSF");
	}

	/*  Read in the number of bonds and then the bonds themselves	*/
	sscanf(buffer, "%d", &numBonds);

	if (numBonds)
		read_bonds(psf_file, params);

	/*  Read until we find the next non-blank line			*/
	ret_code = NAMD_read_line(psf_file, buffer);

	while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
	{
		ret_code = NAMD_read_line(psf_file, buffer);
	}

	/*  Check to make sure we didn't hit the EOF			*/
	if (ret_code != 0)
	{
		NAMD_die("EOF ENCOUNTERED LOOKING FOR NTHETA IN PSF");
	}

	/*  Look for the string "NTHETA"				*/
	if (!NAMD_find_word(buffer, "NTHETA"))
	{
		NAMD_die("DID NOT FIND NTHETA AFTER BOND LIST IN PSF");
	}

	/*  Read in the number of angles and then the angles themselves */
	sscanf(buffer, "%d", &numAngles);

	if (numAngles)
		read_angles(psf_file, params);

	/*  Read until we find the next non-blank line			*/
	ret_code = NAMD_read_line(psf_file, buffer);

	while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
	{
		ret_code = NAMD_read_line(psf_file, buffer);
	}

	/*  Check to make sure we didn't hit the EOF			*/
	if (ret_code != 0)
	{
		NAMD_die("EOF ENCOUNTERED LOOKING FOR NPHI IN PSF");
	}

	/*  Look for the string "NPHI"					*/
	if (!NAMD_find_word(buffer, "NPHI"))
	{
		NAMD_die("DID NOT FIND NPHI AFTER ANGLE LIST IN PSF");
	}

	/*  Read in the number of dihedrals and then the dihedrals      */
	sscanf(buffer, "%d", &numDihedrals);

	if (numDihedrals)
		read_dihedrals(psf_file, params);

	/*  Read until we find the next non-blank line			*/
	ret_code = NAMD_read_line(psf_file, buffer);

	while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
	{
		ret_code = NAMD_read_line(psf_file, buffer);
	}

	/*  Check to make sure we didn't hit the EOF			*/
	if (ret_code != 0)
	{
		NAMD_die("EOF ENCOUNTERED LOOKING FOR NIMPHI IN PSF");
	}

	/*  Look for the string "NIMPHI"				*/
	if (!NAMD_find_word(buffer, "NIMPHI"))
	{
		NAMD_die("DID NOT FIND NIMPHI AFTER ATOM LIST IN PSF");
	}

	/*  Read in the number of Impropers and then the impropers	*/
	sscanf(buffer, "%d", &numImpropers);

	if (numImpropers)
		read_impropers(psf_file, params);

	/*  Read until we find the next non-blank line			*/
	ret_code = NAMD_read_line(psf_file, buffer);

	while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
	{
		ret_code = NAMD_read_line(psf_file, buffer);
	}

	/*  Check to make sure we didn't hit the EOF			*/
	if (ret_code != 0)
	{
		NAMD_die("EOF ENCOUNTERED LOOKING FOR NDON IN PSF");
	}

	/*  Look for the string "NDON"				*/
	if (!NAMD_find_word(buffer, "NDON"))
	{
		NAMD_die("DID NOT FIND NDON AFTER ATOM LIST IN PSF");
	}

	/*  Read in the number of hydrogen bond donors and then the donors */
	sscanf(buffer, "%d", &numDonors);

	if (numDonors)
		read_donors(psf_file);

	/*  Read until we find the next non-blank line			*/
	ret_code = NAMD_read_line(psf_file, buffer);

	while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
	{
		ret_code = NAMD_read_line(psf_file, buffer);
	}

	/*  Check to make sure we didn't hit the EOF			*/
	if (ret_code != 0)
	{
		NAMD_die("EOF ENCOUNTERED LOOKING FOR NACC IN PSF");
	}

	/*  Look for the string "NACC"				*/
	if (!NAMD_find_word(buffer, "NACC"))
	{
		NAMD_die("DID NOT FIND NACC AFTER ATOM LIST IN PSF");
	}

	/*  Read in the number of hydrogen bond donors and then the donors */
	sscanf(buffer, "%d", &numAcceptors);

	if (numAcceptors)
		read_acceptors(psf_file);

	/*  look for the explicit non-bonded exclusion section.     */
	while (!NAMD_find_word(buffer, "NNB"))
	{
		ret_code = NAMD_read_line(psf_file, buffer);

		if (ret_code != 0)
		{
			NAMD_die("EOF ENCOUNTERED LOOKING FOR NNB IN PSF FILE");
		}
	}

	/*  Read in the number of exclusions and then the exclusions    */
	sscanf(buffer, "%d", &numExclusions);

	if (numExclusions)
		read_exclusions(psf_file);

	/*  Close the .psf file.  There is a Group section in the .psf  */
	/*  file after the NNB section, but currently, NAMD does not    */
	/*  use this section for anything.				*/
	fclose(psf_file);

	//  analyze the data and find the status of each atom
	build_atom_status();

	return;
}
/*			END OF FUNCTION read_psf_file			*/

/************************************************************************/
/*									*/
/*				FUNCTION read_atoms			*/
/*									*/
/*   INPUTS:								*/
/*	fd - file pointer to the .psf file				*/
/*	params - Parameters object to use for parameters		*/
/*									*/
/*	this function reads in the Atoms section of the .psf file.      */
/*   This section consists of numAtoms lines that are of the form:	*/
/*     <atom#> <mol> <seg#> <res> <atomname> <atomtype> <charge> <mass> */
/*   Each line is read into the appropriate entry in the atoms array.   */
/*   The parameters object is then used to determine the vdW constants  */
/*   for this atom.							*/
/*									*/
/************************************************************************/

void Molecule::read_atoms(FILE *fd, Parameters *params)

{
	char buffer[512];	// Buffer for reading from file
	int atom_number=0;	// Atom number 
	int last_atom_number=0; // Last atom number, used to assure
				// atoms are in order
	char molecule_name[11]; // Molecule name
	char segment_number[21];// Segment number
	char residue_name[11];	// Residue name
	char atom_name[11];	// Atom name
	char atom_type[11];	// Atom type
	Real charge;		// Charge for the current atom
	Real mass;		// Mass for the current atom
	int read_count;		// Number of fields read by sscanf

	/*  Allocate the atom arrays					*/
	atoms     = new Atom[numAtoms];
	atomNames = new AtomNameInfo[numAtoms];

	if (atoms == NULL || atomNames == NULL )
	{
		NAMD_die("memory allocation failed in Molecule::read_atoms");
	}

	/*  Loop and read in numAtoms atom lines.			*/
	while (atom_number < numAtoms)
	{
		/*  Get the line from the file				*/
		NAMD_read_line(fd, buffer);

		/*  If its blank or a comment, skip it			*/
		if ( (NAMD_blank_string(buffer)) || (buffer[0] == '!') )
			continue;

		/*  Parse up the line					*/
		read_count=sscanf(buffer, "%d %s %s %s %s %s %f %f",
		   &atom_number, molecule_name, segment_number,
		   residue_name, atom_name, atom_type, &charge, &mass);

		/*  Check to make sure we found what we were expecting  */
		if (read_count != 8)
		{
			char err_msg[128];

			sprintf(err_msg, "BAD ATOM LINE FORMAT IN PSF FILE IN ATOM LINE %d\nLINE=%s",
			   last_atom_number+1, buffer);
			NAMD_die(err_msg);
		}

		/*  Make sure the atoms were in sequence		*/
		if (atom_number != last_atom_number+1)
		{
			char err_msg[128];

			sprintf(err_msg, "ATOM NUMBERS OUT OF ORDER AT ATOM #%d OF PSF FILE",
			   last_atom_number+1);
			NAMD_die(err_msg);
		}

		last_atom_number++;

		/*  Dynamically allocate strings for atom name, atom    */
		/*  type, etc so that we only allocate as much space    */
		/*  for these strings as we really need			*/
		atomNames[atom_number-1].resname = new char[strlen(residue_name)+1];
		atomNames[atom_number-1].atomname = new char[strlen(atom_name)+1];
		atomNames[atom_number-1].atomtype = new char[strlen(atom_type)+1];
	
		if (atomNames[atom_number-1].resname == NULL ||
		    atomNames[atom_number-1].atomname == NULL ||
		    atomNames[atom_number-1].atomtype == NULL )
		{
		  NAMD_die("memory allocation failed in Molecule::read_atoms");
		}

		/*  Put the values from this atom into the atoms array  */
		strcpy(atomNames[atom_number-1].resname, residue_name);
		strcpy(atomNames[atom_number-1].atomname, atom_name);
		strcpy(atomNames[atom_number-1].atomtype, atom_type);
		atoms[atom_number-1].mass = mass;
		atoms[atom_number-1].charge = charge;
		atoms[atom_number-1].status = UnknownAtom;

		/*  Determine the type of the atom (H or O) */
		if (atoms[atom_number-1].mass <=3.5) {
		  atoms[atom_number-1].status |= HydrogenAtom;
		} else if ((atomNames[atom_number-1].atomname[0] == 'O') && 
			   (atoms[atom_number-1].mass >= 14.0) && 
			   (atoms[atom_number-1].mass <= 18.0)) {
		  atoms[atom_number-1].status |= OxygenAtom;
		}

		/*  Look up the vdw constants for this atom		*/
		params->assign_vdw_index(atomNames[atom_number-1].atomtype, 
		   &(atoms[atom_number-1]));
        }

	return;
}
/*			END OF FUNCTION read_atoms			*/

/************************************************************************/
/*									*/
/*			FUNCTION read_bonds				*/
/*									*/
/*	read_bonds reads in the bond section of the .psf file.  This    */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are bonded together.  Each atom pair is   */
/*  read in.  Then that parameter object is queried to determine the    */
/*  force constant and rest distance for the bond.			*/
/*									*/
/************************************************************************/

void Molecule::read_bonds(FILE *fd, Parameters *params)

{
	int atom_nums[2];	// Atom indexes for the bonded atoms
	char atom1name[11];	// Atom type for atom #1
	char atom2name[11];	// Atom type for atom #2
	char tmp_string[11];	// Temporary string to read in atom #'s
	int j;			// Loop counter
	int num_read=0;		// Number of bonds read so far

	/*  Allocate the array to hold the bonds			*/
	bonds=new Bond[numBonds];

	if (bonds == NULL)
	{
		NAMD_die("memory allocations failed in Molecule::read_bonds");
	}

	/*  Loop through and read in all the bonds			*/
	while (num_read < numBonds)
	{
		/*  Loop and read in the two atom indexes		*/
		for (j=0; j<2; j++)
		{
			/*  Read the atom number from the file.         */
			/*  Subtract 1 to convert the index from the    */
			/*  1 to NumAtoms used in the file to the       */
			/*  0 to NumAtoms-1 that we need		*/
			atom_nums[j]=NAMD_read_int(fd, "BONDS")-1;

			/*  Check to make sure the index isn't too big  */
			if (atom_nums[j] >= numAtoms)
			{
				char err_msg[128];

				sprintf(err_msg, "BOND INDEX %d GREATER THAN NATOM %d IN BOND # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
				NAMD_die(err_msg);
			}
		}

		/*  Get the atom type for the two atoms.  When we query */
		/*  the parameter object, we need to send the atom type */
		/*  that is alphabetically first as atom 1.		*/
		if (strcasecmp(atomNames[atom_nums[0]].atomtype, 
			   atomNames[atom_nums[1]].atomtype) < 0)
		{
			strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
			strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
		}
		else
		{
			strcpy(atom2name, atomNames[atom_nums[0]].atomtype);
			strcpy(atom1name, atomNames[atom_nums[1]].atomtype);
		}

		/*  Query the parameter object for the constants for    */
		/*  this bond						*/
		Bond *b = &(bonds[num_read]);
		params->assign_bond_index(atom1name, atom2name, b);

		/*  Assign the atom indexes to the array element	*/
		b->atom1=atom_nums[0];
		b->atom2=atom_nums[1];

		num_read++;
	}

	return;
}
/*			END OF FUNCTION read_bonds			*/

/************************************************************************/
/*									*/
/*			FUNCTION read_angles				*/
/*									*/
/*   INPUTS:								*/
/*	fd - File descriptor for .psf file				*/
/*	params - Parameters object to query for parameters		*/
/*									*/
/*	read_angles reads the angle parameters from the .psf file.      */
/*   This section of the .psf file consists of a list of triplets of    */
/*   atom indexes.  Each triplet represents three atoms connected via   */
/*   an angle bond.  The parameter object is queried to obtain the      */
/*   constants for each bond.						*/
/*									*/
/************************************************************************/

void Molecule::read_angles(FILE *fd, Parameters *params)

{
	int atom_nums[3];	//  Atom numbers for the three atoms
	char atom1name[11];	//  Atom type for atom 1
	char atom2name[11];	//  Atom type for atom 2
	char atom3name[11];	//  Atom type for atom 3
	char tmp_string[11];	//  Temporary string for reading atoms
	int j;			//  Loop counter
	int num_read=0;		//  Number of angles read so far

	/*  Alloc the array of angles					*/
	angles=new Angle[numAngles];

	if (angles == NULL)
	{
		NAMD_die("memory allocation failed in Molecule::read_angles");
	}

	/*  Loop through and read all the angles			*/
	while (num_read < numAngles)
	{
		/*  Loop through the 3 atom indexes in the current angle*/
		for (j=0; j<3; j++)
		{
			/*  Read the atom number from the file.         */
			/*  Subtract 1 to convert the index from the    */
			/*  1 to NumAtoms used in the file to the       */
			/*  0 to NumAtoms-1 that we need		*/
			atom_nums[j]=NAMD_read_int(fd, "ANGLES")-1;

			/*  Check to make sure the atom index doesn't   */
			/*  exceed the Number of Atoms			*/
			if (atom_nums[j] >= numAtoms)
			{
				char err_msg[128];

				sprintf(err_msg, "ANGLES INDEX %d GREATER THAN NATOM %d IN ANGLES # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
				NAMD_die(err_msg);
			}
		}

		/*  Place the bond name that is alphabetically first	*/
		/*  in the atom1name.  This is OK since the order of    */
		/*  atom1 and atom3 are interchangable.  And to search  */
		/*  the tree of angle parameters, we need the order     */
		/*  to be predictable.					*/
		if (strcasecmp(atomNames[atom_nums[0]].atomtype, 
			   atomNames[atom_nums[2]].atomtype) < 0)
		{
			strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
			strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
			strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
		}
		else
		{
			strcpy(atom1name, atomNames[atom_nums[2]].atomtype);
			strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
			strcpy(atom3name, atomNames[atom_nums[0]].atomtype);
		}

		/*  Get the constant values for this bond from the	*/
		/*  parameter object					*/
		params->assign_angle_index(atom1name, atom2name, 
		   atom3name, &(angles[num_read]));

		/*  Assign the three atom indices			*/
		angles[num_read].atom1=atom_nums[0];
		angles[num_read].atom2=atom_nums[1];
		angles[num_read].atom3=atom_nums[2];

		num_read++;
	}

	return;
}
/*			END OF FUNCTION read_angles			*/

/************************************************************************/
/*									*/
/*				FUNCTION read_dihedrals			*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor for the .psf file				*/
/*	params - pointer to parameter object				*/
/*									*/
/*	read_dihedreals reads the dihedral section of the .psf file.    */
/*   This section of the file contains a list of quartets of atom       */
/*   numbers.  Each quartet represents a group of atoms that form a     */
/*   dihedral bond.							*/
/*									*/
/************************************************************************/

void Molecule::read_dihedrals(FILE *fd, Parameters *params)

{
	int atom_nums[4];	// The 4 atom indexes
	int last_atom_nums[4];  // Atom numbers from previous bond
	char atom1name[11];	// Atom type for atom 1
	char atom2name[11];	// Atom type for atom 2
	char atom3name[11];	// Atom type for atom 3
	char atom4name[11];	// Atom type for atom 4
	char tmp_string[11];	// Temporary string for reading indexes
	int j;			// loop counter
	int num_read=0;		// number of dihedrals read so far
	int multiplicity=1;	// multiplicity of the current bond
	Bool duplicate_bond;	// Is this a duplicate of the last bond
	int num_unique=0; 	// Number of unique dihedral bonds

	//  Initialize the array used to check for duplicate dihedrals
	for (j=0; j<4; j++)
		last_atom_nums[j] = -1;

	/*  Allocate an array to hold the Dihedrals			*/
	dihedrals = new Dihedral[numDihedrals];

	if (dihedrals == NULL)
	{
		NAMD_die("memory allocation failed in Molecule::read_dihedrals");
	}

	/*  Loop through and read all the dihedrals			*/
	while (num_read < numDihedrals)
	{
		duplicate_bond = TRUE;

		/*  Loop through and read the 4 indexes for this bond   */
		for (j=0; j<4; j++)
		{
			/*  Read the atom number from the file.         */
			/*  Subtract 1 to convert the index from the    */
			/*  1 to NumAtoms used in the file to the       */
			/*  0 to NumAtoms-1 that we need		*/
			atom_nums[j]=NAMD_read_int(fd, "DIHEDRALS")-1;

			/*  Check for an atom index that is too large	*/
			if (atom_nums[j] >= numAtoms)
			{
				char err_msg[128];

				sprintf(err_msg, "DIHEDRALS INDEX %d GREATER THAN NATOM %d IN DIHEDRALS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
				NAMD_die(err_msg);
			}

			//  Check to see if this atom matches the last bond
			if (atom_nums[j] != last_atom_nums[j])
			{
				duplicate_bond = FALSE;
			}

			last_atom_nums[j] = atom_nums[j];
		}

		/*  Get the atom types for the 4 atoms so we can look	*/
		/*  up the constants in the parameter object		*/
		strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
		strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
		strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
		strcpy(atom4name, atomNames[atom_nums[3]].atomtype);

		//  Check to see if this is really a new bond or just
		//  a repeat of the last one
		if (duplicate_bond)
		{
			//  This is a duplicate, so increase the multiplicity
			multiplicity++;

			if (multiplicity == 2)
			{
				numMultipleDihedrals++;
			}
		}
		else
		{
			multiplicity=1;
			num_unique++;
		}

		/*  Get the constants for this dihedral bond		*/
		params->assign_dihedral_index(atom1name, atom2name, 
		   atom3name, atom4name, &(dihedrals[num_unique-1]),
		   multiplicity);

		/*  Assign the atom indexes				*/
		dihedrals[num_unique-1].atom1=atom_nums[0];
		dihedrals[num_unique-1].atom2=atom_nums[1];
		dihedrals[num_unique-1].atom3=atom_nums[2];
		dihedrals[num_unique-1].atom4=atom_nums[3];

		num_read++;
	}

	numDihedrals = num_unique;

	return;
}
/*			END OF FUNCTION read_dihedral			*/

/************************************************************************/
/*									*/
/*				FUNCTION read_impropers			*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor for .psf file				*/
/*	params - parameter object					*/
/*									*/
/*	read_impropers reads the improper section of the .psf file.	*/
/*   This section is identical to the dihedral section in that it is    */
/*   made up of a list of quartets of atom indexes that define the      */
/*   atoms that are bonded together.					*/
/*									*/
/************************************************************************/

void Molecule::read_impropers(FILE *fd, Parameters *params)

{
	int atom_nums[4];	//  Atom indexes for the 4 atoms
	int last_atom_nums[4];	//  Atom indexes from previous bond
	char atom1name[11];	//  Atom type for atom 1
	char atom2name[11];	//  Atom type for atom 2
	char atom3name[11];	//  Atom type for atom 3
	char atom4name[11];	//  Atom type for atom 4
	char tmp_string[11];	//  Temporary string to read in indexes
	int j;			//  Loop counter
	int num_read=0;		//  Number of impropers read so far
	int multiplicity=1;	// multiplicity of the current bond
	Bool duplicate_bond;	// Is this a duplicate of the last bond
	int num_unique=0; 	// Number of unique dihedral bonds

	//  Initialize the array used to look for duplicate improper
	//  entries.  Set them all to -1 so we know nothing will match
	for (j=0; j<4; j++)
		last_atom_nums[j] = -1;

	/*  Allocate the array to hold the impropers			*/
	impropers=new Improper[numImpropers];

	if (impropers == NULL)
	{
	  NAMD_die("memory allocation failed in Molecule::read_impropers");
	}

	/*  Loop through and read all the impropers			*/
	while (num_read < numImpropers)
	{
		duplicate_bond = TRUE;

		/*  Loop through the 4 indexes for this improper	*/
		for (j=0; j<4; j++)
		{
			/*  Read the atom number from the file.         */
			/*  Subtract 1 to convert the index from the    */
			/*  1 to NumAtoms used in the file to the       */
			/*  0 to NumAtoms-1 that we need		*/
			atom_nums[j]=NAMD_read_int(fd, "IMPROPERS")-1;

			/*  Check to make sure the index isn't too big  */
			if (atom_nums[j] >= numAtoms)
			{
				char err_msg[128];

				sprintf(err_msg, "IMPROPERS INDEX %d GREATER THAN NATOM %d IN IMPROPERS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
				NAMD_die(err_msg);
			}

			if (atom_nums[j] != last_atom_nums[j])
			{
				duplicate_bond = FALSE;
			}

			last_atom_nums[j] = atom_nums[j];
		}

		/*  Get the atom types so we can look up the parameters */
		strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
		strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
		strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
		strcpy(atom4name, atomNames[atom_nums[3]].atomtype);

		//  Check to see if this is a duplicate improper
		if (duplicate_bond)
		{
			//  This is a duplicate improper.  So we don't
			//  really count this entry, we just update
			//  the parameters object
			multiplicity++;

			if (multiplicity == 2)
			{
				//  Count the number of multiples.
				numMultipleImpropers++;
			}
		}
		else
		{
			//  Not a duplicate
			multiplicity = 1;
			num_unique++;
		}

		/*  Look up the constants for this bond			*/
		params->assign_improper_index(atom1name, atom2name, 
		   atom3name, atom4name, &(impropers[num_unique-1]),
		   multiplicity);

		/*  Assign the atom indexes				*/
		impropers[num_unique-1].atom1=atom_nums[0];
		impropers[num_unique-1].atom2=atom_nums[1];
		impropers[num_unique-1].atom3=atom_nums[2];
		impropers[num_unique-1].atom4=atom_nums[3];

		num_read++;
	}

	//  Now reset the numImpropers value to the number of UNIQUE
	//  impropers.  Sure, we waste a few entries in the improper_array
	//  on the master node, but it is very little space . . .
	numImpropers = num_unique;

	return;
}
/*			END OF FUNCTION read_impropers			*/

/************************************************************************/
/*									*/
/*			FUNCTION read_donors				*/
/*									*/
/*	read_donors reads in the bond section of the .psf file.  This   */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are part of an H-bond.  Each atom pair is */
/*  read in.                                                            */
/*									*/
/*  Donor atoms are the heavy atoms to which hydrogens are bonded.      */
/*  There will always be a donor atom for each donor pair.  However,    */
/*  for a united-atom model there may not be an explicit hydrogen       */
/*  present, in which case the second atom index in the pair will be    */
/*  given as 0 in the PSF (and stored as -1 in this program's internal  */
/*  storage).                                                           */
/************************************************************************/

void Molecule::read_donors(FILE *fd)

{
	int d[2];               // temporary storage of donor atom index
	int j;			// Loop counter
	int num_read=0;		// Number of bonds read so far
	int num_no_hydr=0;      // Number of bonds with no hydrogen given
	char tmp_string[11];	// Temporary string to read in atom #'s

	/*  Allocate the array to hold the bonds			*/
	donors=new Bond[numDonors];

	if (donors == NULL)
	{
		NAMD_die("memory allocations failed in Molecule::read_donors");
	}

	/*  Loop through and read in all the donors			*/
	while (num_read < numDonors)
	{
		/*  Loop and read in the two atom indexes		*/
		for (j=0; j<2; j++)
		{
			/*  Read the atom number from the file.         */
			/*  Subtract 1 to convert the index from the    */
			/*  1 to NumAtoms used in the file to the       */
			/*  0 to NumAtoms-1 that we need		*/
			d[j]=NAMD_read_int(fd, "DONORS")-1;

			/*  Check to make sure the index isn't too big  */
			if (d[j] >= numAtoms)
			{
				char err_msg[128];

				sprintf(err_msg,
		"DONOR INDEX %d GREATER THAN NATOM %d IN DONOR # %d IN PSF FILE",
					d[j]+1, numAtoms, num_read+1);
				NAMD_die(err_msg);
			}

			/*  Check if there is a hydrogen given */
			if (d[j] < 0)
                          num_no_hydr++;
		}

		/*  Assign the atom indexes to the array element	*/
		Bond *b = &(donors[num_read]);
		b->atom1=d[0];
		b->atom2=d[1];

		num_read++;
	}

	return;
}
/*			END OF FUNCTION read_donors			*/


/************************************************************************/
/*									*/
/*			FUNCTION read_acceptors				*/
/*									*/
/*	read_acceptors reads in the bond section of the .psf file.      */
/*  This section contains a list of pairs of numbers where each pair is */
/*  represents two atoms that are part of an H-bond.  Each atom pair is */
/*  read in.                                                            */
/*									*/
/*  Acceptor atoms are the heavy atoms to which hydrogens directly      */
/*  orient in a hydrogen bond interaction.  There will always be an     */
/*  acceptor atom for each acceptor pair.  The antecedent atom, to      */
/*  which the acceptor is bound, may not be given in the structure,     */
/*  however, in which case the second atom index in the pair will be    */
/*  given as 0 in the PSF (and stored as -1 in this program's internal  */
/*  storage).                                                           */
/************************************************************************/

void Molecule::read_acceptors(FILE *fd)

{
	int d[2];               // temporary storage of atom index
	int j;			// Loop counter
	int num_read=0;		// Number of bonds read so far
        int num_no_ante=0;      // number of pairs with no antecedent
	char tmp_string[11];	// Temporary string to read in atom #'s

	/*  Allocate the array to hold the bonds			*/
	acceptors=new Bond[numAcceptors];

	if (acceptors == NULL)
	{
		NAMD_die("memory allocations failed in Molecule::read_acceptors");
	}

	/*  Loop through and read in all the acceptors			*/
	while (num_read < numAcceptors)
	{
		/*  Loop and read in the two atom indexes		*/
		for (j=0; j<2; j++)
		{
			/*  Read the atom number from the file.         */
			/*  Subtract 1 to convert the index from the    */
			/*  1 to NumAtoms used in the file to the       */
			/*  0 to NumAtoms-1 that we need		*/
			d[j]=NAMD_read_int(fd, "ACCEPTORS")-1;

			/*  Check to make sure the index isn't too big  */
			if (d[j] >= numAtoms)
			{
				char err_msg[128];

				sprintf(err_msg, "ACCEPTOR INDEX %d GREATER THAN NATOM %d IN DONOR # %d IN PSF FILE", d[j]+1, numAtoms, num_read+1);
				NAMD_die(err_msg);
			}

			/*  Check if there is an antecedent given */
			if (d[j] < 0)
                          num_no_ante++;
		}

		/*  Assign the atom indexes to the array element	*/
		Bond *b = &(acceptors[num_read]);
		b->atom1=d[0];
		b->atom2=d[1];

		num_read++;
	}

	return;
}
/*			END OF FUNCTION read_acceptors			*/


/************************************************************************/
/*									*/
/*			FUNCTION read_exclusions			*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor for .psf file				*/
/*									*/
/*	read_exclusions reads in the explicit non-bonded exclusions     */
/*  from the .psf file.  This section is a little funky, so hang on.    */
/*  Ok, first there is a list of atom indexes that is NumExclusions     */
/*  long.  These are in some sense the atoms that will be exlcuded.     */
/*  Following this list is a list of NumAtoms length that is a list     */
/*  of indexes into the list of excluded atoms.  So an example.  Suppose*/
/*  we have a 5 atom simulation with 3 explicit exclusions.  The .psf   */
/*  file could look like:						*/
/*									*/
/*	3!NNB								*/
/*	3 4 5								*/
/*	0 1 3 3 3							*/
/*									*/
/*  This would mean that atom 1 has no explicit exclusions.  Atom 2     */
/*  has an explicit exclusion with atom 3.  Atom 3 has an explicit      */
/*  exclusion with atoms 4 AND 5.  And atoms 4 and 5 have no explicit   */
/*  exclusions.  Got it!?!  I'm not sure who dreamed this up . . .      */
/*									*/
/************************************************************************/

void Molecule::read_exclusions(FILE *fd)

{
	int *exclusion_atoms;	//  Array of indexes of excluded atoms
	char tmp_string[11];	//  temporary string for readin in values
	int num_read=0;		//  Number fo exclusions read in
	int current_index;	//  Current index value
	int last_index;		//  the previous index value
	int insert_index=0;	//  index of where we are in exlcusions array

	/*  Allocate the array of exclusion structures and the array of */
	/*  exlcuded atom indexes					*/
	exclusions      = new Exclusion[numExclusions];
	exclusion_atoms = new int[numExclusions];

	if ( (exclusions == NULL) || (exclusion_atoms == NULL) )
	{
	  NAMD_die("memory allocation failed in Molecule::read_exclusions");
	}

	/*  First, read in the excluded atoms list			*/
	for (num_read=0; num_read<numExclusions; num_read++)
	{
		/*  Read the atom number from the file. Subtract 1 to   */
		/*  convert the index from the 1 to NumAtoms used in the*/
		/*  file to the  0 to NumAtoms-1 that we need		*/
		exclusion_atoms[num_read]=NAMD_read_int(fd, "IMPROPERS")-1;

		/*  Check for an illegal index				*/
		if (exclusion_atoms[num_read] >= numAtoms)
		{
			char err_msg[128];

			sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN PSF FILE", exclusion_atoms[num_read]+1, numAtoms, num_read+1);
			NAMD_die(err_msg);
		}
	}

	/*  Now, go through and read the list of NumAtoms pointers into */
	/*  the array that we just read in				*/
	last_index=0;

	for (num_read=0; num_read<numAtoms; num_read++)
	{
		/*  Read in the current index value			*/
		current_index=NAMD_read_int(fd, "EXCLUSIONS");

		/*  Check for an illegal pointer			*/
		if (current_index>numExclusions)
		{
			char err_msg[128];

			sprintf(err_msg, "EXCLUSION INDEX %d LARGER THAN NUMBER OF EXLCUSIONS %d IN PSF FILE, EXCLUSION #%d\n", 
			   current_index+1, numExclusions, num_read);
			NAMD_die(err_msg);
		}

		/*  Check to see if it matches the last index.  If so   */
		/*  than this atom has no exclusions.  If not, then     */
		/*  we have to build some exclusions			*/
		if (current_index != last_index)
		{
			/*  This atom has some exlcusions.  Loop from   */
			/*  the last_index to the current index.  This  */
			/*  will include how ever many exclusions this  */
			/*  atom has					*/
			for (insert_index=last_index; 
			     insert_index<current_index; insert_index++)
			{
				/*  Assign the two atoms involved.      */
				/*  The first one is our position in    */
				/*  the list, the second is based on    */
				/*  the pointer into the index list     */
				exclusions[insert_index].atom1=num_read;
				exclusions[insert_index].atom2=exclusion_atoms[insert_index];
			}

			last_index=current_index;
		}
	}

	/*  Free our temporary list of indexes				*/
	delete [] exclusion_atoms;

	return;
}
/*			END OF FUNCTION read_exclusions			*/

/************************************************************************/
/*									*/
/*			FUNCTION print_atoms				*/
/*									*/
/*	print_atoms prints out the list of atoms stored in this object. */
/*  It is inteded mainly for debugging purposes.			*/
/*									*/
/************************************************************************/

void Molecule::print_atoms(Parameters *params)

{
	int i;
	Real sigma;
	Real epsilon;
	Real sigma14;
	Real epsilon14;

	DEBUG_MSG("ATOM LIST\n" \
		  << "******************************************\n" \
                  << "NUM  NAME TYPE RES  MASS    CHARGE SIGMA   EPSILON SIGMA14 EPSILON14\n" \
pp		  << endi);

	for (i=0; i<numAtoms; i++)
	{
		params->get_vdw_params(&sigma, &epsilon, &sigma14, &epsilon14, 
				atoms[i].vdw_type);

		DEBUG_MSG(i+1 << " " << atomNames[i].atomname  \
		          << " " << atomNames[i].atomtype << " " \
		          << atomNames[i].resname  << " " << atoms[i].mass  \
			  << " " << atoms[i].charge << " " << sigma \
			  << " " << epsilon << " " << sigma14 \
			  << " " << epsilon14 \
			  << endi);
	}
}
/*			END OF FUNCTION print_atoms			*/

/************************************************************************/
/*									*/
/*			FUNCTION print_bonds				*/
/*									*/
/*	print_bonds prints out the list of bonds stored in this object. */
/*  It is inteded mainly for debugging purposes.			*/
/*									*/
/************************************************************************/

void Molecule::print_bonds(Parameters *params)

{
	int i;
	Real k;
	Real x0;

	DEBUG_MSG("BOND LIST\n" << "********************************\n" \
		  << "ATOM1 ATOM2 TYPE1 TYPE2      k        x0" \
		  << endi);

	for (i=0; i<numBonds; i++)
	{
		params->get_bond_params(&k, &x0, bonds[i].bond_type);

		DEBUG_MSG(bonds[i].atom1+1 << " " \
		   << bonds[i].atom2+1 << " "   \
		   << atomNames[bonds[i].atom1].atomtype << " "  \
		   << atomNames[bonds[i].atom2].atomtype << " " << k \
		   << " " << x0 << endi);
	}
}
/*			END OF FUNCTION print_bonds			*/

/************************************************************************/
/*									*/
/*			FUNCTION print_exclusions			*/
/*									*/
/*	print_exlcusions prints out the list of exlcusions stored in    */
/*  this object.  It is inteded mainly for debugging purposes.		*/
/*									*/
/************************************************************************/

void Molecule::print_exclusions()
{
	int i;

	DEBUG_MSG("EXPLICIT EXCLUSION LIST\n" \
		  << "********************************\n" \
	          << "ATOM1 ATOM2 " \
		  << endi);

	for (i=0; i<numExclusions; i++)
	{
		DEBUG_MSG(exclusions[i].atom1+1 << "  " \
		   << exclusions[i].atom2+1 << endi);
	}
}
/*			END OF FUNCTION print_exclusions		*/

/************************************************************************/
/*									*/
/*			FUNCTION send_Molecule				*/
/*									*/
/*	send_Molecule is used by the Master node to distribute the      */
/*   structural information to all the client nodes.  It is NEVER called*/
/*   by the client nodes.  						*/
/*									*/
/*   NOTE:  Because the structures that are being passed here are not   */
/*  	    directly recognized by the Message and Communication classes*/
/*	    the data is first copied to arrays of basic types (i.e., int*/
/*	    Real, Index, etc.).  While this is a little inefficient,    */
/*	    this communication is only done once per simulation.	*/
/*									*/
/************************************************************************/


void Molecule::send_Molecule(Communicate *com_obj)

{
	Message *msg=new Message; //  Message to send to clients
	Real *a1, *a2;		//  Arrays of reals used to send data
	Index *ind1;		//  Array of Indexes used to send data
	int *i1, *i2, *i3, *i4;	//  Array of ints used to send data
	Vector *v1;		//  Array of vectors used to send data
	int i;			//  Loop counter
	if ( msg == NULL )
	{
	  NAMD_die("Memory allocation failed in Molecule::send_Molecule");
	}

	//  Send the atom information
	a1   = new Real[numAtoms];
	a2   = new Real[numAtoms];
	ind1 = new Index[numAtoms];
	i1   = new int[numAtoms];

	if ( (a1==NULL) || (a2==NULL) || (ind1==NULL) || (i1==NULL) )
	{
	  NAMD_die("memory allocation failed in Molecule::send_Molecule");
	}

	for (i=0; i<numAtoms; i++)
	{
		a1[i]=atoms[i].mass;
		a2[i]=atoms[i].charge;
		ind1[i]=atoms[i].vdw_type;
		i1[i]=atoms[i].status;
	}

	msg->put(numAtoms);
	msg->put(numAtoms, a1).put(numAtoms, a2);
	msg->put(numAtoms, ind1).put(numAtoms, i1);

	delete [] a1;
	delete [] a2;
	delete [] ind1;
	delete [] i1;

	//  Send the bond information
	msg->put(numBonds);

	if (numBonds)
	{
		i1= new int[numBonds];
		i2= new int[numBonds];
		ind1= new Index[numBonds];

		if ( (i1==NULL) || (i2==NULL) || (ind1==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::send_Molecule");
		}

		for (i=0; i<numBonds; i++)
		{
			i1[i]=bonds[i].atom1;
			i2[i]=bonds[i].atom2;
			ind1[i]=bonds[i].bond_type;
		}
	
		msg->put(numBonds, i1).put(numBonds, i2);
		msg->put(numBonds, ind1);
	
		delete [] i1;
		delete [] i2;
		delete [] ind1;
	}

	//  Send the angle information
	msg->put(numAngles);

	if (numAngles)
	{
		i1= new int[numAngles];
		i2= new int[numAngles];
		i3= new int[numAngles];
		ind1= new Index[numAngles];

		if ( (i1==NULL) || (i2==NULL) || (i3==NULL) || (ind1==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::send_Molecule");
		}

		for (i=0; i<numAngles; i++)
		{
			i1[i]=angles[i].atom1;
			i2[i]=angles[i].atom2;
			i3[i]=angles[i].atom3;
			ind1[i]=angles[i].angle_type;
		}

		msg->put(numAngles, i1).put(numAngles, i2);
		msg->put(numAngles, i3).put(numAngles, ind1);

		delete [] i1;
		delete [] i2;
		delete [] i3;
		delete [] ind1;
	}

	//  Send the dihedral information
	msg->put(numDihedrals);

	if (numDihedrals)
	{
		i1= new int[numDihedrals];
		i2= new int[numDihedrals];
		i3= new int[numDihedrals];
		i4= new int[numDihedrals];
		ind1= new Index[numDihedrals];

		if ( (i1==NULL) || (i2==NULL) || (i3==NULL) || (i4==NULL) || 
		     (ind1==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::send_Molecule");
		}

		for (i=0; i<numDihedrals; i++)
		{
			i1[i]=dihedrals[i].atom1;
			i2[i]=dihedrals[i].atom2;
			i3[i]=dihedrals[i].atom3;
			i4[i]=dihedrals[i].atom4;
			ind1[i]=dihedrals[i].dihedral_type;
		}
	
		msg->put(numDihedrals, i1).put(numDihedrals, i2);
		msg->put(numDihedrals, i3).put(numDihedrals, i4);
		msg->put(numDihedrals, ind1);

		delete [] i1;
		delete [] i2;
		delete [] i3;
		delete [] i4;
		delete [] ind1;
	}

	//  Send the improper information
	msg->put(numImpropers);

	if (numImpropers)
	{
		i1= new int[numImpropers];
		i2= new int[numImpropers];
		i3= new int[numImpropers];
		i4= new int[numImpropers];
		ind1= new Index[numImpropers];

		if ( (i1==NULL) || (i2==NULL) || (i3==NULL) || (i4==NULL) || 
		     (ind1==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::send_Molecule");
		}

		for (i=0; i<numImpropers; i++)
		{
			i1[i]=impropers[i].atom1;
			i2[i]=impropers[i].atom2;
			i3[i]=impropers[i].atom3;
			i4[i]=impropers[i].atom4;
			ind1[i]=impropers[i].improper_type;
		}

		msg->put(numImpropers, i1).put(numImpropers, i2);
		msg->put(numImpropers, i3).put(numImpropers, i4);
		msg->put(numImpropers, ind1);

		delete [] i1;
		delete [] i2;
		delete [] i3;
		delete [] i4;
		delete [] ind1;
	}

	// send the hydrogen bond donor information
	msg->put(numDonors);

	if(numDonors)
	{
	        i1= new int[numDonors];
		i2= new int[numDonors];

		if ( (i1==NULL) || (i2==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::send_Molecule");
		}

		for (i=0; i<numDonors; i++)
		{
			i1[i]=donors[i].atom1;
			i2[i]=donors[i].atom2;
		}

		msg->put(numDonors, i1).put(numDonors, i2);

		delete [] i1;
		delete [] i2;
	}

	// send the hydrogen bond acceptor information
	msg->put(numAcceptors);

	if(numAcceptors)
	{
	        i1= new int[numAcceptors];
		i2= new int[numAcceptors];

		if ( (i1==NULL) || (i2==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::send_Molecule");
		}

		for (i=0; i<numAcceptors; i++)
		{
			i1[i]=acceptors[i].atom1;
			i2[i]=acceptors[i].atom2;
		}

		msg->put(numAcceptors, i1).put(numAcceptors, i2);

		delete [] i1;
		delete [] i2;
	}

	//  Send the exclusion information
	msg->put(numExclusions);

	if (numExclusions)
	{
		i1= new int[numExclusions];
		i2= new int[numExclusions];

		if ( (i1==NULL) || (i2==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::send_Molecule");
		}

		for (i=0; i<numExclusions; i++)
		{
			i1[i]=exclusions[i].atom1;
			i2[i]=exclusions[i].atom2;
		}

		msg->put(numExclusions, i1).put(numExclusions, i2);

		delete [] i1;
		delete [] i2;
	}
	
	//  Send the constraint information, if used
	if (simParams->constraintsOn)
	{
	   msg->put(numConstraints);
	   
	   msg->put(numAtoms, consIndexes);
	   
	   if (numConstraints)
	   {
	      a1 = new Real[numConstraints];
	      v1 = new Vector[numConstraints];
	      
	      if ( (a1 == NULL) || (v1 == NULL) )
	      {
		 NAMD_die("Memory allocation failed in Molecule::send_Molecule()");
	      }
	      
	      for (i=0; i<numConstraints; i++)
	      {
		 a1[i] = consParams[i].k;
		 v1[i] = consParams[i].refPos;
	      }
	      
	      msg->put(numConstraints, a1);
	      msg->put(numConstraints, v1);
	      
	      delete [] a1;
	      delete [] v1;
	   }
	}

	//  Send the langevin parameters, if active
	if (simParams->langevinOn || 
	    simParams->tCoupleOn)
	{
		msg->put(numAtoms, langevinParams);
		msg->put(numAtoms, langForceVals);
	}


	// Broadcast the message to the other nodes
	com_obj->broadcast_others(msg, MOLECULETAG);
	
	//  Now build arrays of indexes into these arrays by atom
	build_lists_by_atom();
}
/*			END OF FUNCTION send_Molecule			*/

/************************************************************************/
/*									*/
/*			FUNCTION receive_Molecule			*/
/*									*/
/*	receive_Molecule is used by all the clients to receive the	*/
/*   structural data sent out by the master node.  It is NEVER called   */
/*   by the Master node.						*/
/*									*/
/************************************************************************/

void Molecule::receive_Molecule(Message *msg)

{
	int *i1, *i2, *i3, *i4;		//  Temporary integer arrays
	Real *a1, *a2;			//  Temporary real arrays
	Index *ind1;			//  Temporary array of Indexes
	Vector *v1;			//  Temporary array of Vectors
	int i;				//  Loop counter

	//  Get the atom information
	msg->get(numAtoms);

	delete [] atoms;
	atoms= new Atom[numAtoms];
	a1   = new Real[numAtoms];
	a2   = new Real[numAtoms];
	ind1 = new Index[numAtoms];
	i1   = new int[numAtoms];

	if (atoms==NULL || a1==NULL || a2==NULL || ind1==NULL || i1==NULL)
	{
	  NAMD_die("memory allocation failed in Molecule::receive_Molecule");
	}

	msg->get(a1);
	msg->get(a2);
	msg->get(ind1);
	msg->get(i1);

	for (i=0; i<numAtoms; i++)
	{
		atoms[i].mass = a1[i];
		atoms[i].charge = a2[i];
		atoms[i].vdw_type = ind1[i];
		atoms[i].status = i1[i];
	}

	delete [] a1;
	delete [] a2;
	delete [] ind1;
	delete [] i1;

	//  Get the bond information
	msg->get(numBonds);

	if (numBonds)
	{
		delete [] bonds;
		bonds=new Bond[numBonds];
		i1 = new int[numBonds];
		i2 = new int[numBonds];
		ind1 = new Index[numBonds];

		if ( (bonds==NULL) || (i1==NULL) || (i2==NULL) || (ind1==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::receive_Molecule");
		}

		msg->get(i1);
		msg->get(i2);
		msg->get(ind1);

		for (i=0; i<numBonds; i++)
		{
			bonds[i].atom1 = i1[i];
			bonds[i].atom2 = i2[i];
			bonds[i].bond_type = ind1[i];
		}

		delete [] i1;
		delete [] i2;
		delete [] ind1;
	}

	//  Get the angle information
	msg->get(numAngles);

	if (numAngles)
	{
		delete [] angles;
		angles=new Angle[numAngles];
		i1 = new int[numAngles];
		i2 = new int[numAngles];
		i3 = new int[numAngles];
		ind1 = new Index[numAngles];

		if ( (angles==NULL) || (i1==NULL) || (i2==NULL) || (i3==NULL) ||
		     (ind1==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::receive_Molecule");
		}

		msg->get(i1);
		msg->get(i2);
		msg->get(i3);
		msg->get(ind1);

		for (i=0; i<numAngles; i++)
		{
			angles[i].atom1 = i1[i];
			angles[i].atom2 = i2[i];
			angles[i].atom3 = i3[i];
			angles[i].angle_type = ind1[i];
		}

		delete [] i1;
		delete [] i2;
		delete [] i3;
		delete [] ind1;
	}

	//  Get the dihedral information
	msg->get(numDihedrals);

	if (numDihedrals)
	{
		delete [] dihedrals;
		dihedrals=new Dihedral[numDihedrals];
		i1 = new int[numDihedrals];
		i2 = new int[numDihedrals];
		i3 = new int[numDihedrals];
		i4 = new int[numDihedrals];
		ind1 = new Index[numDihedrals];

		if ( (dihedrals==NULL) || (i1==NULL) || (i2==NULL) || (i3==NULL) ||
		     (i4==NULL) || (ind1==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::receive_Molecule");
		}

		msg->get(i1);
		msg->get(i2);
		msg->get(i3);
		msg->get(i4);
		msg->get(ind1);

		for (i=0; i<numDihedrals; i++)
		{
			dihedrals[i].atom1 = i1[i];
			dihedrals[i].atom2 = i2[i];
			dihedrals[i].atom3 = i3[i];
			dihedrals[i].atom4 = i4[i];
			dihedrals[i].dihedral_type = ind1[i];
		}

		delete [] i1;
		delete [] i2;
		delete [] i3;
		delete [] i4;
		delete [] ind1;
	}

	//  Get the improper information
	msg->get(numImpropers);

	if (numImpropers)
	{
		delete [] impropers;
		impropers=new Improper[numImpropers];
		i1 = new int[numImpropers];
		i2 = new int[numImpropers];
		i3 = new int[numImpropers];
		i4 = new int[numImpropers];
		ind1 = new Index[numImpropers];

		if ( (impropers==NULL) || (i1==NULL) || (i2==NULL) || (i3==NULL) ||
		     (i4==NULL) || (ind1==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::receive_Molecule");
		}

		msg->get(i1);
		msg->get(i2);
		msg->get(i3);
		msg->get(i4);
		msg->get(ind1);

		for (i=0; i<numImpropers; i++)
		{
			impropers[i].atom1 = i1[i];
			impropers[i].atom2 = i2[i];
			impropers[i].atom3 = i3[i];
			impropers[i].atom4 = i4[i];
			impropers[i].improper_type = ind1[i];
		}

		delete [] i1;
		delete [] i2;
		delete [] i3;
		delete [] i4;
		delete [] ind1;
	}

	//  Get the hydrogen bond donors
	msg->get(numDonors);

	if (numDonors)
	{
		delete [] donors;
	        donors=new Bond[numDonors];
		i1 = new int[numDonors];
		i2 = new int[numDonors];

		if ( (donors==NULL) || (i1==NULL) || (i2==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::receive_Molecule");
		}

		msg->get(i1);
		msg->get(i2);

		for (i=0; i<numDonors; i++)
		{
			donors[i].atom1 = i1[i];
			donors[i].atom2 = i2[i];
		}

		delete [] i1;
		delete [] i2;
	}

	//  Get the hydrogen bond acceptors
	msg->get(numAcceptors);

	if (numAcceptors)
	{
		delete [] acceptors;
	        acceptors=new Bond[numAcceptors];
		i1 = new int[numAcceptors];
		i2 = new int[numAcceptors];

		if ( (acceptors==NULL) || (i1==NULL) || (i2==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::receive_Molecule");
		}

		msg->get(i1);
		msg->get(i2);

		for (i=0; i<numAcceptors; i++)
		{
			acceptors[i].atom1 = i1[i];
			acceptors[i].atom2 = i2[i];
		}

		delete [] i1;
		delete [] i2;
	}

	//  Get the exclusion information
	msg->get(numExclusions);

	if (numExclusions)
	{
		delete [] exclusions;
		exclusions=new Exclusion[numExclusions];
		i1 = new int[numExclusions];
		i2 = new int[numExclusions];

		if ( (exclusions==NULL) || (i1==NULL) || (i2==NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::receive_Molecule");
		}

		msg->get(i1);
		msg->get(i2);

		for (i=0; i<numExclusions; i++)
		{
			exclusions[i].atom1 = i1[i];
			exclusions[i].atom2 = i2[i];
		}

		delete [] i1;
		delete [] i2;
	}
	
	//  Get the constraint information, if they are active
	if (simParams->constraintsOn)
	{
	   msg->get(numConstraints);

	   delete [] consIndexes;
	   consIndexes = new int[numAtoms];
	   
	   if (consIndexes == NULL)
	   {
	      NAMD_die("memory allocation failed in Molecule::receive_Molecule");
	   }
	   
	   msg->get(consIndexes);
	   
	   if (numConstraints)
	   {
	      delete [] consParams;
	      consParams = new ConstraintParams[numConstraints];
	      a1         = new Real[numConstraints];
	      v1         = new Vector[numConstraints];
	      
	      if ( (consParams == NULL) || (a1 == NULL) || (v1 == NULL) )
	      {
		 NAMD_die("memory allocation failed in Molecule::receive_Molecule");
	      }
	      
	      msg->get(a1);
	      msg->get(v1);
	      
	      for (i=0; i<numConstraints; i++)
	      {
		 consParams[i].k = a1[i];
		 consParams[i].refPos = v1[i];
	      }
	      
	      delete [] a1;
	      delete [] v1;
	   }
	}
	      

	//  Get the langevin parameters, if they are active
	if (simParams->langevinOn ||
	    simParams->tCoupleOn)
	{
		delete [] langevinParams;
		delete [] langForceVals;
		langevinParams = new Real[numAtoms];
		langForceVals = new Real[numAtoms];

		if ( (langevinParams == NULL) || (langForceVals == NULL) )
		{
			NAMD_die("memory allocation failed in Molecule::receive_Molecule");
		}

		msg->get(langevinParams);
		msg->get(langForceVals);
	}


	//  Now free the message 
	delete msg;
	
	//  analyze the data and find the status of each atom
	build_atom_status();

	//  Now build arrays of indexes into these arrays by atom
	build_lists_by_atom();
}
/*			END OF FUNCTION receive_Molecule		*/

/************************************************************************/
/*									*/
/*			FUNCTION build_lists_by_atom			*/
/*									*/
/*	This function builds O(NumAtoms) arrays that store the bonds,   */
/*  angles, dihedrals, and impropers, that each atom is involved in.    */
/*  This is a space hog, but VERY fast.  This will certainly have to    */
/*  change to make things scalable in memory, but for now, speed is the */
/*  thing!								*/
/*									*/
/************************************************************************/

void Molecule::build_lists_by_atom()
   
{
   int i;			//  Loop counter
   
   bondsByAtom = new LintList[numAtoms];
   anglesByAtom = new LintList[numAtoms];
   dihedralsByAtom = new LintList[numAtoms];
   impropersByAtom = new LintList[numAtoms];
   exclusionsByAtom = new LintList[numAtoms];

   if ( (bondsByAtom == NULL) || (anglesByAtom == NULL) || (dihedralsByAtom == NULL)
       || (impropersByAtom == NULL) || (exclusionsByAtom == NULL) )
   {
      NAMD_die("memory allocation failed in Molecule::build_lists_by_atom");
   }

   DebugM(3,"Building bond lists.\n");
      
   //  Build the bond lists
   for (i=0; i<numBonds; i++)
   {
      bondsByAtom[bonds[i].atom1].add(i);
      bondsByAtom[bonds[i].atom2].add(i);
   }
   
   DebugM(3,"Building angle lists.\n");
      
   //  Build the angle lists
   for (i=0; i<numAngles; i++)
   {  
      anglesByAtom[angles[i].atom1].add(i);
      anglesByAtom[angles[i].atom2].add(i);
      anglesByAtom[angles[i].atom3].add(i);
   }
   
   DebugM(3,"Building improper lists.\n");
      
   //  Build the improper lists
   for (i=0; i<numImpropers; i++)
   {
      impropersByAtom[impropers[i].atom1].add(i);
      impropersByAtom[impropers[i].atom2].add(i);
      impropersByAtom[impropers[i].atom3].add(i);
      impropersByAtom[impropers[i].atom4].add(i);
   }
   
   DebugM(3,"Building dihedral lists.\n");
      
   //  Build the dihedral lists
   for (i=0; i<numDihedrals; i++)
   {
      dihedralsByAtom[dihedrals[i].atom1].add(i);
      dihedralsByAtom[dihedrals[i].atom2].add(i);
      dihedralsByAtom[dihedrals[i].atom3].add(i);
      dihedralsByAtom[dihedrals[i].atom4].add(i);
   }
   
   DebugM(3,"Building exclusion data.\n");
      
   //  Build the arrays of exclusions for each atom
   build_exclusions();

   if (exclusions != NULL)
   	delete [] exclusions;

   // Now, eliminate 1-4 exclusions which are also fully excluded
   exclusionList.sort();
   int numOriginalExclusions = exclusionList.size();
   for(i=1; i<numOriginalExclusions; ++i)
   {
     if ( exclusionList[i].atom1 == exclusionList[i-1].atom1 &&
	  exclusionList[i].atom2 == exclusionList[i-1].atom2 )
     {
       // modified == 0 < modified == 1 so assign first to second
       // this can't do any harm if they are the same but will
       // make the modified one unmodified (full)
       exclusionList[i].modified = exclusionList[i-1].modified;
     }
   }
   exclusionList.uniq();

   DebugM(3,"Building exclusion lists.\n");
      
   int numTotalExclusions = exclusionList.size();
   exclusions = exclusionList.unencap();
   for (i=0; i<numTotalExclusions; i++)
   {
      exclusionsByAtom[exclusions[i].atom1].add(i);
      exclusionsByAtom[exclusions[i].atom2].add(i);
   }

   DebugM(4,numOriginalExclusions << " exclusions, " <<
	    numTotalExclusions << " unique\n");

}
/*		END OF FUNCTION build_lists_by_atom		*/

/****************************************************************/
/*								*/
/*			FUNCTION build_exclusions		*/
/*								*/
/*	This function builds a list of all the exlcusions       */
/*  atoms.  These lists include explicit exclusions as well as  */
/*  exclusions that are calculated based on the bonded structure*/
/*  and the exclusion flag.  For each pair of atoms that are    */
/*  excluded, the larger of the 2 atom indexes is stored in the */
/*  array of the smaller index.  All the arrays are then sorted.*/
/*  Then to determine if two atoms have an exclusion, a binary  */
/*  search is done on the array of the atom with the smaller    */
/*  index for the larger index.					*/
/*	If the exclusion policy is set to scaled1-4, there are  */
/*  actually two lists built.  One contains the pairs of atoms  */
/*  that are to be exlcuded (i.e., explicit exclusions, 1-2,    */
/*  and 1-3 interactions) and the other contains just the 1-4   */
/*  interactions, since they will need to be determined 	*/
/*  independantly of the other exclusions.			*/
/*								*/
/****************************************************************/

void Molecule::build_exclusions()
{
	int i;					//  Loop counter
	ExclusionSettings exclude_flag;		//  Exclusion policy

	//  Allocate an array of intlists to hold the exclusions for
	//  each atom
	all_exclusions = new IntList[numAtoms];

	if (all_exclusions == NULL)
	{
		NAMD_die("Memory allocation died in Molecule::build_exclusions()");
	}

	exclude_flag = simParams->exclude;

	//  If the exclusion policy is scaled 1-4, then allocate
	//  an array of IntList's to hold the 1-4 interactions
	if (exclude_flag == SCALED14)
	{ 
		onefour_exclusions = new IntList[numAtoms];

		if (onefour_exclusions == NULL)
		{
			NAMD_die("Memory allocation died in Molecule::build_exclusions()");
		}
	}

	//  Go through the explicit exclusions and add them to the arrays
	for (i=0; i<numExclusions; i++)
	{
		if (exclusions[i].atom1 < exclusions[i].atom2)
		{
			all_exclusions[exclusions[i].atom1].add(exclusions[i].atom2);
		}
		else
		{
			all_exclusions[exclusions[i].atom2].add(exclusions[i].atom1);
		}
		exclusionList.add(exclusions[i]);
	}

	//  Now calculate the bonded exlcusions based on the exclusion policy
	switch (exclude_flag)
	{
		case NONE:
			break;
		case ONETWO:
			build12excl(all_exclusions);
			break;
		case ONETHREE:
			build12excl(all_exclusions);
			build13excl(all_exclusions);
			break;
		case ONEFOUR:
			build12excl(all_exclusions);
			build13excl(all_exclusions);
			build14excl(all_exclusions,0);
			break;
		case SCALED14:
			build12excl(all_exclusions);
			build13excl(all_exclusions);
			build14excl(onefour_exclusions,1);
			break;
	}

	//  Sort all the arrays so that we can search them using a binary search
	for (i=0; i<numAtoms; i++)
	{
		all_exclusions[i].sort();

		if (exclude_flag == SCALED14)
		{
			onefour_exclusions[i].sort();
		}
	}
}
/*			END OF FUNCTION build_exclusions		*/

/************************************************************************/
/*									*/
/*			FUNCTION build12excl				*/
/*									*/
/*   INPUTS:								*/
/*	lists - An array of IntList objects to put the exclusions into  */
/*									*/
/*	This function determines all the 1-2 exclusions (that is, atoms */
/*   that are bond together by a linear bond) and places these          */
/*   exclusions into the array of IntList's that are passed in.		*/
/*									*/
/************************************************************************/

void Molecule::build12excl(IntList *lists)
   
{
   int current_val;	//  Current value to check
   int i;		//  Loop counter to loop through all atoms
   
   //  Loop through all the atoms marking the bonded interactions for each one
   for (i=0; i<numAtoms; i++)
   {
   	current_val = bondsByAtom[i].head();
   
	//  Loop through all the bonds for this atom
   	while (current_val != LIST_EMPTY)
   	{
	   if (bonds[current_val].atom1 == i)
	   {
	      if (i<bonds[current_val].atom2)
	      {
		 lists[i].add(bonds[current_val].atom2);
		 exclusionList.add(Exclusion(i,bonds[current_val].atom2));
	      }
	   }
	   else
	   {
	      if (i<bonds[current_val].atom1)
	      {
		 lists[i].add(bonds[current_val].atom1);
		 exclusionList.add(Exclusion(i,bonds[current_val].atom1));
	      }
	   }
	    
     	   current_val = bondsByAtom[i].next();
   	}
   }
}
/*			END OF FUNCTION build12excl			*/

/************************************************************************/
/*									*/
/*			FUNCTION build13excl				*/
/*									*/
/*   INPUTS:								*/
/*	lists - Array of IntList objects to put exclusions into		*/
/*									*/
/*	This function calculates the 1-3 exclusions (that is, atoms     */
/*   that are bonded by linear bonds to a common third atom) and places */
/*   these exclusions into the array of IntList objects passed in.	*/
/*									*/
/************************************************************************/

void Molecule::build13excl(IntList *lists)
   
{
   int bond1, bond2;	//  The two bonds being checked
   int middle_atom;	//  Common third atom
   int i;		//  Loop counter to loop through all atoms
   
   //  Loop through all the atoms looking at the bonded connections
   //  for each one
   for (i=0; i<numAtoms; i++)
   {
   	 bond1 = bondsByAtom[i].head();
   
	 //  Loop through all the bonds directly connect to atom i
   	 while (bond1 != LIST_EMPTY)
    	 {
    	 	if (bonds[bond1].atom1 == i)
        	{
			middle_atom=bonds[bond1].atom2;
      		}
      		else
      		{
			middle_atom=bonds[bond1].atom1;
      		}

      		bond2 = bondsByAtom[middle_atom].head();

		//  Now loop through all the bonds connect to the
		//  middle atom
      		while (bond2 != LIST_EMPTY)
      		{
			if (bonds[bond2].atom1 == middle_atom)
			{
				if (i < bonds[bond2].atom2)
				{
					lists[i].add(bonds[bond2].atom2);
					exclusionList.add(Exclusion(i,bonds[bond2].atom2));
				}
			}
			else
			{
				if (i < bonds[bond2].atom1)
				{
					lists[i].add(bonds[bond2].atom1);
					exclusionList.add(Exclusion(i,bonds[bond2].atom1));
				}
			}

     			bond2 = bondsByAtom[middle_atom].next();
		}

      		bond1 = bondsByAtom[i].next();
	}
   }
}
/*			END OF FUNCTION build13excl			*/

/************************************************************************/
/*									*/
/*				FUNCTION build14excl			*/
/*									*/
/*   INPUTS:								*/
/*	lists - Array of IntList objects to put exclusions into		*/
/*									*/
/*	This function calculates all the 1-4 exclusions (that is,	*/
/*   atoms that are connected via a sequence of three linear bonds) and */
/*   places these interactions into the array of IntList object passed  */
/*   in.								*/
/*									*/
/************************************************************************/


void Molecule::build14excl(IntList *lists, int modified)
   
{
   int bond1, bond2, bond3;	//  The two bonds being checked
   int mid1, mid2;		//  Middle atoms
   int i;			//  Counter to loop through all atoms
   
   //  Loop through all the atoms
   for (i=0; i<numAtoms; i++)
   {	
        // Get all the bonds connect directly to atom i
   	bond1 = bondsByAtom[i].head();
   
  	while (bond1 != LIST_EMPTY)
        {
      		if (bonds[bond1].atom1 == i)
      		{
			mid1=bonds[bond1].atom2;
      		}
      		else
      		{
			mid1=bonds[bond1].atom1;
      		}

      		bond2 = bondsByAtom[mid1].head();

		//  Loop through all the bonds connected to atom mid1
      		while (bond2 != LIST_EMPTY)
      		{
        		if (bonds[bond2].atom1 == mid1)
      			{
				mid2 = bonds[bond2].atom2;
      			}
      			else
      			{
				mid2 = bonds[bond2].atom1;
      			}

			//  Make sure that we don't double back to where
			//  we started from.  This causes strange behavior.
			//  Trust me, I've been there . . .
			if (mid2 == i)
			{
     				bond2 = bondsByAtom[mid1].next();
				continue;
			}

			bond3=bondsByAtom[mid2].head();

			//  Loop through all the bonds connected to mid2
      			while (bond3 != LIST_EMPTY)
      			{
    		   		if (bonds[bond3].atom1 == mid2)
      				{
					//  Make sure that we don't double back to where
					//  we started from.  This causes strange behavior.
					//  Trust me, I've been there . . .
					//  I added this!!!  Why wasn't it there before?  -JCP
					if (bonds[bond3].atom2 != mid1)
					if (i<bonds[bond3].atom2)
					{
					   lists[i].add(bonds[bond3].atom2);
					   exclusionList.add(Exclusion(i,bonds[bond3].atom2,modified));
					}
				}
				else
				{
					//  Make sure that we don't double back to where
					//  we started from.  This causes strange behavior.
					//  Trust me, I've been there . . .
					//  I added this!!!  Why wasn't it there before?  -JCP
					if (bonds[bond3].atom1 != mid1)
				   	if (i<bonds[bond3].atom1)
					{
					   lists[i].add(bonds[bond3].atom1);
					   exclusionList.add(Exclusion(i,bonds[bond3].atom1,modified));
					}
				}

				bond3 = bondsByAtom[mid2].next();
			}

     			bond2 = bondsByAtom[mid1].next();
      		}
	    
      		bond1 = bondsByAtom[i].next();
   	}
   }
}
/*			END OF FUNCTION build14excl			*/

/************************************************************************/
/*									*/
/*			FUNCTION build_constraint_params		*/
/*									*/
/*   INPUTS:								*/
/*	consref - Value of consref parameter from config file		*/
/*	conskfile - Value of conskfile from config file			*/
/*	conskcol - Value of conskcol from config file			*/
/*	initial_pdb - PDB object that contains initial positions	*/
/*	cwd - Current working directory					*/
/*									*/
/*	This function builds all the parameters that are necessary	*/
/*   to do harmonic constraints.  This involves looking through		*/
/*   one or more PDB objects to determine which atoms are constrained,  */
/*   and what the force constant and reference position is force each   */
/*   atom that is constrained.  This information is then stored		*/
/*   in the arrays consIndexes and consParams.				*/
/*									*/
/************************************************************************/

void Molecule::build_constraint_params(StringList *consref, 
				       StringList *conskfile, 
				       StringList *conskcol, 
				       PDB *initial_pdb,
				       char *cwd)
   
{
   PDB *refPDB, *kPDB;		//  Pointer to other PDB's if used
   int i;			//  Loop counter
   int current_index=0;		//  Index into values used
   int kcol;			//  Column to look for force constant in
   Real kval;			//  Force constant value retreived
   char filename[129];		//  PDB filename
   
   //  Get the PDB object that contains the reference positions.  If
   //  the user gave another file name, use it.  Otherwise, just use
   //  the PDB file that has the initial coordinates.  i.e., constrain
   //  the atoms around their initial position.  This is the most likely
   //  case anyway
   if (consref == NULL)
   {
      refPDB = initial_pdb;
   }
   else
   {
      if (consref->next != NULL)
      {
	 NAMD_die("Multiple definitions of constraint reference file in configruation file");
      }

      if ( (cwd == NULL) || (consref->data[0] == '/') )
      {
	   strcpy(filename, consref->data);
      }
      else
      {
	   strcpy(filename, cwd);
	   strcat(filename, consref->data);
      }
      
      refPDB = new PDB(filename);
      if ( refPDB == NULL )
      {
	NAMD_die("Memory allocation failed in Molecule::build_constraint_params");
      }
      
      if (refPDB->num_atoms() != numAtoms)
      {
	 NAMD_die("Number of atoms in constraint reference PDB doesn't match coordinate PDB");
      }
   }
   
   //  Get the PDB to read the force constants from.  Again, if the user
   //  gave us another file name, open that one.  Otherwise, just use
   //  the PDB with the initial coordinates
   if (conskfile == NULL)
   {
      kPDB = initial_pdb;
   }
   else
   {
      if (conskfile->next != NULL)
      {
	 NAMD_die("Multiple definitions of constraint constant file in configuration file");
      }

      if ( (consref != NULL) && (strcasecmp(consref->data, conskfile->data) == 0) )
      {
	 //  Same PDB used for reference positions and force constants
         kPDB = refPDB; 
      }
      else
      {
	if ( (cwd == NULL) || (conskfile->data[0] == '/') )
	{
		strcpy(filename, conskfile->data);
	}
	else
	{
		strcpy(filename, cwd);
		strcat(filename, conskfile->data);
	}

      	kPDB = new PDB(filename);
	if ( kPDB == NULL )
	{
	  NAMD_die("Memory allocation failed in Molecule::build_constraint_params");
	}
      
      	if (kPDB->num_atoms() != numAtoms)
      	{
	   NAMD_die("Number of atoms in constraint constant PDB doesn't match coordinate PDB");
        }
      }
   }
   
   //  Get the column that the force constant is going to be in.  It
   //  can be in any of the 5 floating point fields in the PDB, according
   //  to what the user wants.  The allowable fields are X, Y, Z, O, or
   //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
   //  The default is the 4th field, which is the occupancy
   if (conskcol == NULL)
   {
      kcol = 4;
   }
   else
   {
      if (conskcol->next != NULL)
      {
	 NAMD_die("Multiple definitions of harmonic constraint column in config file");
      }
      
      if (strcasecmp(conskcol->data, "X") == 0)
      {
	 kcol=1;
      }
      else if (strcasecmp(conskcol->data, "Y") == 0)
      {
	 kcol=2;
      }
      else if (strcasecmp(conskcol->data, "Z") == 0)
      {
	 kcol=3;
      }
      else if (strcasecmp(conskcol->data, "O") == 0)
      {
	 kcol=4;
      }
      else if (strcasecmp(conskcol->data, "B") == 0)
      {
	 kcol=5;
      }
      else
      {
	 NAMD_die("conskcol must have value of X, Y, Z, O, or B");
      }
   }
   
   //  Allocate an array that will store an index into the constraint
   //  parameters for each atom.  If the atom is not constrained, its
   //  value will be set to -1 in this array.
   consIndexes = new int[numAtoms];
   
   if (consIndexes == NULL)
   {
      NAMD_die("memory allocation failed in Molecule::build_constraint_params()");
   }
   
   //  Loop through all the atoms and find out which ones are constrained
   for (i=0; i<numAtoms; i++)
   {
      //  Get the k value based on where we were told to find it
      switch (kcol)
      {
	 case 1:
	    kval = (kPDB->atom(i))->xcoor();
	    break;
	 case 2:
	    kval = (kPDB->atom(i))->ycoor();
	    break;
	 case 3:
	    kval = (kPDB->atom(i))->zcoor();
	    break;
	 case 4:
	    kval = (kPDB->atom(i))->occupancy();
	    break;
	 case 5:
	    kval = (kPDB->atom(i))->temperaturefactor();
	    break;
      }
      
      if (kval > 0.0)
      {
	 //  This atom is constrained
	 consIndexes[i] = current_index;
	 current_index++;
      }
      else
      {
	 //  This atom is not constrained
	 consIndexes[i] = -1;
      }
   }
   
   if (current_index == 0)
   {
      //  Constraints were turned on, but there weren't really any constrained
      iout << iWARN << "NO CONSTRAINED ATOMS WERE FOUND, BUT CONSTRAINTS ARE ON . . . " << endi;
   }
   else
   {
      //  Allocate an array to hold the constraint parameters
      consParams = new ConstraintParams[current_index];
      
      if (consParams == NULL)
      {
	 NAMD_die("memory allocation failed in Molecule::build_constraint_params");
      }
   }
   
   numConstraints = current_index;
   
   //  Loop through all the atoms and assign the parameters for those
   //  that are constrained
   for (i=0; i<numAtoms; i++)
   {
      if (consIndexes[i] != -1)
      {
	       //  This atom is constrained, so get the k value again
	       switch (kcol)
               {
	          case 1:
	             consParams[consIndexes[i]].k = (kPDB->atom(i))->xcoor();
	             break;
	          case 2:
	             consParams[consIndexes[i]].k = (kPDB->atom(i))->ycoor();
	             break;
	          case 3:
	             consParams[consIndexes[i]].k = (kPDB->atom(i))->zcoor();
	             break;
	          case 4:
	             consParams[consIndexes[i]].k = (kPDB->atom(i))->occupancy();
	             break;
	          case 5:
	             consParams[consIndexes[i]].k = (kPDB->atom(i))->temperaturefactor();
	             break;
	       }
	       
	       //  Get the reference position
	       consParams[consIndexes[i]].refPos.x = (refPDB->atom(i))->xcoor();
	       consParams[consIndexes[i]].refPos.y = (refPDB->atom(i))->ycoor();
	       consParams[consIndexes[i]].refPos.z = (refPDB->atom(i))->zcoor();
      }
   }
   
   //  If we had to create new PDB objects, delete them now
   if (consref != NULL)
   {
      delete refPDB;
   }
   
   if ((conskfile != NULL) &&
       !((consref != NULL) && 
         (strcasecmp(consref->data, conskfile->data) == 0)
        )
      )
   {
      delete kPDB;
   }
}
/*			END OF FUNCTION build_constraint_params		*/

/************************************************************************/
/*									*/
/*			FUNCTION build_langevin_params			*/
/*									*/
/*   INPUTS:								*/
/*	langfile - Value of langevinfile from config file		*/
/*	langcol - Value of langevincol from config file			*/
/*	initial_pdb - PDB object that contains initial positions	*/
/*      cwd - Current working directory					*/
/*									*/
/*	This function builds the array of b values necessary for	*/
/*   Langevin dynamics.  It takes the name of the PDB file and the      */
/*   column in the PDB file that contains the b values.  It then	*/
/*   builds the array langevinParams for use during the program.	*/
/*									*/
/************************************************************************/

void Molecule::build_langevin_params(StringList *langfile, 
				     StringList *langcol, 
				     PDB *initial_pdb,
				     char *cwd)
   
{
   PDB *bPDB;			//  Pointer to PDB object to use
   int bcol;			//  Column that data is in
   Real bval;			//  b value from PDB file
   int i;			//  Loop counter
   BigReal forceConstant;	//  Constant factor in force calc
   char filename[129];		//  Filename
   
   //  Get the PDB object that contains the b values.  If
   //  the user gave another file name, use it.  Otherwise, just use
   //  the PDB file that has the initial coordinates.  
   if (langfile == NULL)
   {
      bPDB = initial_pdb;
   }
   else
   {
      if (langfile->next != NULL)
      {
	 NAMD_die("Multiple definitions of langvein PDB file in configuration file");
      }

      if ( (cwd == NULL) || (langfile->data[0] == '/') )
      {
	   strcpy(filename, langfile->data);
      }
      else
      {
	   strcpy(filename, cwd);
	   strcat(filename, langfile->data);
      }
      
      bPDB = new PDB(filename);
      if ( bPDB == NULL )
      {
	NAMD_die("Memory allocation failed in Molecule::build_langevin_params");
      }
      
      if (bPDB->num_atoms() != numAtoms)
      {
	 NAMD_die("Number of atoms in langevin parameter PDB doesn't match coordinate PDB");
      }
   }
   
   //  Get the column that the b vaules are in.  It
   //  can be in any of the 5 floating point fields in the PDB, according
   //  to what the user wants.  The allowable fields are X, Y, Z, O, or
   //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
   //  The default is the 4th field, which is the occupancy
   if (langcol == NULL)
   {
      bcol = 4;
   }
   else
   {
      if (langcol->next != NULL)
      {
	 NAMD_die("Multiple definitions of langevin parameter column in config file");
      }
      
      if (strcasecmp(langcol->data, "X") == 0)
      {
	 bcol=1;
      }
      else if (strcasecmp(langcol->data, "Y") == 0)
      {
	 bcol=2;
      }
      else if (strcasecmp(langcol->data, "Z") == 0)
      {
	 bcol=3;
      }
      else if (strcasecmp(langcol->data, "O") == 0)
      {
	 bcol=4;
      }
      else if (strcasecmp(langcol->data, "B") == 0)
      {
	 bcol=5;
      }
      else
      {
	 NAMD_die("langevincol must have value of X, Y, Z, O, or B");
      }
   }
   
   //  Allocate the array to hold all the data
   langevinParams = new Real[numAtoms];
   langForceVals = new Real[numAtoms];
   
   if ( (langevinParams == NULL) || (langForceVals == NULL) )
   {
      NAMD_die("memory allocation failed in Molecule::build_langevin_params()");
   }

   //  Calculate the constant portion of the force values.  Note that
   //  because we need to convert from femtoseconds to picoseconds,
   //  the factor of 0.001 is needed.  
   forceConstant = 0.002*TIMEFACTOR*TIMEFACTOR*BOLTZMAN*(simParams->langevinTemp)/(simParams->dt);
   
   //  Loop through all the atoms and get the b value
   for (i=0; i<numAtoms; i++)
   {
      //  Get the k value based on where we were told to find it
      switch (bcol)
      {
	 case 1:
	    bval = (bPDB->atom(i))->xcoor();
	    break;
	 case 2:
	    bval = (bPDB->atom(i))->ycoor();
	    break;
	 case 3:
	    bval = (bPDB->atom(i))->zcoor();
	    break;
	 case 4:
	    bval = (bPDB->atom(i))->occupancy();
	    break;
	 case 5:
	    bval = (bPDB->atom(i))->temperaturefactor();
	    break;
      }
      
      //  Assign the b value
      langevinParams[i] = bval;

      //  Calculate the random force value
      langForceVals[i] = sqrt(forceConstant*atoms[i].mass*bval);
   }
   
   //  If we had to create a PDB object, delete it now
   if (langfile != NULL)
   {
      delete bPDB;
   }
}
/*			END OF FUNCTION build_constraint_params		*/



Bool Molecule::is_hydrogen(int anum)
{
    return ((atoms[anum].status & HydrogenAtom) != 0);
}

Bool Molecule::is_oxygen(int anum)
{
    return ((atoms[anum].status & OxygenAtom) != 0);
}

Bool Molecule::is_hb_donor(int anum)
{
    return ((atoms[anum].status & HBDonorAtom) != 0);
}

Bool Molecule::is_hb_acceptor(int anum)
{
    return ((atoms[anum].status & HBAcceptorAtom) != 0);
}

Bool Molecule::is_hb_antecedent(int anum)
{
    return ((atoms[anum].status & HBAntecedentAtom) != 0);
}

Bool Molecule::is_hb_hydrogen(int anum)
{
    return ((atoms[anum].status & HBHydrogenAtom) != 0);
}


int Molecule::get_mother_atom(int anum)
{
    // for efficiency reasons, we are not checking if anum is already 
    // hydrogen or not. This function must be called for hydrogens only;
    return atoms[anum].partner;
}


// return the list of hbond donor interactions for the Nth's atom.
IntList *Molecule::get_atom_hb_donors(int anum)
{
    return atoms[anum].donorList;
}


// return the list of hbond acceptor interactions for the Nth's atom.
IntList *Molecule::get_atom_hb_acceptors(int anum)
{
    return atoms[anum].acceptorList;
}


// go through the molecular structure, analyze the status of each atom,
// and save the data in the Atom structures stored for each atom.  This
// could be built up incrementally while the molecule is being read in,
// but doing it all in one shot allows us to just send the basic info
// over the network and have each node calculate the rest of the data on
// it's own.
void Molecule::build_atom_status(void) {
  int i, a1, a2;

  // initialize information for each atom (note that the status has
  // already been initialized during the read/receive phase)
  for (i=0; i < numAtoms; i++) {
    atoms[i].partner = (-1);
    atoms[i].donorList = NULL;
    atoms[i].acceptorList = NULL;
  }

  // find which atom each hydrogen is bound to
  for (i=0; i < numBonds; i++) {
    a1 = bonds[i].atom1;
    a2 = bonds[i].atom2;
    if (is_hydrogen(a1))
      atoms[a1].partner = a2;
    if (is_hydrogen(a2))
      atoms[a2].partner = a1;
  }

  // find the hydrogen bond partners for each atom.  For donors and
  // acceptors, the value stored in donorList is actually the index of
  // the donor-hydrogen pair for the bond.
  for (i=0; i < numDonors; i++) {
    a1 = donors[i].atom1;
    a2 = donors[i].atom2;
    atoms[a1].status |= HBDonorAtom;

    // add this donor-hydrogen pair to the list of hb donor pairs for the atom
    if (atoms[a1].donorList == NULL)
      atoms[a1].donorList = new IntList;
    atoms[a1].donorList->add(i);

    // check for explicit hydrogen in this pair
    if (a2 >= 0)
      atoms[a2].status |= HBHydrogenAtom;
  }

  // look for acceptors and antecedents now in the same fashion.
  for (i=0; i < numAcceptors; i++) {
    a1 = acceptors[i].atom1;
    a2 = acceptors[i].atom2;
    atoms[a1].status |= HBAcceptorAtom;

    // add this acceptor-ante pair to the list of hb acc pairs for the atom
    if (atoms[a1].acceptorList == NULL)
      atoms[a1].acceptorList = new IntList;
    atoms[a1].acceptorList->add(i);

    // check for explicit antecedent in this pair
    if (a2 >= 0)
      atoms[a2].status |= HBAntecedentAtom;
  }

}
