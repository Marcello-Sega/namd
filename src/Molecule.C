/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/


/*
   The class Molecule is used to hold all of the structural information
   for a simulation.  This information is read in from a .psf file and
   cross checked with the Parameters object passed in.  All of the structural
   information is then stored in arrays for use.
*/

#include "Molecule.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <converse.h>
#include "strlib.h"
#include "InfoStream.h"
#include "MStream.h"
#include "Communicate.h"
// #include "Node.h"
#include "ObjectArena.h"
#include "Parameters.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Hydrogen.h"
#include "UniqueSetIter.h"


#define MIN_DEBUG_LEVEL 3
// #define DEBUGM
#include "Debug.h"


class ResidueLookupElem
{
public:
  char mySegid[11];
  ResidueLookupElem *next;	// stored as a linked list
  int firstResid;		// valid resid is >= firstResid
  int lastResid;		// valid resid is <= lastResid
  ResizeArray<int> atomIndex;	// 0-based index for first atom in residue

  ResidueLookupElem(void) { next = 0; firstResid = -1; lastResid = -1; }
  ~ResidueLookupElem(void) { delete next; }
  int lookup(const char *segid, int resid, int *begin, int *end) const;
  ResidueLookupElem* append(const char *segid, int resid, int aid);
};

int ResidueLookupElem::lookup(
	const char *segid, int resid, int *begin, int *end) const {
    const ResidueLookupElem *elem = this;
    int rval = -1;  // error
    while ( elem && strcasecmp(elem->mySegid,segid) ) elem = elem->next;
    if ( elem && (resid >= elem->firstResid) && (resid <= elem->lastResid) ) {
      *begin = elem->atomIndex[resid - elem->firstResid];
      *end = elem->atomIndex[resid - elem->firstResid + 1];
      rval = 0;  // no error
    }
    return rval;
}

ResidueLookupElem* ResidueLookupElem::append(
	const char *segid, int resid, int aid) {
    ResidueLookupElem *rval = this;
    if ( firstResid == -1 ) {  // nothing added yet
      strcpy(mySegid,segid);
      firstResid = resid;
      lastResid = resid;
      atomIndex.add(aid);
      atomIndex.add(aid+1);
    } else if ( ! strcasecmp(mySegid,segid) ) {  // same segid
      if ( resid == lastResid ) {  // same resid
        atomIndex[lastResid - firstResid + 1] = aid + 1;
      } else if ( resid < lastResid ) {  // error
	// We can work around this by creating a new segment.
	iout << iWARN << "Residue " << resid <<
	  " out of order in segment " << segid <<
	  ", lookup for additional residues in this segment disabled.\n" << endi;
	rval = next = new(ResidueLookupElem);
	next->append(segid,resid,aid);
      } else {  // new resid
        for ( ; lastResid < resid; ++lastResid ) atomIndex.add(aid);
        atomIndex[lastResid - firstResid + 1] = aid + 1;
      }
    } else {  // new segid
      rval = next = new(ResidueLookupElem);
      next->append(segid,resid,aid);
    }
    return rval;
}


//  Lookup atom id from segment, residue, and name
int Molecule::get_atom_from_name(
	const char *segid, int resid, const char *aname) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }

  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return -1;
  for ( ; i < end; ++i ) {
    if ( ! strcasecmp(aname,atomNames[i].atomname) ) return i;
  }
  return -1;
}

//  Lookup number of atoms in residue from segment and residue
int Molecule::get_residue_size(
	const char *segid, int resid) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }
  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return 0;
  return ( end - i );
}

//  Lookup atom id from segment, residue, and index in residue
int Molecule::get_atom_from_index_in_residue(
	const char *segid, int resid, int index) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }
  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return -1;
  if ( index >= 0 && index < ( end - i ) ) return ( index + i );
  return -1;
}

/************************************************************************/
/*                  */
/*      FUNCTION initialize  */
/*                  */
/*  This is the initializer for the Molecule class.  It simply sets */
/*  the counts for all the various parameters to 0 and sets the pointers*/
/*  to the arrays that will store these parameters to NULL, since they  */
/*  have not been allocated yet.          */
/*                  */
/************************************************************************/

void Molecule::initialize(SimParameters *simParams, Parameters *param)
{
  if ( sizeof(int32) != 4 ) { NAMD_bug("sizeof(int32) != 4"); }
  this->simParams = simParams;
  this->params = param;
  /*  Initialize array pointers to NULL  */
  atoms=NULL;
  atomNames=NULL;
  resLookup=NULL;
  if ( simParams->globalForcesOn ) {
    resLookup = new ResidueLookupElem;
  }
  bonds=NULL;
  angles=NULL;
  dihedrals=NULL;
  impropers=NULL;
  donors=NULL;
  acceptors=NULL;
  exclusions=NULL;
  tmpArena=NULL;
  bondsWithAtom=NULL;
  bondsByAtom=NULL;
  anglesByAtom=NULL;
  dihedralsByAtom=NULL;
  impropersByAtom=NULL;
  exclusionsByAtom=NULL;
  all_exclusions=NULL;
  langevinParams=NULL;
  langForceVals=NULL;
  fixedAtomFlags=NULL;
  cluster=NULL;
  clusterSize=NULL;
  exPressureAtomFlags=NULL;
  rigidBondLengths=NULL;
  consIndexes=NULL;
  consParams=NULL;
  stirIndexes=NULL;
  stirParams=NULL;
  movDragIndexes=NULL;
  movDragParams=NULL;
  rotDragIndexes=NULL;
  rotDragParams=NULL;
  consTorqueIndexes=NULL;
  consTorqueParams=NULL;
  consForceIndexes=NULL;
  consForce=NULL;
//fepb
  fepAtomFlags=NULL;
//fepe

  nameArena = new ObjectArena<char>;
  // nameArena->setAlignment(8);
  arena = new ObjectArena<int32>;
  // arena->setAlignment(32);
  exclArena = new ObjectArena<char>;
  // exclArena->setAlignment(32);

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
  numStirredAtoms=0;
  numMovDrag=0;
  numRotDrag=0;
  numConsTorque=0;
  numConsForce=0;
  numFixedAtoms=0;
  numFixedGroups=0;
  numExPressureAtoms=0;
  numRigidBonds=0;
  numFixedRigidBonds=0;
  numMultipleDihedrals=0;
  numMultipleImpropers=0;
  numCalcBonds=0;
  numCalcAngles=0;
  numCalcDihedrals=0;
  numCalcImpropers=0;
  numCalcExclusions=0;

//fepb
  numFepInitial = 0;
  numFepFinal = 0;
//fepe

}

/*      END OF FUNCTION initialize */

/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for the Molecule class. */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param)
{
  initialize(simParams,param);
}

/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for the Molecule class from CHARMM/XPLOR files. */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param, char *filename)
{
  initialize(simParams,param);

  read_psf_file(filename, param);
}

/*      END OF FUNCTION Molecule      */

/************************************************************************/
/*                  */
/*        FUNCTION Molecule      */
/*                  */
/*  This is the destructor for the class Molecule.  It simply frees */
/*  the memory allocated for each of the arrays used to store the       */
/*  structure information.            */
/*                  */
/************************************************************************/

Molecule::~Molecule()
{
  /*  Check to see if each array was ever allocated.  If it was   */
  /*  then free it            */
  if (atoms != NULL)
    delete [] atoms;

  if (atomNames != NULL)
  {
    // subarrarys allocated from arena - automatically deleted
    delete [] atomNames;
  }
  delete nameArena;

  if (resLookup != NULL)
    delete resLookup;

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
  
  if (fixedAtomFlags != NULL)
       delete [] fixedAtomFlags;

  if (stirIndexes != NULL)
    delete [] stirIndexes;

  
  if (cluster != NULL)
       delete [] cluster;


  if (clusterSize != NULL)
       delete [] clusterSize;

  if (exPressureAtomFlags != NULL)
       delete [] exPressureAtomFlags;

  if (rigidBondLengths != NULL)
       delete [] rigidBondLengths;

//fepb
  if (fepAtomFlags != NULL)
       delete [] fepAtomFlags;
//fepe

  delete arena;
  delete exclArena;
}
/*      END OF FUNCTION Molecule      */

/************************************************************************/
/*                  */
/*        FUNCTION read_psf_file      */
/*                  */
/*   INPUTS:                */
/*  fname - Name of the .psf file to read        */
/*  params - pointer to Parameters object to use to obtain          */
/*     parameters for vdWs, bonds, etc.      */
/*                  */
/*  This function reads a .psf file in.  This is where just about   */
/*   all of the structural information for this class comes from.  The  */
/*   .psf file contains descriptions of the atom, bonds, angles,        */
/*   dihedrals, impropers, and exclusions.  The parameter object is     */
/*   used to look up parameters for each of these entities.    */
/*                  */
/************************************************************************/

void Molecule::read_psf_file(char *fname, Parameters *params)

{
  char err_msg[512];  //  Error message for NAMD_die
  char buffer[512];  //  Buffer for file reading
  int i;      //  Loop counter
  int NumTitle;    //  Number of Title lines in .psf file
  FILE *psf_file;    //  pointer to .psf file
  int ret_code;    //  ret_code from NAMD_read_line calls

  /* Try and open the .psf file           */
  if ( (psf_file = Fopen(fname, "r")) == NULL)
  {
    sprintf(err_msg, "UNABLE TO OPEN .psf FILE %s", fname);
    NAMD_die(err_msg);
  }

  /*  Read till we have the first non-blank line of file    */
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
  /*  If we can't find it, die.               */
  if (!NAMD_find_word(buffer, "psf"))
  {
    sprintf(err_msg, "UNABLE TO FIND \"PSF\" STRING IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to see if we dropped out of the loop because of a     */
  /*  read error.  This shouldn't happen unless there is nothing  */
  /*  but the PSF line in the file        */
  if (ret_code!=0)
  {
    sprintf(err_msg, "MISSING EVERYTHING BUT PSF FROM %s", fname);
    NAMD_die(err_msg);
  }

  /*  This line should have the word "NTITLE" in it specifying    */
  /*  how many title lines there are        */
  if (!NAMD_find_word(buffer, "NTITLE"))
  {
    sprintf(err_msg,"CAN NOT FIND \"NTITLE\" STRING IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  sscanf(buffer, "%d", &NumTitle);

  /*  Now skip the next NTITLE non-blank lines and then read in the*/
  /*  line which should contain NATOM        */
  i=0;

  while ( ((ret_code=NAMD_read_line(psf_file, buffer)) == 0) && 
    (i<NumTitle) )
  {
    if (!NAMD_blank_string(buffer))
      i++;
  }

  /*  Make sure we didn't exit because of a read error    */
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

  /*  Check to make sure we have the line we want      */
  if (!NAMD_find_word(buffer, "NATOM"))
  {
    sprintf(err_msg, "DIDN'T FIND \"NATOM\" IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  /*  Read in the number of atoms, and then the atoms themselves  */
  sscanf(buffer, "%d", &numAtoms);

  read_atoms(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NBONDS IN PSF");
  }

  /*  Look for the string "NBOND"          */
  if (!NAMD_find_word(buffer, "NBOND"))
  {
    NAMD_die("DID NOT FIND NBOND AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of bonds and then the bonds themselves  */
  sscanf(buffer, "%d", &numBonds);

  if (numBonds)
    read_bonds(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NTHETA IN PSF");
  }

  /*  Look for the string "NTHETA"        */
  if (!NAMD_find_word(buffer, "NTHETA"))
  {
    NAMD_die("DID NOT FIND NTHETA AFTER BOND LIST IN PSF");
  }

  /*  Read in the number of angles and then the angles themselves */
  sscanf(buffer, "%d", &numAngles);

  if (numAngles)
    read_angles(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NPHI IN PSF");
  }

  /*  Look for the string "NPHI"          */
  if (!NAMD_find_word(buffer, "NPHI"))
  {
    NAMD_die("DID NOT FIND NPHI AFTER ANGLE LIST IN PSF");
  }

  /*  Read in the number of dihedrals and then the dihedrals      */
  sscanf(buffer, "%d", &numDihedrals);

  if (numDihedrals)
    read_dihedrals(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NIMPHI IN PSF");
  }

  /*  Look for the string "NIMPHI"        */
  if (!NAMD_find_word(buffer, "NIMPHI"))
  {
    NAMD_die("DID NOT FIND NIMPHI AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of Impropers and then the impropers  */
  sscanf(buffer, "%d", &numImpropers);

  if (numImpropers)
    read_impropers(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NDON IN PSF");
  }

  /*  Look for the string "NDON"        */
  if (!NAMD_find_word(buffer, "NDON"))
  {
    NAMD_die("DID NOT FIND NDON AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of hydrogen bond donors and then the donors */
  sscanf(buffer, "%d", &numDonors);

  if (numDonors)
    read_donors(psf_file);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NACC IN PSF");
  }

  /*  Look for the string "NACC"        */
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
  /*  use this section for anything.        */
  Fclose(psf_file);

  //  analyze the data and find the status of each atom
  build_atom_status();

  return;
}
/*      END OF FUNCTION read_psf_file      */

/************************************************************************/
/*                  */
/*        FUNCTION read_atoms      */
/*                  */
/*   INPUTS:                */
/*  fd - file pointer to the .psf file        */
/*  params - Parameters object to use for parameters    */
/*                  */
/*  this function reads in the Atoms section of the .psf file.      */
/*   This section consists of numAtoms lines that are of the form:  */
/*     <atom#> <mol> <seg#> <res> <atomname> <atomtype> <charge> <mass> */
/*   Each line is read into the appropriate entry in the atoms array.   */
/*   The parameters object is then used to determine the vdW constants  */
/*   for this atom.              */
/*                  */
/************************************************************************/

void Molecule::read_atoms(FILE *fd, Parameters *params)

{
  char buffer[512];  // Buffer for reading from file
  int atom_number=0;  // Atom number 
  int last_atom_number=0; // Last atom number, used to assure
        // atoms are in order
  char segment_name[11]; // Segment name
  char residue_number[11]; // Residue number
  char residue_name[11];  // Residue name
  char atom_name[11];  // Atom name
  char atom_type[11];  // Atom type
  Real charge;    // Charge for the current atom
  Real mass;    // Mass for the current atom
  int read_count;    // Number of fields read by sscanf

  /*  Allocate the atom arrays          */
  atoms     = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];
  hydrogenGroup.resize(0);

  if (atoms == NULL || atomNames == NULL )
  {
    NAMD_die("memory allocation failed in Molecule::read_atoms");
  }

  ResidueLookupElem *tmpResLookup = resLookup;

  /*  Loop and read in numAtoms atom lines.      */
  while (atom_number < numAtoms)
  {
    /*  Get the line from the file        */
    NAMD_read_line(fd, buffer);

    /*  If its blank or a comment, skip it      */
    if ( (NAMD_blank_string(buffer)) || (buffer[0] == '!') )
      continue;

    /*  Parse up the line          */
    read_count=sscanf(buffer, "%d %s %s %s %s %s %f %f",
       &atom_number, segment_name, residue_number,
       residue_name, atom_name, atom_type, &charge, &mass);

    /*  Check to make sure we found what we were expecting  */
    if (read_count != 8)
    {
      char err_msg[128];

      sprintf(err_msg, "BAD ATOM LINE FORMAT IN PSF FILE IN ATOM LINE %d\nLINE=%s",
         last_atom_number+1, buffer);
      NAMD_die(err_msg);
    }

    /*  Check if this is in XPLOR format  */
    int atom_type_num;
    if ( sscanf(atom_type, "%d", &atom_type_num) > 0 )
    {
      NAMD_die("Structure (psf) file is in CHARMM format; XPLOR format required.");
    }

    /*  Make sure the atoms were in sequence    */
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
    /*  for these strings as we really need      */
    int reslength = strlen(residue_name)+1;
    int namelength = strlen(atom_name)+1;
    int typelength = strlen(atom_type)+1;

    atomNames[atom_number-1].resname = nameArena->getNewArray(reslength);
    atomNames[atom_number-1].atomname = nameArena->getNewArray(namelength);
    atomNames[atom_number-1].atomtype = nameArena->getNewArray(typelength);
  
    if (atomNames[atom_number-1].resname == NULL)
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

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append(segment_name, atoi(residue_number), atom_number-1);

    /*  Determine the type of the atom (H or O) */
    if (atoms[atom_number-1].mass <=3.5) {
      atoms[atom_number-1].status |= HydrogenAtom;
    } else if ((atomNames[atom_number-1].atomname[0] == 'O') && 
         (atoms[atom_number-1].mass >= 14.0) && 
         (atoms[atom_number-1].mass <= 18.0)) {
      atoms[atom_number-1].status |= OxygenAtom;
    }

    /*  Look up the vdw constants for this atom    */
    params->assign_vdw_index(atomNames[atom_number-1].atomtype, 
       &(atoms[atom_number-1]));
        }

  return;
}
/*      END OF FUNCTION read_atoms      */

/************************************************************************/
/*                  */
/*      FUNCTION read_bonds        */
/*                  */
/*  read_bonds reads in the bond section of the .psf file.  This    */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are bonded together.  Each atom pair is   */
/*  read in.  Then that parameter object is queried to determine the    */
/*  force constant and rest distance for the bond.      */
/*                  */
/************************************************************************/

void Molecule::read_bonds(FILE *fd, Parameters *params)

{
  int atom_nums[2];  // Atom indexes for the bonded atoms
  char atom1name[11];  // Atom type for atom #1
  char atom2name[11];  // Atom type for atom #2
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
  int origNumBonds = numBonds;   // number of bonds in file header

  /*  Allocate the array to hold the bonds      */
  bonds=new Bond[numBonds];

  if (bonds == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_bonds");
  }

  /*  Loop through and read in all the bonds      */
  while (num_read < numBonds)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
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
    /*  that is alphabetically first as atom 1.    */
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
    /*  this bond            */
    Bond *b = &(bonds[num_read]);
    params->assign_bond_index(atom1name, atom2name, b);

    /*  Assign the atom indexes to the array element  */
    b->atom1=atom_nums[0];
    b->atom2=atom_nums[1];

    /*  Make sure this isn't a fake bond meant for shake in x-plor.  */
    Real k, x0;
    params->get_bond_params(&k,&x0,b->bond_type);
    if ( k == 0. ) --numBonds;  // fake bond
    else ++num_read;  // real bond
  }

  /*  Tell user about our subterfuge  */
  if ( numBonds != origNumBonds ) {
    iout << iWARN << "Ignored " << origNumBonds - numBonds <<
            " bonds with zero force constants.\n" << endi;
    iout << iWARN <<
	"Will get H-H distance in rigid H2O from H-O-H angle.\n" << endi;
  }

  return;
}
/*      END OF FUNCTION read_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION read_angles        */
/*                  */
/*   INPUTS:                */
/*  fd - File descriptor for .psf file        */
/*  params - Parameters object to query for parameters    */
/*                  */
/*  read_angles reads the angle parameters from the .psf file.      */
/*   This section of the .psf file consists of a list of triplets of    */
/*   atom indexes.  Each triplet represents three atoms connected via   */
/*   an angle bond.  The parameter object is queried to obtain the      */
/*   constants for each bond.            */
/*                  */
/************************************************************************/

void Molecule::read_angles(FILE *fd, Parameters *params)

{
  int atom_nums[3];  //  Atom numbers for the three atoms
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  register int j;      //  Loop counter
  int num_read=0;    //  Number of angles read so far
  int origNumAngles = numAngles;  // Number of angles in file
  /*  Alloc the array of angles          */
  angles=new Angle[numAngles];

  if (angles == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_angles");
  }

  /*  Loop through and read all the angles      */
  while (num_read < numAngles)
  {
    /*  Loop through the 3 atom indexes in the current angle*/
    for (j=0; j<3; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "ANGLES")-1;

      /*  Check to make sure the atom index doesn't   */
      /*  exceed the Number of Atoms      */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "ANGLES INDEX %d GREATER THAN NATOM %d IN ANGLES # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }
    }

    /*  Place the bond name that is alphabetically first  */
    /*  in the atom1name.  This is OK since the order of    */
    /*  atom1 and atom3 are interchangable.  And to search  */
    /*  the tree of angle parameters, we need the order     */
    /*  to be predictable.          */
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

    /*  Get the constant values for this bond from the  */
    /*  parameter object          */
    params->assign_angle_index(atom1name, atom2name, 
       atom3name, &(angles[num_read]));

    /*  Assign the three atom indices      */
    angles[num_read].atom1=atom_nums[0];
    angles[num_read].atom2=atom_nums[1];
    angles[num_read].atom3=atom_nums[2];

    /*  Make sure this isn't a fake angle meant for shake in x-plor.  */
    Real k, t0, k_ub, r_ub;
    params->get_angle_params(&k,&t0,&k_ub,&r_ub,angles[num_read].angle_type);
    if ( k == 0. && k_ub == 0. ) --numAngles;  // fake angle
    else ++num_read;  // real angle
  }

  /*  Tell user about our subterfuge  */
  if ( numAngles != origNumAngles ) {
    iout << iWARN << "Ignored " << origNumAngles - numAngles <<
            " angles with zero force constants.\n" << endi;
  }

  return;
}
/*      END OF FUNCTION read_angles      */

/************************************************************************/
/*                  */
/*        FUNCTION read_dihedrals      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for the .psf file        */
/*  params - pointer to parameter object        */
/*                  */
/*  read_dihedreals reads the dihedral section of the .psf file.    */
/*   This section of the file contains a list of quartets of atom       */
/*   numbers.  Each quartet represents a group of atoms that form a     */
/*   dihedral bond.              */
/*                  */
/************************************************************************/

void Molecule::read_dihedrals(FILE *fd, Parameters *params)

{
  int atom_nums[4];  // The 4 atom indexes
  int last_atom_nums[4];  // Atom numbers from previous bond
  char atom1name[11];  // Atom type for atom 1
  char atom2name[11];  // Atom type for atom 2
  char atom3name[11];  // Atom type for atom 3
  char atom4name[11];  // Atom type for atom 4
  register int j;      // loop counter
  int num_read=0;    // number of dihedrals read so far
  int multiplicity=1;  // multiplicity of the current bond
  Bool duplicate_bond;  // Is this a duplicate of the last bond
  int num_unique=0;   // Number of unique dihedral bonds

  //  Initialize the array used to check for duplicate dihedrals
  for (j=0; j<4; j++)
    last_atom_nums[j] = -1;

  /*  Allocate an array to hold the Dihedrals      */
  dihedrals = new Dihedral[numDihedrals];

  if (dihedrals == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_dihedrals");
  }

  /*  Loop through and read all the dihedrals      */
  while (num_read < numDihedrals)
  {
    duplicate_bond = TRUE;

    /*  Loop through and read the 4 indexes for this bond   */
    for (j=0; j<4; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "DIHEDRALS")-1;

      /*  Check for an atom index that is too large  */
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

    /*  Get the atom types for the 4 atoms so we can look  */
    /*  up the constants in the parameter object    */
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

    /*  Get the constants for this dihedral bond    */
    params->assign_dihedral_index(atom1name, atom2name, 
       atom3name, atom4name, &(dihedrals[num_unique-1]),
       multiplicity);

    /*  Assign the atom indexes        */
    dihedrals[num_unique-1].atom1=atom_nums[0];
    dihedrals[num_unique-1].atom2=atom_nums[1];
    dihedrals[num_unique-1].atom3=atom_nums[2];
    dihedrals[num_unique-1].atom4=atom_nums[3];

    num_read++;
  }

  numDihedrals = num_unique;

  return;
}
/*      END OF FUNCTION read_dihedral      */

/************************************************************************/
/*                  */
/*        FUNCTION read_impropers      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*  params - parameter object          */
/*                  */
/*  read_impropers reads the improper section of the .psf file.  */
/*   This section is identical to the dihedral section in that it is    */
/*   made up of a list of quartets of atom indexes that define the      */
/*   atoms that are bonded together.          */
/*                  */
/************************************************************************/

void Molecule::read_impropers(FILE *fd, Parameters *params)

{
  int atom_nums[4];  //  Atom indexes for the 4 atoms
  int last_atom_nums[4];  //  Atom indexes from previous bond
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  char atom4name[11];  //  Atom type for atom 4
  register int j;      //  Loop counter
  int num_read=0;    //  Number of impropers read so far
  int multiplicity=1;  // multiplicity of the current bond
  Bool duplicate_bond;  // Is this a duplicate of the last bond
  int num_unique=0;   // Number of unique dihedral bonds

  //  Initialize the array used to look for duplicate improper
  //  entries.  Set them all to -1 so we know nothing will match
  for (j=0; j<4; j++)
    last_atom_nums[j] = -1;

  /*  Allocate the array to hold the impropers      */
  impropers=new Improper[numImpropers];

  if (impropers == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_impropers");
  }

  /*  Loop through and read all the impropers      */
  while (num_read < numImpropers)
  {
    duplicate_bond = TRUE;

    /*  Loop through the 4 indexes for this improper  */
    for (j=0; j<4; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
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

    /*  Look up the constants for this bond      */
    params->assign_improper_index(atom1name, atom2name, 
       atom3name, atom4name, &(impropers[num_unique-1]),
       multiplicity);

    /*  Assign the atom indexes        */
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
/*      END OF FUNCTION read_impropers      */

/************************************************************************/
/*                  */
/*      FUNCTION read_donors        */
/*                  */
/*  read_donors reads in the bond section of the .psf file.  This   */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are part of an H-bond.  Each atom pair is */
/*  read in.                                                            */
/*                  */
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
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
  int num_no_hydr=0;      // Number of bonds with no hydrogen given

  /*  Allocate the array to hold the bonds      */
  donors=new Bond[numDonors];

  if (donors == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_donors");
  }

  /*  Loop through and read in all the donors      */
  while (num_read < numDonors)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
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

    /*  Assign the atom indexes to the array element  */
    Bond *b = &(donors[num_read]);
    b->atom1=d[0];
    b->atom2=d[1];

    num_read++;
  }

  return;
}
/*      END OF FUNCTION read_donors      */


/************************************************************************/
/*                  */
/*      FUNCTION read_acceptors        */
/*                  */
/*  read_acceptors reads in the bond section of the .psf file.      */
/*  This section contains a list of pairs of numbers where each pair is */
/*  represents two atoms that are part of an H-bond.  Each atom pair is */
/*  read in.                                                            */
/*                  */
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
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
        int num_no_ante=0;      // number of pairs with no antecedent

  /*  Allocate the array to hold the bonds      */
  acceptors=new Bond[numAcceptors];

  if (acceptors == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_acceptors");
  }

  /*  Loop through and read in all the acceptors      */
  while (num_read < numAcceptors)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
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

    /*  Assign the atom indexes to the array element  */
    Bond *b = &(acceptors[num_read]);
    b->atom1=d[0];
    b->atom2=d[1];

    num_read++;
  }

  return;
}
/*      END OF FUNCTION read_acceptors      */


/************************************************************************/
/*                  */
/*      FUNCTION read_exclusions      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*                  */
/*  read_exclusions reads in the explicit non-bonded exclusions     */
/*  from the .psf file.  This section is a little funky, so hang on.    */
/*  Ok, first there is a list of atom indexes that is NumExclusions     */
/*  long.  These are in some sense the atoms that will be exlcuded.     */
/*  Following this list is a list of NumAtoms length that is a list     */
/*  of indexes into the list of excluded atoms.  So an example.  Suppose*/
/*  we have a 5 atom simulation with 3 explicit exclusions.  The .psf   */
/*  file could look like:            */
/*                  */
/*  3!NNB                */
/*  3 4 5                */
/*  0 1 3 3 3              */
/*                  */
/*  This would mean that atom 1 has no explicit exclusions.  Atom 2     */
/*  has an explicit exclusion with atom 3.  Atom 3 has an explicit      */
/*  exclusion with atoms 4 AND 5.  And atoms 4 and 5 have no explicit   */
/*  exclusions.  Got it!?!  I'm not sure who dreamed this up . . .      */
/*                  */
/************************************************************************/

void Molecule::read_exclusions(FILE *fd)

{
  int *exclusion_atoms;  //  Array of indexes of excluded atoms
  register int num_read=0;    //  Number fo exclusions read in
  int current_index;  //  Current index value
  int last_index;    //  the previous index value
  register int insert_index=0;  //  index of where we are in exlcusions array

  /*  Allocate the array of exclusion structures and the array of */
  /*  exlcuded atom indexes          */
  exclusions      = new Exclusion[numExclusions];
  exclusion_atoms = new int[numExclusions];

  if ( (exclusions == NULL) || (exclusion_atoms == NULL) )
  {
    NAMD_die("memory allocation failed in Molecule::read_exclusions");
  }

  /*  First, read in the excluded atoms list      */
  for (num_read=0; num_read<numExclusions; num_read++)
  {
    /*  Read the atom number from the file. Subtract 1 to   */
    /*  convert the index from the 1 to NumAtoms used in the*/
    /*  file to the  0 to NumAtoms-1 that we need    */
    exclusion_atoms[num_read]=NAMD_read_int(fd, "IMPROPERS")-1;

    /*  Check for an illegal index        */
    if (exclusion_atoms[num_read] >= numAtoms)
    {
      char err_msg[128];

      sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN PSF FILE", exclusion_atoms[num_read]+1, numAtoms, num_read+1);
      NAMD_die(err_msg);
    }
  }

  /*  Now, go through and read the list of NumAtoms pointers into */
  /*  the array that we just read in        */
  last_index=0;

  for (num_read=0; num_read<numAtoms; num_read++)
  {
    /*  Read in the current index value      */
    current_index=NAMD_read_int(fd, "EXCLUSIONS");

    /*  Check for an illegal pointer      */
    if (current_index>numExclusions)
    {
      char err_msg[128];

      sprintf(err_msg, "EXCLUSION INDEX %d LARGER THAN NUMBER OF EXLCUSIONS %d IN PSF FILE, EXCLUSION #%d\n", 
         current_index+1, numExclusions, num_read);
      NAMD_die(err_msg);
    }

    /*  Check to see if it matches the last index.  If so   */
    /*  than this atom has no exclusions.  If not, then     */
    /*  we have to build some exclusions      */
    if (current_index != last_index)
    {
      /*  This atom has some exlcusions.  Loop from   */
      /*  the last_index to the current index.  This  */
      /*  will include how ever many exclusions this  */
      /*  atom has          */
      for (insert_index=last_index; 
           insert_index<current_index; insert_index++)
      {
        /*  Assign the two atoms involved.      */
        /*  The first one is our position in    */
        /*  the list, the second is based on    */
        /*  the pointer into the index list     */
        int a1 = num_read;
        int a2 = exclusion_atoms[insert_index];
        if ( a1 < a2 ) {
          exclusions[insert_index].atom1 = a1;
          exclusions[insert_index].atom2 = a2;
        } else if ( a2 < a1 ) {
          exclusions[insert_index].atom1 = a2;
          exclusions[insert_index].atom2 = a1;
        } else {
          char err_msg[128];
          sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN PSF FILE\n", a1+1);
          NAMD_die(err_msg);
        }
      }

      last_index=current_index;
    }
  }

  /*  Free our temporary list of indexes        */
  delete [] exclusion_atoms;

  return;
}
/*      END OF FUNCTION read_exclusions      */

/************************************************************************/
/*                  */
/*      FUNCTION print_atoms        */
/*                  */
/*  print_atoms prints out the list of atoms stored in this object. */
/*  It is inteded mainly for debugging purposes.      */
/*                  */
/************************************************************************/

void Molecule::print_atoms(Parameters *params)

{
  register int i;
  Real sigma;
  Real epsilon;
  Real sigma14;
  Real epsilon14;

  DebugM(2,"ATOM LIST\n" \
      << "******************************************\n" \
                  << "NUM  NAME TYPE RES  MASS    CHARGE CHARGE   FEP-CHARGE"  \
      << "SIGMA   EPSILON SIGMA14 EPSILON14\n" \
        << endi);

  for (i=0; i<numAtoms; i++)
  {
    params->get_vdw_params(&sigma, &epsilon, &sigma14, &epsilon14, 
        atoms[i].vdw_type);

    DebugM(2,i+1 << " " << atomNames[i].atomname  \
              << " " << atomNames[i].atomtype << " " \
              << atomNames[i].resname  << " " << atoms[i].mass  \
        << " " << atoms[i].charge << " " << sigma \
        << " " << epsilon << " " << sigma14 \
        << " " << epsilon14 << "\n" \
        << endi);
  }
}
/*      END OF FUNCTION print_atoms      */

/************************************************************************/
/*                  */
/*      FUNCTION print_bonds        */
/*                  */
/*  print_bonds prints out the list of bonds stored in this object. */
/*  It is inteded mainly for debugging purposes.      */
/*                  */
/************************************************************************/

void Molecule::print_bonds(Parameters *params)

{
  register int i;
  Real k;
  Real x0;

  DebugM(2,"BOND LIST\n" << "********************************\n" \
      << "ATOM1 ATOM2 TYPE1 TYPE2      k        x0" \
      << endi);

  for (i=0; i<numBonds; i++)
  {
    params->get_bond_params(&k, &x0, bonds[i].bond_type);

    DebugM(2,bonds[i].atom1+1 << " " \
       << bonds[i].atom2+1 << " "   \
       << atomNames[bonds[i].atom1].atomtype << " "  \
       << atomNames[bonds[i].atom2].atomtype << " " << k \
       << " " << x0 << endi);
  }
}
/*      END OF FUNCTION print_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION print_exclusions      */
/*                  */
/*  print_exlcusions prints out the list of exlcusions stored in    */
/*  this object.  It is inteded mainly for debugging purposes.    */
/*                  */
/************************************************************************/

void Molecule::print_exclusions()
{
  register int i;

  DebugM(2,"EXPLICIT EXCLUSION LIST\n" \
      << "********************************\n" \
            << "ATOM1 ATOM2 " \
      << endi);

  for (i=0; i<numExclusions; i++)
  {
    DebugM(2,exclusions[i].atom1+1 << "  " \
       << exclusions[i].atom2+1 << endi);
  }
}
/*      END OF FUNCTION print_exclusions    */

/************************************************************************/
/*                  */
/*      FUNCTION send_Molecule        */
/*                  */
/*  send_Molecule is used by the Master node to distribute the      */
/*   structural information to all the client nodes.  It is NEVER called*/
/*   by the client nodes.              */
/*                  */
/************************************************************************/


void Molecule::send_Molecule(Communicate *com_obj)

{
  //  Message to send to clients
  MOStream *msg=com_obj->newOutputStream(ALLBUTME, MOLECULETAG, BUFSIZE);
  if ( msg == NULL )
  {
    NAMD_die("Memory allocation failed in Molecule::send_Molecule");
  }

      msg->put(numAtoms);
      msg->put(numAtoms*sizeof(Atom), (char*)atoms);

      //  Send the bond information
      msg->put(numBonds);

      if (numBonds)
      {
        msg->put(numBonds*sizeof(Bond), (char*)bonds);
      }

      //  Send the angle information
      msg->put(numAngles);

      if (numAngles)
      {
        msg->put(numAngles*sizeof(Angle), (char*)angles);
      }

      //  Send the dihedral information
      msg->put(numDihedrals);

      if (numDihedrals)
      {
        msg->put(numDihedrals*sizeof(Dihedral), (char*)dihedrals);
      }

      //  Send the improper information
      msg->put(numImpropers);

      if (numImpropers)
      {
        msg->put(numImpropers*sizeof(Improper), (char*)impropers);
      }

      // send the hydrogen bond donor information
      msg->put(numDonors);

      if(numDonors)
      {
        msg->put(numDonors*sizeof(Bond), (char*)donors);
      }

      // send the hydrogen bond acceptor information
      msg->put(numAcceptors);

      if(numAcceptors)
      {
        msg->put(numAcceptors*sizeof(Bond), (char*)acceptors);
      }

      //  Send the exclusion information
      msg->put(numExclusions);

      if (numExclusions)
      {
        msg->put(numExclusions*sizeof(Exclusion), (char*)exclusions);
      }
      
      //  Send the constraint information, if used
      if (simParams->constraintsOn)
      {
         msg->put(numConstraints);
         
         msg->put(numAtoms, consIndexes);
         
         if (numConstraints)
         {
           msg->put(numConstraints*sizeof(ConstraintParams), (char*)consParams);
         }
      }
      
      //  Send the stirring information, if used
      if (simParams->stirOn)
      {
	       //CkPrintf ("DEBUG: putting numStirredAtoms..\n");
         msg->put(numStirredAtoms);
	       //CkPrintf ("DEBUG: putting numAtoms,stirIndexes.. numAtoms=%d\n",numStirredAtoms);
         msg->put(numAtoms, stirIndexes);
	       //CkPrintf ("DEBUG: if numStirredAtoms..\n");
         if (numStirredAtoms)
         {
           //CkPrintf ("DEBUG: big put, with (char*)stirParams\n");
           msg->put(numStirredAtoms*sizeof(StirParams), (char*)stirParams);
         }
      }
      
      
      //  Send the moving drag information, if used
      if (simParams->movDragOn) {
         msg->put(numMovDrag);
         msg->put(numAtoms, movDragIndexes);
         if (numMovDrag)
         {
           msg->put(numMovDrag*sizeof(MovDragParams), (char*)movDragParams);
         }
      }
      
      //  Send the rotating drag information, if used
      if (simParams->rotDragOn) {
         msg->put(numRotDrag);
         msg->put(numAtoms, rotDragIndexes);
         if (numRotDrag)
         {
           msg->put(numRotDrag*sizeof(RotDragParams), (char*)rotDragParams);
         }
      }
      
      //  Send the "constant" torque information, if used
      if (simParams->consTorqueOn) {
         msg->put(numConsTorque);
         msg->put(numAtoms, consTorqueIndexes);
         if (numConsTorque)
         {
           msg->put(numConsTorque*sizeof(ConsTorqueParams), (char*)consTorqueParams);
         }
      }
      
      // Send the constant force information, if used
      if (simParams->consForceOn)
      { msg->put(numConsForce);
        msg->put(numAtoms, consForceIndexes);
        if (numConsForce)
          msg->put(numConsForce*sizeof(Vector), (char*)consForce);
      }

      //  Send the langevin parameters, if active
      if (simParams->langevinOn || simParams->tCoupleOn)
      {
        msg->put(numAtoms, langevinParams);
        msg->put(numAtoms, langForceVals);
      }

      //  Send fixed atoms, if active
      if (simParams->fixedAtomsOn)
      {
        msg->put(numFixedAtoms);
        msg->put(numAtoms, fixedAtomFlags);
	msg->put(numFixedRigidBonds);
      }

      if (simParams->excludeFromPressure) {
        msg->put(numExPressureAtoms);
        msg->put(numAtoms, exPressureAtomFlags);
      }

//fepb
      // send fep atom info
      if (simParams->fepOn || simParams->lesOn || simParams->pairInteractionOn) {
        msg->put(numFepInitial);
        msg->put(numFepFinal);
        msg->put(numAtoms*sizeof(char), (char*)fepAtomFlags);
      }
//fepe

      // Broadcast the message to the other nodes
      msg->end();
      delete msg;
      
      //  Now build arrays of indexes into these arrays by atom
      build_lists_by_atom();
}
    /*      END OF FUNCTION send_Molecule      */

    /************************************************************************/
    /*                  */
    /*      FUNCTION receive_Molecule      */
    /*                  */
    /*  receive_Molecule is used by all the clients to receive the  */
    /*   structural data sent out by the master node.  It is NEVER called   */
    /*   by the Master node.            */
    /*                  */
    /************************************************************************/

void Molecule::receive_Molecule(MIStream *msg)
{
      //  Get the atom information
      msg->get(numAtoms);

      delete [] atoms;
      atoms= new Atom[numAtoms];

      msg->get(numAtoms*sizeof(Atom), (char*)atoms);

      //  Get the bond information
      msg->get(numBonds);

      if (numBonds)
      {
        delete [] bonds;
        bonds=new Bond[numBonds];

        msg->get(numBonds*sizeof(Bond), (char*)bonds);
      }

      //  Get the angle information
      msg->get(numAngles);

      if (numAngles)
      {
        delete [] angles;
        angles=new Angle[numAngles];

        msg->get(numAngles*sizeof(Angle), (char*)angles);
      }

      //  Get the dihedral information
      msg->get(numDihedrals);

      if (numDihedrals)
      {
        delete [] dihedrals;
        dihedrals=new Dihedral[numDihedrals];

        msg->get(numDihedrals*sizeof(Dihedral), (char*)dihedrals);
      }

      //  Get the improper information
      msg->get(numImpropers);

      if (numImpropers)
      {
        delete [] impropers;
        impropers=new Improper[numImpropers];

        msg->get(numImpropers*sizeof(Improper), (char*)impropers);
      }

      //  Get the hydrogen bond donors
      msg->get(numDonors);

      if (numDonors)
      {
        delete [] donors;
        donors=new Bond[numDonors];

        msg->get(numDonors*sizeof(Bond), (char*)donors);
      }

      //  Get the hydrogen bond acceptors
      msg->get(numAcceptors);

      if (numAcceptors)
      {
        delete [] acceptors;
        acceptors=new Bond[numAcceptors];

        msg->get(numAcceptors*sizeof(Bond), (char*)acceptors);
      }

      //  Get the exclusion information
      msg->get(numExclusions);

      if (numExclusions)
      {
        delete [] exclusions;
        exclusions=new Exclusion[numExclusions];

        msg->get(numExclusions*sizeof(Exclusion), (char*)exclusions);
      }
      
      //  Get the constraint information, if they are active
      if (simParams->constraintsOn)
      {
         msg->get(numConstraints);

         delete [] consIndexes;
         consIndexes = new int32[numAtoms];
         
         msg->get(numAtoms, consIndexes);
         
         if (numConstraints)
         {
           delete [] consParams;
           consParams = new ConstraintParams[numConstraints];
      
           msg->get(numConstraints*sizeof(ConstraintParams), (char*)consParams);
         }
      }

            
      //  Get the stirring information, if stirring is  active
      if (simParams->stirOn)
      {
         msg->get(numStirredAtoms);

         delete [] stirIndexes;
         stirIndexes = new int32[numAtoms];
         
         msg->get(numAtoms, stirIndexes);
         
         if (numStirredAtoms)
         {
           delete [] stirParams;
           stirParams = new StirParams[numStirredAtoms];
      
           msg->get(numStirredAtoms*sizeof(StirParams), (char*)stirParams);
         }
      }
      
      //  Get the moving drag information, if it is active
      if (simParams->movDragOn) {
         msg->get(numMovDrag);
         delete [] movDragIndexes;
         movDragIndexes = new int32[numAtoms];
         msg->get(numAtoms, movDragIndexes);
         if (numMovDrag)
         {
           delete [] movDragParams;
           movDragParams = new MovDragParams[numMovDrag];
           msg->get(numMovDrag*sizeof(MovDragParams), (char*)movDragParams);
         }
      }
      
      //  Get the rotating drag information, if it is active
      if (simParams->rotDragOn) {
         msg->get(numRotDrag);
         delete [] rotDragIndexes;
         rotDragIndexes = new int32[numAtoms];
         msg->get(numAtoms, rotDragIndexes);
         if (numRotDrag)
         {
           delete [] rotDragParams;
           rotDragParams = new RotDragParams[numRotDrag];
           msg->get(numRotDrag*sizeof(RotDragParams), (char*)rotDragParams);
         }
      }
      
      //  Get the "constant" torque information, if it is active
      if (simParams->consTorqueOn) {
         msg->get(numConsTorque);
         delete [] consTorqueIndexes;
         consTorqueIndexes = new int32[numAtoms];
         msg->get(numAtoms, consTorqueIndexes);
         if (numConsTorque)
         {
           delete [] consTorqueParams;
           consTorqueParams = new ConsTorqueParams[numConsTorque];
           msg->get(numConsTorque*sizeof(ConsTorqueParams), (char*)consTorqueParams);
         }
      }
      
      // Get the constant force information, if it's active
      if (simParams->consForceOn)
      { msg->get(numConsForce);
        delete [] consForceIndexes;
        consForceIndexes = new int32[numAtoms];
        msg->get(numAtoms, consForceIndexes);
        if (numConsForce)
        { delete [] consForce;
          consForce = new Vector[numConsForce];
          msg->get(numConsForce*sizeof(Vector), (char*)consForce);
        }
      }

      //  Get the langevin parameters, if they are active
      if (simParams->langevinOn || simParams->tCoupleOn)
      {
        delete [] langevinParams;
        delete [] langForceVals;
        langevinParams = new Real[numAtoms];
        langForceVals = new Real[numAtoms];

        msg->get(numAtoms, langevinParams);
        msg->get(numAtoms, langForceVals);
      }

      //  Get the fixed atoms, if they are active
      if (simParams->fixedAtomsOn)
      {
        delete [] fixedAtomFlags;
        fixedAtomFlags = new int32[numAtoms];

        msg->get(numFixedAtoms);
        msg->get(numAtoms, fixedAtomFlags);
        msg->get(numFixedRigidBonds);
      }

      if (simParams->excludeFromPressure) {
        exPressureAtomFlags = new int32[numAtoms];
        msg->get(numExPressureAtoms);
        msg->get(numAtoms, exPressureAtomFlags);
      }

//fepb
      //receive fep atom info
      if (simParams->fepOn || simParams->lesOn || simParams->pairInteractionOn) {
        delete [] fepAtomFlags;
        fepAtomFlags = new unsigned char[numAtoms];

        msg->get(numFepInitial);
        msg->get(numFepFinal);
        msg->get(numAtoms*sizeof(unsigned char), (char*)fepAtomFlags);
      }
//fepe

      //  Now free the message 
      delete msg;
      
      //  analyze the data and find the status of each atom
      build_atom_status();

      //  Now build arrays of indexes into these arrays by atom
      build_lists_by_atom();
    }
    /*      END OF FUNCTION receive_Molecule    */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_lists_by_atom      */
    /*                  */
    /*  This function builds O(NumAtoms) arrays that store the bonds,   */
    /*  angles, dihedrals, and impropers, that each atom is involved in.    */
    /*  This is a space hog, but VERY fast.  This will certainly have to    */
    /*  change to make things scalable in memory, but for now, speed is the */
    /*  thing!                */
    /*                  */
    /************************************************************************/

    void Molecule::build_lists_by_atom()
       
    {
       register int i;      //  Loop counter

       register int numFixedAtoms = this->numFixedAtoms;
       // if we want forces on fixed atoms then just pretend
       // there are none for the purposes of this routine;
       if ( simParams->fixedAtomsForces ) numFixedAtoms = 0;

//fepb
//     int numFepInitial = this->numFepInitial;
//     int numFepFinal = this->numFepFinal;
//fepe
       tmpArena = new ObjectArena<int32>;
       bondsWithAtom = new int32 *[numAtoms];
       cluster = new int32 [numAtoms];
       clusterSize = new int32 [numAtoms];
       bondsByAtom = new int32 *[numAtoms];
       anglesByAtom = new int32 *[numAtoms];
       dihedralsByAtom = new int32 *[numAtoms];
       impropersByAtom = new int32 *[numAtoms];
       exclusionsByAtom = new int32 *[numAtoms];

       int32 *byAtomSize = new int32[numAtoms];

       const int pair_self = 
         simParams->pairInteractionOn ? simParams->pairInteractionSelf : 0;

       DebugM(3,"Building bond lists.\n");
    
       //  Build the bond lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       for (i=0; i<numBonds; i++)
       {
         byAtomSize[bonds[i].atom1]++;
         byAtomSize[bonds[i].atom2]++;
       }
       for (i=0; i<numAtoms; i++)
       {
         bondsWithAtom[i] = tmpArena->getNewArray(byAtomSize[i]+1);
         bondsWithAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numBonds; i++)
       {
         int a1 = bonds[i].atom1;
         int a2 = bonds[i].atom2;
         bondsWithAtom[a1][byAtomSize[a1]++] = i;
         bondsWithAtom[a2][byAtomSize[a2]++] = i;
       }

       //  Build cluster information (contiguous clusters)
       for (i=0; i<numAtoms; i++) {
         cluster[i] = i;
       }
       for (i=0; i<numAtoms; i++) {
         int ci = i;
         while ( cluster[ci] != ci ) ci = cluster[ci];
         for ( int32 *b = bondsWithAtom[i]; *b != -1; ++b ) {
           int a = bonds[*b].atom1;
           if ( a == i ) a = bonds[*b].atom2;
           if ( a > i ) {
             int ca = a;
             while ( cluster[ca] != ca ) ca = cluster[ca];
             if ( ca > ci ) cluster[ca] = cluster[ci];
             else cluster[ci] = cluster[ca];
           }
         }
       }
       while ( 1 ) {
         int allok = 1;
         for (i=0; i<numAtoms; i++) {
           int ci = cluster[i];
           if ( cluster[ci] != ci ) {
             allok = 0;
             cluster[i] = cluster[ci];
           }
         }
         if ( allok ) break;
       }
       for (i=0; i<numAtoms; i++) {
         clusterSize[i] = 0;
       }
       for (i=0; i<numAtoms; i++) {
         clusterSize[cluster[i]] += 1;
       }
       
       //  Build the bond lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcBonds = 0;
       for (i=0; i<numBonds; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[bonds[i].atom1]
                            && fixedAtomFlags[bonds[i].atom2] ) continue;
   
         if ( pair_self && fepAtomFlags[bonds[i].atom1] != 1) continue;
         byAtomSize[bonds[i].atom1]++;
         numCalcBonds++;
       }
       for (i=0; i<numAtoms; i++)
       {
         bondsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         bondsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numBonds; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[bonds[i].atom1]
                            && fixedAtomFlags[bonds[i].atom2] ) continue;
         if ( pair_self && fepAtomFlags[bonds[i].atom1] != 1) continue;
         int a1 = bonds[i].atom1;
         bondsByAtom[a1][byAtomSize[a1]++] = i;
       }
       
       DebugM(3,"Building angle lists.\n");
    
       //  Build the angle lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcAngles = 0;
       for (i=0; i<numAngles; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[angles[i].atom1]
                            && fixedAtomFlags[angles[i].atom2]
                            && fixedAtomFlags[angles[i].atom3] ) continue;
         if ( pair_self && fepAtomFlags[angles[i].atom1] != 1) continue;
         byAtomSize[angles[i].atom1]++;
         numCalcAngles++;
       }
       for (i=0; i<numAtoms; i++)
       {
         anglesByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         anglesByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numAngles; i++)
       {
         if ( pair_self && fepAtomFlags[angles[i].atom1] != 1) continue;
         if ( numFixedAtoms && fixedAtomFlags[angles[i].atom1]
                            && fixedAtomFlags[angles[i].atom2]
                            && fixedAtomFlags[angles[i].atom3] ) continue;
         int a1 = angles[i].atom1;
         anglesByAtom[a1][byAtomSize[a1]++] = i;
       }
       
       DebugM(3,"Building improper lists.\n");
    
       //  Build the improper lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcImpropers = 0;
       for (i=0; i<numImpropers; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[impropers[i].atom1]
                            && fixedAtomFlags[impropers[i].atom2]
                            && fixedAtomFlags[impropers[i].atom3]
                            && fixedAtomFlags[impropers[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[impropers[i].atom1] != 1) continue;
         byAtomSize[impropers[i].atom1]++;
         numCalcImpropers++;
       }
       for (i=0; i<numAtoms; i++)
       {
         impropersByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         impropersByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numImpropers; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[impropers[i].atom1]
                            && fixedAtomFlags[impropers[i].atom2]
                            && fixedAtomFlags[impropers[i].atom3]
                            && fixedAtomFlags[impropers[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[impropers[i].atom1] != 1) continue;
         int a1 = impropers[i].atom1;
         impropersByAtom[a1][byAtomSize[a1]++] = i;
       }
       
       DebugM(3,"Building dihedral lists.\n");
    
       //  Build the dihedral lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcDihedrals = 0;
       for (i=0; i<numDihedrals; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[dihedrals[i].atom1]
                            && fixedAtomFlags[dihedrals[i].atom2]
                            && fixedAtomFlags[dihedrals[i].atom3]
                            && fixedAtomFlags[dihedrals[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[dihedrals[i].atom1] != 1) continue;
         byAtomSize[dihedrals[i].atom1]++;
         numCalcDihedrals++;
       }
       for (i=0; i<numAtoms; i++)
       {
         dihedralsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         dihedralsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numDihedrals; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[dihedrals[i].atom1]
                            && fixedAtomFlags[dihedrals[i].atom2]
                            && fixedAtomFlags[dihedrals[i].atom3]
                            && fixedAtomFlags[dihedrals[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[dihedrals[i].atom1] != 1) continue;
         int a1 = dihedrals[i].atom1;
         dihedralsByAtom[a1][byAtomSize[a1]++] = i;
       }
    
       DebugM(3,"Building exclusion data.\n");
    
       //  Build the arrays of exclusions for each atom
       build_exclusions();

       //  Remove temporary structures
       delete [] bondsWithAtom;  bondsWithAtom = 0;
       delete tmpArena;  tmpArena = 0;

       if (exclusions != NULL)
      delete [] exclusions;

       // 1-4 exclusions which are also fully excluded were eliminated by hash table
       numTotalExclusions = exclusionSet.size();
       exclusions = new Exclusion[numTotalExclusions];
       UniqueSetIter<Exclusion> exclIter(exclusionSet);
       for ( exclIter=exclIter.begin(),i=0; exclIter != exclIter.end(); exclIter++,i++ )
       {
         exclusions[i] = *exclIter;
       }
       // Free exclusionSet storage
       // exclusionSet.clear(1);
       exclusionSet.clear();

       DebugM(3,"Building exclusion lists.\n");
    
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcExclusions = 0;
       for (i=0; i<numTotalExclusions; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
         byAtomSize[exclusions[i].atom1]++;
         numCalcExclusions++;
       }
       for (i=0; i<numAtoms; i++)
       {
         exclusionsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         exclusionsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
         int a1 = exclusions[i].atom1;
         exclusionsByAtom[a1][byAtomSize[a1]++] = i;
       }

       delete [] byAtomSize;  byAtomSize = 0;


       //  Allocate an array to hold the exclusions for each atom
       all_exclusions = new ExclusionCheck[numAtoms];

       for (i=0; i<numAtoms; i++)
       {
         all_exclusions[i].min = numAtoms;
         all_exclusions[i].max = -1;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         // first atom should alway have lower number!
         int a1 = exclusions[i].atom1;
         int a2 = exclusions[i].atom2;
         if ( numFixedAtoms && fixedAtomFlags[a1]
                            && fixedAtomFlags[a2] ) continue;
         if ( all_exclusions[a1].min > a2 ) all_exclusions[a1].min = a2;
         if ( all_exclusions[a2].min > a1 ) all_exclusions[a2].min = a1;
         if ( all_exclusions[a1].max < a2 ) all_exclusions[a1].max = a2;
         if ( all_exclusions[a2].max < a1 ) all_exclusions[a2].max = a1;
       }
       int exclmem = 0;
       for (i=0; i<numAtoms; i++)
       {
         if ( all_exclusions[i].max != -1 ) {
           int s = all_exclusions[i].max - all_exclusions[i].min + 1;
           char *f = all_exclusions[i].flags = exclArena->getNewArray(s);
           for ( int k=0; k<s; ++k ) f[k] = 0;
           exclmem += s;
         } else {
           all_exclusions[i].flags = 0;
         }
       }
       if ( 0 ) {
         iout << iINFO << numTotalExclusions << " exclusions consume "
            << exclmem << " bytes.\n" << endi;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         int a1 = exclusions[i].atom1;
         int a2 = exclusions[i].atom2;
         if ( numFixedAtoms && fixedAtomFlags[a1]
                            && fixedAtomFlags[a2] ) continue;
         if ( exclusions[i].modified ) {
           all_exclusions[a1].flags[a2-all_exclusions[a1].min] = EXCHCK_MOD;
           all_exclusions[a2].flags[a1-all_exclusions[a2].min] = EXCHCK_MOD;
         } else {
           all_exclusions[a1].flags[a2-all_exclusions[a1].min] = EXCHCK_FULL;
           all_exclusions[a2].flags[a1-all_exclusions[a2].min] = EXCHCK_FULL;
         }
       }

    }
    /*    END OF FUNCTION build_lists_by_atom    */

    /****************************************************************/
    /*                */
    /*      FUNCTION build_exclusions    */
    /*                */
    /*  This function builds a list of all the exlcusions       */
    /*  atoms.  These lists include explicit exclusions as well as  */
    /*  exclusions that are calculated based on the bonded structure*/
    /*  and the exclusion flag.  For each pair of atoms that are    */
    /*  excluded, the larger of the 2 atom indexes is stored in the */
    /*  array of the smaller index.  All the arrays are not sorted. */
    /*  Then to determine if two atoms have an exclusion, a linear  */
    /*  search is done on the array of the atom with the smaller    */
    /*  index for the larger index.          */
    /*  If the exclusion policy is set to scaled1-4, there are  */
    /*  actually two lists built.  One contains the pairs of atoms  */
    /*  that are to be exlcuded (i.e., explicit exclusions, 1-2,    */
    /*  and 1-3 interactions) and the other contains just the 1-4   */
    /*  interactions, since they will need to be determined   */
    /*  independantly of the other exclusions.      */
    /*                */
    /****************************************************************/

    void Molecule::build_exclusions()
    {
      register int i;          //  Loop counter
      ExclusionSettings exclude_flag;    //  Exclusion policy

      exclude_flag = simParams->exclude;
      int stripHGroupExclFlag = (simParams->splitPatch == SPLIT_PATCH_HYDROGEN);

      //  Go through the explicit exclusions and add them to the arrays
      for (i=0; i<numExclusions; i++)
      {
        exclusionSet.add(exclusions[i]);
      }

      // If this is AMBER force field, and readExclusions is TRUE,
      // then all the exclusions were read from parm file, and we
      // shouldn't generate any of them.
      if (!simParams->amberOn || !simParams->readExclusions)
      { //  Now calculate the bonded exlcusions based on the exclusion policy
        switch (exclude_flag)
        {
         case NONE:
           break;
         case ONETWO:
           build12excl();
           break;
          case ONETHREE:
            build12excl();
            build13excl();
	    if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
          case ONEFOUR:
            build12excl();
            build13excl();
            build14excl(0);
	    if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
          case SCALED14:
            build12excl();
            build13excl();
            build14excl(1);
	    if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
        }
      }
      else if (stripHGroupExclFlag && exclude_flag!=NONE && exclude_flag!=ONETWO)
        stripHGroupExcl();

      stripFepExcl();
    }
    /*      END OF FUNCTION build_exclusions    */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build12excl        */
    /*                  */
    /************************************************************************/

    void Molecule::build12excl(void)
       
    {
       int32 *current_val;  //  Current value to check
       register int i;    //  Loop counter to loop through all atoms
       
       //  Loop through all the atoms marking the bonded interactions for each one
       for (i=0; i<numAtoms; i++)
       {
      current_val = bondsWithAtom[i];
       
      //  Loop through all the bonds for this atom
      while (*current_val != -1)
      {
         if (bonds[*current_val].atom1 == i)
         {
      if (i<bonds[*current_val].atom2)
      {
         exclusionSet.add(Exclusion(i,bonds[*current_val].atom2));
      }
         }
         else
         {
      if (i<bonds[*current_val].atom1)
      {
         exclusionSet.add(Exclusion(i,bonds[*current_val].atom1));
      }
         }
    
         ++current_val;
      }
       }
    }
    /*      END OF FUNCTION build12excl      */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build13excl        */
    /*                  */
    /************************************************************************/

    void Molecule::build13excl(void)
       
    {
       int32 *bond1, *bond2;  //  The two bonds being checked
       int middle_atom;  //  Common third atom
       register int i;    //  Loop counter to loop through all atoms
       
       //  Loop through all the atoms looking at the bonded connections
       //  for each one
       for (i=0; i<numAtoms; i++)
       {
       bond1 = bondsWithAtom[i];
       
       //  Loop through all the bonds directly connect to atom i
       while (*bond1 != -1)
       {
        if (bonds[*bond1].atom1 == i)
        {
          middle_atom=bonds[*bond1].atom2;
        }
        else
        {
          middle_atom=bonds[*bond1].atom1;
        }

        bond2 = bondsWithAtom[middle_atom];

        //  Now loop through all the bonds connect to the
        //  middle atom
        while (*bond2 != -1)
        {
          if (bonds[*bond2].atom1 == middle_atom)
          {
            if (i < bonds[*bond2].atom2)
            {
              exclusionSet.add(Exclusion(i,bonds[*bond2].atom2));
            }
          }
          else
          {
            if (i < bonds[*bond2].atom1)
            {
              exclusionSet.add(Exclusion(i,bonds[*bond2].atom1));
            }
          }

          ++bond2;
        }

        ++bond1;
      }
       }
    }
    /*      END OF FUNCTION build13excl      */

    /************************************************************************/
    /*                  */
    /*        FUNCTION build14excl      */
    /*                  */
    /************************************************************************/


    void Molecule::build14excl(int modified)
       
    {
       int32 *bond1, *bond2, *bond3;  //  The two bonds being checked
       int mid1, mid2;    //  Middle atoms
       register int i;      //  Counter to loop through all atoms
       
       //  Loop through all the atoms
       for (i=0; i<numAtoms; i++)
       {  
      // Get all the bonds connect directly to atom i
      bond1 = bondsWithAtom[i];
       
      while (*bond1 != -1)
      {
        if (bonds[*bond1].atom1 == i)
        {
          mid1=bonds[*bond1].atom2;
        }
        else
        {
          mid1=bonds[*bond1].atom1;
        }

        bond2 = bondsWithAtom[mid1];

        //  Loop through all the bonds connected to atom mid1
        while (*bond2 != -1)
        {
          if (bonds[*bond2].atom1 == mid1)
          {
            mid2 = bonds[*bond2].atom2;
          }
          else
          {
            mid2 = bonds[*bond2].atom1;
          }

          //  Make sure that we don't double back to where
          //  we started from.  This causes strange behavior.
          //  Trust me, I've been there . . .
          if (mid2 == i)
          {
            ++bond2;
            continue;
          }

          bond3=bondsWithAtom[mid2];

          //  Loop through all the bonds connected to mid2
          while (*bond3 != -1)
          {
            if (bonds[*bond3].atom1 == mid2)
            {
              //  Make sure that we don't double back to where
              //  we started from.  This causes strange behavior.
              //  Trust me, I've been there . . .
              //  I added this!!!  Why wasn't it there before?  -JCP
              if (bonds[*bond3].atom2 != mid1)
              if (i<bonds[*bond3].atom2)
              {
                 exclusionSet.add(Exclusion(i,bonds[*bond3].atom2,modified));
              }
            }
            else
            {
              //  Make sure that we don't double back to where
              //  we started from.  This causes strange behavior.
              //  Trust me, I've been there . . .
              //  I added this!!!  Why wasn't it there before?  -JCP
              if (bonds[*bond3].atom1 != mid1)
              if (i<bonds[*bond3].atom1)
              {
                 exclusionSet.add(Exclusion(i,bonds[*bond3].atom1,modified));
              }
            }

            ++bond3;
          }

          ++bond2;
        }
    
        ++bond1;
      }
       }
    }
    /*      END OF FUNCTION build14excl      */


    /************************************************************************/
    /*                                                                      */
    /*        FUNCTION stripHGroupExcl                                      */
    /*                                                                      */
    /*  This function removes all exclusions which are entirely             */
    /*  within a single hydrogen group.  This assumes that they all         */
    /*  exist, which should be true for exclusion policy 1-3 or higher.     */
    /*                                                                      */
    /************************************************************************/

  void Molecule::stripHGroupExcl(void)
  {

    HydrogenGroup::iterator h_i, h_e, h_j;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      for ( h_j = h_i + 1; h_j != h_e && ! h_j->isGP; ++h_j ) {
	if ( h_i->atomID < h_j->atomID )
	  exclusionSet.del(Exclusion(h_i->atomID,h_j->atomID));
	else
	  exclusionSet.del(Exclusion(h_j->atomID,h_i->atomID));
      }
    }

  }
    /*      END OF FUNCTION stripHGroupExcl      */

    /************************************************************************/
    /*                                                                      */
    /*        FUNCTION stripFepExcl                                         */
    /*                                                                      */
    /************************************************************************/

  void Molecule::stripFepExcl(void)
  {
    UniqueSet<Exclusion> fepExclusionSet;
    UniqueSetIter<Exclusion> exclIter(exclusionSet);

    if ( simParams->fepOn || simParams->lesOn ) {
       for ( exclIter=exclIter.begin(); exclIter != exclIter.end(); exclIter++ )
       {
         int t1 = get_fep_type(exclIter->atom1);
         int t2 = get_fep_type(exclIter->atom2);
         if ( t1 && t2 && t1 != t2 ) {
           fepExclusionSet.add(*exclIter);
         }
       }
    } else if ( simParams->pairInteractionOn ) {
      for ( exclIter=exclIter.begin(); exclIter != exclIter.end(); exclIter++ )
      {
        int ifep_type = get_fep_type(exclIter->atom1);
        int jfep_type = get_fep_type(exclIter->atom2);
        if ( simParams->pairInteractionSelf ) {
          // for pair-self, both atoms must be in group 1.
          if (ifep_type != 1 || jfep_type != 1) {
            fepExclusionSet.add(*exclIter);
          }
        } else {

          // for pair, must have one from each group.
          if (!(ifep_type == 1 && jfep_type == 2) &&
              !(ifep_type == 2 && jfep_type == 1)) {
            fepExclusionSet.add(*exclIter);
          }
        }
       }
    }

    UniqueSetIter<Exclusion> fepIter(fepExclusionSet);
    for ( fepIter=fepIter.begin(); fepIter != fepIter.end(); fepIter++ )
    {
      exclusionSet.del(*fepIter);
    }

  }
    /*      END OF FUNCTION stripFepExcl      */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_constraint_params    */
    /*                  */
    /*   INPUTS:                */
    /*  consref - Value of consref parameter from config file    */
    /*  conskfile - Value of conskfile from config file      */
    /*  conskcol - Value of conskcol from config file      */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*  cwd - Current working directory          */
    /*                  */
    /*  This function builds all the parameters that are necessary  */
    /*   to do harmonic constraints.  This involves looking through    */
    /*   one or more PDB objects to determine which atoms are constrained,  */
    /*   and what the force constant and reference position is force each   */
    /*   atom that is constrained.  This information is then stored    */
    /*   in the arrays consIndexes and consParams.        */
    /*                  */
    /************************************************************************/

    void Molecule::build_constraint_params(StringList *consref, 
             StringList *conskfile, 
             StringList *conskcol, 
             PDB *initial_pdb,
             char *cwd)
       
    {
       PDB *refPDB, *kPDB;    //  Pointer to other PDB's if used
       register int i;      //  Loop counter
       int current_index=0;    //  Index into values used
       int kcol = 4;      //  Column to look for force constant in
       Real kval = 0;      //  Force constant value retreived
       char filename[129];    //  PDB filename
       
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
       consIndexes = new int32[numAtoms];
       
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
    iout << iWARN << "NO CONSTRAINED ATOMS WERE FOUND, BUT CONSTRAINTS ARE ON . . .\n" << endi;
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
    /*      END OF FUNCTION build_constraint_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_movdrag_params  */
/*                  */
/*   INPUTS:        */
/*  movDragFile - value of movDragFile from the config file */
/*  movDragCol - value of movDragCol from the config file */
/*  movDragVelFile - value of movDragVelFile from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do moving drag. This involves looking through one or more    */
/*  PDB objects to determine which atoms are dragged,  and what the */
/*  drag parameters for each atom are. This information is then stored */
/*  in the arrays movDragIndexes and movDragParams. */
/*                  */
/************************************************************************/

void Molecule::build_movdrag_params(StringList *movDragFile, 
				    StringList *movDragCol, 
				    StringList *movDragVelFile, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *vPDB;        //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for drag tag in
  Real dtval = 0;          //  Drag tag value retreived
  char mainfilename[129];  //  main moving drag PDB filename
  char velfilename[129];   //  moving drag velocity PDB filename
  
  //  Get the PDB to read the moving drag tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (movDragFile == NULL) {
    tPDB = initial_pdb;
    
  } else {

    if (movDragFile->next != NULL) {
      NAMD_die("Multiple definitions of moving drag tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (movDragFile->data[0] == '/') ) {
      strcpy(mainfilename, movDragFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, movDragFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_movdrag_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in moving drag tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom velocities. If no name given, use
  // movDragFile if it is defined. Can NOT use the PDB coordinate
  // file!
  
  if (movDragVelFile == NULL) {
    if (movDragFile == NULL) {
      NAMD_die("Moving drag velocity file can not be same as coordinate PDB file");
    } else {
      if (movDragVelFile->next != NULL) {
	NAMD_die("Multiple definitions of moving drag velocity file in configuration file");
      };
      vPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (movDragVelFile->data[0] == '/') ) {
      strcpy(velfilename, movDragVelFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, movDragVelFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_movdrag_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in moving drag velocity PDB doesn't match coordinate PDB");
    }
  };
  
  
  //  Get the column that the drag tag is going to be in. If
  //  movDragFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If movDragFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (movDragCol == NULL) {
    dtcol = 4;
  } else {
    if (movDragCol->next != NULL) {
      NAMD_die("Multiple definitions of drag column in config file");
    };
    
    if (movDragFile == NULL
	&& strcasecmp(movDragCol->data, "B")
	&& strcasecmp(movDragCol->data, "O")) {
      NAMD_die("Can not read moving drag tags from X, Y, or Z column of the coordinate or velocity file");
    };
    if (!strcasecmp(movDragCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(movDragCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(movDragCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(movDragCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(movDragCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("movDragCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the drag
  //  parameters for each atom.  If the atom is not dragged, its
  //  value will be set to -1 in this array.
  movDragIndexes = new int32[numAtoms];
    if (movDragIndexes == NULL) {
    NAMD_die("memory allocation failed in Molecule::build_movdrag_params()");
  };
  
  //  Loop through all the atoms and find out which ones are dragged
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is dragged
      movDragIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not dragged
      movDragIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    //  Drag was turned on, but there weren't really any dragged
    iout << iWARN << "NO DRAGGED ATOMS WERE FOUND, BUT MOVING DRAG IS ON . . . " << endi;
  } else {
    //  Allocate an array to hold the drag parameters
    movDragParams = new MovDragParams[current_index];
    if (movDragParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_movdrag_params");
    }
  };
  
  numMovDrag = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are dragged
  for (i=0; i<numAtoms; i++) {
    if (movDragIndexes[i] != -1) {
      movDragParams[movDragIndexes[i]].v[0] = (vPDB->atom(i))->xcoor();
      movDragParams[movDragIndexes[i]].v[1] = (vPDB->atom(i))->ycoor();
      movDragParams[movDragIndexes[i]].v[2] = (vPDB->atom(i))->zcoor();
    };
  };
      
  if (movDragFile != NULL) delete tPDB;
  if (movDragVelFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_movdrag_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_rotdrag_params  */
/*                  */
/*   INPUTS:        */
/*  rotDragFile - value of rotDragFile from the config file */
/*  rotDragCol - value of rotDragCol from the config file */
/*  rotDragAxisFile - value of rotDragAxisFile from the config file */
/*  rotDragPivotFile - value of rotDragPivotFile from the config file */
/*  rotDragVelFile - value of rotDragVelFile from the config file */
/*  rotDragVelCol - value of rotDragVelCol from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do moving drag. This involves looking through one or more    */
/*  PDB objects to determine which atoms are dragged,  and what the */
/*  drag parameters for each atom are. This information is then stored */
/*  in the arrays rotDragIndexes and rotDragParams. */
/*                  */
/************************************************************************/

void Molecule::build_rotdrag_params(StringList *rotDragFile, 
				    StringList *rotDragCol, 
				    StringList *rotDragAxisFile, 
				    StringList *rotDragPivotFile, 
				    StringList *rotDragVelFile, 
				    StringList *rotDragVelCol, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *aPDB, *pPDB, *vPDB; //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for drag tag in
  Real dtval = 0;          //  Drag tag value retreived
  int dvcol = 4;           //  Column to look for angular velocity in
  Real dvval = 0;          //  Angular velocity value retreived
  char mainfilename[129];  //  main rotating drag PDB filename
  char axisfilename[129];  //  rotating drag axis PDB filename
  char pivotfilename[129]; //  rotating drag pivot point PDB filename
  char velfilename[129];   //  rotating drag angular velocity PDB filename
  
  //  Get the PDB to read the rotating drag tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (rotDragFile == NULL) {
    tPDB = initial_pdb;
    
  } else {

    if (rotDragFile->next != NULL) {
      NAMD_die("Multiple definitions of rotating drag tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (rotDragFile->data[0] == '/') ) {
      strcpy(mainfilename, rotDragFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, rotDragFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom rotation axes. If no name given, use
  // rotDragFile if both it AND rotDragPivotFile are defined. Can NOT
  // use the PDB coordinate file, nor rotDragPivotFile!

  if (rotDragAxisFile == NULL) {
    if (rotDragFile == NULL) {
      NAMD_die("Rotating drag axis file can not be same as coordinate PDB file");
    } else {
      if (rotDragAxisFile->next != NULL) {
	NAMD_die("Multiple definitions of rotating drag axis file in configuration file");
      };
      if (rotDragPivotFile == NULL) {
	NAMD_die("Need to specify at least one of rotDragAxisFile and rotDragPivotFile; they can not be same");
      };
      aPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (rotDragAxisFile->data[0] == '/') ) {
      strcpy(axisfilename, rotDragAxisFile->data);
    } else {
      strcpy(axisfilename, cwd);
      strcat(axisfilename, rotDragAxisFile->data);
    }
    
    aPDB = new PDB(axisfilename);
    if ( aPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (aPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag axis PDB doesn't match coordinate PDB");
    }
  };
  
  // Get the PDB to read atom rotation pivot points. If no name given,
  // use rotDragFile if both it AND rotDragAxisFile are defined. Can
  // NOT use the PDB coordinate file, nor rotDragAxisFile!

  if (rotDragPivotFile == NULL) {
    if (rotDragFile == NULL) {
      NAMD_die("Rotating drag pivot point file can not be same as coordinate PDB file");
    } else {
      if (rotDragPivotFile->next != NULL) {
	NAMD_die("Multiple definitions of rotating drag pivot point file in configuration file");
      };
      if (rotDragAxisFile == NULL) {
	NAMD_die("Need to specify at least one of rotDragAxisFile and rotDragPivotFile; they can not be same");
      };
      pPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (rotDragPivotFile->data[0] == '/') ) {
      strcpy(pivotfilename, rotDragPivotFile->data);
    } else {
      strcpy(pivotfilename, cwd);
      strcat(pivotfilename, rotDragPivotFile->data);
    }
    
    pPDB = new PDB(pivotfilename);
    if ( pPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (pPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag pivot point PDB doesn't match coordinate PDB");
    }
  };
  
  
  // Get the PDB to read atom angular velocities. If no name given,
  // use rotDragFile (or the coordinate PDB file if rotDragFile is not
  // defined).

  if (rotDragVelFile == NULL) {
    vPDB = tPDB;
  } else {
    if (rotDragVelFile->next != NULL) {
      NAMD_die("Multiple definitions of rotating drag velocity file in configuration file");
    };
    
    if ( (cwd == NULL) || (rotDragVelFile->data[0] == '/') ) {
      strcpy(velfilename, rotDragVelFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, rotDragVelFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag velocity PDB doesn't match coordinate PDB");
    }
  };
  
  //  Get the column that the drag tag is going to be in. If
  //  rotDragFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If rotDragFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (rotDragCol == NULL) {
    dtcol = 4;
  } else {
    if (rotDragCol->next != NULL) {
      NAMD_die("Multiple definitions of drag tag column in config file");
    };
    
    if ( rotDragFile == NULL
	 && (!strcasecmp(rotDragCol->data, "X")
	     || !strcasecmp(rotDragCol->data, "Y")
	     || !strcasecmp(rotDragCol->data, "Z"))) {
      NAMD_die("Can not read rotating drag tags from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(rotDragCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(rotDragCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(rotDragCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(rotDragCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(rotDragCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("rotDragCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Get the column that the drag angular velocity is going to be
  //  in. If rotDragVelFile is defined, it can be in any of the 5
  //  floating point fields in the PDB (X, Y, Z, O, or B) which
  //  correspond to the 1st, 2nd, ... 5th floating point fields. If
  //  NEITHER of rotDragVelFile OR rotDragFile is defined, it can
  //  only be O or B fileds. The default is the O (4th) field, which
  //  is the occupancy.

  if (rotDragVelCol == NULL) {
    dvcol = 4;
  } else {
    if (rotDragVelCol->next != NULL) {
      NAMD_die("Multiple definitions of drag angular velocity column in config file");
    };
    
    if (rotDragVelFile == NULL
	&& rotDragFile == NULL
	&& strcasecmp(rotDragCol->data, "B")
	&& strcasecmp(rotDragCol->data, "O")) {
      NAMD_die("Can not read rotating drag angular velocities from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(rotDragVelCol->data, "X")) {
      dvcol=1;
    } else if (!strcasecmp(rotDragVelCol->data, "Y")) {
      dvcol=2;
    } else if (!strcasecmp(rotDragVelCol->data, "Z")) {
      dvcol=3;
    } else if (!strcasecmp(rotDragVelCol->data, "O")) {
      dvcol=4;
    } else if (!strcasecmp(rotDragVelCol->data, "B")) {
      dvcol=5;
    }
    else {
      NAMD_die("rotDragVelCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the drag
  //  parameters for each atom.  If the atom is not dragged, its
  //  value will be set to -1 in this array.
  rotDragIndexes = new int32[numAtoms];
  if (rotDragIndexes == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_rotdrag_params()");
  };
  
  //  Loop through all the atoms and find out which ones are dragged
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is dragged
      rotDragIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not dragged
      rotDragIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    iout << iWARN << "NO DRAGGED ATOMS WERE FOUND, BUT ROTATING DRAG IS ON . . . " << endi;
  } else {
    rotDragParams = new RotDragParams[current_index];
    if (rotDragParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_rotdrag_params");
    }
  };
  
  numRotDrag = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are dragged
  for (i=0; i<numAtoms; i++) {
    if (rotDragIndexes[i] != -1) {
      rotDragParams[rotDragIndexes[i]].a[0] = (aPDB->atom(i))->xcoor();
      rotDragParams[rotDragIndexes[i]].a[1] = (aPDB->atom(i))->ycoor();
      rotDragParams[rotDragIndexes[i]].a[2] = (aPDB->atom(i))->zcoor();
      rotDragParams[rotDragIndexes[i]].p[0] = (pPDB->atom(i))->xcoor();
      rotDragParams[rotDragIndexes[i]].p[1] = (pPDB->atom(i))->ycoor();
      rotDragParams[rotDragIndexes[i]].p[2] = (pPDB->atom(i))->zcoor();
      switch (dvcol) {
      case 1:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->xcoor();
	break;
      case 2:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->ycoor();
	break;
      case 3:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->zcoor();
	break;
      case 4:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->occupancy();
	break;
      case 5:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->temperaturefactor();
	break;
      };
    };
  };
      
  if (rotDragFile != NULL) delete tPDB;
  if (rotDragAxisFile != NULL) delete aPDB;
  if (rotDragPivotFile != NULL) delete pPDB;
  if (rotDragVelFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_rotdrag_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_constorque_params  */
/*                  */
/*   INPUTS:        */
/*  consTorqueFile - value of consTorqueFile from the config file */
/*  consTorqueCol - value of consTorqueCol from the config file */
/*  consTorqueAxisFile - value of consTorqueAxisFile from the config file */
/*  consTorquePivotFile - value of consTorquePivotFile from the config file */
/*  consTorqueValFile - value of consTorqueValFile from the config file */
/*  consTorqueValCol - value of consTorqueValCol from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do "constant" torque. This involves looking through one or more    */
/*  PDB objects to determine which atoms are torqued,  and what the */
/*  torque parameters for each atom are. This information is then stored */
/*  in the arrays consTorqueIndexes and consTorqueParams. */
/*                  */
/************************************************************************/

void Molecule::build_constorque_params(StringList *consTorqueFile, 
				    StringList *consTorqueCol, 
				    StringList *consTorqueAxisFile, 
				    StringList *consTorquePivotFile, 
				    StringList *consTorqueValFile, 
				    StringList *consTorqueValCol, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *aPDB, *pPDB, *vPDB; //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for torque tag in
  Real dtval = 0;          //  Torque tag value retreived
  int dvcol = 4;           //  Column to look for angular velocity in
  Real dvval = 0;          //  Angular velocity value retreived
  char mainfilename[129];  //  main "constant" torque PDB filename
  char axisfilename[129];  //  "constant" torque axis PDB filename
  char pivotfilename[129]; //  "constant" torque pivot point PDB filename
  char velfilename[129];   //  "constant" torque angular velocity PDB filename
  
  //  Get the PDB to read the "constant" torque tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (consTorqueFile == NULL) {
    tPDB = initial_pdb;
    
  } else {

    if (consTorqueFile->next != NULL) {
      NAMD_die("Multiple definitions of \"constant\" torque tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (consTorqueFile->data[0] == '/') ) {
      strcpy(mainfilename, consTorqueFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, consTorqueFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom rotation axes. If no name given, use
  // consTorqueFile if both it AND consTorquePivotFile are defined. Can NOT
  // use the PDB coordinate file, nor consTorquePivotFile!

  if (consTorqueAxisFile == NULL) {
    if (consTorqueFile == NULL) {
      NAMD_die("\"Constant\" torque axis file can not be same as coordinate PDB file");
    } else {
      if (consTorqueAxisFile->next != NULL) {
	NAMD_die("Multiple definitions of \"constant\" torque axis file in configuration file");
      };
      if (consTorquePivotFile == NULL) {
	NAMD_die("Need to specify at least one of consTorqueAxisFile and consTorquePivotFile; they can not be same");
      };
      aPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (consTorqueAxisFile->data[0] == '/') ) {
      strcpy(axisfilename, consTorqueAxisFile->data);
    } else {
      strcpy(axisfilename, cwd);
      strcat(axisfilename, consTorqueAxisFile->data);
    }
    
    aPDB = new PDB(axisfilename);
    if ( aPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (aPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque axis PDB doesn't match coordinate PDB");
    }
  };
  
  // Get the PDB to read atom rotation pivot points. If no name given,
  // use consTorqueFile if both it AND consTorqueAxisFile are defined. Can
  // NOT use the PDB coordinate file, nor consTorqueAxisFile!

  if (consTorquePivotFile == NULL) {
    if (consTorqueFile == NULL) {
      NAMD_die("\"Constant\" torque pivot point file can not be same as coordinate PDB file");
    } else {
      if (consTorquePivotFile->next != NULL) {
	NAMD_die("Multiple definitions of \"constant\" torque pivot point file in configuration file");
      };
      if (consTorqueAxisFile == NULL) {
	NAMD_die("Need to specify at least one of consTorqueAxisFile and consTorquePivotFile; they can not be same");
      };
      pPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (consTorquePivotFile->data[0] == '/') ) {
      strcpy(pivotfilename, consTorquePivotFile->data);
    } else {
      strcpy(pivotfilename, cwd);
      strcat(pivotfilename, consTorquePivotFile->data);
    }
    
    pPDB = new PDB(pivotfilename);
    if ( pPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (pPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque pivot point PDB doesn't match coordinate PDB");
    }
  };
  
  
  // Get the PDB to read atom angular velocities. If no name given,
  // use consTorqueFile (or the coordinate PDB file if consTorqueFile is not
  // defined).

  if (consTorqueValFile == NULL) {
    vPDB = tPDB;
  } else {
    if (consTorqueValFile->next != NULL) {
      NAMD_die("Multiple definitions of \"constant\" torque velocity file in configuration file");
    };
    
    if ( (cwd == NULL) || (consTorqueValFile->data[0] == '/') ) {
      strcpy(velfilename, consTorqueValFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, consTorqueValFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque velocity PDB doesn't match coordinate PDB");
    }
  };
  
  //  Get the column that the torque tag is going to be in. If
  //  consTorqueFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If consTorqueFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (consTorqueCol == NULL) {
    dtcol = 4;
  } else {
    if (consTorqueCol->next != NULL) {
      NAMD_die("Multiple definitions of torque tag column in config file");
    };
    
    if ( consTorqueFile == NULL
	 && (!strcasecmp(consTorqueCol->data, "X")
	     || !strcasecmp(consTorqueCol->data, "Y")
	     || !strcasecmp(consTorqueCol->data, "Z"))) {
      NAMD_die("Can not read \"constant\" torque tags from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(consTorqueCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(consTorqueCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(consTorqueCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(consTorqueCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(consTorqueCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("consTorqueCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Get the column that the torque value is going to be
  //  in. If consTorqueValFile is defined, it can be in any of the 5
  //  floating point fields in the PDB (X, Y, Z, O, or B) which
  //  correspond to the 1st, 2nd, ... 5th floating point fields. If
  //  NEITHER of consTorqueValFile OR consTorqueFile is defined, it can
  //  only be O or B fileds. The default is the O (4th) field, which
  //  is the occupancy.

  if (consTorqueValCol == NULL) {
    dvcol = 4;
  } else {
    if (consTorqueValCol->next != NULL) {
      NAMD_die("Multiple definitions of torque value column in config file");
    };
    
    if (consTorqueValFile == NULL
	&& consTorqueFile == NULL
	&& strcasecmp(consTorqueCol->data, "B")
	&& strcasecmp(consTorqueCol->data, "O")) {
      NAMD_die("Can not read \"constant\" torque values from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(consTorqueValCol->data, "X")) {
      dvcol=1;
    } else if (!strcasecmp(consTorqueValCol->data, "Y")) {
      dvcol=2;
    } else if (!strcasecmp(consTorqueValCol->data, "Z")) {
      dvcol=3;
    } else if (!strcasecmp(consTorqueValCol->data, "O")) {
      dvcol=4;
    } else if (!strcasecmp(consTorqueValCol->data, "B")) {
      dvcol=5;
    }
    else {
      NAMD_die("consTorqueValCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the torque
  //  parameters for each atom.  If the atom is not torqued, its
  //  value will be set to -1 in this array.
  consTorqueIndexes = new int32[numAtoms];
  if (consTorqueIndexes == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_constorque_params()");
  };
  
  //  Loop through all the atoms and find out which ones are torqued
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is torqued
      consTorqueIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not torqued
      consTorqueIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    iout << iWARN << "NO TORQUED ATOMS WERE FOUND, BUT \"CONSTANT\" TORQUE IS ON . . . " << endi;
  } else {
    consTorqueParams = new ConsTorqueParams[current_index];
    if (consTorqueParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_constorque_params");
    }
  };
  
  numConsTorque = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are torqued
  for (i=0; i<numAtoms; i++) {
    if (consTorqueIndexes[i] != -1) {
      consTorqueParams[consTorqueIndexes[i]].a[0] = (aPDB->atom(i))->xcoor();
      consTorqueParams[consTorqueIndexes[i]].a[1] = (aPDB->atom(i))->ycoor();
      consTorqueParams[consTorqueIndexes[i]].a[2] = (aPDB->atom(i))->zcoor();
      consTorqueParams[consTorqueIndexes[i]].p[0] = (pPDB->atom(i))->xcoor();
      consTorqueParams[consTorqueIndexes[i]].p[1] = (pPDB->atom(i))->ycoor();
      consTorqueParams[consTorqueIndexes[i]].p[2] = (pPDB->atom(i))->zcoor();
      switch (dvcol) {
      case 1:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->xcoor();
	break;
      case 2:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->ycoor();
	break;
      case 3:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->zcoor();
	break;
      case 4:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->occupancy();
	break;
      case 5:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->temperaturefactor();
	break;
      };
    };
  };
      
  if (consTorqueFile != NULL) delete tPDB;
  if (consTorqueAxisFile != NULL) delete aPDB;
  if (consTorquePivotFile != NULL) delete pPDB;
  if (consTorqueValFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_constorque_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_constant_forces    */
/*                  */
/*   INPUTS:                */
/*  filename - PDB file containing the constant forces    */
/*                  */
/*  This function reads the constant forces from the PDB file.  */
/*   The force vector to be applied on each atom is determined by:    */
/*     occupancy*(X,Y,Z)   */
/*   Only non-zero forces are stored   */
/*                  */
/************************************************************************/

void Molecule::build_constant_forces(char *filename)
{ int i, index;
  PDB *forcePDB;
  
  if (!filename) {
    // then all forces are zero to begin with; may be changed by
    // the consforceconfig command.
    iout << iWARN << "NO CONSTANT FORCES SPECIFIED, BUT CONSTANT FORCE IS ON . . .\n" << endi;
    consForceIndexes = new int32[numAtoms];
    for (i=0; i<numAtoms; i++) consForceIndexes[i] = -1;
    return;
  }

  if ((forcePDB=new PDB(filename)) == NULL)
    NAMD_die("Memory allocation failed in Molecule::build_constant_forces");
  if (forcePDB->num_atoms() != numAtoms)
    NAMD_die("Number of atoms in constant force PDB doesn't match coordinate PDB");

  //  Allocate an array that will store an index into the constant force
  //  array for each atom.  If the atom has no constant force applied, its
  //  value will be set to -1 in this array.
  consForceIndexes = new int32[numAtoms];
  if (consForceIndexes == NULL)
    NAMD_die("memory allocation failed in Molecule::build_constant_forces()");

  //  Loop through all the atoms and find out which ones have constant force
  numConsForce = 0;
  for (i=0; i<numAtoms; i++)
    if ((forcePDB->atom(i)->xcoor()==0 && forcePDB->atom(i)->ycoor()==0 &&
         forcePDB->atom(i)->zcoor()==0) || forcePDB->atom(i)->occupancy()==0)
      //  This atom has no constant force
      consForceIndexes[i] = -1;
    else
      //  This atom has constant force
      consForceIndexes[i] = numConsForce++;

  if (numConsForce == 0)
    // Constant force was turned on, but there weren't really any non-zero forces
    iout << iWARN << "NO NON-ZERO FORCES WERE FOUND, BUT CONSTANT FORCE IS ON . . .\n" << endi;
  else
  { // Allocate an array to hold the forces
    consForce = new Vector[numConsForce];
    if (consForce == NULL)
      NAMD_die("memory allocation failed in Molecule::build_constant_forces");
    // Loop through all the atoms and assign the forces
    for (i=0; i<numAtoms; i++)
      if ((index=consForceIndexes[i]) != -1)
      { //  This atom has constant force on it
        consForce[index].x = forcePDB->atom(i)->xcoor() * forcePDB->atom(i)->occupancy();
        consForce[index].y = forcePDB->atom(i)->ycoor() * forcePDB->atom(i)->occupancy();
        consForce[index].z = forcePDB->atom(i)->zcoor() * forcePDB->atom(i)->occupancy();
      }
  }

  delete forcePDB;
}
/*      END OF FUNCTION build_constant_forces    */


void Molecule::build_langevin_params(BigReal coupling, Bool doHydrogen) {

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
  BigReal forceConstant = 0.002*TIMEFACTOR*TIMEFACTOR*BOLTZMAN*
			(simParams->langevinTemp)/(simParams->dt);

  //  Loop through all the atoms and get the b value
  for (int i=0; i<numAtoms; i++)
  {
    BigReal bval = coupling;

    if ( (! doHydrogen) && is_hydrogen(i) ) bval = 0;

    //  Assign the b value
    langevinParams[i] = bval;

    //  Calculate the random force value
    langForceVals[i] = sqrt(forceConstant*atoms[i].mass*bval);
  }

}

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_langevin_params      */
    /*                  */
    /*   INPUTS:                */
    /*  langfile - Value of langevinfile from config file    */
    /*  langcol - Value of langevincol from config file      */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*      cwd - Current working directory          */
    /*                  */
    /*  This function builds the array of b values necessary for  */
    /*   Langevin dynamics.  It takes the name of the PDB file and the      */
    /*   column in the PDB file that contains the b values.  It then  */
    /*   builds the array langevinParams for use during the program.  */
    /*                  */
    /************************************************************************/

    void Molecule::build_langevin_params(StringList *langfile, 
           StringList *langcol, 
           PDB *initial_pdb,
           char *cwd)
       
    {
       PDB *bPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       BigReal forceConstant;  //  Constant factor in force calc
       char filename[129];    //  Filename
       
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
    /*      END OF FUNCTION build_langevin_params    */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_fixed_atoms      */
    /*                  */
    /*   INPUTS:              */
    /*  fixedfile - Value of langevinfile from config file    */
    /*  fixedcol - Value of langevincol from config file    */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*      cwd - Current working directory        */
    /*                  */
    /*  This function builds the list of fixed atoms.      */
    /*   It takes the name of the PDB file and the      */
    /*   column in the PDB file that contains the flags.  It then  */
    /*   builds the array fixedAtomFlags for use during the program.  */
    /*                  */
    /************************************************************************/

    void Molecule::build_fixed_atoms(StringList *fixedfile, 
           StringList *fixedcol, 
           PDB *initial_pdb,
           char *cwd)
       
{
       PDB *bPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       char filename[129];    //  Filename
       
       //  Get the PDB object that contains the b values.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates.  
       if (fixedfile == NULL)
       {
    bPDB = initial_pdb;
       }
       else
       {
    if (fixedfile->next != NULL)
    {
       NAMD_die("Multiple definitions of fixed atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (fixedfile->data[0] == '/') )
    {
         strcpy(filename, fixedfile->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, fixedfile->data);
    }
    
    bPDB = new PDB(filename);
    if ( bPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_fixed_atoms");
    }
    
    if (bPDB->num_atoms() != numAtoms)
    {
       NAMD_die("Number of atoms in fixed atoms PDB doesn't match coordinate PDB");
    }
       }
       
       //  Get the column that the b vaules are in.  It
       //  can be in any of the 5 floating point fields in the PDB, according
       //  to what the user wants.  The allowable fields are X, Y, Z, O, or
       //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
       //  The default is the 4th field, which is the occupancy
       if (fixedcol == NULL)
       {
    bcol = 4;
       }
       else
       {
    if (fixedcol->next != NULL)
    {
       NAMD_die("Multiple definitions of fixed atoms column in config file");
    }
    
    if (strcasecmp(fixedcol->data, "X") == 0)
    {
       bcol=1;
    }
    else if (strcasecmp(fixedcol->data, "Y") == 0)
    {
       bcol=2;
    }
    else if (strcasecmp(fixedcol->data, "Z") == 0)
    {
       bcol=3;
    }
    else if (strcasecmp(fixedcol->data, "O") == 0)
    {
       bcol=4;
    }
    else if (strcasecmp(fixedcol->data, "B") == 0)
    {
       bcol=5;
    }
    else
    {
       NAMD_die("fixedatomscol must have value of X, Y, Z, O, or B");
    }
       }
       
       //  Allocate the array to hold all the data
       fixedAtomFlags = new int32[numAtoms];
       
       if (fixedAtomFlags == NULL)
       {
    NAMD_die("memory allocation failed in Molecule::build_fixed_atoms()");
       }
       
  numFixedAtoms = 0;

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
    if ( bval != 0 ) {
      fixedAtomFlags[i] = 1;
      numFixedAtoms++;
    }
    else {
      fixedAtomFlags[i] = 0;
    }
       }
       
       //  If we had to create a PDB object, delete it now
       if (fixedfile != NULL)
       {
    delete bPDB;
       }

  // now figure out how we interact with rigidBonds 
  // this is mainly for degree of freedom counting
  if ( numRigidBonds ) {
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    int parentIsFixed = 0;
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) {
	parentIsFixed = fixedAtomFlags[h_i->atomID];
	if ( (rigidBondLengths[h_i->atomID] > 0.)  // water
		&& fixedAtomFlags[h_i[1].atomID]
		&& fixedAtomFlags[h_i[2].atomID] ) {
	  ++numFixedRigidBonds;
	}
      } else {
	if ( (rigidBondLengths[h_i->atomID] > 0.)
		&& fixedAtomFlags[h_i->atomID]
		&& parentIsFixed ) {
	  ++numFixedRigidBonds;
	}
      }
    }
  }

  // how many hydrogen groups are completely fixed
  {
    numFixedGroups = 0;
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    int allFixed = 0;
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) {
        if ( allFixed ) ++numFixedGroups;
        allFixed = 1;
      }
      allFixed = allFixed && fixedAtomFlags[h_i->atomID];
    }
    if ( allFixed ) ++numFixedGroups;
  }

}
    /*      END OF FUNCTION build_fixed_atoms    */




/************************************************************************/
/*                                                                      */
/*      FUNCTION build_stirred_atoms                                    */
/*                                                                      */
/*   INPUTS:                                                            */
/*  stirredfile - Value of stirFilename from config file    */
/*  stirredcol - Value of stircol from config file (but B, O only */
/*                    since we need X, Y, Z!     */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory        */
/*                  */
/*  This function builds the list of fixed atoms.      */
/*   It takes the name of the PDB file and the      */
/*   column in the PDB file that contains the flags.  It then  */
/*   builds the array fixedAtomFlags for use during the program.  */
/*                                                                      */
/************************************************************************/

    void Molecule::build_stirred_atoms(StringList *stirredfile, 
           StringList *stirredcol, 
           PDB *initial_pdb,
           char *cwd)
       
{
       PDB *sPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       int current_index=0;    //  Index into values used
       char filename[129];    //  Filename
       
       //  Get the PDB object that contains the b values.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates. 
       //  use posssible only if this is 'optional' in simulation parameters
       // dangerous, since restarted simulations will be different
       if (stirredfile == NULL)
       {
    sPDB = initial_pdb;
    // dangerous, since restarted simulations will be different, so warn
        iout << iWARN << "STIRRING USING INITIAL POSITION FILE FOR REFERENCE POSITIONS" << endi;
       }
       else
       {
    if (stirredfile->next != NULL)
    {
       NAMD_die("Multiple definitions of stirred atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (stirredfile->data[0] == '/') )
    {
         strcpy(filename, stirredfile->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, stirredfile->data);
    }
        
    //CkPrintf ("DEBUG: the stir filename is %s\n",filename);    
    sPDB = new PDB(filename);

    if ( sPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_stirred_atoms");

    }

    if (sPDB->num_atoms() != numAtoms)
      {
	NAMD_die("Number of atoms in stirred atoms PDB doesn't match coordinate PDB");
      }
 
       }

//  Get the column that the b vaules are in.  It
//  can be in any of the 5 floating point fields in the PDB, according
//  to what the user wants.  The allowable fields are X, Y, Z, O, or
//  B which correspond to the 1st, 2nd, ... 5th floating point fields.
//  The default is the 4th field, which is the occupancy

 
      if (stirredcol == NULL)
      {
	 bcol = 4;
      }
      else
    {
      if (stirredcol->next != NULL)
	{
	  NAMD_die("Multiple definitions of stirred atoms column in config file");
	}
      
      if  (strcasecmp(stirredcol->data, "O") == 0)
	{
	  bcol=4;
	}
      else if (strcasecmp(stirredcol->data, "B") == 0)
	{
	  bcol=5;
	}
      else
	{
	  NAMD_die("stirredAtomsCol must have value of O or B");
	}
    }
          
       //  Allocate an array that will store an index into the stir
       //  parameters for each atom.  If the atom is not stirred, its
       //  value will be set to -1 in this array.
       stirIndexes = new int32[numAtoms];
       
       if (stirIndexes == NULL)
       {
    NAMD_die("memory allocation failed in Molecule::build_stirred_params()");
       }
       
       current_index = 0;
       //  Loop through all the atoms and find out which ones are stirred
       for (i=0; i<numAtoms; i++)
       {



	 //  Get the b value based on where we were told to find it
	 switch (bcol)
	   {

	   case 4:
	     bval = (sPDB->atom(i))->occupancy();
	     break;
	   case 5:
	     bval = (sPDB->atom(i))->temperaturefactor();
	     break;
	   }

	// CkPrintf ("DEBUG--for atom i= %d  bcol= %d    bval= %g   occ= %g   realbval= %g  x= %g numAtoms= %d test= %g\n" ,i        ,bcol, bval, (sPDB->atom(i))->occupancy(), (sPDB->atom(i))->temperaturefactor(),(sPDB->atom(i))->xcoor(), sPDB->num_atoms(), 0.123 );
	 //  Assign the b value
	 if ( bval != 0 )
	 {
	   // This atom is stirred 
	   stirIndexes[i] = current_index;
	   current_index++;
	 }
	 else
	 {
           //This atom is not stirred 
	   stirIndexes[i] = -1;
	 }
       }
	 

    

       
       if (current_index == 0)
       {
    //  Stirring was turned on, but there weren't really any stirred atoms found in file
    iout << iWARN << "NO STIRRED ATOMS WERE FOUND, BUT STIRRING TORQUES ARE ON . . .\n" << endi;
       }
       else
       {
    //  Allocate an array to hold the stirring parameters
    stirParams = new StirParams[current_index];
    
    if (stirParams == NULL)
    {
       NAMD_die("memory allocation failed in Molecule::build_stir_params");
    }
       }
       
       numStirredAtoms = current_index;
       
       //  Loop through all the atoms and assign the parameters for those
       //  that are stirred
       for (i=0; i<numAtoms; i++)
       {
    if (stirIndexes[i] != -1)
    {
       
       //  This atom is stirred, so get the reference position
       stirParams[stirIndexes[i]].refPos.x = (sPDB->atom(i))->xcoor();
       stirParams[stirIndexes[i]].refPos.y = (sPDB->atom(i))->ycoor();
       stirParams[stirIndexes[i]].refPos.z = (sPDB->atom(i))->zcoor();
    }
       }
       
       //  If we had to create a PDB object, delete it now
       if (stirredfile != NULL)
       {
	 delete sPDB;
       }
       

    }

    /*      END OF FUNCTION build_stirred_atoms    */

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
   //
   //
   //  FUNCTION build_fep_flags
   //
   // INPUTS:
   // fepfile - Value of fepfile read from config file
   // fepcol - Value of fep column, read from config file
   // initial_pdb - PDB object that contains the initial positions
   //  cwd - current working directory
   //
   // This function builds the array of state values necessary
   // for FEP. It takes the name of the PDB file and column in
   // the PDB file that contains the fep flag. It then builds
   // the array FepParams for use in the program.
   //
   //
   
   void Molecule::build_fep_flags(StringList *fepfile,
         StringList *fepcol,
         PDB *initial_pdb,
         char *cwd)
   {
     PDB *bPDB;  //Pointer to PDB object to use
     int bcol = 5;  //Column that the data is in
     Real bval = 0; //flag from PDB file
     int i;         // loop counter
     char filename[129]; // filename

     // get the pdb object that contains the fep flags.
     // if the user gave another filename, use it, else
     // use the pdb file with the initial coordinates
     if (fepfile == NULL) {
       bPDB = initial_pdb;
     }
     else {
       if (fepfile->next != NULL) {
        NAMD_die("Multiple definitions of fep PDB file in configuration file");
       }
   
       if ((cwd == NULL) || (fepfile->data[0] == '/')) {
         strcpy(filename, fepfile->data);
       }
       else {
         strcpy(filename, cwd);
         strcat(filename, fepfile->data);
       }

       bPDB = new PDB(filename);
       if (bPDB == NULL) {
         NAMD_die("Memory allocation failed in Molecule:build_fep_flags");
       }

       if (bPDB->num_atoms() != numAtoms) {
         NAMD_die("Number of atoms in fep PDB doesnt match coordinate PDB");
       }
    }
   
    // Get the column that the fep flag is in. It can be in any of the 5 
    // floating point fields in the PDB ie X, Y, X, O or B.
    // The default is 4th field ie the occupancy field
    if (fepcol == NULL) {
      bcol = 5;
    }
    else {
      if (fepcol->next != NULL) {
        NAMD_die("Multiple definitions of fep parameter column in config file");
      }

      if (strcasecmp(fepcol->data, "X") == 0) {
       bcol = 1;
      }
      else if (strcasecmp(fepcol->data, "Y") == 0) {
       bcol = 2;
      }
      else if (strcasecmp(fepcol->data, "Z") == 0) {
       bcol = 3;
      }
      else if (strcasecmp(fepcol->data, "O") == 0) {
       bcol = 4;
      }
      else if (strcasecmp(fepcol->data, "B") == 0) {
       bcol = 5;
      }
      else {
       NAMD_die("fepcol must have value of X, Y, Z, O or B");
      }
    }

   cout << "To read fep data from file: " << filename << endl;
   cout << "To read fep flag data from column: " << bcol << endl;
 
   //  Allocate the array to hold all the fep data
   fepAtomFlags = new unsigned char[numAtoms];
       
   if (fepAtomFlags == NULL) {
    NAMD_die("Memory allocation failed in Molecule::build_fep_params()");
   }

   // loop through all the atoms and get the b value
   for (i = 0; i < numAtoms; i++) {
   // Get the fep flag value
    switch (bcol) {
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

    // Assign fep flag value
    if (simParams->lesOn) {
      if ( bval == (int) bval && bval > 0 ) {
        if ( bval > simParams->lesFactor ) 
          NAMD_die("LES flag must be less than or equal to lesFactor.");
        fepAtomFlags[i] = (int) bval;
        numFepFinal++;
        numFepInitial++;
      } else {
        fepAtomFlags[i] = 0;
      }
    } else if (simParams->fepOn) {
      if (bval == 1.0) {
        fepAtomFlags[i] = 1;
        numFepFinal++;
      } else if (bval == -1.0) {
        fepAtomFlags[i] = 2;
        numFepInitial++;
      } else {
        fepAtomFlags[i] = 0;
      }
    } else if (simParams->pairInteractionOn) {
      if (bval == simParams->pairInteractionGroup1) {
        fepAtomFlags[i] = 1;
        numFepInitial++;
      } else if (bval == simParams->pairInteractionGroup2) {
        fepAtomFlags[i] = 2;
        numFepFinal++;
      } else {
        fepAtomFlags[i] = 0;
      }
    }
  }

  // if PDB object was created, delete it
  if (fepfile != NULL) {
    delete bPDB;
  }

}
 // End of function build_fep_flags

//fepe

void Molecule::build_exPressure_atoms(StringList *fixedfile, 
   StringList *fixedcol, PDB *initial_pdb, char *cwd) {
       
  PDB *bPDB;      //  Pointer to PDB object to use
  int bcol = 4;      //  Column that data is in
  Real bval = 0;      //  b value from PDB file
  int i;      //  Loop counter
  char filename[129];    //  Filename

  //  Get the PDB object that contains the b values.  If
  //  the user gave another file name, use it.  Otherwise, just use
  //  the PDB file that has the initial coordinates.
  if (fixedfile == NULL) {
    bPDB = initial_pdb;
  } else {
    if (fixedfile->next != NULL) {
      NAMD_die("Multiple definitions of excluded pressure atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (fixedfile->data[0] == '/') ) {
         strcpy(filename, fixedfile->data);
    } else {
         strcpy(filename, cwd);
         strcat(filename, fixedfile->data);
    }
    bPDB = new PDB(filename);
    if ( bPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_exPressure_atoms");
    }

    if (bPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in excludedPressure atoms PDB doesn't match coordinate PDB");
    }
  }

  //  Get the column that the b vaules are in.  It
  //  can be in any of the 5 floating point fields in the PDB, according
  //  to what the user wants.  The allowable fields are X, Y, Z, O, or
  //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
  //  The default is the 4th field, which is the occupancy
  if (fixedcol == NULL) {
    bcol = 4;
  } else {
    if (fixedcol->next != NULL) {
      NAMD_die("Multiple definitions of excludedPressure atoms column in config file");
    }

    if (strcasecmp(fixedcol->data, "X") == 0) {
       bcol=1;
    } else if (strcasecmp(fixedcol->data, "Y") == 0) {
       bcol=2;
    } else if (strcasecmp(fixedcol->data, "Z") == 0) {
       bcol=3;
    } else if (strcasecmp(fixedcol->data, "O") == 0) {
       bcol=4;
    } else if (strcasecmp(fixedcol->data, "B") == 0) {
       bcol=5;
    } else {
       NAMD_die("excludedPressureFileCol must have value of X, Y, Z, O, or B");
    }
  }

  //  Allocate the array to hold all the data
  exPressureAtomFlags = new int32[numAtoms];

  if (exPressureAtomFlags == NULL) {
    NAMD_die("memory allocation failed in Molecule::build_fixed_atoms()");
  }

  numExPressureAtoms = 0;

  //  Loop through all the atoms and get the b value
  for (i=0; i<numAtoms; i++) {
    //  Get the k value based on where we were told to find it
    switch (bcol) {
       case 1: bval = (bPDB->atom(i))->xcoor(); break;
       case 2: bval = (bPDB->atom(i))->ycoor(); break;
       case 3: bval = (bPDB->atom(i))->zcoor(); break;
       case 4: bval = (bPDB->atom(i))->occupancy(); break;
       case 5: bval = (bPDB->atom(i))->temperaturefactor(); break;
    }

    //  Assign the b value
    if ( bval != 0 ) {
      exPressureAtomFlags[i] = 1;
      numExPressureAtoms++;
    } else {
      exPressureAtomFlags[i] = 0;
    }
  }
  if (fixedfile != NULL) 
    delete bPDB;

  iout << iINFO << "Got " << numExPressureAtoms << " excluded pressure atoms." 
       << endi;
}


    Bool Molecule::is_hydrogen(int anum)
    {
  return ((atoms[anum].status & HydrogenAtom) != 0);
    }

    Bool Molecule::is_oxygen(int anum)
    {
  return ((atoms[anum].status & OxygenAtom) != 0);
    }

    Bool Molecule::is_hydrogenGroupParent(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].isGP);
    }

    Bool Molecule::is_water(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].sortVal == 2);
    }

    int Molecule::get_groupSize(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].atomsInGroup);
    }

    int Molecule::get_mother_atom(int anum)
    {
  // for efficiency reasons, we are not checking if anum is already 
  // hydrogen or not. This function must be called for hydrogens only;
  return atoms[anum].partner;
    }


// go through the molecular structure, analyze the status of each atom,
// and save the data in the Atom structures stored for each atom.  This
// could be built up incrementally while the molecule is being read in,
// but doing it all in one shot allows us to just send the basic info
// over the network and have each node calculate the rest of the data on
// it's own.
void Molecule::build_atom_status(void) {
  register int i;
  int a1, a2, a3;

  // initialize information for each atom (note that the status has
  // already been initialized during the read/receive phase)
  hydrogenGroup.resize(numAtoms);
  HydrogenGroupID *hg = hydrogenGroup.begin();
  for (i=0; i < numAtoms; i++) {
    atoms[i].partner = (-1);
    hg[i].atomID = i;  // currently unsorted
    hg[i].atomsInGroup = 1;  // currently only 1 in group
    hg[i].isGP = 1;  // assume it is a group parent
    hg[i].GPID = i;  // assume it is a group parent
    hg[i].sortVal = 0;  // for group sorting
  }

  // deal with H-H bonds in a sane manner
  // this information will be rewritten later if bonded elsewhere
  int hhbondcount = 0;
  for (i=0; i < numBonds; i++) {
    a1 = bonds[i].atom1;
    a2 = bonds[i].atom2;
    if (is_hydrogen(a1) && is_hydrogen(a2)) {
      ++hhbondcount;
      // make H atoms point at each other for now
      atoms[a1].partner = a2;
      atoms[a2].partner = a1;
      hg[a1].atomsInGroup++;
      hg[a1].GPID = a2;
      hg[a2].atomsInGroup++;
      hg[a2].GPID = a1;
    }
  }
  if ( hhbondcount && ! CmiMyPe() ) {
    iout << iWARN << "Found " << hhbondcount << " H-H bonds.\n" << endi;
  }

  // find which atom each hydrogen is bound to
  // also determine number of atoms in each group
  for (i=0; i < numBonds; i++) {
    a1 = bonds[i].atom1;
    a2 = bonds[i].atom2;
    if (is_hydrogen(a1)) {
      if (is_hydrogen(a2)) continue;
      atoms[a1].partner = a2;
      hg[a2].atomsInGroup++;
      hg[a1].atomsInGroup = 0;
      hg[a1].GPID = a2;
      hg[a1].isGP = 0;
      // check for waters (put them in their own groups: OH or OHH)
      if (is_oxygen(a2))  hg[a2].sortVal++;
    }
    if (is_hydrogen(a2)) {
      atoms[a2].partner = a1;
      hg[a1].atomsInGroup++;
      hg[a2].atomsInGroup = 0;
      hg[a2].GPID = a1;
      hg[a2].isGP = 0;
      // check for waters (put them in their own groups: OH or OHH)
      if (is_oxygen(a1))  hg[a1].sortVal++;
    }
  }

  // check up on our H-H bonds and general sanity check
  int hGPcount = 0;
  for(i=0; i<numAtoms; i++) {
    if ( ! hg[hg[i].GPID].isGP ) {
      char msg[256];
      sprintf(msg, "child atom %d bonded only to child H atoms",i+1);
      NAMD_die(msg);
    }
    if ( hg[i].isGP && is_hydrogen(i) ) {
      if ( hg[i].GPID == i ) continue;  // atomic hydrogen ion
      ++hGPcount;  // molecular hydrogen
      if ( is_hydrogen(hg[i].GPID) && hg[hg[i].GPID].GPID != i ) {
        char msg[256];
        sprintf(msg, "H atom %d bonded only to child H atoms",i+1);
        NAMD_die(msg);
      }
      hg[i].GPID = i;
      hg[hg[i].GPID].atomsInGroup = 0;
      hg[hg[i].GPID].GPID = 0;
      hg[hg[i].GPID].isGP = 0;
      if ( hg[i].atomsInGroup != 2 ) {
        char msg[256];
        sprintf(msg, "H atom %d bonded to multiple H atoms",i+1);
        NAMD_die(msg);
      }
    }
  }
  if ( hGPcount && ! CmiMyPe() ) {
    iout << iWARN << "Found " << hGPcount << " H-H molecules.\n" << endi;
  }

  // sort the hydrogenGroup list and count number of groups
  numHydrogenGroups = 0;
  for(i=0; i<numAtoms; i++)
  {
    // make H follow their group parents.
    if (!hg[i].isGP)  hg[i].sortVal = hg[hg[i].GPID].sortVal;
    else ++numHydrogenGroups;
  }
  hydrogenGroup.sort();

  // finally, add the indexing from atoms[] to hydrogenGroup[]
  waterIndex = numAtoms;
  for(i=0; i<numAtoms; i++)
    {
    atoms[hydrogenGroup[i].atomID].hydrogenList = i;
    // identify where waters start
    if ((hydrogenGroup[i].sortVal == 2) && (i < numAtoms))
  waterIndex = i;
    }

  #if 0
  // debugging code for showing sorted atoms
  for(i=0; i<numAtoms; i++)
    iout << i << " atomID=" << hydrogenGroup[i].atomID
   << " isGP=" << hydrogenGroup[i].isGP
   << " parent=" << hydrogenGroup[i].GPID
   << " #" << hydrogenGroup[i].atomsInGroup
   << " sortVal=" << hydrogenGroup[i].sortVal
   << "\n" << endi;
  #endif

  // now deal with rigidBonds
  if ( simParams->rigidBonds != RIGID_NONE || simParams->mollyOn ) {

    delete [] rigidBondLengths;
    rigidBondLengths = new Real[numAtoms];
    if ( ! rigidBondLengths ) {
      NAMD_die("Memory allocation failed in Molecule::build_atom_status()\n");
    }
    for (i=0; i<numAtoms; ++i) rigidBondLengths[i] = 0;
    int mode = simParams->rigidBonds;
    if ( simParams->mollyOn ) mode = RIGID_ALL;

    // add H-mother lengths or 0 if not constrained
    for (i=0; i < numBonds; i++) {
      a1 = bonds[i].atom1;
      a2 = bonds[i].atom2;
      Real dum, x0;
      params->get_bond_params(&dum,&x0,bonds[i].bond_type);
      if (is_hydrogen(a2)) { int tmp = a1;  a1 = a2;  a2 = tmp; } // swap
      if (is_hydrogen(a1)) {
        if ( is_hydrogen(a2) ) {  // H-H
          if ( ! is_water(a2) ) {  // H-H but not water
	    rigidBondLengths[a1] = ( mode == RIGID_ALL ? x0 : 0. );
	    rigidBondLengths[a2] = ( mode == RIGID_ALL ? x0 : 0. );
          }
        } else if ( is_water(a2) || mode == RIGID_ALL ) {
	  rigidBondLengths[a1] = x0;
	} else {
	  rigidBondLengths[a1] = 0.;
        }
      }
    }

    // zero out H-H lengths - water handled below
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) rigidBondLengths[h_i->atomID] = 0.;
    }

    // fill in H-H lengths for water by searching angles - yuck
    for (i=0; i < numAngles; i++) {
      a2 = angles[i].atom2;
      if ( ! is_water(a2) ) continue;
      a1 = angles[i].atom1;
      a3 = angles[i].atom3;
      if ( rigidBondLengths[a1] != rigidBondLengths[a3] ) {
	NAMD_die("Asymmetric water molecule found???  This can't be right.\n");
      }
      Real dum, t0;
      params->get_angle_params(&dum,&t0,&dum,&dum,angles[i].angle_type);
      rigidBondLengths[a2] = 2. * rigidBondLengths[a1] * sin(0.5*t0);
    }

    // fill in H-H lengths for waters that are missing angles
    int numBondWaters = 0;
    int numFailedWaters = 0;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP && is_water(h_i->atomID) &&
                     rigidBondLengths[h_i->atomID] == 0. ) {
        if ( h_i + 1 == h_e || h_i + 2 == h_e ||
             h_i[1].isGP || h_i[2].isGP || h_i->atomsInGroup != 3 ) {
          NAMD_die("Abnormal water detected.");
        }
        Bond btmp;
        btmp.atom1 = h_i[1].atomID;
        char atom1name[11];
	strcpy(atom1name,get_atomtype(btmp.atom1));
        btmp.atom2 = h_i[2].atomID;
        char atom2name[11];
	strcpy(atom2name,get_atomtype(btmp.atom2));
        params->assign_bond_index(atom1name,atom2name,&btmp);
        Real k, x0;
	x0 = 0.;
        params->get_bond_params(&k,&x0,btmp.bond_type);
	if ( x0 > 0. ) {
          rigidBondLengths[h_i->atomID] = x0;
	  numBondWaters++;
        } else {
	  numFailedWaters++;
        }
      }
    }
    if ( numBondWaters + numFailedWaters ) {
      iout << iWARN << "Missing angles for " <<
	      ( numBondWaters + numFailedWaters ) << " waters.\n" << endi;
    }
    if ( numBondWaters ) {
      iout << iWARN << "Obtained H-H distance from bond parameters for " <<
	      numBondWaters << " waters.\n" << endi;
    }
    if ( numFailedWaters ) {
      iout << iERROR << "Failed to obtain H-H distance from angles or bonds for " <<
	      numFailedWaters << " waters.\n" << endi;
    }

    // in case both molly and rigidBonds are in use make lengths which
    // are molly-only negative and leave lengths which are both alone
    if ( simParams->mollyOn ) {
      mode = simParams->rigidBonds;
      if ( mode == RIGID_NONE ) {
	for (i=0; i<numAtoms; ++i) rigidBondLengths[i] *= -1;
      } else if ( mode == RIGID_WATER ) {
	for (i=0; i<numAtoms; ++i) {
	  if ( ! is_water(i) ) rigidBondLengths[i] *= -1;
	}
      }
    }

    numRigidBonds = 0;
    for (i=0; i<numAtoms; ++i) {
      if ( rigidBondLengths[i] > 0. ) ++numRigidBonds;
    }

  }

}


int Molecule::checkexcl(int atom1, int atom2) const {

  int amin = all_exclusions[atom1].min;
  int amax = all_exclusions[atom1].max;
  if ( atom2 < amin || atom2 > amax ) return 0;
  else return all_exclusions[atom1].flags[atom2-amin];

}


/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for reading AMBER topology data    */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param, Ambertoppar *amber_data)
{
  initialize(simParams,param);

  read_parm(amber_data);
}
/*      END OF FUNCTION Molecule      */


/************************************************************************/
/*                  */
/*      FUNCTION read_parm   */
/*                  */
/*   INPUTS:                */
/*  amber_data - AMBER data structure    */
/*                  */
/*  This function copys AMBER topology data to the corresponding data  */
/*   structures      */
/*                  */
/************************************************************************/

void Molecule::read_parm(Ambertoppar *amber_data)
{
  int i,j,ntheth,nphih,current_index,a1,a2,
      max,min,index,found;

  if (!amber_data->data_read)
    NAMD_die("No data read from parm file yet!");

  // Copy atom informations
  numAtoms = amber_data->Natom;
  atoms = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];
  if (atoms == NULL || atomNames == NULL )
    NAMD_die("memory allocation failed when reading atom information");
  ResidueLookupElem *tmpResLookup = resLookup;
  for (i=0; i<numAtoms; ++i)
  { atomNames[i].resname = nameArena->getNewArray(5);
    atomNames[i].atomname = nameArena->getNewArray(5);
    atomNames[i].atomtype = nameArena->getNewArray(5);
    if (atomNames[i].resname == NULL || atomNames[i].atomname == NULL || atomNames[i].atomtype == NULL)
      NAMD_die("memory allocation failed when reading atom information");
    for (j=0; j<4; ++j)
    { atomNames[i].resname[j] = amber_data->ResNames[amber_data->AtomRes[i]*4+j];
      atomNames[i].atomname[j] = amber_data->AtomNames[i*4+j];
      atomNames[i].atomtype[j] = amber_data->AtomSym[i*4+j];
    }
    atomNames[i].resname[4] = atomNames[i].atomname[4] = atomNames[i].atomtype[4] = '\0';
    strtok(atomNames[i].resname," ");
    strtok(atomNames[i].atomname," ");
    strtok(atomNames[i].atomtype," ");
    atoms[i].mass = amber_data->Masses[i];
    // Divide by 18.2223 to convert to charge in units of the electron charge
    atoms[i].charge = amber_data->Charges[i] / 18.2223;
    atoms[i].vdw_type = amber_data->Iac[i] - 1;

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append("MAIN", amber_data->AtomRes[i]+1, i);

    /*  Determine the type of the atom (H or O) */
    atoms[i].status = UnknownAtom; // the default
    if (atoms[i].mass <=3.5) {
      atoms[i].status |= HydrogenAtom;
    } else if ((atomNames[i].atomname[0] == 'O') && 
         (atoms[i].mass >= 14.0) && 
         (atoms[i].mass <= 18.0)) {
      atoms[i].status |= OxygenAtom;
    }
  }

// Note: In AMBER, the atom numbers in bond, angle and dihedral arrays are in fact
// (3*(atnum-1)). So we divide them by 3 to get the real indices of atom array. Also
// note that NAMD indexes arrays from 0 to NumAtoms-1.

  // Copy bond information
  // Fake bonds (bonds with 0 force constant) are ignored
  Real k, x0;
  numBonds = 0;
  if (amber_data->Nbonh + amber_data->Nbona > 0)
  { bonds = new Bond[amber_data->Nbonh + amber_data->Nbona];
    if (bonds == NULL || amber_data->Nbona < 0)
      NAMD_die("memory allocation failed when reading bond information");
    // Bonds WITH hydrogen
    for (i=0; i<amber_data->Nbonh; ++i)
    { bonds[numBonds].atom1 = amber_data->BondHAt1[i] / 3;
      bonds[numBonds].atom2 = amber_data->BondHAt2[i] / 3;
      bonds[numBonds].bond_type = amber_data->BondHNum[i] - 1;
      if (bonds[numBonds].atom1>=numAtoms || bonds[numBonds].atom2>=numAtoms ||
          bonds[numBonds].bond_type>=amber_data->Numbnd)
      { char err_msg[128];
        sprintf(err_msg, "BOND (WITH H) # %d OVERFLOW IN PARM FILE", i+1);
        NAMD_die(err_msg);
      }
      params->get_bond_params(&k,&x0,bonds[numBonds].bond_type);
      if ( k != 0. ) ++numBonds;  // real bond
    }
    // Bonds WITHOUT hydrogen
    for (i=amber_data->Nbonh; i<amber_data->Nbonh+amber_data->Nbona; ++i)
    { bonds[numBonds].atom1 = amber_data->BondAt1[i-amber_data->Nbonh] / 3;
      bonds[numBonds].atom2 = amber_data->BondAt2[i-amber_data->Nbonh] / 3;
      bonds[numBonds].bond_type = amber_data->BondNum[i-amber_data->Nbonh] - 1;
      if (bonds[i].atom1>=numAtoms || bonds[i].atom2>=numAtoms ||
          bonds[i].bond_type>=amber_data->Numbnd)
      { char err_msg[128];
        sprintf(err_msg, "BOND (WITHOUT H) # %d OVERFLOW IN PARM FILE", i+1-amber_data->Nbonh);
        NAMD_die(err_msg);
      }
      params->get_bond_params(&k,&x0,bonds[numBonds].bond_type);
      if ( k != 0. ) ++numBonds;  // real bond
    }
  }
  /*  Tell user about our subterfuge  */
  if ( numBonds !=  amber_data->Nbonh + amber_data->Nbona) {
    iout << iWARN << "Ignored " << amber_data->Nbonh + amber_data->Nbona - numBonds <<
            " bonds with zero force constants.\n" << endi;
    iout << iWARN <<
	"Will get H-H distance in rigid H2O from H-O-H angle.\n" << endi;
  }

  // Copy angle information
  numAngles = amber_data->Ntheth + amber_data->Ntheta;
  if (numAngles > 0)
  { ntheth = amber_data->Ntheth;
    angles = new Angle[numAngles];
    if (angles == NULL || numAngles < ntheth)
      NAMD_die("memory allocation failed when reading angle information");
    // Angles WITH hydrogen
    for (i=0; i<ntheth; ++i)
    { angles[i].atom1 = amber_data->AngleHAt1[i] / 3;
      angles[i].atom2 = amber_data->AngleHAt2[i] / 3;
      angles[i].atom3 = amber_data->AngleHAt3[i] / 3;
      angles[i].angle_type = amber_data->AngleHNum[i] - 1;
      if (angles[i].atom1>=numAtoms || angles[i].atom2>=numAtoms ||
          angles[i].atom3>=numAtoms || angles[i].angle_type>=amber_data->Numang)
      { char err_msg[128];
        sprintf(err_msg, "ANGLE (WITH H) # %d OVERFLOW IN PARM FILE", i+1);
        NAMD_die(err_msg);
      }
    }
    // Angles WITHOUT hydrogen
    for (i=ntheth; i<numAngles; ++i)
    { angles[i].atom1 = amber_data->AngleAt1[i-ntheth] / 3;
      angles[i].atom2 = amber_data->AngleAt2[i-ntheth] / 3;
      angles[i].atom3 = amber_data->AngleAt3[i-ntheth] / 3;
      angles[i].angle_type = amber_data->AngleNum[i-ntheth] - 1;
      if (angles[i].atom1>=numAtoms || angles[i].atom2>=numAtoms ||
          angles[i].atom3>=numAtoms || angles[i].angle_type>=amber_data->Numang)
      { char err_msg[128];
        sprintf(err_msg, "ANGLE (WITHOUT H) # %d OVERFLOW IN PARM FILE", i+1-ntheth);
        NAMD_die(err_msg);
      }
    }
  }

  numExclusions =  0;
  // If readExclusions is TRUE, then we copy exclusions from parm
  // file; otherwise we skip the exclusions here and generate
  // them later in build_exclusions()
  if (simParams->readExclusions)
  { // Copy exclusion information
    // In Amber data structure, Iblo[] is the number of exclusions
    // for each atom; ExclAt[] is the atom index for the excluded atoms.
    exclusions = new Exclusion[amber_data->Nnb];
    if (exclusions == NULL &&  amber_data->Nnb > 0)
      NAMD_die("memory allocation failed when reading exclusion information");
    current_index = 0;
    for (i=0; i<numAtoms; ++i)
      for (j=0; j<amber_data->Iblo[i]; ++j)
      { if (current_index >= amber_data->Nnb)
        { char err_msg[128];
          sprintf(err_msg, "EXCLUSION INDEX EXCEEDS NUMBER OF EXLCUSIONS %d IN AMBER FILE, AT ATOM #%d\n",
            amber_data->Nnb, i+1);
          NAMD_die(err_msg);
        }
        // There's some 0 in the ExclAt[] list, which is strange
        // and redundant. In this case, I simply ignore such entries.
        if (amber_data->ExclAt[current_index] != 0)
        { // Subtract 1 to convert the index from the 1 to NumAtoms
          // used in the file to the 0 to NumAtoms-1 that we need
          a2 = amber_data->ExclAt[current_index] - 1;
          if (a2 < i)
          { // I assume the latter index be larger than the former
            // one, so that the same exclusion won't be double-counted;
            // if not, give error
            char err_msg[128];
            sprintf(err_msg, "Atom #%d has exclusion with atom #%d, in reverse order.", i+1, a2+1);
            NAMD_die(err_msg);
          }
          else if (a2 == i)
          { char err_msg[128];
            sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN AMBER FILE\n", i+1);
              NAMD_die(err_msg);
          }
          else if (a2 >= numAtoms)
          {  char err_msg[128];
             sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN AMBER FILE",
               a2+1, numAtoms, current_index+1);
             NAMD_die(err_msg);
          }
          exclusions[numExclusions].atom1 = i;
          exclusions[numExclusions].atom2 = a2;
          ++numExclusions;
        }
        ++current_index;
      }
    if (current_index < amber_data->Nnb)
    { char err_msg[128];
      sprintf(err_msg, "Num of exclusions recorded (%d) is smaller than what it's supposed to be (%d)",
        current_index,amber_data->Nnb);
      NAMD_die(err_msg);
    }
  }

  // Copy dihedral information
  numDihedrals = amber_data->Nphih + amber_data->Nphia;
  if (numDihedrals > 0)
  { nphih = amber_data->Nphih;
    dihedrals = new Dihedral[numDihedrals];
    if (dihedrals == NULL || numDihedrals < nphih)
      NAMD_die("memory allocation failed when reading dihedral information");
    // Dihedral WITH hydrogen
    for (i=0; i<nphih; ++i)
    { dihedrals[i].atom1 = amber_data->DihHAt1[i] / 3;
      dihedrals[i].atom2 = amber_data->DihHAt2[i] / 3;
      dihedrals[i].atom3 = amber_data->DihHAt3[i] / 3;
      dihedrals[i].atom4 = amber_data->DihHAt4[i] / 3;
      dihedrals[i].dihedral_type = amber_data->DihHNum[i] - 1;
    }
    // Dihedral WITHOUT hydrogen
    for (i=nphih; i<numDihedrals; ++i)
    { dihedrals[i].atom1 = amber_data->DihAt1[i-nphih] / 3;
      dihedrals[i].atom2 = amber_data->DihAt2[i-nphih] / 3;
      dihedrals[i].atom3 = amber_data->DihAt3[i-nphih] / 3;
      dihedrals[i].atom4 = amber_data->DihAt4[i-nphih] / 3;
      dihedrals[i].dihedral_type = amber_data->DihNum[i-nphih] - 1;
    }
  }
  // In AMBER parm file, dihedrals contain 1-4 exclusion infomation:
  // the 1st and 4th atoms have 1-4 nonbond interation. So we should
  // find them in the exclusion array and change their exclusion to
  // 1-4 type. However, there're two exceptions --
  // 1.If the third atom is negative, it means the end group
  //   interactions are to be ignored;
  // 2.If the fourth atom is negative, it means this is an improper.
  // For the above two cases, the actual atom index is the absolute
  // value of the atom number read; and there's no 1-4 interation
  // for these dihedrals.
  // If readExclusions is not TRUE, then we don't worry about
  // exclusions here.
  for (i=0; i<numDihedrals; ++i)
  { if (dihedrals[i].atom3 < 0 || dihedrals[i].atom4 < 0)
    { dihedrals[i].atom3 = abs(dihedrals[i].atom3);
      dihedrals[i].atom4 = abs(dihedrals[i].atom4);
    }
    else if (simParams->readExclusions)
    { if (dihedrals[i].atom1 < dihedrals[i].atom4)
        a1=dihedrals[i].atom1, a2=dihedrals[i].atom4;
      else
        a1=dihedrals[i].atom4, a2=dihedrals[i].atom1;
      // Since in the exclusion array, atom1 is guaranteed to be
      // ordered, we can do a binary serch to find it first.
      found = 0;
      min=0, max=numExclusions-1;
      while (!found && min<=max)
      { index = (min+max)/2;
        if (exclusions[index].atom1 == a1)
          found = 1;
        else if (exclusions[index].atom1 < a1)
          min = index+1;
        else
          max = index-1;
      }
      if (!found)
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
      // After finding atom1, we do a linear serch to find atom2,
      // in both directions.
      for (j=index-1; j>=0 && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; --j);
      if (j<0 || exclusions[j].atom1!=a1)
        for (j=index; j<numExclusions && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; ++j);
      if (j<numExclusions && exclusions[j].atom1==a1)
        exclusions[j].modified = 1;  // Change the exclusion type to 1-4
      else
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
    }
    if (dihedrals[i].atom1>=numAtoms || dihedrals[i].atom2>=numAtoms ||
        dihedrals[i].atom3>=numAtoms || dihedrals[i].atom4>=numAtoms ||
        dihedrals[i].dihedral_type>=amber_data->Nptra)
    { char err_msg[128];
      sprintf(err_msg, "DIHEDRAL # %d OVERFLOW IN PARM FILE", i+1);
      NAMD_die(err_msg);
    }
  }
  
  //  analyze the data and find the status of each atom
  build_atom_status();

}
/*      END OF FUNCTION read_parm    */


/************************************************************************/
/*                                                                      */
/*      FUNCTION Molecule                                               */
/*                                                                      */
/*  This is the constructor for reading GROMACS topology data           */
/*                                                                      */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param,
		   const GromacsTopFile *gromacsTopFile)
{
  initialize(simParams,param);

  read_parm(gromacsTopFile);
}
/*      END OF FUNCTION Molecule      */

/************************************************************************/
/*                                                                      */
/*      FUNCTION read_parm                                              */
/*                                                                      */
/*   INPUTS:                                                            */
/*  amber_data - AMBER data structure                                   */
/*                                                                      */
/*  This function copys AMBER topology data to the corresponding data   */
/*   structures                                                         */
/*                                                                      */
/************************************************************************/

void Molecule::read_parm(const GromacsTopFile *gf) {
  /*  int i,j,ntheth,nphih,current_index,a1,a2,
      max,min,index,found;*/
  int i;
  
  // Initializes the atom array
  numAtoms = gf->getNumAtoms();
  atoms = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];
  if (atoms == NULL || atomNames == NULL )
    NAMD_die("memory allocation failed when reading atom information");
  ResidueLookupElem *tmpResLookup = resLookup;

  // Copy the individual atoms over
  for (i=0; i<numAtoms; ++i) {
    char *resname = nameArena->getNewArray(11);
    char *atomname = nameArena->getNewArray(11);
    char *atomtype = nameArena->getNewArray(11);
    int resnum,typenum;
    Real charge,mass;

    if (resname == NULL || atomname == NULL || atomtype == NULL)
      NAMD_die("memory allocation failed when reading atom information");

    // get the data out of the GROMACS file
    gf->getAtom(i,&resnum,resname,atomname,atomtype,&typenum,
		&charge,&mass);

    atomNames[i].resname = resname;
    atomNames[i].atomname = atomname;
    atomNames[i].atomtype = atomtype;
    atoms[i].mass = mass;
    atoms[i].charge = charge;
    atoms[i].vdw_type = typenum;

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append("MAIN", resnum+1, i);

    /*  Determine the type of the atom (H or O) */
    // XXX this cannot be done this way in GROMACS
    // For example, in dppc LO2 appears to be an oxygen.
    // And how do the hydrogens in CH3 etc factor in to this?
    atoms[i].status = UnknownAtom; // the default
    if (atoms[i].mass <=3.5) {
      atoms[i].status |= HydrogenAtom;
    } else if ((atomNames[i].atomname[0] == 'O') && 
         (atoms[i].mass >= 14.0) && 
         (atoms[i].mass <= 18.0)) {
      atoms[i].status |= OxygenAtom;
    }
  }

  // Copy bond information
  numBonds = gf->getNumBonds();
  bonds = new Bond[numBonds];
  if (bonds == NULL)
    NAMD_die("memory allocation failed when reading bond information");
  for(i=0;i<numBonds;i++) {
    int type; // to convert the type correctly
    int atom1,atom2;
    gf->getBond(i,&atom1,&atom2,&type);
    bonds[i].atom1 = atom1;
    bonds[i].atom2 = atom2;
    bonds[i].bond_type = (Index)type;
  }

  // Copy angle information
  numAngles = gf->getNumAngles();
  angles = new Angle[numAngles];
  if (angles == NULL)
    NAMD_die("memory allocation failed when reading angle information");
  for(i=0;i<numAngles;i++) {
    int type; // to convert the type correctly
    int atom1,atom2,atom3;
    gf->getAngle(i,&atom1,&atom2,&atom3,&type);

    angles[i].atom1 = atom1;
    angles[i].atom2 = atom2;
    angles[i].atom3 = atom3;

    angles[i].angle_type=type;
  }

  numExclusions =  0;
  exclusions = new Exclusion[numExclusions];

  /*
  // If readExclusions is TRUE, then we copy exclusions from parm
  // file; otherwise we skip the exclusions here and generate
  // them later in build_exclusions()
  if (simParams->readExclusions)
  { // Copy exclusion information
    // In Amber data structure, Iblo[] is the number of exclusions
    // for each atom; ExclAt[] is the atom index for the excluded atoms.
    exclusions = new Exclusion[amber_data->Nnb];
    if (exclusions == NULL &&  amber_data->Nnb > 0)
      NAMD_die("memory allocation failed when reading exclusion information");
    current_index = 0;
    for (i=0; i<numAtoms; ++i)
      for (j=0; j<amber_data->Iblo[i]; ++j)
      { if (current_index >= amber_data->Nnb)
        { char err_msg[128];
          sprintf(err_msg, "EXCLUSION INDEX EXCEEDS NUMBER OF EXLCUSIONS %d IN AMBER FILE, AT ATOM #%d\n",
            amber_data->Nnb, i+1);
          NAMD_die(err_msg);
        }
        // There's some 0 in the ExclAt[] list, which is strange
        // and redundant. In this case, I simply ignore such entries.
        if (amber_data->ExclAt[current_index] != 0)
        { // Subtract 1 to convert the index from the 1 to NumAtoms
          // used in the file to the 0 to NumAtoms-1 that we need
          a2 = amber_data->ExclAt[current_index] - 1;
          if (a2 < i)
          { // I assume the latter index be larger than the former
            // one, so that the same exclusion won't be double-counted;
            // if not, give error
            char err_msg[128];
            sprintf(err_msg, "Atom #%d has exclusion with atom #%d, in reverse order.", i+1, a2+1);
            NAMD_die(err_msg);
          }
          else if (a2 == i)
          { char err_msg[128];
            sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN AMBER FILE\n", i+1);
              NAMD_die(err_msg);
          }
          else if (a2 >= numAtoms)
          {  char err_msg[128];
             sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN AMBER FILE",
               a2+1, numAtoms, current_index+1);
             NAMD_die(err_msg);
          }
          exclusions[numExclusions].atom1 = i;
          exclusions[numExclusions].atom2 = a2;
          ++numExclusions;
        }
        ++current_index;
      }
    if (current_index < amber_data->Nnb)
    { char err_msg[128];
      sprintf(err_msg, "Num of exclusions recorded (%d) is smaller than what it's supposed to be (%d)",
        current_index,amber_data->Nnb);
      NAMD_die(err_msg);
    }
  }
  */

  // Copy dihedral information
  numDihedrals = gf->getNumDihedrals();
  dihedrals = new Dihedral[numDihedrals];
  if (dihedrals == NULL)
    NAMD_die("memory allocation failed when reading dihedral information");
  for(i=0;i<numDihedrals;i++) {
    int type; // to convert the type correctly
    int atom1,atom2,atom3,atom4;
    gf->getDihedral(i,&atom1,&atom2,&atom3,&atom4,&type);
    dihedrals[i].atom1 = atom1;
    dihedrals[i].atom2 = atom2;
    dihedrals[i].atom3 = atom3;
    dihedrals[i].atom4 = atom4;
    dihedrals[i].dihedral_type = type;
  }

  /*
  // In AMBER parm file, dihedrals contain 1-4 exclusion infomation:
  // the 1st and 4th atoms have 1-4 nonbond interation. So we should
  // find them in the exclusion array and change their exclusion to
  // 1-4 type. However, there're two exceptions --
  // 1.If the third atom is negative, it means the end group
  //   interactions are to be ignored;
  // 2.If the fourth atom is negative, it means this is an improper.
  // For the above two cases, the actual atom index is the absolute
  // value of the atom number read; and there's no 1-4 interation
  // for these dihedrals.
  // If readExclusions is not TRUE, then we don't worry about
  // exclusions here.
  for (i=0; i<numDihedrals; ++i)
  { if (dihedrals[i].atom3 < 0 || dihedrals[i].atom4 < 0)
    { dihedrals[i].atom3 = abs(dihedrals[i].atom3);
      dihedrals[i].atom4 = abs(dihedrals[i].atom4);
    }
    else if (simParams->readExclusions)
    { if (dihedrals[i].atom1 < dihedrals[i].atom4)
        a1=dihedrals[i].atom1, a2=dihedrals[i].atom4;
      else
        a1=dihedrals[i].atom4, a2=dihedrals[i].atom1;
      // Since in the exclusion array, atom1 is guaranteed to be
      // ordered, we can do a binary serch to find it first.
      found = 0;
      min=0, max=numExclusions-1;
      while (!found && min<=max)
      { index = (min+max)/2;
        if (exclusions[index].atom1 == a1)
          found = 1;
        else if (exclusions[index].atom1 < a1)
          min = index+1;
        else
          max = index-1;
      }
      if (!found)
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
      // After finding atom1, we do a linear serch to find atom2,
      // in both directions.
      for (j=index-1; j>=0 && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; --j);
      if (j<0 || exclusions[j].atom1!=a1)
        for (j=index; j<numExclusions && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; ++j);
      if (j<numExclusions && exclusions[j].atom1==a1)
        exclusions[j].modified = 1;  // Change the exclusion type to 1-4
      else
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
    }
    if (dihedrals[i].atom1>=numAtoms || dihedrals[i].atom2>=numAtoms ||
        dihedrals[i].atom3>=numAtoms || dihedrals[i].atom4>=numAtoms ||
        dihedrals[i].dihedral_type>=amber_data->Nptra)
    { char err_msg[128];
      sprintf(err_msg, "DIHEDRAL # %d OVERFLOW IN PARM FILE", i+1);
      NAMD_die(err_msg);
    }
  }
  */
  //  analyze the data and find the status of each atom

  build_atom_status();

}
/*      END OF FUNCTION read_parm    */
