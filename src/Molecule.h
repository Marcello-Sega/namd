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
 *	$RCSfile: Molecule.h,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.5 $	$Date: 1996/12/03 17:50:13 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *	This class is used to store all of the structural       
 *   information for a simulation.  It reas in this information 
 *   from a .psf file, cross checks and obtains some information
 *   from the Parameters object that is passed in, and then     
 *   stores all this information for later use.	
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Molecule.h,v $
 * Revision 1.5  1996/12/03 17:50:13  nealk
 * Added nonbonded excl stuff.
 *
 * Revision 1.4  1996/11/21 23:36:04  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/11/05 05:01:23  jim
 * added some consts
 *
 * Revision 1.2  1996/08/16 04:39:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.22  1996/04/23 22:44:48  billh
 * Modified Atom structure to store a list of the hydrogen bond pairs
 * each atom participates in.
 *
 * Revision 1.21  1996/04/18 18:46:18  billh
 * Updated to read hydrogen bond information.
 *
 * Revision 1.20  1996/02/16 16:37:04  gursoy
 * added atomType array and functions to allow is_hydrogen and is_oxygen
 * queries for an atom.
 *
 * Revision 1.19  95/04/06  12:52:03  12:52:03  nelson (Mark T. Nelson)
 * Removed extern class references
 * 
 * Revision 1.18  95/03/08  14:46:22  14:46:22  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.17  95/02/22  14:55:45  14:55:45  nelson (Mark T. Nelson)
 * Made changes for current working directory parameter
 * 
 * Revision 1.16  95/01/26  14:17:38  14:17:38  nelson (Mark T. Nelson)
 * Made additions for Charm22 parameters
 * 
 * Revision 1.15  95/01/19  15:29:03  15:29:03  nelson (Mark T. Nelson)
 * Added langevin dynamics parameters
 * 
 * Revision 1.14  94/10/28  12:47:28  12:47:28  nelson (Mark T. Nelson)
 * Added routine to get atom types for VMD
 * 
 * Revision 1.13  94/10/19  21:43:17  21:43:17  nelson (Mark T. Nelson)
 * Added functions and data structures for harmonic constraints
 * 
 * Revision 1.12  94/10/12  09:38:31  09:38:31  nelson (Mark T. Nelson)
 * Changed the way exclusions are handled so that bonded
 * exclusions are now pre-calculated and all exclusions
 * are stored in arrays by atom that can be searched via
 * binary searches.
 * 
 * Revision 1.11  94/10/08  10:28:35  10:28:35  nelson (Mark T. Nelson)
 * Added routines to get bonded lists by atom
 * 
 * Revision 1.10  94/10/04  10:35:03  10:35:03  nelson (Mark T. Nelson)
 * Changed 12, 13, and 14 search routines
 * 
 * Revision 1.9  94/09/28  10:00:28  10:00:28  nelson (Mark T. Nelson)
 * Added get_improper and get_dihedral functions
 * 
 * Revision 1.8  94/09/24  20:11:37  20:11:37  nelson (Mark T. Nelson)
 * Added routines to check for 1-2, 1-3, and 1-4 interactions
 * 
 * Revision 1.7  94/09/13  14:33:19  14:33:19  gursoy (Attila Gursoy)
 * receive_Molecule gets the message from Node object (charm++ integration)
 * 
 * Revision 1.6  94/09/13  13:47:12  13:47:12  nelson (Mark T. Nelson)
 * Added get_angle function
 * 
 * Revision 1.5  94/09/04  20:20:21  20:20:21  nelson (Mark T. Nelson)
 * Added get_bond function
 * 
 * Revision 1.4  94/08/30  13:56:45  13:56:45  nelson (Mark T. Nelson)
 * added routines atommass and atomcharge
 * 
 * Revision 1.3  94/08/04  08:39:36  08:39:36  nelson (Mark T. Nelson)
 * Added send_Molecule and receive_Molecule functions
 * 
 * Revision 1.2  94/07/07  13:31:23  13:31:23  nelson (Mark T. Nelson)
 * Changed due to changes in Parameters class
 * 
 * Revision 1.1  94/06/22  15:04:16  15:04:16  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef MOLECULE_H

#define MOLECULE_H

#include "common.h"
#include "IntList.h"
#include "LintList.h"
#include "structures.h"
#include "ConfigList.h"
#include "Vector.h"

class SimParameters;
class Parameters;
class PDB;
class Message;

class Molecule
{
typedef struct constraint_params
{
   Real k;		//  Force constant
   Vector refPos;	//  Reference position for restraint
} ConstraintParams;

private:
	Atom *atoms;		//  Array of atom structures
	AtomNameInfo *atomNames;//  Array of atom name info.  Only maintained
				//  on node 0 for VMD interface
	Bond *bonds;		//  Array of bond structures
	Angle *angles;		//  Array of angle structures
	Dihedral *dihedrals;	//  Array of dihedral structures
	Improper *impropers;	//  Array of improper structures
	NonbondedExcl *nonbondedexcls;	// Array of Nonbonded Excl structures
	Bond *donors;	        //  Array of hydrogen bond donor structures
	Bond *acceptors;	//  Array of hydrogen bond acceptor
	Exclusion *exclusions;	//  Array of exclusion structures
	int *consIndexes;	//  Constraint indexes for each atom
	ConstraintParams *consParams;
				//  Parameters for each atom constrained
	Real *langevinParams;   //  b values for langevin dynamics
	Real *langForceVals;    //  Calculated values for langvin random forces

	LintList *bondsByAtom;	//  List of bonds involving each atom
	LintList *anglesByAtom;	//  List of angles involving each atom
	LintList *dihedralsByAtom;
				//  List of dihedrals by atom
	LintList *impropersByAtom;
				//  List of impropers by atom
	LintList *nonbondedexclsByAtom;
				//  List of nonbonded excls involving each atom

	IntList *all_exclusions;
				//  List of all exclusions, including
				//  explicit exclusions and those calculated
				//  from the bonded structure based on the
				//  exclusion policy
	IntList *onefour_exclusions;
				//  List of 1-4 interactions.  This list is
				//  used only if the exclusion policy is 
				//  scaled1-4 to track 1-4 interactions that
				//  need to be handled differently

//	void build_lists_by_atom();
				//  Build the list of structures by atom
	

	void read_atoms(FILE *, Parameters *);
				//  Read in atom info from .psf
	void read_bonds(FILE *, Parameters *);
				//  Read in bond info from .psf
	void read_angles(FILE *, Parameters *);
				//  Read in angle info from .psf
	void read_dihedrals(FILE *, Parameters *);
				//  Read in dihedral info from .psf
	void read_impropers(FILE *, Parameters *);
				//  Read in improper info from .psf
	void read_donors(FILE *);
				//  Read in hydrogen bond donors from .psf
	void read_acceptors(FILE *);
				//  Read in hydrogen bond acceptors from .psf
	void read_exclusions(FILE *);
				//  Read in exclusion info from .psf

	void build12excl(IntList *);
	void build13excl(IntList *);
	void build14excl(IntList *);
	void build_exclusions();

	// analyze the atoms, and determine which are oxygen, hb donors, etc.
	// this is called after a molecule is sent our (or received in)
	void build_atom_status(void);

	// added during the trasition from 1x to 2
	SimParameters *simParams;

public:
	int numAtoms;		//  Number of atoms 
	int numBonds;		//  Number of bonds
	int numAngles;		//  Number of angles
	int numDihedrals;	//  Number of dihedrals
	int numImpropers;	//  Number of impropers
	int numNonbondedExcls;	//  Number of nonbonded excls
	int numDonors;	        //  Number of hydrogen bond donors
	int numAcceptors;	//  Number of hydrogen bond acceptors
	int numExclusions;	//  Number of exclusions
	int numConstraints;	//  Number of atoms constrained

	//  Number of dihedrals with multiple periodicity
	int numMultipleDihedrals; 
	//  Number of impropers with multiple periodicity
	int numMultipleImpropers; 

	Molecule(SimParameters *, Parameters *param=NULL, char *filename=NULL);
	~Molecule();		//  Destructor

	// We Did This - ari and jim
	void build_lists_by_atom();
				//  Build the list of structures by atom

	void read_psf_file(char *, Parameters *);
				//  Read in a .psf file given
				//  the filename and the parameter
				//  object to use
	void send_Molecule(Communicate *);
				//  send the molecular structure 
				//  from the master to the clients
	void receive_Molecule(Message *);
				//  receive the molecular structure
				//  from the master on a client
	
	void build_constraint_params(StringList *, StringList *, StringList *,
				     PDB *, char *);
				//  Build the set of harmonic constraint 
				// parameters

	void build_langevin_params(StringList *, StringList *, PDB *, char *);
				//  Build the set of langevin dynamics parameters

        Bool is_hydrogen(int);     // return true if atom is hydrogen
        Bool is_oxygen(int);       // return true if atom is oxygen
        Bool is_hb_donor(int);     // return true if atom is hbond donor
        Bool is_hb_acceptor(int);  // return true if atom is hbond acceptor
        Bool is_hb_antecedent(int);// return true if atom is hbond antecedent
        Bool is_hb_hydrogen(int);  // return true if atom is hbond hydrogen

        int get_mother_atom(int);  // return mother atom of a hydrogen

	IntList *get_atom_hb_donors(int);    // return list of Nth atom donors
	IntList *get_atom_hb_acceptors(int); // return list of Nth atom acc's

	//  Get the mass of an atom
	Real atommass(int anum) const
	{
		return(atoms[anum].mass);
	}

	//  Get the charge of an atom
	Real atomcharge(int anum) const
	{
		return(atoms[anum].charge);
	}
	
	//  Get the vdw type of an atom
	Index atomvdwtype(int anum) const
	{
	   	return(atoms[anum].vdw_type);
	}

	//  Retrieve a bond structure
	Bond *get_bond(int bnum) const {return (&(bonds[bnum]));}

	//  Retrieve an angle structure
	Angle *get_angle(int anum) const {return (&(angles[anum]));}

	//  Retrieve an improper strutcure
	Improper *get_improper(int inum) const {return (&(impropers[inum]));}

	//  Retrieve a dihedral structure
	Dihedral *get_dihedral(int dnum) const {return (&(dihedrals[dnum]));}

	//  Retrieve a nonbonded excl structure
	NonbondedExcl *get_nonbondedexcl(int nnum) const
		{return (&(nonbondedexcls[nnum]));}

	//  Retrieve a hydrogen bond donor structure
	Bond *get_donor(int dnum) const {return (&(donors[dnum]));}

	//  Retrieve a hydrogen bond acceptor structure
	Bond *get_acceptor(int dnum) const {return (&(acceptors[dnum]));}

	//  Retrieve an atom type
	const char *get_atomtype(int anum) const
	{
		if (atomNames == NULL)
		{
			NAMD_die("Tried to find atom type on node other than node 0");
		}

		return(atomNames[anum].atomtype);
	}
	
	//  The following routines are used to get the list of bonds
	//  for a given atom.  This is used when creating the bond lists
	//  for the force objects
	LintList *get_bonds_for_atom(int anum) {return (&(bondsByAtom[anum]));}
	LintList *get_angles_for_atom(int anum) 
			{return (&(anglesByAtom[anum]));}
	LintList *get_dihedrals_for_atom(int anum) 
			{return (&(dihedralsByAtom[anum]));}
	LintList *get_impropers_for_atom(int anum) 
			{return (&(impropersByAtom[anum]));}
	LintList *get_nonbondedexcls_for_atom(int anum)
			{return (&(nonbondedexclsByAtom[anum]));}
	
	//  Check for exclusions, either explicit or bonded.
	//  Inline this funcion since it is called so often
	Bool checkexcl(int atom1, int atom2) const
        {
	   int check_int;	//  atom whose array we will search
	   int other_int;	//  atom we are looking for

	   //  We want to search the array of the smaller atom
	   if (atom1<atom2)
	   {
		check_int = atom1;
		other_int = atom2;
	   }
	   else
	   {
		check_int = atom2;
		other_int = atom1;
	   }

	   //  Do the search and return the correct value
	   if (all_exclusions[check_int].find(other_int) == INTLIST_NOTFOUND)
	   {
		return(FALSE);
	   }
	   else
	   {
		return(TRUE);
	   }
        }
	
	//  Check for 1-4 exclusions.  This is only valid when the
	//  exclusion policy is set to scaled1-4. Inline this function
	//  since it will be called so often
	Bool check14excl(int atom1, int atom2) const
        {
	   int check_int;
	   int other_int;

	   if (atom1<atom2)
	   {
		check_int = atom1;
		other_int = atom2;
	   }
	   else
	   {
		check_int = atom2;
		other_int = atom1;
	   }

	   if (onefour_exclusions[check_int].find(other_int) == INTLIST_NOTFOUND)
	   {
		return(FALSE);
	   }
	   else
	   {
		return(TRUE);
	   }
	}

	//  Return true or false based on whether the specified atom
	//  is constrained or not.
	Bool is_atom_constrained(int atomnum) const
	{
		if (numConstraints)
		{
			//  Check the index to see if it is constrained
			return(consIndexes[atomnum] != -1);
		}
		else
		{
			//  No constraints at all, so just return FALSE
			return(FALSE);
		}
	}

	//  Get the harmonic constraints for a specific atom
	void get_cons_params(Real &k, Vector &refPos, int atomnum) const
	{
		k = consParams[consIndexes[atomnum]].k;
		refPos = consParams[consIndexes[atomnum]].refPos;
	}

	Real langevin_param(int atomnum) const
	{
		return(langevinParams[atomnum]);
	}

	Real langevin_force_val(int atomnum) const
	{
		return(langForceVals[atomnum]);
	}

	void print_atoms(Parameters *);	
				//  Print out list of atoms
	void print_bonds(Parameters *);	
				//  Print out list of bonds
	void print_exclusions();//  Print out list of exclusions

};

#endif
