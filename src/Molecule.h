//-*-c++-*-
/***************************************************************************/
/*    (C) Copyright 1995,1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *	This class is used to store all of the structural       
 *   information for a simulation.  It reas in this information 
 *   from a .psf file, cross checks and obtains some information
 *   from the Parameters object that is passed in, and then     
 *   stores all this information for later use.	
 *
 ***************************************************************************/


#ifndef MOLECULE_H

#define MOLECULE_H

#include "common.h"
#include "NamdTypes.h"
#include "IntList.h"
#include "LintList.h"
#include "structures.h"
#include "ConfigList.h"
#include "Vector.h"
#include "UniqueSet.h"
#include "ObjectArena.h"
#include "Hydrogen.h"

class SimParameters;
class Parameters;
class PDB;
class MIStream;
class MOStream;

class BondElem;
class AngleElem;
class DihedralElem;
class ImproperElem;
class NonbondedExclElem;
class ResidueLookupElem;

typedef int* intPtr;

// List maintaining the global atom indicies sorted by helix groups.
class Molecule
{
typedef struct constraint_params
{
   Real k;		//  Force constant
   Vector refPos;	//  Reference position for restraint
} ConstraintParams;

friend class BondElem;
friend class AngleElem;
friend class DihedralElem;
friend class ImproperElem;
friend class NonbondedExclElem;

private:
	Atom *atoms;		//  Array of atom structures
	AtomNameInfo *atomNames;//  Array of atom name info.  Only maintained
				//  on node 0 for VMD interface
	ResidueLookupElem *resLookup; // find residues by name
	Bond *bonds;		//  Array of bond structures
	Angle *angles;		//  Array of angle structures
	Dihedral *dihedrals;	//  Array of dihedral structures
	Improper *impropers;	//  Array of improper structures
	Bond *donors;	        //  Array of hydrogen bond donor structures
	Bond *acceptors;	//  Array of hydrogen bond acceptor
	Exclusion *exclusions;	//  Array of exclusion structures
	UniqueSet<Exclusion> exclusionSet;  //  Used for building
	int *consIndexes;	//  Constraint indexes for each atom
	ConstraintParams *consParams;
				//  Parameters for each atom constrained
	Real *langevinParams;   //  b values for langevin dynamics
	Real *langForceVals;    //  Calculated values for langvin random forces
	int *fixedAtomFlags;	//  1 for fixed, -1 for fixed group, else 0
	Real *rigidBondLengths;  //  if H, length to parent or 0. or
				//  if not H, length between children or 0.

	ObjectArena<int> arena;
	int **bondsByAtom;	//  List of bonds involving each atom
	int **anglesByAtom;  //  List of angles involving each atom
	int **dihedralsByAtom;
				//  List of dihedrals by atom
	int **impropersByAtom;
				//  List of impropers by atom
	int **exclusionsByAtom;
				//  List of exclusions by atom

	int **all_exclusions;
				//  List of all exclusions, including
				//  explicit exclusions and those calculated
				//  from the bonded structure based on the
				//  exclusion policy
	int **onefour_exclusions;
				//  List of 1-4 interactions.  This list is
				//  used only if the exclusion policy is 
				//  scaled1-4 to track 1-4 interactions that
				//  need to be handled differently

	void build_lists_by_atom();
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

	void build12excl(void);
	void build13excl(void);
	void build14excl(int);
	void build_exclusions();

	// analyze the atoms, and determine which are oxygen, hb donors, etc.
	// this is called after a molecule is sent our (or received in)
	void build_atom_status(void);

	// added during the trasition from 1x to 2
	SimParameters *simParams;
	Parameters *params;

public:
	int numAtoms;		//  Number of atoms 
	int numBonds;		//  Number of bonds
	int numAngles;		//  Number of angles
	int numDihedrals;	//  Number of dihedrals
	int numImpropers;	//  Number of impropers
	int numDonors;	        //  Number of hydrogen bond donors
	int numAcceptors;	//  Number of hydrogen bond acceptors
	int numExclusions;	//  Number of exclusions
	int numTotalExclusions; //  Real Total Number of Exclusions // hack
	int numConstraints;	//  Number of atoms constrained
	int numFixedAtoms;	//  Number of fixed atoms
	int numHydrogenGroups;	//  Number of hydrogen groups

	//  Number of dihedrals with multiple periodicity
	int numMultipleDihedrals; 
	//  Number of impropers with multiple periodicity
	int numMultipleImpropers; 
	// indexes of "atoms" sorted by hydrogen groups
	HydrogenGroup hydrogenGroup;
	int waterIndex;

	Molecule(SimParameters *, Parameters *param, char *filename=NULL);
	~Molecule();		//  Destructor

	void read_psf_file(char *, Parameters *);
				//  Read in a .psf file given
				//  the filename and the parameter
				//  object to use
	void send_Molecule(Communicate *);
				//  send the molecular structure 
				//  from the master to the clients
	void receive_Molecule(MIStream *);
				//  receive the molecular structure
				//  from the master on a client
	
	void build_constraint_params(StringList *, StringList *, StringList *,
				     PDB *, char *);
				//  Build the set of harmonic constraint 
				// parameters

	void build_langevin_params(StringList *, StringList *, PDB *, char *);
				//  Build the set of langevin dynamics parameters

	void build_fixed_atoms(StringList *, StringList *, PDB *, char *);
				//  Determine which atoms are fixed (if any)

        Bool is_hydrogen(int);     // return true if atom is hydrogen
        Bool is_oxygen(int);       // return true if atom is oxygen
        Bool is_hb_donor(int);     // return true if atom is hbond donor
        Bool is_hb_acceptor(int);  // return true if atom is hbond acceptor
        Bool is_hb_antecedent(int);// return true if atom is hbond antecedent
        Bool is_hb_hydrogen(int);  // return true if atom is hbond hydrogen
	Bool is_hydrogenGroupParent(int); // return true if atom is group parent
	Bool is_water(int);        // return true if atom is part of water 
	int  get_groupSize(int);     // return # atoms in (hydrogen) group
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

	//  Retrieve a hydrogen bond donor structure
	Bond *get_donor(int dnum) const {return (&(donors[dnum]));}

	//  Retrieve a hydrogen bond acceptor structure
	Bond *get_acceptor(int dnum) const {return (&(acceptors[dnum]));}

	//  Retrieve an exclusion structure
	Exclusion *get_exclusion(int ex) const {return (&(exclusions[ex]));}

	//  Retrieve an atom type
	const char *get_atomtype(int anum) const
	{
		if (atomNames == NULL)
		{
			NAMD_die("Tried to find atom type on node other than node 0");
		}

		return(atomNames[anum].atomtype);
	}

	//  Lookup atom id from segment, residue, and name
	int get_atom_from_name(const char *segid, int resid, const char *aname) const;
	
	//  The following routines are used to get the list of bonds
	//  for a given atom.  This is used when creating the bond lists
	//  for the force objects
	int *get_bonds_for_atom(int anum) { return bondsByAtom[anum]; }
	int *get_angles_for_atom(int anum) 
			{ return anglesByAtom[anum]; }
	int *get_dihedrals_for_atom(int anum) 
			{ return dihedralsByAtom[anum]; }
	int *get_impropers_for_atom(int anum) 
			{ return impropersByAtom[anum]; }
	int *get_exclusions_for_atom(int anum)
			{ return exclusionsByAtom[anum]; }
	
	//  Check for exclusions, either explicit or bonded.
	//  Inline this funcion since it is called so often
	Bool checkexcl(int atom1, int atom2) const
        {
	   register int check_int;	//  atom whose array we will search
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
	   register int *list = all_exclusions[check_int];
	   check_int = *list;
	   while( check_int != other_int && check_int != -1 )
	   {
	      check_int = *(++list);
	   }
	   return ( check_int != -1 );
        }
	
	//  Check for 1-4 exclusions.  This is only valid when the
	//  exclusion policy is set to scaled1-4. Inline this function
	//  since it will be called so often
	Bool check14excl(int atom1, int atom2) const
        {
	   register int check_int;	//  atom whose array we will search
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
	   register int *list = onefour_exclusions[check_int];
	   check_int = *list;
	   while( check_int != other_int && check_int != -1 )
	   {
	      check_int = *(++list);
	   }
	   return ( check_int != -1 );
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

	Bool is_atom_fixed(int atomnum) const
	{
		return (numFixedAtoms && fixedAtomFlags[atomnum]);
	}

	Bool is_group_fixed(int atomnum) const
	{
		return (numFixedAtoms && (fixedAtomFlags[atomnum] == -1));
	}

	// 0 if not rigid or length to parent, for parent refers to H-H length
	Real rigid_bond_length(int atomnum) const
	{
		return(rigidBondLengths[atomnum]);
	}

	void print_atoms(Parameters *);	
				//  Print out list of atoms
	void print_bonds(Parameters *);	
				//  Print out list of bonds
	void print_exclusions();//  Print out list of exclusions

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Molecule.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1017 $	$Date: 1998/02/17 06:39:22 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Molecule.h,v $
 * Revision 1.1017  1998/02/17 06:39:22  jim
 * SHAKE/RATTLE (rigidBonds) appears to work!!!  Still needs langevin,
 * proper startup, and degree of freedom tracking.
 *
 * Revision 1.1016  1998/02/11 07:27:24  jim
 * Completed interface to free energy perturbation, added code for
 * determining atomid from segid, resid, and atomname.
 *
 * Revision 1.1015  1997/12/26 23:10:52  milind
 * Made namd2 to compile, link and run under linux. Merged Templates and src
 * directoriies, and removed separate definition and declaration files for
 * templates.
 *
 * Revision 1.1014  1997/10/17 17:16:50  jim
 * Switched from hash tables to checklists, eliminated special exclusion code.
 *
 * Revision 1.1013  1997/10/01 16:46:57  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.1012  1997/09/21 21:58:31  jim
 * Added printing of hydrogen group count.
 *
 * Revision 1.1011  1997/09/19 08:55:33  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1010  1997/04/03 19:59:05  nealk
 * 1) New Fopen() which handles .Z and .gz files.
 * 2) localWaters and localNonWaters lists on each patch.
 *
 * Revision 1.1009  1997/03/31 16:12:54  nealk
 * Atoms now can migrate by hydrogen groups.
 *
 * Revision 1.1008  1997/03/27 17:08:30  nealk
 * Added hydrogen groupings.  Now configuration parameter "splitPatch" determines
 * atom-into-patch distribution.
 *
 * Revision 1.1007  1997/03/19 18:10:15  nealk
 * Added sorted hydrogen group list to molecule.
 *
 * Revision 1.1006  1997/03/19 11:54:33  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1005  1997/03/14 21:40:10  ari
 * Reorganized startup to make possible inital load
 * balancing by changing methods in WorkDistrib.
 * Also made startup more transparent and easier
 * to modify.
 *
 * Revision 1.1004  1997/03/11 23:46:31  ari
 * Improved ComputeNonbondedExcl loadTuples() by overloading the default
 * template method from ComputeHomeTuples and used the checklist suggested
 * by Jim.  Good performance gain.
 *
 * Revision 1.1003  1997/03/11 06:29:14  jim
 * Modified exclusion lookup to use fixed arrays and linear searches.
 *
 * Revision 1.1002  1997/03/11 04:07:55  jim
 * Eliminated use of LintList for by-atom lists.
 * Now uses little arrays managed by ObjectArena<int>.
 *
 * Revision 1.1001  1997/02/10 08:14:36  jim
 * Fixed problem with exclusions where both modified and unmodified
 * versions of the same exclusion could be placed in the list, causing
 * one to be selected more or less randomly.  Also caused different
 * results on different numbers of processors.
 *
 * Revision 1.1000  1997/02/06 15:58:45  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:54  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:25  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:29  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.9  1996/12/06 17:28:15  jim
 * put build_lists_by_atom back where it belongs
 *
 * Revision 1.8  1996/12/06 06:54:41  jim
 * started from 1.4, added support for treating exclusions like bonds
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
