/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   This class is used to store all of the structural       
   information for a simulation.  It reas in this information 
   from a .psf file, cross checks and obtains some information
   from the Parameters object that is passed in, and then     
   stores all this information for later use.	
*/


#ifndef MOLECULE_H

#define MOLECULE_H

#include "common.h"
#include "NamdTypes.h"
#include "IntList.h"
#include "structures.h"
#include "ConfigList.h"
#include "Vector.h"
#include "UniqueSet.h"
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
template<class Type> class ObjectArena;

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
	ObjectArena<char> *nameArena;
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

	ObjectArena<int> *tmpArena;
	int **bondsWithAtom;	//  List of bonds involving each atom

	ObjectArena<int> *arena;
	int **bondsByAtom;	//  List of bonds owned by each atom
	int **anglesByAtom;     //  List of angles owned by each atom
	int **dihedralsByAtom;  //  List of dihedrals owned by each atom
	int **impropersByAtom;  //  List of impropers owned by each atom
	int **exclusionsByAtom; //  List of exclusions owned by each atom

	int **all_exclusions;
				//  List of all exclusions, including
				//  explicit exclusions and those calculated
				//  from the bonded structure based on the
				//  exclusion policy.  Also includes 1-4's.

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
	void stripHGroupExcl(void);
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
	int numRigidBonds;	//  Number of rigid bonds
	int numFixedRigidBonds; //  Number of rigid bonds between fixed atoms

	// The following are needed for error checking because we
	// eliminate bonds, etc. which involve only fixed atoms
	int numCalcBonds;	//  Number of bonds requiring calculation
	int numCalcAngles;	//  Number of angles requiring calculation
	int numCalcDihedrals;	//  Number of dihedrals requiring calculation
	int numCalcImpropers;	//  Number of impropers requiring calculation
	int numCalcExclusions;	//  Number of exclusions requiring calculation

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

	void build_langevin_params(BigReal coupling, Bool doHydrogen);
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

	//  Lookup number of atoms in residue from segment and residue
	int get_residue_size(const char *segid, int resid) const;

	//  Lookup atom id from segment, residue, and index in residue
	int get_atom_from_index_in_residue(const char *segid, int resid, int index) const;

	
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
        //  Returns 1 for full, 2 for 1-4 exclusions.
	int checkexcl(int atom1, int atom2) const;

	int *get_excl_check_for_atom(int anum) const
			 { return all_exclusions[anum]; }

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
	// < 0 implies MOLLY but not SHAKE, > 1 implies both if MOLLY is on
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

