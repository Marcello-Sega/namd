/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   The Parameters class is used to read in and store all of the parameters
   from the parameter files.  These parameters are stored and used to assign
   constants to bonds and atoms as they are read in from the psf file
*/

#ifndef PARAM_H

#define PARAM_H

#include "parm.h"

#include "common.h"
#include "structures.h"
#include "strlib.h"
#include "MStream.h"
//****** BEGIN CHARMM/XPLOR type changes
#include "SimParameters.h"
//****** END CHARMM/XPLOR type changes

#include "GromacsTopFile.h"

class Communicate;
class StringList;

//  The class Parameters is used to store and query the parameters for
//  bonds and atoms.  If the Parameter object resides on the master
//  process, it is responsible for reading in all the parameters and
//  then communicating them to the other processors.  To do this, first
//  the routine read_paramter_file is called to read in each parameter
//  file.  The information that is read in is stored in binary trees
//  for vdw, bond, and angle information and in linked lists for 
//  dihedrals and impropers.  Once all of the files have been read
//  in, the routine done_reading_files is called.  At this point, all
//  of the data that has been read in is copied to arrays.  This is
//  so that each bond and atom can be assigned an index into these
//  arrays to retreive the parameters in constant time.
//
//  Then the psf file is read in.  Each bond and atom is assigned an
//  index into the parameter lists using the functions assign_*_index.
//  Once the psf file has been read in, then the routine 
//  done_reading_structre is called to tell the object that it no
//  longer needs to store the binary trees and linked lists that were
//  used to query the parameters based on atom type.  From this point
//  on, only the indexes will be used.
//
//  The master node then uses send_parameters to send all of these
//  parameters to the other nodes and the objects on all of the other
//  nodes use receive_parameters to accept these parameters.
//
//  From this point on, all of the parameter data is static and the
//  functions get_*_params are used to retrieve the parameter data
//  that is desired.

//  Define the number of Multiples that a Dihedral or Improper
//  bond can have with the Charm22 parameter set
#define MAX_MULTIPLICITY	6

//  Number of characters maximum allowed for storing atom type names
#define MAX_ATOMTYPE_CHARS      6

//****** BEGIN CHARMM/XPLOR type changes
//  Define the numbers associated with each possible parameter-file 
//  type (format) NAMD can handle.
#define paraXplor               0
#define paraCharmm              1
//****** END CHARMM/XPLOR type changes


class BondValue {
public:
	Real k;		//  Force constant for the bond
	Real x0;	//  Rest distance for the bond
};

class AngleValue {
public:
	Real k;		//  Force constant for angle
	Real theta0;	//  Rest angle for angle
	Real k_ub;	//  Urey-Bradley force constant
	Real r_ub;	//  Urey-Bradley distance
};

typedef struct four_body_consts
{
	Real k;		//  Force constant
	int n;		//  Periodicity
	Real delta;	//  Phase shift
} FourBodyConsts;

class DihedralValue {
public:
	int multiplicity;
	FourBodyConsts values[MAX_MULTIPLICITY];
};

class ImproperValue {
public:
	int multiplicity;
	FourBodyConsts values[MAX_MULTIPLICITY];
};

class NonbondedExclValue {
public:
	// need to put parameters here...
	// for now, copy bond
	Real k;		//  Force constant for the bond
	Real x0;	//  Rest distance for the bond
};

typedef struct vdw_val
{
	Real sigma;	//  Sigma value
	Real epsilon;	//  Epsilon value
	Real sigma14;	//  Sigma value for 1-4 interactions
	Real epsilon14; //  Epsilon value for 1-4 interactions
} VdwValue;

//  IndexedVdwPair is used to form a binary search tree that is
//  indexed by vwd_type index.  This is the tree that will be
//  used to search during the actual simulation

typedef struct indexed_vdw_pair
{
   	Index ind1;		//  Index for first atom type
   	Index ind2;		//  Index for second atom type
   	Real A;			//  Parameter A for this pair
	Real A14;		//  Parameter A for 1-4 interactions
	Real B;			//  Parameter B for this pair
	Real B14;		//  Parameter B for 1-4 interactions
	struct indexed_vdw_pair *right;	 //  Right child
   	struct indexed_vdw_pair *left;	 //  Left child
} IndexedVdwPair;

//  Structures that are defined in Parameters.C
struct bond_params;
struct angle_params;
struct improper_params;
struct dihedral_params;
struct vdw_params;
struct vdw_pair_params;

class Parameters
{
private:
        char *atomTypeNames;                    //  Names of atom types
	Bool AllFilesRead;			//  Flag TRUE imples that all
						//  of the parameter files
						//  have been read in and
						//  the arrays have been
						//  created.
//****** BEGIN CHARMM/XPLOR type changes
        int paramType;                          //  Type (format) of parameter-file
//****** END CHARMM/XPLOR type changes
	struct bond_params *bondp;		//  Binary tree of bond params
	struct angle_params *anglep;		//  Binary tree of angle params
	struct improper_params *improperp;	//  Linked list of improper par.
	struct dihedral_params *dihedralp;      //  Linked list of dihedral par.
	struct vdw_params *vdwp;		//  Binary tree of vdw params
	struct vdw_pair_params *vdw_pairp;	//  Binary tree of vdw pairs
public:
	BondValue *bond_array;			//  Array of bond params
	AngleValue *angle_array;		//  Array of angle params
	DihedralValue *dihedral_array;		//  Array of dihedral params
	ImproperValue *improper_array;		//  Array of improper params
	VdwValue *vdw_array;			//  Array of vdw params
private:
	IndexedVdwPair *vdw_pair_tree;		//  Tree of vdw pair params
	int NumBondParams;			//  Number of bond parameters
	int NumAngleParams;			//  Number of angle parameters
	int NumDihedralParams;			//  Number of dihedral params
	int NumImproperParams;			//  Number of improper params
	int NumVdwParams;			//  Number of vdw parameters
        int NumVdwParamsAssigned;               //  Number actually assigned
	int NumVdwPairParams;			//  Number of vdw_pair params
	ResizeArray<char *> error_msgs;		//  Avoids repeating warnings

	int *maxDihedralMults;			//  Max multiplicity for
						//  dihedral bonds
	int *maxImproperMults;			//  Max multiplicity for
						//  improper bonds

	void add_bond_param(char *);		//  Add a bond parameter
	struct bond_params *add_to_bond_tree(struct bond_params * , 
				     struct bond_params *);

	void add_angle_param(char *);		//  Add an angle parameter
	struct angle_params *add_to_angle_tree(struct angle_params * , 
				     struct angle_params *);

	void add_dihedral_param(char *, FILE *); //  Add a dihedral parameter
	void add_to_dihedral_list(struct dihedral_params *);
	void add_to_charmm_dihedral_list(struct dihedral_params *);

	void add_improper_param(char *, FILE *); //  Add an improper parameter
	void add_to_improper_list(struct improper_params *);

	void add_vdw_param(char *);		//  Add a vdw parameter
	struct vdw_params *add_to_vdw_tree(struct vdw_params *, 
				     struct vdw_params *);

	void add_vdw_pair_param(char *);	//  Add a vdw pair parameter
	void add_to_vdw_pair_list(struct vdw_pair_params *);

	void add_hb_pair_param(char *);	//  Add a hydrogen bond pair parameter

	//  All of the traverse routines are used for debugging purposes
	//  to print out the appropriate list of parameters
	void traverse_vdw_pair_params(struct vdw_pair_params *);
	void traverse_vdw_params(struct vdw_params *);
	void traverse_dihedral_params(struct dihedral_params *);
	void traverse_improper_params(struct improper_params *);
	void traverse_angle_params(struct angle_params *);
	void traverse_bond_params(struct bond_params *);

	//  The index_* routines are used to index each of 
	//  the parameters and build the arrays that will be used
	//  for constant time access
	Index index_bonds(struct bond_params *, Index);
	Index index_angles(struct angle_params *, Index);
	Index index_vdw(struct vdw_params *, Index);
	void index_dihedrals();
	void index_impropers();
	
	void convert_vdw_pairs();
	IndexedVdwPair *add_to_indexed_vdw_pairs(IndexedVdwPair *, IndexedVdwPair *);
	
	int vdw_pair_to_arrays(int *, int *, Real *, Real *, Real *, Real *, 
			       int, IndexedVdwPair *);

	//  The free_* routines are used by the destructor to deallocate
	//  memory
	void free_bond_tree(struct bond_params *);
	void free_angle_tree(struct angle_params *);
	void free_dihedral_list(struct dihedral_params *);
	void free_improper_list(struct improper_params *);
	void free_vdw_tree(struct vdw_params *);
	void free_vdw_pair_tree(IndexedVdwPair *);
	void free_vdw_pair_list();

	/* does the actual initialization, once the variables have all
	   been given default values.  See Parameters() below */
        void read_parm(const GromacsTopFile *, Bool min);
public:
        //****** BEGIN CHARMM/XPLOR type changes
        //// added SimParameters to argument list
	Parameters(SimParameters *, StringList *f=NULL);
        //****** END CHARMM/XPLOR type changes
        
        Parameters(Ambertoppar *, BigReal);
        void read_parm(Ambertoppar *, BigReal);

        /* initializes this to hold the set of parameters in the
	   GROMACS topology file <gf>.  If the parameter <min> is on,
	   this assumes that we are going to do minimization and
	   therefore can't have atoms with zero VDW - it will add a
	   small repulsive term to these. */
        Parameters(const GromacsTopFile *gf, Bool min);
        
	~Parameters();				//  Destructor

        // return a string for the Nth atom type.  This can only be
        // called after all the param files have been read and the type
        // names have been indexed.  The Nth atom type refers to the same
        // index of the Nth vdw parameter (i.e. there are NumVdwParams names).
	char *atom_type_name(Index a) {
	  return (atomTypeNames + (a * (MAX_ATOMTYPE_CHARS + 1)));
        }

	//  Read a parameter file
	void read_parameter_file(char *);

        //****** BEGIN CHARMM/XPLOR type changes
	void read_charmm_parameter_file(char *);
        //****** END CHARMM/XPLOR type changes

	//  Signal the parameter object that all of
	//  the parameter files have been read in
	void done_reading_files();

	//  Signal the parameter object that the
	//  structure file has been read in
	void done_reading_structure();

	//  The assign_*_index routines are used to assign
	//  an index to atoms or bonds.  If an specific atom
	//  or bond type can't be found, then the program 
	//  terminates
	void assign_vdw_index(char *, Atom *);	//  Assign a vdw index to
						//  an atom
	void assign_bond_index(char *, char *, Bond *); 
						//  Assign a bond index
						//  to a bond
	void assign_angle_index(char *, char *, char *, Angle *);
						//  Assign an angle index
						//  to an angle
	void assign_dihedral_index(char *, char*, char*, char *, Dihedral *, int);
						//  Assign a dihedral index
						//  to a dihedral
	void assign_improper_index(char *, char*, char*, char *, Improper *, int);
						//  Assign an improper index
						//  to an improper

	//  send_parameters is used by the master process to
	//  communicate the paramters to all the other processors
	void send_Parameters(Communicate *);

	//  receive_parameters is used by all the child processes
	//  to receive the parameters from the master process
	void receive_Parameters(MIStream *);

	//  The get_*_params routines are the routines that really
	//  do all the work for this object.  Given an index, they
	//  access the parameters and return the relevant information
	void get_bond_params(Real *k, Real *x0, Index index)
	{
		*k = bond_array[index].k;
		*x0 = bond_array[index].x0;
	}

	void get_angle_params(Real *k, Real *theta0, Real *k_ub, Real *r_ub,
			      Index index)
	{
		*k = angle_array[index].k;
		*theta0 = angle_array[index].theta0;
		*k_ub = angle_array[index].k_ub;
		*r_ub = angle_array[index].r_ub;
	}

	int get_improper_multiplicity(Index index)
	{
		return(improper_array[index].multiplicity);
	}

	int get_dihedral_multiplicity(Index index)
	{
		return(dihedral_array[index].multiplicity);
	}

	void get_improper_params(Real *k, int *n, Real *delta, 
				 Index index, int mult)
	{
		if ( (mult<0) || (mult>MAX_MULTIPLICITY) )
		{
			NAMD_die("Bad mult index in Parameters::get_improper_params");
		}

		*k = improper_array[index].values[mult].k;
		*n = improper_array[index].values[mult].n;
		*delta = improper_array[index].values[mult].delta;
	}

	void get_dihedral_params(Real *k, int *n, Real *delta, 
				 Index index, int mult)
	{
		if ( (mult<0) || (mult>MAX_MULTIPLICITY) )
		{
			NAMD_die("Bad mult index in Parameters::get_dihedral_params");
		}

		*k = dihedral_array[index].values[mult].k;
		*n = dihedral_array[index].values[mult].n;
		*delta = dihedral_array[index].values[mult].delta;
	}

	void get_vdw_params(Real *sigma, Real *epsilon, Real *sigma14, 
			    Real *epsilon14, Index index)
	{
	    if ( vdw_array ) {
		*sigma = vdw_array[index].sigma;
		*epsilon = vdw_array[index].epsilon;
		*sigma14 = vdw_array[index].sigma14;
		*epsilon14 = vdw_array[index].epsilon14;
	    } else {
		*sigma = 0.;
		*epsilon = 0.;
		*sigma14 = 0.;
		*epsilon14 = 0.;
	    }
	}

	int get_vdw_pair_params(Index ind1, Index ind2, Real *, Real *, Real *, Real *);
						//  Find a vwd_pair parameter

        int get_num_vdw_params(void) { return NumVdwParamsAssigned; }

	//  The print_*_params are provided for debugging purposes
	void print_bond_params();		//  Print bonds params
	void print_angle_params();		//  Print angle params
	void print_dihedral_params();		//  Print dihedral params
	void print_improper_params();		//  Print improper params
	void print_vdw_params();		//  Print vdw params
	void print_vdw_pair_params();		//  Print vdw_pair params
	void print_param_summary();		//  Print a summary of params
};

#endif

