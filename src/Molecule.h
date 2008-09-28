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

#include "parm.h"

#include "common.h"
#include "NamdTypes.h"
#include "structures.h"
#include "ConfigList.h"
#include "Vector.h"
#include "UniqueSet.h"
#include "Hydrogen.h"
#include "GromacsTopFile.h"
/* BEGIN gf */
#include "GridForceGrid.h"
#include "Tensor.h"
/* END gf */

#include <vector>
using namespace std;

class SimParameters;
class Parameters;
class ConfigList;
class PDB;
class MIStream;
class MOStream;

class BondElem;
class AngleElem;
class DihedralElem;
class ImproperElem;
class CrosstermElem;
class ResidueLookupElem;
template<class Type> class ObjectArena;

class ExclusionCheck {
public:
  int32 min,max;
  char *flags;

  ExclusionCheck(){
      min=0;
      max=-1;
      flags = NULL;
  }
  ExclusionCheck(const ExclusionCheck& chk){
      min = chk.min;
      max = chk.max;
      if(max>min){ 
	flags = new char[max-min+1];
	memcpy(flags, chk.flags, sizeof(char)*(max-min+1));
      }
  }
  ExclusionCheck &operator=(const ExclusionCheck& chk){
    min = chk.min;
    max = chk.max;
    if(flags) delete [] flags;
    flags = NULL;
    if(max>min){
	flags = new char[max-min+1];
	memcpy(flags, chk.flags, sizeof(char)*(max-min+1));
    }
    return *this;
  }
  ~ExclusionCheck(){
      if(flags) delete [] flags;
  }
};
#define EXCHCK_FULL 1
#define EXCHCK_MOD 2

// List maintaining the global atom indicies sorted by helix groups.
class Molecule
{
typedef struct constraint_params
{
   Real k;    //  Force constant
   Vector refPos; //  Reference position for restraint
} ConstraintParams;



/* BEGIN gf */
typedef struct gridfrc_params
{
    Real k;	// force multiplier
    Charge q;	// charge
} GridforceParams;
/* END gf */


typedef struct stir_params
{
  Real startTheta;
  Vector refPos;  //  Reference position for stirring
} StirParams;

typedef struct movdrag_params
{
  Vector v;            //  Linear velocity (A/step)
} MovDragParams;


typedef struct rotdrag_params
{
   Real v;              //  Angular velocity (deg/step)
   Vector a;            //  Rotation axis
   Vector p;            //  Rotation pivot point
} RotDragParams;

typedef struct constorque_params
{
   Real v;              //  "Torque" value (Kcal/(mol*A^2))
   Vector a;            //  Rotation axis
   Vector p;            //  Rotation pivot point
} ConsTorqueParams;

friend class BondElem;
friend class AngleElem;
friend class DihedralElem;
friend class ImproperElem;
friend class CrosstermElem;

private:
  void initialize(SimParameters *, Parameters *param);
            // Sets most data to zero

  #ifdef MEM_OPT_VERSION 
  //Indexing to constant pools to save space
  AtomCstInfo *atoms;
  Index *eachAtomMass; //after initialization, this could be freed (possibly)
  Index *eachAtomCharge; //after initialization, this could be freed (possibly)
  AtomNameIdx *atomNames;
  ObjectArena<char> *nameArena; //the space for names to be allocated  
  #else
  Atom *atoms;    //  Array of atom structures
  ObjectArena<char> *nameArena;
  AtomNameInfo *atomNames;//  Array of atom name info.  Only maintained on node 0 for VMD interface
  #endif

  ResidueLookupElem *resLookup; // find residues by name

  #ifndef MEM_OPT_VERSION
  //replaced by atom signatures
  Bond *bonds;    //  Array of bond structures
  Angle *angles;    //  Array of angle structures
  Dihedral *dihedrals;  //  Array of dihedral structures
  Improper *impropers;  //  Array of improper structures                          
  Crossterm *crossterms;  //  Array of cross-term structures
  #endif
  
  Bond *donors;         //  Array of hydrogen bond donor structures
  Bond *acceptors;  //  Array of hydrogen bond acceptor

  #ifndef MEM_OPT_VERSION
  //These will be replaced by exclusion signatures
  Exclusion *exclusions;  //  Array of exclusion structures
  UniqueSet<Exclusion> exclusionSet;  //  Used for building
  #endif

  int32 *consIndexes; //  Constraint indexes for each atom
  ConstraintParams *consParams;

/* BEGIN gf */
  int32 *gridfrcIndexes;
  GridforceParams *gridfrcParams;
  GridforceGrid *gridfrcGrid;
/* END gf */

        //  Parameters for each atom constrained
  int32 *stirIndexes; //Stirring indexes for each atoms
  StirParams *stirParams;
                          // Paramters for each atoms stirred
  int32 *movDragIndexes;  //  Moving drag indexes for each atom
  MovDragParams *movDragParams;
                                //  Parameters for each atom moving-dragged
  int32 *rotDragIndexes;  //  Rotating drag indexes for each atom
  RotDragParams *rotDragParams;
                                //  Parameters for each atom rotation-dragged

  Real *langevinParams;   //  b values for langevin dynamics
  int32 *fixedAtomFlags;  //  1 for fixed, -1 for fixed group, else 0
  int32 *exPressureAtomFlags; // 1 for excluded, -1 for excluded group.

  #ifdef MEM_OPT_VERSION  
  //Assumption: all atoms are arranged in the order of clusters. In other words,
  //all atoms for a cluster must appear before/after any atoms in other clusters
  //The first atom in the cluster (which has the lowest atom id) stores the cluster size
  //other atoms in the cluster stores -1
  int32 *clusterSigs;
  #else
        int32 *cluster;   //  first atom of connected cluster
        int32 *clusterSize; //  size of connected cluster or 0
  #endif
                            // 
        Real *rigidBondLengths;  //  if H, length to parent or 0. or
        //  if not H, length between children or 0.
//fepb
        unsigned char *fepAtomFlags; 
//fepe

#ifndef MEM_OPT_VERSION
  ObjectArena<int32> *tmpArena;
  int32 **bondsWithAtom;  //  List of bonds involving each atom
  ObjectArena<int32> *arena;
#endif

#ifdef MEM_OPT_VERSION
  Index *eachAtomSig;
  Index *eachAtomExclSig;
#else
//function is replaced by atom signatures
  int32 **bondsByAtom;  //  List of bonds owned by each atom
  int32 **anglesByAtom;     //  List of angles owned by each atom
  int32 **dihedralsByAtom;  //  List of dihedrals owned by each atom
  int32 **impropersByAtom;  //  List of impropers owned by each atom
  int32 **crosstermsByAtom;  //  List of crossterms owned by each atom
    
  int32 **exclusionsByAtom; //  List of exclusions owned by each atom
  int32 **fullExclusionsByAtom; //  List of atoms excluded for each atom
  int32 **modExclusionsByAtom; //  List of atoms modified for each atom
  ObjectArena<char> *exclArena;
  ExclusionCheck *all_exclusions;
#endif
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
  void read_crossterms(FILE *, Parameters *);
        //  Read in cross-term info from .psf
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
  #ifdef MEM_OPT_VERSION
  void stripFepFixedExcl(void);
  #else
  void stripFepExcl(void);
  #endif

  void build_exclusions();

  // analyze the atoms, and determine which are oxygen, hb donors, etc.
  // this is called after a molecule is sent our (or received in)
  void build_atom_status(void);

  // added during the trasition from 1x to 2
  SimParameters *simParams;
  Parameters *params;

  void read_parm(const GromacsTopFile *);  

public:
  // data for TIP4P
  Real r_om;
  Real r_ohc;



#ifndef MEM_OPT_VERSION
  Atom *getAtoms () const { return atoms; }
#endif

  int numAtoms;   //  Number of atoms                   

  int numRealBonds;   //  Number of bonds for exclusion determination
  int numBonds;   //  Number of bonds calculated, including extras
  int numAngles;    //  Number of angles
  int numDihedrals; //  Number of dihedrals
  int numImpropers; //  Number of impropers
  int numCrossterms; //  Number of cross-terms
  int numDonors;          //  Number of hydrogen bond donors
  int numAcceptors; //  Number of hydrogen bond acceptors
  int numExclusions;  //  Number of exclusions
  int numTotalExclusions; //  Real Total Number of Exclusions // hack
  int numLP; // Number of lone pairs
  
  int numConstraints; //  Number of atoms constrained
/* BEGIN gf */
  int numGridforces;  //  Number of atoms in gridforce file
/* END gf */
  int numMovDrag;         //  Number of atoms moving-dragged
  int numRotDrag;         //  Number of atoms rotating-dragged
  int numConsTorque;  //  Number of atoms "constant"-torqued
  int numFixedAtoms;  //  Number of fixed atoms
  int numStirredAtoms;  //  Number of stirred atoms
  int numExPressureAtoms; //  Number of atoms excluded from pressure
  int numHydrogenGroups;  //  Number of hydrogen groups
  int numFixedGroups; //  Number of totally fixed hydrogen groups
  int numRigidBonds;  //  Number of rigid bonds
  int numFixedRigidBonds; //  Number of rigid bonds between fixed atoms
//fepb
        int numFepInitial;  // no. of fep atoms with initial flag
        int numFepFinal;  // no. of fep atoms with final flag
//fepe

  int numConsForce; //  Number of atoms that have constant force applied
  int32 *consForceIndexes;//  Constant force indexes for each atom
  Vector *consForce;  //  Constant force array

  int32 *consTorqueIndexes; //  "Constant" torque indexes for each atom
  ConsTorqueParams *consTorqueParams;
                                //  Parameters for each atom "constant"-torqued

  // The following are needed for error checking because we
  // eliminate bonds, etc. which involve only fixed atoms
  int numCalcBonds; //  Number of bonds requiring calculation
  int numCalcAngles;  //  Number of angles requiring calculation
  int numCalcDihedrals; //  Number of dihedrals requiring calculation
  int numCalcImpropers; //  Number of impropers requiring calculation
  int numCalcCrossterms; //  Number of cross-terms requiring calculation
  int numCalcExclusions;  //  Number of exclusions requiring calculation

  //  Number of dihedrals with multiple periodicity
  int numMultipleDihedrals; 
  //  Number of impropers with multiple periodicity
  int numMultipleImpropers; 
  // indexes of "atoms" sorted by hydrogen groups
  HydrogenGroup hydrogenGroup;
  int waterIndex;

  Molecule(SimParameters *, Parameters *param);
  Molecule(SimParameters *, Parameters *param, char *filename, ConfigList *cfgList=NULL);  
  
  Molecule(SimParameters *, Parameters *, Ambertoppar *);
  void read_parm(Ambertoppar *);

  Molecule(SimParameters *, Parameters *, const GromacsTopFile *);

  ~Molecule();    //  Destructor

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

/* BEGIN gf */
  void build_gridforce_params(StringList *, StringList *, StringList *, StringList *, PDB *, char *);
        //  Build the set of gridForce-style force pars
/* END gf */

  void build_movdrag_params(StringList *, StringList *, StringList *, 
          PDB *, char *);
        //  Build the set of moving drag pars

  void build_rotdrag_params(StringList *, StringList *, StringList *,
          StringList *, StringList *, StringList *,
          PDB *, char *);
        //  Build the set of rotating drag pars

  void build_constorque_params(StringList *, StringList *, StringList *,
             StringList *, StringList *, StringList *,
             PDB *, char *);
        //  Build the set of "constant" torque pars


  void build_constant_forces(char *);
        //  Build the set of constant forces

  void build_langevin_params(BigReal coupling, Bool doHydrogen);
  void build_langevin_params(StringList *, StringList *, PDB *, char *);
        //  Build the set of langevin dynamics parameters

  void build_fixed_atoms(StringList *, StringList *, PDB *, char *);
        //  Determine which atoms are fixed (if any)

  void build_stirred_atoms(StringList *, StringList *, PDB *, char *);
        //  Determine which atoms are stirred (if any)

  void build_extra_bonds(Parameters *parameters, StringList *file);

//fepb
        void build_fep_flags(StringList *, StringList *, PDB *, char *);
                               // selection of the mutant atoms
//fepe

  void build_exPressure_atoms(StringList *, StringList *, PDB *, char *);
        //  Determine which atoms are excluded from
                                //  pressure (if any)

  void reloadCharges(float charge[], int n);

        Bool is_lp(int);     // return true if atom is a lone pair
        Bool is_drude(int);     // return true if atom is a Drude particle
        Bool is_hydrogen(int);     // return true if atom is hydrogen
        Bool is_oxygen(int);       // return true if atom is oxygen
  Bool is_hydrogenGroupParent(int); // return true if atom is group parent
  Bool is_water(int);        // return true if atom is part of water 
  int  get_groupSize(int);     // return # atoms in (hydrogen) group
        int get_mother_atom(int);  // return mother atom of a hydrogen

  #ifdef MEM_OPT_VERSION
  int get_cluster_size(int idx) const { return clusterSigs[idx]; }  
  #else
  int get_cluster(int anum) const { return cluster[anum]; }
  int get_clusterSize(int anum) const { return clusterSize[anum]; }
  #endif

  #ifdef CHARMIZE_NAMD
  Atom *getAllAtoms() {
      return atoms;
  }
  #endif

  //  Get the mass of an atom
  Real atommass(int anum) const
  {
    #ifdef MEM_OPT_VERSION
    return atomMassPool[eachAtomMass[anum]];
    #else
    return(atoms[anum].mass);
    #endif
  }

  //  Get the charge of an atom
  Real atomcharge(int anum) const
  {
    #ifdef MEM_OPT_VERSION
    return atomChargePool[eachAtomCharge[anum]];
    #else
    return(atoms[anum].charge);
    #endif
  }
  
  //  Get the vdw type of an atom
  Index atomvdwtype(int anum) const
  {      
      return(atoms[anum].vdw_type);
  }

  #ifndef MEM_OPT_VERSION
  //  Retrieve a bond structure
  Bond *get_bond(int bnum) const {return (&(bonds[bnum]));}

  //  Retrieve an angle structure
  Angle *get_angle(int anum) const {return (&(angles[anum]));}

  //  Retrieve an improper strutcure
  Improper *get_improper(int inum) const {return (&(impropers[inum]));}

  //  Retrieve a dihedral structure
  Dihedral *get_dihedral(int dnum) const {return (&(dihedrals[dnum]));}

  //  Retrieve a cross-term strutcure
  Crossterm *get_crossterm(int inum) const {return (&(crossterms[inum]));}
  #endif

  //  Retrieve a hydrogen bond donor structure
  Bond *get_donor(int dnum) const {return (&(donors[dnum]));}

  //  Retrieve a hydrogen bond acceptor structure
  Bond *get_acceptor(int dnum) const {return (&(acceptors[dnum]));}

  //  Retrieve an exclusion structure
  #ifndef MEM_OPT_VERSION
  Exclusion *get_exclusion(int ex) const {return (&(exclusions[ex]));}
  #endif

  //  Retrieve an atom type
  const char *get_atomtype(int anum) const
  {
    if (atomNames == NULL)
    {
      NAMD_die("Tried to find atom type on node other than node 0");
    }

    #ifdef MEM_OPT_VERSION    
    return atomTypePool[atomNames[anum].atomtypeIdx];
    #else
    return(atomNames[anum].atomtype);
    #endif
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

  #ifndef MEM_OPT_VERSION
  int32 *get_bonds_for_atom(int anum)
      { return bondsByAtom[anum]; } 
  int32 *get_angles_for_atom(int anum) 
      { return anglesByAtom[anum]; }
  int32 *get_dihedrals_for_atom(int anum) 
      { return dihedralsByAtom[anum]; }
  int32 *get_impropers_for_atom(int anum) 
      { return impropersByAtom[anum]; }  
  int32 *get_crossterms_for_atom(int anum) 
      { return crosstermsByAtom[anum]; }  
  int32 *get_exclusions_for_atom(int anum)
      { return exclusionsByAtom[anum]; }
  const int32 *get_full_exclusions_for_atom(int anum) const
      { return fullExclusionsByAtom[anum]; }
  const int32 *get_mod_exclusions_for_atom(int anum) const
      { return modExclusionsByAtom[anum]; }
  #endif
  
  //  Check for exclusions, either explicit or bonded.
        //  Returns 1 for full, 2 for 1-4 exclusions.
  #ifdef MEM_OPT_VERSION
  int checkExclByIdx(int idx1, int atom1, int atom2) const;
  const ExclusionCheck *get_excl_check_for_idx(int idx) const{      
      return &exclChkSigPool[idx];
  }
  #else
  int checkexcl(int atom1, int atom2) const;

  const ExclusionCheck *get_excl_check_for_atom(int anum) const{      
      return &all_exclusions[anum];             
  }
  #endif

/* BEGIN gf */
  // Return true or false based on whether or not the atom
  // is subject to grid force
  Bool is_atom_gridforced(int atomnum) const
  {
      if (numGridforces)
      {
	  return(gridfrcIndexes[atomnum] != -1);
      }
      else
      {
	  return(FALSE);
      }
  }
/* END gf */

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

  //  Return true or false based on whether the specified atom
  //  is moving-dragged or not.
  Bool is_atom_movdragged(int atomnum) const
  {
    if (numMovDrag)
    {
      //  Check the index to see if it is constrained
      return(movDragIndexes[atomnum] != -1);
    }
    else
    {
      //  No constraints at all, so just return FALSE
      return(FALSE);
    }
  }

  //  Return true or false based on whether the specified atom
  //  is rotating-dragged or not.
  Bool is_atom_rotdragged(int atomnum) const
  {
    if (numRotDrag)
    {
      //  Check the index to see if it is constrained
      return(rotDragIndexes[atomnum] != -1);
    }
    else
    {
      //  No constraints at all, so just return FALSE
      return(FALSE);
    }
  }

  //  Return true or false based on whether the specified atom
  //  is "constant"-torqued or not.
  Bool is_atom_constorqued(int atomnum) const
  {
    if (numConsTorque)
    {
      //  Check the index to see if it is constrained
      return(consTorqueIndexes[atomnum] != -1);
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

/* BEGIN gf */
  void get_gridfrc_params(Real &k, Charge &q, int atomnum) const
  {
      k = gridfrcParams[gridfrcIndexes[atomnum]].k;
      q = gridfrcParams[gridfrcIndexes[atomnum]].q;
  }
  
  const GridforceGrid* get_gridfrc_grid(void) const
  {
      return gridfrcGrid;
  }
/* END gf */

  Real langevin_param(int atomnum) const
  {
    return(langevinParams[atomnum]);
  }

  //  Get the stirring constraints for a specific atom
  void get_stir_refPos(Vector &refPos, int atomnum) const
  {
    refPos = stirParams[stirIndexes[atomnum]].refPos;
  }


  void put_stir_startTheta(Real theta, int atomnum) const
  {
    stirParams[stirIndexes[atomnum]].startTheta = theta;
  }


  Real get_stir_startTheta(int atomnum) const
  {
    return stirParams[stirIndexes[atomnum]].startTheta;
  }
 

  //  Get the moving drag factor for a specific atom
  void get_movdrag_params(Vector &v, int atomnum) const
  {
    v = movDragParams[movDragIndexes[atomnum]].v;
  }

  //  Get the rotating drag pars for a specific atom
  void get_rotdrag_params(BigReal &v, Vector &a, Vector &p, 
        int atomnum) const
  {
    v = rotDragParams[rotDragIndexes[atomnum]].v;
    a = rotDragParams[rotDragIndexes[atomnum]].a;
    p = rotDragParams[rotDragIndexes[atomnum]].p;
  }

  //  Get the "constant" torque pars for a specific atom
  void get_constorque_params(BigReal &v, Vector &a, Vector &p, 
        int atomnum) const
  {
    v = consTorqueParams[consTorqueIndexes[atomnum]].v;
    a = consTorqueParams[consTorqueIndexes[atomnum]].a;
    p = consTorqueParams[consTorqueIndexes[atomnum]].p;
  }

//fepb
        unsigned char get_fep_type(int anum) const
        {
                return(fepAtomFlags[anum]);
        }
//fepe

  Bool is_atom_fixed(int atomnum) const
  {
    return (numFixedAtoms && fixedAtomFlags[atomnum]);
  }

        
  Bool is_atom_stirred(int atomnum) const
  {
    if (numStirredAtoms)
    {
      //  Check the index to see if it is constrained
      return(stirIndexes[atomnum] != -1);
    }
    else
    {
      //  No constraints at all, so just return FALSE
      return(FALSE);
    }
  }
  

  Bool is_group_fixed(int atomnum) const
  {
    return (numFixedAtoms && (fixedAtomFlags[atomnum] == -1));
  }
  Bool is_atom_exPressure(int atomnum) const
  {
    return (numExPressureAtoms && exPressureAtomFlags[atomnum]);
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

private:
  //Read in a compressed .psf file (whose format is not standard) for
  //the sake of very very large simulations (say, 100M atoms system)
  void read_compressed_psf_file(char *, Parameters *);   

#ifdef MEM_OPT_VERSION
public:  
  int atomSigPoolSize;
  AtomSignature *atomSigPool;

  /* All the following are temporary variables for reading the compressed psf file */
  //declarations for atoms' constant information  
  int segNamePoolSize; //Its value is usually less than 5
  char **segNamePool; //This seems not to be important, but it only occupied very little space.

  int resNamePoolSize;
  char **resNamePool;

  int atomNamePoolSize;
  char **atomNamePool;

  int atomTypePoolSize;
  char **atomTypePool;

  int chargePoolSize;
  Real *atomChargePool;

  int massPoolSize;
  Real *atomMassPool;

  AtomSigID getAtomSigId(int aid) {
      return eachAtomSig[aid]; 
  }
  ExclSigID getAtomExclSigId(int aid) const {
      return eachAtomExclSig[aid];
  }

  //Indicates the size of both exclSigPool and exclChkSigPool
  int exclSigPoolSize;
  //this will be deleted after build_lists_by_atom
  ExclusionSignature *exclSigPool;
  //This is the final data structure we want to store  
  ExclusionCheck *exclChkSigPool;

  void addNewExclSigPool(const vector<ExclusionSignature>&);

  void build_excl_check_signatures();

  void delEachAtomSigs(){    
      //for NAMD-smp version, only one Molecule object is held
      //on each node, therefore, only one deletion operation should
      //be taken on a node, otherwise, there possibly would be some
      //wierd memory problems. The same reason applies to other deletion
      //operations inside the Molecule object.   
      if(CmiMyRank()) return;

      delete [] eachAtomSig;
      delete [] eachAtomExclSig;
      eachAtomSig = NULL;
      eachAtomExclSig = NULL;
  }

  void delChargeSpace(){
      if(CmiMyRank()) return;

      delete [] atomChargePool;
      delete [] eachAtomCharge;
      atomChargePool = NULL;
      eachAtomCharge = NULL;
  }
  
  void delMassSpace(){
      if(CmiMyRank()) return;

      delete [] atomMassPool;
      delete [] eachAtomMass;
      atomMassPool = NULL;
      eachAtomMass = NULL;
  }
  
  void delClusterSigs() {
      if(CmiMyRank()) return;      

      delete [] clusterSigs;
      clusterSigs = NULL;
  }

  void delOtherEachAtomStructs();


private:
  Index insert_new_mass(Real newMass);
  //void markClusterIdx(int curClusterIdx, int startAtomID);  

  //This function checks whether the exclusion (atom1 and atom2) will
  //be stripped by the original stripFepExcl function and by the fixed atoms.
  //When stripping exclusions in the stripHGroupExcl function (where we actually
  //not performing stripping but only deducting numCalcExclusions), we use
  //this function to avoid the repeated deduction. Only exclusions that are only
  //stripped due to hydrogen group is counted
  int exclStrippedByFepOrFixedAtoms(int atom1, int atom2);
#endif

};

#endif

