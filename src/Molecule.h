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

#include "molfile_plugin.h"
#include "ParallelIOMgr.h"

#include <vector>
using namespace std;

class SimParameters;
class Parameters;
class ConfigList;
class PDB;
class MIStream;
class MOStream;

class ExclElem;
class BondElem;
class AngleElem;
class DihedralElem;
class ImproperElem;
class TholeElem;  // Drude model
class AnisoElem;  // Drude model
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

//only used for compressing the molecule information
typedef struct seg_resid
{
    char segname[11];
    int resid;
}AtomSegResInfo;

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

friend class ExclElem;
friend class BondElem;
friend class AngleElem;
friend class DihedralElem;
friend class ImproperElem;
friend class TholeElem;  // Drude model
friend class AnisoElem;  // Drude model
friend class CrosstermElem;
friend class WorkDistrib;

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

  AtomSegResInfo *atomSegResids; //only used for generating compressed molecule info

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

  // DRUDE
  DrudeConst *drudeConsts;  // supplement Atom data (length of Atom array)
  Thole *tholes;            // Thole (screened Coulomb) interactions
  Aniso *anisos;            // anisotropic terms
  Lphost *lphosts;          // lone pair hosts
  int32 *lphostIndexes;     // index for each LP into lphosts array
  // DRUDE

  int32 *consIndexes; //  Constraint indexes for each atom
  ConstraintParams *consParams;

/* BEGIN gf */
  int32 **gridfrcIndexes;
  GridforceParams **gridfrcParams;
  GridforceGrid **gridfrcGrid;
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
  //A generally true assumption: all atoms are arranged in the order of clusters. In other words,
  //all atoms for a cluster must appear before/after any atoms in other clusters
  //The first atom in the cluster (which has the lowest atom id) stores the cluster size
  //other atoms in the cluster stores -1
  int32 *clusterSigs;
  int32 numClusters;
  //indicate whether the atom ids of a cluster are contiguous. If not, then
  //the general true assumption above mentioned will not hold. 
  int32 isClusterContiguous;
  #else
        int32 *cluster;   //  first atom of connected cluster
  #endif
  //In the memory optimized version: it will be NULL if the general
  //true assumption mentioned above holds. If not, its size is numClusters.
  //In the ordinary version: its size is numAtoms, and indicates the size
  //of connected cluster or 0.
  int32 *clusterSize;
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
  AtomSigID *eachAtomSig;
  ExclSigID *eachAtomExclSig;
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

  // DRUDE
  int32 **tholesByAtom;  // list of Thole correction terms owned by each atom
  int32 **anisosByAtom;  // list of anisotropic terms owned by each atom
  // DRUDE
#endif


  //occupancy and bfactor data from plugin-based IO implementation of loading structures
  float *occupancy;
  float *bfactor;

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

  // DRUDE: PSF reading
  void read_lphosts(FILE *);
        //  Read in lone pair hosts from Drude PSF
  void read_anisos(FILE *);
        //  Read in anisotropic terms from Drude PSF
  // DRUDE

  //pluginIO-based loading atoms' structure
  void plgLoadAtomBasics(molfile_atom_t *atomarray);
  void plgLoadBonds(int *from, int *to); //atom index is 1-based in the parameters
  void plgLoadAngles(int *plgAngles);
  void plgLoadDihedrals(int *plgDihedrals);
  void plgLoadImpropers(int *plgImpropers);
  void plgLoadCrossterms(int *plgCterms);


  void build12excl(void);
  void build13excl(void);
  void build14excl(int);

  // DRUDE: extend exclusions for Drude and LP
  void build_inherited_excl(int);
  // DRUDE
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
  // DRUDE
  int is_drude_psf;      // flag for reading Drude PSF
  int is_lonepairs_psf;  // flag for reading lone pairs from PSF
  // DRUDE

  // data for TIP4P
  Real r_om;
  Real r_ohc;

  // data for tail corrections
  BigReal tail_corr_ener;
  BigReal tail_corr_virial;



#ifdef MEM_OPT_VERSION
  AtomCstInfo *getAtoms() const { return atoms; }
  AtomNameIdx *getAtomNames() const { return atomNames; }
#else
  Atom *getAtoms () const { return atoms; }
  AtomNameInfo *getAtomNames() const { return atomNames; }
#endif

  AtomSegResInfo *getAtomSegResInfo() const { return atomSegResids; }
  
  // return number of fixed atoms, guarded by SimParameters
  int num_fixed_atoms() const {
    // local variables prefixed by s_
    int s_NumFixedAtoms = (simParams->fixedAtomsOn ? numFixedAtoms : 0);
    return s_NumFixedAtoms;  // value is "turned on" SimParameters
  }

  int num_fixed_groups() const {
    // local variables prefixed by s_
    int s_NumFixedAtoms = num_fixed_atoms();
    int s_NumFixedGroups = (s_NumFixedAtoms ? numFixedGroups : 0);
    return s_NumFixedGroups;  // value is "turned on" SimParameters
  }

  int num_group_deg_freedom() const {
    // local variables prefixed by s_
    int s_NumGroupDegFreedom = 3 * numHydrogenGroups;
    int s_NumFixedAtoms = num_fixed_atoms();
    int s_NumFixedGroups = num_fixed_groups();
    if (s_NumFixedGroups) s_NumGroupDegFreedom -= 3 * s_NumFixedGroups;
    if ( ! (s_NumFixedAtoms || numConstraints
          || simParams->comMove || simParams->langevinOn) ) {
      s_NumGroupDegFreedom -= 3;
    }
    return s_NumGroupDegFreedom;
  }

  int num_deg_freedom(int isInitialReport = 0) const {
    // local variables prefixed by s_
    int s_NumDegFreedom = 3 * numAtoms;
    int s_NumFixedAtoms = num_fixed_atoms();
    if (s_NumFixedAtoms) s_NumDegFreedom -= 3 * s_NumFixedAtoms;
    if (numLonepairs) s_NumDegFreedom -= 3 * numLonepairs;
    if ( ! (s_NumFixedAtoms || numConstraints
          || simParams->comMove || simParams->langevinOn) ) {
      s_NumDegFreedom -= 3;
    }
    if ( ! isInitialReport && simParams->pairInteractionOn) {
      //
      // DJH: a kludge?  We want to initially report system degrees of freedom
      //
      // this doesn't attempt to deal with fixed atoms or constraints
      s_NumDegFreedom = 3 * numFepInitial;
    }
    int s_NumFixedRigidBonds = 
      (simParams->fixedAtomsOn ? numFixedRigidBonds : 0);
    if (simParams->watmodel == WAT_TIP4) {
      // numLonepairs is subtracted here because all lonepairs have a rigid bond
      // to oxygen, but all of the LP degrees of freedom are dealt with above
      s_NumDegFreedom -= (numRigidBonds - s_NumFixedRigidBonds - numLonepairs);
    }
    else {
      // Non-polarized systems don't have LPs.
      // For Drude model, bonds that attach LPs are not counted as rigid;
      // LPs have already been subtracted from degrees of freedom.
      s_NumDegFreedom -= (numRigidBonds - s_NumFixedRigidBonds);
    }
    return s_NumDegFreedom;
  }

  int numAtoms;   //  Number of atoms                   

  int numRealBonds;   //  Number of bonds for exclusion determination
  int numBonds;   //  Number of bonds calculated, including extras
  int numAngles;    //  Number of angles
  int numDihedrals; //  Number of dihedrals
  int suspiciousAlchBonds;    //  angles dropped due to crossing FEP partitions
  int alchDroppedAngles;    //  angles dropped due to crossing FEP partitions
  int alchDroppedDihedrals; //  dihedrals dropped due to crossing FEP partitions
  int alchDroppedImpropers; //  impropers dropped due to crossing FEP partitions
  int numImpropers; //  Number of impropers
  int numCrossterms; //  Number of cross-terms
  int numDonors;          //  Number of hydrogen bond donors
  int numAcceptors; //  Number of hydrogen bond acceptors
  int numExclusions;  //  Number of exclusions
  int numTotalExclusions; //  Real Total Number of Exclusions // hack

  // DRUDE
  int numLonepairs; // Number of lone pairs
  int numDrudeAtoms;  // Number of Drude particles
  int numTholes;  // Number of Thole terms
  int numAnisos;  // Number of anisotropic terms
  int numLphosts;  // Number of lone pair hosts
  // DRUDE
  
  int numConstraints; //  Number of atoms constrained
/* BEGIN gf */
  int numGridforceGrids;//  Number of gridforce grids
  int *numGridforces;	//  Number of atoms in gridforce file (array, one per grid)
/* END gf */
  int numMovDrag;         //  Number of atoms moving-dragged
  int numRotDrag;         //  Number of atoms rotating-dragged
  int numConsTorque;  //  Number of atoms "constant"-torqued
  int numFixedAtoms;  //  Number of fixed atoms
  int numStirredAtoms;  //  Number of stirred atoms
  int numExPressureAtoms; //  Number of atoms excluded from pressure
  int numHydrogenGroups;  //  Number of hydrogen groups
  int maxHydrogenGroupSize;  //  Max atoms per hydrogen group
  int numMigrationGroups;  //  Number of migration groups
  int maxMigrationGroupSize;  //  Max atoms per migration group
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

  // DRUDE
  int numCalcTholes;  // Number of Thole correction terms requiring calculation
  int numCalcAnisos;  // Number of anisotropic terms requiring calculation
  // DRUDE

  //  Number of dihedrals with multiple periodicity
  int numMultipleDihedrals; 
  //  Number of impropers with multiple periodicity
  int numMultipleImpropers; 
  // indexes of "atoms" sorted by hydrogen groups
  HydrogenGroup hydrogenGroup;

  Molecule(SimParameters *, Parameters *param);
  Molecule(SimParameters *, Parameters *param, char *filename, ConfigList *cfgList=NULL);  
  Molecule(SimParameters *simParams, Parameters *param, molfile_plugin_t *pIOHdl, void *pIOFileHdl, int natoms);
  
  Molecule(SimParameters *, Parameters *, Ambertoppar *);
  void read_parm(Ambertoppar *);

  Molecule(SimParameters *, Parameters *, const GromacsTopFile *);

  ~Molecule();    //  Destructor

  void read_psf_file(char *, Parameters *);
        //  Read in a .psf file given
        //  the filename and the parameter
        //  object to use  

  void send_Molecule(MOStream *);
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

  void build_langevin_params(BigReal coupling, BigReal drudeCoupling,
      Bool doHydrogen);
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
        void delete_alch_bonded(void);
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
  //the way to get the cluster size if the atom ids of the cluster are
  //contiguous. The input parameter is the atom id that leads the cluster.
  int get_cluster_size_con(int aid) const { return clusterSigs[aid]; }  
  //the way to get the cluster size if the atoms ids of the cluster are
  //not contiguous. The input parameter is the cluster index.
  int get_cluster_size_uncon(int cIdx) const { return clusterSize[cIdx]; }
  int get_cluster_idx(int aid) const { return clusterSigs[aid]; }
  int get_num_clusters() const { return numClusters; }
  int get_cluster_contiguity() { return isClusterContiguous; }
  #else
  int get_cluster(int anum) const { return cluster[anum]; }
  int get_clusterSize(int anum) const { return clusterSize[anum]; }
  #endif

  const float *getOccupancyData() { return (const float *)occupancy; }
  void setOccupancyData(molfile_atom_t *atomarray);
  void freeOccupancyData() { delete [] occupancy; occupancy=NULL; }

  const float *getBFactorData() { return (const float *)bfactor; }
  void setBFactorData(molfile_atom_t *atomarray);
  void freeBFactorData() { delete [] bfactor; bfactor=NULL; }

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

  // DRUDE: retrieve lphost structure
  Lphost *get_lphost(int atomid) const {
    // don't call unless simParams->drudeOn == TRUE
    // otherwise lphostIndexes array doesn't exist!
    int index = lphostIndexes[atomid];
    return (index != -1 ? &(lphosts[index]) : NULL);
  }
  // DRUDE

  #ifndef MEM_OPT_VERSION
  Bond *getAllBonds() const {return bonds;}
  Angle *getAllAngles() const {return angles;}
  Improper *getAllImpropers() const {return impropers;}
  Dihedral *getAllDihedrals() const {return dihedrals;}
  Crossterm *getAllCrossterms() const {return crossterms;}
  #endif

  // DRUDE: retrieve entire lphosts array
  Lphost *getAllLphosts() const { return lphosts; }
  // DRUDE

  //  Retrieve a hydrogen bond donor structure
  Bond *get_donor(int dnum) const {return (&(donors[dnum]));}  

  //  Retrieve a hydrogen bond acceptor structure
  Bond *get_acceptor(int dnum) const {return (&(acceptors[dnum]));} 

  Bond *getAllDonors() const {return donors;}
  Bond *getAllAcceptors() const {return acceptors;}

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
  Bool is_atom_gridforced(int atomnum, int gridnum) const
  {
      if (numGridforceGrids)
      {
	  return(gridfrcIndexes[gridnum][atomnum] != -1);
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
  void get_gridfrc_params(Real &k, Charge &q, int atomnum, int gridnum) const
  {
      k = gridfrcParams[gridnum][gridfrcIndexes[gridnum][atomnum]].k;
      q = gridfrcParams[gridnum][gridfrcIndexes[gridnum][atomnum]].q;
  }
  
  GridforceGrid* get_gridfrc_grid(int gridnum) const
  {
      return gridfrcGrid[gridnum];
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

public:  
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Osman Sarood
// Parallel Input
// comment out the read_compressed_psf_file from above and make it public. Declare the following new methods as well.
  int isOccupancyValid, isBFactorValid;

  void read_hdr_info(char *fname, Parameters *params, ConfigList *cfgList=0);
  int getNumCalcExclusions(){return numCalcExclusions;}
  void setNumCalcExclusions(int x){numCalcExclusions= x;}
  void read_compressed_psf_file(char *, Parameters *, ConfigList *cfgList);
  void read_basic_info(char *, Parameters *, ConfigList *cfgList);
  void getCountsToMaster();
#ifdef MEM_OPT_VERSION
  Index getEachAtomMass(int i){return eachAtomMass[i];}
  Index getEachAtomCharge(int i){return eachAtomCharge[i];}

  ExclSigID getAtomExclSigId(int aid) const {
      return eachAtomExclSig[aid];
  }

  void read_compressed_psf_file_parallelIO(char *fname, Parameters *params);
  Real *getAtomMassPool(){return atomMassPool;}
  Real *getAtomChargePool(){return atomChargePool;}
  AtomCstInfo *getAtoms(){return atoms;}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

#endif

};

#endif

