/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

#include "ConfigList.h"
#include "common.h"
#include "Communicate.h"
#include "Vector.h"
#include "Lattice.h"

class ParseOptions;

//  The class SimParameters is really just a glorified structure used to
//  maintain the global simulation parameters.  The only functions
//  associated with the class are used to get the parameters from the
//  ConfigList object, to send that Parameters from the master node 
//  to the other nodes, and to receive the Parameters on the other nodes.


//  The following definitions are used to distinguish between possible
//  bonded exclusion settings
typedef int  ExclusionSettings;

#define NONE		0
#define	ONETWO  	1
#define	ONETHREE	2
#define	ONEFOUR		3
#define	SCALED14 	4

//  The following definitions are used to distinguish between multiple
//  timestep integration schemes
typedef int  MTSChoices;

#define NAIVE		0
#define VERLETI		1

//  The following definitions are used to distinuish between multiple
//  long-short range force splittings
#define SHARP		0
#define XPLOR		1
#define C1		2

//  The following definitions are used to distinguish among load
//  balancing strategies
#define LDBSTRAT_NONE    0
#define LDBSTRAT_REFINEONLY 1
#define LDBSTRAT_ALG7 2
#define LDBSTRAT_OTHER  99

// The following definitions are used to distinguish between patch-splitting
// strategies
#define SPLIT_PATCH_POSITION	0	// atom position determines patch
#define SPLIT_PATCH_HYDROGEN	1	// hydrogen groups are not broken up

// The following definitions are used to distinguish the range of rigid
// bond calculations: none, all bonds to hydrogen, or only water
#define RIGID_NONE    0
#define RIGID_ALL     1
#define RIGID_WATER   2

class SimParameters
{
private:
public:

//  MAKE SURE THAT THIS CLASS CAN BE BIT COPIED OR YOU WILL HAVE TO
//  ADD SPECIAL CODE TO send_SimParameters() and receive_SimParameters()

  char dummy;
	BigReal dt;	   		//  Timestep size
	int N;		   		//  Number of steps to be performed
	int stepsPerCycle;		//  Number of timesteps per cycle

	zVector cellBasisVector1;	//  Basis vector for periodic cell
	zVector cellBasisVector2;	//  Basis vector for periodic cell
	zVector cellBasisVector3;	//  Basis vector for periodic cell
	zVector cellOrigin;		//  Fixed center of periodic cell
	Lattice lattice;		//  All data for periodic cell
	
	int nonbondedFrequency;		//  Number of timesteps between
					//  nonbonded evaluation
	int fullElectFrequency;		//  Number of timesteps between
					//  full electrostatic evaluation
        BigReal fmaTheta;	        //  DPMTA theta value
	int ldbStrategy;                //  What type of load balancing
	int ldbPeriod;                  //  How often to do load balancing
	int firstLdbStep;		//  What step to do the first 
                                        //  load-balance on.
	BigReal initialTemp;   		//  Initial temperature for the 
					//  simulation
	Bool comMove;     		//  Should the center of mass be 
					//  able to move
	Bool wrapWater;			//  Wrap water around on output
	BigReal dielectric;   		//  Dielectric constant
	ExclusionSettings exclude;      //  What electrostatic exclusions should
					//  be made
	BigReal scale14;		//  Scaling factor for 1-4 
					//  electrostatics
	BigReal nonbondedScaling;	//  Scaling factor for nonbonded forces
	int dcdFrequency;		//  How often (in timesteps) should
					//  a DCD trajectory file be updated
	int velDcdFrequency;		//  How often (in timesteps) should
					//  a velocity DCD file be updated
	int xstFrequency;		//  How often (in timesteps) should
					//  a XST trajectory file be updated
	char dcdFilename[129];		//  DCD filename
	char velDcdFilename[129];       //  Velocity DCD filename
	char xstFilename[129];		//  Extended system trajectory filename
	char outputFilename[129];	//  Output file name.  This name will
					//  have .coor appended to it 
					//  for the coordinates and 
					//  .vel appended to
					//  it for the velocities
	char restartFilename[129];	//  Base name of the restart file
	int restartFrequency;		//  How often (in timesteps) shoud the
					//  restart files be updated
	Bool binaryRestart;		//  should restart files be
					//  binary format rather than PDB
	Bool binaryOutput;		//  should output files be
					//  binary format rather than PDB
	BigReal cutoff;			//  Cutoff distance
	BigReal eleccutoff;		//  electrostatic Cutoff distance
	BigReal vdwcutoff;		//  vdw Cutoff distance
	BigReal margin;			//  Fudge factor on patch size
	BigReal patchDimension;		//  Dimension of each side of a patch
					//  This is either cutoff+margin or
					//  pairlistDist+margin depending on
					//  whether or not switching is on
					//  or not
	Bool switchingActive;		//  Flag TRUE->using switching function
					//  for electrostatics and vdw
	BigReal switchingDist;		//  Distance at which switching
					//  becomes active
	BigReal vdwswitchDist;		//  Distance at which vdw switching
					//  becomes active
	BigReal elecswitchDist;		//  Distance at which electrostatic
					//  switching becomes active
	BigReal pairlistDist;		//  Distance within which atom pairs 
					//  should be added to pairlist

        Bool plMarginCheckOn;           //  Should atom movement be checked
	  				//  each tstep to see if an atom
					//  has moved too far for the current
					//  pairlistdist

	Bool constraintsOn;		//  Flag TRUE-> harmonic constraints 
					//  active
	int constraintExp;		//  Exponent for harmonic constraints

        //****** BEGIN selective restraints (X,Y,Z) changes 
        Bool selectConstraintsOn;       //  Flag TRUE-> selective restraints  
                                        //  active
        Bool constrXOn, constrYOn,       
             constrZOn;                 //  Flag TRUE-> select which Cartesian 
                                        //  component to restrain
        //****** END selective restraints (X,Y,Z) changes 
	BigReal constraintScaling;	//  Scaling factor for constraint forces

        //****** BEGIN CHARMM/XPLOR type changes
        Bool paraTypeXplorOn;           //  FLAG TRUE-> parametrs are XPLOR format (default)
        Bool paraTypeCharmmOn;          //  FLAG TRUE-> parametrs are CHARMM format
        //****** END CHARMM/XPLOR type changes

        //****** BEGIN moving constraints changes 
        Bool movingConstraintsOn;       //  Flag TRUE-> moving constraints 
                                        //  active
        zVector movingConsVel;           //  Velocity of the movement, A/timestep
        //****** END moving constraints changes 
        //****** BEGIN rotating constraints changes 
        Bool rotConstraintsOn;          //  Flag TRUE-> rotating constraints 
                                        //  active
        zVector rotConsAxis;             //  Axis of rotation
        zVector rotConsPivot;            //  Pivot point of rotation
        BigReal rotConsVel;             //  Velocity of rotation, Deg/timestep
        //****** END rotating constraints changes 
        //****** BEGIN SMD constraints changes   
        Bool SMDOn;                     //  Flag TRUE-> SMD constraints active
        BigReal SMDVel;                 //  Velocity of the movement, A/timestep
        zVector SMDDir;                  //  Direction of the movement
        BigReal SMDk; 			//  Elastic constant for SMD
 	char SMDFile[128];		//  File for SMD information
        int SMDOutputFreq;              //  Output frequency for SMD constr.
        //****** END SMD constraints changes 

	Bool globalForcesOn;		//  Are global forces present?
	Bool tclForcesOn;		//  Are Tcl forces present?
	Bool freeEnergyOn;		//  Doing free energy perturbation?
	Bool tclOn;			//  Are Tcl scripts present?
	Bool miscForcesOn;		//  Using misc forces?

	Bool fixedAtomsOn;		//  Are there fixed atoms?

	Bool langevinOn;		//  Flag TRUE-> langevin dynamics active
	BigReal langevinTemp;		//  Temperature for Langevin dynamics
	BigReal langevinDamping;	//  Damping coefficient (1/ps)
	Bool langevinHydrogen;		//  Flag TRUE-> apply to hydrogens

	Bool globalOn;			//  Flag TRUE-> use global integrator
	Bool dihedralOn;		//  Flag TRUE-> dihedral dynamics active
	Bool COLDOn;			//  Flag TRUE-> constrained overdamped
					//  langevin dynamics active
	BigReal COLDRate;		//  Damping coefficient for COLD.
	BigReal COLDTemp;		//  Temperature for COLD.

	Bool tCoupleOn;			//  Flag TRUE-> Temperature coupling 
					//  active
	BigReal tCoupleTemp;		//  Temperature for temp coupling

	int rescaleFreq;		//  Velocity rescale frequency
	BigReal rescaleTemp;		//  Temperature to rescale to

	int reassignFreq;		//  Velocity reassignment frequency
	BigReal reassignTemp;		//  Temperature to reassign to
	BigReal reassignIncr;		//  Added to reassignTemp each time
	BigReal reassignHold;		//  Hold reassignTemp at this value

	Bool useGroupPressure;		//  Use group rather than atomic
					//  quantities for pressure calc

	Bool useFlexibleCell;		//  Use anisotropic cell fluctuations

	Bool berendsenPressureOn;	//  Berendsen pressure bath
	BigReal berendsenPressureTarget;
	BigReal berendsenPressureCompressibility;
	BigReal berendsenPressureRelaxationTime;
	int berendsenPressureFreq;

	Bool langevinPistonOn;		//  Langevin piston pressure control
	BigReal langevinPistonTarget;
	BigReal langevinPistonPeriod;
	BigReal langevinPistonDecay;
	BigReal langevinPistonTemp;

        BigReal surfaceTensionTarget;

	zVector strainRate;
	zVector strainRate2; // off diagonal elements (xy, xz, yz)

	unsigned int randomSeed;	//  Seed for random number generator

	Bool FMAOn;                     //  Flag TRUE-> FMA active
	int FMALevels;			//  Number of Levels for FMA
	int FMAMp;			//  Number of multipole terms for FMA
	Bool FMAFFTOn;			//  FFT on/off flag for FMA
	int FMAFFTBlock;		//  FFT blocking factor for FMA

	Bool fullDirectOn;		//  Should direct calculations of
					//  full electrostatics be performed?

	Bool PMEOn;			//  Flag TRUE -> PME active
	BigReal PMETolerance;		//  Direct space tolerance
        BigReal PMEEwaldCoefficient;    //  From tolerance and cutoff
	int PMEInterpOrder;		//  Order of interpolation
	int PMEGridSizeX;		//  No. of grid points in x dim
	int PMEGridSizeY;		//  No. of grid points in y dim
	int PMEGridSizeZ;		//  No. of grid points in z dim

	Bool useDPME;			//  Flag TRUE -> old DPME code

	Bool minimizeOn;		//  Flag TRUE-> minimization active
	BigReal maximumMove;		//  Maximum movement per timestep 
					//  during minimization

	Bool sphericalBCOn;		//  Flag TRUE-> spherical boundary 
					//  conditions are active
	zVector sphericalCenter;		//  Center specified by user
	BigReal sphericalBCk1;		//  First force constant for 
					//  spherical BC
	BigReal sphericalBCk2;		//  Second force constant for 
					//  spherical BC
	BigReal sphericalBCr1;		//  First radius for spherical BC
	BigReal sphericalBCr2;		//  Second radius for spherical BC
	int sphericalBCexp1;		//  First radius for spherical BC
	int sphericalBCexp2;		//  Second radius for spherical BC

        Bool cylindricalBCOn;           //  Flag TRUE->cylindrical boundary
                                        //  conditions are active
        zVector cylindricalCenter;
	char cylindricalBCAxis;		//  'x', 'y', or 'z'
        BigReal cylindricalBCr1;
        BigReal cylindricalBCr2;
        BigReal cylindricalBCl1;
        BigReal cylindricalBCl2;
        int cylindricalBCexp1;
        int cylindricalBCexp2;
        BigReal cylindricalBCk1;
        BigReal cylindricalBCk2;

	Bool eFieldOn;                  //  Should a electric field be applied
	zVector eField;                  //  Electric field vector to be applied

	int outputEnergies;		//  Number of timesteps between energy
					//  outputs

	int outputMomenta;		//  Number of timesteps between momentum
					//  outputs

	int outputTiming;		//  Number of timesteps between timing
					//  outputs

	int outputPressure;		//  Number of timesteps between pressure
					//  tensor outputs

	int firstTimestep;		//  Starting timestep.  Will be 0 unless
					//  restarting a simulation

	MTSChoices MTSAlgorithm;	//  What multiple timestep algorithm
					//  to use

	int longSplitting;		//  What electrostatic splitting 	
					//  to use

	int splitPatch;			// How are patches determined?
	BigReal hgroupCutoff;		// what is the added hydrogen margin?

	int mollyOn;			// mollify long range forces?
	BigReal mollyTol;		// error tolerance for molly
	int mollyIter;			// max number of iterations for molly

        int rigidBonds;                 // what type of rigid bonds to hydrogens
                                        // none, all, or only water

        BigReal rigidTol;               // error tolerance for rigid bonds
        int rigidIter;                  // Number of NR iterations 
	int rigidDie;			// die if rigidTol not achieved

	Bool testOn;			//  Do tests rather than simulation
	Bool commOnly;			//  Don't do any force evaluations

	int totalAtoms;			//  Total Number of atoms in simulation

	//
        // hydrogen bond simulation parameters
        //

        // should the hydrogen bond term be used?  If FALSE, all other
	// hydrogen bond parameters are unnecessary in simulation.
	Bool HydrogenBonds;

	// should the antecedent atom be used in the calculation of hbonds?
	Bool useAntecedent;

	// exponents used in hydrogen bond energy function:
	//   aaAngleExp = exp for H-A-AA angle term (n)
	//   haAngleExp = exp for D-H-A angle term (m)
	//   distAttExp = exp for attractive A-D distance term (j)
	//   distRepExp = exp for repulsive A-D distance term (i)
	int aaAngleExp, haAngleExp, distAttExp, distRepExp;

	// cutoff D-H-A angle, and on/off angles for switch fcn (in degrees)
	BigReal dhaCutoffAngle, dhaOnAngle, dhaOffAngle;

	// cutoff distance for D-A separation in hbonds (in Angstroms), and
	// on/off distances for hbond radial term switching function
	BigReal daCutoffDist, daOnDist, daOffDist;

	// IMD parameters
	int IMDon;
	int IMDport;
 	int IMDfreq;

public:

	SimParameters() {};
	SimParameters(ConfigList *c, char *&cwd) {
	  initialize_config_data(c,cwd);
	};
	~SimParameters() {};

	void initialize_config_data(ConfigList *, char *&cwd);
					//  Initialize SimParameters data
					//  from the ConfigList object
	void send_SimParameters(Communicate *);	
					//  Used by the master process
					//  to send the paramters to
					//  the other processors
	void receive_SimParameters(MIStream *);  
					//  Used by the other processors
					//  to receive the data from the
					//  master process
	void scriptSet(const char *, const char *);
					//  Set parameters at run time

private:

	void config_parser(ParseOptions &opts);

	void config_parser_basic(ParseOptions &opts);
	void config_parser_fileio(ParseOptions &opts);
	void config_parser_fullelect(ParseOptions &opts);
	void config_parser_methods(ParseOptions &opts);
	void config_parser_constraints(ParseOptions &opts);
	void config_parser_boundary(ParseOptions &opts);
	void config_parser_misc(ParseOptions &opts);

	void check_config(ParseOptions &opts, ConfigList *config, char *&cwd);

	void print_config(ParseOptions &opts, ConfigList *config, char *&cwd);

	int fmaFrequency;		//  outdated parameter name
	char loadStrategy[65];		//  Load balancing strategy
};

#endif

