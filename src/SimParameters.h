//-*-c++-*-
/***************************************************************************/
/*        (C) Copyright 1995,1996,1997 The Board of Trustees of the        */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/


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
#define VERLETII	2
#define VERLETX		3

//  The following definitions are used to distinuish between multiple
//  long-short range force splittings
#define SHARP		0
#define XPLOR		1
#define C1		2
#define SKEEL		3

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

//****** BEGIN SMD constraints changes 

// The following definitions are used to distinguish between
// changing moving constraint direction strategies
#define SMD_UNIFORM   0
#define SMD_GAUSSIAN  1
//****** END SMD constraints changes 


class SimParameters
{
private:
public:
  char dummy;
	BigReal dt;	   		//  Timestep size
	int N;		   		//  Number of steps to be performed
	int stepsPerCycle;		//  Number of timesteps per cycle

	Vector cellBasisVector1;	//  Basis vector for periodic cell
	Vector cellBasisVector2;	//  Basis vector for periodic cell
	Vector cellBasisVector3;	//  Basis vector for periodic cell
	Vector cellOrigin;		//  Fixed center of periodic cell
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
	int electForceDcdFrequency;	//  How often (in timesteps) should
					//  long and short range electrostatic
					//  force DCD files be updated
	int allForceDcdFrequency;	//  How often (in timesteps) should
					//  DCD file with total force be updated
	int xstFrequency;		//  How often (in timesteps) should
					//  a XST trajectory file be updated
	int vmdFrequency;		//  How often (in timesteps) should
					//  a VMD trajectory to be collected
	char dcdFilename[129];		//  DCD filename
	char velDcdFilename[129];       //  Velocity DCD filename
	char xstFilename[129];		//  Extended system trajectory filename
	char electForceDcdFilename[129];//  long and short range elect force 
					//  DCD filename
	char allForceDcdFilename[129];  //  File name for DCD file containing
					//  total forces
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

        //****** BEGIN CHARMM/XPLOR type changes
        Bool paraTypeXplorOn;           //  FLAG TRUE-> parametrs are XPLOR format (default)
        Bool paraTypeCharmmOn;          //  FLAG TRUE-> parametrs are CHARMM format
        //****** END CHARMM/XPLOR type changes

        //****** BEGIN moving constraints changes 
        Bool movingConstraintsOn;       //  Flag TRUE-> moving constraints 
                                        //  active
        Vector movingConsVel;           //  Velocity of the movement, A/timestep
        //****** END moving constraints changes 
        //****** BEGIN rotating constraints changes 
        Bool rotConstraintsOn;          //  Flag TRUE-> rotating constraints 
                                        //  active
        Vector rotConsAxis;             //  Axis of rotation
        Vector rotConsPivot;            //  Pivot point of rotation
        BigReal rotConsVel;             //  Velocity of rotation, Deg/timestep
        //****** END rotating constraints changes 
        //****** BEGIN SMD constraints changes   
        Bool SMDOn;                     //  Flag TRUE-> SMD constraints 
                                        //  active
        int SMDExp;                     //  Exponent for SMD constraints
        Vector SMDRefPos;               //  Initial pos. of SMD restraint point
        BigReal SMDk;                   //  SMD force constant (kcal/mol)
        BigReal SMDVel;                 //  Velocity of the movement, A/timestep
        Vector SMDDir;                  //  Direction of the movement
        int SMDAtom;                    //  Index of the atom to be moved
        int SMDOutputFreq;              //  Output frequency for SMD constr.

        int SMDTStamp;                  //  Time of the last refPos change
        Bool SMDChDirOn;                //  Is changing directions on?
        BigReal SMDVmin;                //  Min allowed ave velocity 
        int SMDVminTave;                //  Time for averaging Vmin (timesteps)
        BigReal SMDConeAngle;           //  Angle of the cone to choose 
                                        //    directions from
        int SMDChDirMethod;             //  Distribution to use (uniform/gauss)
        BigReal SMDGaussW;              //  width of gaussian distribution

        Bool SMDChForceOn;              //  Is changing forces on?
        BigReal SMDVmax;                //  Max allowed ave velocity
        int SMDVmaxTave;                //  Time for averaging Vmax (timesteps)
        BigReal SMDFmin;                //  Force (pN) to which to reset
        //****** END SMD constraints changes 

	Bool globalForcesOn;		//  Are global forces present?
	Bool tclForcesOn;		//  Are Tcl forces present?
	Bool freeEnergyOn;		//  Doing free energy perturbation?
	Bool tclOn;			//  Are Tcl scripts present?

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
	Vector strainRate;

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
	int PMEInterpOrder;		//  Order of interpolation
	int PMEGridSizeX;		//  No. of grid points in x dim
	int PMEGridSizeY;		//  No. of grid points in y dim
	int PMEGridSizeZ;		//  No. of grid points in z dim

	Bool minimizeOn;		//  Flag TRUE-> minimization active
	BigReal maximumMove;		//  Maximum movement per timestep 
					//  during minimization

	Bool sphericalBCOn;		//  Flag TRUE-> spherical boundary 
					//  conditions are active
	BigReal sphericalBCk1;		//  First force constant for 
					//  spherical BC
	BigReal sphericalBCk2;		//  Second force constant for 
					//  spherical BC
	BigReal sphericalBCr1;		//  First radius for spherical BC
	BigReal sphericalBCr2;		//  Second radius for spherical BC
	int sphericalBCexp1;		//  First radius for spherical BC
	int sphericalBCexp2;		//  Second radius for spherical BC
	Bool sphericalCenterCOM;	//  Are the spherical boundaries centered
					//  around the center of mass?
	Vector sphericalCenter;		//  Center specified by user

        Bool cylindricalCenterCOM;
        Vector cylindricalCenter;
        Bool cylindricalBCOn;           //  Flag TRUE->cylindrical boundary
                                        //  conditions are active
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
	Vector eField;                  //  Electric field vector to be applied

	int outputEnergies;		//  Number of timesteps between energy
					//  outputs

	int outputMomenta;		//  Number of timesteps between momentum
					//  outputs

	int outputTiming;		//  Number of timesteps between timing
					//  outputs

	int firstTimestep;		//  Starting timestep.  Will be 0 unless
					//  restarting a simulation

	MTSChoices MTSAlgorithm;	//  What multiple timestep algorithm
					//  to use

	int longSplitting;		//  What electrostatic splitting 	
					//  to use

	int splitPatch;			// How are patches determined?
	BigReal hgroupCutoff;		// what is the added hydrogen margin?

        int rigidBonds;                 // what type of rigid bonds to hydrogens
                                        // none, all, or only water


        BigReal rigidTol;               // error tolerance for rigid bonds
        int rigidIter;                  // Number of NR iterations 

	Bool testOn;			//  Do tests rather than simulation

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
	//****** BEGIN SMD constraints changes
	char chDirMethod[65];           // SMD changing direction method
	//****** END SMD constraints changes



};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: SimParameters.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1036 $	$Date: 1999/05/27 19:00:46 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: SimParameters.h,v $
 * Revision 1.1036  1999/05/27 19:00:46  jim
 * Added nonbondedScaling parameter and fixed Tcl scripting bug.
 *
 * Revision 1.1035  1999/05/26 22:23:57  jim
 * Added basic Tcl scripting, fixed bugs in broadcasts.
 *
 * Revision 1.1034  1999/05/14 18:47:27  jim
 * Split up huge functions into several smaller ones to help compilers.
 *
 * Revision 1.1033  1999/03/17 21:26:35  jim
 * Switching internal nomenclature from fmaFrequency to fullElectFrequency.
 *
 * Revision 1.1032  1999/03/10 05:11:35  jim
 * Added reassignHold parameter.
 *
 * Revision 1.1031  1999/03/09 01:44:17  jim
 * Added langevinDamping and langevinHydrogen parameters.
 *
 * Revision 1.1030  1999/02/02 08:02:37  ferenc
 * Added support for CHARMM parameter format in parameter files.
 *
 * Revision 1.1029  1999/01/08 23:24:54  ferenc
 * added selective position restraints for specific Cartesian components
 *
 * Revision 1.1028  1999/01/06 22:50:33  jim
 * Anisotropic (flexible cell) Langevin Piston pressure control finished.
 *
 * Revision 1.1027  1998/11/29 22:01:00  jim
 * Added group-based pressure control to work with rigidBonds.
 * New option useGroupPressure, turned on automatically if needed.
 *
 * Revision 1.1026  1998/10/26 17:14:51  krishnan
 * Added binaryOutput
 *
 * Revision 1.1025  1998/10/01 00:31:30  sergei
 * added rotating restraints feature;
 * changed the moving restraints from only moving one atom to moving all
 * atoms that are restrained. One-atom pulling is available in SMD feature.
 *
 * Revision 1.1024  1998/09/13 21:06:15  jim
 * Cleaned up output, defaults, etc.
 *
 * Revision 1.1023  1998/08/17 23:34:30  jim
 * Added options for Langevin piston and wrapWater.
 *
 * Revision 1.1022  1998/08/04 04:07:23  jim
 * Added extended system file support and fixed lack of endi in SimParameters.
 *
 * Revision 1.1021  1998/08/03 15:31:22  jim
 * Added temperature reassignment.
 *
 * Revision 1.1020  1998/04/06 16:34:11  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1019  1998/03/31 04:55:48  jim
 * Added test mode, fixed errors in virial with full electrostatics.
 *
 * Revision 1.1018  1998/02/17 06:39:25  jim
 * SHAKE/RATTLE (rigidBonds) appears to work!!!  Still needs langevin,
 * proper startup, and degree of freedom tracking.
 *
 * Revision 1.1017  1998/02/10 05:35:06  jim
 * Split ComputeGlobal into different classes and files.
 * Switched globalForces and globalForcesTcl to tclForces and tclForcesScript.
 * Added (soon to be used) freeEnergy and freeEnergyConfig.
 *
 * Revision 1.1016  1998/01/05 20:28:20  sergei
 * Introduced SMD parameters.
 *
 * Revision 1.1015  1997/12/19 23:48:51  jim
 * Added Tcl interface for calculating forces.
 *
 * Revision 1.1014  1997/10/01 16:47:03  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.1013  1997/09/19 08:55:38  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1012  1997/08/18 17:45:11  sergei
 * added moving restraint capability with input from config file
 * (for one atom only)
 *
 * Revision 1.1011  1997/05/29 19:13:16  nealk
 * Added hgroupCutoff
 *
 * Revision 1.1010  1997/04/16 23:44:05  brunner
 * Put ldbStrategy={none|refineonly|alg7}, ldbPeriod, and firstLdbStep
 * in SimParameters.
 *
 * Revision 1.1009  1997/04/08 21:08:51  jim
 * Contant pressure now correct on multiple nodes, should work with MTS.
 *
 * Revision 1.1008  1997/04/04 17:31:43  brunner
 * New charm fixes for CommunicateConverse, and LdbCoordinator data file
 * output, required proxies, and idle time.
 *
 * Revision 1.1007  1997/03/27 17:08:32  nealk
 * Added hydrogen groupings.  Now configuration parameter "splitPatch" determines
 * atom-into-patch distribution.
 *
 * Revision 1.1006  1997/03/27 08:04:26  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1005  1997/03/25 23:01:05  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
 * Revision 1.1004  1997/03/21 23:05:47  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1003  1997/03/16 19:44:08  jim
 * Added cylindricalBCAxis option to cylindrical boundary conditions.
 *
 * Revision 1.1002  1997/03/15 22:15:33  jim
 * Added ComputeCylindricalBC.  Doesn't break anything but untested and
 * cylinder is along x axis (will fix soon).
 *
 * Revision 1.1001  1997/03/04 22:37:19  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1000  1997/02/06 15:59:21  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:31  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/01/28 17:28:52  jim
 * First top-down changes for periodic boundary conditions, added now to
 * avoid conflicts with Ari's migration system.
 *
 * Revision 1.778  1997/01/28 00:31:27  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:42  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:37:02  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.3  1996/12/11 00:05:20  milind
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 04:39:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.42  1996/05/08 20:54:53  gursoy
 * put rigid bonds tolerance parameter
 *
 * Revision 1.41  1996/04/18 18:46:18  billh
 * Updated to read hydrogen bond information.
 *
 * Revision 1.40  1996/03/28 11:09:05  gursoy
 * declared rigidBonds
 *
 * Revision 1.39  96/03/28  10:42:18  10:42:18  gursoy (Attila Gursoy)
 * added rigid options
 * 
 * Revision 1.38  96/01/28  21:50:08  21:50:08  jean (Jean Wolfgang)
 * Attempting to make stable RCS without Mark Nelson's
 * fma/pairlist decoupling
 * 
 * Revision 1.41  1996/01/20 22:13:49  nelson
 * Added SKEEL definition for switching function types
 *
 * Revision 1.40  95/12/06  22:14:12  22:14:12  brunner (Robert Brunner)
 * Took out eswitchingDist, and added eleccutoff, elecswitchDist,
 * vdwcutoff, and vdwswitchDist, for separate cutoffs when switching
 * active
 * 
 * Revision 1.39  1995/12/06 17:12:28  brunner
 * Added fmaTheta, the DPMTA theta parameter.  Previously, it was
 * hard-coded in as 0.715
 *
 * Revision 1.38  1995/12/05 21:12:34  brunner
 * Added eswitchdist, for separate electrostatic and vdw switching
 *
 * Revision 1.37  1995/11/19 21:04:55  nelson
 * Added fmaFrequency parameter
 *
 * Revision 1.36  95/10/27  21:37:16  21:37:16  jim (Jim Phillips)
 * Added global integration methods.
 * Specifically, globalTest, dihedral, COLD, COLDTemp, and COLDRate.
 * 
 * Revision 1.35  95/09/28  13:19:43  13:19:43  brunner (Robert Brunner)
 * Added a new lbdStrategy possible value, "bisection", which causes
 * the recursive bisection distribution algorithm to be executed periodically
 * during the run, to account for atom migration.
 * 
 * Revision 1.34  95/09/26  15:15:08  15:15:08  nelson (Mark T. Nelson)
 * Added all force DCD files and cleaned up symantics of long and short
 * range electrostatic force DCD files.
 * 
 * Revision 1.33  95/09/26  13:27:44  13:27:44  nelson (Mark T. Nelson)
 * Added temperature coupling
 * 
 * Revision 1.32  95/08/30  14:05:04  14:05:04  nelson (Mark T. Nelson)
 * Added options for short range force DCD files and binary coordinate
 * files
 * 
 * Revision 1.31  95/08/28  13:19:31  13:19:31  nelson (Mark T. Nelson)
 * Added option for specifying the long range force splitting to use and
 * also fixed *stupid* bug in forceDCDfile parameter initialization
 * 
 * Revision 1.30  95/08/16  12:47:07  12:47:07  nelson (Mark T. Nelson)
 * Added parameters for specifing a center for the spherical boundary conditions
 * 
 * Revision 1.29  95/08/11  14:52:06  14:52:06  nelson (Mark T. Nelson)
 * Added options for force DCD files and MTS algorithms
 * 
 * Revision 1.28  95/07/25  11:42:23  11:42:23  brunner (Robert Brunner)
 * Added plMarginCheck, variable simParam->plMarginCheckOn, to check
 * to see if an atom has moved in such a way that it may be closer than the
 * cutoff distance.  The default setting is off.
 * 
 * Revision 1.27  95/07/14  14:19:11  14:19:11  nelson (Mark T. Nelson)
 * Added binary velocity restart files
 * 
 * Revision 1.26  95/05/23  14:11:14  14:11:14  nelson (Mark T. Nelson)
 * Added options for electric field force
 * 
 * Revision 1.25  95/04/06  15:48:04  15:48:04  nelson (Mark T. Nelson)
 * Added fullDirectOn option
 * 
 * Revision 1.24  95/04/06  14:38:12  14:38:12  nelson (Mark T. Nelson)
 * Added firstTimestep parameter for restarting
 * 
 * Revision 1.23  95/03/22  11:21:36  11:21:36  nelson (Mark T. Nelson)
 * Added parameters for spherical boundary conditions and outputEnergies,
 * added functions check_duplicate() and check_read_format(), and some
 * other cleanup
 * 
 * Revision 1.22  95/03/08  14:47:40  14:47:40  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.21  95/02/22  14:55:21  14:55:21  nelson (Mark T. Nelson)
 * Made changes for the cwd parameter
 * 
 * Revision 1.20  95/02/03  14:01:53  14:01:53  brunner (Robert Brunner)
 * Added LDBSTRAT_NOLOCAL, as a possible value for ldbStrategy
 * 
 * Revision 1.19  95/01/31  19:47:49  19:47:49  nelson (Mark T. Nelson)
 * Added parameters for velocity rescaling
 * 
 * Revision 1.18  95/01/19  17:18:16  17:18:16  brunner (Robert Brunner)
 * Added ldbStrategy, which can be either none or random, and
 * changed stepsperldbcycle to ldbstepspercycle and sendldbstep to 
 * ldbsendstep.
 * 
 * Revision 1.17  95/01/19  15:28:23  15:28:23  nelson (Mark T. Nelson)
 * Added langveinOn and langevinTemp parameters
 * 
 * Revision 1.16  95/01/12  14:30:20  14:30:20  nelson (Mark T. Nelson)
 * Added randomSeed parameters
 * 
 * Revision 1.15  94/12/19  15:24:30  15:24:30  nelson (Mark T. Nelson)
 * Added options for minization
 * 
 * Revision 1.14  94/12/19  10:10:05  10:10:05  nelson (Mark T. Nelson)
 * Added FMA parameters
 * 
 * Revision 1.13  94/12/01  14:15:59  14:15:59  brunner (Robert Brunner)
 * Added stepsPerLdbCycle and sendLdbStep
 * 
 * Revision 1.12  94/11/29  13:34:31  13:34:31  nelson (Mark T. Nelson)
 * Added FMAOn
 * 
 * Revision 1.11  94/10/19  21:40:49  21:40:49  nelson (Mark T. Nelson)
 * Added constraintExp for harmonic constraints
 * 
 * Revision 1.10  94/10/18  11:42:38  11:42:38  nelson (Mark T. Nelson)
 * Added new parameters for velocity dcd files, constraints, and
 * switching function
 * 
 * Revision 1.9  94/10/12  12:34:30  12:34:30  nelson (Mark T. Nelson)
 * Added option for restartFilename
 * 
 * Revision 1.8  94/10/04  11:53:57  11:53:57  nelson (Mark T. Nelson)
 * Changed allowable values for the exclude parameter
 * 
 * Revision 1.7  94/09/28  16:42:19  16:42:19  gursoy (Attila Gursoy)
 * added vmdFrequency as a public member
 * 
 * Revision 1.6  94/09/12  17:47:50  17:47:50  gursoy (Attila Gursoy)
 * took out the receive-message to Node module for charm++ integration
 * 
 * Revision 1.5  94/08/30  13:59:39  13:59:39  nelson (Mark T. Nelson)
 * Added total atoms parameter
 * 
 * Revision 1.4  94/08/08  15:35:06  15:35:06  nelson (Mark T. Nelson)
 * added necessary header files
 * 
 * Revision 1.3  94/08/04  13:15:35  13:15:35  nelson (Mark T. Nelson)
 * Added stepsPerCycle
 * 
 * Revision 1.2  94/08/01  10:36:54  10:36:54  nelson (Mark T. Nelson)
 * Added send and receive routines
 * 
 * Revision 1.1  94/07/08  13:04:23  13:04:23  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/
