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
 *	$RCSfile: SimParameters.C,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1017 $	$Date: 1997/05/29 19:12:11 $
 *
 ***************************************************************************
 * DESCRIPTION:
 * 	SimParameters is just a glorified structure to hold the global
 * static simulation parameters such as timestep size, cutoff, etc. that
 * are read in from the configuration file.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: SimParameters.C,v $
 * Revision 1.1017  1997/05/29 19:12:11  nealk
 * Modified so hydrogen grouping margin offset is a user-defined option.
 *
 * Revision 1.1016  1997/05/09 18:24:24  nealk
 * 1. Added hydrogen grouping code to improve performance in ComputeNonbondedBase
 *    CODE ONLY WORKS WITH HYDROGEN GROUPING!
 * 2. Increased the hydrogen group cutoff side from 2A to 2.5A -- 2A gave
 *    fractionally different values after 100 iterations.  2.5A gives same numbers.
 * 3. Made migration by hydrogen grouping the default in SimParameters.
 *
 * Revision 1.1015  1997/04/24 18:51:45  nealk
 * Corrected parameter bug with DPMTA options.
 * In particular: FFT Block must be 4 and FFTMp must be a multiple
 * of FFT Block.
 *
 * Revision 1.1014  1997/04/16 23:44:04  brunner
 * Put ldbStrategy={none|refineonly|alg7}, ldbPeriod, and firstLdbStep
 * in SimParameters.
 *
 * Revision 1.1013  1997/04/16 22:12:20  brunner
 * Fixed an LdbCoordinator bug, and cleaned up timing and Ldb output some.
 *
 * Revision 1.1012  1997/04/08 21:08:49  jim
 * Contant pressure now correct on multiple nodes, should work with MTS.
 *
 * Revision 1.1011  1997/03/31 19:23:04  nealk
 * Hard-coded +2 angstrom increase in patch margin when using hydrogen grouping.
 *
 * Revision 1.1010  1997/03/27 17:08:31  nealk
 * Added hydrogen groupings.  Now configuration parameter "splitPatch" determines
 * atom-into-patch distribution.
 *
 * Revision 1.1009  1997/03/27 08:04:24  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1008  1997/03/27 03:16:57  jim
 * Added code to check virial calculation, fixed problems with DPMTA and PBC's.
 *
 * Revision 1.1007  1997/03/25 23:01:03  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
 * Revision 1.1006  1997/03/21 23:05:45  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1005  1997/03/19 18:10:17  nealk
 * Added sorted hydrogen group list to molecule.
 *
 * Revision 1.1004  1997/03/19 11:54:57  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1003  1997/03/16 19:44:06  jim
 * Added cylindricalBCAxis option to cylindrical boundary conditions.
 *
 * Revision 1.1002  1997/03/15 22:15:31  jim
 * Added ComputeCylindricalBC.  Doesn't break anything but untested and
 * cylinder is along x axis (will fix soon).
 *
 * Revision 1.1001  1997/03/04 22:37:18  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1000  1997/02/06 15:59:20  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:30  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.2  1997/02/06 02:35:34  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778.2.1  1997/01/28 17:28:51  jim
 * First top-down changes for periodic boundary conditions, added now to
 * avoid conflicts with Ari's migration system.
 *
 * Revision 1.778  1997/01/28 00:31:26  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:37:01  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.8  1996/12/11 00:04:23  milind
 * *** empty log message ***
 *
 * Revision 1.7  1996/12/10 00:13:12  ari
 * *** empty log message ***
 *
 * Revision 1.6  1996/12/04 17:16:32  jim
 * ComputeNonbondedUtil::select() now caches simulation parameters
 *
 * Revision 1.5  1996/11/21 01:00:49  jim
 * added calls to ComputeNonbondedUtil::select()
 *
 * Revision 1.4  1996/11/13 22:42:43  nealk
 * One more /n and iINFO minor change to make the output look nice.
 *
 * Revision 1.3  1996/11/13 21:35:29  nealk
 * Corrected location of CR and iINFO statements.
 *
 * Revision 1.2  1996/11/11 19:54:09  nealk
 * Modified to use InfoStream instead of Inform.
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.67  1996/05/08 20:54:53  gursoy
 * put rigid bonds tolerance parameter
 *
 * Revision 1.66  1996/04/26 16:05:11  gursoy
 * added rigidBonds to the message sent to other nodes
 *
 * Revision 1.65  1996/04/18 18:46:18  billh
 * Updated to read hydrogen bond information.
 *
 * Revision 1.64  1996/03/28 11:41:31  gursoy
 *  printing rigid options if chosen
 *
 * Revision 1.63  96/03/28  11:09:23  11:09:23  gursoy (Attila Gursoy)
 * added rigid bond options
 * 
 * Revision 1.62  96/03/27  03:09:06  03:09:06  jim (Jim Phillips)
 * fixed bug from unset COLDOn flag
 * 
 * Revision 1.61  1996/03/25 13:56:17  brunner
 * Fixed calculation of switchdist from electswitchdist and vdwswitchdist.
 * This bug had no effect if only switchdist was defined.
 *
 * Revision 1.60  1996/02/22 14:22:45  gursoy
 * the integration scheme were reported wrong. fixed it.
 *
 * Revision 1.59  96/01/28  21:50:08  21:50:08  jean (Jean Wolfgang)
 * Attempting to make stable RCS without Mark Nelson's
 * fma/pairlist decoupling
 * 
 * Revision 1.61  1995/12/06 22:16:16  brunner
 * Added separate vdw and electric force cutoffs and switchDists
 * when switching is on.  These default to cutoff and switchDist,
 * but cutoff will be the max of the two cutoffs, and switchDist
 * will be the min of the two switchDists, if separate values are
 * defined.
 * This makes sense if you look at the switching setup code after
 * the if (switchingActive)
 *
 * Revision 1.60  1995/12/06 17:12:28  brunner
 * Added fmaTheta, the DPMTA theta parameter.  Previously, it was
 * hard-coded in as 0.715
 *
 * Revision 1.59  1995/12/05 21:12:34  brunner
 * Added eswitchdist, for separate electrostatic and vdw switching
 *
 * Revision 1.58  1995/11/19 21:04:41  nelson
 * Added fmaFrequency parameter
 *
 * Revision 1.57  95/10/27  21:36:33  21:36:33  jim (Jim Phillips)
 * Added global integration methods.
 * Specifically, globalTest, dihedral, COLD, COLDTemp, and COLDRate.
 * 
 * Revision 1.56  95/10/17  15:14:47  15:14:47  nelson (Mark T. Nelson)
 * Fixed some minor allocation bugs
 * 
 * Revision 1.55  95/10/13  16:28:28  16:28:28  hazen (Brett Hazen)
 * Updated memory allocation to use C++ new/delete
 * 
 * Revision 1.54  1995/09/28  13:18:39  brunner
 * Added a new lbdStrategy possible value, "bisection", which causes
 * the recursive bisection distribution algorithm to be executed periodically
 * during the run, to account for atom migration.
 *
 * Revision 1.53  95/09/26  15:15:07  15:15:07  nelson (Mark T. Nelson)
 * Added all force DCD files and cleaned up symantics of long and short
 * range electrostatic force DCD files.
 * 
 * Revision 1.52  95/09/26  13:27:43  13:27:43  nelson (Mark T. Nelson)
 * Added temperature coupling
 * 
 * Revision 1.51  95/09/21  17:45:47  17:45:47  billh (Bill Humphrey)
 * Use `\0' instead of NULL when terminating a string (or comparing a
 * character to the string terminator character ... yes they're both 0,
 * but NULL is a void *, and '\0' is a char)
 * 
 * Revision 1.50  95/08/31  13:12:33  13:12:33  nelson (Mark T. Nelson)
 * Changed default values for a few parameters
 * 
 * Revision 1.49  95/08/30  14:04:44  14:04:44  nelson (Mark T. Nelson)
 * Added options for short range force DCD files and binary coordinate
 * files
 * 
 * Revision 1.48  95/08/28  13:18:50  13:18:50  nelson (Mark T. Nelson)
 * Added option for specifying the long range force splitting to use and
 * also fixed *stupid* bug in forceDCDfile parameter initialization
 * 
 * Revision 1.47  95/08/23  16:04:04  16:04:04  nelson (Mark T. Nelson)
 * Fixed problem with lack of initialization of forceDcdFilename variable
 * 
 * Revision 1.46  95/08/16  12:46:49  12:46:49  nelson (Mark T. Nelson)
 * Added parameters for specifing a center for the spherical boundary condition
 * 
 * Revision 1.45  95/08/11  14:51:53  14:51:53  nelson (Mark T. Nelson)
 * Added options for force DCD files and MTS algorithms
 * 
 * Revision 1.44  95/07/26  14:45:59  14:45:59  nelson (Mark T. Nelson)
 * Fixed STUPID error in consexp value
 * 
 * Revision 1.43  95/07/25  11:40:09  11:40:09  brunner (Robert Brunner)
 * Added plMarginCheck, variable simParam->plMarginCheckOn, to check
 * to see if an atom has moved in such a way that it may be closer than the
 * cutoff distance.  The default setting is off.
 * 
 * Revision 1.42  95/07/14  14:19:10  14:19:10  nelson (Mark T. Nelson)
 * Added binary velocity restart files
 * 
 * Revision 1.41  95/05/23  14:11:06  14:11:06  nelson (Mark T. Nelson)
 * Added options for electric field force
 * 
 * Revision 1.40  95/05/21  19:46:38  19:46:38  nelson (Mark T. Nelson)
 * Fixed STUPID mistake in the SphericalBC parameters
 * 
 * Revision 1.39  95/05/19  20:32:05  20:32:05  nelson (Mark T. Nelson)
 * Fixed problem where restart file options were ignored
 * 
 * Revision 1.38  95/05/12  15:24:16  15:24:16  nelson (Mark T. Nelson)
 * Fixed bug in way that minimization and initial velocities were 
 * parsed
 * 
 * Revision 1.37  95/04/26  16:31:30  16:31:30  nelson (Mark T. Nelson)
 * Corrected error in dcdfilename processing
 * 
 * Revision 1.36  95/04/10  15:41:56  15:41:56  nelson (Mark T. Nelson)
 * Fixed uninitized values to make ObjectCenter happy
 * 
 * Revision 1.35  95/04/10  15:13:11  15:13:11  nelson (Mark T. Nelson)
 * Fixed several miscellaneous problems with the incorporation of ParseOptions
 * 
 * Revision 1.34  95/04/10  11:28:35  11:28:35  nelson (Mark T. Nelson)
 * Added include of strlib for AIX to resolve strcasecmp and strncasecmp
 * 
 * Revision 1.33  95/04/10  10:17:16  10:17:16  nelson (Mark T. Nelson)
 * Added ifdef's to handle FMA and MDComm parameters correctly
 * 
 * Revision 1.32  95/04/06  15:47:54  15:47:54  nelson (Mark T. Nelson)
 * Added fullDirectOn option
 * 
 * Revision 1.31  95/04/06  14:38:10  14:38:10  nelson (Mark T. Nelson)
 * Added firstTimestep parameter for restarting
 * 
 * Revision 1.30  95/04/05  15:31:13  15:31:13  nelson (Mark T. Nelson)
 * Added ParseOptions object for MUCH better integrity checking
 * 
 * Revision 1.29  95/03/30  21:57:30  21:57:30  nelson (Mark T. Nelson)
 * Made minor changes for pairlistdist and switching
 * 
 * Revision 1.28  95/03/22  11:22:20  11:22:20  nelson (Mark T. Nelson)
 * Added parameters for spherical boundary conditions and outputEnergies,
 * added functions check_duplicate() and check_read_format(), and some
 * other cleanup
 * 
 * Revision 1.27  95/03/08  14:47:34  14:47:34  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.26  95/02/22  14:55:11  14:55:11  nelson (Mark T. Nelson)
 * Made changes for the cwd parameter
 * 
 * Revision 1.25  95/02/03  14:00:55  14:00:55  brunner (Robert Brunner)
 * Added ldbStrategy=nolocality, which sets the ldbStrategy variable to
 * LDBSTRAT_NOLOCAL
 * 
 * Revision 1.24  95/01/31  19:47:37  19:47:37  nelson (Mark T. Nelson)
 * Added parameters for velocity rescaling
 * 
 * Revision 1.23  95/01/30  15:08:51  15:08:51  nelson (Mark T. Nelson)
 * Set langevinTemp even if langevin dynamics weren't active to
 * get rid of ObjectCenter warnings
 * 
 * Revision 1.22  95/01/19  17:17:11  17:17:11  brunner (Robert Brunner)
 * Added ldbStrategy, which can be either none or random, and
 * changed stepsperldbcycle to ldbstepspercycle and sendldbstep to 
 * ldbsendstep.
 * 
 * Revision 1.21  95/01/19  15:28:12  15:28:12  nelson (Mark T. Nelson)
 * Added langveinOn and langevinTemp parameters
 * 
 * Revision 1.20  95/01/12  14:30:09  14:30:09  nelson (Mark T. Nelson)
 * Added randomSeed parameters
 * 
 * Revision 1.19  94/12/19  15:24:22  15:24:22  nelson (Mark T. Nelson)
 * Added options for minization
 * 
 * Revision 1.18  94/12/19  14:25:59  14:25:59  nelson (Mark T. Nelson)
 * Changed pairlistdist value when FMA is active
 * 
 * Revision 1.17  94/12/19  10:09:58  10:09:58  nelson (Mark T. Nelson)
 * Added FMA parameters
 * 
 * Revision 1.16  94/12/14  16:08:17  16:08:17  brunner (Robert Brunner)
 * Added stepsPerLdbCycle and sendLdbStep correctly this time.
 * 
 * Revision 1.15  94/12/01  14:15:03  14:15:03  brunner (Robert Brunner)
 * Added stepsPerLdbCycle and sendLdbStep
 * 
 * Revision 1.14  94/11/29  13:34:25  13:34:25  nelson (Mark T. Nelson)
 * Added FMAOn
 * 
 * Revision 1.13  94/11/09  10:43:48  10:43:48  nelson (Mark T. Nelson)
 * Added vmdFreq configuration option
 * 
 * Revision 1.12  94/10/31  14:37:39  14:37:39  nelson (Mark T. Nelson)
 * Added messages indicating activation of VMD interface
 * 
 * Revision 1.11  94/10/28  12:50:15  12:50:15  nelson (Mark T. Nelson)
 * Added check for mdcomm option for VMD
 * 
 * Revision 1.10  94/10/19  21:41:04  21:41:04  nelson (Mark T. Nelson)
 * Added constraintExp for harmonic constraints
 * 
 * Revision 1.9  94/10/18  11:42:14  11:42:14  nelson (Mark T. Nelson)
 * Added new parameters for velocity dcd files, constraints, and
 * switching function
 * 
 * Revision 1.8  94/10/12  12:43:30  12:43:30  nelson (Mark T. Nelson)
 * Added check for restartFrequency with no restartFilename
 * 
 * Revision 1.7  94/10/12  12:34:19  12:34:19  nelson (Mark T. Nelson)
 * Added option for restartFilename
 * 
 * Revision 1.6  94/10/04  11:54:09  11:54:09  nelson (Mark T. Nelson)
 * Changed allowable values for the exclude parameter
 * 
 * Revision 1.5  94/09/12  17:48:45  17:48:45  gursoy (Attila Gursoy)
 * took out the receive-message to Node module for charm++ integration
 * 
 * Revision 1.4  94/08/30  13:59:32  13:59:32  nelson (Mark T. Nelson)
 * Added total atoms parameter
 * 
 * Revision 1.3  94/08/04  13:15:54  13:15:54  nelson (Mark T. Nelson)
 * Added stepsPerCycle
 * 
 * Revision 1.2  94/08/01  10:36:39  10:36:39  nelson (Mark T. Nelson)
 * Added send and receive routines
 * 
 * Revision 1.1  94/07/08  13:04:22  13:04:22  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/SimParameters.C,v 1.1017 1997/05/29 19:12:11 nealk Exp $";


#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "ConfigList.h"
#include "SimParameters.h"
#include "ParseOptions.h"
#include "structures.h"
#include "Communicate.h"
#include "Message.h"
#include <stdio.h>
#include "InfoStream.h"
#include <time.h>
#include <unistd.h>

#ifdef AIX
#include "strlib.h"		//  For strcasecmp and strncasecmp
#endif

// #define DEBUGM
#include "Debug.h"

/************************************************************************/
/*									*/
/*			FUNCTION initialize_config_data			*/
/*									*/
/*	This function is used by the master process to populate the     */
/*   simulation parameters from a ConfigList object that is passed to   */
/*   it.  Each parameter is checked to make sure that it has a value    */
/*   that makes sense, and that it doesn't conflict with any other      */
/*   values that have been given.					*/
/*									*/
/************************************************************************/

void SimParameters::initialize_config_data(ConfigList *config, char *&cwd)

{
   ParseOptions opts;   //  Object to check consistency of config file
   int len;		//  String length
   char tmpstr[257];	//  Temporary string
   StringList *current; //  Pointer to config option list
   char loadStrategy[65];//  Load balancing strategy
   char filename[129];  //  Temporary file name

   //  Set all variable to fallback default values.  This is not really
   //  necessary, as we give default values when we set up the ParseOptions
   //  object, but it helps the debuggers figure out we've initialized the
   //  variables.
   HydrogenBonds = FALSE;
   useAntecedent = TRUE;
   aaAngleExp = 2;
   haAngleExp = 4;
   distAttExp = 4;
   distRepExp = 6;
   dhaCutoffAngle = 100.0;
   dhaOnAngle = 60.0;
   dhaOffAngle = 80.0;
   daCutoffDist = 7.5;
   daOnDist = 5.5;
   daOffDist = 6.5;
   
   //  So first we set up the ParseOptions objects so that it will check
   //  all of the logical rules that the configuration file must follow.

   ////// basic options
   opts.require("main", "timestep", "size of the timestep, in fs",
		&dt, 1.0);
   opts.range("timestep", POSITIVE);
   opts.units("timestep", FSEC);

   opts.require("main", "numsteps", "number of timesteps to perform",
		&N);
   opts.range("numsteps", POSITIVE);

   opts.optional("main", "stepspercycle",
      "Number of steps between pairlist generation and FMA execution, if active", 
      &stepsPerCycle, 15);
   opts.range("stepspercycle", POSITIVE);

   opts.optional("main", "cutoff", "local electrostatic and Vdw distance", 
      &cutoff);
   opts.range("cutoff", POSITIVE);
   opts.units("cutoff", ANGSTROM);
   
   opts.require("main", "exclude", "Electrostatic exclusion policy",
		PARSE_STRING);

   opts.optional("exclude", "1-4scaling", "1-4 electrostatic scaling factor",
		 &scale14, 1.0);
   opts.range("1-4scaling", POSITIVE);

   opts.optionalB("main", "switching",
		 "Should a smoothing function be used?", &switchingActive, TRUE);
   
   opts.optional("switching", "switchdist",
		 "Distance for switching function activation",
		 &switchingDist);
   opts.range("switchdist", POSITIVE);
   opts.units("switchdist", ANGSTROM);

   opts.optional("switching", "elecswitchdist",
		 "Distance for electrostatic switching function activation",
		 &elecswitchDist);
   opts.range("elecswitchdist", POSITIVE);
   opts.units("elecswitchdist", ANGSTROM);

   opts.optional("switching", "eleccutoff",
		 "Distance for electrostatic cutoff",
		 &eleccutoff);
   opts.range("eleccutoff", POSITIVE);
   opts.units("eleccutoff", ANGSTROM);

   opts.optional("switching", "vdwswitchdist",
		 "Distance for vdw switching function activation",
		 &vdwswitchDist);
   opts.range("vdwswitchdist", POSITIVE);
   opts.units("vdwswitchdist", ANGSTROM);

   opts.optional("switching", "vdwcutoff",
		 "Distance for vdw cutoff",
		 &vdwcutoff);
   opts.range("vdwcutoff", POSITIVE);
   opts.units("vdwcutoff", ANGSTROM);

   opts.optional("switching", "pairlistdist",  "Pairlist inclusion distance",
		 &pairlistDist);
   opts.range("pairlistdist", POSITIVE);
   opts.units("pairlistdist", ANGSTROM);

   opts.optionalB("main", "plMarginCheck", 
		  "Check atom movement since pairlist made?",
		  &plMarginCheckOn, FALSE);

   opts.optional("main", "temperature", "initial temperature",
		 &initialTemp);
   opts.range("temperature", NOT_NEGATIVE);
   opts.units("temperature", KELVIN);

   opts.optionalB("main", "COMmotion", "should the center of mass move?",
		  &comMove, FALSE);

   opts.optional("main", "dielectric", "dielectric constant",
		 &dielectric, 1.0);
   opts.range("dielectric", POSITIVE); // Hmmm, dielectric < 1 ...

   opts.optional("main", "margin", "Patch width margin", &margin, 1.0);
   opts.range("margin", NOT_NEGATIVE);
   opts.units("margin", ANGSTROM);

   opts.optional("main", "seed", "Initial random number seed", &randomSeed);
   opts.range("seed", POSITIVE);

   opts.optional("main", "outputEnergies", "How often to print energies in timesteps",
		 &outputEnergies, 1);
   opts.range("outputEnergies", POSITIVE);
		 
   opts.optional("main", "MTSAlgorithm", "Multiple timestep algorithm",
		PARSE_STRING);

   opts.optional("main", "longSplitting", "Long range force splitting option",
		PARSE_STRING);

   opts.optional("main", "splitPatch", "Atom into patch splitting option",
		PARSE_STRING);
   opts.optional("main", "hgroupCutoff", "Hydrogen margin", &hgroupCutoff);

   opts.optional("main", "cellBasisVector1", "Basis vector for periodic cell",
		&cellBasisVector1);
   opts.optional("main", "cellBasisVector2", "Basis vector for periodic cell",
		&cellBasisVector2);
   opts.optional("main", "cellBasisVector3", "Basis vector for periodic cell",
		&cellBasisVector3);
   opts.optional("main", "cellOrigin", "Fixed center of periodic cell",
		&cellOrigin);

   opts.optional("main", "rigidBonds", "Rigid bonds to hydrogen",PARSE_STRING);
   opts.optional("main", "rigidTolerance", 
                  "Error tolerance for rigid bonds to hydrogen",&rigidTol);

   opts.optional("main", "nonbondedFreq", "Nonbonded evaluation frequency",
		&nonbondedFrequency, 1);
   opts.range("nonbondedFreq", POSITIVE);

   /////////////// file I/O
   opts.optional("main", "cwd", "current working directory", PARSE_STRING);

   opts.require("main", "coordinates", "initial PDB coordinate file",
		PARSE_STRING);
   opts.optional("main", "velocities",
		 "initial velocities, given as a PDB file", PARSE_STRING);
   opts.optional("main", "binvelocities",
		 "initial velocities, given as a binary restart", PARSE_STRING);
   opts.optional("main", "bincoordinates",
		 "initial coordinates in a binary restart file", PARSE_STRING);
   opts.require("main", "structure", "initial PSF structure file",
		PARSE_STRING);
   opts.require("main", "parameters",
"CHARMm 19 or CHARMm 22 compatable force field file (multiple "
"inputs allowed)", PARSE_MULTIPLES);
   
   opts.require("main", "outputname",
		"prefix for the final PDB position and velocity filenames", 
		outputFilename);

   opts.optional("main", "DCDfile", "DCD trajectory output file name",
		 dcdFilename);
   opts.require("DCDfile", "DCDfreq", "Frequency of DCD trajectory output, in "
		"timesteps", &dcdFrequency);
   opts.range("DCDfreq", POSITIVE);

   opts.optional("main", "velDCDfile", "velocity DCD output file name",
		 velDcdFilename);
   opts.require("velDCDfile", "velDCDfreq", "Frequency of velocity "
		"DCD output, in timesteps", &velDcdFrequency);
   opts.range("velDCDfreq", POSITIVE);
   
   opts.optional("main", "electForceDCDfile", "Electrostatic force DCD output file name",
		 electForceDcdFilename);
   opts.require("electForceDCDfile", "electForceDCDfreq", "Frequency of electorstatic force "
		"DCD output, in timesteps", &electForceDcdFrequency);
   opts.range("electForceDCDfreq", POSITIVE);
   
   opts.optional("main", "allForceDCDfile", "Total force DCD output file name",
		 allForceDcdFilename);
   opts.require("allForceDCDfile", "allForceDCDfreq", "Frequency of total force "
		"DCD output, in timesteps", &allForceDcdFrequency);
   opts.range("allForceDCDfreq", POSITIVE);
   
   opts.optional("main", "restartname", "Prefix for the position and velocity "
		 "PDB files used for restarting", restartFilename);
   opts.require("restartname", "restartfreq", "Frequency of restart file "
		"generation", &restartFrequency);
   opts.range("restartfreq", POSITIVE);
   opts.optionalB("restartname", "binaryrestart", "Specify use of binary restart files ", 
		   &binaryRestart, TRUE);
   
   /////////// FMA options
#ifdef DPMTA
   DebugM(1,"DPMTA setup start\n");
   //  PMTA is included, so really get these values
   opts.optionalB("main", "FMA", "Should FMA be used?", &FMAOn, FALSE);
   opts.optional("FMA", "FMALevels", "Tree levels to use in FMA", &FMALevels,
		 5);
   opts.range("FMALevels", POSITIVE);
   opts.optional("FMA", "FMAMp", "Number of FMA multipoles", &FMAMp, 8);
   opts.range("FMAMp", POSITIVE);
   opts.optionalB("FMA", "FMAFFT", "Use FFT enhancement in FMA?", &FMAFFTOn, TRUE);
   opts.optional("FMAFFT", "FMAFFTBlock", "FFT blocking factor",
		&FMAFFTBlock, 4);
   opts.range("FMAFFTBlock", POSITIVE);
   DebugM(1,"DPMTA setup end\n");
#else
   //  PMTA is NOT included.  So just set all the values to 0.
   FMAOn = FALSE;
   FMALevels = 0;
   FMAMp = 0;
   FMAFFTOn = FALSE;
   FMAFFTBlock = 0;
#endif

   opts.optional("main", "fmaFrequency",
      "Number of steps between full electrostatic executions", 
      &fmaFrequency);
   opts.range("fmaFrequency", POSITIVE);

   opts.optional("main", "fmaTheta",
      "FMA theta parameter value", 
      &fmaTheta,0.715);
   opts.range("fmaTheta", POSITIVE);

   opts.optionalB("main", "FullDirect", "Should direct calculations of full electrostatics be performed?",
		  &fullDirectOn, FALSE);

   /////////// Special Dynamics Methods
   opts.optionalB("main", "minimization", "Should minimization be performed?",
		  &minimizeOn, FALSE);
   opts.optional("minimization", "maximumMove", "Maximum atom movement per step "
		 "during minimization", &maximumMove);
   opts.range("maximumMove", POSITIVE);
   opts.units("maximumMove", ANGSTROM);

   opts.optionalB("main", "Langevin", "Should Langevin dynamics be performed?",
		  &langevinOn, FALSE);
   opts.require("Langevin", "langevinTemp", "Temperature for heat bath in Langevin "
		 "dynamics", &langevinTemp);
   opts.range("langevinTemp", NOT_NEGATIVE);
   opts.units("langevinTemp", KELVIN);
   opts.optional("Langevin", "langevinFile", "PDB file with temperature "
		 "coupling terms (B(i)) (default is the PDB input file)",
		 PARSE_STRING);
   opts.optional("Langevin", "langevinCol", "Column in the langevinFile "
		 "containing the temperature coupling term B(i);\n"
		 "default is 'O'", PARSE_STRING);

   //  Dihedral angle dynamics
   opts.optionalB("main", "globalTest", "Should global integration (for development) be used?",
		&globalOn, FALSE);
   opts.optionalB("main", "dihedral", "Should dihedral angle dynamics be performed?",
		&dihedralOn, FALSE);
   COLDOn = FALSE;
   opts.optionalB("dihedral", "COLD", "Should overdamped Langevin dynamics be performed?",
		&COLDOn, FALSE);
   opts.require("COLD", "COLDTemp", "Temperature for heat bath in COLD",
		&COLDTemp);
   opts.range("COLDTemp", NOT_NEGATIVE);
   opts.units("COLDTemp", KELVIN);
   opts.require("COLD", "COLDRate", "Damping rate for COLD",
		&COLDRate, 3000.0);
   opts.range("COLDRate", NOT_NEGATIVE);

   //  Get the parameters for temperature coupling
   opts.optionalB("main", "tcouple", 
		  "Should temperature coupling be performed?",
		  &tCoupleOn, FALSE);
   opts.require("tcouple", "tCoupleTemp", 
		"Temperature for temperature coupling", &tCoupleTemp);
   opts.range("tCoupleTemp", NOT_NEGATIVE);
   opts.units("tCoupleTemp", KELVIN);
   opts.optional("tCouple", "tCoupleFile", "PDB file with temperature "
		 "coupling terms (B(i)) (default is the PDB input file)",
		 PARSE_STRING);
   opts.optional("tCouple", "tCoupleCol", "Column in the tCoupleFile "
		 "containing the temperature coupling term B(i);\n"
		 "default is 'O'", PARSE_STRING);

   opts.optional("main", "rescaleFreq", "Number of steps between "
		"velocity rescaling", &rescaleFreq);
   opts.range("rescaleFreq", POSITIVE);
   opts.optional("main", "rescaleTemp", "Target temperature for velocity rescaling",
		&rescaleTemp);
   opts.range("rescaleTemp", NOT_NEGATIVE);
   opts.units("rescaleTemp", KELVIN);

   ////  Berendsen pressure bath coupling
   opts.optionalB("main", "BerendsenPressure", 
		  "Should Berendsen pressure bath coupling be performed?",
		  &berendsenPressureOn, FALSE);
   opts.require("BerendsenPressure", "BerendsenPressureTarget",
		"Target pressure for pressure coupling",
		&berendsenPressureTarget);
   // opts.units("BerendsenPressureTarget",);
   opts.require("BerendsenPressure", "BerendsenPressureCompressibility",
		"Isothermal compressibility for pressure coupling",
		&berendsenPressureCompressibility);
   // opts.units("BerendsenPressureCompressibility",);
   opts.require("BerendsenPressure", "BerendsenPressureRelaxationTime",
		"Relaxation time for pressure coupling",
		&berendsenPressureRelaxationTime);
   opts.range("BerendsenPressureRelaxationTime", POSITIVE);
   opts.units("BerendsenPressureRelaxationTime", FSEC);
   opts.optional("BerendsenPressure", "BerendsenPressureFreq",
		"Number of steps between volume rescaling",
		&berendsenPressureFreq, 1);
   opts.range("BerendsenPressureFreq", POSITIVE);

   ////  Harmonic Constraints
   opts.optionalB("main", "constraints", "Are harmonic constraints active?",
		 &constraintsOn, FALSE);
   opts.require("constraints", "consexp", "Exponent for harmonic potential",
		&constraintExp, 2);
   opts.range("consexp", POSITIVE);
   opts.require("constraints", "consref", "PDB file containing reference "
		"positions (defaults to the initial coordinates)",
		PARSE_STRING);
   opts.require("constraints", "conskfile", "PDB file containing force "
		"constaints in one of the columns (defaults to "
		"consref)", PARSE_STRING);
   opts.require("constraints", "conskcol", "Column of conskfile to use "
		"for the force constants (defaults to O)", PARSE_STRING);

   //// Spherical Boundary Conditions
   opts.optionalB("main", "sphericalBC", "Are spherical boundary counditions "
		  "active?", &sphericalBCOn, FALSE);
   opts.require("sphericalBC", "sphericalBCr1", "Radius for first sphere "
		 "potential", &sphericalBCr1);
   opts.range("sphericalBCr1", POSITIVE);
   opts.units("sphericalBCr1", ANGSTROM);
   opts.require("sphericalBC", "sphericalBCk1", "Force constant for first "
		"sphere potential (+ is an inward force, - outward)",
		&sphericalBCk1);
   opts.units("sphericalBCk1", KCAL);
   opts.optional("sphericalBC", "sphericalBCexp1", "Exponent for first "
		"sphere potential", &sphericalBCexp1, 2);
   opts.range("sphericalBCexp1", POSITIVE);
   
   opts.optional("sphericalBCr1", "sphericalBCr2", "Radius for second sphere "
		 "potential", &sphericalBCr2);
   opts.range("sphericalBCr2", POSITIVE);
   opts.units("sphericalBCr2", ANGSTROM);
   opts.require("sphericalBCr2", "sphericalBCk2", "Force constant for second "
		"sphere potential (+ is an inward force, - outward)",
		&sphericalBCk2);
   opts.units("sphericalBCk2", KCAL);
   opts.optional("sphericalBCr2", "sphericalBCexp2", "Exponent for second "
		"sphere potential", &sphericalBCexp2, 2);
   opts.range("sphericalBCexp2", POSITIVE);
   opts.optional("sphericalBC", "sphericalBCCenter", "Center of spherical boundaries",
		&sphericalCenter);

   /////////////// Cylindrical Boundary Conditions
   opts.optionalB("main", "cylindricalBC", "Are cylindrical boundary counditions "
                  "active?", &cylindricalBCOn, FALSE);
   opts.require("cylindricalBC", "cylindricalBCr1", "Radius for first cylinder "
                 "potential", &cylindricalBCr1);
   opts.range("cylindricalBCr1", POSITIVE);
   opts.units("cylindricalBCr1", ANGSTROM);
   opts.require("cylindricalBC", "cylindricalBCk1", "Force constant for first "
                "cylinder potential (+ is an inward force, - outward)",
                &cylindricalBCk1);
   opts.units("cylindricalBCk1", KCAL);
   opts.optional("cylindricalBC", "cylindricalBCexp1", "Exponent for first "
                "cylinder potential", &cylindricalBCexp1, 2);
   opts.range("cylindricalBCexp1", POSITIVE);


// additions beyond those already found in spherical parameters    JJU
   opts.optional("cylindricalBC", "cylindricalBCAxis", "Cylinder axis (defaults to x)",
		PARSE_STRING);
   opts.require ("cylindricalBC", "cylindricalBCl1", "Length of first cylinder",
                 &cylindricalBCl1);
   opts.range("cylindricalBCl1", POSITIVE);
   opts.units("cylindricalBCl1", ANGSTROM);
   opts.optional ("cylindricalBCl1", "cylindricalBCl2", "Length of second cylinder",
                  &cylindricalBCl2);
   opts.range ("cylindricalBCl2", POSITIVE);
   opts.units ("cylindricalBCl2", ANGSTROM);
// end  additions

   opts.optional("cylindricalBCr1", "cylindricalBCr2", "Radius for second cylinder "
                 "potential", &cylindricalBCr2);
   opts.range("cylindricalBCr2", POSITIVE);
   opts.units("cylindricalBCr2", ANGSTROM);
   opts.require("cylindricalBCr2", "cylindricalBCk2", "Force constant for second "
                "cylinder potential (+ is an inward force, - outward)",
                &cylindricalBCk2);
   opts.units("cylindricalBCk2", KCAL);
   opts.optional("cylindricalBCr2", "cylindricalBCexp2", "Exponent for second "
                "cylinder potential", &cylindricalBCexp2, 2);
   opts.range("cylindricalBCexp2", POSITIVE);
   opts.optional("cylindricalBC", "cylindricalBCCenter", "Center of cylindrical boundaries", &cylindricalCenter);

   ///////////////  Electric field options
   opts.optionalB("main", "eFieldOn", "Should and electric field be applied",
                 &eFieldOn, FALSE);
   opts.require("eFieldOn", "eField", "Electric field vector", &eField);
   
   ///////////////  Load balance options
   opts.optional("main", "ldbStrategy", "Load balancing strategy",
		 loadStrategy);
   opts.optional("ldbStrategy", "ldbPeriod",
		 "steps between load balancing", &ldbPeriod);
   opts.range("ldbPeriod", POSITIVE);
   opts.optional("ldbStrategy", "firstLdbStep", 
		 "when to start load balancing",
		 &firstLdbStep);
   opts.range("firstLdbStep", POSITIVE);

   //////  MDComm options
#ifdef MDCOMM
   //  MDComm is included, so really get this value
   opts.optionalB("main", "mdcomm", "Enable external communication?",
		 PARSE_BOOL);
   opts.optional("mdcomm", "vmdfreq", "VMD update frequency", &vmdFrequency, 1);
   opts.range("vmdfreq", POSITIVE);
#else
   //  MDComm is NOT included, so just set the vmdFrequency value
   vmdFrequency = -1;
#endif

   /////  Restart timestep option
   opts.optional("main", "firsttimestep", "Timestep to start simulation at",
		 &firstTimestep, 0);
   opts.range("firsttimestep", NOT_NEGATIVE);
 

   ///////////////  hydrogen bond computation options
   opts.optionalB("main", "hbonds", "Use explicit hydrogen bond term",
                 &HydrogenBonds, FALSE);
   opts.optionalB("hbonds","hbAntecedents","Include Antecedent in hbond term",
                 &useAntecedent, TRUE);
   opts.optional("hbonds","hbAAexp","Hbond AA-A-H angle cos exponential",
                 &aaAngleExp, 2);
   opts.optional("hbonds","hbHAexp","Hbond D-H-A angle cos exponential",
                 &haAngleExp, 4);
   opts.optional("hbonds","hbDistAexp","Hbond A-D dist attractive exponential",
                 &distAttExp, 4);
   opts.optional("hbonds","hbDistRexp","Hbond A-D dist repulstive exponential",
                 &distRepExp, 6);
   opts.optional("hbonds","hbCutoffAngle","Hbond D-H-A cutoff angle",
                 &dhaCutoffAngle, 100.0);
   opts.range("hbCutoffAngle", NOT_NEGATIVE);
   opts.optional("hbonds","hbOnAngle","Hbond D-H-A switch function on angle",
                 &dhaOnAngle, 60.0);
   opts.range("hbOnAngle", NOT_NEGATIVE);
   opts.optional("hbonds","hbOffAngle","Hbond D-H-A switch function off angle",
                 &dhaOffAngle, 80.0);
   opts.range("hbOffAngle", NOT_NEGATIVE);
   opts.optional("hbonds","hbCutoffDist","Hbond A-D cutoff distance",
                 &daCutoffDist, 7.5);
   opts.range("hbCutoffDist", POSITIVE);
   opts.units("hbCutoffDist", ANGSTROM);
   opts.optional("hbonds","hbOnDist","Hbond A-D switch function on distance",
                 &daOnDist, 5.5);
   opts.range("hbOnDist", POSITIVE);
   opts.units("hbOnDist", ANGSTROM);
   opts.optional("hbonds","hbOffDist","Hbond A-D switch function off distance",
                 &daOffDist, 6.5);
   opts.range("hbOffDist", POSITIVE);
   opts.units("hbOffDist", ANGSTROM);


   ///////////////////////////////// check the internal consistancy
   if (!opts.check_consistancy()) 
   {
      NAMD_die("Internal error in configuration file parser");
   }

   // Now, feed the object with the actual configuration options through the
   // ParseOptions file and make sure everything is OK
   if (!opts.set(*config)) 
   {
      NAMD_die("ERROR(S) IN THE CONFIGURATION FILE");
   }

   //// now do more logic stuff that can't be done by the ParseOptions object


   //  Make sure that both a temperature and a velocity PDB were
   //  specified
   if (opts.defined("temperature") &&
       (opts.defined("velocities") || opts.defined("binvelocities")) ) 
   {
      NAMD_die("Cannot specify both an initial temperature and a velocity file");
   }

   //  Check for frequencies
   if (!opts.defined("dcdfreq"))
   {
	dcdFrequency = -1;
   }

   if (!opts.defined("veldcdfreq"))
   {
	velDcdFrequency = -1;
   }
   
   if (!opts.defined("electforcedcdfreq"))
   {
        electForceDcdFrequency = -1;
   }

   if (!opts.defined("allforcedcdfreq"))
   {
        allForceDcdFrequency = -1;
   }

   if (!opts.defined("restartname"))
   {
	restartFrequency = -1;
	binaryRestart = FALSE;
   }

   //  Take care of filenames
   if (!opts.defined("dcdfile"))
   {
	dcdFilename[0] = STRINGNULL;
   }

   if (!opts.defined("veldcdfile"))
   {
	velDcdFilename[0] = STRINGNULL;
   }

   if (!opts.defined("restartname"))
   {
	restartFilename[0] = STRINGNULL;
   }

   if (!opts.defined("electForceDCDfile"))
   {
	electForceDcdFilename[0] = STRINGNULL;
   }

   if (!opts.defined("allForceDCDfile"))
   {
	allForceDcdFilename[0] = STRINGNULL;
   }
   
   //  If minimization isn't on, must have a temp or velocity
   if (!minimizeOn && !opts.defined("temperature") && 
       !opts.defined("velocities") && !opts.defined("binvelocities") ) 
   {
      NAMD_die("Must have either an initial temperature or a velocity file");
   }

   if (opts.defined("velocities") || opts.defined("binvelocities") )
   {
	initialTemp = -1.0;
   }

   ///// periodic cell parameters

   /* Save for more flexible PBCs in future
   if ( opts.defined("cellBasisVector3") &&
	! opts.defined("cellBasisVector2") )
   {
	NAMD_die("Used cellBasisVector3 without cellBasisVector2!");
   }

   if ( opts.defined("cellBasisVector2") &&
	! opts.defined("cellBasisVector1") )
   {
	NAMD_die("Used cellBasisVector2 without cellBasisVector1!");
   }
   */

   if (  cellBasisVector1.y || cellBasisVector1.z
      || cellBasisVector2.x || cellBasisVector2.z
      || cellBasisVector3.x || cellBasisVector3.y )
   {
	NAMD_die("Basis vectors 1/2/3 must align with x/y/z axes.");
   }

   if (  cellBasisVector1.x < 0
      || cellBasisVector2.y < 0
      || cellBasisVector3.z < 0 )
   {
	NAMD_die("Basis vector elements must be positive.");
   }

   lattice.set(cellBasisVector1,cellBasisVector2,cellBasisVector3,cellOrigin);

   char s[129];

   ///// cylindricalBC stuff
   if ( ! opts.defined("cylindricalBCAxis") )
   {
      cylindricalBCAxis = 'x';
   }
   else
   {
     opts.get("cylindricalBCAxis", s);

     if (!strcasecmp(s, "x"))
     {
      cylindricalBCAxis = 'x';
     }
     else if (!strcasecmp(s, "y"))
     {
      cylindricalBCAxis = 'y';
     }
     else if (!strcasecmp(s, "z"))
     {
      cylindricalBCAxis = 'z';
     }
   else
     {
      char err_msg[128];

      sprintf(err_msg, "Illegal value '%s' for 'cylindricalBCAxis' in configuration file", s);
      NAMD_die(err_msg);
     }
   }

   ///// exclude stuff
   opts.get("exclude", s);

   if (!strcasecmp(s, "none"))
   {
      exclude = NONE;
   }
   else if (!strcasecmp(s, "1-2"))
   {
      exclude = ONETWO;
   }
   else if (!strcasecmp(s, "1-3"))
   {
      exclude = ONETHREE;
   }
   else if (!strcasecmp(s, "1-4"))
   {
      exclude = ONEFOUR;
   }
   else if (!strcasecmp(s, "scaled1-4"))
   {
      exclude = SCALED14;
   }
   else
   {
      char err_msg[128];

      sprintf(err_msg, "Illegal value '%s' for 'exclude' in configuration file",
	 s);
      NAMD_die(err_msg);
   }

   if (opts.defined("1-4scaling") &&
       exclude != SCALED14)
   {
      iout << iWARN << "Value given for '1-4scaling', but 1-4 scaling "
	      << "not in effect.\n This value will be ignored!" << endi;
   }

   //  Get multiple timestep integration scheme
   if (!opts.defined("MTSAlgorithm"))
   {
	MTSAlgorithm = NAIVE;
   }
   else
   {
	opts.get("MTSAlgorithm", s);

	if (!strcasecmp(s, "naive"))
	{
		MTSAlgorithm = NAIVE;
	}
	else if (!strcasecmp(s, "verleti"))
	{
		MTSAlgorithm = VERLETI;
	}
	else if (!strcasecmp(s, "verletii"))
	{
		MTSAlgorithm = VERLETII;
	}
	else if (!strcasecmp(s, "verletx"))
	{
		MTSAlgorithm = VERLETX;
	}
	else
	{
		char err_msg[129];

		sprintf(err_msg, 
		   "Illegal value '%s' for 'MTSAlgorithm' in configuration file", 
		   s);
		NAMD_die(err_msg);
	}
   }

   //  Get the atom-into-patch splitting specification
   if (!opts.defined("hgroupCutoff"))
   {
	hgroupCutoff = 2.5;  // add 2.5 angstroms
   }
   if (!opts.defined("splitPatch"))
   {
	splitPatch = SPLIT_PATCH_HYDROGEN;
   }
   else
   {
	opts.get("splitPatch", s);
	if (!strcasecmp(s, "position"))
		splitPatch = SPLIT_PATCH_POSITION;
	else if (!strcasecmp(s,"hydrogen"))
		{
		splitPatch = SPLIT_PATCH_HYDROGEN;
		// increase margin by 1 hydrogen bond length
		margin += hgroupCutoff;	// assume no greater than 2 angstroms
		}
	else
	{
		char err_msg[129];
		sprintf(err_msg, 
		   "Illegal value '%s' for 'splitPatch' in configuration file", 
		   s);
		NAMD_die(err_msg);
	}
   }

   //  Get the long range force splitting specification
   if (!opts.defined("longSplitting"))
   {
	longSplitting = SHARP;
   }
   else
   {
	opts.get("longSplitting", s);
	if (!strcasecmp(s, "sharp"))
		longSplitting = SHARP;
	else if (!strcasecmp(s, "xplor"))
		longSplitting = XPLOR;
	else if (!strcasecmp(s, "c1"))
		longSplitting = C1;
	else
	{
		char err_msg[129];

		sprintf(err_msg, 
		   "Illegal value '%s' for 'longSplitting' in configuration file", 
		   s);
		NAMD_die(err_msg);
	}
   }


   // take care of rigid bond options
   if (!opts.defined("rigidBonds"))
   {
      rigidBonds = RIGID_NONE;
      rigidTol   = 0.0;
   }
   else
   {
      opts.get("rigidBonds", s); 
      if (!strcasecmp(s, "all"))
      {
          rigidBonds = RIGID_ALL;
      }
      else  if (!strcasecmp(s, "water"))
      {
           rigidBonds = RIGID_WATER;
      } 
      else
      {
           rigidBonds = RIGID_NONE;
      }

      if (!opts.defined("rigidTolerance"))
      {
         rigidTol = 0.00001;
         iout << iWARN << "rigidBonds is on but rigidTol is not provided:" ;
         iout << iWARN << " using default value" << endi;
      }
   }
   
   //  Take care of switching stuff
   if (switchingActive)
   {

     if (opts.defined("cutoff") && opts.defined("eleccutoff") &&
	 opts.defined("vdwcutoff")) {
       iout << iWARN
	 << "PARAMETERS cutoff, eleccutoff, and vdwcutoff\n"
	 << "all defined -- cutoff ignored\n" 
	 << endi;
     }

     if (opts.defined("switchDist") && opts.defined("elecswitchDist") &&
	 opts.defined("vdwswitchDist")) {
       iout << iWARN
	 << "PARAMETERS switchDist, elecswitchDist, and vdwswitchDist\n"
	 << "all defined -- switchDist ignored\n" 
	 << endi;
     }

     // Check if enough things are defined

     if (!opts.defined("eleccutoff")) {
       if (opts.defined("cutoff"))
	 eleccutoff=cutoff;
       else {
	 char err_msg[129];

	 sprintf(err_msg,
		 "Either cutoff or eleccutoff must be defined");
	 NAMD_die(err_msg);
       }
     }

     if (!opts.defined("vdwcutoff")) {
       if (opts.defined("cutoff"))
	 vdwcutoff=cutoff;
       else {
	 char err_msg[129];

	 sprintf(err_msg,
		 "Either cutoff or vdwcutoff must be defined");
	 NAMD_die(err_msg);
       }
     }

     if (!opts.defined("elecswitchDist")) {
       if (opts.defined("switchDist"))
	 elecswitchDist=switchingDist;
       else {
	 char err_msg[129];

	 sprintf(err_msg,
		 "Either switchDist or elecswitchDist must be defined");
	 NAMD_die(err_msg);
       }
     }

     if (!opts.defined("vdwswitchDist")) {
       if (opts.defined("switchDist"))
	 vdwswitchDist=switchingDist;
       else {
	 char err_msg[129];

	 sprintf(err_msg,
		 "Either switchDist or vdwswitchDist must be defined");
	 NAMD_die(err_msg);
       }
     }

     // Let cutoff be the max of vdw and elec.
     if (eleccutoff > vdwcutoff)
       cutoff = eleccutoff;
     else cutoff = vdwcutoff;

     // Let switchDist be the min of vdw and elec.
     if (elecswitchDist < vdwswitchDist)
       switchingDist = elecswitchDist;
     else switchingDist = vdwswitchDist;

     if ( (elecswitchDist>cutoff) || (elecswitchDist<0) )
     {
       char err_msg[129];

       sprintf(err_msg, 
	       "Switching distance 'elecswitchDist' muct be between 0 and the eleccutoff, which is %f", eleccutoff);
       NAMD_die(err_msg);
     }

     if ( (vdwswitchDist>cutoff) || (vdwswitchDist<0) )
     {
       char err_msg[129];

       sprintf(err_msg, 
	       "Switching distance 'vdwswitchDist' muct be between 0 and the vdwcutoff, which is %f", vdwcutoff);
       NAMD_die(err_msg);
     }
   }

   if (!opts.defined("pairlistDist"))
   {
	pairlistDist = cutoff;
   }
   else if (pairlistDist < cutoff)
   {
	NAMD_die("pairlistDist must be >= cutoff distance");
   }

   patchDimension = pairlistDist + margin;

   //  Turn on global integration if not explicitly specified

   if ( dihedralOn ) globalOn = TRUE;

   //  Make sure modes don't conflict
   if (minimizeOn && langevinOn) 
   {
      NAMD_die("Minimization and Langevin dynamics are mutually exclusive dynamics modes");
   }

   if (minimizeOn && tCoupleOn) 
   {
      NAMD_die("Minimization and temperature coupling are mutually exclusive dynamics modes");
   }

   if (langevinOn && tCoupleOn)
   {
      NAMD_die("Langevin dynamics and temperature coupling are mutually exclusive dynamics modes");
   }

   if (tCoupleOn && opts.defined("rescaleFreq") )
   {
      NAMD_die("Temperature coupling and temperature rescaling are mutually exclusive");
   }

   if (globalOn && comm->nodes() > 1)
   {
      NAMD_die("Global integration does not run in parallel (yet).");
   }

   if (COLDOn && langevinOn)
   {
      NAMD_die("COLD and Langevin dynamics are mutually exclusive dynamics modes");
   }
   if (COLDOn && minimizeOn)
   {
      NAMD_die("COLD and minimization are mutually exclusive dynamics modes");
   }
   if (COLDOn && tCoupleOn)
   {
      NAMD_die("COLD and temperature coupling are mutually exclusive dynamics modes");
   }
   if (COLDOn && opts.defined("rescaleFreq"))
   {
      NAMD_die("COLD and velocity rescaling are mutually exclusive dynamics modes");
   }

   //  Set the default value for the maximum movement parameter
   //  for minimization
   if (opts.defined("minimization") && !opts.defined("maximumMove")) 
   {
      maximumMove = 0.75 * pairlistDist/stepsPerCycle;
   }

   if (!opts.defined("minimization"))
   {
	maximumMove = 0;
   }

   if (opts.defined("rescaleFreq"))
   {
	if (!opts.defined("rescaleTemp"))
	{
		if (opts.defined("temperature"))
		{
			rescaleTemp = initialTemp;
		}
		else
		{
			NAMD_die("Must give a rescale temperature if rescaleFreq is defined");
		}
	}
   }
   else
   {
	rescaleFreq = -1;
	rescaleTemp = 0.0;
   }

   if (opts.defined("rescaleTemp"))
   {
	if (!opts.defined("rescaleFreq"))
	{
		NAMD_die("Must give a rescale temperature if rescaleFreq is given");
	}
   }

   if (minimizeOn && rescaleFreq > 0) 
   {
  	NAMD_die("Minimization and temperature rescaling are mutually exclusive dynamics modes");
   }

   if (!opts.defined("seed")) 
   {
      randomSeed = (unsigned int) time(NULL);
   }

   //  Now take care of cwd processing
   if (opts.defined("cwd"))
   {
	//  First allocate and get the cwd value
	current = config->find("cwd");

	len = strlen(current->data);

	if ( chdir(current->data) )
 	{
		NAMD_die("chdir() to given cwd failed!");
	}

	if (current->data[len-1] != '/')
		len++;

	cwd = new char[len+1];

	if (cwd == NULL)
	{
	   NAMD_die("memory allocation failed in SimParameters::initialize_config_data");
	}

	strcpy(cwd, current->data);

	if (current->data[strlen(current->data)-1] != '/')
		strcat(cwd, "/");

	//  The following checks and prepends should be unnecessary
	//  if we just change directories to "cwd"!

	//  Now check the file names given and see if we should
	//  prepend the cwd value

	if (opts.defined("dcdfile") && (dcdFilename[0] != '/') )
	{
		sprintf(tmpstr, "%s%s", cwd, dcdFilename);
		strcpy(dcdFilename, tmpstr);
	}

	if (opts.defined("veldcdfile") && (velDcdFilename[0] != '/') )
	{
		sprintf(tmpstr, "%s%s", cwd, velDcdFilename);
		strcpy(velDcdFilename, tmpstr);
	}
	
	if (opts.defined("electforcedcdfile") && (electForceDcdFilename[0] != '/') )
	{
		sprintf(tmpstr, "%s%s", cwd, electForceDcdFilename);
		strcpy(electForceDcdFilename, tmpstr);
	}

	if (opts.defined("allforcedcdfile") && (allForceDcdFilename[0] != '/') )
	{
		sprintf(tmpstr, "%s%s", cwd, allForceDcdFilename);
		strcpy(allForceDcdFilename, tmpstr);
	}

	if (opts.defined("restartname") && (restartFilename[0] != '/') )
	{
		sprintf(tmpstr, "%s%s", cwd, restartFilename);
		strcpy(restartFilename, tmpstr);
	}

	//  output name must ALWAYS be given, so don't need to check if it
	//  is defined or not
	if (outputFilename[0] != '/')
	{
		sprintf(tmpstr, "%s%s", cwd, outputFilename);
		strcpy(outputFilename, tmpstr);
	}
   }

   //  Set up load balancing variables
   if (opts.defined("ldbStrategy"))
   {
	//  Assign the load balancing strategy
	if (strcasecmp(loadStrategy, "none") == 0)
	{
	    ldbStrategy=LDBSTRAT_NONE;
	}
	else if (strcasecmp(loadStrategy, "refineonly") == 0)
	{
	    ldbStrategy=LDBSTRAT_REFINEONLY;
	}
	else if (strcasecmp(loadStrategy, "alg7") == 0)
	{
	    ldbStrategy=LDBSTRAT_ALG7;
	}
	else if (strcasecmp(loadStrategy, "other") == 0)
	{
	    ldbStrategy=LDBSTRAT_OTHER;
	}
	else
	{
	    NAMD_die("Unknown ldbStrategy selected");
	}

	if (!opts.defined("ldbPeriod"))
	{
		ldbPeriod=10*stepsPerCycle;
	}

	//  Set default values
	if (!opts.defined("firstLdbStep"))
	{
		firstLdbStep=ldbPeriod;
	}
   }
   else
   {
	ldbStrategy=LDBSTRAT_NONE;
	ldbPeriod=stepsPerCycle;
	firstLdbStep=ldbPeriod;
   }

#ifdef MDCOMM
   if (!opts.defined("mdcomm"))
   {
	vmdFrequency = -1;
   }
#endif

   if (firstTimestep >= N)
   {
	NAMD_die("First timestep must be less than number of steps to perform!!!");
   }

   if ( (firstTimestep%stepsPerCycle) != 0)
   {
	NAMD_die("First timestep must be a multiple of stepsPerCycle!!");
   }

   if (FMAOn && fullDirectOn)
   {
	NAMD_die("Can't do FMA and full electrostatic calculations!!!");
   }

   //  Take care of initializing FMA values to something if FMA is not
   //  active
   if (!FMAOn)
   {
   	FMALevels = 0;
   	FMAMp = 0;
   	FMAFFTOn = FALSE;
   	FMAFFTBlock = 0;
   }
   else
   {
	// idiot checking: frm bug reported by Tom Bishop.
	// DPMTA requires: (#terms in fma)/(fft blocking factor) = integer.
	if (FMAFFTBlock != 4)
		NAMD_die("FMAFFTBlock: Block length must be 4 for short FFT's");
	if (FMAMp % FMAFFTBlock != 0)
		NAMD_die("FMAMp: multipole term must be multiple of block length (FMAFFTBlock)");
    }

   if (!FMAOn && !fullDirectOn)
   {
	fmaFrequency = 0;
   }
   else
   {
	if (!opts.defined("fmaFrequency"))
	{
		fmaFrequency = stepsPerCycle;
	}
	else
	{
		if ( (fmaFrequency > stepsPerCycle) || ( (stepsPerCycle % fmaFrequency) != 0) )
		{
			NAMD_die("stepsPerCycle must be a multiple of fmaFrequency");
		}
	}

	if (!opts.defined("fmaTheta"))
	  fmaTheta=0.715;  /* Suggested by Duke developers */
   }

   if (!opts.defined("constraints"))
   {
	constraintExp = 0;
   }

   if (!sphericalBCOn)
   {
	sphericalBCr1 = 0.0;
	sphericalBCk1 = 0.0;
	sphericalBCexp1 = 0;
	sphericalBCr2 = 0.0;
	sphericalBCk2 = 0.0;
	sphericalBCexp2 = 0;
   }
   else if (!opts.defined("sphericalBCr2"))
   {
      sphericalBCr2 = -1.0;
      sphericalBCk2 = 0.0;
      sphericalBCexp2 = 0;
   }


   if (opts.defined("sphericalBCCenter"))
   {
	sphericalCenterCOM = FALSE;
   }
   else
   {
	sphericalCenterCOM = TRUE;
	sphericalCenter.x = 0.0;
	sphericalCenter.y = 0.0;
	sphericalCenter.z = 0.0;
   }

   if (!cylindricalBCOn)
   {
    cylindricalBCr1 = 0.0;
    cylindricalBCk1 = 0.0;
    cylindricalBCexp1 = 0;
    cylindricalBCr2 = 0.0;
    cylindricalBCk2 = 0.0;
    cylindricalBCexp2 = 0;
    cylindricalBCl1 = 0.0;
    cylindricalBCl2 = 0.0;
   }
   else if (!opts.defined("cylindricalBCr2"))
   {
    cylindricalBCr2 = -1.0;
    cylindricalBCk2 = 0.0;
    cylindricalBCexp2 = 0;
    cylindricalBCl2 = 0.0;
   }
   if (opts.defined ("cylindricalBCCenter"))
   {
    cylindricalCenterCOM = FALSE;
    cylindricalCenter.x = 0.0;
    cylindricalCenter.y = 0.0;
    cylindricalCenter.z = 0.0;
   }
   else
   {
    cylindricalCenterCOM = TRUE;
   }

   if (!eFieldOn)
   {
        eField.x = 0.0;
        eField.y = 0.0;
        eField.z = 0.0;
   }

   if (!opts.defined("langevin"))
   {
	langevinTemp = 0.0;
   }

   if (!opts.defined("tcouple"))
   {
	tCoupleTemp = 0.0;
   }

   if (HydrogenBonds)
   {
     if (daCutoffDist > pairlistDist)
       NAMD_die("Hydrogen bond cutoff distance must be <= pairlist distance");
   }


   //  Now that we have read everything, print it out so that
   //  the user knows what is going on
   iout << iINFO << "SIMULATION PARAMETERS:\n";
   iout << iINFO << "TIMESTEP               " << dt << "\n" << endi;
   iout << iINFO << "NUMBER OF STEPS        " << N << "\n";
   iout << iINFO << "STEPS PER CYCLE        " << stepsPerCycle << "\n";

   iout << iINFO << "PERIODIC CELL          " << lattice.dimension() << "\n";
   iout << iINFO << "PERIODIC CELL CENTER   " << lattice.origin() << "\n";

   if (ldbStrategy==LDBSTRAT_NONE)  {
     iout << iINFO << "LOAD BALANCE STRATEGY  none\n";
   } else {
     if (ldbStrategy==LDBSTRAT_REFINEONLY) {
       iout << iINFO << "LOAD BALANCE STRATEGY  Refine-only\n";
     } else if (ldbStrategy==LDBSTRAT_ALG7)  {
       iout << iINFO << "LOAD BALANCE STRATEGY  Alg7\n";
     } else if (ldbStrategy==LDBSTRAT_OTHER)  {
       iout << iINFO << "LOAD BALANCE STRATEGY  Other\n";
     }
     iout << iINFO << "LDB PERIOD             " << ldbPeriod << " steps\n";
     iout << iINFO << "FIRST LDB TIMESTEP     " << firstLdbStep 
	  << "\n" << endi;
   }
   
   if (initialTemp < 0)
   {
	current = config->find("velocities");

	if (current == NULL)
	{
		current = config->find("binvelocities");
	}

	if ( (cwd == NULL) || (current->data[0] == '/') )
	{
		strcpy(filename, current->data);
	}
	else
	{
		strcpy(filename, cwd);
		strcat(filename, current->data);
	}

	iout << iINFO << "VELOCITY FILE          " << filename << endi;
   }
   else
   {
	iout << iINFO << "INITIAL TEMPERATURE    " 
		 << initialTemp << "\n";
   }

   iout << iINFO << "CENTER OF MASS MOVING? ";

   if (comMove)
   {
   	iout << "YES\n";
   }
   else
   {
   	iout << "NO\n";
   }

   iout << iINFO << "DIELECTRIC             " 
   	 << dielectric << "\n";

   iout << iINFO << "EXCLUDE                ";

   switch (exclude)
   {
   	case NONE:
   		iout << "NONE\n";
   		break;
   	case ONETWO:
   		iout << "ONETWO\n";
   		break;
   	case ONETHREE:
   		iout << "ONETHREE\n";
   		break;
   	case ONEFOUR:
   		iout << "ONE-FOUR\n";
   		break;
   	default:
   		iout << "SCALED ONE-FOUR\n";
   		break;
   }

   if (exclude == SCALED14)
   {
   	iout << iINFO << "1-4 SCALE FACTOR       " << scale14 << "\n" << endi;
   }

   if (dcdFrequency > 0)
   {
   	iout << iINFO << "DCD FILENAME           " 
   		 << dcdFilename << "\n";
   	iout << iINFO << "DCD FREQUENCY          " 
   		 << dcdFrequency << "\n";
   }
   else
   {
   	iout << iINFO << "NO DCD TRAJECTORY OUTPUT\n";
   }
   
   if (velDcdFrequency > 0)
   {
   	iout << iINFO << "VELOCITY DCD FILENAME  " 
   		 << velDcdFilename << "\n";
   	iout << iINFO << "VELOCITY DCD FREQUENCY " 
   		 << velDcdFrequency << "\n";
   }
   else
   {
   	iout << iINFO << "NO VELOCITY DCD OUTPUT\n";
   }
   
   if (electForceDcdFrequency > 0)
   {
   	iout << iINFO << "ELECT FORCE DCD NAME   " 
   		 << electForceDcdFilename << "\n";
   	iout << iINFO << "ELECT FORCE DCD FREQ   " 
   		 << electForceDcdFrequency << endi;
   }

   if (allForceDcdFrequency > 0)
   {
   	iout << iINFO << "TOTAL FORCE DCD NAME   " 
   		 << allForceDcdFilename << "\n";
   	iout << iINFO << "TOTAL FORCE DCD FREQ   " 
   		 << allForceDcdFrequency << "\n";
   }
   	
   iout << iINFO << "OUTPUT FILENAME        " 
   	 << outputFilename << "\n";

   if (restartFrequency == -1)
   {
   	iout << iINFO << "NO RESTART FILE\n";
   }
   else
   {
   	iout << iINFO << "RESTART FILENAME       "
   		 << restartFilename << "\n";
   	iout << iINFO << "RESTART FREQUENCY      " 
   		 << restartFrequency << "\n";

	if (binaryRestart)
	{
		iout << iINFO << "BINARY RESTART FILES WILL BE USED\n";
	}
   }
   
   if (switchingActive)
   {
      iout << iINFO << "SWITCHING ACTIVE\n";
      iout << iINFO << "SWITCHING ON           "
               << switchingDist << "\n";
      iout << iINFO << "SWITCHING OFF          "
         	    << cutoff << "\n";
      iout << iINFO << "E-SWITCHING ON         "
               << elecswitchDist << "\n";
      iout << iINFO << "E-SWITCHING OFF        "
               << eleccutoff << "\n";
      iout << iINFO << "VDW-SWITCHING ON       "
               << vdwswitchDist << "\n";
      iout << iINFO << "VDW-SWITCHING OFF      "
               << vdwcutoff << "\n";
      iout << iINFO << "PAIRLIST DISTANCE      "
               << pairlistDist << "\n";
   }
   else
   {
      iout << iINFO << "CUTOFF                 " 
   	    << cutoff << "\n";
   }

   if (plMarginCheckOn)
     iout << iINFO << "PAIRLIST CHECK ON\n";
   else 
     iout << iINFO << "PAIRLIST CHECK OFF\n";
            
   iout << iINFO << "MARGIN                 " 
   	 << margin << "\n";
   
   iout << iINFO << "PATCH DIMENSION        "
            << patchDimension << "\n";

   if (outputEnergies != 1)
   {
      iout << iINFO << "ENERGY OUTPUT STEPS    "
   	    << outputEnergies << "\n";
   }
   
   if (constraintsOn)
   {
      iout << iINFO << "HARMONIC CONSTRAINTS ACTIVE\n";

      iout << iINFO << "HARMONIC CONS EXP      "
   	    << constraintExp << "\n";
   }

   if (globalOn && ! dihedralOn)
   {
      iout << iINFO << "GLOBAL INTEGRATION TEST MODE ACTIVE\n";
   }

   if (dihedralOn)
   {
      iout << iINFO << "DIHEDRAL ANGLE DYNAMICS ACTIVE\n";
      if (!COLDOn)
      {
         iout << iINFO << "*** DIHEDRAL ANGLE DYNAMICS IS HIGHLY EXPERIMENTAL ***\n";
         iout << iINFO << "PLEASE CONSIDER USING THE COLD OPTION AS WELL\n";
      }
   }

   if (COLDOn)
   {
      iout << iINFO << "COLD (CONSTRAINED OVERDAMPED LANGEVIN DYNAMICS) ACTIVE\n";

      iout << iINFO << "COLD TARGET TEMP       "
   	    << COLDTemp << "\n";

      iout << iINFO << "COLD COLLISION RATE    "
   	    << COLDRate << "\n";
   }

   if (cylindricalBCOn)
   {
    iout << iINFO << "CYLINDRICAL BOUNDARY CONDITIONS ACTIVE\n";
    iout << iINFO << "AXIS                     " << cylindricalBCAxis << "\n";
    iout << iINFO << "RADIUS #1                " << cylindricalBCr1 << "\n";
    iout << iINFO << "FORCE CONSTANT #1        " << cylindricalBCk1 << "\n";
    iout << iINFO << "EXPONENT #1              " << cylindricalBCexp1 << "\n";
    iout << iINFO << "LENGTH #1                " << cylindricalBCl1 << "\n";
    if (cylindricalBCr2 > 0.0)
    {
     iout << iINFO << "RADIUS #2               " << cylindricalBCr2 << "\n";
     iout << iINFO << "FORCE CONSTANT #2       " << cylindricalBCk2 << "\n";
     iout << iINFO << "EXPONENT #2             " << cylindricalBCexp2 << "\n";
     iout << iINFO << "LENGTH #2               " << cylindricalBCl2 << "\n";
    }
    if (cylindricalCenterCOM)
    {
      iout << iINFO << "CYLINDRICAL BOUNDARIES CENTERED AROUND COM\n";
    }
    else
    {
    iout << iINFO << "CYLINDER BOUNDARY CENTER(" << cylindricalCenter.x << ", "
             << cylindricalCenter.y << ", " << cylindricalCenter.z << ")\n";
    }
  }

   if (sphericalBCOn)
   {
      iout << iINFO << "SPHERICAL BOUNDARY CONDITIONS ACTIVE\n";

      iout << iINFO << "RADIUS #1              "
   	    << sphericalBCr1 << "\n";
      iout << iINFO << "FORCE CONSTANT #1      "
   	    << sphericalBCk1 << "\n";
      iout << iINFO << "EXPONENT #1            "
   	    << sphericalBCexp1 << "\n";

      if (sphericalBCr2 > 0)
      {
      	iout << iINFO << "RADIUS #2              "
   	         << sphericalBCr2 << "\n";
      	iout << iINFO << "FORCE CONSTANT #2      "
   	    	 << sphericalBCk2 << "\n";
      	iout << iINFO << "EXPONENT #2            "
   		 << sphericalBCexp2 << "\n";
      }

      if (sphericalCenterCOM)
      {
	iout << iINFO << "SPHERICAL BOUNDARIES CENTERED AROUND COM\n";
      }
      else
      {
	iout << iINFO << "SPHERE BOUNDARY CENTER(" << sphericalCenter.x << ", "
		 << sphericalCenter.y << ", " << sphericalCenter.z << ")\n";
      }
   }
   
   if (eFieldOn)
   {
      iout << iINFO << "ELECTRIC FIELD ACTIVE\n";
      
      iout << iINFO << "E-FIELD VECTOR         ("
	       << eField.x << ", " << eField.y
	       << ", " << eField.z << ")\n";
   }

   if (langevinOn)
   {
      iout << iINFO << "LANGEVIN DYNAMICS ACTIVE\n";
      iout << iINFO << "LANGEVIN TEMPERATURE   "
   	    << langevinTemp << "\n";
   }

   if (tCoupleOn)
   {
      iout << iINFO << "TEMPERATURE COUPLING ACTIVE\n";
      iout << iINFO << "COUPLING TEMPERATURE   "
   	    << tCoupleTemp << "\n";
   }

   if (minimizeOn)
   {
      iout << iINFO << "MINIMIZATION ACTIVE\n";

      iout << iINFO << "MAXIMUM MOVEMENT       "
   	    << maximumMove << "\n";
   }

   if (rescaleFreq > 0)
   {
   	iout << iINFO << "VELOCITY RESCALE FREQ  "
   		 << rescaleFreq << "\n";
   	iout << iINFO << "VELOCITY RESCALE TEMP  "
   		 << rescaleTemp << "\n";
   }

   if (berendsenPressureOn)
   {
   	iout << iINFO << "BERENDSEN PRESSURE COUPLING ACTIVE\n";
   	iout << iINFO << "    TARGET PRESSURE IS "
   		 << berendsenPressureTarget << " BAR\n";
   	iout << iINFO << "    COMPRESSIBILITY ESTIMATE IS "
   		 << berendsenPressureCompressibility << " BAR^(-1)\n";
   	iout << iINFO << "    RELAXATION TIME IS "
   		 << berendsenPressureRelaxationTime << " FS\n";
	berendsenPressureTarget /= PRESSUREFACTOR;
	berendsenPressureCompressibility *= PRESSUREFACTOR;
   	iout << iINFO << "    APPLIED EVERY "
   		 << berendsenPressureFreq << " STEPS\n";
   }

   if (vmdFrequency > 0)
   {
   	iout << iINFO << "VMD INTERFACE ON\n"
   		 << "VMD FRREQUENCY    "
   		 << vmdFrequency << "\n";
   }

   if (FMAOn)
   {
   	iout << iINFO << "FMA ACTIVE\n";
   	iout << iINFO << "FMA EXECUTION FREQ     "
   		 << fmaFrequency << "\n";
   	iout << iINFO << "FMA THETA              "
   		 << fmaTheta << "\n";
   }

   if (fullDirectOn)
   {
	iout << iINFO << "DIRECT FULL ELECTROSTATIC CALCULATIONS ACTIVE\n";
   }

   if (MTSAlgorithm != NAIVE)
   {
	if (MTSAlgorithm == VERLETI)
	{
		iout << iINFO << "VERLET I MTS SCHEME\n";
	}
	else if (MTSAlgorithm == VERLETII )
        {
		iout << iINFO << "VERLET II MTS SCHEME\n";
        }
        else
	{
		iout << iINFO << "VERLET X MTS SCHEME\n";
	}
   }

   if (longSplitting == SHARP)
	iout << iINFO << "SHARP SPLITTING OF LONG RANGE ELECTROSTATICS\n";
   else if (longSplitting == XPLOR)
	iout << iINFO << "XPLOR SPLITTING OF LONG RANGE ELECTROSTATICS\n";
   else if (longSplitting == C1)
	iout << iINFO << "C1 SPLITTING OF LONG RANGE ELECTROSTATICS\n";

   if (splitPatch == SPLIT_PATCH_POSITION)
	iout << iINFO << "PLACING ATOMS IN PATCHES BY POSITION\n";
   else if (splitPatch == SPLIT_PATCH_HYDROGEN)
	iout << iINFO << "PLACING ATOMS IN PATCHES BY HYDROGEN GROUPS\n";

   if (rigidBonds == RIGID_ALL)
   {
     iout << iINFO <<"RIGID BONDS TO HYDROGEN : ALL, TOLERANCE=" << rigidTol << "\n";
   }
   else if (rigidBonds == RIGID_WATER)
   {
    iout << iINFO<<"RIGID BONDS TO HYDROGEN :  WATER, TOLERANCE="<< rigidTol << "\n";
   }
   

   if (nonbondedFrequency != 1)
   {
     iout << iINFO << "NONBONDED FORCES EVALUATED EVERY " << nonbondedFrequency << " STEPS\n";
   }

   iout << iINFO << "RANDOM NUMBER SEED     "
   	 << randomSeed << "\n";


   iout << iINFO << "USE HYDROGEN BONDS?    ";
   if (HydrogenBonds)
   {
	iout << "YES\n" << endi;
	iout << iINFO << "USE ANTECEDENT ATOMS?  ";
	iout << (useAntecedent ? "YES" : "NO");
        iout << "\nHB DIST CUT, ON, OFF   ";
	iout << daCutoffDist << " , " << daOnDist << " , " << daOffDist;
        iout << "\nHB ANGLE CUT, ON, OFF  ";
	iout << dhaCutoffAngle << " , " << dhaOnAngle << " , ";
	iout << dhaOffAngle;
        iout << "\nHB ATT, REP exponents  ";
	iout << distAttExp << " , " << distRepExp;
        iout << "\nHB AA, HA exponents    ";
	iout << aaAngleExp << " , " << haAngleExp;
	iout << "\n" << endi;
   }
   else
   {
	iout << "NO\n" << endi;
   }


   iout << iINFO << "Here we go config->find\n" << endi;
   current = config->find("coordinates");
   iout << iINFO << "Here done config->find\n" << endi;

   if ( (cwd == NULL) || (current->data[0] == '/') )
   {
        iout << iINFO << "Here cwd==NULL and current is "
	<< current->data << '\n' << endi;
   	strcpy(filename, current->data);
   }
   else
   {
        iout << iINFO << "cwd != NULL and not abs\n" << endi;
   	strcpy(filename, cwd);
   	strcat(filename, current->data);
   }


   iout << iINFO << "COORDINATE PDB         " << filename << '\n' << endi;

   if (opts.defined("bincoordinates"))
   {
	current = config->find("bincoordinates");

   	if ( (cwd == NULL) || (current->data[0] == '/') )
   	{
   		strcpy(filename, current->data);
   	}
   	else
   	{
   		strcpy(filename, cwd);
   		strcat(filename, current->data);
   	}

   	iout << iINFO << "BINARY COORDINATES     " 
   	         << filename << "\n";
   }

   current = config->find("structure");

   if ( (cwd == NULL) || (current->data[0] == '/') )
   {
   	strcpy(filename, current->data);
   }
   else
   {
   	strcpy(filename, cwd);
   	strcat(filename, current->data);
   }

   iout << iINFO << "STRUCTURE FILE         " 
   	 << filename << "\n" << endi;

   current = config->find("parameters");

   while (current != NULL)
   {
   	if ( (cwd == NULL) || (current->data[0] == '/') )
   	{
   		strcpy(filename, current->data);
   	}
   	else
   	{
   		strcpy(filename, cwd);
   		strcat(filename, current->data);
   	}

   	iout << iINFO << "PARAMETERS             " 
   		 << filename << "\n" << endi;
   	current = current->next;
   }

   if (firstTimestep)
   {
	iout << iINFO << "FIRST TIMESTEP         "
		 << firstTimestep << "\n" << endi;
   }

}
/*		END OF FUNCTION initialize_config_data		*/

/****************************************************************/
/*								*/
/*		FUNCTION send_SimParameters			*/
/*								*/
/*	This function is used by the master process to broadcast*/
/*  the parameter data to all the other nodes.  It just builds  */
/*  a message with all the relevant data and broadcasts it to   */
/*  the other nodes.  The routine receive_SimParameters is used */
/*  by all the other nodes to receive this message.		*/
/*								*/
/****************************************************************/

void SimParameters::send_SimParameters(Communicate *com_obj)

{
	Message *msg = new Message;	//  Message to send
	if ( msg == NULL )
	{
	  NAMD_die("memory allocation failed in SimParameters::send_SimParameters");
	}

	msg->put(dt).put(N).put(stepsPerCycle);
	msg->put(ldbStrategy).put(ldbPeriod).put(firstLdbStep);
	msg->put(initialTemp).put(comMove);
	msg->put(dielectric).put(exclude).put(scale14);
	msg->put(dcdFrequency).put(velDcdFrequency).put(vmdFrequency);
	msg->put(dcdFilename).put(velDcdFilename).put(outputFilename);
	msg->put(restartFilename).put(restartFrequency).put(cutoff);
	msg->put(eleccutoff).put(vdwcutoff);
	msg->put(margin).put(patchDimension).put(switchingActive);
	msg->put(switchingDist).put(elecswitchDist).put(vdwswitchDist);
	msg->put(pairlistDist).put(plMarginCheckOn).put(constraintsOn);
	msg->put(constraintExp).put(FMAOn).put(FMALevels).put(FMAMp);
	msg->put(FMAFFTOn).put(FMAFFTBlock).put(minimizeOn);
	msg->put(maximumMove).put(totalAtoms).put(randomSeed);
	msg->put(langevinOn).put(langevinTemp).put(globalOn);
	msg->put(dihedralOn).put(COLDOn).put(COLDRate).put(COLDTemp);
	msg->put(rescaleFreq).put(rescaleTemp);
	msg->put(sphericalBCOn).put(sphericalBCr1);
	msg->put(sphericalBCr2).put(sphericalBCk1).put(sphericalBCk2);
	msg->put(sphericalBCexp1).put(sphericalBCexp2);
	msg->put(firstTimestep).put(fullDirectOn);
	msg->put(eFieldOn).put(&eField).put(binaryRestart);
	msg->put(electForceDcdFilename).put(electForceDcdFrequency);
	msg->put(allForceDcdFilename).put(allForceDcdFrequency);
	msg->put(MTSAlgorithm).put(sphericalCenterCOM).put(&sphericalCenter);
	msg->put(longSplitting).put(tCoupleOn).put(tCoupleTemp);
	msg->put(splitPatch);
	msg->put(fmaFrequency).put(fmaTheta);
        msg->put(rigidBonds);
        msg->put(rigidTol);
	msg->put(nonbondedFrequency);

	// send hydrogen bond data
	msg->put(HydrogenBonds).put(useAntecedent);
	msg->put(aaAngleExp).put(haAngleExp).put(distAttExp).put(distRepExp);
	msg->put(dhaCutoffAngle).put(dhaOnAngle).put(dhaOffAngle);
	msg->put(daCutoffDist).put(daOnDist).put(daOffDist);

	// send cylindrical boundary conditions data
        msg->put(cylindricalBCOn).put(cylindricalBCr1);
        msg->put(cylindricalBCr2).put(cylindricalBCk1);
        msg->put(cylindricalBCk2).put(cylindricalBCl1);
        msg->put(cylindricalBCl2);
        msg->put(&cylindricalCenter);
        msg->put(cylindricalBCexp1).put(cylindricalBCexp2);
        msg->put(cylindricalBCAxis);

	// send periodic box data
	msg->put(cellBasisVector1.x);
	msg->put(cellBasisVector2.y);
	msg->put(cellBasisVector3.z);
	msg->put(&cellOrigin);

	// send pressure control data
	msg->put(berendsenPressureOn);
	msg->put(berendsenPressureTarget);
	msg->put(berendsenPressureCompressibility);
	msg->put(berendsenPressureRelaxationTime);
	msg->put(berendsenPressureFreq);

	// now broadcast this info to all other nodes
	com_obj->broadcast_others(msg, SIMPARAMSTAG);
}
/*		END OF FUNCITON send_SimParameters		*/

/****************************************************************/
/*								*/
/*			FUNCTION receive_SimParameters		*/
/*								*/
/*	This function is used by all the child nodes to 	*/
/*  receive the simulation parameters from the master node.	*/
/*								*/
/****************************************************************/

void SimParameters::receive_SimParameters(Message *msg)

{
	//  Get each of the parameters from the message
	msg->get(dt);
	msg->get(N);
	msg->get(stepsPerCycle);
	msg->get(ldbStrategy);
	msg->get(ldbPeriod);
	msg->get(firstLdbStep);
	msg->get(initialTemp);
	msg->get(comMove);
	msg->get(dielectric);
	msg->get(exclude);
	msg->get(scale14);
	msg->get(dcdFrequency);
	msg->get(velDcdFrequency);
	msg->get(vmdFrequency);
	msg->get(dcdFilename);
	msg->get(velDcdFilename);
	msg->get(outputFilename);
	msg->get(restartFilename);
	msg->get(restartFrequency);
	msg->get(cutoff);
	msg->get(eleccutoff);
	msg->get(vdwcutoff);
	msg->get(margin);
	msg->get(patchDimension);
	msg->get(switchingActive);
	msg->get(switchingDist);
	msg->get(elecswitchDist);
	msg->get(vdwswitchDist);
	msg->get(pairlistDist);
	msg->get(plMarginCheckOn);
	msg->get(constraintsOn);
	msg->get(constraintExp);
	msg->get(FMAOn);
	msg->get(FMALevels);
	msg->get(FMAMp);
	msg->get(FMAFFTOn);
	msg->get(FMAFFTBlock);
	msg->get(minimizeOn);
	msg->get(maximumMove);
	msg->get(totalAtoms);
	msg->get(randomSeed);
	msg->get(langevinOn);
	msg->get(langevinTemp);
	msg->get(globalOn);
	msg->get(dihedralOn);
	msg->get(COLDOn);
	msg->get(COLDRate);
	msg->get(COLDTemp);
	msg->get(rescaleFreq);
	msg->get(rescaleTemp);
	msg->get(sphericalBCOn);
	msg->get(sphericalBCr1);
	msg->get(sphericalBCr2);
	msg->get(sphericalBCk1);
	msg->get(sphericalBCk2);
	msg->get(sphericalBCexp1);
	msg->get(sphericalBCexp2);
	msg->get(firstTimestep);
	msg->get(fullDirectOn);
	msg->get(eFieldOn);
	msg->get(&eField);
	msg->get(binaryRestart);
	msg->get(electForceDcdFilename);
	msg->get(electForceDcdFrequency);
	msg->get(allForceDcdFilename);
	msg->get(allForceDcdFrequency);
	msg->get(MTSAlgorithm);
	msg->get(sphericalCenterCOM);
	msg->get(&sphericalCenter);
	msg->get(longSplitting);
	msg->get(tCoupleOn);
	msg->get(tCoupleTemp);
	msg->get(splitPatch);
	msg->get(fmaFrequency);
	msg->get(fmaTheta);
	msg->get(rigidBonds);
	msg->get(rigidTol);
	msg->get(nonbondedFrequency);

	// receive hydrogen bond data
	msg->get(HydrogenBonds);
	msg->get(useAntecedent);
	msg->get(aaAngleExp);
	msg->get(haAngleExp);
	msg->get(distAttExp);
	msg->get(distRepExp);
	msg->get(dhaCutoffAngle);
	msg->get(dhaOnAngle);
	msg->get(dhaOffAngle);
	msg->get(daCutoffDist);
	msg->get(daOnDist);
	msg->get(daOffDist);

	// receive cylindrical boundary conditions data
        msg->get(cylindricalBCOn);
        msg->get(cylindricalBCr1);
        msg->get(cylindricalBCr2);
        msg->get(cylindricalBCk1);
        msg->get(cylindricalBCk2);
        msg->get(cylindricalBCl1);
        msg->get(cylindricalBCl2);
        msg->get(&cylindricalCenter);
        msg->get(cylindricalBCexp1);
        msg->get(cylindricalBCexp2);
        msg->get(cylindricalBCAxis);

	// receive periodic box data
	msg->get(cellBasisVector1.x);
	msg->get(cellBasisVector2.y);
	msg->get(cellBasisVector3.z);
	msg->get(&cellOrigin);
	lattice.set(cellBasisVector1,cellBasisVector2,cellBasisVector3,cellOrigin);

	// receive pressure control data
	msg->get(berendsenPressureOn);
	msg->get(berendsenPressureTarget);
	msg->get(berendsenPressureCompressibility);
	msg->get(berendsenPressureRelaxationTime);
	msg->get(berendsenPressureFreq);

	//  Free the message
	delete msg;

}
/*			END OF FUNCTION receive_SimParameters	*/


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1017 $	$Date: 1997/05/29 19:12:11 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: SimParameters.C,v $
 * Revision 1.1017  1997/05/29 19:12:11  nealk
 * Modified so hydrogen grouping margin offset is a user-defined option.
 *
 * Revision 1.1016  1997/05/09 18:24:24  nealk
 * 1. Added hydrogen grouping code to improve performance in ComputeNonbondedBase
 *    CODE ONLY WORKS WITH HYDROGEN GROUPING!
 * 2. Increased the hydrogen group cutoff side from 2A to 2.5A -- 2A gave
 *    fractionally different values after 100 iterations.  2.5A gives same numbers.
 * 3. Made migration by hydrogen grouping the default in SimParameters.
 *
 * Revision 1.1015  1997/04/24 18:51:45  nealk
 * Corrected parameter bug with DPMTA options.
 * In particular: FFT Block must be 4 and FFTMp must be a multiple
 * of FFT Block.
 *
 * Revision 1.1014  1997/04/16 23:44:04  brunner
 * Put ldbStrategy={none|refineonly|alg7}, ldbPeriod, and firstLdbStep
 * in SimParameters.
 *
 * Revision 1.1013  1997/04/16 22:12:20  brunner
 * Fixed an LdbCoordinator bug, and cleaned up timing and Ldb output some.
 *
 * Revision 1.1012  1997/04/08 21:08:49  jim
 * Contant pressure now correct on multiple nodes, should work with MTS.
 *
 * Revision 1.1011  1997/03/31 19:23:04  nealk
 * Hard-coded +2 angstrom increase in patch margin when using hydrogen grouping.
 *
 * Revision 1.1010  1997/03/27 17:08:31  nealk
 * Added hydrogen groupings.  Now configuration parameter "splitPatch" determines
 * atom-into-patch distribution.
 *
 * Revision 1.1009  1997/03/27 08:04:24  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1008  1997/03/27 03:16:57  jim
 * Added code to check virial calculation, fixed problems with DPMTA and PBC's.
 *
 * Revision 1.1007  1997/03/25 23:01:03  jim
 * Added nonbondedFrequency parameter and multiple time-stepping
 *
 * Revision 1.1006  1997/03/21 23:05:45  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1005  1997/03/19 18:10:17  nealk
 * Added sorted hydrogen group list to molecule.
 *
 * Revision 1.1004  1997/03/19 11:54:57  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
