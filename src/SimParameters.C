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
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/06 20:38:38 $
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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/SimParameters.C,v 1.1 1996/08/06 20:38:38 ari Exp $";

#include "ConfigList.h"
#include "SimParameters.h"
#include "ParseOptions.h"
#include "structures.h"
#include "Communicate.h"
#include "Message.h"
#include <stdio.h>
#include "Inform.h"
#include <time.h>

#ifdef AIX
#include "strlib.h"		//  For strcasecmp and strncasecmp
#endif

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


   opts.optional("main", "rigidBonds", "Rigid bonds to hydrogen",PARSE_STRING);
   opts.optional("main", "rigidTolerance", 
                  "Error tolerance for rigid bonds to hydrogen",&rigidTol);
   
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

   ///////////////  Electric field options
   opts.optionalB("main", "eFieldOn", "Should and electric field be applied",
                 &eFieldOn, FALSE);
   opts.require("eFieldOn", "eField", "Electric field vector", &eField);
   
   ///////////////  Load balance options
   opts.optional("main", "ldbstrategy", "Load balancing strategy",
		 loadStrategy);
   opts.optional("ldbstrategy", "ldbstepspercycle",
		 "steps between load balancing", &ldbStepsPerCycle);
   opts.range("ldbstepspercycle", POSITIVE);
   opts.optional("ldbstrategy", "ldbsendstep", "when to send load stats",
		 &ldbSendStep);
   opts.range("ldbsendstep", POSITIVE);

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

   ///// exclude stuff
   char s[129];
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
      namdWarn << "Value given for '1-4scaling', but 1-4 scaling "
	      << "not in effect.\n This value will be ignored!" << sendmsg;
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

   //  Get the long range force splitting specification
   if (!opts.defined("longSplitting"))
   {
	longSplitting = SHARP;
   }
   else
   {
	opts.get("longSplitting", s);
	if (!strcasecmp(s, "sharp"))
	{
		longSplitting = SHARP;
	}
	else if (!strcasecmp(s, "xplor"))
	{
		longSplitting = XPLOR;
	}
	else if (!strcasecmp(s, "c1"))
	{
		longSplitting = C1;
	}
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
         namdWarn << "rigidBonds is on but rigidTol is not provided:" ;
         namdWarn << " using default value" << sendmsg;
      }
   }
   
   //  Take care of switching stuff
   if (switchingActive)
   {

     if (opts.defined("cutoff") && opts.defined("eleccutoff") &&
	 opts.defined("vdwcutoff")) {
       namdWarn
	 << "PARAMETERS cutoff, eleccutoff, and vdwcutoff\n"
	 << "all defined -- cutoff ignored\n" 
	 << sendmsg;
     }

     if (opts.defined("switchDist") && opts.defined("elecswitchDist") &&
	 opts.defined("vdwswitchDist")) {
       namdWarn
	 << "PARAMETERS switchDist, elecswitchDist, and vdwswitchDist\n"
	 << "all defined -- switchDist ignored\n" 
	 << sendmsg;
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
   if (opts.defined("ldbstrategy"))
   {
	//  Set default values
	if (!opts.defined("ldbsendstep"))
	{
		ldbSendStep=ldbStepsPerCycle-1;
	}

	if (!opts.defined("ldbStepsPerCycle"))
	{
		ldbStepsPerCycle=4*stepsPerCycle;
	}

	//  Check to make things look OK
	if (ldbStepsPerCycle % stepsPerCycle != 0)
	{
		NAMD_die("Number of steps per load balance cycle must be a multiple of stepsPerCycle");
	}

	if (ldbStepsPerCycle/stepsPerCycle < 2)
	{
		NAMD_die("Number of steps per load balance cycle must be at least 2 * stepsPerCycle");
	}

	if (ldbSendStep >= ldbStepsPerCycle )
	{
		NAMD_die("ldbSendStep must be less than ldbStepsPerCycle");
	}

	//  Assign the load balancing strategy
	if (strcasecmp(loadStrategy, "none") == 0)
	{
	    ldbStrategy=LDBSTRAT_NONE;
	}
	else if (strcasecmp(loadStrategy, "random") == 0)
	{
	    ldbStrategy=LDBSTRAT_RANDOM;
	}
	else if (strcasecmp(loadStrategy, "nolocality") == 0)
	{
	    ldbStrategy=LDBSTRAT_NOLOCAL;
	}
	else if (strcasecmp(loadStrategy, "bisection") == 0)
	{
	    ldbStrategy=LDBSTRAT_RBISEC;
	}
	else if (strcasecmp(loadStrategy, "other") == 0)
	{
	    ldbStrategy=LDBSTRAT_OTHER;
	}
	else
	{
	    NAMD_die("Unknown ldbStrategy selected");
	}
   }
   else
   {
	ldbStrategy=LDBSTRAT_NONE;
	ldbStepsPerCycle=4*stepsPerCycle;
	ldbSendStep=ldbStepsPerCycle-1;
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
   namdInfo << "SIMULATION PARAMETERS:\n";
   namdInfo << "TIMESTEP               " << dt << sendmsg;
   namdInfo << "NUMBER OF STEPS        " << N << "\n";
   namdInfo << "STEPS PER CYCLE        " << stepsPerCycle << "\n";
   if (ldbStrategy==LDBSTRAT_NONE)  {
     namdInfo << "LOAD BALANCE STRATEGY  none\n";
   } else if (ldbStrategy==LDBSTRAT_RANDOM)  {
     namdInfo << "LOAD BALANCE STRATEGY  random\n";
     namdInfo << "STEPS PER LDB CYCLE    " << ldbStepsPerCycle << "\n";
     namdInfo << "SEND LDB TIMESTEP      " << ldbSendStep << sendmsg;
   } else if (ldbStrategy==LDBSTRAT_NOLOCAL)  {
     namdInfo << "LOAD BALANCE STRATEGY  nolocality\n";
     namdInfo << "STEPS PER LDB CYCLE    " << ldbStepsPerCycle << "\n";
     namdInfo << "SEND LDB TIMESTEP      " << ldbSendStep << sendmsg;
   } else if (ldbStrategy==LDBSTRAT_RBISEC)  {
     namdInfo << "LOAD BALANCE STRATEGY  bisection\n";
     namdInfo << "STEPS PER LDB CYCLE    " << ldbStepsPerCycle << "\n";
     namdInfo << "SEND LDB TIMESTEP      " << ldbSendStep << "\n";
   } else if (ldbStrategy==LDBSTRAT_OTHER)  {
     namdInfo << "LOAD BALANCE STRATEGY  other\n";
     namdInfo << "STEPS PER LDB CYCLE    " << ldbStepsPerCycle << "\n";
     namdInfo << "SEND LDB TIMESTEP      " << ldbSendStep << sendmsg;
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

	namdInfo << "VELOCITY FILE          " << filename << sendmsg;
   }
   else
   {
	namdInfo << "INITIAL TEMPERATURE    " 
		 << initialTemp << "\n";
   }

   namdInfo << "CENTER OF MASS MOVING? ";

   if (comMove)
   {
   	namdInfo << "YES\n";
   }
   else
   {
   	namdInfo << "NO\n";
   }

   namdInfo << "DIELECTRIC             " 
   	 << dielectric << "\n";

   namdInfo << "EXCLUDE                ";

   switch (exclude)
   {
   	case NONE:
   		namdInfo << "NONE\n";
   		break;
   	case ONETWO:
   		namdInfo << "ONETWO\n";
   		break;
   	case ONETHREE:
   		namdInfo << "ONETHREE\n";
   		break;
   	case ONEFOUR:
   		namdInfo << "ONE-FOUR\n";
   		break;
   	default:
   		namdInfo << "SCALED ONE-FOUR\n";
   		break;
   }

   if (exclude == SCALED14)
   {
   	namdInfo << "1-4 SCALE FACTOR       " 
   		 << scale14 << sendmsg;
   }

   if (dcdFrequency > 0)
   {
   	namdInfo << "DCD FILENAME           " 
   		 << dcdFilename << "\n";
   	namdInfo << "DCD FREQUENCY          " 
   		 << dcdFrequency << "\n";
   }
   else
   {
   	namdInfo << "NO DCD TRAJECTORY OUTPUT\n";
   }
   
   if (velDcdFrequency > 0)
   {
   	namdInfo << "VELOCITY DCD FILENAME  " 
   		 << velDcdFilename << "\n";
   	namdInfo << "VELOCITY DCD FREQUENCY " 
   		 << velDcdFrequency << "\n";
   }
   else
   {
   	namdInfo << "NO VELOCITY DCD OUTPUT\n";
   }
   
   if (electForceDcdFrequency > 0)
   {
   	namdInfo << "ELECT FORCE DCD NAME   " 
   		 << electForceDcdFilename << "\n";
   	namdInfo << "ELECT FORCE DCD FREQ   " 
   		 << electForceDcdFrequency << sendmsg;
   }

   if (allForceDcdFrequency > 0)
   {
   	namdInfo << "TOTAL FORCE DCD NAME   " 
   		 << allForceDcdFilename << "\n";
   	namdInfo << "TOTAL FORCE DCD FREQ   " 
   		 << allForceDcdFrequency << "\n";
   }
   	
   namdInfo << "OUTPUT FILENAME        " 
   	 << outputFilename << "\n";

   if (restartFrequency == -1)
   {
   	namdInfo << "NO RESTART FILE\n";
   }
   else
   {
   	namdInfo << "RESTART FILENAME       "
   		 << restartFilename << "\n";
   	namdInfo << "RESTART FREQUENCY      " 
   		 << restartFrequency << "\n";

	if (binaryRestart)
	{
		namdInfo << "BINARY RESTART FILES WILL BE USED\n";
	}
   }
   
   if (switchingActive)
   {
      namdInfo << "SWITCHING ACTIVE\n";
      namdInfo << "SWITCHING ON           "
               << switchingDist << "\n";
      namdInfo << "SWITCHING OFF          "
         	    << cutoff << "\n";
      namdInfo << "E-SWITCHING ON         "
               << elecswitchDist << "\n";
      namdInfo << "E-SWITCHING OFF        "
               << eleccutoff << "\n";
      namdInfo << "VDW-SWITCHING ON       "
               << vdwswitchDist << "\n";
      namdInfo << "VDW-SWITCHING OFF      "
               << vdwcutoff << "\n";
      namdInfo << "PAIRLIST DISTANCE      "
               << pairlistDist << "\n";
   }
   else
   {
      namdInfo << "CUTOFF                 " 
   	    << cutoff << "\n";
   }

   if (plMarginCheckOn)
     namdInfo << "PAIRLIST CHECK ON\n";
   else 
     namdInfo << "PAIRLIST CHECK OFF\n";
            
   namdInfo << "MARGIN                 " 
   	 << margin << "\n";
   
   namdInfo << "PATCH DIMENSION        "
            << patchDimension << "\n";

   if (outputEnergies != 1)
   {
      namdInfo << "ENERGY OUTPUT STEPS    "
   	    << outputEnergies << "\n";
   }
   
   if (constraintsOn)
   {
      namdInfo << "HARMONIC CONSTRAINTS ACTIVE\n";

      namdInfo << "HARMONIC CONS EXP      "
   	    << constraintExp << "\n";
   }

   if (globalOn && ! dihedralOn)
   {
      namdInfo << "GLOBAL INTEGRATION TEST MODE ACTIVE\n";
   }

   if (dihedralOn)
   {
      namdInfo << "DIHEDRAL ANGLE DYNAMICS ACTIVE\n";
      if (!COLDOn)
      {
         namdInfo << "*** DIHEDRAL ANGLE DYNAMICS IS HIGHLY EXPERIMENTAL ***\n";
         namdInfo << "PLEASE CONSIDER USING THE COLD OPTION AS WELL\n";
      }
   }

   if (COLDOn)
   {
      namdInfo << "COLD (CONSTRAINED OVERDAMPED LANGEVIN DYNAMICS) ACTIVE\n";

      namdInfo << "COLD TARGET TEMP       "
   	    << COLDTemp << "\n";

      namdInfo << "COLD COLLISION RATE    "
   	    << COLDRate << "\n";
   }

   if (sphericalBCOn)
   {
      namdInfo << "SPHERICAL BOUNDARY CONDITIONS ACTIVE\n";

      namdInfo << "RADIUS #1              "
   	    << sphericalBCr1 << "\n";
      namdInfo << "FORCE CONSTANT #1      "
   	    << sphericalBCk1 << "\n";
      namdInfo << "EXPONENT #1            "
   	    << sphericalBCexp1 << "\n";

      if (sphericalBCr2 > 0)
      {
      	namdInfo << "RADIUS #2              "
   	         << sphericalBCr2 << "\n";
      	namdInfo << "FORCE CONSTANT #2      "
   	    	 << sphericalBCk2 << "\n";
      	namdInfo << "EXPONENT #2            "
   		 << sphericalBCexp2 << "\n";
      }

      if (sphericalCenterCOM)
      {
	namdInfo << "SPHERICAL BOUNDARIES CENTERED AROUND COM\n";
      }
      else
      {
	namdInfo << "SPHERE BOUNDARY CENTER(" << sphericalCenter.x << ", "
		 << sphericalCenter.y << ", " << sphericalCenter.z << ")\n";
      }
   }
   
   if (eFieldOn)
   {
      namdInfo << "ELECTRIC FIELD ACTIVE\n";
      
      namdInfo << "E-FIELD VECTOR         ("
	       << eField.x << ", " << eField.y
	       << ", " << eField.z << ")\n";
   }

   if (langevinOn)
   {
      namdInfo << "LANGEVIN DYNAMICS ACTIVE\n";
      namdInfo << "LANGEVIN TEMPERATURE   "
   	    << langevinTemp << "\n";
   }

   if (tCoupleOn)
   {
      namdInfo << "TEMPERATURE COUPLING ACTIVE\n";
      namdInfo << "COUPLING TEMPERATURE   "
   	    << tCoupleTemp << "\n";
   }

   if (minimizeOn)
   {
      namdInfo << "MINIMIZATION ACTIVE\n";

      namdInfo << "MAXIMUM MOVEMENT       "
   	    << maximumMove << "\n";
   }

   if (rescaleFreq > 0)
   {
   	namdInfo << "VELOCITY RESCALE FREQ  "
   		 << rescaleFreq << "\n";
   	namdInfo << "VELOCITY RESCALE TEMP  "
   		 << rescaleTemp << "\n";
   }

   if (vmdFrequency > 0)
   {
   	namdInfo << "VMD INTERFACE ON\n"
   		 << "VMD FRREQUENCY    "
   		 << vmdFrequency << "\n";
   }

   if (FMAOn)
   {
   	namdInfo << "FMA ACTIVE\n";
   	namdInfo << "FMA EXECUTION FREQ     "
   		 << fmaFrequency << "\n";
   	namdInfo << "FMA THETA              "
   		 << fmaTheta << "\n";
   }

   if (fullDirectOn)
   {
	namdInfo << "DIRECT FULL ELECTROSTATIC CALCULATIONS ACTIVE\n";
   }

   if (MTSAlgorithm != NAIVE)
   {
	if (MTSAlgorithm == VERLETI)
	{
		namdInfo << "VERLET I MTS SCHEME\n";
	}
	else if (MTSAlgorithm == VERLETII )
        {
		namdInfo << "VERLET II MTS SCHEME\n";
        }
        else
	{
		namdInfo << "VERLET X MTS SCHEME\n";
	}
   }

   if (longSplitting == SHARP)
   {
	namdInfo << "SHARP SPLITTING OF LONG RANGE ELECTROSTATICS\n";
   }
   else if (longSplitting == XPLOR)
   {
	namdInfo << "XPLOR SPLITTING OF LONG RANGE ELECTROSTATICS\n";
   }
   else if (longSplitting == C1)
   {
	namdInfo << "C1 SPLITTING OF LONG RANGE ELECTROSTATICS\n";
   }

   if (rigidBonds == RIGID_ALL)
   {
     namdInfo <<"RIGID BONDS TO HYDROGEN : ALL, TOLERANCE=" << rigidTol << "\n";
   }
   else if (rigidBonds == RIGID_WATER)
   {
    namdInfo<<"RIGID BONDS TO HYDROGEN :  WATER, TOLERANCE="<< rigidTol << "\n";
   }
   

   namdInfo << "RANDOM NUMBER SEED     "
   	 << randomSeed << "\n";


   namdInfo << "USE HYDROGEN BONDS?    ";
   if (HydrogenBonds)
   {
	namdInfo << "YES" << sendmsg;
	namdInfo << "USE ANTECEDENT ATOMS?  ";
	namdInfo << (useAntecedent ? "YES" : "NO");
        namdInfo << "\nHB DIST CUT, ON, OFF   ";
	namdInfo << daCutoffDist << " , " << daOnDist << " , " << daOffDist;
        namdInfo << "\nHB ANGLE CUT, ON, OFF  ";
	namdInfo << dhaCutoffAngle << " , " << dhaOnAngle << " , ";
	namdInfo << dhaOffAngle;
        namdInfo << "\nHB ATT, REP exponents  ";
	namdInfo << distAttExp << " , " << distRepExp;
        namdInfo << "\nHB AA, HA exponents    ";
	namdInfo << aaAngleExp << " , " << haAngleExp;
	namdInfo << sendmsg;
   }
   else
   {
	namdInfo << "NO" << sendmsg;
   }


   namdInfo << "Here we go config->find " << sendmsg;
   current = config->find("coordinates");
   namdInfo << "Here done config->find " << sendmsg;

   if ( (cwd == NULL) || (current->data[0] == '/') )
   {
        namdInfo << "Here cwd==NULL and current is "
	<< current->data << sendmsg;
   	strcpy(filename, current->data);
   }
   else
   {
        namdInfo << "cwd != NULL and not abs" << sendmsg;
   	strcpy(filename, cwd);
   	strcat(filename, current->data);
   }


   namdInfo << "COORDINATE PDB         " 
   	 << filename << sendmsg;

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

   	namdInfo << "BINARY COORDINATES     " 
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

   namdInfo << "STRUCTURE FILE         " 
   	 << filename << sendmsg;

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

   	namdInfo << "PARAMETERS             " 
   		 << filename << sendmsg;
   	current = current->next;
   }

   if (firstTimestep)
   {
	namdInfo << "FIRST TIMESTEP         "
		 << firstTimestep << sendmsg;
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
	msg->put(ldbStrategy).put(ldbStepsPerCycle).put(ldbSendStep);
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
	msg->put(fmaFrequency).put(fmaTheta);
        msg->put(rigidBonds);
        msg->put(rigidTol);

	// send hydrogen bond data
	msg->put(HydrogenBonds).put(useAntecedent);
	msg->put(aaAngleExp).put(haAngleExp).put(distAttExp).put(distRepExp);
	msg->put(dhaCutoffAngle).put(dhaOnAngle).put(dhaOffAngle);
	msg->put(daCutoffDist).put(daOnDist).put(daOffDist);

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
	msg->get(ldbStepsPerCycle);
	msg->get(ldbSendStep);
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
	msg->get(fmaFrequency);
	msg->get(fmaTheta);
	msg->get(rigidBonds);
	msg->get(rigidTol);

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

	//  Free the message
	delete msg;
}
/*			END OF FUNCTION receive_SimParameters	*/



