/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   SimParameters is just a glorified structure to hold the global
   static simulation parameters such as timestep size, cutoff, etc. that
   are read in from the configuration file.
*/

#include "InfoStream.h"
#include "ComputeNonbondedUtil.h"
#include "ConfigList.h"
#include "SimParameters.h"
#include "ParseOptions.h"
#include "structures.h"
#include "Communicate.h"
#include "MStream.h"
#include <stdio.h>
#include <time.h>
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_NO_TYPE_PREFIX
#include <fftw.h>
#include <rfftw.h>
#else
#include <sfftw.h>
#include <srfftw.h>
#endif
#endif
#if defined(WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define CHDIR _chdir
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#else
#include <unistd.h>
#define CHDIR chdir
#define PATHSEP '/'
#define PATHSEPSTR "/"
#endif
#include <fstream>
using namespace std;

#ifdef WIN32
extern "C" {
  double erfc(double);
}
#endif

#include "strlib.h"    //  For strcasecmp and strncasecmp

// #define DEBUGM
#include "Debug.h"

#define XXXBIGREAL 1.0e32

/************************************************************************/
/*                  */
/*      FUNCTION initialize_config_data      */
/*                  */
/*  This function is used by the master process to populate the     */
/*   simulation parameters from a ConfigList object that is passed to   */
/*   it.  Each parameter is checked to make sure that it has a value    */
/*   that makes sense, and that it doesn't conflict with any other      */
/*   values that have been given.          */
/*                  */
/************************************************************************/

void SimParameters::initialize_config_data(ConfigList *config, char *&cwd)

{

   ParseOptions opts;   //  Object to check consistency of config file

   config_parser(opts);

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

   check_config(opts,config,cwd);

   print_config(opts,config,cwd);

}

/************************************************************************/
/*                                                                      */
/*      FUNCTION scriptSet                                              */
/*                                                                      */
/************************************************************************/

int atobool(const char *s) {
  return ( (! strncasecmp(s,"yes",8)) ||
           (! strncasecmp(s,"on",8)) || 
           (! strncasecmp(s,"true",8)) );
};
                         
void SimParameters::scriptSet(const char *param, const char *value) {

#define MAX_SCRIPT_PARAM_SIZE 128
#define SCRIPT_PARSE_BOOL(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR) = atobool(value); return; } }
#define SCRIPT_PARSE_INT(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR) = atoi(value); return; } }
#define SCRIPT_PARSE_FLOAT(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR) = atof(value); return; } }
#define SCRIPT_PARSE_VECTOR(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR).set(value); return; } }
#define SCRIPT_PARSE_STRING(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { strcpy(VAR,value); return; } }

  SCRIPT_PARSE_FLOAT("scriptArg1",scriptArg1)
  SCRIPT_PARSE_FLOAT("scriptArg2",scriptArg2)
  SCRIPT_PARSE_FLOAT("scriptArg3",scriptArg3)
  SCRIPT_PARSE_FLOAT("scriptArg4",scriptArg4)
  SCRIPT_PARSE_FLOAT("scriptArg5",scriptArg5)
  SCRIPT_PARSE_INT("numsteps",N)
  SCRIPT_PARSE_INT("firsttimestep",firstTimestep)
  SCRIPT_PARSE_FLOAT("reassignTemp",reassignTemp)
  SCRIPT_PARSE_FLOAT("rescaleTemp",rescaleTemp)
  // SCRIPT_PARSE_BOOL("Langevin",langevinOn)
  SCRIPT_PARSE_FLOAT("langevinTemp",langevinTemp)
  SCRIPT_PARSE_FLOAT("initialTemp",initialTemp)
  SCRIPT_PARSE_BOOL("useGroupPressure",useGroupPressure)
  SCRIPT_PARSE_BOOL("useFlexibleCell",useFlexibleCell)
  SCRIPT_PARSE_BOOL("useConstantArea",useConstantArea)
  SCRIPT_PARSE_BOOL("useConstantRatio",useConstantRatio)
  SCRIPT_PARSE_BOOL("LangevinPiston",langevinPistonOn)
  SCRIPT_PARSE_FLOAT("LangevinPistonTarget",langevinPistonTarget)
  SCRIPT_PARSE_FLOAT("LangevinPistonPeriod",langevinPistonPeriod)
  SCRIPT_PARSE_FLOAT("LangevinPistonDecay",langevinPistonDecay)
  SCRIPT_PARSE_FLOAT("LangevinPistonTemp",langevinPistonTemp)
  SCRIPT_PARSE_FLOAT("SurfaceTensionTarget",surfaceTensionTarget)
  SCRIPT_PARSE_BOOL("BerendsenPressure",berendsenPressureOn)
  SCRIPT_PARSE_FLOAT("BerendsenPressureTarget",berendsenPressureTarget)
  SCRIPT_PARSE_FLOAT("BerendsenPressureCompressibility",
				berendsenPressureCompressibility)
  SCRIPT_PARSE_FLOAT("BerendsenPressureRelaxationTime",
				berendsenPressureRelaxationTime)
  SCRIPT_PARSE_FLOAT("constraintScaling",constraintScaling)
  SCRIPT_PARSE_STRING("outputname",outputFilename)
  SCRIPT_PARSE_STRING("tclBCArgs",tclBCArgs)
  SCRIPT_PARSE_VECTOR("eField",eField)
  SCRIPT_PARSE_VECTOR("stirAxis",stirAxis)
  SCRIPT_PARSE_VECTOR("stirPivot",stirPivot)

  if ( ! strncasecmp(param,"fixedatoms",MAX_SCRIPT_PARAM_SIZE) ) {
    if ( ! fixedAtomsOn )
      NAMD_die("FixedAtoms may not be enabled in a script.");
    if ( ! fixedAtomsForces )
      NAMD_die("To use fixedAtoms in script first use fixedAtomsForces yes.");
    fixedAtomsOn = atobool(value);
    return;
  }

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
  SCRIPT_PARSE_INT("fepEquilSteps",fepEquilSteps)

  if ( ! strncasecmp(param,"lambda",MAX_SCRIPT_PARAM_SIZE) ) {
    lambda = atof(value);
    ComputeNonbondedUtil::select();
    return;
  }

  if ( ! strncasecmp(param,"lambda2",MAX_SCRIPT_PARAM_SIZE) ) {
    lambda2 = atof(value);
    ComputeNonbondedUtil::select();
    return;
  }
//fepe

  if ( ! strncasecmp(param,"nonbondedScaling",MAX_SCRIPT_PARAM_SIZE) ) {
    nonbondedScaling = atof(value);
    ComputeNonbondedUtil::select();
    return;
  }
  if ( ! strncasecmp(param,"commOnly",MAX_SCRIPT_PARAM_SIZE) ) {
    commOnly = atobool(value);
    ComputeNonbondedUtil::select();
    return;
  }

  char *error = new char[2 * MAX_SCRIPT_PARAM_SIZE + 100];
  sprintf(error,"Setting parameter %s from script failed!\n",param);
  NAMD_die(error);

}

/************************************************************************/
/*                                                                      */
/*      FUNCTION config_parser                                          */
/*                                                                      */
/************************************************************************/
                         
void SimParameters::config_parser(ParseOptions &opts) {

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

   config_parser_basic(opts);
   config_parser_fileio(opts);
   config_parser_fullelect(opts);
   config_parser_methods(opts);
   config_parser_constraints(opts);
   config_parser_movdrag(opts);
   config_parser_rotdrag(opts);
   config_parser_constorque(opts);
   config_parser_boundary(opts);
   config_parser_misc(opts);

}

void SimParameters::config_parser_basic(ParseOptions &opts) {
   
   //  So first we set up the ParseOptions objects so that it will check
   //  all of the logical rules that the configuration file must follow.

   ////// basic options
   opts.require("main", "timestep", "size of the timestep, in fs",
    &dt, 1.0);
   opts.range("timestep", NOT_NEGATIVE);
   opts.units("timestep", N_FSEC);

   opts.optional("main", "numsteps", "number of timesteps to perform",
    &N,0);
   opts.range("numsteps", NOT_NEGATIVE);

   opts.optional("main", "stepspercycle",
      "Number of steps between atom migrations", 
      &stepsPerCycle, 20);
   opts.range("stepspercycle", POSITIVE);

   opts.optional("main", "cutoff", "local electrostatic and Vdw distance", 
      &cutoff);
   opts.range("cutoff", POSITIVE);
   opts.units("cutoff", N_ANGSTROM);
   
   opts.optional("main", "nonbondedScaling", "nonbonded scaling factor",
     &nonbondedScaling, 1.0);
   opts.range("nonbondedScaling", NOT_NEGATIVE);

   opts.optional("main", "limitDist", "limit nonbonded below this distance",
     &limitDist, 0.0);
   opts.range("limitDist", NOT_NEGATIVE);

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
   opts.units("switchdist", N_ANGSTROM);

   opts.optional("switching", "elecswitchdist",
     "Distance for electrostatic switching function activation",
     &elecswitchDist);
   opts.range("elecswitchdist", POSITIVE);
   opts.units("elecswitchdist", N_ANGSTROM);

   opts.optional("switching", "eleccutoff",
     "Distance for electrostatic cutoff",
     &eleccutoff);
   opts.range("eleccutoff", POSITIVE);
   opts.units("eleccutoff", N_ANGSTROM);

   opts.optional("switching", "vdwswitchdist",
     "Distance for vdw switching function activation",
     &vdwswitchDist);
   opts.range("vdwswitchdist", POSITIVE);
   opts.units("vdwswitchdist", N_ANGSTROM);

   opts.optional("switching", "vdwcutoff",
     "Distance for vdw cutoff",
     &vdwcutoff);
   opts.range("vdwcutoff", POSITIVE);
   opts.units("vdwcutoff", N_ANGSTROM);

   opts.optional("main", "pairlistdist",  "Pairlist inclusion distance",
     &pairlistDist);
   opts.range("pairlistdist", POSITIVE);
   opts.units("pairlistdist", N_ANGSTROM);

   opts.optional("main", "pairlistMinProcs",  "Min procs for pairlists",
     &pairlistMinProcs,1);
   opts.range("pairlistMinProcs", POSITIVE);

   opts.optional("main", "pairlistsPerCycle",  "regenerate x times per cycle",
     &pairlistsPerCycle,2);
   opts.range("pairlistsPerCycle", POSITIVE);

   opts.optional("main", "outputPairlists", "how often to print warnings",
     &outputPairlists, 0);
   opts.range("outputPairlists", NOT_NEGATIVE);

   opts.optional("main", "pairlistShrink",  "tol *= (1 - x) on regeneration",
     &pairlistShrink,0.01);
   opts.range("pairlistShrink", NOT_NEGATIVE);

   opts.optional("main", "pairlistGrow",  "tol *= (1 + x) on trigger",
     &pairlistGrow, 0.01);
   opts.range("pairlistGrow", NOT_NEGATIVE);

   opts.optional("main", "pairlistTrigger",  "trigger is atom > (1 - x) * tol",
     &pairlistTrigger, 0.3);
   opts.range("pairlistTrigger", NOT_NEGATIVE);

   opts.optional("main", "temperature", "initial temperature",
     &initialTemp);
   opts.range("temperature", NOT_NEGATIVE);
   opts.units("temperature", N_KELVIN);

   opts.optionalB("main", "COMmotion", "allow initial center of mass movement",
      &comMove, FALSE);

   opts.optionalB("main", "zeroMomentum", "constrain center of mass",
      &zeroMomentum, FALSE);
   opts.optionalB("zeroMomentum", "zeroMomentumAlt", "constrain center of mass",
      &zeroMomentumAlt, FALSE);

   opts.optionalB("main", "wrapWater", "wrap waters around periodic boundaries on output",
      &wrapWater, FALSE);
   opts.optionalB("main", "wrapAll", "wrap all clusters around periodic boundaries on output",
      &wrapAll, FALSE);
   opts.optionalB("main", "wrapNearest", "wrap to nearest image to cell origin",
      &wrapNearest, FALSE);

   opts.optional("main", "dielectric", "dielectric constant",
     &dielectric, 1.0);
   opts.range("dielectric", POSITIVE); // Hmmm, dielectric < 1 ...

   opts.optional("main", "margin", "Patch width margin", &margin, XXXBIGREAL);
   opts.range("margin", NOT_NEGATIVE);
   opts.units("margin", N_ANGSTROM);

   opts.optional("main", "seed", "Initial random number seed", &randomSeed);
   opts.range("seed", POSITIVE);

   opts.optional("main", "outputEnergies", "How often to print energies in timesteps",
     &outputEnergies, 1);
   opts.range("outputEnergies", POSITIVE);
     
   opts.optional("main", "outputMomenta", "How often to print linear and angular momenta in timesteps",
     &outputMomenta, 0);
   opts.range("outputMomenta", NOT_NEGATIVE);
     
   opts.optional("main", "outputTiming", "How often to print timing data in timesteps",
     &outputTiming);
   opts.range("outputTiming", NOT_NEGATIVE);
     
   opts.optional("main", "outputPressure", "How often to print pressure data in timesteps",
     &outputPressure, 0);
   opts.range("outputPressure", NOT_NEGATIVE);
     
   opts.optionalB("main", "mergeCrossterms", "merge crossterm energy with dihedral when printing?",
      &mergeCrossterms, TRUE);

   opts.optional("main", "MTSAlgorithm", "Multiple timestep algorithm",
    PARSE_STRING);

   opts.optional("main", "longSplitting", "Long range force splitting option",
    PARSE_STRING);

   opts.optional("main", "splitPatch", "Atom into patch splitting option",
    PARSE_STRING);
   opts.optional("main", "hgroupCutoff", "Hydrogen margin", &hgroupCutoff, 2.5);

   opts.optional("main", "extendedSystem",
    "Initial configuration of extended system variables and periodic cell",
    PARSE_STRING);

   opts.optional("main", "cellBasisVector1", "Basis vector for periodic cell",
    &cellBasisVector1);
   opts.optional("main", "cellBasisVector2", "Basis vector for periodic cell",
    &cellBasisVector2);
   opts.optional("main", "cellBasisVector3", "Basis vector for periodic cell",
    &cellBasisVector3);
   opts.optional("main", "cellOrigin", "Fixed center of periodic cell",
    &cellOrigin);

   opts.optionalB("main", "molly", "Rigid bonds to hydrogen",&mollyOn,FALSE);
   opts.optional("main", "mollyTolerance", "Error tolerance for MOLLY",
                 &mollyTol, 0.00001);
   opts.optional("main", "mollyIterations", 
		 "Max number of iterations for MOLLY", &mollyIter, 100);

   opts.optional("main", "rigidBonds", "Rigid bonds to hydrogen",PARSE_STRING);
   opts.optional("main", "rigidTolerance", 
                 "Error tolerance for rigid bonds to hydrogen",
                 &rigidTol, 1.0e-8);
   opts.optional("main", "rigidIterations", 
		 "Max number of SHAKE iterations for rigid bonds to hydrogen",
		 &rigidIter, 100);
   opts.optionalB("main", "rigidDieOnError", 
		 "Die if rigidTolerance is not achieved after rigidIterations",
		 &rigidDie, TRUE);
   opts.optionalB("main", "useSettle",
                  "Use the SETTLE algorithm for rigid waters",
                 &useSettle, TRUE);

   opts.optional("main", "nonbondedFreq", "Nonbonded evaluation frequency",
    &nonbondedFrequency, 1);
   opts.range("nonbondedFreq", POSITIVE);

   opts.optionalB("main", "outputPatchDetails", "print number of atoms in each patch",
      &outputPatchDetails, FALSE);
}

void SimParameters::config_parser_fileio(ParseOptions &opts) {
   
   /////////////// file I/O

   opts.optional("main", "cwd", "current working directory", PARSE_STRING);

// In order to include AMBER options, "coordinates", "structure"
// and "parameters" are now optional, not required. The presence
// of them will be checked later in check_config()

//   opts.require("main", "coordinates", "initial PDB coordinate file",
//    PARSE_STRING);
   opts.optional("main", "coordinates", "initial PDB coordinate file",
    PARSE_STRING);

   opts.optional("main", "velocities",
     "initial velocities, given as a PDB file", PARSE_STRING);
   opts.optional("main", "binvelocities",
     "initial velocities, given as a binary restart", PARSE_STRING);
   opts.optional("main", "bincoordinates",
     "initial coordinates in a binary restart file", PARSE_STRING);

//   opts.require("main", "structure", "initial PSF structure file",
//    PARSE_STRING);
   opts.optional("main", "structure", "initial PSF structure file",
    PARSE_STRING);

//   opts.require("main", "parameters",
//"CHARMm 19 or CHARMm 22 compatable force field file (multiple "
//"inputs allowed)", PARSE_MULTIPLES);
   opts.optional("main", "parameters",
"CHARMm 19 or CHARMm 22 compatable force field file (multiple "
"inputs allowed)", PARSE_MULTIPLES);


   //****** BEGIN CHARMM/XPLOR type changes
   //// enable XPLOR as well as CHARMM input files for parameters
   opts.optionalB("parameters", "paraTypeXplor", "Parameter file in Xplor format?", &paraTypeXplorOn, FALSE);
   opts.optionalB("parameters", "paraTypeCharmm", "Parameter file in Charmm format?", &paraTypeCharmmOn, FALSE); 
   //****** END CHARMM/XPLOR type changes
   
   opts.require("main", "outputname",
    "prefix for the final PDB position and velocity filenames", 
    outputFilename);

   opts.optional("main", "auxFile", "Filename for data stream output",
     auxFilename);

   opts.optional("main", "DCDfreq", "Frequency of DCD trajectory output, in "
    "timesteps", &dcdFrequency, 0);
   opts.range("DCDfreq", NOT_NEGATIVE);
   opts.optional("DCDfreq", "DCDfile", "DCD trajectory output file name",
     dcdFilename);
   opts.optionalB("DCDfreq", "DCDunitcell", "Store unit cell in dcd timesteps?",
       &dcdUnitCell);

   opts.optional("main", "velDCDfreq", "Frequency of velocity "
    "DCD output, in timesteps", &velDcdFrequency, 0);
   opts.range("velDCDfreq", NOT_NEGATIVE);
   opts.optional("velDCDfreq", "velDCDfile", "velocity DCD output file name",
     velDcdFilename);
   
   opts.optional("main", "XSTfreq", "Frequency of XST trajectory output, in "
    "timesteps", &xstFrequency, 0);
   opts.range("XSTfreq", NOT_NEGATIVE);
   opts.optional("XSTfreq", "XSTfile", "Extended sytem trajectory output "
    "file name", xstFilename);

   opts.optional("main", "restartfreq", "Frequency of restart file "
    "generation", &restartFrequency, 0);
   opts.range("restartfreq", NOT_NEGATIVE);
   opts.optional("restartfreq", "restartname", "Prefix for the position and "
     "velocity PDB files used for restarting", restartFilename);
   opts.optionalB("restartfreq", "restartsave", "Save restart files with "
     "unique filenames rather than overwriting", &restartSave, FALSE);

   opts.optionalB("restartfreq", "binaryrestart", "Specify use of binary restart files ", 
       &binaryRestart, TRUE);

   opts.optionalB("outputname", "binaryoutput", "Specify use of binary output files ", 
       &binaryOutput, TRUE);

   opts.optionalB("main", "amber", "Is it AMBER force field?",
       &amberOn, FALSE);
   opts.optionalB("amber", "readexclusions", "Read exclusions from parm file?",
       &readExclusions, TRUE);
   opts.require("amber", "scnb", "1-4 VDW interactions are divided by scnb",
       &vdwscale14, 2.0);
   opts.require("amber", "parmfile", "AMBER parm file", PARSE_STRING);
   opts.optional("amber", "ambercoor", "AMBER coordinate file", PARSE_STRING);

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
// begin fep output options    
   opts.optional("fep", "fepoutfreq", "Frequency of FEP energy output in "
     "timesteps", &fepOutFreq, 5);
   opts.range("fepoutfreq", NOT_NEGATIVE);
   opts.optional("fepoutfreq", "fepoutfile", "FEP energy output filename",
     fepOutFile);
// end fep output options
//fepe

   /* GROMACS options */
   opts.optionalB("main", "gromacs", "Use GROMACS-like force field?",
       &gromacsOn, FALSE);
   opts.require("gromacs", "grotopfile", "GROMACS topology file",
		PARSE_STRING);
   opts.optional("gromacs", "grocoorfile","GROMACS coordinate file",
		 PARSE_STRING);

  // OPLS options
   opts.optionalB("main", "vdwGeometricSigma",
       "Use geometric mean to combine L-J sigmas, as for OPLS",
       &vdwGeometricSigma, FALSE);

}


void SimParameters::config_parser_fullelect(ParseOptions &opts) {
   
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

   opts.optional("main", "fullElectFrequency",
      "Number of steps between full electrostatic executions", 
      &fullElectFrequency);
   opts.range("fullElectFrequency", POSITIVE);

   //  USE OF THIS PARAMETER DISCOURAGED
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


   ///////////  Particle Mesh Ewald

   opts.optionalB("main", "PME", "Use particle mesh Ewald for electrostatics?",
	&PMEOn, FALSE);
   opts.optional("PME", "PMETolerance", "PME direct space tolerance",
	&PMETolerance, 1.e-6);
   opts.optional("PME", "PMEInterpOrder", "PME interpolation order",
	&PMEInterpOrder, 4);  // cubic interpolation is default
   opts.require("PME", "PMEGridSizeX", "PME grid in x dimension",
	&PMEGridSizeX);
   opts.require("PME", "PMEGridSizeY", "PME grid in y dimension",
	&PMEGridSizeY);
   opts.require("PME", "PMEGridSizeZ", "PME grid in z dimension",
	&PMEGridSizeZ);
   opts.optional("PME", "PMEProcessors",
	"PME FFT and reciprocal sum processor count", &PMEProcessors, 0);
   opts.optionalB("main", "PMEBarrier", "Use barrier in PME?",
	&PMEBarrier, FALSE);

#ifdef DPME
   opts.optionalB("PME", "useDPME", "Use old DPME code?", &useDPME, FALSE);
#else
   useDPME = 0;
#endif

   opts.optionalB("main", "FFTWEstimate", "Use estimates to optimize FFTW?",
	&FFTWEstimate, FALSE);
   opts.optionalB("main", "FFTWUseWisdom", "Read/save wisdom file for FFTW?",
	&FFTWUseWisdom, TRUE);
   opts.optional("FFTWUseWisdom", "FFTWWisdomFile", "File for FFTW wisdom",
	FFTWWisdomFile);

}

void SimParameters::config_parser_methods(ParseOptions &opts) {
   
   /////////// Special Dynamics Methods
   opts.optionalB("main", "minimization", "Should minimization be performed?",
      &minimizeCGOn, FALSE);
   opts.optional("main", "minTinyStep", "very first minimization steps",
      &minTinyStep, 1.0e-6);
   opts.range("minTinyStep", POSITIVE);
   opts.optional("main", "minBabyStep", "initial minimization steps",
      &minBabyStep, 1.0e-2);
   opts.range("minBabyStep", POSITIVE);
   opts.optional("main", "minLineGoal", "line minimization gradient reduction",
      &minLineGoal, 1.0e-4);
   opts.range("minLineGoal", POSITIVE);

   opts.optionalB("main", "velocityQuenching",
      "Should old-style minimization be performed?", &minimizeOn, FALSE);

   opts.optional("main", "maximumMove", "Maximum atom movement per step", &maximumMove, 0.0);
   opts.range("maximumMove", NOT_NEGATIVE);
   opts.units("maximumMove", N_ANGSTROM);

   opts.optionalB("main", "Langevin", "Should Langevin dynamics be performed?",
      &langevinOn, FALSE);
   opts.require("Langevin", "langevinTemp", "Temperature for heat bath in Langevin "
     "dynamics", &langevinTemp);
   opts.range("langevinTemp", NOT_NEGATIVE);
   opts.units("langevinTemp", N_KELVIN);
   opts.optional("Langevin", "langevinDamping", "Damping coefficient (1/ps)",
      &langevinDamping);
   opts.range("langevinDamping", POSITIVE);
   opts.optionalB("Langevin", "langevinHydrogen", "Should Langevin dynamics be applied to hydrogen atoms?",
      &langevinHydrogen);
   opts.optional("Langevin", "langevinFile", "PDB file with temperature "
     "coupling terms (B(i)) (default is the PDB input file)",
     PARSE_STRING);
   opts.optional("Langevin", "langevinCol", "Column in the langevinFile "
     "containing the temperature coupling term B(i);\n"
     "default is 'O'", PARSE_STRING);

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
//  alchemical fep options
   opts.optionalB("main", "fep", "Is chemical fep being performed?",
     &fepOn, FALSE);
   opts.require("fep", "lambda", "Coupling parameter value", &lambda);
   opts.require("fep", "lambda2", "Coupling comparison value", &lambda2);
   opts.optional("fep", "fepFile", "PDB file with perturbation flags "
     "default is the input PDB file", PARSE_STRING); 
   opts.optional("fep", "fepCol", "Column in the fepFile with the "
     "perturbation flag", PARSE_STRING);
   opts.optional("fep", "fepEquilSteps", "Equilibration steps, before "
     "data collection in the fep window", &fepEquilSteps, 0);
   opts.range("fepEquilSteps", NOT_NEGATIVE);
// end FEP options
//fepe

   opts.optionalB("main", "les", "Is locally enhanced sampling enabled?",
     &lesOn, FALSE);
   opts.require("les", "lesFactor", "Local enhancement factor", &lesFactor);
   opts.optional("les", "lesFile", "PDB file with enhancement flags "
     "default is the input PDB file", PARSE_STRING); 
   opts.optional("les", "lesCol", "Column in the lesFile with the "
     "enhancement flag", PARSE_STRING);
   opts.optionalB("les", "lesReduceTemp", "Reduce enhanced atom temperature?",
     &lesReduceTemp, FALSE);
   opts.optionalB("les", "lesReduceMass", "Reduce enhanced atom mass?",
     &lesReduceMass, FALSE);

   // Pair interaction calculations
    opts.optionalB("main", "pairInteraction", 
	"Are pair interactions calculated?", &pairInteractionOn, FALSE);
    opts.optional("pairInteraction", "pairInteractionFile", 
	"PDB files with interaction flags " "default is the input PDB file", 
	PARSE_STRING);
    opts.optional("pairInteraction", "pairInteractionCol", 
	"Column in the pairInteractionFile with the interaction flags",
	PARSE_STRING);
    opts.require("pairInteraction", "pairInteractionGroup1",
        "Flag for interaction group 1", &pairInteractionGroup1);
    opts.optional("pairInteraction", "pairInteractionGroup2",
        "Flag for interaction group 2", &pairInteractionGroup2, -1);
    opts.optionalB("pairInteraction", "pairInteractionSelf",
        "Compute only within-group interactions?", &pairInteractionSelf, 
        FALSE);

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
   opts.units("COLDTemp", N_KELVIN);
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
   opts.units("tCoupleTemp", N_KELVIN);
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
   opts.units("rescaleTemp", N_KELVIN);

   opts.optional("main", "reassignFreq", "Number of steps between "
    "velocity reassignment", &reassignFreq);
   opts.range("reassignFreq", POSITIVE);
   opts.optional("main", "reassignTemp", "Target temperature for velocity reassignment",
    &reassignTemp);
   opts.range("reassignTemp", NOT_NEGATIVE);
   opts.units("reassignTemp", N_KELVIN);
   opts.optional("main", "reassignIncr", "Temperature increment for velocity reassignment",
    &reassignIncr);
   opts.units("reassignIncr", N_KELVIN);
   opts.optional("main", "reassignHold", "Final holding temperature for velocity reassignment",
    &reassignHold);
   opts.range("reassignHold", NOT_NEGATIVE);
   opts.units("reassignHold", N_KELVIN);

   ////  Group rather than atomic pressure
   opts.optionalB("main", "useGroupPressure", 
      "Use group rather than atomic quantities for pressure control?",
      &useGroupPressure, FALSE);

   ////  Anisotropic cell fluctuations
   opts.optionalB("main", "useFlexibleCell",
      "Use anisotropic cell fluctuation for pressure control?",
      &useFlexibleCell, FALSE);

   ////  Constant dimension ratio in X-Y plane
   opts.optionalB("main", "useConstantRatio",
      "Use constant X-Y ratio for pressure control?",
      &useConstantRatio, FALSE);

   ////  Constant area and normal pressure conditions
   opts.optionalB("main", "useConstantArea",
      "Use constant area for pressure control?",
      &useConstantArea, FALSE);

   //// Exclude atoms from pressure
   opts.optionalB("main", "excludeFromPressure",
	"Should some atoms be excluded from pressure rescaling?",
	&excludeFromPressure, FALSE);
   opts.optional("excludeFromPressure", "excludeFromPressureFile",
	"PDB file for atoms to be excluded from pressure",
        PARSE_STRING);
   opts.optional("excludeFromPressure", "excludeFromPressureCol", 
        "Column in the excludeFromPressureFile"
        "containing the flags (nonzero means excluded);\n"
        "default is 'O'", PARSE_STRING);

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
   opts.units("BerendsenPressureRelaxationTime", N_FSEC);
   opts.optional("BerendsenPressure", "BerendsenPressureFreq",
    "Number of steps between volume rescaling",
    &berendsenPressureFreq, 1);
   opts.range("BerendsenPressureFreq", POSITIVE);

   ////  Langevin Piston pressure control
   opts.optionalB("main", "LangevinPiston",
      "Should Langevin piston pressure control be used?",
      &langevinPistonOn, FALSE);
   opts.require("LangevinPiston", "LangevinPistonTarget",
      "Target pressure for pressure control",
      &langevinPistonTarget);
   opts.require("LangevinPiston", "LangevinPistonPeriod",
      "Oscillation period for pressure control",
      &langevinPistonPeriod);
   opts.range("LangevinPistonPeriod", POSITIVE);
   opts.units("LangevinPistonPeriod", N_FSEC);
   opts.require("LangevinPiston", "LangevinPistonDecay",
      "Decay time for pressure control",
      &langevinPistonDecay);
   opts.range("LangevinPistonDecay", POSITIVE);
   opts.units("LangevinPistonDecay", N_FSEC);
   opts.require("LangevinPiston", "LangevinPistonTemp",
      "Temperature for pressure control piston",
      &langevinPistonTemp);
   opts.range("LangevinPistonTemp", POSITIVE);
   opts.units("LangevinPistonTemp", N_KELVIN);
   opts.optional("LangevinPiston", "StrainRate",
      "Initial strain rate for pressure control (x y z)",
      &strainRate);

   //// Surface tension 
   opts.optional("main", "SurfaceTensionTarget",
      "Surface tension in the x-y plane",
      &surfaceTensionTarget, 0);

   //// Pressure Profile calculations
   opts.optionalB("main", "pressureprofile", "Compute pressure profile?",
     &pressureProfileOn, FALSE);
   opts.require("pressureprofile", "pressureprofileslabs", 
     "Number of pressure profile slabs", &pressureProfileSlabs, 10);
   opts.optional("pressureprofile", "pressureprofilefreq",
     "How often to store profile data", &pressureProfileFreq, 1);
   opts.optional("pressureprofile", "pressureProfileAtomTypes", 
     "Number of pressure profile atom types", &pressureProfileAtomTypes, 1);
   opts.range("pressureProfileAtomTypes", POSITIVE);
   opts.optional("pressureProfile", "pressureProfileAtomTypesFile", 
	"PDB files with pressure profile atom types" "default is the input PDB file", 
	PARSE_STRING);
   opts.optional("pressureProfile", "pressureProfileAtomTypesCol", 
	"Column in the pressureProfileAtomTypesFile with the atom types ",
	PARSE_STRING);
   opts.optionalB("pressureProfile", "pressureProfileEwald", 
       "Compute Ewald contribution to pressure profile",
       &pressureProfileEwaldOn, FALSE);
   opts.optional("pressureProfile", "pressureProfileEwaldX",
       "Ewald grid size X", &pressureProfileEwaldX, 10);
   opts.range("pressureProfileEwaldX", POSITIVE);
   opts.optional("pressureProfile", "pressureProfileEwaldY",
       "Ewald grid size Y", &pressureProfileEwaldY, 10);
   opts.range("pressureProfileEwaldY", POSITIVE);
   opts.optional("pressureProfile", "pressureProfileEwaldZ",
       "Ewald grid size Z", &pressureProfileEwaldZ, 10);
   opts.range("pressureProfileEwaldZ", POSITIVE);
}

void SimParameters::config_parser_constraints(ParseOptions &opts) {
   
   ////  Fixed Atoms
   opts.optionalB("main", "fixedatoms", "Are there fixed atoms?",
    &fixedAtomsOn, FALSE);
   opts.optionalB("fixedatoms", "fixedAtomsForces",
     "Calculate forces between fixed atoms?  (Required to unfix during run.)",
     &fixedAtomsForces, FALSE);
   opts.optional("fixedatoms", "fixedAtomsFile", "PDB file with flags for "
     "fixed atoms (default is the PDB input file)",
     PARSE_STRING);
   opts.optional("fixedatoms", "fixedAtomsCol", "Column in the fixedAtomsFile "
     "containing the flags (nonzero means fixed);\n"
     "default is 'O'", PARSE_STRING);

   ////  Harmonic Constraints
   opts.optionalB("main", "constraints", "Are harmonic constraints active?",
     &constraintsOn, FALSE);
   opts.require("constraints", "consexp", "Exponent for harmonic potential",
    &constraintExp, 2);
   opts.range("consexp", POSITIVE);
   opts.require("constraints", "consref", "PDB file containing reference "
    "positions",
    PARSE_STRING);
   opts.require("constraints", "conskfile", "PDB file containing force "
    "constaints in one of the columns", PARSE_STRING);
   opts.require("constraints", "conskcol", "Column of conskfile to use "
    "for the force constants", PARSE_STRING);
   opts.require("constraints", "constraintScaling", "constraint scaling factor",
     &constraintScaling, 1.0);
   opts.range("constraintScaling", NOT_NEGATIVE);



   //****** BEGIN selective restraints (X,Y,Z) changes

   //// selective restraints (X,Y,Z) 
   opts.optionalB("constraints", "selectConstraints", 
   "Restrain only selected Cartesian components of the coordinates?",
     &selectConstraintsOn, FALSE);
   opts.optionalB("selectConstraints", "selectConstrX",  
   "Restrain X components of coordinates ", &constrXOn, FALSE);
   opts.optionalB("selectConstraints", "selectConstrY",  
   "Restrain Y components of coordinates ", &constrYOn, FALSE);
   opts.optionalB("selectConstraints", "selectConstrZ",  
   "Restrain Z components of coordinates ", &constrZOn, FALSE);
   //****** END selective restraints (X,Y,Z) changes
 

   //****** BEGIN moving constraints changes 

   //// Moving Harmonic Constraints
   opts.optionalB("constraints", "movingConstraints",
      "Are some of the constraints moving?", 
      &movingConstraintsOn, FALSE);
   opts.require("movingConstraints", "movingConsVel",
    "Velocity of the movement, A/timestep", &movingConsVel);
   //****** END moving constraints changes 

   // BEGIN rotating constraints changes
   opts.optionalB("constraints", "rotConstraints",
      "Are the constraints rotating?", 
      &rotConstraintsOn, FALSE);
   opts.require("rotConstraints", "rotConsAxis",
    "Axis of rotation", &rotConsAxis);
   opts.require("rotConstraints", "rotConsPivot",
    "Pivot point of rotation", 
    &rotConsPivot);
   opts.require("rotConstraints", "rotConsVel",
    "Velocity of rotation, deg/timestep", &rotConsVel);

   // END rotating constraints changes

   // external command forces
   opts.optionalB("main", "extForces", "External command forces?",
      &extForcesOn, FALSE);
   opts.require("extForces", "extForcesCommand",
      "External forces command", extForcesCommand);
   opts.require("extForces", "extCoordFilename",
      "External forces coordinate filename", extCoordFilename);
   opts.require("extForces", "extForceFilename",
      "External forces force filename", extForceFilename);

   //****** BEGIN SMD constraints changes 

   // SMD constraints
   opts.optionalB("main", "SMD",
      "Do we use SMD option?", 
      &SMDOn, FALSE);
   opts.require("SMD", "SMDVel",
		"Velocity of the movement, A/timestep", &SMDVel);
   opts.range("SMDVel", NOT_NEGATIVE);
   opts.require("SMD", "SMDDir",
		"Direction of movement", &SMDDir);
   opts.require("SMD", "SMDk",
                "Elastic constant for SMD", &SMDk);
   opts.range("SMDk", NOT_NEGATIVE);
   opts.require("SMD", "SMDFile",
		"File for SMD information",
                 SMDFile);
   opts.optional("SMD", "SMDOutputFreq",
		 "Frequency of output",
		 &SMDOutputFreq, 1);
   opts.range("SMDOutputFreq", POSITIVE);
   
   //****** END SMD constraints changes 

   // TMD parameters
   opts.optionalB("main", "TMD", "Perform Targeted MD?", &TMDOn, FALSE);
   opts.require("TMD", "TMDk", "Elastic constant for TMD", &TMDk);
   opts.range("TMDk", NOT_NEGATIVE);
   opts.require("TMD", "TMDFile", "File for TMD information", TMDFile);
   opts.optional("TMD", "TMDOutputFreq", "Frequency of TMD output", 
       &TMDOutputFreq, 1);
   opts.range("TMDOutputFreq", POSITIVE);
   opts.require("TMD", "TMDLastStep", "Last TMD timestep", &TMDLastStep);
   opts.range("TMDLastStep", POSITIVE);
   opts.optional("TMD", "TMDFirstStep", "First TMD step (default 0)", &TMDFirstStep, 0);
   opts.optional("TMD", "TMDInitialRMSD", "Target RMSD at first TMD step (default 0 to use initial coordinates)", &TMDInitialRMSD);
   opts.optional("TMD", "TMDFinalRMSD", "Target RMSD at last TMD step (default 0 )", &TMDFinalRMSD, 0);
   opts.range("TMDInitialRMSD", NOT_NEGATIVE);

   // End of TMD parameters

   ////  Global Forces / Tcl
   opts.optionalB("main", "tclForces", "Are Tcl global forces active?",
     &tclForcesOn, FALSE);
   opts.require("tclForces", "tclForcesScript",
     "Tcl script for global forces", PARSE_MULTIPLES);

   ////  Boundary Forces / Tcl
   opts.optionalB("main", "tclBC", "Are Tcl boundary forces active?",
     &tclBCOn, FALSE);
   opts.require("tclBC", "tclBCScript",
     "Tcl script defining calcforces for boundary forces", PARSE_STRING);
   tclBCScript = 0;
   opts.optional("tclBC", "tclBCArgs", "Extra args for calcforces command",
     tclBCArgs);
   tclBCArgs[0] = 0;

   ////  Global Forces / Misc
   opts.optionalB("main", "miscForces", "Are misc global forces active?",
     &miscForcesOn, FALSE);
   opts.optional("miscForces", "miscForcesScript",
     "script for misc forces", PARSE_MULTIPLES);

   ////  Free Energy Perturbation
   opts.optionalB("main", "freeEnergy", "Perform free energy perturbation?",
     &freeEnergyOn, FALSE);
   opts.require("freeEnergy", "freeEnergyConfig",
     "Configuration file for free energy perturbation", PARSE_MULTIPLES);

   ////  Constant Force
   opts.optionalB("main", "constantforce", "Apply constant force?",
     &consForceOn, FALSE);
   opts.optional("constantforce", "consForceFile",
       "Configuration file for constant forces", PARSE_STRING);
}

void SimParameters::config_parser_movdrag(ParseOptions &opts) {
   //// moving drag
   opts.optionalB("main", "movDragOn", "Do we apply moving drag?",
      &movDragOn, FALSE);
   opts.require("movDragOn", "movDragFile",
      "Main moving drag PDB file", movDragFile);
   opts.require("movDragOn", "movDragCol",
      "Main moving drag PDB column", PARSE_STRING);
   opts.require("movDragOn", "movDragGlobVel",
      "Global moving drag velocity (A/step)", &movDragGlobVel);
   opts.require("movDragOn", "movDragVelFile",
      "Moving drag linear velocity file", movDragVelFile);
}

void SimParameters::config_parser_rotdrag(ParseOptions &opts) {
   //// rotating drag
   opts.optionalB("main", "rotDragOn", "Do we apply rotating drag?",
      &rotDragOn, FALSE);
   opts.require("rotDragOn", "rotDragFile",
      "Main rotating drag PDB file", rotDragFile);
   opts.require("rotDragOn", "rotDragCol",
      "Main rotating drag PDB column", PARSE_STRING);
   opts.require("rotDragOn", "rotDragAxisFile",
      "Rotating drag axis file", rotDragAxisFile);
   opts.require("rotDragOn", "rotDragPivotFile",
      "Rotating drag pivot point file", rotDragPivotFile);
   opts.require("rotDragOn", "rotDragGlobVel",
      "Global rotating drag angular velocity (deg/step)", &rotDragGlobVel);
   opts.require("rotDragOn", "rotDragVelFile",
      "Rotating drag angular velocity file", rotDragVelFile);
   opts.require("rotDragOn", "rotDragVelCol",
      "Rotating drag angular velocity column", PARSE_STRING);
}

void SimParameters::config_parser_constorque(ParseOptions &opts) {
   //// "constant" torque
   opts.optionalB("main", "consTorqueOn", "Do we apply \"constant\" torque?",
      &consTorqueOn, FALSE);
   opts.require("consTorqueOn", "consTorqueFile",
      "Main \"constant\" torque PDB file", consTorqueFile);
   opts.require("consTorqueOn", "consTorqueCol",
      "Main \"constant\" torque PDB column", PARSE_STRING);
   opts.require("consTorqueOn", "consTorqueAxisFile",
      "\"Constant\" torque axis file", consTorqueAxisFile);
   opts.require("consTorqueOn", "consTorquePivotFile",
      "\"Constant\" torque pivot point file", consTorquePivotFile);
   opts.require("consTorqueOn", "consTorqueGlobVal",
      "Global \"constant\" torque value (Kcal/(mol*A^2))", &consTorqueGlobVal);
   opts.require("consTorqueOn", "consTorqueValFile",
      "\"constant\" torque factors file", consTorqueValFile);
   opts.require("consTorqueOn", "consTorqueValCol",
      "\"constant\" torque factors column", PARSE_STRING);
}

void SimParameters::config_parser_boundary(ParseOptions &opts) {
    
   //// Spherical Boundary Conditions
   opts.optionalB("main", "sphericalBC", "Are spherical boundary counditions "
      "active?", &sphericalBCOn, FALSE);
   opts.require("sphericalBC", "sphericalBCCenter",
     "Center of spherical boundaries", &sphericalCenter);
   opts.require("sphericalBC", "sphericalBCr1", "Radius for first sphere "
     "potential", &sphericalBCr1);
   opts.range("sphericalBCr1", POSITIVE);
   opts.units("sphericalBCr1", N_ANGSTROM);
   opts.require("sphericalBC", "sphericalBCk1", "Force constant for first "
    "sphere potential (+ is an inward force, - outward)",
    &sphericalBCk1);
   opts.units("sphericalBCk1", N_KCAL);
   opts.optional("sphericalBC", "sphericalBCexp1", "Exponent for first "
    "sphere potential", &sphericalBCexp1, 2);
   opts.range("sphericalBCexp1", POSITIVE);
   
   opts.optional("sphericalBCr1", "sphericalBCr2", "Radius for second sphere "
     "potential", &sphericalBCr2);
   opts.range("sphericalBCr2", POSITIVE);
   opts.units("sphericalBCr2", N_ANGSTROM);
   opts.require("sphericalBCr2", "sphericalBCk2", "Force constant for second "
    "sphere potential (+ is an inward force, - outward)",
    &sphericalBCk2);
   opts.units("sphericalBCk2", N_KCAL);
   opts.optional("sphericalBCr2", "sphericalBCexp2", "Exponent for second "
    "sphere potential", &sphericalBCexp2, 2);
   opts.range("sphericalBCexp2", POSITIVE);

   /////////////// Cylindrical Boundary Conditions
   opts.optionalB("main", "cylindricalBC", "Are cylindrical boundary counditions "
                  "active?", &cylindricalBCOn, FALSE);
   opts.require("cylindricalBC", "cylindricalBCr1", "Radius for first cylinder "
                 "potential", &cylindricalBCr1);
   opts.range("cylindricalBCr1", POSITIVE);
   opts.units("cylindricalBCr1", N_ANGSTROM);
   opts.require("cylindricalBC", "cylindricalBCk1", "Force constant for first "
                "cylinder potential (+ is an inward force, - outward)",
                &cylindricalBCk1);
   opts.units("cylindricalBCk1", N_KCAL);
   opts.optional("cylindricalBC", "cylindricalBCexp1", "Exponent for first "
                "cylinder potential", &cylindricalBCexp1, 2);
   opts.range("cylindricalBCexp1", POSITIVE);


// additions beyond those already found in spherical parameters    JJU
   opts.optional("cylindricalBC", "cylindricalBCAxis", "Cylinder axis (defaults to x)",
    PARSE_STRING);
   opts.require("cylindricalBC", "cylindricalBCCenter",
     "Center of cylindrical boundaries", &cylindricalCenter);
   opts.require ("cylindricalBC", "cylindricalBCl1", "Length of first cylinder",
                 &cylindricalBCl1);
   opts.range("cylindricalBCl1", POSITIVE);
   opts.units("cylindricalBCl1", N_ANGSTROM);
   opts.optional ("cylindricalBCl1", "cylindricalBCl2", "Length of second cylinder",
                  &cylindricalBCl2);
   opts.range ("cylindricalBCl2", POSITIVE);
   opts.units ("cylindricalBCl2", N_ANGSTROM);
// end  additions

   opts.optional("cylindricalBCr1", "cylindricalBCr2", "Radius for second cylinder "
                 "potential", &cylindricalBCr2);
   opts.range("cylindricalBCr2", POSITIVE);
   opts.units("cylindricalBCr2", N_ANGSTROM);
   opts.require("cylindricalBCr2", "cylindricalBCk2", "Force constant for second "
                "cylinder potential (+ is an inward force, - outward)",
                &cylindricalBCk2);
   opts.units("cylindricalBCk2", N_KCAL);
   opts.optional("cylindricalBCr2", "cylindricalBCexp2", "Exponent for second "
                "cylinder potential", &cylindricalBCexp2, 2);
   opts.range("cylindricalBCexp2", POSITIVE);

   ///////////////  Electric field options
   opts.optionalB("main", "eFieldOn", "Should and electric field be applied",
                 &eFieldOn, FALSE);
   opts.require("eFieldOn", "eField", "Electric field vector", &eField);

      ///////////////  Stir options
   opts.optionalB("main", "stirOn", "Should stirring torque be applied",
                 &stirOn, FALSE);
   opts.optional("stirOn", "stirFilename", "PDB file with flags for "
     "stirred atoms (default is the PDB input file)",
		 PARSE_STRING);
   opts.optional("stirOn", "stirredAtomsCol", "Column in the stirredAtomsFile "
		 "containing the flags (nonzero means fixed);\n"
		 "default is 'O'", PARSE_STRING);
   opts.require("stirOn", "stirStartingTheta", "Stir starting theta offset", &stirStartingTheta);
   opts.require("stirOn", "stirK", "Stir force harmonic spring constant", &stirK);
   //should make this optional, compute from firsttimestep * stirVel
   opts.require("stirOn", "stirVel", "Stir angular velocity (deg/timestep)", &stirVel);
   opts.require("stirOn", "stirAxis", "Stir axis (direction vector)", &stirAxis);
   opts.require("stirOn", "stirPivot", "Stir pivot point (coordinate)", &stirPivot);
}

void SimParameters::config_parser_misc(ParseOptions &opts) {
   
   ///////////////  Load balance options
   opts.optional("main", "ldbStrategy", "Load balancing strategy",
     loadStrategy);
   opts.optional("main", "ldbPeriod",
     "steps between load balancing", &ldbPeriod);
   opts.range("ldbPeriod", POSITIVE);
   opts.optional("main", "firstLdbStep", 
     "when to start load balancing",
     &firstLdbStep);
   opts.range("firstLdbStep", POSITIVE);
   opts.optional("main", "ldbBackgroundScaling",
     "background load scaling", &ldbBackgroundScaling);
   opts.range("ldbBackgroundScaling", NOT_NEGATIVE);
   opts.optional("main", "ldbPMEBackgroundScaling",
     "PME node background load scaling", &ldbPMEBackgroundScaling);
   opts.range("ldbPMEBackgroundScaling", NOT_NEGATIVE);
   opts.optional("main", "ldbHomeBackgroundScaling",
     "home node background load scaling", &ldbHomeBackgroundScaling);
   opts.range("ldbHomeBackgroundScaling", NOT_NEGATIVE);
   opts.optionalB("main", "ldbUnloadPME", "no load on PME nodes",
     &ldbUnloadPME, FALSE);
   opts.optionalB("main", "ldbUnloadSMP", "no load on one pe of SMP node",
     &ldbUnloadSMP, FALSE);
   opts.optionalB("main", "ldbUnloadZero", "no load on pe zero",
     &ldbUnloadZero, FALSE);
   opts.optionalB("main", "ldbUnloadRankZero", "no load on rank zero",
     &ldbUnloadRankZero, FALSE);
   opts.optionalB("main", "noPatchesOnZero", "no patches on pe zero",
     &noPatchesOnZero, FALSE);
   opts.optional("main", "procsPerNode", "Number of Processor per node",
     &procsPerNode);
   opts.range("procsPerNode", POSITIVE);
   opts.optional("main", "ldbUnloadRank", "no load on rank pe of a node",
     &ldbUnloadRank);
   opts.range("ldbUnloadRank", POSITIVE);
   opts.optionalB("main", "twoAwayX", "half-size patches in 1st dimension",
     &twoAwayX, FALSE);
   opts.optionalB("main", "twoAwayY", "half-size patches in 1st dimension",
     &twoAwayY, FALSE);
   opts.optionalB("main", "twoAwayZ", "half-size patches in 1st dimension",
     &twoAwayZ, FALSE);

   /////  Restart timestep option
   opts.optional("main", "firsttimestep", "Timestep to start simulation at",
     &firstTimestep, 0);
   opts.range("firsttimestep", NOT_NEGATIVE);
 
   /////  Test mode options
   opts.optionalB("main", "test", "Perform self-tests rather than simulation",
		&testOn, FALSE);
   opts.optionalB("main", "commOnly", "Do not evaluate forces or integrate",
		&commOnly, FALSE);

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
   opts.units("hbCutoffDist", N_ANGSTROM);
   opts.optional("hbonds","hbOnDist","Hbond A-D switch function on distance",
                 &daOnDist, 5.5);
   opts.range("hbOnDist", POSITIVE);
   opts.units("hbOnDist", N_ANGSTROM);
   opts.optional("hbonds","hbOffDist","Hbond A-D switch function off distance",
                 &daOffDist, 6.5);
   opts.range("hbOffDist", POSITIVE);
   opts.units("hbOffDist", N_ANGSTROM);

   // IMD options
   opts.optionalB("main","IMDon","Connect using IMD?",&IMDon, FALSE);
   opts.require("IMDon","IMDport", "Port to which to bind", &IMDport);
   opts.range("IMDport",POSITIVE);
   opts.require("IMDon","IMDfreq", "Frequency at which to report", &IMDfreq);
   opts.range("IMDfreq",POSITIVE);
   opts.optionalB("IMDon","IMDwait","Pause until IMD connection?",&IMDwait,
     FALSE);
   opts.optionalB("IMDon","IMDignore","Ignore forces, etc.?",&IMDignore,
     FALSE);

   // Maximum Partition options
   opts.optional("main", "maxSelfPart", 
     "maximum number of self partitions in one patch", &maxSelfPart, 50);
   opts.range("maxSelfPart",POSITIVE);
   opts.optional("main", "maxPairPart", 
     "maximum number of pair partitions in one patch", &maxPairPart, 20);
   opts.range("maxPairPart",POSITIVE);
   opts.optional("main", "numAtomsSelf", 
               "maximum number of atoms in one self compute distribution", 
               &numAtomsSelf, 125);
   opts.range("numAtomsSelf",NOT_NEGATIVE);
   opts.optional("main", "numAtomsPair", 
               "maximum number of atoms in one pair compute distribution", 
               &numAtomsPair, 200);
   opts.range("numAtomsPair",NOT_NEGATIVE);
   opts.optional("main", "numAtomsPair2", 
               "maximum number of atoms in one pair compute distribution", 
               &numAtomsPair2, 400);
   opts.range("numAtomsPair2",NOT_NEGATIVE);
   opts.optional("main", "minAtomsPerPatch", 
               "minimum average atoms per patch", 
               &minAtomsPerPatch, 100);
   opts.range("minAtomsPerPatch",NOT_NEGATIVE);

   // Maximum exclusion flags per atom
   opts.optional("main", "maxExclusionFlags", 
     "maximum number of exclusion flags per atom", &maxExclusionFlags, 100);
   opts.range("maxExclusionFlags",POSITIVE);
}

void SimParameters::check_config(ParseOptions &opts, ConfigList *config, char *&cwd) {
   
   int len;    //  String length
   StringList *current; //  Pointer to config option list

   //  Take care of cwd processing
   if (opts.defined("cwd"))
   {
    //  First allocate and get the cwd value
    current = config->find("cwd");

    len = strlen(current->data);

    if ( CHDIR(current->data) )
    {
      NAMD_die("chdir() to given cwd failed!");
    } else {
      iout << iINFO << "Changed directory to " << current->data << "\n" << endi;
    }

    if (current->data[len-1] != PATHSEP)
      len++;

    cwd = new char[len+1];

    strcpy(cwd, current->data);

    if (current->data[strlen(current->data)-1] != PATHSEP)
      strcat(cwd, PATHSEPSTR);
   }

   // If it's not AMBER||GROMACS, then "coordinates", "structure"
   // and "parameters" must be specified.
   if (!amberOn && !gromacsOn)
   { if (!opts.defined("coordinates"))
       NAMD_die("coordinates not found in the configuration file!");
     if (!opts.defined("structure"))
       NAMD_die("structure not found in the configuration file!");
     if (!opts.defined("parameters"))
       NAMD_die("parameters not found in the configuration file!");
   }
   
   // In any case, there should be either "coordinates" or
   // "ambercoor", but not both
   if (opts.defined("coordinates") && opts.defined("ambercoor"))
     NAMD_die("Cannot specify both coordinates and ambercoor!");
   if (!opts.defined("coordinates") && !opts.defined("ambercoor")
       && !opts.defined("grocoorfile"))
     NAMD_die("Coordinate file not found!");

   //  Make sure that both a temperature and a velocity PDB were
   //  specified
   if (opts.defined("temperature") &&
       (opts.defined("velocities") || opts.defined("binvelocities")) ) 
   {
      NAMD_die("Cannot specify both an initial temperature and a velocity file");
   }

   if (! opts.defined("auxFile")) {
     strcpy(auxFilename,outputFilename);
     strcat(auxFilename,".aux");
   }

   //  Check for frequencies
   if (dcdFrequency) {
     if (! opts.defined("dcdfile")) {
       strcpy(dcdFilename,outputFilename);
       strcat(dcdFilename,".dcd");
     }
   } else {
     dcdFilename[0] = STRINGNULL;
   }

   if (velDcdFrequency) {
     if (! opts.defined("veldcdfile")) {
       strcpy(velDcdFilename,outputFilename);
       strcat(velDcdFilename,".veldcd");
     }
   } else {
     velDcdFilename[0] = STRINGNULL;
   }
   
   if (xstFrequency) {
     if (! opts.defined("xstfile")) {
       strcpy(xstFilename,outputFilename);
       strcat(xstFilename,".xst");
     }
   } else {
     xstFilename[0] = STRINGNULL;
   }

   if (restartFrequency) {
     if (! opts.defined("restartname")) {
       strcpy(restartFilename,outputFilename);
       strcat(restartFilename,".restart");
     }
   } else {
     restartFilename[0] = STRINGNULL;
     restartSave = FALSE;
     binaryRestart = FALSE;
   }


   if (!amberOn)
   { //****** BEGIN CHARMM/XPLOR type changes
     //// set default
     if (!paraTypeXplorOn && !paraTypeCharmmOn) 
     {
       paraTypeXplorOn = TRUE;
     }
     //// make sure that there is just one type of input parameters specified
     if (paraTypeXplorOn && paraTypeCharmmOn) 
     {
       NAMD_die("Please specify either XPLOR or CHARMM format for parameters!");
     }
     //****** END CHARMM/XPLOR type changes
   }

   
   //  If minimization isn't on, must have a temp or velocity
   if (!(minimizeOn||minimizeCGOn) && !opts.defined("temperature") && 
       !opts.defined("velocities") && !opts.defined("binvelocities") ) 
   {
      NAMD_die("Must have either an initial temperature or a velocity file");
   }

   if (minimizeOn||minimizeCGOn) { initialTemp = 0.0; }
   if (opts.defined("velocities") || opts.defined("binvelocities") )
   {
  initialTemp = -1.0;
   }

   ///// periodic cell parameters

   if ( opts.defined("extendedSystem") )
   {
     current = config->find("extendedSystem");

     iout << iINFO << "EXTENDED SYSTEM FILE   "
        << current->data << "\n" << endi;

     ifstream xscFile(current->data);
     if ( ! xscFile ) NAMD_die("Unable to open extended system file.\n");

     char labels[1024];
     do {
       if ( ! xscFile ) NAMD_die("Error reading extended system file.\n");
       xscFile.getline(labels,1023);
     } while ( strncmp(labels,"#$LABELS ",9) );

     int a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z;
     a_x = a_y = a_z = b_x = b_y = b_z = c_x = c_y = c_z = -1;
     int o_x, o_y, o_z, s_u, s_v, s_w, s_x, s_y, s_z;
     o_x = o_y = o_z = s_u = s_v = s_w = s_x = s_y = s_z = -1;

     int pos = 0;
     char *l_i = labels + 8;
     while ( *l_i ) {
       if ( *l_i == ' ' ) { ++l_i; continue; }
       char *l_i2;
       for ( l_i2 = l_i; *l_i2 && *l_i2 != ' '; ++l_i2 );
       if ( (l_i2 - l_i) == 3 && (l_i[1] == '_') ) {
	 if (l_i[0] == 'a' && l_i[2] == 'x') a_x = pos;
	 if (l_i[0] == 'a' && l_i[2] == 'y') a_y = pos;
	 if (l_i[0] == 'a' && l_i[2] == 'z') a_z = pos;
	 if (l_i[0] == 'b' && l_i[2] == 'x') b_x = pos;
	 if (l_i[0] == 'b' && l_i[2] == 'y') b_y = pos;
	 if (l_i[0] == 'b' && l_i[2] == 'z') b_z = pos;
	 if (l_i[0] == 'c' && l_i[2] == 'x') c_x = pos;
	 if (l_i[0] == 'c' && l_i[2] == 'y') c_y = pos;
	 if (l_i[0] == 'c' && l_i[2] == 'z') c_z = pos;
	 if (l_i[0] == 'o' && l_i[2] == 'x') o_x = pos;
	 if (l_i[0] == 'o' && l_i[2] == 'y') o_y = pos;
	 if (l_i[0] == 'o' && l_i[2] == 'z') o_z = pos;
	 if (l_i[0] == 's' && l_i[2] == 'u') s_u = pos;
	 if (l_i[0] == 's' && l_i[2] == 'v') s_v = pos;
	 if (l_i[0] == 's' && l_i[2] == 'w') s_w = pos;
	 if (l_i[0] == 's' && l_i[2] == 'x') s_x = pos;
	 if (l_i[0] == 's' && l_i[2] == 'y') s_y = pos;
	 if (l_i[0] == 's' && l_i[2] == 'z') s_z = pos;
       }
       ++pos;
       l_i = l_i2;
     }
     int numpos = pos;

     for ( pos = 0; pos < numpos; ++pos ) {
       double tmp;
       xscFile >> tmp;
       if ( ! xscFile ) NAMD_die("Error reading extended system file.\n");
       if ( pos == a_x ) cellBasisVector1.x = tmp;
       if ( pos == a_y ) cellBasisVector1.y = tmp;
       if ( pos == a_z ) cellBasisVector1.z = tmp;
       if ( pos == b_x ) cellBasisVector2.x = tmp;
       if ( pos == b_y ) cellBasisVector2.y = tmp;
       if ( pos == b_z ) cellBasisVector2.z = tmp;
       if ( pos == c_x ) cellBasisVector3.x = tmp;
       if ( pos == c_y ) cellBasisVector3.y = tmp;
       if ( pos == c_z ) cellBasisVector3.z = tmp;
       if ( pos == o_x ) cellOrigin.x = tmp;
       if ( pos == o_y ) cellOrigin.y = tmp;
       if ( pos == o_z ) cellOrigin.z = tmp;
       if ( pos == s_u ) strainRate2.x = tmp;
       if ( pos == s_v ) strainRate2.y = tmp;
       if ( pos == s_w ) strainRate2.z = tmp;
       if ( pos == s_x ) strainRate.x = tmp;
       if ( pos == s_y ) strainRate.y = tmp;
       if ( pos == s_z ) strainRate.z = tmp;
     }

   }

   if ( cellBasisVector3.length2() && ! cellBasisVector2.length2() ) {
     NAMD_die("Used cellBasisVector3 without cellBasisVector2!");
   }

   if ( cellBasisVector2.length2() && ! cellBasisVector1.length2() ) {
     NAMD_die("Used cellBasisVector2 without cellBasisVector1!");
   }

   if ( cellOrigin.length2() && ! cellBasisVector1.length2() ) {
     NAMD_die("Used cellOrigin without cellBasisVector1!");
   }

   lattice.set(cellBasisVector1,cellBasisVector2,cellBasisVector3,cellOrigin);

   if (! opts.defined("DCDunitcell")) {
      dcdUnitCell = lattice.a_p() && lattice.b_p() && lattice.c_p();
   }

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
       splitPatch = SPLIT_PATCH_HYDROGEN;
     else
     {
       char err_msg[129];
       sprintf(err_msg, 
          "Illegal value '%s' for 'splitPatch' in configuration file", 
       s);
       NAMD_die(err_msg);
     }
   }

   ///// exclude stuff
   opts.get("exclude", s);

   if (!strcasecmp(s, "none"))
   {
      exclude = NONE;
      splitPatch = SPLIT_PATCH_POSITION;
   }
   else if (!strcasecmp(s, "1-2"))
   {
      exclude = ONETWO;
      splitPatch = SPLIT_PATCH_POSITION;
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

   if (scale14 != 1.0 && exclude != SCALED14)
   {
      iout << iWARN << "Exclude is not scaled1-4; 1-4scaling ignored.\n" << endi;
   }

   //  Get multiple timestep integration scheme
   if (!opts.defined("MTSAlgorithm"))
   {
  MTSAlgorithm = VERLETI;
   }
   else
   {
  opts.get("MTSAlgorithm", s);

  if (!strcasecmp(s, "naive"))
  {
    MTSAlgorithm = NAIVE;
  }
  if (!strcasecmp(s, "constant"))
  {
    MTSAlgorithm = NAIVE;
  }
  else if (!strcasecmp(s, "impulse"))
  {
    MTSAlgorithm = VERLETI;
  }
  else if (!strcasecmp(s, "verleti"))
  {
    MTSAlgorithm = VERLETI;
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
  longSplitting = C1;
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

   patchDimension = pairlistDist;

   if ( splitPatch == SPLIT_PATCH_HYDROGEN ) {
     patchDimension += hgroupCutoff;
   }

   BigReal defaultMargin = 0.0;
   if (berendsenPressureOn || langevinPistonOn) {
      defaultMargin = ( useFlexibleCell ? 0.06 : 0.03 ) * patchDimension;
   }
   if ( margin == XXXBIGREAL ) {
     margin = defaultMargin;
   }
   if ( defaultMargin != 0.0 && margin == 0.0 ) {
     margin = defaultMargin;
     iout << iWARN << "ALWAYS USE NON-ZERO MARGIN WITH CONSTANT PRESSURE!\n";
     iout << iWARN << "CHANGING MARGIN FROM 0 to " << margin << "\n" << endi;
   }

   patchDimension += margin;

   //  Turn on global integration if not explicitly specified

   if ( dihedralOn ) globalOn = TRUE;

   //  Make sure modes don't conflict
   if ((minimizeOn||minimizeCGOn) && langevinOn) 
   {
      NAMD_die("Minimization and Langevin dynamics are mutually exclusive dynamics modes");
   }

   if ((minimizeOn||minimizeCGOn) && tCoupleOn) 
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

   if (globalOn && CkNumPes() > 1)
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

   if (splitPatch == SPLIT_PATCH_POSITION && mollyOn )
   {
      NAMD_die("splitPatch hydrogen is required for MOLLY");
   }

   if (splitPatch == SPLIT_PATCH_POSITION && rigidBonds != RIGID_NONE)
   {
      NAMD_die("splitPatch hydrogen is required for rigidBonds");
   }

   //  Set the default value for the maximum movement parameter
   //  for minimization
   if (minimizeOn && (maximumMove == 0.0)) 
   {
      maximumMove = 0.75 * pairlistDist/stepsPerCycle;
   }

   if (langevinOn) {
     if ( ! opts.defined("langevinDamping") ) langevinDamping = 0.0;
     if ( ! opts.defined("langevinHydrogen") ) langevinHydrogen = TRUE;
     if ( (opts.defined("langevinDamping") || opts.defined("langevinHydrogen"))
       && (opts.defined("langevinFile") || opts.defined("langevinCol")) )
       NAMD_die("To specify Langevin dynamics parameters, use either langevinDamping and langevinHydrogen or langevinFile and langevinCol.  Do not combine them.");
     if ( opts.defined("langevinHydrogen") && langevinDamping == 0.0 )
       NAMD_die("langevinHydrogen requires langevinDamping to be set.");
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
    NAMD_die("Must give a rescale freqency if rescaleTemp is given");
  }
   }

   if ((minimizeOn||minimizeCGOn) && rescaleFreq > 0) 
   {
    NAMD_die("Minimization and temperature rescaling are mutually exclusive dynamics modes");
   }

   if (opts.defined("reassignFreq"))
   {
  if (!opts.defined("reassignTemp"))
  {
    if (opts.defined("temperature"))
    {
      reassignTemp = initialTemp;
    }
    else
    {
      NAMD_die("Must give a reassign temperature if reassignFreq is defined");
    }
  }
   }
   else
   {
  reassignFreq = -1;
  reassignTemp = 0.0;
   }

   if (opts.defined("reassignTemp"))
   {
  if (!opts.defined("reassignFreq"))
  {
    NAMD_die("Must give a reassignment freqency if reassignTemp is given");
  }
   }

   if (opts.defined("reassignIncr"))
   {
  if (!opts.defined("reassignFreq"))
  {
    NAMD_die("Must give a reassignment freqency if reassignIncr is given");
  }
   }
   else
   {
  reassignIncr = 0.0;
   }

   if (opts.defined("reassignHold"))
   {
  if (!opts.defined("reassignIncr"))
  {
    NAMD_die("Must give a reassignment increment if reassignHold is given");
  }
   }
   else
   {
  reassignHold = 0.0;
   }

   if ((minimizeOn||minimizeCGOn) && reassignFreq > 0) 
   {
    NAMD_die("Minimization and temperature reassignment are mutually exclusive dynamics modes");
   }

   if (!opts.defined("seed")) 
   {
      randomSeed = (unsigned int) time(NULL);
   }

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy

   if (fepOn) {
     if	     (rescaleFreq > 0) 	fepTemp = rescaleTemp;
     else if (reassignFreq > 0)	fepTemp = reassignTemp;
     else if (langevinOn) 	fepTemp = langevinTemp;
     else if (tCoupleOn) 	fepTemp = tCoupleTemp;
     else NAMD_die("Alchemical FEP can be performed only in constant temperature simulations");

     if (reassignFreq > 0 && reassignIncr != 0)
	NAMD_die("reassignIncr cannot be used in alchemical FEP runs");

     if (lambda < 0.0 || lambda > 1.0 || lambda2 < 0.0 || lambda2 > 1.0)
        NAMD_die("lambda values should be in the range [0.0, 1.0]");
  
     if (!opts.defined("fepoutfile")) {
       strcpy(fepOutFile, outputFilename);
       strcat(fepOutFile, ".fep");
     }
   } else {
     lambda = lambda2 = 0;
     fepOutFile[0] = STRINGNULL;
   }
//fepe

   if ( fepOn && lesOn )
     NAMD_die("Sorry, combined LES and FEP is not implemented.");
   if ( lesOn && ( lesFactor < 1 || lesFactor > 15 ) ) {
     NAMD_die("lesFactor must be positive and less than 16");
   }
   if ((pairInteractionOn && fepOn) || (pairInteractionOn && lesOn)) 
     NAMD_die("Sorry, pair interactions may not be calculated when LES or FEP is enabled.");

   //  Set up load balancing variables
   if (opts.defined("ldbStrategy"))
   {
     //  Assign the load balancing strategy
     if (strcasecmp(loadStrategy, "none") == 0)
       ldbStrategy=LDBSTRAT_NONE;
     else if (strcasecmp(loadStrategy, "refineonly") == 0)
       ldbStrategy=LDBSTRAT_REFINEONLY;
     else if (strcasecmp(loadStrategy, "alg7") == 0)
       ldbStrategy=LDBSTRAT_ALG7;
     else if (strcasecmp(loadStrategy, "orb") == 0)
       ldbStrategy=LDBSTRAT_ALGORB;
     else if (strcasecmp(loadStrategy, "neighbor") == 0)
       ldbStrategy=LDBSTRAT_ALGNBOR;
     else if (strcasecmp(loadStrategy, "other") == 0)
       ldbStrategy=LDBSTRAT_OTHER;
     else
       NAMD_die("Unknown ldbStrategy selected");
   } else {
     ldbStrategy=LDBSTRAT_OTHER;
   }

  if (!opts.defined("ldbPeriod"))
  {
    ldbPeriod=200*stepsPerCycle;
  }

  //  Set default values
  if (!opts.defined("firstLdbStep"))
  {
    firstLdbStep=5*stepsPerCycle;
  }

  if (ldbPeriod <= firstLdbStep) {
    NAMD_die("ldbPeriod must greater than firstLdbStep.");
  }

   if (N < firstTimestep) { N = firstTimestep; }

   if ( (firstTimestep%stepsPerCycle) != 0)
   {
  NAMD_die("First timestep must be a multiple of stepsPerCycle!!");
   }

   //  Make sure only one full electrostatics algorithm is selected
   {
     int i = 0;
     if ( FMAOn ) ++i;
     if ( PMEOn ) ++i;
     if ( fullDirectOn ) ++i;
     if ( i > 1 )
	NAMD_die("More than one full electrostatics algorithm selected!!!");
   }

   if (!opts.defined("ldbBackgroundScaling")) {
     ldbBackgroundScaling = 1.0;
   }
   if (!opts.defined("ldbPMEBackgroundScaling")) {
     ldbPMEBackgroundScaling = ldbBackgroundScaling;
   }
   if (!opts.defined("ldbHomeBackgroundScaling")) {
     ldbHomeBackgroundScaling = ldbBackgroundScaling;
   }
   if (!opts.defined("procsPerNode"))
   {
     procsPerNode=1;
   }
   if (!opts.defined("ldbUnloadRank"))
   {
     ldbUnloadRank = procsPerNode-1;
   }
   if (ldbUnloadRank >= procsPerNode)
	NAMD_die("Invalid unload rank of SMP node!!!");



   //  Check on PME parameters
   if (PMEOn) {  // idiot checking
     if ( lattice.volume() == 0. ) {
	NAMD_die("PME requires periodic boundary conditions.");
     }
     if ( PMEGridSizeX < 5 ) {
	NAMD_die("PMEGridSizeX (number of grid points) is very small.");
     }
     if ( PMEGridSizeY < 5 ) {
	NAMD_die("PMEGridSizeY (number of grid points) is very small.");
     }
     if ( PMEGridSizeZ < 5 ) {
	NAMD_die("PMEGridSizeZ (number of grid points) is very small.");
     }
     BigReal tolerance = PMETolerance;
     BigReal ewaldcof = 1.0;
     while ( erfc(ewaldcof*cutoff)/cutoff >= tolerance ) ewaldcof *= 2.0;
     BigReal ewaldcof_lo = 0.;
     BigReal ewaldcof_hi = ewaldcof;
     for ( int i = 0; i < 100; ++i ) {
       ewaldcof = 0.5 * ( ewaldcof_lo + ewaldcof_hi );
       if ( erfc(ewaldcof*cutoff)/cutoff >= tolerance ) {
         ewaldcof_lo = ewaldcof;
       } else {
         ewaldcof_hi = ewaldcof;
       }
     }
     PMEEwaldCoefficient = ewaldcof;
   } else {  // initialize anyway
     useDPME = 0;
     PMEGridSizeX = 0;
     PMEGridSizeY = 0;
     PMEGridSizeZ = 0;
     PMEEwaldCoefficient = 0;
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

   if ( (nonbondedFrequency > stepsPerCycle) || ( (stepsPerCycle % nonbondedFrequency) != 0) )
   {
     NAMD_die("stepsPerCycle must be a multiple of nonbondedFreq");
   }

   if (!FMAOn && !PMEOn && !fullDirectOn)
   {
     fullElectFrequency = 0;
   }
   else
   {
     if (!opts.defined("fullElectFrequency"))
     {
       if (opts.defined("fmaFrequency")) {
         iout << iWARN << "The parameter fmaFrequency has been renamed fullElectFrequency.\n" << endi;
         fullElectFrequency = fmaFrequency;
       } else {
         iout << iWARN << "The parameter fullElectFrequency now defaults to nonbondedFreq (" << nonbondedFrequency << ") rather than stepsPerCycle.\n" << endi;
         fullElectFrequency = nonbondedFrequency;
       }
     }
     else
     {
       if (opts.defined("fmaFrequency")) {
         iout << iWARN << "Ignoring redundant parameter fmaFrequency in favor of fullElectFrequency.\n" << endi;
       }
       if ( (fullElectFrequency > stepsPerCycle) || ( (stepsPerCycle % fullElectFrequency) != 0) )
       {
         NAMD_die("stepsPerCycle must be a multiple of fullElectFrequency");
       }
     }

     if ( (nonbondedFrequency > fullElectFrequency) || ( (fullElectFrequency % nonbondedFrequency) != 0) )
     {
       NAMD_die("fullElectFrequency must be a multiple of nonbondedFreq");
     }


     if (!opts.defined("fmaTheta"))
     fmaTheta=0.715;  /* Suggested by Duke developers */
   }

   if ( lesOn && ( FMAOn || useDPME || fullDirectOn ) ) {
     NAMD_die("Sorry, LES is only implemented for PME full electrostatics.");
   }
   if ( fepOn && ( FMAOn || useDPME || fullDirectOn ) ) {
     NAMD_die("Sorry, FEP is only implemented for PME full electrostatics.");
   }
   if ( pairInteractionOn && FMAOn ) {
     NAMD_die("Sorry, pairInteraction not implemented for FMA.");
   }
   if ( pairInteractionOn && useDPME ) {
     NAMD_die("Sorry, pairInteraction not implemented for DPME.");
   }
   if ( pairInteractionOn && fullDirectOn ) {
     NAMD_die("Sorry, pairInteraction not implemented for full direct electrostatics.");
   }
   if ( ! pairInteractionOn ) {
     pairInteractionSelf = 0;
   }
   if ( pairInteractionOn && !pairInteractionSelf && !config->find("pairInteractionGroup2")) 
     NAMD_die("pairInteractionGroup2 must be specified");

   if ( ! fixedAtomsOn ) {
     fixedAtomsForces = 0;
   }

   if (!opts.defined("constraints"))
   {
     constraintExp = 0;     
     constraintScaling = 1.0;     

     //****** BEGIN selective restraints (X,Y,Z) changes
     selectConstraintsOn = FALSE;
     //****** END selective restraints (X,Y,Z) changes
 
     //****** BEGIN moving constraints changes 
     movingConstraintsOn = FALSE;
     //****** END moving constraints changes 
     //****** BEGIN rotating constraints changes 
     rotConstraintsOn = FALSE;
    //****** END rotating constraints changes 
   } 
   //****** BEGIN rotating constraints changes 
   else {
     if (rotConstraintsOn) {
       rotConsAxis = rotConsAxis.unit();
     }
   }
   if(opts.defined("rotConstraints") 
      && opts.defined("movingConstraints")) {
     NAMD_die("Rotating and moving constraints are mutually exclusive!");
   }
   //****** END rotating constraints changes 

   //****** BEGIN selective restraints (X,Y,Z) changes
   if(opts.defined("selectConstraints") && !opts.defined("selectConstrX")
      && !opts.defined("selectConstrY") && !opts.defined("selectConstrZ")) {
     NAMD_die("selectConstraints was specified, but no Cartesian components were defined!");
   }
   if (!opts.defined("selectConstraints")) {
       constrXOn = FALSE;
       constrYOn = FALSE;
       constrZOn = FALSE;
   }
   //****** END selective restraints (X,Y,Z) changes


   //****** BEGIN SMD constraints changes 
   
   if (!opts.defined("SMD")) {     
     SMDOn = FALSE;
   }

   if (SMDOn) {
     // normalize direction
     if (SMDDir.length2() == 0) {
       NAMD_die("SMD direction vector must be non-zero");
     }
     else {
       SMDDir = SMDDir.unit();
     }

     if (SMDOutputFreq > 0 && SMDOutputFreq < stepsPerCycle
	 || SMDOutputFreq % stepsPerCycle != 0) {
       NAMD_die("SMDOutputFreq must be a multiple of stepsPerCycle");
     }
   }
     
   //****** END SMD constraints changes 
   
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

   if (!eFieldOn)
   {
        eField.x = 0.0;
        eField.y = 0.0;
        eField.z = 0.0;
   }

   if (!stirOn)
   { 
     stirFilename[0] = STRINGNULL;
     stirStartingTheta = 0.0;
     stirVel = 0.0;
     stirK = 0.0;
     stirAxis.x = 0.0;
     stirAxis.y = 0.0;
     stirAxis.z = 0.0;
     stirPivot.x = 0.0;
     stirPivot.y = 0.0;
     stirPivot.z = 0.0;
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

   // If we're doing pair interaction, set 
   // outputEnergies to 1 to make NAMD not die (the other nonbonded code paths 
   // aren't defined when these options are enabled), and set nonbondedFreq to 
   // 1 to avoid getting erroneous output.  Warn the user of what we're doing.
   if (pairInteractionOn) {
	   if (outputEnergies != 1) {
		   iout << iWARN << "Setting outputEnergies to 1 due to\n";
		   iout << iWARN << "pairInteraction calculations\n" << endi;
		   outputEnergies  = 1;
	   }
   }
   if (pairInteractionOn || pressureProfileOn) {
	   if (nonbondedFrequency != 1) {
		   iout << iWARN << "Setting nonbondedFreq to 1 due to\n";
		   iout << iWARN << "pairInteraction or pressure profile calculations\n" << endi;
	   }
   }

   // print timing at a reasonable interval by default
   if (!opts.defined("outputTiming"))
   {
      outputTiming = firstLdbStep;
      int ot2 = 10 * outputEnergies;
      if ( outputTiming < ot2 ) outputTiming = ot2;
   }

}

void SimParameters::print_config(ParseOptions &opts, ConfigList *config, char *&cwd) {

   StringList *current; //  Pointer to config option list

   //  Now that we have read everything, print it out so that
   //  the user knows what is going on
   iout << iINFO << "SIMULATION PARAMETERS:\n";
   iout << iINFO << "TIMESTEP               " << dt << "\n" << endi;
   iout << iINFO << "NUMBER OF STEPS        " << N << "\n";
   iout << iINFO << "STEPS PER CYCLE        " << stepsPerCycle << "\n";
   iout << endi;

   if ( lattice.a_p() || lattice.b_p() || lattice.c_p() ) {
     if ( lattice.a_p() )
       iout << iINFO << "PERIODIC CELL BASIS 1  " << lattice.a() << "\n";
     if ( lattice.b_p() )
       iout << iINFO << "PERIODIC CELL BASIS 2  " << lattice.b() << "\n";
     if ( lattice.c_p() )
       iout << iINFO << "PERIODIC CELL BASIS 3  " << lattice.c() << "\n";
     iout << iINFO << "PERIODIC CELL CENTER   " << lattice.origin() << "\n";
     if (wrapWater) {
       iout << iINFO << "WRAPPING WATERS AROUND PERIODIC BOUNDARIES ON OUTPUT.\n";
     }
     if (wrapAll) {
       iout << iINFO << "WRAPPING ALL CLUSTERS AROUND PERIODIC BOUNDARIES ON OUTPUT.\n";
     }
     if (wrapNearest) {
       iout << iINFO << "WRAPPING TO IMAGE NEAREST TO PERIODIC CELL CENTER.\n";
     }
     iout << endi;
   }

   if (ldbStrategy==LDBSTRAT_NONE)  {
     iout << iINFO << "LOAD BALANCE STRATEGY  none\n" << endi;
   } else {
     if (ldbStrategy==LDBSTRAT_REFINEONLY) {
       iout << iINFO << "LOAD BALANCE STRATEGY  Refine-only\n";
     } else if (ldbStrategy==LDBSTRAT_ALG7)  {
       iout << iINFO << "LOAD BALANCE STRATEGY  Alg7\n";
     } else if (ldbStrategy==LDBSTRAT_ALGORB)  {
       iout << iINFO << "LOAD BALANCE STRATEGY  RecBisection\n";
     } else if (ldbStrategy==LDBSTRAT_ALGNBOR)  {
       iout << iINFO << "LOAD BALANCE STRATEGY  Neighborhood\n";
     } else if (ldbStrategy==LDBSTRAT_OTHER)  {
       iout << iINFO << "LOAD BALANCE STRATEGY  Other\n";
     }
     iout << iINFO << "LDB PERIOD             " << ldbPeriod << " steps\n";
     iout << iINFO << "FIRST LDB TIMESTEP     " << firstLdbStep << "\n";
     iout << iINFO << "LDB BACKGROUND SCALING " << ldbBackgroundScaling << "\n";
     iout << iINFO << "HOM BACKGROUND SCALING " << ldbHomeBackgroundScaling << "\n";
     if ( PMEOn ) {
       iout << iINFO << "PME BACKGROUND SCALING "
				<< ldbPMEBackgroundScaling << "\n";
     if ( ldbUnloadPME )
     iout << iINFO << "REMOVING LOAD FROM PME NODES" << "\n";
     }
     if ( ldbUnloadZero ) iout << iINFO << "REMOVING LOAD FROM NODE 0\n";
     iout << endi;
     if ( CkNumPes() > 64 || ( IMDon && CkNumPes() > 8 ) ) {
       noPatchesOnZero = TRUE;
     }
     if ( noPatchesOnZero ) iout << iINFO << "REMOVING PATCHES FROM NODE 0\n";
     iout << endi;
     if ( ldbUnloadRankZero ) iout << iINFO << "REMOVING LOAD FROM RANK 0\n";
     iout << endi;
   }

   iout << iINFO << "MAX SELF PARTITIONS    " << maxSelfPart << "\n"
        << iINFO << "MAX PAIR PARTITIONS    " << maxPairPart << "\n"
        << iINFO << "SELF PARTITION ATOMS   " << numAtomsSelf << "\n"
        << iINFO << "PAIR PARTITION ATOMS   " << numAtomsPair << "\n"
        << iINFO << "PAIR2 PARTITION ATOMS  " << numAtomsPair2 << "\n"
        << iINFO << "MIN ATOMS PER PATCH    " << minAtomsPerPatch << "\n"
        << endi;
   
   if (initialTemp < 0)
   {
  current = config->find("velocities");

  if (current == NULL)
  {
    current = config->find("binvelocities");
  }

  iout << iINFO << "VELOCITY FILE          " << current->data << "\n";
   }
   else
   {
  iout << iINFO << "INITIAL TEMPERATURE    " 
     << initialTemp << "\n";
   }
   iout << endi;

   iout << iINFO << "CENTER OF MASS MOVING INITIALLY? ";

   if (comMove)
   {
     iout << "YES\n";
   }
   else
   {
     iout << "NO\n";
   }
   iout << endi;

   if ( zeroMomentum ) {
     iout << iINFO << "REMOVING CENTER OF MASS DRIFT DURING SIMULATION";
     if ( zeroMomentumAlt ) iout << " (ALT METHOD)";
     iout << "\n" << endi;
   }

   iout << iINFO << "DIELECTRIC             " 
      << dielectric << "\n";

   if ( nonbondedScaling != 1.0 )
   {
     iout << iINFO << "NONBONDED SCALING    " << nonbondedScaling << "\n" << endi;
   }

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
   iout << endi;

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
     iout << iINFO << "DCD FIRST STEP         " 
        << ( firstTimestep + dcdFrequency ) << "\n";
     if ( dcdUnitCell ) {
       iout << iINFO << "DCD FILE WILL CONTAIN UNIT CELL DATA\n";
     }
   }
   else
   {
     iout << iINFO << "NO DCD TRAJECTORY OUTPUT\n";
   }
   iout << endi;
   
   if (xstFrequency > 0)
   {
     iout << iINFO << "XST FILENAME           " 
        << xstFilename << "\n";
     iout << iINFO << "XST FREQUENCY          " 
        << xstFrequency << "\n";
   }
   else
   {
     iout << iINFO << "NO EXTENDED SYSTEM TRAJECTORY OUTPUT\n";
   }
   iout << endi;
   
   if (velDcdFrequency > 0)
   {
     iout << iINFO << "VELOCITY DCD FILENAME  " 
        << velDcdFilename << "\n";
     iout << iINFO << "VELOCITY DCD FREQUENCY " 
        << velDcdFrequency << "\n";
     iout << iINFO << "VELOCITY DCD FIRST STEP         " 
        << ( firstTimestep + velDcdFrequency ) << "\n";
   }
   else
   {
     iout << iINFO << "NO VELOCITY DCD OUTPUT\n";
   }
   iout << endi;
   
   iout << iINFO << "OUTPUT FILENAME        " 
      << outputFilename << "\n" << endi;
   if (binaryOutput)
   {
     iout << iINFO << "BINARY OUTPUT FILES WILL BE USED\n" << endi;
   }

   if (! restartFrequency)
   {
     iout << iINFO << "NO RESTART FILE\n";
   }
   else
   {
     iout << iINFO << "RESTART FILENAME       "
        << restartFilename << "\n";
     iout << iINFO << "RESTART FREQUENCY      " 
        << restartFrequency << "\n";
  if (restartSave) {
    iout << iINFO << "RESTART FILES WILL NOT BE OVERWRITTEN\n";
  }

  if (binaryRestart)
  {
    iout << iINFO << "BINARY RESTART FILES WILL BE USED\n";
  }
   }
   iout << endi;
   
   if (switchingActive)
   {
      iout << iINFO << "SWITCHING ACTIVE\n";
      iout << iINFO << "SWITCHING ON           "
               << switchingDist << "\n";
      iout << iINFO << "SWITCHING OFF          "
               << cutoff << "\n";
      if ( elecswitchDist != switchingDist ||
           vdwswitchDist  != switchingDist ||
           eleccutoff     != cutoff        ||
           vdwcutoff      != cutoff         ) {
        NAMD_die("Separate elect and vdw switching parameters not supported\n");
      }
    /*
      iout << iINFO << "E-SWITCHING ON         "
               << elecswitchDist << "\n";
      iout << iINFO << "E-SWITCHING OFF        "
               << eleccutoff << "\n";
      iout << iINFO << "VDW-SWITCHING ON       "
               << vdwswitchDist << "\n";
      iout << iINFO << "VDW-SWITCHING OFF      "
               << vdwcutoff << "\n";
    */
   }
   else
   {
      iout << iINFO << "CUTOFF                 " 
         << cutoff << "\n";
   }
   iout << iINFO << "PAIRLIST DISTANCE      " << pairlistDist << "\n";
   iout << iINFO << "PAIRLIST SHRINK RATE   " << pairlistShrink << "\n";
   iout << iINFO << "PAIRLIST GROW RATE     " << pairlistGrow << "\n";
   iout << iINFO << "PAIRLIST TRIGGER       " << pairlistTrigger << "\n";
   iout << iINFO << "PAIRLISTS PER CYCLE    " << pairlistsPerCycle << "\n";
   if ( outputPairlists )
     iout << iINFO << "PAIRLIST OUTPUT STEPS  " << outputPairlists << "\n";
   iout << endi;

   if ( pairlistMinProcs > 1 )
     iout << iINFO << "REQUIRING " << pairlistMinProcs << " PROCESSORS FOR PAIRLISTS\n";
   usePairlists = ( CkNumPes() >= pairlistMinProcs );
   iout << iINFO << "PAIRLISTS " << ( usePairlists ? "ENABLED" : "DISABLED" )
							<< "\n" << endi;

   iout << iINFO << "MARGIN                 " << margin << "\n";

   if ( splitPatch == SPLIT_PATCH_HYDROGEN ) {
      iout << iINFO << "HYDROGEN GROUP CUTOFF  " << hgroupCutoff << "\n";
   }
   
   iout << iINFO << "PATCH DIMENSION        "
            << patchDimension << "\n";

   iout << endi;

   if (outputEnergies != 1)
   {
      iout << iINFO << "ENERGY OUTPUT STEPS    "
         << outputEnergies << "\n";
      iout << endi;
   }

   if (mergeCrossterms) {
      iout << iINFO << "CROSSTERM ENERGY INCLUDED IN DIHEDRAL\n" << endi;
   }
   
   if (outputMomenta != 0)
   {
      iout << iINFO << "MOMENTUM OUTPUT STEPS  "
         << outputMomenta << "\n";
      iout << endi;
   }
   
   if (outputTiming != 0)
   {
      iout << iINFO << "TIMING OUTPUT STEPS    "
         << outputTiming << "\n";
      iout << endi;
   }
   
   if (outputPressure != 0)
   {
      iout << iINFO << "PRESSURE OUTPUT STEPS  "
         << outputPressure << "\n";
      iout << endi;
   }
   
   if (fixedAtomsOn)
   {
      iout << iINFO << "FIXED ATOMS ACTIVE\n";
      if ( fixedAtomsForces )
	iout << iINFO << "FORCES BETWEEN FIXED ATOMS ARE CALCULATED\n";
      iout << endi;
   }

   if (constraintsOn)
   {
      iout << iINFO << "HARMONIC CONSTRAINTS ACTIVE\n";

      iout << iINFO << "HARMONIC CONS EXP      "
         << constraintExp << "\n";

      if (constraintScaling != 1.0) {
        iout << iINFO << "HARMONIC CONS SCALING  "
         << constraintScaling << "\n";
      }

      //****** BEGIN selective restraints (X,Y,Z) changes 

      if (selectConstraintsOn) {
	iout << iINFO << "SELECTED CARTESIAN COMPONENTS OF HARMONIC RESTRAINTS ACTIVE\n";

        if (constrXOn)
	iout << iINFO << "RESTRAINING X-COMPONENTS OF CARTESIAN COORDINATES!\n";

        if (constrYOn)
	iout << iINFO << "RESTRAINING Y-COMPONENTS OF CARTESIAN COORDINATES!\n";

        if (constrZOn)
	iout << iINFO << "RESTRAINING Z-COMPONENTS OF CARTESIAN COORDINATES!\n";
      }
      iout << endi;
      //****** END selective restraints (X,Y,Z) changes 

      //****** BEGIN moving constraints changes 

      if (movingConstraintsOn) {
	iout << iINFO << "MOVING HARMONIC CONSTRAINTS ACTIVE\n";

	iout << iINFO << "MOVING CONSTRAINT VELOCITY    "
	     << movingConsVel << " ANGSTROM/TIMESTEP\n";
	
	iout << iINFO << "ALL CONSTRAINED ATOMS WILL MOVE\n";
      }
      //****** END moving constraints changes 
      iout << endi;

      //****** BEGIN rotating constraints changes 

      if (rotConstraintsOn) {
	iout << iINFO << "ROTATING HARMONIC CONSTRAINTS ACTIVE\n";

	iout << iINFO << "AXIS OF ROTATION    "
	     << rotConsAxis << "\n";
	
	iout << iINFO << "PIVOT OF ROTATION   "
	     << rotConsPivot << "\n";

	iout << iINFO << "ROTATING CONSTRAINT VELOCITY    "
	     << rotConsVel << " DEGREES/TIMESTEP\n";
      }
      iout << endi;
      //****** END rotating constraints changes 
   }


   // moving drag
   if (movDragOn) {
     iout << iINFO << "MOVING DRAG ACTIVE.\n";
     
     iout << iINFO << "MOVING DRAG MAIN PDB FILE "
	  << movDragFile << "\n";
     
     iout << iINFO << "MOVING DRAG GLOBAL VELOCITY (A/step) "
	  << movDragGlobVel << "\n";
     
     iout << iINFO << "MOVING DRAG LINEAR VELOCITY FILE " 
	  << movDragVelFile << "\n";
     
     iout << endi;
   }
   
   // rotating drag
   if (rotDragOn) {
     iout << iINFO << "ROTATING DRAG ACTIVE.\n";
     
     iout << iINFO << "ROTATING DRAG MAIN PDB FILE "
	  << rotDragFile << "\n";
     
     iout << iINFO << "ROTATING DRAG AXIS FILE " 
	  << rotDragAxisFile << "\n";
     
     iout << iINFO << "ROTATING DRAG PIVOT POINT FILE " 
	  << rotDragPivotFile << "\n";
     
     iout << iINFO << "ROTATING DRAG GLOBAL ANGULAR VELOCITY (deg/step) "
	  << rotDragGlobVel << "\n";
     
     iout << iINFO << "ROTATING DRAG ANGULAR VELOCITY FILE " 
	  << rotDragVelFile << "\n";
     
     iout << endi;
   }
   

   // "constant" torque
   if (consTorqueOn) {
     iout << iINFO << "\"CONSTANT\" TORQUE ACTIVE.\n";
     
     iout << iINFO << "\"CONSTANT\" TORQUE MAIN PDB FILE "
	  << consTorqueFile << "\n";
     
     iout << iINFO << "\"CONSTANT\" TORQUE AXIS FILE " 
	  << consTorqueAxisFile << "\n";
     
     iout << iINFO << "\"CONSTANT\" TORQUE PIVOT POINT FILE " 
	  << consTorquePivotFile << "\n";
     
     iout << iINFO << "\"CONSTANT\" TORQUE GLOBAL VALUE (Kcal/(mol*A^2)) "
	  << consTorqueGlobVal << "\n";
     
     iout << iINFO << "\"CONSTANT\" TORQUE DACTORS FILE " 
	  << consTorqueValFile << "\n";
     
     iout << endi;
   }
   

   //****** BEGIN SMD constraints changes 
   
   if (SMDOn) {
     iout << iINFO << "SMD ACTIVE\n";
     
     iout << iINFO << "SMD VELOCITY    "
	  << SMDVel << " ANGSTROM/TIMESTEP\n";
	
     iout << iINFO << "SMD DIRECTION   "
	  << SMDDir << "\n";
 
     iout << iINFO << "SMD K   " 
          << SMDk << "\n";

     iout << iINFO << "SMD OUTPUT FREQUENCY   "
	  << SMDOutputFreq << " TIMESTEPS\n";
    
     iout << iINFO << "SMD FILE " << SMDFile << "\n"; 

     iout << endi;
   }
   
   //****** END SMD constraints changes 

   if (TMDOn) {
     iout << iINFO << "TMD ACTIVE BETWEEN STEPS " << TMDFirstStep 
          << " and " << TMDLastStep << "\n";
     iout << iINFO << "TMD K  " << TMDk << "\n";
     iout << iINFO << "TMD FILE  " << TMDFile << "\n";
     iout << iINFO << "TMD OUTPUT FREQUENCY  " << TMDOutputFreq << "\n";
     if (TMDInitialRMSD) {
       iout << iINFO << "TMD TARGET RMSD AT FIRST STEP  " << TMDInitialRMSD << "\n";
     } else {
       iout << iINFO << "TMD TARGET RMSD AT FIRST STEP COMPUTED FROM INITIAL COORDINATES\n";
     }
     iout << iINFO << "TMD TARGET RMSD AT FINAL STEP  " << TMDFinalRMSD << "\n";
     iout << endi;
   }

   
//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
//  Alchemical FEP status

//   current = config->find("fepOutFile");
   if (fepOn)
   {
     iout << iINFO << "ALCHEMICAL FEP ON\n";
     iout << iINFO << "CURRENT LAMBDA VALUE     "
          << lambda << "\n";
     iout << iINFO << "COMPARISON LAMBDA VALUE  "
          << lambda2 << "\n";
   }

//fepe

   if ( lesOn ) {
     iout << iINFO << "LOCALLY ENHANCED SAMPLING ACTIVE\n";
     iout << iINFO << "LOCAL ENHANCEMENT FACTOR IS "
          << lesFactor << "\n";
     if ( lesReduceTemp ) iout << iINFO
       << "SCALING ENHANCED ATOM TEMPERATURE BY 1/" << lesFactor << "\n";
     if ( lesReduceMass ) iout << iINFO
       << "SCALING ENHANCED ATOM MASS BY 1/" << lesFactor << "\n";
   }
   
   if ( pairInteractionOn ) {
     iout << iINFO << "PAIR INTERACTION CALCULATIONS ACTIVE\n";
     iout << iINFO << "USING FLAG " << pairInteractionGroup1 
          << " FOR GROUP 1\n";
     if (pairInteractionSelf) {
       iout << iINFO << "COMPUTING ONLY SELF INTERACTIONS FOR GROUP 1 ATOMS\n";
     } else {
       iout << iINFO << "USING FLAG " << pairInteractionGroup2 
            << " FOR GROUP 2\n";
     }
   }

   if (consForceOn)
     iout << iINFO << "CONSTANT FORCE ACTIVE\n";

   // external command forces

   if (extForcesOn) {
     iout << iINFO << "EXTERNAL COMMAND FORCES ACTIVE\n";
     iout << iINFO << "EXT FORCES COMMAND: " << extForcesCommand << "\n";
     iout << iINFO << "EXT COORD FILENAME: " << extCoordFilename << "\n";
     iout << iINFO << "EXT FORCE FILENAME: " << extForceFilename << "\n";
     iout << endi;
   }

   tclBCScript = 0;
   if (tclBCOn) {
     iout << iINFO << "TCL BOUNDARY FORCES ACTIVE\n";
     current = config->find("tclBCScript");
     if ( current ) {
       tclBCScript = current->data;
       iout << iINFO << "TCL BOUNDARY FORCES SCRIPT   " << current->data << "\n";
     }
       iout << iINFO << "TCL BOUNDARY FORCES ARGS     " << tclBCArgs << "\n";
     iout << endi;
   }
   
   // Global forces configuration

   globalForcesOn = ( tclForcesOn || freeEnergyOn || miscForcesOn ||
                      (IMDon) || SMDOn || TMDOn);

   if (tclForcesOn)
   {
     iout << iINFO << "TCL GLOBAL FORCES ACTIVE\n";

     current = config->find("tclForcesScript");

     for ( ; current; current = current->next ) {

     if ( strstr(current->data,"\n") ) {
       iout << iINFO << "TCL GLOBAL FORCES SCRIPT INLINED IN CONFIG FILE\n";
       continue;
     }

     iout << iINFO << "TCL GLOBAL FORCES SCRIPT   " << current->data << "\n";

     }
     iout << endi;
   }

   if (miscForcesOn)
   {
     iout << iINFO << "MISC FORCES ACTIVE\n";

     current = config->find("miscForcesScript");

     for ( ; current; current = current->next ) {

     if ( strstr(current->data,"\n") ) {
       iout << iINFO << "MISC FORCES SCRIPT INLINED IN CONFIG FILE\n";
       continue;
     }

     iout << iINFO << "MISC FORCES SCRIPT   " << current->data << "\n";

     }
     iout << endi;
   }

   if (freeEnergyOn)
   {
     iout << iINFO << "FREE ENERGY PERTURBATION ACTIVE\n";

     current = config->find("freeEnergyConfig");

     for ( ; current; current = current->next ) {

     if ( strstr(current->data,"\n") ) {
       iout << iINFO << "FREE ENERGY PERTURBATION SCRIPT INLINED IN CONFIG FILE\n";
       continue;
     }

     iout << iINFO << "FREE ENERGY PERTURBATION SCRIPT   " << current->data << "\n";

     }
     iout << endi;
   }

   if (IMDon)
   {
     iout << iINFO << "INTERACTIVE MD ACTIVE\n";
     iout << iINFO << "INTERACTIVE MD PORT    " << IMDport << "\n";
     iout << iINFO << "INTERACTIVE MD FREQ    " << IMDfreq << "\n";
     if (IMDignore) {
        iout << iINFO << "INTERACTIVE MD WILL NOT INFLUENCE SIMULATION\n";
     } else {
       if (IMDwait) iout << iINFO << "WILL AWAIT INTERACTIVE MD CONNECTION\n";
     }
     iout << endi;
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
    iout << iINFO << "CYLINDER BOUNDARY CENTER(" << cylindricalCenter.x << ", "
             << cylindricalCenter.y << ", " << cylindricalCenter.z << ")\n";
    iout << endi;
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

      iout << iINFO << "SPHERE BOUNDARY CENTER(" << sphericalCenter.x << ", "
               << sphericalCenter.y << ", " << sphericalCenter.z << ")\n";
      iout << endi;
   }
   
   if (eFieldOn)
   {
      iout << iINFO << "ELECTRIC FIELD ACTIVE\n";
      
      iout << iINFO << "E-FIELD VECTOR         ("
         << eField.x << ", " << eField.y
         << ", " << eField.z << ")\n";
      iout << endi;
   }

      if (stirOn)
   {
      iout << iINFO << "STIRRING TORQUES ACTIVE\n";
      
      iout << iINFO << "STIR STARTING THETA   (deg)  "<< stirStartingTheta << "\n";
      iout << iINFO << "STIR ANGULAR VELOCITY (deg/ts)   " << stirVel <<"\n";
      iout << iINFO << "STIR FORCE HARMONIC SPRING CONSTANT "<< stirK << "\n";
      iout << iINFO << "STIR AXIS OF ROTATION (DIRECTION)      ("
         << stirAxis.x << ", " << stirAxis.y
         << ", " << stirAxis.z << ")\n";
	          iout << iINFO << "STIR PIVOT POINT (COORDINATE)           ("
         << stirPivot.x << ", " << stirPivot.y
         << ", " << stirPivot.z << ")\n";
      current = config->find("stirFilename");
		  
      iout << iINFO << "STIR ATOMS AND ORIGINAL POSITIONS FROM FILE    " <<current ->data << '\n';
      current = config->find("stirredAtomsCol");
      iout << iINFO <<"STIR FILE COLUMN " << current ->data << '\n';
      iout << endi;
   }

   if (langevinOn)
   {
      iout << iINFO << "LANGEVIN DYNAMICS ACTIVE\n";
      iout << iINFO << "LANGEVIN TEMPERATURE   "
         << langevinTemp << "\n";
      if (langevinDamping > 0.0) {
	iout << iINFO << "LANGEVIN DAMPING COEFFICIENT IS "
		<< langevinDamping << " INVERSE PS\n";
	if (langevinHydrogen)
		iout << iINFO << "LANGEVIN DYNAMICS APPLIED TO HYDROGENS\n";
	else
		iout << iINFO << "LANGEVIN DYNAMICS NOT APPLIED TO HYDROGENS\n";
      } else {
	iout << iINFO << "LANGEVIN DAMPING COEFFICIENTS DETERMINED FROM FILES\n";
        current = config->find("langevinFile");
	if ( current ) iout << iINFO << "LANGEVIN DAMPING FILE:  " <<
          current->data << "\n";
        else iout << iINFO << "LANGEVIN DAMPING FILE IS COORDINATE PDB\n";
        current = config->find("langevinCol");
	if ( current ) iout << iINFO << "LANGEVIN DAMPING COLUMN:  " <<
          current->data << "\n";
        else iout << iINFO << "LANGEVIN DAMPING COLUMN:  DEFAULT (4TH, O)\n";
      }
      iout << endi;
   }

   if (tCoupleOn)
   {
      iout << iINFO << "TEMPERATURE COUPLING ACTIVE\n";
      iout << iINFO << "COUPLING TEMPERATURE   "
         << tCoupleTemp << "\n";
      iout << endi;
   }

   if (minimizeOn)
   {
      iout << iINFO << "OLD STYLE MINIMIZATION ACTIVE\n";
      iout << endi;
   }

   if (minimizeCGOn)
   {
      iout << iINFO << "CONJUGATE GRADIENT MINIMIZATION ACTIVE\n";
      iout << iINFO << "LINE MINIMIZATION GOAL = " << minLineGoal << "\n";
      iout << iINFO << "BABY STEP SIZE = " << minBabyStep << "\n";
      iout << iINFO << "TINY STEP SIZE = " << minTinyStep << "\n";
      iout << endi;
   }

   if (maximumMove)
   {
      iout << iINFO << "MAXIMUM MOVEMENT       "
         << maximumMove << "\n";
      iout << endi;
   }

   if (rescaleFreq > 0)
   {
     iout << iINFO << "VELOCITY RESCALE FREQ  "
        << rescaleFreq << "\n";
     iout << iINFO << "VELOCITY RESCALE TEMP  "
        << rescaleTemp << "\n";
     iout << endi;
   }

   if (reassignFreq > 0)
   {
     iout << iINFO << "VELOCITY REASSIGNMENT FREQ  "
        << reassignFreq << "\n";
     iout << iINFO << "VELOCITY REASSIGNMENT TEMP  "
        << reassignTemp << "\n";
     if ( reassignIncr != 0. )
       iout << iINFO << "VELOCITY REASSIGNMENT INCR  "
        << reassignIncr << "\n";
     if ( reassignHold != 0. )
       iout << iINFO << "VELOCITY REASSIGNMENT HOLD  "
        << reassignHold << "\n";
     iout << endi;
   }

   if (berendsenPressureOn && langevinPistonOn)
   {
      NAMD_die("Multiple pressure control algorithms selected!\n");
   }

   if (excludeFromPressure) {
     iout << iINFO << "EXCLUDE FROM PRESSURE ACTIVE\n";
   }
   if (useConstantArea && useConstantRatio) {
     NAMD_die("useConstantArea and useConstantRatio are mutually exclusive.\n");
   }
   if (useConstantRatio && !useFlexibleCell) {
     NAMD_die("useConstantRatio requires useFlexibleCell.\n");
   }
   if (useConstantArea && surfaceTensionTarget) {
     NAMD_die("surfaceTensionTarget and useConstantArea are mutually exclusive.\n");
   }
   if (useConstantArea && !useFlexibleCell) {
     NAMD_die("useConstantArea requires useFlexibleCell.\n");
   }

   if (berendsenPressureOn || langevinPistonOn) {
     if (rigidBonds != RIGID_NONE && useGroupPressure == FALSE) {
       useGroupPressure = TRUE;
       iout << iWARN << "Option useGroupPressure is being enabled "
            << "due to pressure control with rigidBonds.\n" << endi;
     }
   }

   if (berendsenPressureOn)
   {
     if ( ! opts.defined("BerendsenPressureFreq") ) {
	berendsenPressureFreq = nonbondedFrequency;
	if ( fullElectFrequency )
		berendsenPressureFreq = fullElectFrequency;
     }
     if ( (berendsenPressureFreq % nonbondedFrequency) || ( fullElectFrequency
		&& (berendsenPressureFreq % fullElectFrequency) ) )
	NAMD_die("berendsenPressureFreq must be a multiple of both fullElectFrequency and nonbondedFrequency\n");
     iout << iINFO << "BERENDSEN PRESSURE COUPLING ACTIVE\n";
     iout << iINFO << "    TARGET PRESSURE IS "
        << berendsenPressureTarget << " BAR\n";
     iout << iINFO << "    COMPRESSIBILITY ESTIMATE IS "
        << berendsenPressureCompressibility << " BAR^(-1)\n";
     iout << iINFO << "    RELAXATION TIME IS "
        << berendsenPressureRelaxationTime << " FS\n";
     iout << iINFO << "    APPLIED EVERY "
        << berendsenPressureFreq << " STEPS\n";
     iout << iINFO << "    PRESSURE CONTROL IS "
	<< (useGroupPressure?"GROUP":"ATOM") << "-BASED\n";
     iout << endi;
     berendsenPressureTarget /= PRESSUREFACTOR;
     berendsenPressureCompressibility *= PRESSUREFACTOR;
   }

   if (langevinPistonOn)
   {
     iout << iINFO << "LANGEVIN PISTON PRESSURE CONTROL ACTIVE\n";
     iout << iINFO << "       TARGET PRESSURE IS "
        << langevinPistonTarget << " BAR\n";
     iout << iINFO << "    OSCILLATION PERIOD IS "
        << langevinPistonPeriod << " FS\n";
     iout << iINFO << "            DECAY TIME IS "
        << langevinPistonDecay << " FS\n";
     iout << iINFO << "    PISTON TEMPERATURE IS "
        << langevinPistonTemp << " K\n";
     iout << iINFO << "      PRESSURE CONTROL IS "
	<< (useGroupPressure?"GROUP":"ATOM") << "-BASED\n";
     iout << iINFO << "   INITIAL STRAIN RATE IS "
        << strainRate << "\n";
     iout << endi;
     langevinPistonTarget /= PRESSUREFACTOR;
   }

   if (berendsenPressureOn || langevinPistonOn) {
     iout << iINFO << "      CELL FLUCTUATION IS "
	    << (useFlexibleCell?"AN":"") << "ISOTROPIC\n";
     if (useConstantRatio) 
       iout << iINFO << "    SHAPE OF CELL IS CONSTRAINED IN X-Y PLANE\n";
     if (useConstantArea) 
       iout << iINFO << "    CONSTANT AREA PRESSURE CONTROL ACTIVE\n";
   }

   if (surfaceTensionTarget != 0)
   {
     iout << iINFO << "SURFACE TENSION CONTROL ACTIVE\n";
     iout << iINFO << "      TARGET SURFACE TENSION IS "
          << surfaceTensionTarget << " DYN/CM\n";
     iout << endi;
     // multiply by 100 to convert from dyn/cm to bar-Angstroms, then divide
     // by PRESSURE factor to convert bar to NAMD internal pressure units. 
     surfaceTensionTarget *= 100.0 / PRESSUREFACTOR;
   }

   if (pressureProfileOn) {
     if ((berendsenPressureOn || langevinPistonOn) && !dcdUnitCell) {
#if 1
       iout << iWARN << "Turning on dcdUnitCell so that trajectory files contain unit cell data.\n" << endi;
       dcdUnitCell = 1;
#else
       NAMD_die("Sorry, pressure profile not implemented for constant pressure.");
#endif
     }
     // if Ewald is on, only calculate Ewald
     if (pressureProfileEwaldOn)
       pressureProfileOn = 0;

     if (pressureProfileSlabs < 1) 
       NAMD_die("pressureProfileSlabs must be positive.");
     iout << iINFO << "PRESSURE PROFILE CALCULATIONS ACTIVE\n";
     iout << iINFO << "      NUMBER OF SLABS: " << pressureProfileSlabs << "\n";
     iout << iINFO << "      SLAB THICKNESS: " << cellBasisVector3.z / pressureProfileSlabs
                   << "\n";
     iout << iINFO << "      TIMESTEPS BETWEEN DATA OUTPUT: " 
                   << pressureProfileFreq << "\n";
     iout << iINFO << "      NUMBER OF ATOM TYPES: " << pressureProfileAtomTypes << "\n";
     iout << endi;
   } else {
     pressureProfileEwaldOn = 0;
     pressureProfileAtomTypes = 1;
   }

   if (FMAOn)
   {
     iout << iINFO << "FMA ACTIVE\n";
     iout << iINFO << "FMA THETA              "
        << fmaTheta << "\n";
     iout << endi;
   }

   FFTWWisdomString = 0;
   if (PMEOn)
   {
     iout << iINFO << "PARTICLE MESH EWALD (PME) ACTIVE\n";
     iout << iINFO << "PME TOLERANCE               "
	<< PMETolerance << "\n";
     iout << iINFO << "PME EWALD COEFFICIENT       "
	<< PMEEwaldCoefficient << "\n";
     iout << iINFO << "PME INTERPOLATION ORDER     "
	<< PMEInterpOrder << "\n";
     iout << iINFO << "PME GRID DIMENSIONS         "
	<< PMEGridSizeX << " "
	<< PMEGridSizeY << " "
	<< PMEGridSizeZ << "\n";
     if ( PMEBarrier ) {
       iout << iINFO << "PME BARRIER ENABLED\n";
     }
     iout << endi;
     if ( useDPME ) iout << iINFO << "USING OLD DPME CODE\n";
#ifdef NAMD_FFTW
     else if ( FFTWUseWisdom ) {  // handle FFTW wisdom
       if (! opts.defined("FFTWWisdomFile")) {
         strcpy(FFTWWisdomFile,"FFTW_NAMD_");
         strcat(FFTWWisdomFile,NAMD_VERSION);
	 strcat(FFTWWisdomFile,"_");
	 strcat(FFTWWisdomFile,NAMD_PLATFORM);
	 strcat(FFTWWisdomFile,".txt");
       }
       iout << iINFO << "Attempting to read FFTW data from "
		<< FFTWWisdomFile << "\n" << endi;
       FILE *wisdom_file = fopen(FFTWWisdomFile,"r");
       if ( wisdom_file ) {
	 fftw_import_wisdom_from_file(wisdom_file);
	 fclose(wisdom_file);
       }

       int nrp = 1;

       // rules based on work available
       int minslices = 1;
       int dimx = PMEGridSizeX;
       int nrpx = ( dimx + minslices - 1 ) / minslices;
       if ( nrpx > nrp ) nrp = nrpx;
       int dimy = PMEGridSizeY;
       int nrpy = ( dimy + minslices - 1 ) / minslices;
       if ( nrpy > nrp ) nrp = nrpy;

       // rules based on processors available
       int nrpp = CkNumPes();
       // if ( nrpp > 32 ) nrpp = 32;  // cap to limit messages
       if ( nrpp < nrp ) nrp = nrpp;

       // user override
       int nrps = PMEProcessors;
       if ( nrps > CkNumPes() ) nrps = CkNumPes();
       if ( nrps > 0 ) nrp = nrps;

       // make sure there aren't any totally empty processors
       int bx = ( dimx + nrp - 1 ) / nrp;
       int nrpbx = ( dimx + bx - 1 ) / bx;
       int by = ( dimy + nrp - 1 ) / nrp;
       int nrpby = ( dimy + by - 1 ) / by;
       nrp = ( nrpby > nrpbx ? nrpby : nrpbx );
       if ( bx != ( dimx + nrp - 1 ) / nrp )
         NAMD_bug("Error in selecting number of PME processors.");
       if ( by != ( dimy + nrp - 1 ) / nrp )
         NAMD_bug("Error in selecting number of PME processors.");

       // numGridPes = nrpbx;
       // numTransPes = nrpby;
       // numRecipPes = nrp;
       int block2 = (PMEGridSizeY + nrp - 1) / nrp;
       int block2_min = PMEGridSizeY % block2;
       if ( ! block2_min ) block2_min = block2;
       int dim3 = 2 * (PMEGridSizeZ/2 + 1);

       int n[3]; n[0] = PMEGridSizeX; n[1] = PMEGridSizeY; n[2] = PMEGridSizeZ;
       fftw_complex *work = new fftw_complex[n[0]];
       float *grid1 = new float[n[1]*dim3];
       float *grid2 = new float[n[0]*block2*dim3];
       iout << iINFO << "Optimizing 6 FFT steps.  1..." << endi;
       rfftwnd_destroy_plan( rfftwnd_create_plan_specific(
	 2, n+1, FFTW_REAL_TO_COMPLEX,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, grid1, 1, 0, 0) );
       iout << " 2..." << endi;
       fftw_destroy_plan( fftw_create_plan_specific(n[0], FFTW_REAL_TO_COMPLEX,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) grid2,
	 block2*dim3/2, work, 1) );
       iout << " 3..." << endi;
       fftw_destroy_plan( fftw_create_plan_specific(n[0], FFTW_REAL_TO_COMPLEX,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) grid2,
	 block2_min*dim3/2, work, 1) );
       iout << " 4..." << endi;
       fftw_destroy_plan( fftw_create_plan_specific(n[0], FFTW_COMPLEX_TO_REAL,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) grid2,
	 block2*dim3/2, work, 1) );
       iout << " 5..." << endi;
       fftw_destroy_plan( fftw_create_plan_specific(n[0], FFTW_COMPLEX_TO_REAL,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) grid2,
	 block2_min*dim3/2, work, 1) );
       iout << " 6..." << endi;
       rfftwnd_destroy_plan( rfftwnd_create_plan_specific(
	 2, n+1, FFTW_COMPLEX_TO_REAL,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, grid1, 1, 0, 0) );
       iout << "   Done.\n" << endi;
       delete [] work;

       iout << iINFO << "Writing FFTW data to "
		<< FFTWWisdomFile << "\n" << endi;
       wisdom_file = fopen(FFTWWisdomFile,"w");
       if ( wisdom_file ) {
	 fftw_export_wisdom_to_file(wisdom_file);
	 fclose(wisdom_file);
       }

       FFTWWisdomString = fftw_export_wisdom_to_string();
     }
#endif
     iout << endi;
   }

   if (fullDirectOn)
   {
     iout << iINFO << "DIRECT FULL ELECTROSTATIC CALCULATIONS ACTIVE\n";
     iout << endi;
   }

   if ( FMAOn || PMEOn || fullDirectOn )
   {
     iout << iINFO << "FULL ELECTROSTATIC EVALUATION FREQUENCY      "
	<< fullElectFrequency << "\n";
     iout << endi;
   }

  if (MTSAlgorithm == NAIVE)
  {
    iout << iINFO << "USING NAIVE (CONSTANT FORCE) MTS SCHEME.\n" << endi;
  }
  if (MTSAlgorithm == VERLETI )
  {
    iout << iINFO << "USING VERLET I (r-RESPA) MTS SCHEME.\n" << endi;
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

   iout << endi;

   if (mollyOn)
   {
     iout << iINFO << "SLOW FORCE MOLLIFICATION : \n";
     iout << iINFO << "         ERROR TOLERANCE : " << mollyTol << "\n";
     iout << iINFO << "          MAX ITERATIONS : " << mollyIter << "\n";
     iout << endi;
   }

   if (rigidBonds != RIGID_NONE)
   {
     iout << iINFO << "RIGID BONDS TO HYDROGEN : ";
     if (rigidBonds == RIGID_ALL)    iout << "ALL\n";
     if (rigidBonds == RIGID_WATER)  iout << "WATER\n";
     iout << iINFO << "        ERROR TOLERANCE : " << rigidTol << "\n";
     iout << iINFO << "         MAX ITERATIONS : " << rigidIter << "\n";
     if (useSettle) iout << iINFO << "RIGID WATER USING SETTLE ALGORITHM\n";
     iout << endi;
   }
   

   if (nonbondedFrequency != 1)
   {
     iout << iINFO << "NONBONDED FORCES EVALUATED EVERY " << nonbondedFrequency << " STEPS\n";
   }

   iout << iINFO << "RANDOM NUMBER SEED     "
      << randomSeed << "\n";

   iout << endi;

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

// If this is AMBER, then print AMBER options

   if (amberOn)
   { iout << iINFO << "Using AMBER format force field!\n";
     current = config->find("parmfile");
     iout << iINFO << "AMBER PARM FILE        " << current->data << '\n';
     if (opts.defined("coordinates"))
     { current = config->find("coordinates");
       iout << iINFO << "COORDINATE PDB         " << current->data << '\n';
     }
     else
     { current = config->find("ambercoor");
       iout << iINFO << "AMBER COORDINATE FILE  " << current->data << '\n';
     }
     if (readExclusions)
       iout << iINFO << "Exclusions will be read from PARM file!\n";
     else
       iout << iINFO << "Exclusions in PARM file will be ignored!\n";
     iout << iINFO << "SCNB (VDW SCALING)     " << vdwscale14 << "\n" << endi;
   }
   else if(gromacsOn)
   {
     iout << iINFO << "Using GROMACS format force field!\n";

     current = config->find("grotopfile");
     // it should be defined, but, just in case...
     if (current == NULL)
       NAMD_die("no GROMACS topology file defined!?");
     iout << iINFO << "GROMACS TOPO FILE        " << current->data << '\n';

     // XXX handle the two types of coordinates more gracefully
     current = config->find("grocoorfile");
     if (current == NULL) {
       current = config->find("coordinates");
       if (current == NULL) {
	 NAMD_die("no coordinate file defined!?");
       }
     }
     iout << iINFO << "GROMACS COOR FILE        " << current->data << '\n' 
	  << endi;

   }
   else {
     current = config->find("coordinates");
     iout << iINFO << "COORDINATE PDB         " << current->data << '\n' << endi;

     current = config->find("structure");

     iout << iINFO << "STRUCTURE FILE         " 
        << current->data << "\n" << endi;

     //****** BEGIN CHARMM/XPLOR type changes
     if (paraTypeXplorOn)
     {
       iout << iINFO << "PARAMETER file: XPLOR format! (default) \n" << endi;
     }
     else if (paraTypeCharmmOn)
     {
       iout << iINFO << "PARAMETER file: CHARMM format! \n" << endi;
     }
     //****** END CHARMM/XPLOR type changes

     current = config->find("parameters");

     while (current != NULL)
     {
       iout << iINFO << "PARAMETERS             " 
          << current->data << "\n" << endi;
       current = current->next;
     }
   }

     iout << iINFO << "USING " <<
        ( vdwGeometricSigma ? "GEOMETRIC" : "ARITHMETIC" ) <<
        " MEAN TO COMBINE L-J SIGMA PARAMETERS\n" << endi;

   if (opts.defined("bincoordinates"))
   {
     current = config->find("bincoordinates");

     iout << iINFO << "BINARY COORDINATES     " 
              << current->data << "\n";
   }


   if (firstTimestep)
   {
  iout << iINFO << "FIRST TIMESTEP         "
     << firstTimestep << "\n" << endi;
   }
}
/*    END OF FUNCTION initialize_config_data    */

/****************************************************************/
/*                */
/*    FUNCTION send_SimParameters      */
/*                */
/*  This function is used by the master process to broadcast*/
/*  the parameter data to all the other nodes.  It just builds  */
/*  a message with all the relevant data and broadcasts it to   */
/*  the other nodes.  The routine receive_SimParameters is used */
/*  by all the other nodes to receive this message.    */
/*                */
/****************************************************************/

void SimParameters::send_SimParameters(Communicate *com_obj)

{
  MOStream *msg = com_obj->newOutputStream(ALLBUTME, SIMPARAMSTAG, BUFSIZE);
  if ( msg == NULL )
  {
    NAMD_die("memory allocation failed in SimParameters::send_SimParameters");
  }

  msg->put(sizeof(SimParameters),(char*)this);
  if ( FFTWWisdomString ) {
    int fftwlen = strlen(FFTWWisdomString) + 1;
    msg->put(fftwlen);
    msg->put(fftwlen,FFTWWisdomString);
  }
  if ( tclBCScript ) {
    int tcllen = strlen(tclBCScript) + 1;
    msg->put(tcllen);
    msg->put(tcllen,tclBCScript);
  }

  msg->end();
}
/*    END OF FUNCITON send_SimParameters    */

/****************************************************************/
/*                */
/*      FUNCTION receive_SimParameters    */
/*                */
/*  This function is used by all the child nodes to   */
/*  receive the simulation parameters from the master node.  */
/*                */
/****************************************************************/

void SimParameters::receive_SimParameters(MIStream *msg)

{
  msg->get(sizeof(SimParameters),(char*)this);
  if ( FFTWWisdomString ) {
    int fftwlen;
    msg->get(fftwlen);
    FFTWWisdomString = new char[fftwlen];
    msg->get(fftwlen,FFTWWisdomString);
#ifdef NAMD_FFTW
    fftw_import_wisdom_from_string(FFTWWisdomString);
#endif
  }
  if ( tclBCScript ) {
    int tcllen;
    msg->get(tcllen);
    tclBCScript = new char[tcllen];
    msg->get(tcllen,tclBCScript);
  }

  delete msg;
}
/*      END OF FUNCTION receive_SimParameters  */

