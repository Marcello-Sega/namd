/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Holds pointers to large molecule data structure, simulation parameters...
*/

#include "charm++.h"

#include "Inform.h"
#include "common.h"
#include "InfoStream.h"
#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "ConfigList.h"
#include "PDB.h"
#include "NamdState.h"
#include "Controller.h"
#include "SMD.h"
#include <unistd.h>
#include <sys/stat.h>

//#define DEBUGM
#include "Debug.h"

NamdState::NamdState()
{
    configList = NULL;
    simParameters = NULL;
    parameters = NULL;
    molecule = NULL;
    pdb = NULL;
    smdData = NULL;
}

int
NamdState::status()
{
    int ret=0;
    if (configList != NULL) {
      DebugM(1, "Config List exists\n");
    } else ret++;

    if (simParameters != NULL) {
      DebugM(1, "SimParameters exists\n");
    }
    else ret++;

    if (parameters != NULL) {
      DebugM(1,"Parameters exists\n");
    }
    else ret++;

    if (molecule != NULL) {
      DebugM(1, "Molecule exists\n");
    }
    else ret++;

    if (pdb != NULL) {
      DebugM(1,"PDB exists \n");
    }
    else ret++;
    
    if (smdData != NULL) {
      DebugM(1, "SMDData exists\n");
    }
    else ret++;

    return(ret);
}

void NamdState:: useController(Controller *controllerPtr)
{
  controller=controllerPtr;
}

void NamdState:: runController(void)
{
  controller->run();
}

extern void read_binary_coors(char *fname, PDB *pdbobj);

int
NamdState::configFileInit(char *confFile)
{
  struct stat statBuf;

  char *currentdir=confFile;
  char *tmp;
  for(tmp=confFile;*tmp;++tmp); // find final null
  for( ; tmp != confFile && *tmp != '/'; --tmp); // find last '/'
  if ( tmp != confFile )
  {
    *tmp = 0; confFile = tmp + 1; 
    if ( chdir(currentdir) ) NAMD_die("chdir() failed!");
    iout << iINFO << "Changed directory to " << currentdir << "\n" << endi;
  }
  else if ( *tmp == '/' ) // config file in / is odd, but it might happen
    if ( chdir("/") ) NAMD_die("chdir() failed!");
  currentdir = NULL;

  iout << iINFO << "Configuration file is " << confFile << "\n" << endi;
  if (stat(confFile, &statBuf)) {
    NAMD_die("Simulation config file is not accessible.");
  }

  if ( NULL == confFile || NULL == (configList = new ConfigList(confFile)) ) {
    NAMD_die("Simulation config file is empty.");
  }
  if (!configList->okay()) {
    NAMD_die("Simulation config file is incomplete or contains errors.");
  }
  DebugM(1,"NamdState::configFileInit configList okay\n");

  StringList *moleculeFilename = configList->find("structure");
  StringList *parameterFilename = configList->find("parameters");
  StringList *coordinateFilename = configList->find("coordinates");
  simParameters =  new SimParameters(configList,currentdir);
  lattice = simParameters->lattice;
  //****** BEGIN CHARMM/XPLOR type changes
  parameters = new Parameters(simParameters, parameterFilename);
  //****** END CHARMM/XPLOR type changes
  parameters->print_param_summary();

  molecule = new Molecule(simParameters, parameters, moleculeFilename->data);
  pdb = new PDB(coordinateFilename->data);
  if (pdb->num_atoms() != molecule->numAtoms) {
    iout << iERRORF 
      << "Number of pdb and psf atoms are not the same!" << "\n" << endi;
    return(1);
  }
  smdData = new SMDData(simParameters);
  smdData->init(pdb);

  StringList *binCoordinateFilename = configList->find("bincoordinates");
  if ( binCoordinateFilename ) {
    read_binary_coors(binCoordinateFilename->data, pdb);
  }

	//  If constraints are active, build the parameters necessary
	if (simParameters->constraintsOn)
	{
	   molecule->build_constraint_params(configList->find("consref"),
					      configList->find("conskfile"),
					      configList->find("conskcol"),
					      pdb,
					      NULL);
	}

	if (simParameters->fixedAtomsOn)
	{
	   molecule->build_fixed_atoms(configList->find("fixedatomsfile"),
					configList->find("fixedatomscol"),
					pdb,
					NULL);
	}

	//  If langevin dynamics or temperature coupling are active, build 
	//  the parameters necessary
	if (simParameters->langevinOn)
	{
	  if (simParameters->langevinDamping == 0.0) {
	    molecule->build_langevin_params(configList->find("langevinfile"),
					    configList->find("langevincol"),
					    pdb,
					    NULL);
	  } else {
	    molecule->build_langevin_params(simParameters->langevinDamping,
					    simParameters->langevinHydrogen);
	  }
	}
	else if (simParameters->tCoupleOn)
	{
	   //  Temperature coupling uses the same parameters, but with different
	   //  names . . .
	   molecule->build_langevin_params(configList->find("tcouplefile"),
					    configList->find("tcouplecol"),
					    pdb,
					    NULL);
	}
	
	iout << iINFO << "****************************\n";
	iout << iINFO << "STRUCTURE SUMMARY:\n";
	iout << iINFO << molecule->numAtoms << " ATOMS\n";
	iout << iINFO << molecule->numBonds << " BONDS\n";
	iout << iINFO << molecule->numAngles << " ANGLES\n";
	iout << iINFO << molecule->numDihedrals << " DIHEDRALS\n";
	iout << iINFO << molecule->numImpropers << " IMPROPERS\n";
	iout << iINFO << molecule->numExclusions << " EXCLUSIONS\n";

        //****** BEGIN CHARMM/XPLOR type changes
	if ((molecule->numMultipleDihedrals) && (simParameters->paraTypeXplorOn))
	{
		iout << iINFO << molecule->numMultipleDihedrals 
	     << " DIHEDRALS WITH MULTIPLE PERIODICITY (BASED ON PSF FILE)\n";
	}
	if ((molecule->numMultipleDihedrals) && (simParameters->paraTypeCharmmOn))
	{
		iout << iINFO << molecule->numMultipleDihedrals 
	 << " DIHEDRALS WITH MULTIPLE PERIODICITY IGNORED (BASED ON PSF FILE) \n";
		iout << iINFO  
	 << " CHARMM MULTIPLICITIES BASED ON PARAMETER FILE INFO! \n";
	}
        //****** END CHARMM/XPLOR type changes

	if (molecule->numMultipleImpropers)
	{
		iout << iINFO << molecule->numMultipleImpropers 
			 << " IMPROPERS WITH MULTIPLE PERIODICITY\n";
	}
	
	if (simParameters->constraintsOn)
	{
	   iout << iINFO << molecule->numConstraints << " CONSTRAINTS\n";
	}

	if (simParameters->fixedAtomsOn)
	{
	   iout << iINFO << molecule->numFixedAtoms << " FIXED ATOMS\n";
	}

	if (simParameters->rigidBonds)
	{
	   iout << iINFO << molecule->numRigidBonds << " RIGID BONDS\n";
	}

	if (simParameters->fixedAtomsOn && simParameters->rigidBonds)
	{
	   iout << iINFO << molecule->numFixedRigidBonds <<
			" RIGID BONDS BETWEEN FIXED ATOMS\n";
	}

	{
	  // Copied from Controller::printEnergies()
	  int numAtoms = molecule->numAtoms;
	  int numDegFreedom = 3 * numAtoms;
	  int numFixedAtoms = molecule->numFixedAtoms;
	  if ( numFixedAtoms ) numDegFreedom -= 3 * numFixedAtoms;
	  if ( ! ( numFixedAtoms || molecule->numConstraints
		|| simParameters->comMove || simParameters->langevinOn ) ) {
	    numDegFreedom -= 3;
	  }
	  int numRigidBonds = molecule->numRigidBonds;
	  int numFixedRigidBonds = molecule->numFixedRigidBonds;
	  numDegFreedom -= ( numRigidBonds - numFixedRigidBonds );
	  iout << iINFO << numDegFreedom << " DEGREES OF FREEDOM\n";
	}

	iout << iINFO << molecule->numHydrogenGroups << " HYDROGEN GROUPS\n";

	iout << iINFO << "*****************************\n";
	iout << endi;

  DebugM(4, "::configFileInit() - printing Molecule Information\n");

  namdDebug.on(1);

  molecule->print_atoms(parameters);
  molecule->print_bonds(parameters);
  molecule->print_exclusions();

  namdDebug.on(0);

  DebugM(4, "::configFileInit() - done printing Molecule Information\n");
  DebugM(1, "::configFileInit() - done\n");

  return(0);
}

