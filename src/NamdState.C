/*
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Holds pointers to large molecule data structure, simulation parameters...
*/

#include "charm++.h"

#include "common.h"
#include "InfoStream.h"
#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "ConfigList.h"
#include "PDB.h"
#include "NamdState.h"
#include "Controller.h"
#include "ScriptTcl.h"
#ifndef WIN32
#include <unistd.h>
#endif
#include <sys/stat.h>
#include "parm.h"

//#define DEBUGM
#include "Debug.h"

NamdState::NamdState()
{
    configList = NULL;
    simParameters = NULL;
    parameters = NULL;
    molecule = NULL;
    pdb = NULL;
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

int NamdState::configListInit(ConfigList *cfgList) {
  configList = cfgList;
  if (!configList->okay()) {
    NAMD_die("Simulation config file is incomplete or contains errors.");
  }
  DebugM(1,"NamdState::configFileInit configList okay\n");

  char *currentdir = 0;
  simParameters =  new SimParameters(configList,currentdir);
  lattice = simParameters->lattice;
  
  // If it's AMBER force field, read the AMBER style files;
  // Otherwise read the CHARMM style files

  if (simParameters->amberOn)
  { StringList *parmFilename = configList->find("parmfile");
    StringList *coorFilename = configList->find("ambercoor");
    // "amber" is a temporary data structure, which records all
    // the data from the parm file. After copying them into
    // molecule, parameter and pdb structures, it will be deleted.
    Ambertoppar *amber;
    amber = new Ambertoppar;
    if (amber->readparm(parmFilename->data))
    { parameters = new Parameters(amber, simParameters->vdwscale14);
      molecule = new Molecule(simParameters, parameters, amber);
      if (coorFilename != NULL)
        pdb = new PDB(coorFilename->data,amber);
      delete amber;
    }
    else
      NAMD_die("Failed to read AMBER parm file!");
    parameters->print_param_summary();
  }
  else
  { StringList *moleculeFilename = configList->find("structure");
    StringList *parameterFilename = configList->find("parameters");
    //****** BEGIN CHARMM/XPLOR type changes
    // For AMBER use different constructor based on parm_struct!!!  -JCP
    parameters = new Parameters(simParameters, parameterFilename);
    //****** END CHARMM/XPLOR type changes
    parameters->print_param_summary();

    molecule = new Molecule(simParameters, parameters, moleculeFilename->data);
  }
  
  StringList *coordinateFilename = configList->find("coordinates");
  if (coordinateFilename != NULL)
    pdb = new PDB(coordinateFilename->data);
  if (pdb->num_atoms() != molecule->numAtoms) {
    iout << iERRORF 
      << "Number of pdb and psf atoms are not the same!" << "\n" << endi;
    return(1);
  }

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

	// If constant forces are active, read the forces necessary
	if (simParameters->consForceOn)
	  molecule->build_constant_forces(configList->find("consforcefile")->data);

        if (simParameters->excludeFromPressure) {
           molecule->build_exPressure_atoms(
             configList->find("excludeFromPressureFile"),
             configList->find("excludeFromPressureCol"),
	     pdb, NULL);
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

	if (simParameters->consForceOn)
	  iout << iINFO << molecule->numConsForce << " CONSTANT FORCES\n";

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

        {
          BigReal totalMass = 0;
          BigReal totalCharge = 0;
          int i;
          for ( i = 0; i < molecule->numAtoms; ++i ) {
            totalMass += molecule->atommass(i);
            totalCharge += molecule->atomcharge(i);
          }
          iout << iINFO << "TOTAL MASS = " << totalMass << " amu\n"; 
          iout << iINFO << "TOTAL CHARGE = " << totalCharge << " e\n"; 
        }

	iout << iINFO << "*****************************\n";
	iout << endi;

  DebugM(4, "::configFileInit() - printing Molecule Information\n");

  molecule->print_atoms(parameters);
  molecule->print_bonds(parameters);
  molecule->print_exclusions();

  DebugM(4, "::configFileInit() - done printing Molecule Information\n");
  DebugM(1, "::configFileInit() - done\n");

  return(0);
}

