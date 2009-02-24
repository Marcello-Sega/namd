/*
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Holds pointers to large molecule data structure, simulation parameters...
*/

#include "InfoStream.h"
#include "common.h"
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

#include "CompressPsf.h"
#include "PluginIOMgr.h"

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
  fflush(stdout);
  lattice = simParameters->lattice;

  StringList *molInfoFilename;
  // If it's AMBER force field, read the AMBER style files;
  // if it's GROMACS, read the GROMACS files;
  // Otherwise read the CHARMM style files

  if (simParameters->amberOn)
  { StringList *parmFilename = configList->find("parmfile");
    molInfoFilename = parmFilename;
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
  else if (simParameters->gromacsOn) {
    StringList *topFilename = configList->find("grotopfile");
    molInfoFilename = topFilename;
    StringList *coorFilename = configList->find("grocoorfile");
    // "gromacsFile" is a temporary data structure, which records all
    // the data from the topology file. After copying it into the
    // molecule and parameter and pdb, it will be deleted.
    GromacsTopFile *gromacsFile;
    gromacsFile = new GromacsTopFile(topFilename->data);
    parameters = new Parameters(gromacsFile,simParameters->minimizeCGOn);
    if (coorFilename != NULL)
      pdb = new PDB(coorFilename->data,gromacsFile);

    molecule = new Molecule(simParameters, parameters, gromacsFile);
    // XXX does Molecule(needAll,these,arguments)?

    delete gromacsFile; // XXX unimplemented

    // XXX add error handling when the file doesn't exist
    // XXX make sure the right things happen when the parameters are
    // not even specified.
    // NAMD_die("Failed to read AMBER parm file!");
    parameters->print_param_summary();
  }
  else if (simParameters->usePluginIO){
    if(simParameters->useCompressedPsf) {
        NAMD_die("Compressed psf file is not supported currently when using plugin IO");
    }
    PluginIOMgr *pIOMgr = new PluginIOMgr();
    
    iout << iWARN << "Plugin-based I/O is still in development and may still have bugs" << endl;

    molfile_plugin_t *pIOHandle = pIOMgr->getPlugin();
    if (pIOHandle == NULL) {
        NAMD_die("ERROR: Failed to match requested plugin type");
    }
    if ( pIOHandle->open_file_read == NULL )
       NAMD_die("ERROR: Selected plugin type cannot open files"); 
    if ( pIOHandle->read_structure == NULL )
       NAMD_die("ERROR: Selected plugin type cannot read structures"); 
    if ( pIOHandle->read_next_timestep == NULL )
       NAMD_die("ERROR: Selected plugin type cannot read coordinates"); 

    StringList *moleculeFilename = configList->find("structure");
    molInfoFilename = moleculeFilename;
    StringList *parameterFilename = configList->find("parameters");
    //****** BEGIN CHARMM/XPLOR type changes
    // For AMBER use different constructor based on parm_struct!!!  -JCP
    parameters = new Parameters(simParameters, parameterFilename);
    parameters->print_param_summary();

    int numAtoms = 0;
    //TODO: not sure about the name field in the handler
    void *plgFile = pIOHandle->open_file_read(moleculeFilename->data, 
                                              pIOHandle->name, &numAtoms);
    if(plgFile ==  NULL) {
        NAMD_die("ERROR: Opening structure file failed!");
    }

    double fileReadTime = CmiWallTimer();
    molecule = new Molecule(simParameters, parameters, pIOHandle, plgFile, numAtoms);
    iout << iINFO << "TIME FOR LOAD MOLECULE STRUCTURE INFORMATION: " << CmiWallTimer() - fileReadTime << "\n" << endi;

    /* If we are only generating compressed molecule information, the PDB object is not needed */
    if(!simParameters->genCompressedPsf) {        
        fileReadTime = CmiWallTimer();
        //get the occupancy data from the Molecule object and then free it
        //as it is stored in the Molecule object.
        pdb = new PDB(pIOHandle, plgFile, molecule->numAtoms, molecule->getOccupancyData(), molecule->getBFactorData());
        molecule->freeOccupancyData();
        molecule->freeBFactorData();
        iout << iINFO << "TIME FOR LOADING ATOMS' COORDINATES INFORMATION: " << CmiWallTimer() - fileReadTime << "\n" << endi;
    }

    pIOHandle->close_file_read(plgFile);
    delete pIOMgr;
  }
  else
  { StringList *moleculeFilename = configList->find("structure");
    molInfoFilename = moleculeFilename; 
    StringList *parameterFilename = configList->find("parameters");
    //****** BEGIN CHARMM/XPLOR type changes
    // For AMBER use different constructor based on parm_struct!!!  -JCP
    parameters = new Parameters(simParameters, parameterFilename);
    //****** END CHARMM/XPLOR type changes    

    parameters->print_param_summary();

    double fileReadTime = CmiWallTimer();
    molecule = new Molecule(simParameters, parameters, moleculeFilename->data);
    iout << iINFO << "TIME FOR READING PSF FILE: " << CmiWallTimer() - fileReadTime << "\n" << endi;    
  }
  fflush(stdout);

  if (simParameters->extraBondsOn) {        
      molecule->build_extra_bonds(parameters, configList->find("extraBondsFile"));         
  }
  if(simParameters->genCompressedPsf) {
    #ifdef MEM_OPT_VERSION
      iout << "In the memory optimized version, the compression of molecule information is not supported! "
              << "Please use the non-memory optimized version.\n" <<endi;
    #else
      double fileReadTime = CmiWallTimer();
      compress_molecule_info(molecule, molInfoFilename->data, parameters, simParameters, configList);
      iout << "Finished compressing molecule information, which takes " << CmiWallTimer()-fileReadTime <<"(s)\n"<<endi;
    #endif
      CkExit();
  }

  //If using plugin-based IO, the PDB object is already created!
  StringList *coordinateFilename = NULL;
  if(!simParameters->usePluginIO) {
      coordinateFilename = configList->find("coordinates");
    
      double fileReadTime = CmiWallTimer();
    #ifdef MEM_OPT_VERSION
      if (coordinateFilename != NULL)
        pdb = new PDB(coordinateFilename->data, molecule->numAtoms);
    #else
      if (coordinateFilename != NULL)
        pdb = new PDB(coordinateFilename->data);
      if (pdb->num_atoms() != molecule->numAtoms) {
        NAMD_die("Number of pdb and psf atoms are not the same!");
      }
    #endif
      iout << iINFO << "TIME FOR READING PDB FILE: " << CmiWallTimer() - fileReadTime << "\n" << endi;
      iout << iINFO << "\n" << endi;
  }

  StringList *binCoordinateFilename = configList->find("bincoordinates");
  if ( binCoordinateFilename ) {
    read_binary_coors(binCoordinateFilename->data, pdb);
  }
  

	//  If constraints are active, build the parameters necessary
	if (simParameters->constraintsOn)
	{
           StringList *consRefFile = configList->find("consref");
           StringList *consKFile = configList->find("conskfile");

          if (coordinateFilename != NULL) {
           if(strcasecmp(coordinateFilename->data, consRefFile->data)==0)
                consRefFile = NULL;
           if(strcasecmp(coordinateFilename->data, consKFile->data)==0)
                consKFile = NULL;
          }

           molecule->build_constraint_params(consRefFile, consKFile,
                                             configList->find("conskcol"),
                                             pdb,
                                             NULL);
	}
	//CkPrintf ("DEBUG--check if StirOn to build stir params..\n");
	if (simParameters->stirOn)
	{	
	//CkPrintf ("DEBUG--now to build stir params..\n");
	  
	   molecule->build_stirred_atoms(configList->find("stirFilename"),
				       configList->find("stirredAtomsCol"),
				       pdb,
				       NULL);
	}

	/*if (simParameters->extraBondsOn) {
       #ifdef MEM_OPT_VERSION
       iout << "Information about the extra bonds have been already integrated into the compressed psf file!\n" <<endi;
       #else
	   molecule->build_extra_bonds(parameters,
				configList->find("extraBondsFile"));
       #endif
	}*/
	
	if (simParameters->fixedAtomsOn)
	{
	   molecule->build_fixed_atoms(configList->find("fixedatomsfile"),
					configList->find("fixedatomscol"),
					pdb,
					NULL);
	}
	
	/* BEGIN gf */
	if (simParameters->mgridforceOn)
	{
	    molecule->build_gridforce_params(configList->find("gridforcefile"),
					     configList->find("gridforcecol"),
					     configList->find("gridforceqcol"),
					     configList->find("gridforcevfile"),
					     pdb,
					     NULL);
	}
	/* END gf */

	// If constant forces are active, read the forces necessary
	if (simParameters->consForceOn) {
    char *filename = NULL;
    if (configList->find("consforcefile"))
      filename = configList->find("consforcefile")->data;
    molecule->build_constant_forces(filename);
  }

        if (simParameters->excludeFromPressure) {
           molecule->build_exPressure_atoms(
             configList->find("excludeFromPressureFile"),
             configList->find("excludeFromPressureCol"),
	     pdb, NULL);
        }

	// If moving drag is active, build the parameters necessary
	if (simParameters->movDragOn) {
	  molecule->build_movdrag_params(configList->find("movDragFile"),
					 configList->find("movDragCol"),
					 configList->find("movDragVelFile"),
					 pdb,
					 NULL);
	}

	// If rotating drag is active, build the parameters necessary
	if (simParameters->rotDragOn) {
	  molecule->build_rotdrag_params(configList->find("rotDragFile"),
					 configList->find("rotDragCol"),
					 configList->find("rotDragAxisFile"),
					 configList->find("rotDragPivotFile"),
					 configList->find("rotDragVelFile"),
					 configList->find("rotDragVelCol"),
					 pdb,
					 NULL);
	}

	// If "constant" torque is active, build the parameters necessary
	if (simParameters->consTorqueOn) {
	  molecule->build_constorque_params(configList->find("consTorqueFile"),
				       configList->find("consTorqueCol"),
				       configList->find("consTorqueAxisFile"),
				       configList->find("consTorquePivotFile"),
				       configList->find("consTorqueValFile"),
				       configList->find("consTorqueValCol"),
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

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
     //identify the mutant atoms for fep simulation
        if (simParameters->fepOn) {
           molecule->build_fep_flags(configList->find("fepfile"),
                configList->find("fepcol"), pdb, NULL);
        }
        if (simParameters->thermInt) {
           molecule->build_fep_flags(configList->find("tifile"),
                configList->find("ticol"), pdb, NULL);
        }
//fepe
        if (simParameters->lesOn) {
	   if (simParameters->fepOn || simParameters->thermInt) NAMD_bug("FEP/TI and LES are incompatible!");
           molecule->build_fep_flags(configList->find("lesfile"),
                configList->find("lescol"), pdb, NULL);
        }
        if (simParameters->pairInteractionOn) {
           molecule->build_fep_flags(configList->find("pairInteractionFile"),
                configList->find("pairInteractionCol"), pdb, NULL);
        }      
        if (simParameters->pressureProfileAtomTypes > 1) {
          molecule->build_fep_flags(configList->find("pressureProfileAtomTypesFile"),
                configList->find("pressureProfileAtomTypesCol"), pdb, NULL);
        }

	iout << iINFO << "****************************\n";
	iout << iINFO << "STRUCTURE SUMMARY:\n";
	iout << iINFO << molecule->numAtoms << " ATOMS\n";
	iout << iINFO << molecule->numBonds << " BONDS\n";
	iout << iINFO << molecule->numAngles << " ANGLES\n";
	iout << iINFO << molecule->numDihedrals << " DIHEDRALS\n";
	iout << iINFO << molecule->numImpropers << " IMPROPERS\n";
	iout << iINFO << molecule->numCrossterms << " CROSSTERMS\n";
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

        if (simParameters->stirOn)
	  iout << iINFO << molecule->numStirredAtoms << " STIRRED ATOMS\n";

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
	
	/* BEGIN gf */
	if (simParameters->mgridforceOn)
	{
            int i;
	    iout << iINFO << molecule->numGridforceGrids 
	         << " GRIDS ACTIVE\n";
	}
	/* END gf */

//Modifications for alchemical fep
//SD & CC, CNRS - LCTN, Nancy
        if (simParameters->fepOn || simParameters->thermInt) {
           iout << iINFO << molecule->numFepInitial <<
               " ATOMS TO DISAPPEAR IN FINAL STATE\n";
           iout << iINFO << molecule->numFepFinal <<
               " ATOMS TO APPEAR IN FINAL STATE\n";
        }
//fepe

        if (simParameters->lesOn) {
           iout << iINFO << molecule->numFepInitial <<
               " LOCALLY ENHANCED ATOMS ENABLED\n";
        }
       
        if (simParameters->pairInteractionOn) {
           iout << iINFO << "PAIR INTERACTION GROUP 1 CONTAINS "
                <<  molecule->numFepInitial << " ATOMS\n";
           if (!simParameters->pairInteractionSelf) {
             iout << iINFO << "PAIR INTERACTION GROUP 2 CONTAINS "
                  <<  molecule->numFepFinal << " ATOMS\n";
           }
        }
           
#if 1
        if (molecule->numLonepairs != 0) {
          iout << iINFO << molecule->numLonepairs << " LONE PAIRS\n";
        }
        if (molecule->numDrudeAtoms != 0) {
          iout << iINFO << molecule->numDrudeAtoms << " DRUDE ATOMS\n";
        }
        iout << iINFO << molecule->num_deg_freedom(1)
             << " DEGREES OF FREEDOM\n";
        if (simParameters->drudeOn) {
          int g_bond = 3 * molecule->numDrudeAtoms;
          int g_com = molecule->num_deg_freedom(1) - g_bond;
          iout << iINFO << g_com << " DRUDE COM DEGREES OF FREEDOM\n";
          iout << iINFO << g_bond << " DRUDE BOND DEGREES OF FREEDOM\n";
        }
#endif
#if 0
	{
	  // Copied from Controller::printEnergies()
	  int numAtoms = molecule->numAtoms;
	  int numDegFreedom = 3 * numAtoms;
    int numLonepairs = molecule->numLonepairs;
	  int numFixedAtoms = molecule->numFixedAtoms;
	  if ( numFixedAtoms ) numDegFreedom -= 3 * numFixedAtoms;
	  if ( ! ( numFixedAtoms || molecule->numConstraints
		|| simParameters->comMove || simParameters->langevinOn ) ) {
	    numDegFreedom -= 3;
	  }
    if (numLonepairs) numDegFreedom -= 3 * numLonepairs;
	  int numRigidBonds = molecule->numRigidBonds;
	  int numFixedRigidBonds = molecule->numFixedRigidBonds;
    // numLonepairs is subtracted here because all lonepairs have a rigid bond
    // to oxygen, but all of the LP degrees of freedom are dealt with above
    numDegFreedom -= ( numRigidBonds - numFixedRigidBonds - numLonepairs);
	  iout << iINFO << numDegFreedom << " DEGREES OF FREEDOM\n";
	}
#endif

	iout << iINFO << molecule->numHydrogenGroups << " HYDROGEN GROUPS\n";
	if (simParameters->fixedAtomsOn)
	{
	   iout << iINFO << molecule->numFixedGroups <<
			" HYDROGEN GROUPS WITH ALL ATOMS FIXED\n";
	}

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
        fflush(stdout);

  DebugM(4, "::configFileInit() - printing Molecule Information\n");

  molecule->print_atoms(parameters);
  molecule->print_bonds(parameters);
  molecule->print_exclusions();
  fflush(stdout);

  DebugM(4, "::configFileInit() - done printing Molecule Information\n");
  DebugM(1, "::configFileInit() - done\n");

  return(0);
}

