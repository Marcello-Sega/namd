/***************************************************************************/
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Holds pointers to large molecule data structure, simulation
 *		Parameters...
 ***************************************************************************/

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
    iout << iERRORF << "Configuration File is not accessible." << "\n" << endi;
    return(1);
  }

  if ( NULL == confFile || NULL == (configList = new ConfigList(confFile)) ) {
    iout << iERRORF << "Configuration List is NULL." << "\n" << endi;
    return(1);
  }
  if (!configList->okay()) {
    iout << iERRORF << " ConfigList is bad." << "\n" << endi;
    return(1);
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

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdState.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1019 $	$Date: 1999/07/08 21:26:54 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdState.C,v $
 * Revision 1.1019  1999/07/08 21:26:54  jim
 * Eliminated compiler warnings.
 *
 * Revision 1.1018  1999/03/09 01:44:15  jim
 * Added langevinDamping and langevinHydrogen parameters.
 *
 * Revision 1.1017  1999/02/02 08:02:31  ferenc
 * Added support for CHARMM parameter format in parameter files.
 *
 * Revision 1.1016  1998/10/24 19:57:43  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1015  1998/09/13 21:06:09  jim
 * Cleaned up output, defaults, etc.
 *
 * Revision 1.1014  1998/03/03 23:05:18  brunner
 * Changed include files for new simplified Charm++ include file structure.
 *
 * Revision 1.1013  1998/02/18 05:38:31  jim
 * RigidBonds mainly finished.  Now temperature is correct and a form
 * of Langevin dynamics works with constraints.
 *
 * Revision 1.1012  1998/01/13 23:51:45  sergei
 * *** empty log message ***
 *
 * Revision 1.1010  1997/09/21 21:58:31  jim
 * Added printing of hydrogen group count.
 *
 * Revision 1.1009  1997/09/19 09:39:05  jim
 * Small tweaks for fixed atoms.
 *
 * Revision 1.1008  1997/09/19 08:55:34  jim
 * Added rudimentary but relatively efficient fixed atoms.  New options
 * are fixedatoms, fixedatomsfile, and fixedatomscol (nonzero means fixed).
 * Energies will be affected, although this can be fixed with a little work.
 *
 * Revision 1.1007  1997/03/24 01:43:59  jim
 * Added Langevin dynamics.
 *
 * Revision 1.1006  1997/03/21 23:05:38  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1005  1997/03/18 18:09:06  jim
 * Revamped collection system to ensure ordering and eliminate
 * unnecessary collections.  Also reduced make dependencies.
 *
 * Revision 1.1004  1997/03/10 17:40:13  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1003  1997/03/04 22:37:12  ari
 * Clean up of code.  Debug statements removal, dead code removal.
 * Minor fixes, output fixes.
 * Commented some code from the top->down.  Mainly reworked Namd, Node, main.
 *
 * Revision 1.1002  1997/02/26 23:18:44  jim
 * Now should read binary coordinate files - untested.
 *
 * Revision 1.1001  1997/02/11 18:51:48  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1000  1997/02/06 15:58:48  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.779  1997/02/06 15:53:16  ari
 * Updating Revision Line, getting rid of branches
 *
 * Revision 1.778.2.1  1997/02/06 02:35:26  jim
 * Implemented periodic boundary conditions - may not work with
 * atom migration yet, but doesn't seem to alter calculation,
 * appears to work correctly when turned on.
 * NamdState chdir's to same directory as config file in argument.
 *
 * Revision 1.778  1997/01/28 00:30:56  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/24 02:29:52  jim
 * Fixed bug where only first parameter file was read!
 * Added files for hydrogen bond parameter reading.
 *
 * Revision 1.777  1997/01/17 19:36:31  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.6  1996/12/10 00:13:12  ari
 * *** empty log message ***
 *
 * Revision 1.5  1996/12/06 17:28:15  jim
 * put build_lists_by_atom back where it belongs
 *
 * Revision 1.4  1996/11/21 23:34:24  ari
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 04:55:30  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 04:39:46  ari
 * Initial revision
 *
 ***************************************************************************/

