/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/
#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "ConfigList.h"
#include "PDB.h"
#include "NamdState.h"
#include <unistd.h>

#define DEBUGM
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
	iout << "Config List exists\n" << endi;
    } else {
	ret++;
    }

    if (simParameters != NULL) {
	iout << "SimParameters exists\n" << endi;
    } else {
	ret++;
    }

    if (parameters != NULL) {
	iout << "Parameters exists\n" << endi;
    } else {
	ret++;
    }

    if (molecule != NULL) {
	iout << "Molecule exists\n" << endi;
    } else {
	ret++;
    }

    if (pdb != NULL) {
	iout << "PDB exists \n" << endi;
    } else {
	ret++;
    }

    return(ret);
}

int
NamdState::configFileInit(char *confFile)
{
  char *currentdir=confFile;
  char *tmp;
  for(tmp=confFile;*tmp;++tmp); // find final null
  for( ; tmp != confFile && *tmp != '/'; --tmp); // find last '/'
  if ( tmp != confFile )
  {
    *tmp = 0; confFile = tmp + 1; 
    if ( chdir(currentdir) ) NAMD_die("chdir() failed!");
    CPrintf("NamdState::configFileInit running in %s\n", currentdir);
  }
  else if ( *tmp == '/' ) // config file in / is odd, but it might happen
    if ( chdir("/") ) NAMD_die("chdir() failed!");
  currentdir = NULL;

  CPrintf("NamdState::configFileInit running %s\n",confFile);

  if ( NULL == confFile || NULL == (configList = new ConfigList(confFile)) ) {
    CPrintf("NamdState::configFileInit() Config File is NULL\n");
    return(1);
  }
  if (!configList->okay()) {
    CPrintf("NamdState::configFileInit() ConfigList is bad\n");
    return(1);
  }

  CPrintf("NamdState::configFileInit configList okay\n");

  StringList *moleculeFilename = configList->find("structure");

  CPrintf("NamdState::configFileInit got moleculeFilename \n");

  StringList *parameterFilename = configList->find("parameters");

  StringList *coordinateFilename = configList->find("coordinates");

  CPrintf("NamdState::configFileInit got Coordinates \n");

  iout << "files are : " << moleculeFilename->data << " and "
      << parameterFilename->data << " and " << coordinateFilename->data
      << "\n" << endi;

  simParameters =  new SimParameters(configList,currentdir);

  parameters = new Parameters(parameterFilename);
  parameters->print_param_summary();

  molecule = new Molecule(simParameters, parameters, moleculeFilename->data);
  DebugM(1, "Done Reading Molecule file " << moleculeFilename->data << endl );

  pdb = new PDB(coordinateFilename->data);
  DebugM(1, "Done Reading Coordinate file " << coordinateFilename->data << endl );

  if (pdb->num_atoms() != molecule->numAtoms)
  {
    CPrintf("Number of pdb and psf atoms are not the same!");
    return(1);
  }

  molecule->print_atoms(parameters);
  molecule->print_bonds(parameters);
  molecule->print_exclusions();


  DebugM(1, "::configFileInit() - done\n");

  return(0);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdState.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.779 $	$Date: 1997/02/06 15:53:16 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdState.C,v $
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

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/NamdState.C,v 1.779 1997/02/06 15:53:16 ari Exp $";
