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
	iout << iINFO << "Config List exists\n" << endi;
    } else {
	ret++;
    }

    if (simParameters != NULL) {
	iout << iINFO << "SimParameters exists\n" << endi;
    } else {
	ret++;
    }

    if (parameters != NULL) {
	iout << iINFO << "Parameters exists\n" << endi;
    } else {
	ret++;
    }

    if (molecule != NULL) {
	iout << iINFO << "Molecule exists\n" << endi;
    } else {
	ret++;
    }

    if (pdb != NULL) {
	iout << iINFO << "PDB exists\n" << endi;
    } else {
	ret++;
    }

    return(ret);
}

int
NamdState::configFileInit(char *confFile)
{
  char *currentdir=NULL;

  DebugM(1,"NamdState::configFileInit running " << confFile << "\n");

  if ( NULL == confFile || NULL == (configList = new ConfigList(confFile)) ) {
    DebugM(1,"NamdState::configFileInit() Config File is NULL\n");
    return(1);
  }
  if (!configList->okay()) {
    DebugM(1,"NamdState::configFileInit() ConfigList is bad\n");
    return(1);
  }

  StringList *moleculeFilename = configList->find("structure");
  StringList *parameterFilename = configList->find("parameters");
  StringList *coordinateFilename = configList->find("coordinates");
  iout << iINFO << "files are : " << moleculeFilename->data << " and "
      << parameterFilename->data << " and " << coordinateFilename->data
      << "\n" << endi;

  simParameters =  new SimParameters(configList,currentdir);

  parameters = new Parameters(parameterFilename->data);
  parameters->print_param_summary();

  molecule = new Molecule(simParameters, parameters, moleculeFilename->data);
  iout << iINFO << "Done Reading Molecule file\n" << endi;

  pdb = new PDB(coordinateFilename->data);
  iout << iINFO << "Done Reading Coordinate file\n" << endi;

  if (pdb->num_atoms() != molecule->numAtoms)
  {
    iout << iWARN << "Number of pdb and psf atoms are not the same!\n" << endi;
    return(1);
  }

  return(0);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdState.C,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/11/14 21:00:16 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdState.C,v $
 * Revision 1.3  1996/11/14 21:00:16  nealk
 * Now uses iout stream and DebugM for output.
 *
 * Revision 1.2  1996/08/16 04:55:30  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 04:39:46  ari
 * Initial revision
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/NamdState.C,v 1.3 1996/11/14 21:00:16 nealk Exp $";
