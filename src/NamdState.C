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
	cout << "Config List exists" << endl;
    } else {
	ret++;
    }

    if (simParameters != NULL) {
	cout << "SimParameters exists" << endl;
    } else {
	ret++;
    }

    if (parameters != NULL) {
	cout << "Parameters exists" << endl;
    } else {
	ret++;
    }

    if (molecule != NULL) {
	cout << "Molecule exists" << endl;
    } else {
	ret++;
    }

    if (pdb != NULL) {
	cout << "PDB exists" << endl;
    } else {
	ret++;
    }

    return(ret);
}

int
NamdState::configFileInit(char *confFile)
{
  char *currentdir=NULL;

  CPrintf("NamdState::configFileInit running %s\n",confFile);

  if ( NULL == confFile || NULL == (configList = new ConfigList(confFile)) ) {
    CPrintf("NamdState::configFileInit() Config File is NULL\n");
    return(1);
  }
  if (!configList->okay()) {
    CPrintf("NamdState::configFileInit() ConfigList is bad\n");
    return(1);
  }

  StringList *moleculeFilename = configList->find("structure");
  StringList *parameterFilename = configList->find("parameters");
  StringList *coordinateFilename = configList->find("coordinates");
  cout << "files are : " << moleculeFilename->data << " and "
      << parameterFilename->data << " and " << coordinateFilename->data
      << endl;

  simParameters =  new SimParameters(configList,currentdir);

  parameters = new Parameters(parameterFilename->data);
  parameters->print_param_summary();

  molecule = new Molecule(simParameters, parameters, moleculeFilename->data);
  cout << "Done Reading Molecule file" << endl;

  pdb = new PDB(coordinateFilename->data);
  cout << "Done Reading Coordinate file" << endl;

  if (pdb->num_atoms() != molecule->numAtoms)
  {
    CPrintf("Number of pdb and psf atoms are not the same!");
    return(1);
  }

  return(0);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: NamdState.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/08/16 04:55:30 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: NamdState.C,v $
 * Revision 1.2  1996/08/16 04:55:30  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 04:39:46  ari
 * Initial revision
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/NamdState.C,v 1.2 1996/08/16 04:55:30 ari Exp $";
