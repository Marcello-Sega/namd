/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/Namd.C,v 1.1 1996/08/16 01:54:59 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "Namd.h"
#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "ConfigList.h"
//#include "Node.h"

// Namd(void ) is the constructor for the startup node.  It needs to
// read in file data,
Namd::Namd(void)
{
}


// ~Namd(void) just needs to tell all the slave nodes to die.
Namd::~Namd(void)
{
}


// startup(char *) 
void Namd::startup(char *confFile)
{
  Molecule *molecule;
  Parameters *parameter;
  SimParameters *simParams;
  ConfigList *configList;

    char *currentdir=NULL;
  // this function is pretty much similar to master_startupp of namd.1.4
  // for this skeleton version it is empty
  CPrintf("Namd::startup running %s\n",confFile);

  if ( confFile != NULL )
    configList = new ConfigList(confFile);
  else
    NAMD_die("Namd::startup() Config File ptr is NULL");

  StringList *moleculeFilename = configList->find("structure");
  StringList *parameterFilename = configList->find("parameters");
  StringList *coordinateFilename = configList->find("coordinates");

  cout << "files are : " << moleculeFilename->data << " and "
      << parameterFilename->data << " and " << coordinateFilename->data
      << endl;

  simParams =  new SimParameters();
  molecule = new Molecule(simParams);
  parameter = new Parameters();

  simParams->initialize_config_data(configList,currentdir);
  parameter->read_parameter_file(parameterFilename->data);
  parameter->done_reading_files();
  parameter->print_param_summary();

  cout << "Reading Molecule file" << endl;
  molecule->read_psf_file(moleculeFilename->data, parameter);
  cout << "Done Reading Molecule file" << endl;
}


// run(void) runs the specified simulation to completion
void Namd::run(void)
{
  // node->run();
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Namd.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/16 01:54:59 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Namd.C,v $
 * Revision 1.1  1996/08/16 01:54:59  ari
 * Initial revision
 *
 *
 *
 ***************************************************************************/
