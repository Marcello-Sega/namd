/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: main.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/08/06 20:38:38 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: main.C,v $
 * Revision 1.2  1996/08/06 20:38:38  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/02 19:20:13  gursoy
 * Initial revision
 *
 *
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/main.C,v 1.2 1996/08/06 20:38:38 ari Exp $";

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.top.h"
#include "main.h"

#include "Molecule.h"
#include "Parameters.h"
#include "SimParameters.h"
#include "Communicate.h"
#include "Inform.h"

Communicate *comm;
Inform namdErr("ERROR");
Inform namdWarn("Warning");
Inform namdInfo("Info");
Inform namdDebug("** DEBUG **");

class main : public chare_object
{
public:
  main(int argc, char **argv)
  {

    Molecule *molecule;
    Parameters *parameter;
    SimParameters *simParams;
    ConfigList *configList;

    char *currentdir=NULL;

    comm = new Communicate();

    configList = new ConfigList(argv[1]);

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

    CharmExit();
  }
};

#include "main.bot.h"
