/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   This object outputs the data collected on the master node
*/

#include <string.h>
#include <stdlib.h>

#include "InfoStream.h"
#include "IMDOutput.h"
#include "Output.h"
#include "dcdlib.h"
#include "strlib.h"
#include "Molecule.h"
#include "Node.h"
#include "Parameters.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Vector.h"
#include "structures.h"
#include "MStream.h"
#include "Communicate.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "ScriptTcl.h"
#include "Lattice.h"

// These make the NAMD 1 names work in NAMD 2
#define namdMyNode Node::Object()
#define simParams simParameters
#define pdbData pdb

/************************************************************************/
/*                  */
/*      FUNCTION Output          */
/*                  */
/*  This is the constructor for the Ouput class.  It just sets   */
/*  up some values for the VMD connection.        */
/*                  */
/************************************************************************/

Output::Output() { }

/*      END OF FUNCTION Output        */

/************************************************************************/
/*                  */
/*      FUNCTION ~Output        */
/*                  */
/************************************************************************/

Output::~Output() { }

/*      END OF FUNCTION ~Output        */

/************************************************************************/
/*                  */
/*      FUNCTION coordinate        */
/*                  */
/*   INPUTS:                */
/*  timestep - Timestep coordinates were accumulated for    */
/*  n - This is the number of coordinates accumulated.    */
/*  vel - Array of Vectors containing the velocities    */
/*                  */
/*  This function receives the coordinates accumulated for a given  */
/*   timestep from the Collect object and calls the appropriate output  */
/*   functions.  ALL routines used to output coordinates information    */
/*   should be called from here.          */
/*                  */
/************************************************************************/

int Output::coordinateNeeded(int timestep)
{
  SimParameters *simParams = Node::Object()->simParameters;

  int positionsNeeded = 0;

  if ( timestep >= 0 ) {

    //  Output a DCD trajectory 
    if ( simParams->dcdFrequency &&
       ((timestep % simParams->dcdFrequency) == 0) )
    { positionsNeeded |= 1; }

    //  Output a restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
    { positionsNeeded |= 2; }

    //  Iteractive MD
    if ( simParams->IMDon &&
       ( ((timestep % simParams->IMDfreq) == 0) ||
         (timestep == simParams->firstTimestep) ) )
      { positionsNeeded |= 1; }

  }

  //  Output final coordinates
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN ||
					timestep == EVAL_MEASURE)
  {
    positionsNeeded |= 2;
  }

  return positionsNeeded;
}

template <class xVector, class xDone>
void wrap_coor_int(xVector *coor, Lattice &lattice, xDone *done) {
  SimParameters *simParams = Node::Object()->simParameters;
  if ( *done ) return;
  *done = 1;
  if ( ! ( simParams->wrapAll || simParams->wrapWater ) ) return;
  const int wrapNearest = simParams->wrapNearest;
  const int wrapAll = simParams->wrapAll;
  Molecule *molecule = Node::Object()->molecule;
  int n = molecule->numAtoms;
  int i;
  Position *con = new Position[n];
  for ( i = 0; i < n; ++i ) {
    con[i] = 0;
    int ci = molecule->get_cluster(i);
    con[ci] += coor[i];
  }
  for ( i = 0; i < n; ++i ) {
    if ( ! wrapAll && ! molecule->is_water(i) ) continue;
    int ci = molecule->get_cluster(i);
    if ( ci == i ) {
      Vector coni = con[i] / molecule->get_clusterSize(i);
      Vector trans = ( wrapNearest ?
	lattice.wrap_nearest_delta(coni) : lattice.wrap_delta(coni) );
      con[i] = trans;
    }
    coor[i] = coor[i] + con[ci];
  }
  delete [] con;
}

void wrap_coor(Vector *coor, Lattice &lattice, double *done) {
  wrap_coor_int(coor,lattice,done);
};

void wrap_coor(FloatVector *coor, Lattice &lattice, float *done) {
  wrap_coor_int(coor,lattice,done);
};

void Output::coordinate(int timestep, int n, Vector *coor, FloatVector *fcoor,
							Lattice &lattice)
{
  SimParameters *simParams = Node::Object()->simParameters;
  double coor_wrapped = 0;
  float fcoor_wrapped = 0;

  if ( timestep >= 0 ) {

    //  Output a DCD trajectory 
    if ( simParams->dcdFrequency &&
       ((timestep % simParams->dcdFrequency) == 0) )
    {
      wrap_coor(fcoor,lattice,&fcoor_wrapped);
      output_dcdfile(timestep, n, fcoor, 
          simParams->dcdUnitCell ? &lattice : NULL);
    }

    //  Output a restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
    {
      iout << "WRITING COORDINATES TO RESTART FILE AT STEP "
				<< timestep << "\n" << endi;
      wrap_coor(coor,lattice,&coor_wrapped);
      output_restart_coordinates(coor, n, timestep);
      iout << "FINISHED WRITING RESTART COORDINATES\n" <<endi;
    }

    //  Interactive MD
    if ( simParams->IMDon &&
       ( ((timestep % simParams->IMDfreq) == 0) ||
         (timestep == simParams->firstTimestep) ) )
    {
      IMDOutput *imd = Node::Object()->imd;
      wrap_coor(fcoor,lattice,&fcoor_wrapped);
      if (imd != NULL) imd->gather_coordinates(timestep, n, fcoor);
    }

  }

  if (timestep == EVAL_MEASURE)
  {
#ifdef NAMD_TCL
    wrap_coor(coor,lattice,&coor_wrapped);
    Node::Object()->getScript()->measure(coor);
#endif
  }

  //  Output final coordinates
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
  {
    iout << "WRITING COORDINATES TO OUTPUT FILE AT STEP "
				<< simParams->N << "\n" << endi;
    wrap_coor(coor,lattice,&coor_wrapped);
    output_final_coordinates(coor, n, simParams->N);
  }

  //  Close trajectory files
  if (timestep == END_OF_RUN)
  {
    if (simParams->dcdFrequency) output_dcdfile(END_OF_RUN,0,0, 
        simParams->dcdUnitCell ? &lattice : NULL);
  }

}
/*    END OF FUNCTION coordinate        */

/************************************************************************/
/*                  */
/*      FUNCTION velocity        */
/*                  */
/*   INPUTS:                */
/*  timestep - Timestep velocities were accumulated for    */
/*  n - This is the number of velocities accumulated.    */
/*  vel - Array of Vectors containing the velocities    */
/*                  */
/*  This function receives the velocities accumulated for a given   */
/*   timestep from the Collect object and calls the appropriate output  */
/*   functions.  ALL routines used to output velocity information should*/
/*   be called from here.            */
/*                  */
/************************************************************************/

int Output::velocityNeeded(int timestep)
{
  SimParameters *simParams = Node::Object()->simParameters;

  int velocitiesNeeded = 0;

  if ( timestep >= 0 ) {

    //  Output a velocity DCD trajectory
    if ( simParams->velDcdFrequency &&
       ((timestep % simParams->velDcdFrequency) == 0) )
      { velocitiesNeeded |= 1; }

    //  Output a restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
      { velocitiesNeeded |= 2; }

  }

  //  Output final velocities
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
  {
    velocitiesNeeded |= 2;
  }

  return velocitiesNeeded;
}

void Output::velocity(int timestep, int n, Vector *vel)
{
  SimParameters *simParams = Node::Object()->simParameters;

  if ( timestep >= 0 ) {

    //  Output velocity DCD trajectory
    if ( simParams->velDcdFrequency &&
       ((timestep % simParams->velDcdFrequency) == 0) )
    {
      output_veldcdfile(timestep, n, vel);
    }

  //  Output restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
    {
      iout << "WRITING VELOCITIES TO RESTART FILE AT STEP "
				<< timestep << "\n" << endi;
      output_restart_velocities(timestep, n, vel);
      iout << "FINISHED WRITING RESTART VELOCITIES\n" <<endi;
    }

  }

  //  Output final velocities
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
  {
    iout << "WRITING VELOCITIES TO OUTPUT FILE AT STEP "
				<< simParams->N << "\n" << endi;
    output_final_velocities(simParams->N, n, vel);
  }

  //  Close trajectory files
  if (timestep == END_OF_RUN)
  {
    if (simParams->velDcdFrequency) output_veldcdfile(END_OF_RUN,0,0);
  }

}
/*      END OF FUNCTION velocity      */

/************************************************************************/
/*                  */
/*      FUNCTION output_restart_coordinates    */
/*                  */
/*   INPUTS:                */
/*  coor - Array of vectors containing current positions    */
/*  n - Number of coordinates to output        */
/*  timestep - Timestep for which the coordinates are being written */
/*                  */
/*  This function writes out the current positions of all the atoms */
/*   in PDB format to the restart file specified by the user in the   */
/*   configuration file.            */
/*                  */
/************************************************************************/

void Output::output_restart_coordinates(Vector *coor, int n, int timestep)

{
  char comment[128];    //  Comment for header of PDB file
  char timestepstr[20];

  int baselen = strlen(namdMyNode->simParams->restartFilename);
  char *restart_name = new char[baselen+26];

  strcpy(restart_name, namdMyNode->simParams->restartFilename);
  if ( namdMyNode->simParams->restartSave ) {
    sprintf(timestepstr,".%d",timestep);
    strcat(restart_name, timestepstr);
  }
  strcat(restart_name, ".coor");

  NAMD_backup_file(restart_name,".old");

  //  Check to see if we should generate a binary or PDB file
  if (!namdMyNode->simParams->binaryRestart)
  {
    //  Generate a PDB restart file
    sprintf(comment, "RESTART COORDINATES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    namdMyNode->pdbData->set_all_positions(coor);
    namdMyNode->pdbData->write(restart_name, comment);
  }
  else
  {
    //  Generate a binary restart file
    write_binary_file(restart_name, n, coor);
  }

  delete [] restart_name;
}
/*      END OF FUNCTION output_restart_coordinates  */

/************************************************************************/
/*                  */
/*      FUNCTION output_restart_velocities    */
/*                  */
/*   INPUTS:                */
/*  vel - Array of vectors containing current velocities    */
/*  timestep - Timestep for which the velocities are being written  */
/*                  */
/*  This function writes out the current velocites of all the atoms */
/*   in PDB format to the restart file specified by the user in the   */
/*   configuration file.            */
/*                  */
/************************************************************************/

void Output::output_restart_velocities(int timestep, int n, Vector *vel)

{
  char comment[128];    //  comment for the header of PDB file
  char timestepstr[20];

  int baselen = strlen(namdMyNode->simParams->restartFilename);
  char *restart_name = new char[baselen+26];

  strcpy(restart_name, namdMyNode->simParams->restartFilename);
  if ( namdMyNode->simParams->restartSave ) {
    sprintf(timestepstr,".%d",timestep);
    strcat(restart_name, timestepstr);
  }
  strcat(restart_name, ".vel");

  NAMD_backup_file(restart_name,".old");

  //  Check to see if we should write out a PDB or a binary file
  if (!namdMyNode->simParams->binaryRestart)
  {
    //  Write the coordinates to a PDB file.  Multiple them by 20
    //  first to make the numbers bigger
    sprintf(comment, "RESTART VELOCITIES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    scale_vels(vel, n, PDBVELFACTOR);
    namdMyNode->pdbData->set_all_positions(vel);
    namdMyNode->pdbData->write(restart_name, comment);
    scale_vels(vel, n, PDBVELINVFACTOR);
  }
  else
  {
    //  Write the velocities to a binary file
    write_binary_file(restart_name, n, vel);
  }

  delete [] restart_name;
}
/*      END OF FUNCTION output_restart_velocities  */

/************************************************************************/
/*                  */
/*      FUNCTION output_dcdfile        */
/*                  */
/*   INPUTS:                */
/*  timestep - Current timestep          */
/*  n - Number of atoms in simulation        */
/*  coor - Coordinate vectors for all atoms        */
/*  lattice - periodic cell data; NULL if not to be written */
/*                  */
/*  This function maintains the interface between the Output object */
/*   and the dcd writing routines contained in dcdlib.      */
/*                  */
/************************************************************************/

#define RAD2DEG 180.0/3.14159265359

void Output::output_dcdfile(int timestep, int n, FloatVector *coor,
    const Lattice *lattice)

{
  static Bool first=TRUE;  //  Flag indicating first call
  static int fileid;  //  File id for the dcd file
  static float *x, *y, *z; // Arrays to hold x, y, and z arrays
  int i;      //  Loop counter
  int ret_code;    //  Return code from DCD calls
  SimParameters *simParams = namdMyNode->simParams;

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( timestep == END_OF_RUN ) {
    if ( ! first ) {
      iout << "CLOSING COORDINATE DCD FILE\n" << endi;
      close_dcd_write(fileid);
    } else {
      iout << "COORDINATE DCD FILE WAS NOT CREATED\n" << endi;
    }
    return;
  }

  if (first)
  {
    //  Allocate x, y, and z arrays since the DCD file routines
    //  need them passed as three independant arrays to be
    //  efficient
    x = new float[n];
    y = new float[n];
    z = new float[n];

    if ( (x==NULL) || (y==NULL) || (z==NULL) )
    {
      NAMD_err("memory allocation failed in Output::output_dcdfile");
    }

    //  Open the DCD file
    iout << "OPENING COORDINATE DCD FILE\n" << endi;

    fileid=open_dcd_write(simParams->dcdFilename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "DCD file %s already exists!!",
        simParams->dcdFilename);

      NAMD_err(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open DCD file %s",
        simParams->dcdFilename);

      NAMD_err(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->dcdFrequency;
    NSTEP = NSAVC * (simParams->N/NSAVC);
    NPRIV = simParams->firstTimestep+NSAVC;
    NPRIV = NSAVC * (NPRIV/NSAVC);
    NFILE = (NSTEP-NPRIV)/NSAVC + 1;

    //  Write out the header
    ret_code = write_dcdheader(fileid, 
        simParams->dcdFilename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR, lattice != NULL);


    if (ret_code<0)
    {
      NAMD_err("Writing of DCD header failed!!");
    }

    first = FALSE;
  }

  //  Copy the coordinates for output
  for (i=0; i<n; i++)
  {
    x[i] = coor[i].x;
    y[i] = coor[i].y;
    z[i] = coor[i].z;
  }

  //  Write out the values for this timestep
  iout << "WRITING COORDINATES TO DCD FILE AT STEP "
	<< timestep << "\n" << endi;
  if (lattice) {
    double unitcell[6];
    if (lattice->a_p() && lattice->b_p() && lattice->c_p()) {
      const Vector &a=lattice->a();
      const Vector &b=lattice->b();
      const Vector &c=lattice->c();
      unitcell[0] = a.length();
      unitcell[2] = b.length();
      unitcell[5] = c.length();
      double cosAB = (a*b)/(unitcell[0]*unitcell[2]);
      double cosAC = (a*c)/(unitcell[0]*unitcell[5]);
      double cosBC = (b*c)/(unitcell[2]*unitcell[5]);
      if (cosAB > 1.0) cosAB = 1.0; else if (cosAB < -1.0) cosAB = -1.0;
      if (cosAC > 1.0) cosAC = 1.0; else if (cosAC < -1.0) cosAC = -1.0;
      if (cosBC > 1.0) cosBC = 1.0; else if (cosBC < -1.0) cosBC = -1.0;
      unitcell[1] = cosAB;
      unitcell[3] = cosAC;
      unitcell[4] = cosBC;
    } else {
      unitcell[0] = unitcell[2] = unitcell[5] = 1.0;
      unitcell[1] = unitcell[3] = unitcell[4] = 0.0;
    }
    ret_code = write_dcdstep(fileid, n, x, y, z, unitcell);
  } else {
    ret_code = write_dcdstep(fileid, n, x, y, z, NULL);
  }
  if (ret_code < 0)
  {
    NAMD_err("Writing of DCD step failed!!");
  }

}
/*      END OF FUNCTION output_dcdfile      */

/************************************************************************/
/*                  */
/*      FUNCTION output_final_coordinates    */
/*                  */
/*   INPUTS:                */
/*  coor - Array of vectors containing final coordinates    */
/*  n - Number of coordinates to output        */
/*  timestep - Timestep that coordinates are being written in  */
/*                  */
/*  This function writes out the final coordinates for the    */
/*   simulation in PDB format to the file specified in the config  */
/*   file.                */
/*                  */
/************************************************************************/

void Output::output_final_coordinates(Vector *coor, int n, int timestep)

{
  char output_name[140];  //  Output filename
  char comment[128];    //  comment for PDB header

  //  Built the output filename
  strcpy(output_name, namdMyNode->simParams->outputFilename);
  strcat(output_name, ".coor");

  NAMD_backup_file(output_name);

  //  Check to see if we should write out a binary file or a
  //  PDB file
  if (!namdMyNode->simParams->binaryOutput)
  {
    sprintf(comment, "FINAL COORDINATES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    namdMyNode->pdbData->set_all_positions(coor);
    namdMyNode->pdbData->write(output_name, comment);
  }
  else
  {
    //  Write the velocities to a binary file
    write_binary_file(output_name, n, coor);
  }
}
/*    END OF FUNCTION output_final_coordinates    */

/************************************************************************/
/*                  */
/*      FUNCTION output_final_velocities    */
/*                  */
/*   INPUTS:                */
/*  vel - Array of vectors containing final velocities    */
/*  timestep - Timestep that vleocities are being written in  */
/*                  */
/*  This function writes out the final vleocities for the    */
/*   simulation in PDB format to the file specified in the config  */
/*   file.                */
/*                  */
/************************************************************************/

void Output::output_final_velocities(int timestep, int n, Vector *vel)

{
  char output_name[140];  //  Output filename
  char comment[128];    //  Comment for PDB header

  //  Build the output filename
  strcpy(output_name, namdMyNode->simParams->outputFilename);
  strcat(output_name, ".vel");

  NAMD_backup_file(output_name);

  //  Check to see if we should write a PDB or binary file
  if (!(namdMyNode->simParams->binaryOutput))
  {
    //  Write the final velocities to a PDB file
    sprintf(comment, "FINAL VELOCITIES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

    scale_vels(vel, n, PDBVELFACTOR);
    namdMyNode->pdbData->set_all_positions(vel);
    namdMyNode->pdbData->write(output_name, comment);
    scale_vels(vel, n, PDBVELINVFACTOR);
  }
  else
  {
    //  Write the coordinates to a binary file
    write_binary_file(output_name, n, vel);
  }

}
/*      END OF FUNCTION output_final_velocities    */

/************************************************************************/
/*                  */
/*      FUNCTION output_veldcdfile      */
/*                  */
/*   INPUTS:                */
/*  timestep - Current timestep          */
/*  n - Number of atoms in simulation        */
/*  coor - velocity vectors for all atoms        */
/*                  */
/*  This function maintains the interface between the Output object */
/*   and the dcd writing routines contained in dcdlib.  This fucntion   */
/*   writes out the velocity vectors in DCD format.      */
/*                  */
/************************************************************************/

void Output::output_veldcdfile(int timestep, int n, Vector *vel)

{
  static Bool first=TRUE;  //  Flag indicating first call
  static int fileid;  //  File id for the dcd file
  static float *x, *y, *z; // Arrays to hold x, y, and z arrays
  int i;      //  Loop counter
  int ret_code;    //  Return code from DCD calls
  SimParameters *simParams = Node::Object()->simParameters;

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( timestep == END_OF_RUN ) {
    if ( ! first ) {
      iout << "CLOSING VELOCITY DCD FILE\n" << endi;
      close_dcd_write(fileid);
    } else {
      iout << "VELOCITY DCD FILE WAS NOT CREATED\n" << endi;
    }
    return;
  }

  if (first)
  {
    //  Allocate x, y, and z arrays since the DCD file routines
    //  need them passed as three independant arrays to be
    //  efficient
    x = new float[n];
    y = new float[n];
    z = new float[n];

    if ( (x==NULL) || (y==NULL) || (z==NULL) )
    {
      NAMD_err("memory allocation failed in Output::output_veldcdfile");
    }

    //  Open the DCD file
    iout << "OPENING VELOCITY DCD FILE\n" << endi;

    fileid=open_dcd_write(namdMyNode->simParams->velDcdFilename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "Velocity DCD file %s already exists!!",
        namdMyNode->simParams->velDcdFilename);

      NAMD_err(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open velocity DCD file %s",
        namdMyNode->simParams->velDcdFilename);

      NAMD_err(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->velDcdFrequency;
    NSTEP = NSAVC * (simParams->N/NSAVC);
    NPRIV = simParams->firstTimestep+NSAVC;
    NPRIV = NSAVC * (NPRIV/NSAVC);
    NFILE = (NSTEP-NPRIV)/NSAVC + 1;

    //  Write out the header
    const int with_unitcell = 0;
    ret_code = write_dcdheader(fileid, 
        simParams->velDcdFilename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR, with_unitcell);


    if (ret_code<0)
    {
      NAMD_err("Writing of velocity DCD header failed!!");
    }

    first = FALSE;
  }

  //  Copy the coordinates for output
  for (i=0; i<n; i++)
  {
    x[i] = vel[i].x;
    y[i] = vel[i].y;
    z[i] = vel[i].z;
  }

  //  Write out the values for this timestep
  iout << "WRITING VELOCITIES TO DCD FILE AT STEP "
	<< timestep << "\n" << endi;
  ret_code = write_dcdstep(fileid, n, x, y, z, NULL);

  if (ret_code < 0)
  {
    NAMD_err("Writing of velocity DCD step failed!!");
  }

}
/*      END OF FUNCTION output_veldcdfile    */

/************************************************************************/
/*                  */
/*      FUNCTION write_binary_file      */
/*                  */
/*   INPUTS:                */
/*  fname - file name to write velocities to      */
/*  n - Number of atoms in system          */
/*  vels - Array of vectors            */
/*                  */
/*  This function writes out vectors in binary format to    */
/*   the specified file.            */
/*                  */
/************************************************************************/

void Output::write_binary_file(char *fname, int n, Vector *vecs)

{
  FILE *fp;    //  File descriptor
  int32 n32 = n;

  //  open the file and die if the open fails
  if ( (fp = fopen(fname, "wb")) == NULL)
  {
    char errmsg[256];

    sprintf(errmsg, "Unable to open binary file %s", fname);
    NAMD_err(errmsg);
  }

  //  Write out the number of atoms and the vectors
  fwrite(&n32, sizeof(int32), 1, fp);
  fwrite(vecs, sizeof(Vector), n, fp);

  if ( ferror(fp) || fclose(fp) )
  {
    char errmsg[256];

    sprintf(errmsg, "Error on write to binary file %s", fname);
    NAMD_err(errmsg);
  }
}
/*      END OF FUNCTION write_binary_file    */

/************************************************************************/
/*                  */
/*      FUNCTION scale_vels        */
/*                  */
/*   INPUTS:                */
/*  v - Array of velocity vectors          */
/*  n - Number of atoms in system          */
/*  fact - Scaling factor            */
/*                  */
/*  This function scales all the vectors passed in by a constant  */
/*   factor.  This is used before writing out velocity vectors that  */
/*   need to be resized.            */
/*                  */
/************************************************************************/

void Output::scale_vels(Vector *v, int n, Real fact)

{
  int i;

  for (i=0; i<n; i++)
  {
    v[i].x *= fact;
    v[i].y *= fact;
    v[i].z *= fact;
  }
}
/*      END OF FUNCTION scale_vels      */

