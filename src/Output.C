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
#include "InfoStream.h"

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

void Output::coordinate(int timestep, int n, Vector *coor, FloatVector *fcoor)
{
  SimParameters *simParams = Node::Object()->simParameters;

  if ( timestep >= 0 ) {

    //  Output a DCD trajectory 
    if ( simParams->dcdFrequency &&
       ((timestep % simParams->dcdFrequency) == 0) )
    {
      output_dcdfile(timestep, n, fcoor);
    }

    //  Output a restart file
    if ( simParams->restartFrequency &&
       ((timestep % simParams->restartFrequency) == 0) )
    {
      iout << "WRITING COORDINATES TO RESTART FILE AT STEP "
				<< timestep << "\n" << endi;
      output_restart_coordinates(coor, n, timestep);
    }

    //  Interactive MD
    if ( simParams->IMDon &&
       ( ((timestep % simParams->IMDfreq) == 0) ||
         (timestep == simParams->firstTimestep) ) )
    {
      IMDOutput *imd = Node::Object()->imd;
      if (imd != NULL) imd->gather_coordinates(timestep, n, fcoor);
    }

  }

  if (timestep == EVAL_MEASURE)
  {
#ifdef NAMD_TCL
    Node::Object()->getScript()->measure(coor);
#endif
  }

  //  Output final coordinates
  if (timestep == FILE_OUTPUT || timestep == END_OF_RUN)
  {
    iout << "WRITING COORDINATES TO OUTPUT FILE AT STEP "
				<< simParams->N << "\n" << endi;
    output_final_coordinates(coor, n, simParams->N);
  }

  //  Close trajectory files
  if (timestep == END_OF_RUN)
  {
    if (simParams->dcdFrequency) output_dcdfile(END_OF_RUN,0,0);
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
  static Bool first=TRUE;    //  Is this the first time function
          //  has been called?
  static char *restart_name,      //  filenames for restart and backup
        *restart_bak;  //  restart files
  char comment[128];    //  Comment for header of PDB file

  //  Check to see if this is the first time the function has been
  //  called
  if (first)
  {
    //  This is the first invocation of the function, so build
    //  the restart and backup filenames
    restart_name = new char[strlen(namdMyNode->simParams->restartFilename)+6];
    restart_bak = new char[strlen(namdMyNode->simParams->restartFilename)+10];

    if ( (restart_name == NULL) || (restart_bak == NULL) )
    {
      NAMD_die("memory allocation failed in Output::coordinate");
    }

    strcpy(restart_name, namdMyNode->simParams->restartFilename);
    strcat(restart_name, ".coor");
    strcpy(restart_bak, namdMyNode->simParams->restartFilename);
    strcat(restart_bak, ".coor.old");

    first=FALSE;
  }
#ifdef WIN32
  remove(restart_bak);
#endif
  rename(restart_name, restart_bak);

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
  static Bool first=TRUE;    //  is this the first invocation
          //  of the function?
  static char *restart_name,      //  filenames for the restart and
        *restart_bak;  //  backup files
  char comment[128];    //  comment for the header of the
          //  PDB file

  //  Check to see if this is the first invocation
  if (first)
  {
    //  This is the first call to this function, so set up
    //  the names that are to be used for the restart and
    //  backup files
    restart_name = new char[strlen(namdMyNode->simParams->restartFilename)+5];
    restart_bak = new char[strlen(namdMyNode->simParams->restartFilename)+9];

    if ( (restart_name == NULL) || (restart_bak == NULL) )
    {
      NAMD_die("memory allocation failed in Output::output_restart_velocites");
    }

    strcpy(restart_name, namdMyNode->simParams->restartFilename);
    strcat(restart_name, ".vel");
    strcpy(restart_bak, namdMyNode->simParams->restartFilename);
    strcat(restart_bak, ".vel.old");

    first=FALSE;
  }
#ifdef WIN32
  remove(restart_bak);
#endif
  rename(restart_name, restart_bak);

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
/*                  */
/*  This function maintains the interface between the Output object */
/*   and the dcd writing routines contained in dcdlib.      */
/*                  */
/************************************************************************/

void Output::output_dcdfile(int timestep, int n, FloatVector *coor)

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
      NAMD_die("memory allocation failed in Output::output_dcdfile");
    }

    //  Open the DCD file
    iout << "OPENING COORDINATE DCD FILE\n" << endi;

    fileid=open_dcd_write(simParams->dcdFilename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "DCD file %s already exists!!",
        simParams->dcdFilename);

      NAMD_die(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open DCD file %s",
        simParams->dcdFilename);

      NAMD_die(err_msg);
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
        simParams->dt/TIMEFACTOR);


    if (ret_code<0)
    {
      NAMD_die("Writing of DCD header failed!!");
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
  ret_code = write_dcdstep(fileid, n, x, y, z);

  if (ret_code < 0)
  {
    NAMD_die("Writing of DCD step failed!!");
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

  char bfname[140];
  strcpy(bfname,output_name);
  strcat(bfname,".BAK");
#ifdef WIN32
  remove(bfname);
#endif
  rename(output_name,bfname);

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

  char bfname[140];
  strcpy(bfname,output_name);
  strcat(bfname,".BAK");
#ifdef WIN32
  remove(bfname);
#endif
  rename(output_name,bfname);

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
      NAMD_die("memory allocation failed in Output::output_veldcdfile");
    }

    //  Open the DCD file
    iout << "OPENING VELOCITY DCD FILE\n" << endi;

    fileid=open_dcd_write(namdMyNode->simParams->velDcdFilename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "Velocity DCD file %s already exists!!",
        namdMyNode->simParams->velDcdFilename);

      NAMD_die(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open velocity DCD file %s",
        namdMyNode->simParams->velDcdFilename);

      NAMD_die(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->velDcdFrequency;
    NSTEP = NSAVC * (simParams->N/NSAVC);
    NPRIV = simParams->firstTimestep+NSAVC;
    NPRIV = NSAVC * (NPRIV/NSAVC);
    NFILE = (NSTEP-NPRIV)/NSAVC + 1;

    //  Write out the header
    ret_code = write_dcdheader(fileid, 
        simParams->velDcdFilename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR);


    if (ret_code<0)
    {
      NAMD_die("Writing of velocity DCD header failed!!");
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
  ret_code = write_dcdstep(fileid, n, x, y, z);

  if (ret_code < 0)
  {
    NAMD_die("Writing of velocity DCD step failed!!");
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
    NAMD_die(errmsg);
  }

  //  Write out the number of atoms and the vectors
  fwrite(&n32, sizeof(int32), 1, fp);
  fwrite(vecs, sizeof(Vector), n, fp);

  if ( ferror(fp) || fclose(fp) )
  {
    char errmsg[256];

    sprintf(errmsg, "Error on write to binary file %s", fname);
    NAMD_die(errmsg);
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

