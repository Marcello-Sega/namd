/***************************************************************************/
/*      (C) Copyright 1995,1996,1997 The Board of Trustees of the          */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 * This object outputs the data collected on the master node
 ***************************************************************************/

#include <string.h>
#include <stdlib.h>

#include "IMDOutput.h"
#include "Output.h"
#include "dcdlib.h"
#include "strlib.h"
#include "Inform.h"
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

  //  Output a DCD trajectory 
  if ( (simParams->dcdFrequency != -1) &&
       ((timestep % simParams->dcdFrequency) == 0) &&
       (timestep != simParams->firstTimestep) )
  {
    positionsNeeded |= 1;
  }

  //  Output a restart file
  if ( (simParams->restartFrequency != -1) &&
       ((timestep % simParams->restartFrequency) == 0) &&
       (timestep != simParams->firstTimestep) )
  {
    positionsNeeded |= 2;
  }

  //  Output final coordinates
  if (timestep == 0)
  {
    positionsNeeded |= 2;
  }

  if (simParams->IMDon) 
    if (( (timestep % simParams->IMDfreq)== 0) ||
        ( timestep == simParams->firstTimestep)) 
      positionsNeeded |= 1;
 
  return positionsNeeded;
}

void Output::coordinate(int timestep, int n, Vector *coor, FloatVector *fcoor)
{
  SimParameters *simParams = Node::Object()->simParameters;

  //  Output a DCD trajectory 
  if ( (simParams->dcdFrequency != -1) &&
       ((timestep % simParams->dcdFrequency) == 0) &&
       (timestep != simParams->firstTimestep) )
  {
    iout << "WRITING COORDINATES TO DCD FILE AT STEP "
				<< timestep << "\n" << endi;
    output_dcdfile(timestep, n, fcoor);
  }

  //  Output a restart file
  if ( (simParams->restartFrequency != -1) &&
       ((timestep % simParams->restartFrequency) == 0) &&
       (timestep != simParams->firstTimestep) )
  {
    iout << "WRITING COORDINATES TO RESTART FILE AT STEP "
				<< timestep << "\n" << endi;
    output_restart_coordinates(coor, n, timestep);
  }

  //  Output final coordinates
  if (timestep == 0)
  {
    iout << "WRITING COORDINATES TO OUTPUT FILE AT STEP "
				<< simParams->N << "\n" << endi;
    output_final_coordinates(coor, n, simParams->N);
  }

  if (simParams->IMDon) {
    if (( (timestep % simParams->IMDfreq)==0) ||
        ( timestep == simParams->firstTimestep)) {
      IMDOutput *imd = Node::Object()->imd;
      if (imd != NULL)
        imd->gather_coordinates(timestep, n, fcoor);
      }
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

  //  Output a velocity DCD trajectory
  if ( (simParams->velDcdFrequency != -1) &&
       ((timestep % simParams->velDcdFrequency) == 0) &&
       (timestep != simParams->firstTimestep) )
  {
    velocitiesNeeded |= 1;
  }

  //  Output a restart file
  if ( (simParams->restartFrequency != -1) &&
       ((timestep % simParams->restartFrequency) == 0) &&
       (timestep != simParams->firstTimestep) )
  {
    velocitiesNeeded |= 2;
  }

  //  Output final velocities
  if (timestep == 0)
  {
    velocitiesNeeded |= 2;
  }

  return velocitiesNeeded;
}

void Output::velocity(int timestep, int n, Vector *vel)
{
  SimParameters *simParams = Node::Object()->simParameters;

  //  Output velocity DCD trajectory
  if ( (simParams->velDcdFrequency != -1) &&
       ((timestep % simParams->velDcdFrequency) == 0) &&
       (timestep != simParams->firstTimestep) )
  {
    iout << "WRITING VELOCITIES TO DCD FILE AT STEP "
				<< timestep << "\n" << endi;
    output_veldcdfile(timestep, n, vel);
  }

  //  Output restart file
  if ( (simParams->restartFrequency != -1) &&
       ((timestep % simParams->restartFrequency) == 0) &&
       (timestep != simParams->firstTimestep) )
  {
    iout << "WRITING VELOCITIES TO RESTART FILE AT STEP "
				<< timestep << "\n" << endi;
    output_restart_velocities(timestep, n, vel);
  }

  //  Output final velocities
  if (timestep == 0)
  {
    iout << "WRITING VELOCITIES TO OUTPUT FILE AT STEP "
				<< simParams->N << "\n" << endi;
    output_final_velocities(simParams->N, n, vel);
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
    strcat(restart_bak, ".coor.BAK");

    first=FALSE;
  }
  else
  {
    //  This is not the first invocation, so move the current version
    //  of the file to the backup filename
    rename(restart_name, restart_bak);
  }

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
    strcat(restart_bak, ".vel.BAK");

    first=FALSE;
  }
  else
  {
    //  This is *not* the first call, so move the current version
    //  of the restart file to the backup file before we write out
    //  the new restart file
    rename(restart_name, restart_bak);
  }

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
  ret_code = write_dcdstep(fileid, n, x, y, z);

  if (ret_code < 0)
  {
    NAMD_die("Writing of DCD step failed!!");
  }

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( (timestep+simParams->dcdFrequency) >
        simParams->N)
  {
    close_dcd_write(fileid);
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
  static char output_name[140];  //  Output filename
  char comment[128];    //  comment for PDB header

  //  Built the output filename
  strcpy(output_name, namdMyNode->simParams->outputFilename);
  strcat(output_name, ".coor");

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
  static char output_name[140];  //  Output filename
  char comment[128];    //  Comment for PDB header

  //  Build the output filename
  strcpy(output_name, namdMyNode->simParams->outputFilename);
  strcat(output_name, ".vel");

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
  ret_code = write_dcdstep(fileid, n, x, y, z);

  if (ret_code < 0)
  {
    NAMD_die("Writing of velocity DCD step failed!!");
  }

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( (timestep+namdMyNode->simParams->velDcdFrequency) >
        namdMyNode->simParams->N)
  {
    close_dcd_write(fileid);
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
  if ( (fp = fopen(fname, "w")) == NULL)
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

/********************************************************************************/
/*                    */
/*        FUNCTION long_force        */
/*                    */
/*   INPUTS:                  */
/*  tiemstep - Current timestep value          */
/*  n - Number of forces              */
/*  f - Array of force vectors            */
/*                    */
/*  This is the public method that is used to take the forces from the  */
/*   collect object and output them into a force DCD.  It does nothing more  */
/*   than call the force DCD generation routine.        */
/*                    */
/********************************************************************************/

void Output::long_force(int timestep, int n, Vector *f)
   
{
   output_longforcedcdfile(timestep, n, f);
}
/*      END OF FUNCTION long_force        */

/********************************************************************************/
/*                    */
/*        FUNCTION short_force        */
/*                    */
/*   INPUTS:                  */
/*  tiemstep - Current timestep value          */
/*  n - Number of forces              */
/*  f - Array of force vectors            */
/*                    */
/*  This is the public method that is used to take the forces from the  */
/*   collect object and output them into a force DCD.  It does nothing more  */
/*   than call the force DCD generation routine.        */
/*                    */
/********************************************************************************/

void Output::short_force(int timestep, int n, Vector *f)
   
{
   output_shortforcedcdfile(timestep, n, f);
}
/*      END OF FUNCTION force          */

/********************************************************************************/
/*                    */
/*        FUNCTION all_force        */
/*                    */
/*   INPUTS:                  */
/*  tiemstep - Current timestep value          */
/*  n - Number of forces              */
/*  f - Array of force vectors            */
/*                    */
/*  This is the public method that is used to take the forces from the  */
/*   collect object and output them into a force DCD.  It does nothing more  */
/*   than call the force DCD generation routine.        */
/*                    */
/********************************************************************************/

void Output::all_force(int timestep, int n, Vector *f)
   
{
   output_allforcedcdfile(timestep, n, f);
}
/*      END OF FUNCTION all_force        */

/********************************************************************************/
/*                    */
/*      FUNCTION output_longforcedcdfile      */
/*                    */
/*   INPUTS:                  */
/*  timestep - Current timestep value          */
/*  n - Number of forces              */
/*  forces - Array of force vectors            */
/*                    */
/*  Output forces to a force DCD file.          */
/*                    */
/********************************************************************************/

void Output::output_longforcedcdfile(int timestep, int n, Vector *forces)

{
  static Bool first=TRUE;  //  Flag indicating first call
  static int fileid;  //  File id for the dcd file
  static float *x, *y, *z; // Arrays to hold x, y, and z arrays
  int i;      //  Loop counter
  int ret_code;    //  Return code from DCD calls
  char filename[257];  //  DCD filename
  SimParameters *simParams = Node::Object()->simParameters;

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
      NAMD_die("memory allocation failed in Output::output_longforcedcdfile");
    }

    //  Open the DCD file
    sprintf(filename, "%s.l", namdMyNode->simParams->electForceDcdFilename);
    fileid=open_dcd_write(filename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "Force DCD file %s already exists!!",
        filename);

      NAMD_die(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open force DCD file %s",
        filename);

      NAMD_die(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->electForceDcdFrequency;
    NSTEP = NSAVC * (simParams->N/NSAVC);
    NPRIV = simParams->firstTimestep+NSAVC;
    NPRIV = NSAVC * (NPRIV/NSAVC);
    NFILE = (NSTEP-NPRIV)/NSAVC + 1;

    //  Write out the header
    ret_code = write_dcdheader(fileid, 
        filename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR);


    if (ret_code<0)
    {
      NAMD_die("Writing of long range elect force DCD header failed!!");
    }

    first = FALSE;
  }

  //  Copy the coordinates for output
  for (i=0; i<n; i++)
  {
    x[i] = forces[i].x;
    y[i] = forces[i].y;
    z[i] = forces[i].z;
  }

  //  Write out the values for this timestep
  ret_code = write_dcdstep(fileid, n, x, y, z);
  
  
  if (ret_code < 0)
  {
    NAMD_die("Writing of long range elect force DCD step failed!!");
  }

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( (timestep+namdMyNode->simParams->electForceDcdFrequency) >
        namdMyNode->simParams->N)
  {
    close_dcd_write(fileid);
  }
}
/*      END OF FUNCTION output_longforcedcdfile      */

/********************************************************************************/
/*                    */
/*      FUNCTION output_shortforcedcdfile      */
/*                    */
/*   INPUTS:                  */
/*  timestep - Current timestep value          */
/*  n - Number of forces              */
/*  forces - Array of force vectors            */
/*                    */
/*  Output forces to a short force DCD file.        */
/*                    */
/********************************************************************************/

void Output::output_shortforcedcdfile(int timestep, int n, Vector *forces)

{
  static Bool first=TRUE;  //  Flag indicating first call
  static int fileid;  //  File id for the dcd file
  static float *x, *y, *z; // Arrays to hold x, y, and z arrays
  int i;      //  Loop counter
  int ret_code;    //  Return code from DCD calls
  char filename[257];
  SimParameters *simParams = Node::Object()->simParameters;

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
      NAMD_die("memory allocation failed in Output::output_shortforcedcdfile");
    }

    //  Open the DCD file
    sprintf(filename, "%s.s", namdMyNode->simParams->electForceDcdFilename);
    fileid=open_dcd_write(filename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "Short range force DCD file %s already exists!!",
        filename);

      NAMD_die(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open short range force DCD file %s",
        filename);

      NAMD_die(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->electForceDcdFrequency;
    NSTEP = NSAVC * (simParams->N/NSAVC);
    NPRIV = simParams->firstTimestep+NSAVC;
    NPRIV = NSAVC * (NPRIV/NSAVC);
    NFILE = (NSTEP-NPRIV)/NSAVC + 1;

    //  Write out the header
    ret_code = write_dcdheader(fileid, 
        filename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR);


    if (ret_code<0)
    {
      NAMD_die("Writing of short range force DCD header failed!!");
    }

    first = FALSE;
  }

  //  Copy the coordinates for output
  for (i=0; i<n; i++)
  {
    x[i] = forces[i].x;
    y[i] = forces[i].y;
    z[i] = forces[i].z;
  }

  //  Write out the values for this timestep
  ret_code = write_dcdstep(fileid, n, x, y, z);
  
  if (ret_code < 0)
  {
    NAMD_die("Writing of short range force DCD step failed!!");
  }

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( (timestep+namdMyNode->simParams->electForceDcdFrequency) >
        namdMyNode->simParams->N)
  {
    close_dcd_write(fileid);
  }
}
/*      END OF FUNCTION output_shortforcedcdfile    */

/********************************************************************************/
/*                    */
/*      FUNCTION output_allforcedcdfile        */
/*                    */
/*   INPUTS:                  */
/*  timestep - Current timestep value          */
/*  n - Number of forces              */
/*  forces - Array of force vectors            */
/*                    */
/*  Output forces to a force DCD file.          */
/*                    */
/********************************************************************************/

void Output::output_allforcedcdfile(int timestep, int n, Vector *forces)

{
  static Bool first=TRUE;  //  Flag indicating first call
  static int fileid;  //  File id for the dcd file
  static float *x, *y, *z; // Arrays to hold x, y, and z arrays
  int i;      //  Loop counter
  int ret_code;    //  Return code from DCD calls
  SimParameters *simParams = Node::Object()->simParameters;

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
      NAMD_die("memory allocation failed in Output::output_allforcedcdfile");
    }

    //  Open the DCD file
    fileid=open_dcd_write(namdMyNode->simParams->allForceDcdFilename);

    if (fileid == DCD_FILEEXISTS)
    {
      char err_msg[257];

      sprintf(err_msg, "Force DCD file %s already exists!!",
        namdMyNode->simParams->allForceDcdFilename);

      NAMD_die(err_msg);
    }
    else if (fileid < 0)
    {
      char err_msg[257];

      sprintf(err_msg, "Couldn't open force DCD file %s",
        namdMyNode->simParams->allForceDcdFilename);

      NAMD_die(err_msg);
    }

    int NSAVC, NFILE, NPRIV, NSTEP;
    NSAVC = simParams->allForceDcdFrequency;
    NSTEP = NSAVC * (simParams->N/NSAVC);
    NPRIV = simParams->firstTimestep+NSAVC;
    NPRIV = NSAVC * (NPRIV/NSAVC);
    NFILE = (NSTEP-NPRIV)/NSAVC + 1;

    //  Write out the header
    ret_code = write_dcdheader(fileid, 
        simParams->allForceDcdFilename,
        n, NFILE, NPRIV, NSAVC, NSTEP,
        simParams->dt/TIMEFACTOR);


    if (ret_code<0)
    {
      NAMD_die("Writing of force DCD header failed!!");
    }

    first = FALSE;
  }

  //  Copy the coordinates for output
  for (i=0; i<n; i++)
  {
    x[i] = forces[i].x;
    y[i] = forces[i].y;
    z[i] = forces[i].z;
  }

  //  Write out the values for this timestep
  ret_code = write_dcdstep(fileid, n, x, y, z);
  
  
  if (ret_code < 0)
  {
    NAMD_die("Writing of all force DCD step failed!!");
  }

  //  If this is the last time we will be writing coordinates,
  //  close the file before exiting
  if ( (timestep+namdMyNode->simParams->allForceDcdFrequency) >
        namdMyNode->simParams->N)
  {
    close_dcd_write(fileid);
  }
}
/*      END OF FUNCTION output_allforcedcdfile      */


/***************************************************************************
 * RCS INFORMATION:
 *
 *  $RCSfile: Output.C,v $
 *  $Author: jim $  $Locker:  $    $State: Exp $
 *  $Revision: 1.24 $  $Date: 1999/09/12 19:33:18 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Output.C,v $
 * Revision 1.24  1999/09/12 19:33:18  jim
 * Collections now use floats when possible.
 *
 * Revision 1.23  1999/09/02 23:04:50  justin
 * Eliminated MDComm from all files and Makefiles
 *
 * Revision 1.22  1999/08/20 21:58:18  jim
 * Changes from Justin for interactive molecular dynamics code.
 *
 * Revision 1.21  1999/08/16 22:19:40  jim
 * Incorporated Justin's interactive MD code.
 *
 * Revision 1.20  1999/06/21 16:15:33  jim
 * Improved scripting, run now ends and generates output.
 *
 * Revision 1.19  1999/05/25 21:48:48  jim
 * Modified DCD code to be compatible with Quanta.
 *
 * Revision 1.18  1999/05/11 23:56:40  brunner
 * Changes for new charm version
 *
 * Revision 1.17  1998/10/25 04:42:15  krishnan
 * The format of the output files now depend on the binaryOutput parameter
 * instead of binaryRestart
 *
 * Revision 1.16  1998/10/24 19:57:45  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.15  1998/09/02 20:38:31  jim
 * Improved error checking on file output.
 *
 * Revision 1.14  1998/04/30 04:53:28  jim
 * Added forces from MDComm and other improvements to ComputeGlobal.
 *
 * Revision 1.13  1998/04/14 03:19:22  jim
 * Fixed up MDCOMM code.
 *
 * Revision 1.12  1997/11/10 16:45:56  milind
 * Made comm a Cpv Variable.
 *
 * Revision 1.11  1997/10/01 16:46:59  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.10  1997/09/22 20:24:48  brunner
 * Corrected velocity PDB output by replacing conversion factor of 20 with
 * conversion for Ang/ps.  This now shoud match XPLOR.
 *
 * Revision 1.9  1997/09/18 22:22:21  jim
 * Cleaned up coordinate ane velocity output code a little.
 * Won't create restart file on first time step.
 *
 * Revision 1.8  1997/09/18 21:05:12  brunner
 * Made DCD files no update on first timestep, if firstTimestep
 * is non-zero.
 *
 * Revision 1.7  1997/08/13 21:00:18  brunner
 * Made binary files always use 32 bits for the number of atoms, so that it
 * works on both 64 and 32-bit machines.  Also, I made Inform.C use CkPrintf,
 * to fix the I/O buffering.
 *
 * Revision 1.6  1997/04/10 17:28:45  brunner
 * Made time output in energy actually work.  The one in output doesn't do
 * anything.
 *
 * Revision 1.4  1997/03/19 11:54:36  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.2  1997/02/26 18:38:24  jim
 * Eliminated +1 from output timestep checks, now makes sense.
 * This should match changes being made in NAMD 1.X.
 *
 * Revision 1.1  1997/02/11 22:56:15  jim
 * Added dcd file writing.
 *
 * Revision 1.44  1997/01/03 22:13:55  brunner
 * Added momentum output
 *
 * Revision 1.43  1996/12/03 18:46:48  brunner
 * Put cylindrical BC into Constraint output field
 *
 * Revision 1.42  1996/10/24 19:00:22  brunner
 * Made timestep go from 0->N, inclusive
 *
 * Revision 1.41  1996/09/13 15:14:52  nealk
 * Changed structures to use array atom[] rather than explicit atom1, 
 * atom2, etc.
 *
 * Revision 1.40  1996/05/08 20:55:23  gursoy
 * modified calculation of temperature
 * to acocound number of degrees of freedom
 *
 * Revision 1.39  1996/04/27 23:40:59  billh
 * Modified to calculate and print hydrogen bond forces and energies.
 *
 * Revision 1.38  1996/01/28 21:50:08  jean
 * Attempting to make stable RCS without Mark Nelson's
 * fma/pairlist decoupling
 *
 * Revision 1.38  1995/12/05 22:13:34  brunner
 * Rearranged #ifdef MDCOMM 's so it will compile without MDCOMM defined
 *
 * Revision 1.37  1995/10/27 12:24:14  dalke
 * Removed bad delete's
 *
 * Revision 1.36  95/10/26  16:24:27  16:24:27  dalke (Andrew Dalke)
 * added new mdcomm code
 * 
 * Revision 1.35  1995/10/09  18:32:31  hazen
 * Updated memory allocation to use C++ new/delete
 * and freed temporary memory allocated in arrays X, Y, Z in
 * functions output_* and strings restart_name and restart_bak.
 * Also freed vmdData struct members X, Y, Z and others in destructor.
 *
 * Revision 1.34  1995/09/30  19:55:14  billh
 * Added ETITLE and ENERGY strings  in energy output records.
 * Added tiny Elvis warning if temperature exceeds 1000 degrees.
 *
 * Revision 1.33  95/09/26  15:15:04  15:15:04  nelson (Mark T. Nelson)
 * Added all force DCD files and cleaned up symantics of long and short
 * range electrostatic force DCD files.
 * 
 * Revision 1.32  95/09/26  13:28:21  13:28:21  nelson (Mark T. Nelson)
 * Added broadcast of temperature for temperature coupling
 * 
 * Revision 1.31  95/08/31  10:45:41  10:45:41  nelson (Mark T. Nelson)
 * Changed so that DCD output write during both first and last timestep
 * 
 * Revision 1.30  95/08/30  14:03:09  14:03:09  nelson (Mark T. Nelson)
 * Added short range DCD file output and binary coordinate restart files
 * 
 * Revision 1.29  95/08/11  14:52:17  14:52:17  nelson (Mark T. Nelson)
 * Added routines for the generation of force DCD files
 * 
 * Revision 1.28  95/07/14  14:19:08  14:19:08  nelson (Mark T. Nelson)
 * Added binary velocity restart files
 * 
 * Revision 1.27  95/07/11  09:18:15  09:18:15  nelson (Mark T. Nelson)
 * Changed scaling factors for velocity restart files so that it matches
 * X-PLOR files
 * 
 * Revision 1.26  95/07/05  11:48:56  11:48:56  nelson (Mark T. Nelson)
 * Changed scaling factor on velocity restart files to match X-PLOR
 * 
 * Revision 1.25  95/06/20  11:34:33  11:34:33  nelson (Mark T. Nelson)
 * Added explicit int declaration to get rid of warning messages
 * 
 * Revision 1.24  95/05/23  14:11:51  14:11:51  nelson (Mark T. Nelson)
 * Added electric field energy values
 * 
 * Revision 1.23  95/04/10  13:35:43  13:35:43  nelson (Mark T. Nelson)
 * Removed some nulls to silence gcc
 * 
 * Revision 1.22  95/03/22  11:22:59  11:22:59  nelson (Mark T. Nelson)
 * Added output opition for sphereical boundary conditions and 
 * added outputEnergies option
 * 
 * Revision 1.21  95/03/08  14:47:32  14:47:32  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.20  95/02/17  12:37:16  12:37:16  nelson (Mark T. Nelson)
 * Added shifting factor for velocity PDB files
 * 
 * Revision 1.19  95/01/31  19:46:07  19:46:07  nelson (Mark T. Nelson)
 * Added call for rescaling velocities
 * 
 * Revision 1.18  94/11/23  09:38:04  09:38:04  nelson (Mark T. Nelson)
 * Added functions for outputing VMD patch load data
 * 
 * Revision 1.17  94/11/10  20:31:43  20:31:43  nelson (Mark T. Nelson)
 * Added code to send Patch data to VMD
 * 
 * Revision 1.16  94/11/09  13:30:31  13:30:31  nelson (Mark T. Nelson)
 * Removed debugging statement in vmd output function
 * 
 * Revision 1.15  94/11/09  13:28:18  13:28:18  nelson (Mark T. Nelson)
 * Changed elapsed_cpu to elapsed_time
 * 
 * Revision 1.14  94/10/31  19:25:21  19:25:21  nelson (Mark T. Nelson)
 * Added functions for VMD interface
 * 
 * Revision 1.13  94/10/19  21:39:59  21:39:59  nelson (Mark T. Nelson)
 * Added harmonic constraint energy output
 * 
 * Revision 1.12  94/10/18  20:32:42  20:32:42  nelson (Mark T. Nelson)
 * Added functions to output velocity DCD files
 * 
 * Revision 1.11  94/10/14  09:36:49  09:36:49  nelson (Mark T. Nelson)
 * Added output of final coordinate and velocity files as well
 * as rearranging and commenting other output routines
 * 
 * Revision 1.10  94/10/12  15:27:57  15:27:57  nelson (Mark T. Nelson)
 * Added comment line to restart files
 * 
 * Revision 1.9  94/10/12  15:18:44  15:18:44  nelson (Mark T. Nelson)
 * Added restart file capability
 * 
 * Revision 1.8  94/10/07  12:19:25  12:19:25  nelson (Mark T. Nelson)
 * Added temperature output
 * 
 * Revision 1.7  94/10/06  21:37:15  21:37:15  nelson (Mark T. Nelson)
 * added dcd file output
 * 
 * Revision 1.6  94/10/06  15:51:36  15:51:36  gursoy (Attila Gursoy)
 * coordinates and velocities are passed as an array(it was a message)
 * 
 * Revision 1.5  94/09/30  11:09:41  11:09:41  nelson (Mark T. Nelson)
 * added Kinetic energy
 * 
 * Revision 1.4  94/09/30  09:10:07  09:10:07  nelson (Mark T. Nelson)
 * Enhanced output format
 * 
 * Revision 1.3  94/09/29  16:35:20  16:35:20  nelson (Mark T. Nelson)
 * Added total energy to output
 * 
 * Revision 1.2  94/09/29  12:32:11  12:32:11  nelson (Mark T. Nelson)
 * Changed format of energy output
 * 
 * Revision 1.1  94/09/29  12:04:57  12:04:57  gursoy (Attila Gursoy)
 * Initial revision
 * 
 ***************************************************************************/
