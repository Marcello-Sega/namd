/***************************************************************************/
/*      (C) Copyright 1995,1996,1997 The Board of Trustees of the          */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 * This object outputs the data collected on the master node
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Output.C,v 1.9 1997/09/18 22:22:21 jim Exp $";

#include <string.h>
#include <stdlib.h>

#include "Output.h"
#include "dcdlib.h"
#include "strlib.h"
//#include "vmd_tags.h"
#include "Inform.h"
#include "Molecule.h"
#include "Node.h"
#include "Parameters.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Vector.h"
#include "structures.h"
#include "Message.h"
#include "Communicate.h"

// These make the NAMD 1 names work in NAMD 2
#define namdMyNode Node::Object()
#define simParams simParameters
#define pdbData pdb

#ifdef OLD_MDCOMM
//  The function vmd_send_static() and vmd_send_dynamic are basically just
//  wrapper functions for the member function send_vmd_static() and
//  send_vmd_dyn().  But since the mdcomm routines want to be passed function
//  pointers, we need these functions since C++ won't allow the passing
//  of member functions as pointers
int vmd_send_static(mdc_app_arena *);
int vmd_send_dyn(mdc_app_arena *);
#endif

/************************************************************************/
/*									*/
/*			FUNCTION Output					*/
/*									*/
/*	This is the constructor for the Ouput class.  It just sets 	*/
/*  up some values for the VMD connection.				*/
/*									*/
/************************************************************************/

Output::Output()
{
        // initialize flags for warning messages
        teTempWarning = TRUE;

#ifdef MDCOMM
	//  Allocate a structure for the VMD connection.  This will 
	//  only be used if there actually is a VMD connection, but
	//  it doesn't occupy much space so we'll just allocate it here
	// vmdData = new VmdDynData;
	vmdData = (VmdDynData *) calloc(1,sizeof(VmdDynData));
	if ( vmdData == NULL )
	{
	  NAMD_die("memory allocation failed in Output::Output");
	}
	vmdData->timestep     = (int *)   calloc(1, sizeof(int));
	vmdData->elapsed_time = (float *) calloc(1, sizeof(float));
	vmdData->elapsed_cpu  = (float *) calloc(1, sizeof(float));
	vmdData->T            = (float *) calloc(1, sizeof(float));
	vmdData->Epot         = (float *) calloc(1, sizeof(float));
	vmdData->Etot         = (float *) calloc(1, sizeof(float));
	vmdData->Evdw         = (float *) calloc(1, sizeof(float));
	vmdData->Eelec        = (float *) calloc(1, sizeof(float));
	vmdData->Esup         = (float *) calloc(1, sizeof(float));
	vmdData->Ehbo         = (float *) calloc(1, sizeof(float));
	vmdData->Ebond        = (float *) calloc(1, sizeof(float));
	vmdData->Eangle       = (float *) calloc(1, sizeof(float));
	vmdData->Edihe        = (float *) calloc(1, sizeof(float));
	vmdData->Eimpr        = (float *) calloc(1, sizeof(float));
	vmdData->numPatches   = (int *)   calloc(1, sizeof(int));

	vmdStaticData = (VmdStaticData *) malloc(sizeof(VmdStaticData));
	if (vmdStaticData == NULL) 
	{
	  NAMD_die("memory allocation failed in Output::Output");
	}
	vmdHandle = NULL;
#else
	/*  Not needed, because they are not declared
	vmdData = NULL;
	vmdStaticData = NUL;
	vmdHandle = NULL;
	*/
#endif
}
/*			END OF FUNCTION Output				*/

/************************************************************************/
/*									*/
/*			FUNCTION ~Output				*/
/*									*/
/*	This is the destructor for the Output class, just free the	*/
/*   VMD data if there is any . . .					*/
/*									*/
/************************************************************************/

Output::~Output()
{
  //  Free the vmdData structure	
#ifdef MDCOMM
  if (vmdData != NULL)
  {
    free(vmdData);
  }
  //  Also the vmdStaticData static data structure;
  if (vmdStaticData != NULL)
  {
    free(vmdStaticData);    
  }
#endif

}
/*			END OF FUNCTION ~Output				*/

#if 0

/************************************************************************/
/*									*/
/*			FUNCTION energy					*/
/*									*/
/*   INPUTS:								*/
/*	timestep - timestep the energies are for			*/
/*	energy - Array containing individual energy components		*/
/*									*/
/*	This function recieves the accumulated energies from the	*/
/*   Collect object and outputs them.  It also sums them all and reports*/
/*   the total energy and also computes and reports the Temperature.	*/
/*									*/
/************************************************************************/

void Output::energy(int timestep, BigReal *energy)
{
   BigReal totalEnergy=0.0;
   int i;
   char out_string[512];
   char tmp_string[21];
   BigReal temperature;
   BigReal momentum;
   static int counter = 0;	//  Number of output lines printed

   for (i=0; i<maxEnergyCount; i++)
   {
     if (i != ENERGY_HBIN)
       totalEnergy += energy[i];
   }

   // temperature = 4.0 * energy[ENERGY_KINE] /
   //  (6.0 * namdMyNode->structure->numAtoms * BOLTZMAN);
 
   // temperature depends on teh number of degrees of freedon.
   // calculation is as follows

   temperature =  2.0*energy[ENERGY_KINE] / 
                  (namdMyNode->numDegFreedom*BOLTZMAN);

   momentum = sqrt(energy[MOMENTUM_X]*energy[MOMENTUM_X] +
		   energy[MOMENTUM_Y]*energy[MOMENTUM_Y] +
		   energy[MOMENTUM_Z]*energy[MOMENTUM_Z]);

   if ( (timestep%namdMyNode->simParams->outputEnergies) == 0)
   {
      if ( (counter % 10) == 0)
      {
	   namdInfo << "ETITLE:     TS    BOND        ANGLE       DIHED       IMPRP       ELECT       VDW       ";

	   if (namdMyNode->simParams->constraintsOn 
	       || namdMyNode->simParams->sphericalBCOn
	       || namdMyNode->simParams->cylindricalBCOn)
	   {
	   	namdInfo << "CONS      ";
	   }

	   if (namdMyNode->simParams->eFieldOn)
	   {
	   	namdInfo << "EFIELD    ";
	   }

	   if (namdMyNode->simParams->HydrogenBonds)
	   {
	   	namdInfo << "HBONDS     ";
	   	namdInfo << "Num HB     ";
	   }

	   namdInfo << "KINETIC       TOTAL      TEMP       |MOMENTUM|\n";
      }

      sprintf(out_string, "ENERGY: %6d ",timestep);

      for (i=0; i<maxEnergyCount; i++)
      {
	if (( (i==ENERGY_CONS) &&
	      !( (namdMyNode->simParams->constraintsOn) 
		|| (namdMyNode->simParams->sphericalBCOn)
		|| (namdMyNode->simParams->cylindricalBCOn)
	      )
	    ) 
	    || ( (i==ENERGY_EFLD) && !(namdMyNode->simParams->eFieldOn))
	    || ( (i==ENERGY_HBEN) && !(namdMyNode->simParams->HydrogenBonds))
	    || ( (i==ENERGY_HBIN) && !(namdMyNode->simParams->HydrogenBonds))
	    || (i==MOMENTUM_X) || (i==MOMENTUM_Y) || (i==MOMENTUM_Z) 
	   )
	{
	  continue;
	}
	sprintf(tmp_string, "%.4f ", energy[i]);
	NAMD_pad(tmp_string, 12);
	strcat(out_string, tmp_string);
      }

      sprintf(tmp_string, "%.4f ", totalEnergy);
      NAMD_pad(tmp_string, 11);
      strcat(out_string, tmp_string);

      sprintf(tmp_string, "%.4f", temperature);
      NAMD_pad(tmp_string, 11);
      strcat(out_string, tmp_string);
	
      sprintf(tmp_string, "%.4e", momentum);
      NAMD_pad(tmp_string, 11);
      strcat(out_string, tmp_string);
	
      namdInfo << out_string << sendmsg;

      // print warning message if temperature exceeds certain value
      if(temperature > OUTPUT_TEMPERATURE_WARNING && teTempWarning)
      {
	namdWarn << temperature << " degrees!!  ";
	namdWarn << "Man, tiny Elvis, it is hot in here!" << sendmsg;
	teTempWarning = FALSE;
      }

      counter++;
   }

   //  If temperature rescaling is active, deposit the current temperature
   //  for later use
   if (namdMyNode->simParams->rescaleFreq > 0)
   {
	namdMyNode->deposit_temperature(timestep, temperature);
   }

   //  If temperature coupling is active, broadcast the temperature
   //  to all the nodes
   if (namdMyNode->simParams->tCoupleOn)
   {
	Message *tempMsg;

	tempMsg = new Message;
	if ( tempMsg == NULL )
	{
	  NAMD_die("memory allocation failed in Output::energy");
	}

	tempMsg->put(timestep);
	tempMsg->put(temperature);

	comm->broadcast_all(tempMsg, LASTTEMPTAG);
   }

#ifdef MDCOMM
   //  If we have a VMD connection, and it is either the 0th timestep
   //  or a multiple of the vmdFrequency, send the energies to the
   //  VMD connection as well
   if ( (namdMyNode->simParams->vmdFrequency != -1) && 
	( ( ( timestep % namdMyNode->simParams->vmdFrequency) == 0) || 
	  (timestep==0) ) )
   {
	gather_vmd_energies(timestep, energy, temperature, totalEnergy);
   }
#endif
}
/*			END OF FUNCTION energy				*/

#endif

/************************************************************************/
/*									*/
/*			FUNCTION coordinate				*/
/*									*/
/*   INPUTS:								*/
/*	timestep - Timestep coordinates were accumulated for		*/
/*	n - This is the number of coordinates accumulated.		*/
/*	vel - Array of Vectors containing the velocities		*/
/*									*/
/*	This function receives the coordinates accumulated for a given  */
/*   timestep from the Collect object and calls the appropriate output  */
/*   functions.  ALL routines used to output coordinates information    */
/*   should be called from here.					*/
/*									*/
/************************************************************************/

int Output::coordinateNeeded(int timestep)
{
	SimParameters *simParams = Node::Object()->simParameters;

	int positionsNeeded = 0;

	//  Output a DCD trajectory 
	if ( (simParams->dcdFrequency != -1) &&
	     ((timestep % simParams->dcdFrequency) == 0) &&
	     ((timestep != simParams->firstTimestep)||(timestep == 0)) )
	{
		positionsNeeded = 1;
	}

	//  Output a restart file
	if ( (simParams->restartFrequency != -1) &&
	     ((timestep % simParams->restartFrequency) == 0) &&
	     (timestep != simParams->firstTimestep) )
	{
		positionsNeeded = 1;
	}

	//  Output final coordinates
	if (timestep == simParams->N)
	{
		positionsNeeded = 1;
	}

#ifdef MDCOMM
	//  If there is an active VMD connection and it is either the
	//  0th timestep or a mutiple of the vmdFrequency, then send
	//  the coordinates to VMD as well
        if ( (simParams->vmdFrequency != -1) && 
	   ( ( (timestep % simParams->vmdFrequency) == 0) ||
	     (timestep == simParams->firstTimestep) ) )
        {
		positionsNeeded = 1;
        }
#endif

	return positionsNeeded;
}

void Output::coordinate(int timestep, int n, Vector *coor)
{
	SimParameters *simParams = Node::Object()->simParameters;

	//  Output a DCD trajectory 
	if ( (simParams->dcdFrequency != -1) &&
	     ((timestep % simParams->dcdFrequency) == 0) &&
	     ((timestep != simParams->firstTimestep)||(timestep == 0)) )
	{
		output_dcdfile(timestep, n, coor);
	}

	//  Output a restart file
	if ( (simParams->restartFrequency != -1) &&
	     ((timestep % simParams->restartFrequency) == 0) &&
	     (timestep != simParams->firstTimestep) )
	{
		output_restart_coordinates(coor, n, timestep);
	}

	//  Output final coordinates
	if (timestep == simParams->N)
	{
		output_final_coordinates(coor, n, timestep);
	}

#ifdef MDCOMM
	//  If there is an active VMD connection and it is either the
	//  0th timestep or a mutiple of the vmdFrequency, then send
	//  the coordinates to VMD as well
        if ( (simParams->vmdFrequency != -1) && 
	   ( ( (timestep % simParams->vmdFrequency) == 0) ||
	     (timestep == simParams->firstTimestep) ) )
        {
		gather_vmd_coords(timestep, n, coor);
        }
#endif
}
/*		END OF FUNCTION coordinate				*/

/************************************************************************/
/*									*/
/*			FUNCTION velocity				*/
/*									*/
/*   INPUTS:								*/
/*	timestep - Timestep velocities were accumulated for		*/
/*	n - This is the number of velocities accumulated.		*/
/*	vel - Array of Vectors containing the velocities		*/
/*									*/
/*	This function receives the velocities accumulated for a given   */
/*   timestep from the Collect object and calls the appropriate output  */
/*   functions.  ALL routines used to output velocity information should*/
/*   be called from here.						*/
/*									*/
/************************************************************************/

int Output::velocityNeeded(int timestep)
{
	SimParameters *simParams = Node::Object()->simParameters;

	int velocitiesNeeded = 0;

	//  Output a velocity DCD trajectory
	if ( (simParams->velDcdFrequency != -1) &&
	     ((timestep % simParams->velDcdFrequency) == 0) &&
	     ((timestep != simParams->firstTimestep)||(timestep == 0)) )
	{
		velocitiesNeeded = 1;
	}

	//  Output a restart file
	if ( (simParams->restartFrequency != -1) &&
	     ((timestep % simParams->restartFrequency) == 0) &&
	     (timestep != simParams->firstTimestep) )
	{
		velocitiesNeeded = 1;
	}

	//  Output final velocities
	if (timestep == simParams->N)
	{
		velocitiesNeeded = 1;
	}

	return velocitiesNeeded;
}

void Output::velocity(int timestep, int n, Vector *vel)
{
	SimParameters *simParams = Node::Object()->simParameters;

	//  Output velocity DCD trajectory
	if ( (simParams->velDcdFrequency != -1) &&
	     ((timestep % simParams->velDcdFrequency) == 0) &&
	     ((timestep != simParams->firstTimestep)||(timestep == 0)) )
	{
		output_veldcdfile(timestep, n, vel);
	}

	//  Output restart file
	if ( (simParams->restartFrequency != -1) &&
	     ((timestep % simParams->restartFrequency) == 0) &&
	     (timestep != simParams->firstTimestep) )
	{
		output_restart_velocities(timestep, n, vel);
	}

	//  Output final velocities
	if (timestep == simParams->N)
	{
		output_final_velocities(timestep, n, vel);
	}
}
/*			END OF FUNCTION velocity			*/

/************************************************************************/
/*									*/
/*			FUNCTION output_restart_coordinates		*/
/*									*/
/*   INPUTS:								*/
/*	coor - Array of vectors containing current positions		*/
/*	n - Number of coordinates to output				*/
/*	timestep - Timestep for which the coordinates are being written */
/*									*/
/*	This function writes out the current positions of all the atoms */
/*   in PDB format to the restart file specified by the user in the 	*/
/*   configuration file.						*/
/*									*/
/************************************************************************/

void Output::output_restart_coordinates(Vector *coor, int n, int timestep)

{
	static Bool first=TRUE;		//  Is this the first time function
					//  has been called?
	static char *restart_name,      //  filenames for restart and backup
		    *restart_bak;	//  restart files
	char comment[128];		//  Comment for header of PDB file

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
/*			END OF FUNCTION output_restart_coordinates	*/

/************************************************************************/
/*									*/
/*			FUNCTION output_restart_velocities		*/
/*									*/
/*   INPUTS:								*/
/*	vel - Array of vectors containing current velocities		*/
/*	timestep - Timestep for which the velocities are being written  */
/*									*/
/*	This function writes out the current velocites of all the atoms */
/*   in PDB format to the restart file specified by the user in the 	*/
/*   configuration file.						*/
/*									*/
/************************************************************************/

void Output::output_restart_velocities(int timestep, int n, Vector *vel)

{
	static Bool first=TRUE;		//  is this the first invocation
					//  of the function?
	static char *restart_name,    	//  filenames for the restart and
		    *restart_bak;	//  backup files
	char comment[128];		//  comment for the header of the
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

		scale_vels(vel, n, 20);
		namdMyNode->pdbData->set_all_positions(vel);
		namdMyNode->pdbData->write(restart_name, comment);
		scale_vels(vel, n, 0.05);
	}
	else
	{
		//  Write the velocities to a binary file
		write_binary_file(restart_name, n, vel);
	}
}
/*			END OF FUNCTION output_restart_velocities	*/

/************************************************************************/
/*									*/
/*			FUNCTION output_dcdfile				*/
/*									*/
/*   INPUTS:								*/
/*	timestep - Current timestep					*/
/*	n - Number of atoms in simulation				*/
/*	coor - Coordinate vectors for all atoms				*/
/*									*/
/*	This function maintains the interface between the Output object */
/*   and the dcd writing routines contained in dcdlib.			*/
/*									*/
/************************************************************************/

void Output::output_dcdfile(int timestep, int n, Vector *coor)

{
	static Bool first=TRUE;	//  Flag indicating first call
	static int fileid;	//  File id for the dcd file
	static float *x, *y, *z; // Arrays to hold x, y, and z arrays
	int i;			//  Loop counter
	int ret_code;		//  Return code from DCD calls

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
		fileid=open_dcd_write(namdMyNode->simParams->dcdFilename);

		if (fileid == DCD_FILEEXISTS)
		{
			char err_msg[257];

			sprintf(err_msg, "DCD file %s already exists!!",
				namdMyNode->simParams->dcdFilename);

			NAMD_die(err_msg);
		}
		else if (fileid < 0)
		{
			char err_msg[257];

			sprintf(err_msg, "Couldn't open DCD file %s",
				namdMyNode->simParams->dcdFilename);

			NAMD_die(err_msg);
		}


		//  Write out the header
		ret_code = write_dcdheader(fileid, 
				namdMyNode->simParams->dcdFilename,
				n, 
				(int) ((namdMyNode->simParams->N+1)/namdMyNode->simParams->dcdFrequency)+1,
				0, 
				namdMyNode->simParams->dcdFrequency, 
				namdMyNode->simParams->dt/TIMEFACTOR);


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
	if ( (timestep+namdMyNode->simParams->dcdFrequency) >
	      namdMyNode->simParams->N)
	{
		close_dcd_write(fileid);
	}

}
/*			END OF FUNCTION output_dcdfile			*/

/************************************************************************/
/*									*/
/*			FUNCTION output_final_coordinates		*/
/*									*/
/*   INPUTS:								*/
/*	coor - Array of vectors containing final coordinates		*/
/*	n - Number of coordinates to output				*/
/*	timestep - Timestep that coordinates are being written in	*/
/*									*/
/*	This function writes out the final coordinates for the		*/
/*   simulation in PDB format to the file specified in the config	*/
/*   file.								*/
/*									*/
/************************************************************************/

void Output::output_final_coordinates(Vector *coor, int n, int timestep)

{
	static char output_name[140];	//  Output filename
	char comment[128];		//  comment for PDB header

	//  Built the output filename
	strcpy(output_name, namdMyNode->simParams->outputFilename);
	strcat(output_name, ".coor");

	//  Check to see if we should write out a binary file or a
	//  PDB file
	if (!namdMyNode->simParams->binaryRestart)
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
/*		END OF FUNCTION output_final_coordinates		*/

/************************************************************************/
/*									*/
/*			FUNCTION output_final_velocities		*/
/*									*/
/*   INPUTS:								*/
/*	vel - Array of vectors containing final velocities		*/
/*	timestep - Timestep that vleocities are being written in	*/
/*									*/
/*	This function writes out the final vleocities for the		*/
/*   simulation in PDB format to the file specified in the config	*/
/*   file.								*/
/*									*/
/************************************************************************/

void Output::output_final_velocities(int timestep, int n, Vector *vel)

{
	static char output_name[140];	//  Output filename
	char comment[128];		//  Comment for PDB header

	//  Build the output filename
	strcpy(output_name, namdMyNode->simParams->outputFilename);
	strcat(output_name, ".vel");

	//  Check to see if we should write a PDB or binary file
	if (!(namdMyNode->simParams->binaryRestart))
	{
		//  Write the final velocities to a PDB file
		sprintf(comment, "FINAL VELOCITIES WRITTEN BY NAMD AT TIMESTEP %d", timestep);

		scale_vels(vel, n, 20);
		namdMyNode->pdbData->set_all_positions(vel);
		namdMyNode->pdbData->write(output_name, comment);
		scale_vels(vel, n, 0.05);
	}
	else
	{
		//  Write the coordinates to a binary file
		write_binary_file(output_name, n, vel);
	}

}
/*			END OF FUNCTION output_final_velocities		*/

/************************************************************************/
/*									*/
/*			FUNCTION output_veldcdfile			*/
/*									*/
/*   INPUTS:								*/
/*	timestep - Current timestep					*/
/*	n - Number of atoms in simulation				*/
/*	coor - velocity vectors for all atoms				*/
/*									*/
/*	This function maintains the interface between the Output object */
/*   and the dcd writing routines contained in dcdlib.  This fucntion   */
/*   writes out the velocity vectors in DCD format.			*/
/*									*/
/************************************************************************/

void Output::output_veldcdfile(int timestep, int n, Vector *vel)

{
	static Bool first=TRUE;	//  Flag indicating first call
	static int fileid;	//  File id for the dcd file
	static float *x, *y, *z; // Arrays to hold x, y, and z arrays
	int i;			//  Loop counter
	int ret_code;		//  Return code from DCD calls

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


		//  Write out the header
		ret_code = write_dcdheader(fileid, namdMyNode->simParams->velDcdFilename, n, 
				(int) ((namdMyNode->simParams->N+1)/namdMyNode->simParams->velDcdFrequency)+1,
				0, namdMyNode->simParams->velDcdFrequency, 
				namdMyNode->simParams->dt/TIMEFACTOR);


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
/*			END OF FUNCTION output_veldcdfile		*/

/************************************************************************/
/*									*/
/*			FUNCTION write_binary_file			*/
/*									*/
/*   INPUTS:								*/
/*	fname - file name to write velocities to			*/
/*	n - Number of atoms in system					*/
/*	vels - Array of vectors						*/
/*									*/
/*	This function writes out vectors in binary format to		*/
/*   the specified file.						*/
/*									*/
/************************************************************************/

void Output::write_binary_file(char *fname, int n, Vector *vecs)

{
	FILE *fp;		//  File descriptor
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

	fclose(fp);
}
/*			END OF FUNCTION write_binary_file		*/

/************************************************************************/
/*									*/
/*			FUNCTION scale_vels				*/
/*									*/
/*   INPUTS:								*/
/*	v - Array of velocity vectors					*/
/*	n - Number of atoms in system					*/
/*	fact - Scaling factor						*/
/*									*/
/*	This function scales all the vectors passed in by a constant	*/
/*   factor.  This is used before writing out velocity vectors that	*/
/*   need to be resized.						*/
/*									*/
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
/*			END OF FUNCTION scale_vels			*/

#ifdef MDCOMM

/************************************************************************/
/*									*/
/*			FUNCTION build_vmd_static_data			*/
/*									*/
/*	This function builds and returns a structure that contains      */
/*  all the static data that is send to VMD using the MDComm library    */
/*  routines.								*/
/*									*/
/************************************************************************/

VmdStaticData *Output::build_vmd_static_data(VmdStaticData *sdata)
{
	int i, j, k;			//  Loop counters
	Bond *bondPtr;			//  Pointer to bonded structure
	PDB *pdbData;			//  PDB object ot use
	int strLen;			//  Length of the current string
	char strBuf[10];		//  Temporary string buffer
	Real sigma, epsilon, s14, e14;  //  Vdw wall parameters
	BigReal two_6;			//  2^(1/6)

	//  Get the pdb object to use
	pdbData=namdMyNode->pdbData;

	//  Calculate 2^(1/6) to use in vdw radii calculation
	two_6 = pow(2.0, 0.166666666667);

	//  Allocate a static data structure
	if (sdata == NULL)
        {
	  sdata = (VmdStaticData *) malloc(sizeof(VmdStaticData));
	}

	if (sdata == NULL)
	{
	  NAMD_die("memory allocation failed in Output::build_vmd_static_data()");
	}

	//  Set the atom name length
	sdata->nameLen = VMD_NAME_LEN;

	//  Get the number of atoms
	sdata->numAtoms = namdMyNode->structure->numAtoms;

	//  Allocate data for structures that require it
#ifdef OLD_MDCOMM
	sdata->atomNames =  new char[sdata->numAtoms*VMD_NAME_LEN]; 
	sdata->atomNameIndexes =  new int[sdata->numAtoms]; 
	sdata->atomTypes =  new char[sdata->numAtoms*VMD_NAME_LEN]; 
	sdata->atomTypeIndexes =  new int[sdata->numAtoms]; 
	sdata->bonds = new int[namdMyNode->structure->numBonds*2];
	sdata->resIds = new char[sdata->numAtoms*VMD_NAME_LEN];
	sdata->resNames = new char[sdata->numAtoms*VMD_NAME_LEN];
	sdata->radii =  new float[sdata->numAtoms]; 
	sdata->charge =  new float[sdata->numAtoms]; 
	sdata->mass =  new float[sdata->numAtoms]; 
	sdata->occupancy =  new float[sdata->numAtoms]; 
	sdata->beta =  new float[sdata->numAtoms]; 
	sdata->segIds = new char[sdata->numAtoms*VMD_NAME_LEN];
#endif
	sdata->atomNames =  (char *) calloc(sdata->numAtoms*VMD_NAME_LEN, 
					    sizeof(char)); 
	sdata->atomNameIndexes =  (int *) calloc(sdata->numAtoms, sizeof(int));
	sdata->atomTypes =  (char *) calloc(sdata->numAtoms*VMD_NAME_LEN, 
					    sizeof(char)); 
	sdata->atomTypeIndexes =  (int *) calloc(sdata->numAtoms, 
						 sizeof(int)); 
	sdata->bonds = (int *) calloc(namdMyNode->structure->numBonds*2, 
				      sizeof(int));
	sdata->resIds = (char *) calloc(sdata->numAtoms*VMD_NAME_LEN, 
					sizeof(char));
	sdata->resNames = (char *) calloc(sdata->numAtoms*VMD_NAME_LEN, 
					  sizeof(char));
	sdata->radii =  (float *) calloc(sdata->numAtoms, sizeof(float)); 
	sdata->charge =  (float *) calloc(sdata->numAtoms, sizeof(float)); 
	sdata->mass =  (float *) calloc(sdata->numAtoms, sizeof(float)); 
	sdata->occupancy =  (float *) calloc(sdata->numAtoms, sizeof(float)); 
	sdata->beta =  (float *) calloc(sdata->numAtoms, sizeof(float)); 
	sdata->segIds = (char *) calloc(sdata->numAtoms*VMD_NAME_LEN, 
					sizeof(char));


	if ( (sdata->atomNames == NULL) ||
	     (sdata->atomNameIndexes == NULL) ||
	     (sdata->atomTypes == NULL) ||
	     (sdata->atomTypeIndexes == NULL) ||
	     (sdata->bonds == NULL) ||
	     (sdata->resIds == NULL) ||
	     (sdata->resNames == NULL) ||
	     (sdata->radii == NULL) ||
	     (sdata->charge == NULL) ||
	     (sdata->mass == NULL) ||
	     (sdata->occupancy == NULL) ||
	     (sdata->beta == NULL) ||
	     (sdata->segIds == NULL) )
	{
		NAMD_die("memory allocation failed in Output::build_vmd_static_data()");
	}

	//  Initialize unknown quantities
	sdata->numAtomNames=0;
	sdata->numAtomTypes=0;

	//  Get the patch info
	sdata->maxNumPatches = namdMyNode->patchMap->maxPatchNum;

	//  Get the bond info
	sdata->numBonds = namdMyNode->structure->numBonds;

	for (i=0, j=0; i<sdata->numBonds; i++, j+=2)
	{
		bondPtr = namdMyNode->structure->get_bond(i);

		sdata->bonds[j] = bondPtr->atom[0];
		sdata->bonds[j+1] = bondPtr->atom[1];
	}
	//  Set the other (unused) static data members
	sdata->nHBondDonors = 0;
	sdata->donor = (int *) NULL;
	sdata->donorh = (int *) NULL;
	sdata->nHBondAcceptors = 0;
	sdata->acceptor = (int *) NULL;
	sdata->acceptora = (int *) NULL;
       

	//  Make a pass through all the atoms.  During this pass, we will
	//  get atom information that we can directly and build the list
	//  of unique atom names and types
	for (i=0, j=0; i<sdata->numAtoms; i++, j+=VMD_NAME_LEN)
	{
		sdata->mass[i] = namdMyNode->structure->atommass(i);
		sdata->charge[i] = namdMyNode->structure->atomcharge(i);
		sdata->occupancy[i] = (pdbData->atom(i))->occupancy();
		sdata->beta[i] = (pdbData->atom(i))->temperaturefactor();

		//  Copy the residue name
		strcpy(strBuf, (pdbData->atom(i))->residuename());
		strLen = strlen(strBuf);

		for (k=0; k<VMD_NAME_LEN; k++)
		{
			if (k<strLen)
			{
				sdata->resNames[j+k] = strBuf[k];
			}
			else
			{
				sdata->resNames[j+k] = ' ';
			}
		}

		//  Copy the residue id
		sprintf(strBuf, "%d", (pdbData->atom(i))->residueseq());
		strLen = strlen(strBuf);

		for (k=0; k<VMD_NAME_LEN; k++)
		{
			if (k<strLen)
			{
				sdata->resIds[j+k] = strBuf[k];
			}
			else
			{
				sdata->resIds[j+k] = ' ';
			}
		}

		//  Copy the Segment id
		strcpy(strBuf,(pdbData->atom(i))->segmentname());
		strLen = strlen(strBuf);

		for (k=0; k<VMD_NAME_LEN; k++)
		{
			if (k<strLen)
			{
				sdata->segIds[j+k] = strBuf[k];
			}
			else
			{
				sdata->segIds[j+k] = ' ';
			}
		}

		//  Calculate the vdw radii
		namdMyNode->params->get_vdw_params(&sigma, &epsilon, &s14, &e14,
						   namdMyNode->structure->atomvdwtype(i));

		sdata->radii[i] = sigma*two_6*0.5;  //  Note that the extra factor
						    //  of 0.5 is purely a fudge
						    //  factor to make VMD look
						    //  better!!

		//  Insert the atom name list
		insert_into_unique_list((pdbData->atom(i))->name(), 
					sdata->atomNames,
					sdata->numAtomNames);

		//  Insert into atom type list
		insert_into_unique_list(namdMyNode->structure->get_atomtype(i), 
					sdata->atomTypes, sdata->numAtomTypes);
		
	}

	//  Now make another pass and get the atom name and type indexes 
	//  for each atom
	for (i=0; i<sdata->numAtoms; i++)
	{
		if ( (sdata->atomNameIndexes[i]=search_list((pdbData->atom(i))->name(), 
							    sdata->atomNames, 
							    sdata->numAtomNames)) == -1)
		{
			NAMD_die("Couldn't find atom name in unique list");
		}

		if ((sdata->atomTypeIndexes[i]=search_list(namdMyNode->structure->get_atomtype(i), 
						sdata->atomTypes, sdata->numAtomTypes)) == -1)
		{
			NAMD_die("Couldn't find atom type in unique list");
		}
	}

	return(sdata);
}
/*			END OF FUNCTION build_vmd_static_data		*/

/************************************************************************/
/*									*/
/*			FUNCTION insert_into_unique_list		*/
/*									*/
/*   INPUTS:								*/
/*	item - String to insert into list				*/
/*	list - List to insert string into				*/
/*	size - Current number of entries in list			*/
/*									*/
/*	This function inserts a value into the atom name or atom type   */
/*   list.  This is kind of an ugly operation, since these lists are    */
/*   maintained as one long string, where each VMD_NAME_LEN characters  */
/*   is a different name.  Since these lists should contain only unique */
/*   values, if the value already exists in the list, then nothing is   */
/*   done.								*/
/*									*/
/************************************************************************/

void Output::insert_into_unique_list(const char *item, char *list, int &size)

{
     int position;		//  Position to insert value at			
     int i;			//  Loop counter
     int strLen;		//  Length of the string

     if ( search_list(item, list, size) == -1 )
     {
	//  The item doesn't already exist in the list, so insert it

	//  First, find where it goes in the list
	position = 0;

	while ( (position < size) && 
		(name_cmp(list+position*VMD_NAME_LEN, item) < 0) )
	{
		position++;
	}

	//  Next, if its not going at the end of the list, shift data to the
	//  right
	if (position != size)
	{
		for(i=(size*VMD_NAME_LEN)-1; i>=(position*VMD_NAME_LEN); i--)
		{
			list[i+VMD_NAME_LEN] = list[i];
		}
	}

	strLen = strlen(item);

	//  Now insert the new string in the right place
	for (i=0; i<VMD_NAME_LEN; i++)
	{
		if (i<strLen)
		{
			list[position*VMD_NAME_LEN + i] = item[i];
		}
		else
		{
			list[position*VMD_NAME_LEN + i] = ' ';
		}
	}

	size++;
    }
}
/*			END OF FUNCTION insert_into_unique_list		*/

/************************************************************************/
/*									*/
/*			FUNCTION namd_cmp				*/
/*									*/
/*   INPUTS:								*/
/*	str1 - 1st string to compare					*/
/*	str2 - 2nd string to compare					*/
/*									*/
/*	This fucntion is just a specialized string compare.  Because    */
/*   of the structure of the unique lists that we are taking, the	*/
/*   strings passed may be null terminated or just pointers into	*/
/*   the list itself.  So this function deals with either case and	*/
/*   returns the same comparison codes as strcmp.			*/
/*									*/
/************************************************************************/

int Output::name_cmp(const char str1[], const char str2[])

{
	char cp1[10], cp2[10];	  	//  Copies of strings 1 and 2
	int i;				//  Loop counter

	//  First, make copies of strings one and two.  If either is
	//  currently not null terminated, then make the copies so
	for (i=0; ( (i<VMD_NAME_LEN) && (str1[i]!='\0') ); i++)
	{
		cp1[i] = str1[i];
	}

	cp1[i]='\0';

	for (i=0; ( (i<VMD_NAME_LEN) && (str2[i]!='\0') ); i++)
	{
		cp2[i] = str2[i];
	}

	cp2[i]='\0';

	//  Pad them with spaces to be the full length
	NAMD_pad(cp1, VMD_NAME_LEN);
	NAMD_pad(cp2, VMD_NAME_LEN);

	//  Now we can just compare them, since they are now NULL terminated
	//  and space padded to be the same format
	return(strcmp(cp1, cp2));
}
/*			END OF FUNCTION name_cmp			*/

/************************************************************************/
/*									*/
/*			FUNCTION search_list				*/
/*									*/
/*   INPUTS:								*/
/*	item - String to look for in list				*/
/*	list - Unique list to search through				*/
/*	size - Size of the list to search				*/
/*									*/
/*	This function does a binary search through the sorted		*/
/*   unique list of atom names or atom types for the string indicated   */
/*   by item.								*/
/*	If the item is found, it position in the list is returned.	*/
/*   otherwise, a -1 is returned.					*/
/*									*/
/************************************************************************/


int Output::search_list(const char *item, const char *list, int size)

{
	int low, high;		//  Low and high search boundaries
	int ret_code;		//  Return code from namd_cmp
	int pos;		//  Current search position

	if (size == 0)
	{
		//  If the list is empty, return item not found
		return(-1);
	}

	low=0;
	high=size-1;

	pos=size/2;

	//  Keep searching until we find it, or the search is over
	while ( ((ret_code = name_cmp(item, (list+pos*VMD_NAME_LEN))) != 0) &&
		(low != high) && (low<high) )
	{
		if (ret_code < 0)
		{
			high = pos-1;
		}
		else
		{
			low = pos+1;
		}

		pos = (low+high)/2;
	}

	if (ret_code == 0)
	{
		return(pos);
	}
	else
	{
		return(-1);
	}
}
/*			END OF FUNCTION search_list			*/

/************************************************************************/
/*									*/
/*			FUNCTION free_vmd_static_data			*/
/*									*/
/*   INPUTS:								*/
/*	ptr - Pointer to the structure to be freed			*/
/*									*/
/*	This function frees a structure containing the static data	*/
/*   for VMD along with all of its components.				*/
/*									*/
/************************************************************************/

void Output::free_vmd_static_data(VmdStaticData *ptr)

{
 /*   NOTE (RK) - actually, the whole thing doesn't get freed because    */
 /*   we will need info from the static data structure (specifically,    */
 /*   the number of atoms) for use in each call to rapp_update_client.   */
 /*   For now, just free the memory-hogging pieces, and it can be done   */
 /*   more cleanly later.                                                */

       free(ptr->atomNames);
       free(ptr->atomNameIndexes);
       free(ptr->atomTypes);
       free(ptr->atomTypeIndexes);
       free(ptr->bonds);
       free(ptr->resIds);
       free(ptr->resNames);
       free(ptr->radii);
       free(ptr->charge);
       free(ptr->mass);
       free(ptr->occupancy);
       free(ptr->beta);
       free(ptr->segIds);
       /* free(ptr);*/

}
/*			END OF FUNCTION free_vmd_static			*/

/************************************************************************/
/*									*/
/*			FUNCTION gather_vmd_energies			*/
/*									*/
/*   INPUTS:								*/
/*	timestep - Timestep that the energies are for			*/
/*	energies - Array containing individual energy components	*/
/*	T - Temperature							*/
/*	totalEnergy - Epot + Ekin					*/
/*									*/
/*	This function gathers the energy information that is sent	*/
/*   to VMD as part of the dynamic information.  The energy information */
/*   is palced into the structure pointed to by vmdData, and if 	*/
/*   the corrseponding coordinate information has also arrived, then    */
/*   the data is sent							*/
/*									*/
/************************************************************************/

void Output::gather_vmd_energies(int timestep, BigReal *energies, BigReal T, 
				 BigReal totalEnergy)

{
	if (vmdData == NULL)
	{
		NAMD_die("Gathered VMD energy data without allocating structure!");
	}

	//  Check the timestep
	if (vmdData->coordsArrived && (*vmdData->timestep != timestep) )
	{
		//  The coordinates have already arrived, and the timestep
		//  doesn't match our timestep.  Something is hosed!
		NAMD_die("Mismatch between timestep on VMD energies and coordinates");
	}
	else
	{
		//  Set the timestep
		*vmdData->timestep = timestep;
	}

	//  Populate the approiate data values
	*vmdData->T = T;

	*vmdData->elapsed_time = (float) timestep*namdMyNode->simParams->dt*0.001;
	*vmdData->elapsed_cpu = 0.0;	//  Currently not implemented

	*vmdData->Epot = totalEnergy-energies[7];  //  Epot = Etot - Ekin
	*vmdData->Etot = totalEnergy;
	*vmdData->Evdw=energies[5];
	*vmdData->Eelec=energies[4];
	*vmdData->Ebond=energies[0];
	*vmdData->Eangle=energies[1];
	*vmdData->Edihe=energies[2];
	*vmdData->Eimpr=energies[3];

	if (vmdData->coordsArrived)
	{
		//  Coordinates are already here, so send the data
		if (timestep == 0)
		{
			//  If this is the first collection, initialize
			//  the connection.  During the initialization,
			//  the dynamic data is sent
			initialize_vmd_connection();
		}
		else
		{
			//  This is a normal update, so just call
			//  the update function
			if (timestep+namdMyNode->simParams->vmdFrequency<namdMyNode->simParams->N)
			{
				//  This is NOT the last update
				rapp_app_update_client(vmdHandle, 1);
//				vmd_process_events();
			}
			else
			{
				rapp_app_update_client(vmdHandle, 0);
				close_vmd_connection();
			}
		}

		//  Reset the arrived flags
		vmdData->energiesArrived = FALSE;
		vmdData->coordsArrived = FALSE;
	}
	else
	{
		//  Coordinates aren't here yet, so just set the flag
		//  The data will be sent when the coordinates arrive
		vmdData->energiesArrived = TRUE;
	}
}
/*			END OF FUNCTION gather_vmd_energies		*/

/************************************************************************/
/*									*/
/*			FUNCTION gather_vmd_coords			*/
/*									*/
/*   INPUTS:								*/
/*	timestep - Timestep that the coordinates are for		*/
/*	N - Number of coordinates					*/
/*	coords - Array of Vectors containing positions			*/
/*									*/
/*	This function gathers the current positions of all the atoms	*/
/*   in the simulation to send to VMD.					*/
/*									*/
/************************************************************************/

void Output::gather_vmd_coords(int timestep, int N, Vector *coords)

{
	static Bool first=TRUE;		//  Flag TRUE -> first call to function
	int i;				//  Loop counter

	if (vmdData == NULL)
	{
		NAMD_die("Tried to collect namd coords without allocating structure");
	}

	//  Check the timestep
	if (vmdData->energiesArrived && (*vmdData->timestep != timestep) )
	{
		//  The energies have already been received, and the timestep
		//  from the energy values doesn't match our value.  Something
		//  is hosed.
		NAMD_die("Mismatch between timestep on VMD energies and coordinates");
	}
	else
	{
		//  Set the timestep value
		*vmdData->timestep = timestep;
	}

	if (first)
	{
		Vector pOrigin;		//  Origin of the patch

		//  This is the first call to the function, so allocate
		//  arrays to hold the X, Y, and Z coordinates
		vmdData->X = (float *) calloc(N, sizeof(float));
		vmdData->Y = (float *) calloc(N, sizeof(float));
		vmdData->Z = (float *) calloc(N, sizeof(float));

		if ( (vmdData->X == NULL) || (vmdData->Y == NULL) ||
		     (vmdData->Z == NULL) )
		{
			NAMD_die("malloc failed in Output::gather_vmd_coords");
		}

		vmdStaticData->numAtoms = N;

		//  Now get the patch information
		*vmdData->numPatches = namdMyNode->patchMap->maxPatchNum;

		//  Allocate arrays to hold all the data
#ifdef OLD_MDCOMM
		vmdData->pXOrigins = new float[vmdData->numPatches];
		vmdData->pYOrigins = new float[vmdData->numPatches];
		vmdData->pZOrigins = new float[vmdData->numPatches];
		vmdData->patchLength = new float[vmdData->numPatches];
		vmdData->patchWidth = new float[vmdData->numPatches];
		vmdData->patchHeight = new float[vmdData->numPatches];
		vmdData->patchAtomNums = new float[vmdData->numPatches];
		vmdData->patchLoads = new float[vmdData->numPatches];
		vmdData->patchNode = new float[vmdData->numPatches];
#endif

		vmdData->pXOrigins = (float *) calloc(*vmdData->numPatches, 
						      sizeof(float));
		vmdData->pYOrigins = (float *) calloc(*vmdData->numPatches, 
						      sizeof(float));
		vmdData->pZOrigins = (float *) calloc(*vmdData->numPatches, 
						      sizeof(float));
		vmdData->patchLength = (float *) calloc(*vmdData->numPatches, 
							sizeof(float));
		vmdData->patchWidth = (float *) calloc(*vmdData->numPatches, 
						       sizeof(float));
		vmdData->patchHeight = (float *) calloc(*vmdData->numPatches, 
							sizeof(float));
		vmdData->patchAtomNums = (float *) calloc(*vmdData->numPatches,
							  sizeof(float));
		vmdData->patchLoads = (float *) calloc(*vmdData->numPatches, 
						       sizeof(float));
		vmdData->patchNode = (float *) calloc(*vmdData->numPatches, 
						      sizeof(float));
		
		if ( (vmdData->pXOrigins == NULL) ||
		     (vmdData->pYOrigins == NULL) ||
		     (vmdData->pZOrigins == NULL) ||
		     (vmdData->patchLength == NULL) ||
		     (vmdData->patchWidth == NULL) ||
		     (vmdData->patchHeight == NULL) ||
		     (vmdData->patchAtomNums == NULL) ||
		     (vmdData->patchLoads == NULL) ||
		     (vmdData->patchNode == NULL) )
		{
			NAMD_die("memory allocation failed in Output::gather_vmd_coords");
		}

		//  Now populate the patch data
		for (i=0; i<*vmdData->numPatches; i++)
		{
			namdMyNode->patchMap->get_patch_origin(i, pOrigin);

			vmdData->pXOrigins[i] = pOrigin.x;
			vmdData->pYOrigins[i] = pOrigin.y;
			vmdData->pZOrigins[i] = pOrigin.z;

			vmdData->patchWidth[i] = namdMyNode->simParams->patchDimension;
			vmdData->patchHeight[i] = namdMyNode->simParams->patchDimension;
			vmdData->patchLength[i] = namdMyNode->simParams->patchDimension;

			vmdData->patchLoads[i] = 1.0;
			vmdData->patchNode[i] = namdMyNode->patchMap->patch_node(i);
			vmdData->patchAtomNums[i] = namdMyNode->patchMap->get_num_atoms(i);
		}

		first=FALSE;
	}

	//  Copy the current coordinates to the vmdData structure
	for (i=0; i<N; i++)
	{
		vmdData->X[i] = coords[i].x;
		vmdData->Y[i] = coords[i].y;
		vmdData->Z[i] = coords[i].z;
	}

	if (vmdData->energiesArrived)
	{
		//  The energies have already arrived, send
		//  the data now
		if (timestep == 0)
		{
			//  This is the first timestep, so initialize
			//  the connection.  During this initialization
			//  the dynamic data is sent for the first time
		        //  Actually, this is no longer true (RJK)
			initialize_vmd_connection();
		}
		//  This is a normal timestep, so just call
		//  the update function
		if (timestep+namdMyNode->simParams->vmdFrequency <
		    namdMyNode->simParams->N)
		{
		  //  This is NOT the last update
		  rapp_app_update_client(vmdHandle, 1);
		}
		else
		{
		  rapp_app_update_client(vmdHandle, 0);
		  close_vmd_connection();
		}

		//  Reset the arrived flags
		vmdData->energiesArrived = FALSE;
		vmdData->coordsArrived = FALSE;
	}
	else
	{
	   //  Energies haven't arrived, so just set the flag and
	   //  go on.  The data will be sent when the coordinates arrive
	   vmdData->coordsArrived=TRUE;
	}
}
/*			END OF FUNCTION gather_vmd_coords		*/

/************************************************************************/
/*									*/
/*			FUNCTION print_vmd_static_data			*/
/*									*/
/*	This is a debugging function that dumps the static data gathered*/
/*   for VMD to the screen.  It is never called during normal runs	*/
/*   and is provided only for testing purposes.				*/
/*									*/
/************************************************************************/

void Output::print_vmd_static_data()
{
	VmdStaticData *data = NULL;
	int i, j;
	char buf[10];

	data = build_vmd_static_data(data);

	namdInfo << "VMD STATIC DATA\n";
	namdInfo << "NumAtoms          " << data->numAtoms << "\n";
	namdInfo << "NumAtomNames      " << data->numAtomNames << "\n";
	namdInfo << "NumAtomTypes      " << data->numAtomTypes << "\n";
	namdInfo << "NumBonds          " << data->numBonds << "\n";
	namdInfo << "MaxNumPatches     " << data->maxNumPatches << sendmsg;

	namdInfo << "ATOM NAMES:\n";

	for (i=0; i<data->numAtomNames; i++)
	{
		for (j=0; j<VMD_NAME_LEN; j++)
		{
			buf[j] = data->atomNames[i*VMD_NAME_LEN+j];
		}

		buf[j] = '\0';

		namdInfo << buf << sendmsg;
	}
		
	namdInfo << "ATOM TYPES:\n";

	for (i=0; i<data->numAtomTypes; i++)
	{
		for (j=0; j<VMD_NAME_LEN; j++)
		{
			buf[j] = data->atomTypes[i*VMD_NAME_LEN+j];
		}

		buf[j] = '\0';

		namdInfo << buf << sendmsg;
	}

	namdInfo << "BONDS\n";

	for (i=0; i<data->numBonds; i++)
	{
		namdInfo << data->bonds[2*i] << " bonded to " << data->bonds[2*i+1] << sendmsg;
	}

	namdInfo << "ATOM INFORMATION\n";

	for (i=0; i<data->numAtoms; i++)
	{
		namdInfo << "ATOM NUMBER " << i << "\n";

		for (j=0; j<VMD_NAME_LEN; j++)
		{
			buf[j] = data->resIds[i*VMD_NAME_LEN+j];
		}

		buf[j] = '\0';

		namdInfo << "RES ID:        " << buf << "\n";

		for (j=0; j<VMD_NAME_LEN; j++)
		{
			buf[j] = data->resNames[i*VMD_NAME_LEN+j];
		}

		buf[j] = '\0';

		namdInfo << "RES NAME:      " << buf << "\n";
		
		for (j=0; j<VMD_NAME_LEN; j++)
		{
			buf[j] = data->segIds[i*VMD_NAME_LEN+j];
		}

		buf[j] = '\0';

		namdInfo << "SEG ID:        " << buf << "\n";

		namdInfo << "RADII:         " << data->radii[i] << "\n";
		namdInfo << "CHARGE:        " << data->charge[i] << "\n";
		namdInfo << "MASS:          " << data->mass[i] << "\n";
		namdInfo << "OCCUPANCY:     " << data->occupancy[i] << "\n";
		namdInfo << "BETA:          " << data->beta[i] << "\n";
		namdInfo << "ATOM NAME IND: " << data->atomNameIndexes[i] << "\n";
		namdInfo << "ATOM TYPE IND: " << data->atomTypeIndexes[i] << "\n";
		namdInfo << sendmsg;
	}

	free_vmd_static_data(data);
}
/*			END OF FUNCTION print_vmd_static_data		*/

/************************************************************************/
/*									*/
/*			FUNCTION initialize_vmd_connection		*/
/*									*/
/*	This function initializes the VMD connection.  It calls		*/
/*   the setup and initialization routines from the mdcomm library.	*/
/*   During this process, the static data and the first set of		*/
/*   dynamic data is sent.						*/
/*									*/
/************************************************************************/

void Output::initialize_vmd_connection()
{
	vmdHandle = rapp_app_setup();

	if (vmdHandle == NULL)
	{
		NAMD_die("rapp_app_setup failed to initialize VMD connection");
	}

        vmdStaticData = build_vmd_static_data(vmdStaticData);

        rapp_app_init(vmdHandle, mdcomm_send_static, mdcomm_send_dynamic,
                      vmdStaticData, vmdData, NULL);

        /*free_vmd_static_data(vmdStaticData);*/

}
/*			END OF FUNCTION initialize_vmd_connection	*/

/************************************************************************/
/*									*/
/*			FUNCTION close_vmd_connection			*/
/*									*/
/*	This function terminates the vmd_connection.			*/
/*									*/
/************************************************************************/

void Output::close_vmd_connection()
{
        if (vmdHandle == NULL)
        {
                NAMD_die("Tried to close VMD connection before it was opened");
        }

        rapp_app_exit(vmdHandle);
}
/*			END OF FUNCTION close_vmd_connection		*/

#ifdef OLD_MDCOMM
/************************************************************************/
/*									*/
/*			FUNCTION vmd_process_events			*/
/*									*/
/*	This function tells the mdcomm to process any event messages.   */
/*  This could include internal events such as RATE, KILL, etc. or      */
/*  user defined events setup using mdc_app_set_handler.		*/
/*									*/
/************************************************************************/

void Output::vmd_process_events()
{
	if (vmdArena == NULL)
	{
		return;
	}

	mdc_app_process_events(vmdArena);
}
/*			END OF FUNCTION vmd_process_events		*/

/************************************************************************/
/*									*/
/*			FUNCTION vmd_send_ascii_int			*/
/*									*/
/*   INPUTS:								*/
/*	arena - Arena structure for mdcomm				*/
/*	val - integer value to be sent					*/
/*	tag - Tag to send value with					*/
/*									*/
/*	This function sends a single integer to VMD in ascii format.    */
/*   It does this by converting it into a string and then sending	*/
/*   the string.							*/
/*									*/
/************************************************************************/

void Output::vmd_send_ascii_int(mdc_app_arena *arena, int val, int tag)
{
	char buf[64];

	sprintf(buf, "%d", val);
	mdc_app_send(arena, buf, MDC_BYTE, strlen(buf)+1, tag);
}
/*		END OF FUNCTION vmd_send_ascii_int			*/

/************************************************************************/
/*									*/
/*			FUNCTION vmd_send_ascii_float			*/
/*									*/
/*   INPUTS:								*/
/*	arena - Arena structure for mdcomm				*/
/*	val - floating point value to be sent				*/
/*	tag - Tag to send value with					*/
/*									*/
/*	This function sends a single floating point value to VMD in	*/
/*   ascii format.  It does this by converting it into a string and     */
/*   then sending the string.						*/
/*									*/
/************************************************************************/

void Output::vmd_send_ascii_float(mdc_app_arena *arena, float val, int tag)
{
	char buf[64];

	sprintf(buf, "%4.1f", val);
	mdc_app_send(arena, buf, MDC_BYTE, strlen(buf)+1, tag);
}
/*			END OF FUNCTION vmd_send_ascii_float		*/

/************************************************************************/
/*									*/
/*			FUNCTION send_vmd_static			*/
/*									*/
/*   INPUTS:								*/
/*	arena - Arena function for mdcomm calls				*/
/*									*/
/*	This function sends the static data to VMD.  It first builds	*/
/*   a structure with all the relevant information and then sends       */
/*   all the data involved using the appropriate mdcomm calls.  This    */
/*   function is called by the wrapper function vmd_send_static()	*/
/*   which is called by the mdcomm routines.				*/
/*									*/
/************************************************************************/

void Output::send_vmd_static(mdc_app_arena *arena)
{
	VmdStaticData *vdata;		//  Structure holding the data

	//  Populate the structure
	vdata = build_vmd_static_data();

	//  Send the information
	vmd_send_ascii_int(arena, VMD_NAME_LEN, XFER_NAMELEN);
	vmd_send_ascii_int(arena, vdata->numAtoms, XFER_NATOMS);
	vmd_send_ascii_int(arena, vdata->numAtomNames, XFER_NATOMNAMES);
	vmd_send_ascii_int(arena, vdata->numAtomTypes, XFER_NATOMTYPES);
	vmd_send_ascii_int(arena, 0, XFER_NHBONDDONORS);		//  NOTE: H-bonds not implemented
	vmd_send_ascii_int(arena, 0, XFER_NHBONDACCEPTORS);		//   so 0's sent instead
	vmd_send_ascii_int(arena, vdata->numBonds, XFER_NBONDS);

	mdc_app_send(arena, vdata->atomNames, MDC_BYTE, 
		     vdata->numAtomNames*VMD_NAME_LEN, XFER_ATOMNAMELIST);
	mdc_app_send(arena, vdata->atomNameIndexes, MDC_INT, 
		     vdata->numAtoms, XFER_ATOMNAMEIDX);
	mdc_app_send(arena, vdata->atomTypes, MDC_BYTE, 
		     vdata->numAtomTypes*VMD_NAME_LEN, XFER_ATOMTYPELIST);
	mdc_app_send(arena, vdata->atomTypeIndexes, MDC_INT, 
		     vdata->numAtoms, XFER_ATOMTYPEIDX);
	mdc_app_send(arena, NULL, MDC_INT, 0, XFER_DONOR);		//  NOTE:  hydrogen bonding
	mdc_app_send(arena, NULL, MDC_INT, 0, XFER_DONORH);		//   currently not implemented
	mdc_app_send(arena, NULL, MDC_INT, 0, XFER_ACCEPTOR);		//   in namd so all 0's are
	mdc_app_send(arena, NULL, MDC_INT, 0, XFER_ACCEPTORA);		//   sent here
	mdc_app_send(arena, vdata->resIds, MDC_BYTE, 
		     vdata->numAtoms*VMD_NAME_LEN, XFER_RESIDS);
	mdc_app_send(arena, vdata->resNames, MDC_BYTE, 
		     vdata->numAtoms*VMD_NAME_LEN, XFER_RESNAMES);
	mdc_app_send(arena, vdata->radii, MDC_FLOAT, 
		     vdata->numAtoms, XFER_RADII);
	mdc_app_send(arena, vdata->bonds, MDC_INT, 
		     vdata->numBonds*2, XFER_BONDS);
	mdc_app_send(arena, vdata->charge, MDC_FLOAT, 
		     vdata->numAtoms, XFER_CHARGE);
	mdc_app_send(arena, vdata->mass, MDC_FLOAT, 
		     vdata->numAtoms, XFER_MASS);
	mdc_app_send(arena, vdata->occupancy, MDC_FLOAT, 
		     vdata->numAtoms, XFER_G);
	mdc_app_send(arena, vdata->beta, MDC_FLOAT, 
		     vdata->numAtoms, XFER_K);
	mdc_app_send(arena, vdata->segIds, MDC_BYTE, 
		     vdata->numAtoms*VMD_NAME_LEN, XFER_SEGIDS);
	vmd_send_ascii_int(arena, vdata->maxNumPatches, XFER_MAXPATCHES);

	//  Free the data structure
	free_vmd_static_data(vdata);
}
/*			END OF FUNCTION send_vmd_static				*/

/********************************************************************************/
/*										*/
/*			FUNCTION send_vmd_dyn					*/
/*										*/
/*   INPUTS:									*/
/*	arena - Arena structure for mdcomm calls				*/
/*										*/
/*	This function sends the dynamic data to VMD.  It is called by the	*/
/*   wrapper function vmd_send_dyn() is called by the mdcomm functions.		*/
/*										*/
/********************************************************************************/

void Output::send_vmd_dyn(mdc_app_arena *arena)
{
	vmd_send_ascii_int(arena, vmdData->timestep, XFER_STEP);
	vmd_send_ascii_float(arena, vmdData->elapsed_time, XFER_DT);
	vmd_send_ascii_float(arena, vmdData->elapsed_cpu, XFER_CPU);
	vmd_send_ascii_float(arena, vmdData->T, XFER_TEMP);
	vmd_send_ascii_float(arena, vmdData->Epot, XFER_EPOT);
	vmd_send_ascii_float(arena, vmdData->Etot, XFER_ETOT);
	vmd_send_ascii_float(arena, vmdData->Evdw, XFER_EVDW);
	vmd_send_ascii_float(arena, vmdData->Eelec, XFER_EELEC);
	vmd_send_ascii_float(arena, 0.0, XFER_ESUP);
	vmd_send_ascii_float(arena, 0.0, XFER_EHBO);
	vmd_send_ascii_float(arena, vmdData->Ebond, XFER_EBOND);
	vmd_send_ascii_float(arena, vmdData->Eangle, XFER_EANGLE);
	vmd_send_ascii_float(arena, vmdData->Edihe, XFER_EDIHED);
	vmd_send_ascii_float(arena, vmdData->Eimpr, XFER_EIMPR);
	mdc_app_send(arena, vmdData->X, MDC_FLOAT, vmdData->numAtoms, 
		     XFER_XCOORDS);
	mdc_app_send(arena, vmdData->Y, MDC_FLOAT, vmdData->numAtoms, 
		     XFER_YCOORDS);
	mdc_app_send(arena, vmdData->Z, MDC_FLOAT, vmdData->numAtoms, 
		     XFER_ZCOORDS);

	vmd_send_ascii_int(arena, vmdData->numPatches, XFER_NPATCHES);

	mdc_app_send(arena, vmdData->pXOrigins, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHX);
	mdc_app_send(arena, vmdData->pYOrigins, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHY);
	mdc_app_send(arena, vmdData->pZOrigins, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHZ);

	mdc_app_send(arena, vmdData->patchLength, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHXLEN);
	mdc_app_send(arena, vmdData->patchWidth, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHYLEN);
	mdc_app_send(arena, vmdData->patchHeight, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHZLEN);

	mdc_app_send(arena, vmdData->patchLoads, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHLOAD);
	mdc_app_send(arena, vmdData->patchAtomNums, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHATOMS);
	mdc_app_send(arena, vmdData->patchNode, MDC_FLOAT, vmdData->numPatches,
		     XFER_PATCHNODE);
}
/*			END OF FUNCITON send_vmd_dyn				*/

/********************************************************************************/
/*										*/
/*				FUNCTION vmd_send_static			*/
/*										*/
/*   INPUTS:									*/
/*	arena  - Arena structure used by mdcomm functions			*/
/*										*/
/*	This is a wrapper function for the function Output::send_vmd_static().  */
/*   It only exists since the mdcomm routines want a function pointer passed    */
/*   to them, and it is impossible to send them a member function of the 	*/
/*   Output object.  So instead, we give them this function that does nothing   */
/*   more than make the appropriate call to the Output object.			*/
/*										*/
/********************************************************************************/

int vmd_send_static(mdc_app_arena *arena)
{
	namdMyNode->output->send_vmd_static(arena);

	return(0);
}
/*				END OF FUNCTION vmd_send_static			*/

/************************************************************************/
/*									*/
/*				FUNCTION vmd_send_dyn			*/
/*									*/
/*   INPUTS:								*/
/*	arena  - Arena structure used by mdcomm functions		*/
/*									*/
/*	This is a wrapper function for the function 			*/
/*   Output::send_vmd_dyn().  It only exists since the mdcomm routines  */
/*   want a function pointer passed to them, and it is impossible to    */
/*   send them a member function of the Output object.  So instead, we  */
/*   give them this function that does nothing more than make the       */
/*   appropriate call to the Output object.				*/
/*									*/
/************************************************************************/

int vmd_send_dyn(mdc_app_arena *arena)
{
	namdMyNode->output->send_vmd_dyn(arena);

	return(0);
}
/*			END OF FUNCTION vmd_send_dyn			*/

#endif

/************************************************************************/
/*									*/
/*			FUNCTION recv_vmd_patch_loads			*/
/*									*/
/*	This function receives the load information from each patch for */
/*  output to VMD.							*/
/*									*/
/************************************************************************/

void Output::recv_vmd_patch_loads()
{
	int notReceived=namdNumNodes;	//  Number of nodes we are waiting for
	int tag=PATCHLOADTAG;		//  tag for load messages
	int node;			//  Node we received message from
	Message *msg;			//  Message received
	int i;				//  Loop counter
	int num_patches;		//  Number of patches in the current
					//  message
	Real *loads;			//  Array of load information
	int *patches;			//  Array of patch numbers
	int *atomsPerPatch;		//  Number of atoms per patch

	//  Loop while we are still waiting for at least one node
	while (notReceived)
	{
		//  Receive a patch message from any node
		node = -1;

		while ( (msg=comm->receive(node, tag)) != NULL)
		{
			//  Get the number of patches in this message
			msg->get(num_patches);

			//  If there is any real load information, get it
			if (num_patches)
			{
				patches = msg->get_intlist_by_ref();
				loads = msg->get_reallist_by_ref();
				atomsPerPatch = msg->get_intlist_by_ref();

				//  Loop through and put the information into
				//  the arrays of data for all patches
				for (i=0; i<num_patches; i++)
				{
					vmdData->patchLoads[patches[i]] = loads[i];
					vmdData->patchAtomNums[patches[i]] = atomsPerPatch[i];
				}
			}

			delete msg;

			notReceived--;

			node=-1;
		}
	}
}
/*			END OF FUNCTION recv_vmd_patch_loads			*/

#endif 	/*  MDCOMM  */

/********************************************************************************/
/*										*/
/*				FUNCTION long_force				*/
/*										*/
/*   INPUTS:									*/
/*	tiemstep - Current timestep value					*/
/*	n - Number of forces							*/
/*	f - Array of force vectors						*/
/*										*/
/*	This is the public method that is used to take the forces from the	*/
/*   collect object and output them into a force DCD.  It does nothing more	*/
/*   than call the force DCD generation routine.				*/
/*										*/
/********************************************************************************/

void Output::long_force(int timestep, int n, Vector *f)
   
{
   output_longforcedcdfile(timestep, n, f);
}
/*			END OF FUNCTION long_force				*/

/********************************************************************************/
/*										*/
/*				FUNCTION short_force				*/
/*										*/
/*   INPUTS:									*/
/*	tiemstep - Current timestep value					*/
/*	n - Number of forces							*/
/*	f - Array of force vectors						*/
/*										*/
/*	This is the public method that is used to take the forces from the	*/
/*   collect object and output them into a force DCD.  It does nothing more	*/
/*   than call the force DCD generation routine.				*/
/*										*/
/********************************************************************************/

void Output::short_force(int timestep, int n, Vector *f)
   
{
   output_shortforcedcdfile(timestep, n, f);
}
/*			END OF FUNCTION force					*/

/********************************************************************************/
/*										*/
/*				FUNCTION all_force				*/
/*										*/
/*   INPUTS:									*/
/*	tiemstep - Current timestep value					*/
/*	n - Number of forces							*/
/*	f - Array of force vectors						*/
/*										*/
/*	This is the public method that is used to take the forces from the	*/
/*   collect object and output them into a force DCD.  It does nothing more	*/
/*   than call the force DCD generation routine.				*/
/*										*/
/********************************************************************************/

void Output::all_force(int timestep, int n, Vector *f)
   
{
   output_allforcedcdfile(timestep, n, f);
}
/*			END OF FUNCTION all_force				*/

/********************************************************************************/
/*										*/
/*			FUNCTION output_longforcedcdfile			*/
/*										*/
/*   INPUTS:									*/
/*	timestep - Current timestep value					*/
/*	n - Number of forces							*/
/*	forces - Array of force vectors						*/
/*										*/
/*	Output forces to a force DCD file.					*/
/*										*/
/********************************************************************************/

void Output::output_longforcedcdfile(int timestep, int n, Vector *forces)

{
	static Bool first=TRUE;	//  Flag indicating first call
	static int fileid;	//  File id for the dcd file
	static float *x, *y, *z; // Arrays to hold x, y, and z arrays
	int i;			//  Loop counter
	int ret_code;		//  Return code from DCD calls
	char filename[257];	//  DCD filename

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


		//  Write out the header
		ret_code = write_dcdheader(fileid, filename, n, 
				(int) ((namdMyNode->simParams->N+1)/namdMyNode->simParams->electForceDcdFrequency)+1,
				0, namdMyNode->simParams->electForceDcdFrequency, 
				namdMyNode->simParams->dt/TIMEFACTOR);


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
/*			END OF FUNCTION output_longforcedcdfile			*/

/********************************************************************************/
/*										*/
/*			FUNCTION output_shortforcedcdfile			*/
/*										*/
/*   INPUTS:									*/
/*	timestep - Current timestep value					*/
/*	n - Number of forces							*/
/*	forces - Array of force vectors						*/
/*										*/
/*	Output forces to a short force DCD file.				*/
/*										*/
/********************************************************************************/

void Output::output_shortforcedcdfile(int timestep, int n, Vector *forces)

{
	static Bool first=TRUE;	//  Flag indicating first call
	static int fileid;	//  File id for the dcd file
	static float *x, *y, *z; // Arrays to hold x, y, and z arrays
	int i;			//  Loop counter
	int ret_code;		//  Return code from DCD calls
	char filename[257];

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

		//  Write out the header
		ret_code = write_dcdheader(fileid, filename, n, 
				(int) ((namdMyNode->simParams->N+1)/namdMyNode->simParams->electForceDcdFrequency)+1,
				0, namdMyNode->simParams->electForceDcdFrequency, 
				namdMyNode->simParams->dt/TIMEFACTOR);


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
/*			END OF FUNCTION output_shortforcedcdfile		*/

/********************************************************************************/
/*										*/
/*			FUNCTION output_allforcedcdfile				*/
/*										*/
/*   INPUTS:									*/
/*	timestep - Current timestep value					*/
/*	n - Number of forces							*/
/*	forces - Array of force vectors						*/
/*										*/
/*	Output forces to a force DCD file.					*/
/*										*/
/********************************************************************************/

void Output::output_allforcedcdfile(int timestep, int n, Vector *forces)

{
	static Bool first=TRUE;	//  Flag indicating first call
	static int fileid;	//  File id for the dcd file
	static float *x, *y, *z; // Arrays to hold x, y, and z arrays
	int i;			//  Loop counter
	int ret_code;		//  Return code from DCD calls

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


		//  Write out the header
		ret_code = write_dcdheader(fileid, namdMyNode->simParams->allForceDcdFilename, n, 
				(int) ((namdMyNode->simParams->N+1)/namdMyNode->simParams->allForceDcdFrequency)+1,
				0, namdMyNode->simParams->allForceDcdFrequency, 
				namdMyNode->simParams->dt/TIMEFACTOR);


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
/*			END OF FUNCTION output_allforcedcdfile			*/


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Output.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.9 $	$Date: 1997/09/18 22:22:21 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Output.C,v $
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
 * works on both 64 and 32-bit machines.  Also, I made Inform.C use CPrintf,
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
