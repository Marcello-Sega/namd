/***************************************************************************/
/*                                                                         */
/*    (C) Copyright 1995,1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 * This object outputs the data collected on the master  node
 ***************************************************************************/

#ifndef OUTPUT_H

#define OUTPUT_H

#include "common.h"

#ifdef MDCOMM
#include <rapp_app.h>                           // RAPP application defns
#include "mdcomm.h"                             // NAMD-specific RAPP function
#endif

class Vector;

#define OUTPUT_TEMPERATURE_WARNING 1000.0  /* warn if temp > this value */
#define VMD_NAME_LEN	4		   /* Length of name fields for VMD */

#ifdef OLD_VMD
//  This is all the static data that needs to be transferred to VMD
typedef struct vmd_static_data
{
	int numAtoms;			//  Total number of atoms
	int numAtomNames;		//  Number of unique atom names
	int numAtomTypes;		//  Number of unique atom types
	int numBonds;			//  Number of linear bonds
	int maxNumPatches;		//  Max number of patches
	char *atomNames;		//  List of unique atom names
	int *atomNameIndexes;		//  Indexes for each atom into the
					//  atomNames array
	char *atomTypes;		//  List of unique atom types
	int *atomTypeIndexes;		//  Indexes for each atom into the
					//  atomTypes array
	int *bonds;			//  List of bonded atoms
	char *resIds;			//  Residue ids for each atom
	char *resNames;			//  Residue names for each atom
	float *radii;			//  Calculated vdw radii
	float *charge;			//  Charge for each atom
	float *mass;			//  Mass for each atom
	float *occupancy;		//  Occupancy for each atom
	float *beta;			//  Beta coupling for each atom
	char *segIds;			//  Segment IDs for each atom
} VmdStaticData;

// This structure is used to transfer dynamic data to VMD
typedef struct vmd_dyn_data
{
	Bool energiesArrived;	//  Flag TRUE->energies have been set
	Bool coordsArrived;	//  Flag TRUE->coords have been set
	int numAtoms;		//  Number of atoms in simulation
	int timestep;		//  Current timestep
	float elapsed_time;	//  Total elapsed simulation time
	float elapsed_cpu;	//  Total elapsed cpu time - CURRENTLY NOT IMPLEMENTED
	float T;		//  Temperature
	float Epot;		//  Potential energy
	float Etot;		//  Total energy
	float Evdw;		//  Van der Waals energy
	float Eelec;		//  Electrostatic energy
	float Ebond;		//  Linear bond energy
	float Eangle;		//  Angle bond energy
	float Edihe;		//  Dihedral bond energy
	float Eimpr;		//  Improper bond energy
	float *X;		//  Array of X coordinates
	float *Y;		//  Array of Y coordinates
	float *Z;		//  Array of Z coordinates
	int numPatches;		//  Number of patches in simulation
	float *pXOrigins;	//  Origins of patches
	float *pYOrigins;	//  Origins of patches
	float *pZOrigins;	//  Origins of patches
	float *patchLength;	//  Length of each patch
	float *patchWidth;	//  Length of each patch
	float *patchHeight;	//  Length of each patch
	float *patchAtomNums;	//  Number of atoms on each patch
	float *patchLoads;	//  Number of atoms on each patch
	float *patchNode;	//  Number of atoms on each patch
} VmdDynData;
#endif // OLD_VMD

class Output 
{

private:
   int teTempWarning;                           //  flag for large temp warn.

#ifdef MDCOMM
   VmdDynData *vmdData;				//  Dynamic data to pass to VMD
   VmdStaticData *vmdStaticData;                //  Static data to pass to VM
   rapp_app_handle *vmdHandle;                  //  Structure used to 
						//  communicate with mdcomm 
						//  routines
#endif 

   void output_dcdfile(int, int, Vector *);  	//  output coords to dcd file
   void output_veldcdfile(int, int, Vector *); 	//  output velocities to
						//  dcd file
   void output_longforcedcdfile(int, int, Vector *); // output long range elect force 
						//  to dcd file
   void output_shortforcedcdfile(int, int, Vector *); // output short range elect force 
						//  to dcd file
   void output_allforcedcdfile(int, int, Vector *); // output total forces to dcd file

   void output_restart_coordinates(Vector *, int, int);
						//  output coords to 
						//  restart file
   void output_restart_velocities(int, int, Vector *);
						//  output velocities to 
						//  restart file
   void output_final_coordinates(Vector *, int, int);//  output final coordinates
   void output_final_velocities(int, int, Vector *);	//  output final coordinates

   void scale_vels(Vector *, int, Real);	//  scale velocity vectors before output
   void write_binary_file(char *, int, Vector *); // Write a binary restart file with
						//  coordinates or velocities

#ifdef MDCOMM
   VmdStaticData *build_vmd_static_data(VmdStaticData *);
                                                //  Build the static data 
						//  structures for VMD

   void free_vmd_static_data(VmdStaticData *);	//  Free the static data 
						//  structures for VMD

   void insert_into_unique_list(const char *, char *, int&); 
						//  Maintain a unique list of
						//  strings

   int search_list(const char *, const char *, int);  
						//  Do a binary search on
						//  a list of strings

   int name_cmp(const char *, const char *);	//  Do a compare of name strings


   void gather_vmd_energies(int, BigReal *, BigReal, BigReal);
						//  Collect energy values for VMD
   void gather_vmd_coords(int, int, Vector *);
						//  Collect coordinates for VMD

#ifdef VMD_OLD
   void vmd_send_ascii_int(mdc_app_arena *, int, int);
						//  Send an integer value in ascii
						//  to VMD
   void vmd_send_ascii_float(mdc_app_arena *, float, int);
						//  Send a floating point value
						//  to VMD in ascii
#endif // VMD_OLD
#endif

public :
   Output();					//  Constructor
   ~Output();					//  Destructor
   void energy(int, BigReal *);			//  Output energies

   static int coordinateNeeded(int);
   void coordinate(int, int, Vector *);		//  Produce appropriate 
						//  coordinate output for 
						//  the current timestep
   static int velocityNeeded(int);
   void velocity(int, int, Vector *);		//  Produce appropriate velocity
						//  output for the current 
						//  timestep
   void long_force(int, int, Vector *);		//  Produce a long range elect force 
						//  output for the current timestep
   void short_force(int, int, Vector *);	//  Produce a short range elect force 
						//  output for the current timestep
   void all_force(int, int, Vector *);		//  Produce a total force output for 
						//  the current timestep

#ifdef MDCOMM
   void initialize_vmd_connection();		//  Startup the VMD connection
   void close_vmd_connection();			//  Close the VMD connection
   void vmd_process_events();			//  Process events from VMD

#ifdef VMD_OLD
   void send_vmd_static(mdc_app_arena *);	//  Send static data to VMD
   void send_vmd_dyn(mdc_app_arena *);		//  Send dynamic data to VMD
#endif
   void send_vmd_static(rapp_active_socket_t *, void *);
                                                //  Send static data to VMD
   void send_vmd_dyn(rapp_active_socket_t *, void *, void *);
                                                //  Send dynamic data to VMD


   void print_vmd_static_data();		//  Debugging routine to print
						//  VMD static info to screen
   void recv_vmd_patch_loads();			//  Collect Patch load messages
#endif /* MDCOMM */
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Output.h,v $
 *	$Author: milind $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1997/04/04 23:34:23 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Output.h,v $
 * Revision 1.4  1997/04/04 23:34:23  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.3  1997/03/19 11:54:38  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1  1997/02/11 22:56:16  jim
 * Added dcd file writing.
 *
 * Revision 1.19  1996/01/28 21:50:08  jean
 * Attempting to make stable RCS without Mark Nelson's
 * fma/pairlist decoupling
 *
 * Revision 1.19  1995/12/05 22:13:34  brunner
 * Rearranged #ifdef MDCOMM 's so it will compile without MDCOMM defined
 *
 * Revision 1.18  1995/10/26 16:24:28  dalke
 * added new mdcomm code
 *
 * Revision 1.17  1995/09/30  19:55:53  billh
 * Added ETITLE and ENERGY strings  in energy output records.
 * Added tiny Elvis warning if temperature exceeds 1000 degrees.
 *
 * Revision 1.16  95/09/26  15:15:05  15:15:05  nelson (Mark T. Nelson)
 * Added all force DCD files and cleaned up symantics of long and short
 * range electrostatic force DCD files.
 * 
 * Revision 1.15  95/08/30  14:03:31  14:03:31  nelson (Mark T. Nelson)
 * Added short range DCD file output and binary coordinate restart files
 * 
 * Revision 1.14  95/08/11  14:52:29  14:52:29  nelson (Mark T. Nelson)
 * Added routines for the generation of force DCD files
 * 
 * Revision 1.13  95/07/14  14:19:09  14:19:09  nelson (Mark T. Nelson)
 * Added binary velocity restart files
 * 
 * Revision 1.12  95/04/06  12:52:05  12:52:05  nelson (Mark T. Nelson)
 * Removed extern class references
 * 
 * Revision 1.11  95/03/08  14:47:39  14:47:39  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.10  95/02/17  12:37:36  12:37:36  nelson (Mark T. Nelson)
 * Added shifting factor for velocity PDB files
 * 
 * Revision 1.9  94/11/23  09:37:09  09:37:09  nelson (Mark T. Nelson)
 * Added recv_vmd_patch_loads()
 * 
 * Revision 1.8  94/11/10  20:27:55  20:27:55  nelson (Mark T. Nelson)
 * Added data for Patch information being sent to VMD
 * 
 * Revision 1.7  94/11/09  13:28:05  13:28:05  nelson (Mark T. Nelson)
 * Changed elapsed_cpu to elapsed_time
 * 
 * Revision 1.6  94/10/31  19:25:31  19:25:31  nelson (Mark T. Nelson)
 * Added functions for VMD interface
 * 
 * Revision 1.5  94/10/18  20:32:55  20:32:55  nelson (Mark T. Nelson)
 * Added functions to output velocity DCD files
 * 
 * Revision 1.4  94/10/14  09:37:21  09:37:21  nelson (Mark T. Nelson)
 * Added output of final coordinate and velocity files as well
 * as rearranging and commenting other output routines
 * 
 * Revision 1.3  94/10/06  21:37:25  21:37:25  nelson (Mark T. Nelson)
 * added DCD file interface
 * 
 * Revision 1.2  94/10/06  15:50:47  15:50:47  gursoy (Attila Gursoy)
 * coordinates and velocities are passed as an array (it was a message previously
 * 
 * Revision 1.1  94/09/29  12:04:56  12:04:56  gursoy (Attila Gursoy)
 * Initial revision
 * 
 ***************************************************************************/
