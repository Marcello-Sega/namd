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

class Vector;

class Output 
{

private:
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
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Output.h,v $
 *	$Author: justin $	$Locker:  $		$State: Exp $
 *	$Revision: 1.7 $	$Date: 1999/09/02 23:04:51 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Output.h,v $
 * Revision 1.7  1999/09/02 23:04:51  justin
 * Eliminated MDComm from all files and Makefiles
 *
 * Revision 1.6  1998/04/15 22:13:52  jim
 * Make depends returns same results regardless of DPME, DPMTA, TCL or MDCOMM.
 *
 * Revision 1.5  1998/04/14 03:19:23  jim
 * Fixed up MDCOMM code.
 *
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
