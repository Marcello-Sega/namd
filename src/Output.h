/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   This object outputs the data collected on the master  node
*/

#ifndef OUTPUT_H
#define OUTPUT_H

#include "common.h"

class Vector;
class FloatVector;

class Output 
{

private:
   void output_dcdfile(int, int, FloatVector *);  //  output coords to dcd file
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
   void coordinate(int, int, Vector *, FloatVector *);
						//  Produce appropriate 
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

