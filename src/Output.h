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
class Lattice;

// semaphore "steps", must be negative
#define FILE_OUTPUT -1
#define END_OF_RUN -2
#define EVAL_MEASURE -3

class Output 
{

private:

   //  output coords to dcd file
   //  Pass non-NULL Lattice to include unit cell in the timesteps.
   void output_dcdfile(int, int, FloatVector *, const Lattice *); 
   void output_veldcdfile(int, int, Vector *); 	//  output velocities to
						//  dcd file

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
   void coordinate(int, int, Vector *, FloatVector *, Lattice &);
						//  Produce appropriate 
						//  coordinate output for 
						//  the current timestep
   static int velocityNeeded(int);
   void velocity(int, int, Vector *);		//  Produce appropriate velocity
						//  output for the current 
						//  timestep
};

#endif

