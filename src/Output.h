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

#ifdef MEM_OPT_VERSION
class ParOutput{
private:
    void output_veldcdfile_master(int timestep, int n);
    void output_veldcdfile_slave(int timestep, int fID, int tID, Vector *vecs);
    void output_restart_velocities_master(int timestep, int n);
    void output_restart_velocities_slave(int timestep, int parN, Vector *vecs, int64 offset);
    void output_final_velocities_master(int n);
    void output_final_velocities_slave(int parN, Vector *vecs, int64 offset);

    void output_dcdfile_master(int timestep, int n, const Lattice *lat);
    void output_dcdfile_slave(int timestep, int fID, int tID, FloatVector *fvecs);
    void output_restart_coordinates_master(int timestep, int n);
    void output_restart_coordinates_slave(int timestep, int parN, Vector *vecs, int64 offset);
    void output_final_coordinates_master(int n);
    void output_final_coordinates_slave(int parN, Vector *vecs, int64 offset);

    void write_binary_file_master(char *fname, int n);
    void write_binary_file_slave(char *fname, int parN, Vector *vecs, int64 offset);

    int dcdFileID;
    Bool dcdFirst;
    float *dcdX, *dcdY, *dcdZ;

    int veldcdFileID;
    Bool veldcdFirst;    
    float *veldcdX, *veldcdY, *veldcdZ;

public:
    ParOutput(){
        dcdFileID=veldcdFileID=-99999;
        dcdFirst=veldcdFirst=TRUE;
        dcdX=dcdY=dcdZ=veldcdX=veldcdY=veldcdZ=NULL;        
    }
    ~ParOutput() {}

    void velocityMaster(int timestep, int n);
    void velocitySlave(int timestep, int fID, int tID, Vector *vecs);

    void coordinateMaster(int timestep, int n, Lattice &lat);
    void coordinateSlave(int timestep, int fID, int tID, Vector *vecs, FloatVector *fvecs);
};
#endif

#endif

