//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/
 
/***************************************************************************
 * DESCRIPTION:
 *	SMD provides moving atomic restraints (harmonic constraints).  See the
 * Programmer's Guide for details.
 *
 ***************************************************************************/

#ifndef SMD_H
#define SMD_H

#include "Vector.h"
#include "common.h"

class SMDDataMsg;
class SimParameters;
class MIStream;
class MOStream;

class SMDData {
public:
  SMDData(SimParameters *simP);
  ~SMDData() {};

  // these functions are for sending a few numbers every ?? timesteps
  void sendData(int t);  //  receive the SMD data  
  void recvData(SMDDataMsg *msg);
                         //  receive the SMD data  

  void init(PDB *); // initialize on Node 0

  // these functions are for sending the initial SMDData obj
  // from Node 0 to other Nodes
  void send_SMDData(Communicate *com_obj);
  void receive_SMDData(MIStream *msg);


  void output(int t, const Vector &p, const Vector &f);
                         //  output the data

  void update (int t, const Vector &p);   
                         // update the data if necessary

  // get the data members for force calculation
  void get_smd_params(int &t, Vector &dir, Vector &refpos) 
    const {
    t = timeStamp;
    dir = direction;
    refpos = refPos;
  }

 
private:
  int timeStamp;          // time the refPos was last changed
  Vector direction;       // direction of restraint point movement
  Vector refPos;          // restraint point position
  Vector atomPosVmin;     // restrained atom position for Vmin averaging
  Vector atomPosVmax;     // restrained atom position for Vmax averaging

  SimParameters *simParams; // for convinience

  void choose_direction(int t, const Vector &p);
                         // choose new direction of the movement

  void reset_force(int t, const Vector &p);
                         // reset the force to Fmin
  
  void output_new_refpos(int t, const Vector &p);
                         // output info about changed direction and/or refPos

};


#endif /* SMD_H */
