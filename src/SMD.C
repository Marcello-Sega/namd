/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/
       
#include "common.h"
#include "NamdOneTools.h"
#include "Node.h"
#include "Output.h"
#include "PDB.h"
#include "PDBData.h"
#include "SimParameters.h"
#include "SMD.h"
#include "SMDMsgs.h"
#include "Node.top.h"
#include "Debug.h"
#include <math.h>

/************************************************************************/
/*									*/
/*			FUNCTION SMDData      			        */
/*									*/
/************************************************************************/
SMDData::SMDData(SimParameters *simP) {

  simParams = simP;

}

// initialize on node 0
SMDData::init(PDB *pdb) {

  timeStamp = simParams->firstTimestep;
  direction = simParams->SMDDir;
  refPos = simParams->SMDRefPos;

  int atomNum = simParams->SMDAtom;

  PDBAtom *atom =  pdb->atom(atomNum);

  atomPosVmin.x = atomPosVmax.x = atom->xcoor();
  atomPosVmin.y = atomPosVmax.y = atom->ycoor();
  atomPosVmin.z = atomPosVmax.z = atom->zcoor();

}
 
// update the data if necessary; 
// t = current timestep
// p = current atom position
void SMDData::update (int t, const Vector &p) {
  BigReal aveVel;
  Vector vel;
  int doSend = 0;

  // don't do anything on the very first timestep
  if (t == simParams->firstTimestep) return;

  // check if we need to change direction
  if (simParams->SMDChDirOn) {
    if (t % simParams->SMDVminTave == 0) {
      vel = (p - atomPosVmin)/simParams->SMDVminTave;
      aveVel = vel.dot(direction);
      if (aveVel < simParams->SMDVmin) {
	choose_direction(t, p);
      }
      atomPosVmin = p;
      doSend = 1;
    }
  }

  // check if we need to reset force
  if (simParams->SMDChForceOn) {
    if (t % simParams->SMDVmaxTave == 0) {
      vel = (p - atomPosVmax)/simParams->SMDVmaxTave;
      aveVel = vel.dot(direction);
      if (aveVel > simParams->SMDVmax) {
	reset_force(t, p);
      }
      atomPosVmax = p;
      doSend = 1;
    }
  }

  if (doSend) {
    sendData(t);
  }
  
} // update()


// choose a new direction for restraint point movement
// this also resets force to Fmin
void SMDData::choose_direction(int t, const Vector &p) {
  int method = simParams->SMDChDirMethod;
  BigReal coneAngle = simParams->SMDConeAngle;
  Vector initDir = simParams->SMDDir;
  Vector a;  // a vector orthogonal to initDir
  Vector b;  // a vector forming the right coord. system with initDir and a
  BigReal phi;   // phi angle [0, 2*Pi] with respect to vector a in a-b plane
  BigReal theta; // theta angle with respect to init direction 
  BigReal Pi = 4.*atan(1.);

  // get the random values for theta and phi
  phi = NAMD_random() * 2. * Pi;

  switch(method) {
  case SMD_GAUSSIAN:
    theta = coneAngle + 1.0;
    while (theta > coneAngle) { // reject everything that's beyond coneAngle
      theta = simParams->SMDGaussW * fabs(gaussian_random_number());
    }
    break;
  case SMD_UNIFORM:
  default:
    theta = NAMD_random() * coneAngle;
    break;
  }
  theta *= Pi/180.0;   // transform to radians

  // now we need to produce a vector that forms angle theta with initDir
  
  // first, construct a vector that is perpendicular to initDir.
  Vector tmp(0.0, 1.0, 1.0); // temporary vector
  a = cross(initDir, tmp);
  if (a.length2() == 0) { // tmp is collinear to initDir
    tmp.x = 1.0; // now it is not, so try again
    a = cross(initDir, tmp);
  }
  a = a.unit();  // normalize the vector

  // complete the coordinate system with the third vector
  b = cross(initDir, a);

  // and the new direction is...
  direction = a * cos(phi) * sin(theta) 
    + b * sin(phi) * sin(theta) + initDir * cos(theta);
  direction = direction.unit();
  
  // Reset the force in the new direction.
  reset_force(t, p);
  
}

// reset the force to Fmin in the current direction
void SMDData::reset_force(int t, const Vector &p) {
  // A   = A +       pN          / (kcal/mol/A^2) / (kcal/mol/(pN*A)) * A
  refPos = p + (simParams->SMDFmin / simParams->SMDk / 69.479) * direction;
  timeStamp = t;
  
  // print some output 
  output_new_refpos(t, p);

}

// output the data
void SMDData::output(int t, const Vector &p, const Vector &f) {
  if ( t % (10 * simParams->SMDOutputFreq) == 0 ) {
    iout <<  "SMDTITLE: TS     CURRENT_POSITION"
	 <<  "                       FORCE"
	 <<  "    RESTRAINT_REF_POSITION"
	 <<  "              CURRENT_DIRECTION" 
	 <<  "        TSTAMP"
	 <<  "\n" << endi;
  }

  Vector fpN = f*PNPERKCALMOL;
  iout << "SMD " 
       << t << "    " 
       << p.x << " " << p.y << " " << p.z << "     "
       << fpN.x << " " << fpN.y << " " << fpN.z << "    "
       << refPos.x << " " << refPos.y << " " << refPos.z << "    "
       << direction.x << " " << direction.y << " " << direction.z << "    "
       << timeStamp << "\n" << endi;
}

// output info about changed direction and/or refPos
void SMDData::output_new_refpos(int t, const Vector &p) {
  
  static int titledone = 0;

  if (!titledone) {
    iout << "SMDCTITLE: TS          "
	 << "CURRENT_POSITION                "
	 << "RESTRAINT_REF_POS               "
	 << "DIRECTION               "
	 << "TSTAMP"  
	 << "\n" << endi;
    titledone = 1;
  } 

  iout << "SMDChange   " 
       << t << "    " 
       << p.x << " " << p.y << " " << p.z << "     "
       << refPos.x << " " << refPos.y << " " << refPos.z << "    "
       << direction.x << " " << direction.y << " " << direction.z << "    "
       << "\n" << endi;
}


// send the SMD data to other nodes
void SMDData::sendData(int t) {
  SMDDataMsg *msg = new (MsgIndex(SMDDataMsg)) SMDDataMsg;

  msg->curTime = t;
  msg->timeStamp = timeStamp;
  msg->direction = direction;
  msg->refPos = refPos;
  msg->atomPosVmin = atomPosVmin;
  msg->atomPosVmax = atomPosVmax;
  
  int pe = CMyPe();
  Node::Object()->sendSMDData(msg);
}


// receive the SMD data from other Nodes 
void SMDData::recvData(SMDDataMsg *msg) {
  int t;
  t = msg->curTime;  // t is wasted, but can be used for error checking
  timeStamp = msg->timeStamp;
  direction = msg->direction;
  refPos = msg->refPos;
  atomPosVmin = msg->atomPosVmin;
  atomPosVmax = msg->atomPosVmax;

  int pe = CMyPe();

  delete msg;
}



// send the SMD data from Node 0
void SMDData::send_SMDData(Communicate *com_obj) {
   //  Message to send to clients
  MOStream *msg=com_obj->newOutputStream(ALLBUTME, MOLECULETAG, BUFSIZE);
  
  if ( msg == NULL )
  {
    NAMD_die("Memory allocation failed in SMDData::send_SMDData");
  }

  msg->put(timeStamp);
  msg->put(&direction);
  msg->put(&refPos);
  msg->put(&atomPosVmin);
  msg->put(&atomPosVmax);


}  


// receive the SMD data from Node 0
void SMDData::receive_SMDData(MIStream *msg) {

  msg->get(timeStamp);
  msg->get(&direction);
  msg->get(&refPos);
  msg->get(&atomPosVmin);
  msg->get(&atomPosVmax);

  delete msg;
}
