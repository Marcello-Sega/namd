/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/
       
#include "ComputeSMD.h"
#include "Node.h"
#include "SMD.h"
#include "SimParameters.h"
#include "Patch.h"

/************************************************************************/
/*									*/
/*			FUNCTION ComputeSMD	   		        */
/*									*/
/************************************************************************/

ComputeSMD::ComputeSMD(ComputeID c, PatchID pid)
  : ComputePatch(c,pid)
{
	reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

	SimParameters *simParams = Node::Object()->simParameters;

	//  Get parameters from the SimParameters object
	consExp = simParams->SMDExp;
	k = simParams->SMDk;
	moveVel = simParams->SMDVel;
	moveAtom = simParams->SMDAtom;
	outputFreq = simParams->SMDOutputFreq;
	chDirOn = simParams->SMDChDirOn;
	chForceOn =  simParams->SMDChForceOn;
	

}
/*			END OF FUNCTION ComputeSMD		*/

/************************************************************************/
/*									*/
/*			FUNCTION ~ComputeSMD			        */
/*									*/
/*	This is the destructor for the ComputeSMD force object.	        */
/*									*/
/************************************************************************/

ComputeSMD::~ComputeSMD()

{
  delete reduction;
}
/*			END OF FUNCTION ~ComputeSMD		       */

/************************************************************************/
/*									*/
/*				FUNCTION force				*/
/*									*/
/************************************************************************/

void ComputeSMD::doForce(Position* p, Results* res, AtomProperties* a)

{
	SMDData *smdData = Node::Object()->smdData;
	Vector refPos;		//  Reference position
	BigReal r, r2; 	//  r=distance between atom position and the
			//  reference position, r2 = r^2
	Vector Rij;	//  vector between current position and reference pos
	BigReal value;	//  Current calculated value

	// aliases to work with old code
	Force *f = res->f[Results::normal];
	BigReal energy = 0;

	int currentTime = patch->flags.step;
	int timeStamp; 
	Vector moveDir;

	for (int localID=0; localID<numAtoms; ++localID) {	  
	  if (a[localID].id == moveAtom) { 
	    // this call will do all the necessary changes
	    smdData->update(currentTime, p[localID]);

	    // now calculate the force and energy
	    smdData->get_smd_params(timeStamp, moveDir, refPos);
	    Rij = refPos + (currentTime - timeStamp) * moveVel * moveDir
	      - p[localID];

	    //  Calculate the distance and the distance squared
	    r2 = Rij.length2();
	    r = sqrt(r2);

	    //  Only calculate the energy and the force if the distance is
	    //  non-zero.   Otherwise, you could end up dividing by 0, which
	    //  is bad
	    if (r>0.0)
	    {
	      value=k;
      
	      //  Loop through and multiple k by r consExp times.
	      //  i.e.  calculate kr^e/e
	      //  I know, I could use pow(), but I don't trust it.
	      for (int k=0; k<consExp; ++k)
	      {
		value *= r;
	      }

	      //  Add to the energy total
	      energy += value/consExp;
      
	      //  Now calculate the force, which is kr^(e-1).  Also, divide
	      //  by another factor of r to normalize the vector before we
	      //  multiple it by the magnitude
	      value /= r2;
      
	      Rij.mult(value);
      
	      f[localID] += Rij;

	    }
	    
	    if (currentTime % outputFreq == 0) {
	      smdData->output(currentTime, p[localID], Rij);
	    }
	    
	    break;  // change this when there are many atoms restrained.

	  } // if found atom
	} // loop

	reduction->item(REDUCTION_SMD_ENERGY) += energy;
	reduction->submit();

}
/*			END OF FUNCTION force				*/

