/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

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
  : ComputeHomePatch(c,pid)
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
 	projectForce = simParams->SMDProjectForce;
	

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

void ComputeSMD::doForce(Position* p, Results* res, AtomProperties* a, Transform* t)

{
	SMDData *smdData = Node::Object()->smdData;
	Lattice &lattice = patch->lattice;
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
	    Position atomPos = lattice.reverse_transform(p[localID],t[localID]);

	    // this call will do all the necessary changes
	    smdData->update(currentTime, atomPos);

	    // now calculate the force and energy
	    smdData->get_smd_params(timeStamp, moveDir, refPos);
	    Rij = refPos + (currentTime - timeStamp) * moveVel * moveDir
	      - atomPos;

	    //  Calculate the distance and the distance squared
            if (projectForce) {
              r = Rij*moveDir;   // Take dot product of R along pulling dir.
            } else {
	      r2 = Rij.length2();
	      r = sqrt(r2);
            }

	    //  Only calculate the energy and the force if the distance is
	    //  non-zero.   Otherwise, you could end up dividing by 0, which
	    //  is bad
            Vector smdForce;
	    if (r>0.0)
	    {
	      value=k;
      
	      //  Loop through and multiple k by r consExp times.
	      //  i.e.  calculate kr^e/e
	      //  I know, I could use pow(), but I don't trust it.
	      for (int i=0; i<consExp; ++i)
	      {
		value *= r;
	      }

	      //  Add to the energy total
	      energy += value/consExp;
      
	      //  Now calculate the force, which is kr^(e-1).  Also, divide
	      //  by another factor of r to normalize the vector before we
	      //  multiple it by the magnitude
	      if (projectForce) {
                value /= r;
                smdForce = value*moveDir;
	      } else {
	        value /= r2;
	        Rij.mult(value);
	        smdForce = Rij;
              }
	    }
	    f[localID] += smdForce;
 
	    if (currentTime % outputFreq == 0) {
	      smdData->output(currentTime, atomPos, smdForce);
	    }
	    
	    break;  // change this when there are many atoms restrained.

	  } // if found atom
	} // loop

	reduction->item(REDUCTION_SMD_ENERGY) += energy;
	reduction->submit();

}
/*			END OF FUNCTION force				*/

