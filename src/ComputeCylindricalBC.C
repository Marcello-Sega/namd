/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/
       
#include "ComputeCylindricalBC.h"
#include "Node.h"
#include "SimParameters.h"

/************************************************************************/
/*									*/
/*			FUNCTION ComputeCylindricalBC			*/
/*									*/
/*	This is the constructor for the ComputeCylindricalBC force object.*/
/*   It is responsible for getting all the parameters from the 		*/
/*   SimParameters object and then determining if this object needs	*/
/*   to perform any computation.  It only needs to do so if there is	*/
/*   some portion of the patch that lays outside of the cylindrical	*/
/*   boundaries.							*/
/*									*/
/************************************************************************/

ComputeCylindricalBC::ComputeCylindricalBC(ComputeID c, PatchID pid)
  : ComputePatch(c,pid)
{
	reduction = ReductionMgr::Object();
	reduction->Register(REDUCTION_BC_ENERGY);
	fake_seq = 0;

	SimParameters *simParams = Node::Object()->simParameters;

	//  Get parameters from the SimParameters object
	r1 = simParams->cylindricalBCr1;
	r2 = simParams->cylindricalBCr2;
	r1_2 = r1*r1;
	r2_2 = r2*r2;
	k1 = simParams->cylindricalBCk1;
	k2 = simParams->cylindricalBCk2;
	exp1 = simParams->cylindricalBCexp1;
	exp2 = simParams->cylindricalBCexp2;
	//Additions for ends
	l1 = simParams->cylindricalBCl1;
	l2 = simParams->cylindricalBCl2;
	l1_2 = l1*l1;
	l2_2 = l2*l2;

	//  Check to see if this is one set of parameters or two
	if (r2 > -1.0)
	{
		twoForces = TRUE;
	}
	else
	{
		twoForces = FALSE;
	}

	//  Get the center of the cylinder, either the center of mass or
	//  a user-defined center
						/* I will want to have
						   a user defined center */
	if (simParams->cylindricalCenterCOM)
	{
		NAMD_die("Sorry, can't center about center of mass yet.\n");
		// center = namdMyNode->com;
	}
	else
	{
		center = simParams->cylindricalCenter;
	}

// Removed work needed checks for simplicity -JCP
//	doAnything = TRUE;
//	doLateral = TRUE;
//	doFaces = TRUE;

/*
	Vector comparePoint;	//  Point on patch boundary to compare
	Vector midPoint;	//  Middle of the patch
	Vector diff;		//  Difference between compare point
				//  and center of cylinder
	BigReal dist2;		//  Distance between point and center
				//  of cylinder squared
	BigReal dist3;		//  Distance between point and center

	//  Now, determine what point on the patch is furthest to
	//  the center of the cylinder.  This is done using the following
	//  algorithm.  For each dimension, compare the midpoint of the
	//  patch to the center of the cylinder.  If the center of the
	//  cylinder is less than the midpoint, take the origin plus the
	//  patch size as the cooridinate in that dimension to check.
	//  If it is greater than the midpoint, choose the origin's
	//  coordinate in this dimension as the cooridinate to use.

	//  Determine the midpoint
	//  size is dimension of patch.  Apparently patches are square.
	//  Are square patches necessarily optimal for my example
	//  of an elongated cylinder?
	midPoint.x = origin->x+0.5*size;
	midPoint.y = origin->y+0.5*size;
	midPoint.z = origin->z+0.5*size;

	//  One only needs to check y  and z dimensions.
	//  Z-dimension
	if (center.z < midPoint.z)
	{
		comparePoint.z = origin->z+size;
	}
	else 
	{
		comparePoint.z = origin->z;
	}

	//  Y-dimension
	if (center.y < midPoint.y)
	{
		comparePoint.y = origin->y+size;
	}
	else 
	{
		comparePoint.y = origin->y;
	}

	//  X-dimension
		comparePoint.x= center.x;


	//  Now, find a vector from the center to the comparison point
	diff = center-comparePoint;

	//  Calculate the distance squared
	dist2 = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;

        // Now perform computations not for the lateral sides, but
        // for the faces.  This involves the x-direction only.
        if (center.x < midPoint.x)
        {
                comparePoint.x = origin->x+size;
        }
        else
        {
                comparePoint.x = origin->x;
        }
        dist3 = fabs (comparePoint.x-center.x);

	//  Compare it to the radius values and determine if there
	//  is really any work for this patch to do
	if (dist2 > r1_2)
	{
		//  Outside of radius 1, we really have work to do
		doAnything = TRUE;
		doLateral = TRUE;
	}
	else if ( twoForces && (dist2 > r2_2) )
	{
		//  Two sets of parameters, and we are inside of
		//  the second radius, we really have work to do
		doAnything = TRUE;
		doLateral = TRUE;
	}
	else
	{
		//  Inside of all active radii, do nothing
		doAnything = FALSE;
		doLateral = FALSE;
	}

        // Handle cases for the ends
        if (dist3 > l1)
        {
                // beyond cylinder length 1, work to do
                doAnything = TRUE;
                doFaces = TRUE;
        }
        else if ( twoForces && (dist2 > l2))
        {
                // Two sets of parameters, work to do
                doAnything = TRUE;
                doFaces = TRUE;
        }
	else
	{
	doFaces = FALSE;
	}
*/


}
/*			END OF FUNCTION ComputeCylindricalBC		*/

/************************************************************************/
/*									*/
/*			FUNCTION ~ComputeCylindricalBC			*/
/*									*/
/*	This is the destructor for the ComputeCylindricalBC force object.	*/
/*   It currently does ABSOLUTELY NOTHING!!				*/
/*									*/
/************************************************************************/

ComputeCylindricalBC::~ComputeCylindricalBC()

{
	reduction->unRegister(REDUCTION_BC_ENERGY);
}
/*			END OF FUNCTION ~ComputeCylindricalBC		*/

/************************************************************************/
/*									*/
/*				FUNCTION force				*/
/*									*/
/*   INPUTS:								*/
/*	numAtoms - Number of coordinates being passed			*/
/*	x - Array of atom coordinates					*/
/*	forces - Array of force vectors					*/
/*									*/
/*	This function calculates the force and energy for the cylindri. */
/*   boundary conditions for this patch.				*/
/*									*/
/************************************************************************/

void ComputeCylindricalBC::doForce(Position* p, Results* r, AtomProperties* a)

{
	Vector diff;		//  Distance from atom to center of cylinder
	Vector f;		//  Calculated force vector
	int i, j;		//  Loop counters
	BigReal dist, dist_2;	//  Distance from atom to center, and this
				//  distance squared
	BigReal rval;		//  Difference between distance from atom
				//  to center and radius of cylinder
	BigReal eval;		//  Energy value for this atom
	BigReal fval;		//  Force magnitude for this atom

	// aliases to work with old code
	Position *x = p;
	Force *forces = r->f[Results::normal];
	BigReal energy = 0;

	// There are a couple of possibilities.  We could require
	// boundary conditions to be applied across the lateral
	// area, or across the faces, or both.
//        if (doLateral)
        {
	//  Loop through and check each atom
	for (i=0; i<numAtoms; i++)
	{
		//  Calculate the vector from the atom to the center of the
		//  cylinder
		diff.z = x[i].z - center.z;
		diff.y = x[i].y - center.y;
		diff.x = 0.0;
		
		//  Calculate the distance squared
		dist_2 = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;

		//  Look to see if we are outside either radius
		if ( (dist_2 > r1_2) || (twoForces && (dist_2 > r2_2)) )
		{
			//  Calculate the distance to the center
			dist = sqrt(dist_2);

			//  Normalize the direction vector
			diff.div(dist);

			//  Check to see if we are outside radius 1
			if (dist > r1)
			{
//printf ("In first force dist=%f\n", dist);
				//  Assign the force vector to the
				//  unit direction vector
				f.x = diff.x;
				f.y = diff.y;
				f.z = diff.z;

				//  Calculate the energy which is
				//  e = k1*(r_i-r_center)^exp1
				eval = k1;
				rval = fabs(dist - r1);

				for (j=0; j<exp1; j++)
				{
					eval *= rval;
				}

				energy += eval;

				//  Now calculate the force which is
				//  e = -k1*exp1*(r_i-r_center1)^(exp1-1)
				fval = -exp1*k1;

				for (j=0; j<exp1-1; j++)
				{
					fval *= rval;
				}

				//  Multiply the force magnitude to the
				//  unit direction vector to get the
				//  resulting force
				f.mult(fval);

				//  Add the force to the force vectors
				forces[i].x += f.x;
				forces[i].y += f.y;
				forces[i].z += f.z;
			}

			//  Check to see if there is a second radius
			//  and if we are outside of it
			if (twoForces && (dist > r2) )
			{
//printf ("In two forces dist = %f\n", dist);
				//  Assign the force vector to the
				//  unit direction vector
				f.x = diff.x;
				f.y = diff.y;
				f.z = diff.z;

				//  Calculate the energy which is
				//  e = k2*(r_i-r_center2)^exp2
				eval = k2;
				rval = fabs(dist - r2);

				for (j=0; j<exp2; j++)
				{
					eval *= rval;
				}

				energy += eval;

				//  Now calculate the force which is
				//  e = -k2*exp2*(r_i-r_center2)^(exp2-1)
				fval = -exp2*k2;

				for (j=0; j<exp2-1; j++)
				{
					fval *= rval;
				}

				//  Multiply the force magnitude to the
				//  unit direction vector to get the
				//  resulting force
				f.mult(fval);

				//  Add the force to the force vectors
				forces[i].x += f.x;
				forces[i].y += f.y;
				forces[i].z += f.z;
			}
		}
	}
  }    // End lateral condition








// Handle ends
// if (doFaces)
     {
       //  Loop through and check each atom
        for (i=0; i<numAtoms; i++)
        {
                //  Calculate the vector from the atom to the center of the
                //  cylinder
		diff.z = 0.0;
		diff.y = 0.0;
                diff.x = x[i].x - center.x;

                //  Calculate the distance squared
                dist_2 = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;

                //  Look to see if we are outside either radius
                if ( (dist_2 > l1_2)  || (twoForces && (dist_2 > l2_2)) )
                {
                        //  Calculate the distance to the center
                        dist = sqrt(dist_2);

                        //  Normalize the direction vector
                        diff.div(dist);

                        //  Check to see if we are outside radius 1
                        if (dist > l1)
                        {
//printf ("Do faces one force z=%f\n", dist);
                                //  Assign the force vector to the
                                //  unit direction vector
                                f.x = diff.x;
                                f.y = diff.y;
                                f.z = diff.z;

                                //  Calculate the energy which is
                                //  e = k1*(r_i-r_center)^exp1
                                eval = k1;
                                rval = fabs(dist - l1);

                                for (j=0; j<exp1; j++)
                                {
                                        eval *= rval;
                                }

                                energy += eval;

                                //  Now calculate the force which is
                                //  e = -k1*exp1*(r_i-r_center1)^(exp1-1)
                                fval = -exp1*k1;

                                for (j=0; j<exp1-1; j++)
                                {
                                        fval *= rval;
                                }

                                //  Multiply the force magnitude to the
                                //  unit direction vector to get the
                                //  resulting force
                                f.mult(fval);

                                //  Add the force to the force vectors
                                forces[i].x += f.x;
                                forces[i].y += f.y;
                                forces[i].z += f.z;
                        }

                        //  Check to see if there is a second radius
                        //  and if we are outside of it
                        if (twoForces && (dist > l2) )
                        {
//printf ("Do faces two force z=%f\n", dist);
                                //  Assign the force vector to the
                                //  unit direction vector
                                f.x = diff.x;
                                f.y = diff.y;
                                f.z = diff.z;

                                //  Calculate the energy which is
                                //  e = k2*(r_i-r_center2)^exp2
                                eval = k2;
                                rval = fabs(dist - l2);

                                for (j=0; j<exp2; j++)
                                {
                                        eval *= rval;
                                }

                                energy += eval;

                                //  Now calculate the force which is
                                //  e = -k2*exp2*(r_i-r_center2)^(exp2-1)
                                fval = -exp2*k2;

                                for (j=0; j<exp2-1; j++)
                                {
                                        fval *= rval;
                                }

                                //  Multiply the force magnitude to the
                                //  unit direction vector to get the
                                //  resulting force
                                f.mult(fval);

                                //  Add the force to the force vectors
                                forces[i].x += f.x;
                                forces[i].y += f.y;
                                forces[i].z += f.z;
                        }
                }
        }

    reduction->submit(fake_seq, REDUCTION_BC_ENERGY, energy);

    }
// END Additions for ends



}
/*			END OF FUNCTION force				*/
