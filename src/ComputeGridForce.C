/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeGridForce.h"
#include "Node.h"
#include "SimParameters.h"
#include "HomePatch.h"
#include "Molecule.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"


ComputeGridForce::ComputeGridForce(ComputeID c, PatchID pid)
  : ComputeHomePatch(c,pid)
{

    reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

}
/*			END OF FUNCTION ComputeGridForce	*/


ComputeGridForce::~ComputeGridForce()

{
    delete reduction;
}
/*			END OF FUNCTION ~ComputeGridForce	*/


void ComputeGridForce::doForce(FullAtom* p, Results* r)
{
    SimParameters *simParams = Node::Object()->simParameters;
    Molecule *mol = Node::Object()->molecule;
    //Lattice lattice = homePatch->lattice;
    
    Force *forces = r->f[Results::normal];
    BigReal energy = 0;
    Force extForce = 0.;
    Tensor extVirial;
    
    Real scale;				// Scaling factor
    Charge charge;			// Charge
    Molecule::GridforceGridbox gbox;	// Structure with potential info
    Vector loc;				// Fractional location within box
    Vector inv;				// Inverse vector
    Force ftemp[2][2], f;
    float v;
    
    //  Loop through and check each atom
    for (int i = 0; i < numAtoms; i++) {
	if (mol->is_atom_gridforced(p[i].id)) {
	    mol->get_gridfrc_params(scale, charge, p[i].id);
	    
	    // Get wrapped position
	    Vector pos = p[i].position;
	    pos += homePatch->lattice.wrap_delta(pos);
	    
	    int err = mol->get_gridfrc_grid(gbox, loc, pos);
	    if (err) {
		//iout << "This atom is out!\n" << endi;
		continue;  // This means the current atom is outside the potential
	    }
	    
	    f = 0;
	    v = 0;
	    for (int j = 0; j < 8; j++) {
		// This is simply the shortest possible way of writing
		// this that I could think of, sorry if it's
		// confusing. (j & 4) is either nonzero or zero for
		// being on either the high or low k3 (high or low z)
		// side of the cube, respectively. Similarly, (j & 2)
		// talks about the k2 (y) direction, and (j & 1) the
		// k1 (x).  The factors are multiplied such that, if
		// we're very near a cube face, then the forces from
		// the opposite face hardly contribute to f. This way
		// the forces we apply are nice and continuous cube to
		// cube.
		
		float factor = ((j & 4) ? loc[0] : 1-loc[0])
		    * ((j & 2) ? loc[1] : 1-loc[1])
		    * ((j & 1) ? loc[2] : 1-loc[2]);
		
		f += factor * gbox.f[j];
		v += factor * gbox.v[j];
	    }
	    
	    Force force = charge * scale * simParams->gridforceScale;
	    for (int j = 0; j < 3; j++) force[j] *= f.dot(gbox.inv[j]);
	    
	    DebugM(4, "force = " << force[0] << " " <<  force[1] << " " <<  force[2] << "\n" << endi);
	    
	    forces[i] += force;
	    extForce += force;
	    Position vpos = homePatch->lattice.reverse_transform(
		p[i].position, p[i].transform );
	    //energy -= force * (vpos - homePatch->lattice.origin());
	    //energy += v * scale * simParams->gridforceScale.x		// only makes sense when scaling is isotropic
	    extVirial += outer(force,vpos);
	}
    }
    
    reduction->item(REDUCTION_MISC_ENERGY) += energy;
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
    reduction->submit();
    
}
/*			END OF FUNCTION force				*/
