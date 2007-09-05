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
    Position loc;			// Fractional location within box
    Force f;
    float v;
    Force f2;
    float v2;
    Force f3;
    float v3;
    BigReal xphi[4], yphi[4], zphi[4];
    BigReal dxphi[4], dyphi[4], dzphi[4];
    BigReal t;
    
    Position center;
    Tensor inv;
    mol->get_gridfrc_info(center, inv);
    
    //  Loop through and check each atom
    for (int i = 0; i < numAtoms; i++) {
	if (mol->is_atom_gridforced(p[i].id)) {
	    mol->get_gridfrc_params(scale, charge, p[i].id);
	    
	    // Wrap coordinates using grid center
	    Position pos = p[i].position;
	    pos += homePatch->lattice.wrap_delta(p[i].position);
	    pos += homePatch->lattice.delta(pos, center) - (pos - center);
	    
	    int err = mol->get_gridfrc_grid(gbox, loc, pos);
	    mol->get_gridfrc_mgradv(p[i], f3, v3);
	    
	    if (err) {
		DebugM(4, "force = 0 0 0\n" << endi);
		//DebugM(4, "force2 = 0 0 0\n" << endi);
		//continue;  // This means the current atom is outside the potential
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
	    
	    // TESTING
	    // k1 direction
	    t = loc[0] + 1;
	    xphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t);
	    dxphi[0] = (1.5 * t - 2) * (2 - t);
	    t--;
	    xphi[1] = (1 - t) * (1 + t - 1.5 * t * t);
	    dxphi[1] = (-5 + 4.5 * t) * t;
	    t--;
	    xphi[2] = (1 + t) * (1 - t - 1.5 * t * t);
	    dxphi[2] = (-5 - 4.5 * t) * t;
	    t--;
	    xphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t);
	    dxphi[3] = (1.5 * t + 2) * (2 + t);
	    
	    // k2 direction
	    t = loc[1] + 1;
	    yphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t);
	    dyphi[0] = (1.5 * t - 2) * (2 - t);
	    t--;
	    yphi[1] = (1 - t) * (1 + t - 1.5 * t * t);
	    dyphi[1] = (-5 + 4.5 * t) * t;
	    t--;
	    yphi[2] = (1 + t) * (1 - t - 1.5 * t * t);
	    dyphi[2] = (-5 - 4.5 * t) * t;
	    t--;
	    yphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t);
	    dyphi[3] = (1.5 * t + 2) * (2 + t);
	    
	    // k3 direction
	    t = loc[2] + 1;
	    zphi[0] = 0.5 * (1 - t) * (2 - t) * (2 - t);
	    dzphi[0] = (1.5 * t - 2) * (2 - t);
	    t--;
	    zphi[1] = (1 - t) * (1 + t - 1.5 * t * t);
	    dzphi[1] = (-5 + 4.5 * t) * t;
	    t--;
	    zphi[2] = (1 + t) * (1 - t - 1.5 * t * t);
	    dzphi[2] = (-5 - 4.5 * t) * t;
	    t--;
	    zphi[3] = 0.5 * (1 + t) * (2 + t) * (2 + t);
	    dzphi[3] = (1.5 * t + 2) * (2 + t);
	    
	    f2 = 0;
	    v2 = 0;
	    for (int ii = 0; ii < 4; ii++) {
		for (int jj = 0; jj < 4; jj++) {
		    for(int kk = 0; kk < 4; kk++) {
			// Subtract since f_x = -pV/px
			f2.x -= gbox.v2[ii][jj][kk] * dxphi[ii] * yphi[jj] * zphi[kk];
			f2.y -= gbox.v2[ii][jj][kk] * xphi[ii] * dyphi[jj] * zphi[kk];
			f2.z -= gbox.v2[ii][jj][kk] * xphi[ii] * yphi[jj] * dzphi[kk];
			v2 += gbox.v2[ii][jj][kk] * xphi[ii] * yphi[jj] * zphi[kk];
		    }
		}
	    }
	    
	    Force force = scale * Tensor::diagonal(simParams->gridforceScale) * charge * (inv * f);
	    Force force2 = scale * Tensor::diagonal(simParams->gridforceScale) * charge * (inv * f2);
	    Force force3 = scale * Tensor::diagonal(simParams->gridforceScale) * charge * f3; // inv mult already done
	    
	    if (!err) {
		DebugM(4, "force = " << force << "\n" << endi);
		DebugM(4, "v = " << v << "\n" << endi);
	    }
	    DebugM(4, "force2 = " << force2 << "\n" << endi);
	    DebugM(4, "v2 = " << v2 << "\n" << endi);
	    DebugM(4, "force3 = " << force3 << "\n" << endi);
	    DebugM(4, "v3 = " << v3 << "\n" << endi);
	    
	    force = force3;
	    v = v3;
	    
	    forces[i] += force;
	    extForce += force;
	    Position vpos = homePatch->lattice.reverse_transform(
	      p[i].position, p[i].transform );
	    
	    DebugM(4, "transform = " << (int)p[i].transform.i << " "
		   << (int)p[i].transform.j << " " << (int)p[i].transform.k << "\n" << endi);
	    
	    //energy -= force * (vpos - homePatch->lattice.origin());
	    if (simParams->gridforceScale.x == simParams->gridforceScale.y 
		&& simParams->gridforceScale.x == simParams->gridforceScale.z)
	    {
		// only makes sense when scaling is isotropic
		energy += v * scale * simParams->gridforceScale.x;
	    }
	    extVirial += outer(force,vpos);
	}
    }
    
    reduction->item(REDUCTION_MISC_ENERGY) += energy;
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
    reduction->submit();
    
}
/*			END OF FUNCTION force				*/
