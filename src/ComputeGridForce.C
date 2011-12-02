/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeGridForce.h"
#include "GridForceGrid.h"
#include "Node.h"
#include "SimParameters.h"
#include "HomePatch.h"
#include "Molecule.h"

#define MIN_DEBUG_LEVEL 2
//#define DEBUGM
#include "Debug.h"

#include "MGridforceParams.h"

#define GF_FORCE_OUTPUT
#define GF_FORCE_OUTPUT_FREQ 100


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
    
    Real scale;			// Scaling factor
    Charge charge;		// Charge
    //GridforceGrid::Box box;	// Structure with potential info
    Vector dV;
    float V;
    int numAtoms = homePatch->getNumAtoms();
    
    for (int gridnum = 0; gridnum < mol->numGridforceGrids; gridnum++) {
	const GridforceGrid *grid = mol->get_gridfrc_grid(gridnum);
	
// #ifdef DEBUGM
// 	if (homePatch->flags.step % 100 == 1) {
// 	    float* grid_flat = NULL;
// 	    int sz = grid->get_all_gridvals(&grid_flat);
// 	    for (int i = 0; i < sz; i++) {
// 		DebugM(4, "(step " << homePatch->flags.step << ") grid_flat[" << gridnum << "][" << i << "] = " << grid_flat[i] << "\n" << endi);
// 	    }
// 	    delete [] grid_flat;
// 	}
// #endif
	
	Position center = grid->get_center();
//	Tensor inv = grid->get_inv();
	Vector gfScale = grid->get_scale();
	
	DebugM(1, "center = " << center << "\n" << endi);
	
	//  Loop through and check each atom
	for (int i = 0; i < numAtoms; i++) {
//	    iout << iINFO << "Gridforced " << p[i].id << "--" << gridnum << "\n" << endi;
	    
	    if (mol->is_atom_gridforced(p[i].id, gridnum))
	    {
		DebugM(1, "Atom " << p[i].id << " is gridforced\n" << endi);
		
		mol->get_gridfrc_params(scale, charge, p[i].id, gridnum);
		
		// Wrap coordinates using grid center
		Position pos = p[i].position;
		pos += homePatch->lattice.wrap_delta(p[i].position);
		pos += homePatch->lattice.delta(pos, center) - (pos - center);
		
		DebugM(1, "pos = " << pos << "\n" << endi);
		
		// Here's where the action happens
		int err = grid->compute_VdV(pos, V, dV);
		if (err) {
		    DebugM(2, "V = 0\n" << endi);
		    DebugM(2, "dV = 0 0 0\n" << endi);
		    continue;  // This means the current atom is outside the potential
		}
		
		Force force = scale * Tensor::diagonal(gfScale) * (-charge * dV);
		
		DebugM(2, "scale = " << scale << " gfScale = " << gfScale << " charge = " << charge << "\n" << endi);
		
		DebugM(2, "V = " << V << "\n" << endi);
		DebugM(2, "dV = " << dV << "\n" << endi);
		DebugM(2, "grid = " << gridnum << " force = " << force << " pos = " << pos << " V = " << V << " dV = " << dV << " step = " << homePatch->flags.step << " index = " << p[i].id << "\n" << endi);
		
		if (V != V) {
		    iout << iWARN << "V is NaN!\natomid = " << p[i].id << " loc = " << p[i].position << " V = " << V << "\n" << endi;
		}
		
		//#ifdef GF_FORCE_OUTPUT
		//	    int t = patch->flags.step;
		//	    if (t % GF_FORCE_OUTPUT_FREQ == 0) {
		//		iout << "GRIDFORCE (TS ATOM FORCE)  " << t << ' ' << i << ' ' << force*PNPERKCALMOL << '\n' << endi;
		//	    }
		//#endif
	    
		forces[i] += force;
		extForce += force;
		Position vpos = homePatch->lattice.reverse_transform(p[i].position, p[i].transform);
		
		DebugM(1, "transform = " << (int)p[i].transform.i << " "
		       << (int)p[i].transform.j << " " << (int)p[i].transform.k << "\n" << endi);
		
		//energy -= force * (vpos - homePatch->lattice.origin());
		if (gfScale.x == gfScale.y && gfScale.x == gfScale.z)
		{
		    // only makes sense when scaling is isotropic
		    energy += scale * gfScale.x * (charge * V);
		    
		    // add something when we're off the grid? I'm thinking no
		}
		extVirial += outer(force,vpos);
	    }
	}
    }
//    for(int i=0; i < numAtoms; i++) {
//      iout << "ZZZAtom " << p[i].id << " Force " << forces[i].x << "," 
//           << forces[i].y << "," << forces[i].z
//           << endl;
//    }
    reduction->item(REDUCTION_MISC_ENERGY) += energy;
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
    reduction->submit();
}
/*			END OF FUNCTION force				*/
