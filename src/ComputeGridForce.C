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
	    if (err) {
		DebugM(4, "force4 = 0 0 0\n" << endi);
		DebugM(4, "v4 = 0\n" << endi);
		continue;  // This means the current atom is outside the potential
	    }
	    
	    float a[64];
// 	    for (int j = 0; j < 64; j++) {
// 		a[j] = 0;
// 		for (int k = 0; k < 64; k++) {
// 		    a[j] += A[j][k] * gbox.b[k];
// 		}
// 	    }
	    
	    // Multiply 'gbox.b' vector by matrix ... ugly but efficient (?)
	    a[0] = gbox.b[0];
	    a[1] = gbox.b[8];
	    a[2] = -3*gbox.b[0] + 3*gbox.b[1] - 2*gbox.b[8] - gbox.b[9];
	    a[3] = 2*gbox.b[0] - 2*gbox.b[1] + gbox.b[8] + gbox.b[9];
	    a[4] = gbox.b[16];
	    a[5] = gbox.b[32];
	    a[6] = -3*gbox.b[16] + 3*gbox.b[17] - 2*gbox.b[32] - gbox.b[33];
	    a[7] = 2*gbox.b[16] - 2*gbox.b[17] + gbox.b[32] + gbox.b[33];
	    a[8] = -3*gbox.b[0] + 3*gbox.b[2] - 2*gbox.b[16] - gbox.b[18];
	    a[9] = -3*gbox.b[8] + 3*gbox.b[10] - 2*gbox.b[32] - gbox.b[34];
	    a[10] = 9*gbox.b[0] - 9*gbox.b[1] - 9*gbox.b[2] + 9*gbox.b[3] + 6*gbox.b[8] + 3*gbox.b[9] - 6*gbox.b[10] - 3*gbox.b[11]
		+ 6*gbox.b[16] - 6*gbox.b[17] + 3*gbox.b[18] - 3*gbox.b[19] + 4*gbox.b[32] + 2*gbox.b[33] + 2*gbox.b[34] + gbox.b[35];
	    a[11] = -6*gbox.b[0] + 6*gbox.b[1] + 6*gbox.b[2] - 6*gbox.b[3] - 3*gbox.b[8] - 3*gbox.b[9] + 3*gbox.b[10] + 3*gbox.b[11]
		- 4*gbox.b[16] + 4*gbox.b[17] - 2*gbox.b[18] + 2*gbox.b[19] - 2*gbox.b[32] - 2*gbox.b[33] - gbox.b[34] - gbox.b[35];
	    a[12] = 2*gbox.b[0] - 2*gbox.b[2] + gbox.b[16] + gbox.b[18];
	    a[13] = 2*gbox.b[8] - 2*gbox.b[10] + gbox.b[32] + gbox.b[34];
	    a[14] = -6*gbox.b[0] + 6*gbox.b[1] + 6*gbox.b[2] - 6*gbox.b[3] - 4*gbox.b[8] - 2*gbox.b[9] + 4*gbox.b[10] + 2*gbox.b[11]
		- 3*gbox.b[16] + 3*gbox.b[17] - 3*gbox.b[18] + 3*gbox.b[19] - 2*gbox.b[32] - gbox.b[33] - 2*gbox.b[34] - gbox.b[35];
	    a[15] = 4*gbox.b[0] - 4*gbox.b[1] - 4*gbox.b[2] + 4*gbox.b[3] + 2*gbox.b[8] + 2*gbox.b[9] - 2*gbox.b[10] - 2*gbox.b[11]
		+ 2*gbox.b[16] - 2*gbox.b[17] + 2*gbox.b[18] - 2*gbox.b[19] + gbox.b[32] + gbox.b[33] + gbox.b[34] + gbox.b[35];
	    a[16] = gbox.b[24];
	    a[17] = gbox.b[40];
	    a[18] = -3*gbox.b[24] + 3*gbox.b[25] - 2*gbox.b[40] - gbox.b[41];
	    a[19] = 2*gbox.b[24] - 2*gbox.b[25] + gbox.b[40] + gbox.b[41];
	    a[20] = gbox.b[48];
	    a[21] = gbox.b[56];
	    a[22] = -3*gbox.b[48] + 3*gbox.b[49] - 2*gbox.b[56] - gbox.b[57];
	    a[23] = 2*gbox.b[48] - 2*gbox.b[49] + gbox.b[56] + gbox.b[57];
	    a[24] = -3*gbox.b[24] + 3*gbox.b[26] - 2*gbox.b[48] - gbox.b[50];
	    a[25] = -3*gbox.b[40] + 3*gbox.b[42] - 2*gbox.b[56] - gbox.b[58];
	    a[26] = 9*gbox.b[24] - 9*gbox.b[25] - 9*gbox.b[26] + 9*gbox.b[27] + 6*gbox.b[40] + 3*gbox.b[41] - 6*gbox.b[42] - 3*gbox.b[43]
		+ 6*gbox.b[48] - 6*gbox.b[49] + 3*gbox.b[50] - 3*gbox.b[51] + 4*gbox.b[56] + 2*gbox.b[57] + 2*gbox.b[58] + gbox.b[59];
	    a[27] = -6*gbox.b[24] + 6*gbox.b[25] + 6*gbox.b[26] - 6*gbox.b[27] - 3*gbox.b[40] - 3*gbox.b[41] + 3*gbox.b[42] + 3*gbox.b[43]
		- 4*gbox.b[48] + 4*gbox.b[49] - 2*gbox.b[50] + 2*gbox.b[51] - 2*gbox.b[56] - 2*gbox.b[57] - gbox.b[58] - gbox.b[59];
	    a[28] = 2*gbox.b[24] - 2*gbox.b[26] + gbox.b[48] + gbox.b[50];
	    a[29] = 2*gbox.b[40] - 2*gbox.b[42] + gbox.b[56] + gbox.b[58];
	    a[30] = -6*gbox.b[24] + 6*gbox.b[25] + 6*gbox.b[26] - 6*gbox.b[27] - 4*gbox.b[40] - 2*gbox.b[41] + 4*gbox.b[42] + 2*gbox.b[43]
		- 3*gbox.b[48] + 3*gbox.b[49] - 3*gbox.b[50] + 3*gbox.b[51] - 2*gbox.b[56] - gbox.b[57] - 2*gbox.b[58] - gbox.b[59];
	    a[31] = 4*gbox.b[24] - 4*gbox.b[25] - 4*gbox.b[26] + 4*gbox.b[27] + 2*gbox.b[40] + 2*gbox.b[41] - 2*gbox.b[42] - 2*gbox.b[43]
		+ 2*gbox.b[48] - 2*gbox.b[49] + 2*gbox.b[50] - 2*gbox.b[51] + gbox.b[56] + gbox.b[57] + gbox.b[58] + gbox.b[59];
	    a[32] = -3*gbox.b[0] + 3*gbox.b[4] - 2*gbox.b[24] - gbox.b[28];
	    a[33] = -3*gbox.b[8] + 3*gbox.b[12] - 2*gbox.b[40] - gbox.b[44];
	    a[34] = 9*gbox.b[0] - 9*gbox.b[1] - 9*gbox.b[4] + 9*gbox.b[5] + 6*gbox.b[8] + 3*gbox.b[9] - 6*gbox.b[12] - 3*gbox.b[13]
		+ 6*gbox.b[24] - 6*gbox.b[25] + 3*gbox.b[28] - 3*gbox.b[29] + 4*gbox.b[40] + 2*gbox.b[41] + 2*gbox.b[44] + gbox.b[45];
	    a[35] = -6*gbox.b[0] + 6*gbox.b[1] + 6*gbox.b[4] - 6*gbox.b[5] - 3*gbox.b[8] - 3*gbox.b[9] + 3*gbox.b[12] + 3*gbox.b[13]
		- 4*gbox.b[24] + 4*gbox.b[25] - 2*gbox.b[28] + 2*gbox.b[29] - 2*gbox.b[40] - 2*gbox.b[41] - gbox.b[44] - gbox.b[45];
	    a[36] = -3*gbox.b[16] + 3*gbox.b[20] - 2*gbox.b[48] - gbox.b[52];
	    a[37] = -3*gbox.b[32] + 3*gbox.b[36] - 2*gbox.b[56] - gbox.b[60];
	    a[38] = 9*gbox.b[16] - 9*gbox.b[17] - 9*gbox.b[20] + 9*gbox.b[21] + 6*gbox.b[32] + 3*gbox.b[33] - 6*gbox.b[36] - 3*gbox.b[37]
		+ 6*gbox.b[48] - 6*gbox.b[49] + 3*gbox.b[52] - 3*gbox.b[53] + 4*gbox.b[56] + 2*gbox.b[57] + 2*gbox.b[60] + gbox.b[61];
	    a[39] = -6*gbox.b[16] + 6*gbox.b[17] + 6*gbox.b[20] - 6*gbox.b[21] - 3*gbox.b[32] - 3*gbox.b[33] + 3*gbox.b[36] + 3*gbox.b[37]
		- 4*gbox.b[48] + 4*gbox.b[49] - 2*gbox.b[52] + 2*gbox.b[53] - 2*gbox.b[56] - 2*gbox.b[57] - gbox.b[60] - gbox.b[61];
	    a[40] = 9*gbox.b[0] - 9*gbox.b[2] - 9*gbox.b[4] + 9*gbox.b[6] + 6*gbox.b[16] + 3*gbox.b[18] - 6*gbox.b[20] - 3*gbox.b[22]
		+ 6*gbox.b[24] - 6*gbox.b[26] + 3*gbox.b[28] - 3*gbox.b[30] + 4*gbox.b[48] + 2*gbox.b[50] + 2*gbox.b[52] + gbox.b[54];
	    a[41] = 9*gbox.b[8] - 9*gbox.b[10] - 9*gbox.b[12] + 9*gbox.b[14] + 6*gbox.b[32] + 3*gbox.b[34] - 6*gbox.b[36] - 3*gbox.b[38]
		+ 6*gbox.b[40] - 6*gbox.b[42] + 3*gbox.b[44] - 3*gbox.b[46] + 4*gbox.b[56] + 2*gbox.b[58] + 2*gbox.b[60] + gbox.b[62];
	    a[42] = -27*gbox.b[0] + 27*gbox.b[1] + 27*gbox.b[2] - 27*gbox.b[3] + 27*gbox.b[4] - 27*gbox.b[5] - 27*gbox.b[6] + 27*gbox.b[7]
		- 18*gbox.b[8] - 9*gbox.b[9] + 18*gbox.b[10] + 9*gbox.b[11] + 18*gbox.b[12] + 9*gbox.b[13] - 18*gbox.b[14] - 9*gbox.b[15]
		- 18*gbox.b[16] + 18*gbox.b[17] - 9*gbox.b[18] + 9*gbox.b[19] + 18*gbox.b[20] - 18*gbox.b[21] + 9*gbox.b[22] - 9*gbox.b[23]
		- 18*gbox.b[24] + 18*gbox.b[25] + 18*gbox.b[26] - 18*gbox.b[27] - 9*gbox.b[28] + 9*gbox.b[29] + 9*gbox.b[30] - 9*gbox.b[31]
		- 12*gbox.b[32] - 6*gbox.b[33] - 6*gbox.b[34] - 3*gbox.b[35] + 12*gbox.b[36] + 6*gbox.b[37] + 6*gbox.b[38] + 3*gbox.b[39]
		- 12*gbox.b[40] - 6*gbox.b[41] + 12*gbox.b[42] + 6*gbox.b[43] - 6*gbox.b[44] - 3*gbox.b[45] + 6*gbox.b[46] + 3*gbox.b[47]
		- 12*gbox.b[48] + 12*gbox.b[49] - 6*gbox.b[50] + 6*gbox.b[51] - 6*gbox.b[52] + 6*gbox.b[53] - 3*gbox.b[54] + 3*gbox.b[55]
		- 8*gbox.b[56] - 4*gbox.b[57] - 4*gbox.b[58] - 2*gbox.b[59] - 4*gbox.b[60] - 2*gbox.b[61] - 2*gbox.b[62] - gbox.b[63];
	    a[43] = 18*gbox.b[0] - 18*gbox.b[1] - 18*gbox.b[2] + 18*gbox.b[3] - 18*gbox.b[4] + 18*gbox.b[5] + 18*gbox.b[6] - 18*gbox.b[7]
		+ 9*gbox.b[8] + 9*gbox.b[9] - 9*gbox.b[10] - 9*gbox.b[11] - 9*gbox.b[12] - 9*gbox.b[13] + 9*gbox.b[14] + 9*gbox.b[15]
		+ 12*gbox.b[16] - 12*gbox.b[17] + 6*gbox.b[18] - 6*gbox.b[19] - 12*gbox.b[20] + 12*gbox.b[21] - 6*gbox.b[22] + 6*gbox.b[23]
		+ 12*gbox.b[24] - 12*gbox.b[25] - 12*gbox.b[26] + 12*gbox.b[27] + 6*gbox.b[28] - 6*gbox.b[29] - 6*gbox.b[30] + 6*gbox.b[31]
		+ 6*gbox.b[32] + 6*gbox.b[33] + 3*gbox.b[34] + 3*gbox.b[35] - 6*gbox.b[36] - 6*gbox.b[37] - 3*gbox.b[38] - 3*gbox.b[39]
		+ 6*gbox.b[40] + 6*gbox.b[41] - 6*gbox.b[42] - 6*gbox.b[43] + 3*gbox.b[44] + 3*gbox.b[45] - 3*gbox.b[46] - 3*gbox.b[47]
		+ 8*gbox.b[48] - 8*gbox.b[49] + 4*gbox.b[50] - 4*gbox.b[51] + 4*gbox.b[52] - 4*gbox.b[53] + 2*gbox.b[54] - 2*gbox.b[55]
		+ 4*gbox.b[56] + 4*gbox.b[57] + 2*gbox.b[58] + 2*gbox.b[59] + 2*gbox.b[60] + 2*gbox.b[61] + gbox.b[62] + gbox.b[63];
	    a[44] = -6*gbox.b[0] + 6*gbox.b[2] + 6*gbox.b[4] - 6*gbox.b[6] - 3*gbox.b[16] - 3*gbox.b[18] + 3*gbox.b[20] + 3*gbox.b[22]
		- 4*gbox.b[24] + 4*gbox.b[26] - 2*gbox.b[28] + 2*gbox.b[30] - 2*gbox.b[48] - 2*gbox.b[50] - gbox.b[52] - gbox.b[54];
	    a[45] = -6*gbox.b[8] + 6*gbox.b[10] + 6*gbox.b[12] - 6*gbox.b[14] - 3*gbox.b[32] - 3*gbox.b[34] + 3*gbox.b[36] + 3*gbox.b[38]
		- 4*gbox.b[40] + 4*gbox.b[42] - 2*gbox.b[44] + 2*gbox.b[46] - 2*gbox.b[56] - 2*gbox.b[58] - gbox.b[60] - gbox.b[62];
	    a[46] = 18*gbox.b[0] - 18*gbox.b[1] - 18*gbox.b[2] + 18*gbox.b[3] - 18*gbox.b[4] + 18*gbox.b[5] + 18*gbox.b[6] - 18*gbox.b[7]
		+ 12*gbox.b[8] + 6*gbox.b[9] - 12*gbox.b[10] - 6*gbox.b[11] - 12*gbox.b[12] - 6*gbox.b[13] + 12*gbox.b[14] + 6*gbox.b[15]
		+ 9*gbox.b[16] - 9*gbox.b[17] + 9*gbox.b[18] - 9*gbox.b[19] - 9*gbox.b[20] + 9*gbox.b[21] - 9*gbox.b[22] + 9*gbox.b[23]
		+ 12*gbox.b[24] - 12*gbox.b[25] - 12*gbox.b[26] + 12*gbox.b[27] + 6*gbox.b[28] - 6*gbox.b[29] - 6*gbox.b[30] + 6*gbox.b[31]
		+ 6*gbox.b[32] + 3*gbox.b[33] + 6*gbox.b[34] + 3*gbox.b[35] - 6*gbox.b[36] - 3*gbox.b[37] - 6*gbox.b[38] - 3*gbox.b[39]
		+ 8*gbox.b[40] + 4*gbox.b[41] - 8*gbox.b[42] - 4*gbox.b[43] + 4*gbox.b[44] + 2*gbox.b[45] - 4*gbox.b[46] - 2*gbox.b[47]
		+ 6*gbox.b[48] - 6*gbox.b[49] + 6*gbox.b[50] - 6*gbox.b[51] + 3*gbox.b[52] - 3*gbox.b[53] + 3*gbox.b[54] - 3*gbox.b[55]
		+ 4*gbox.b[56] + 2*gbox.b[57] + 4*gbox.b[58] + 2*gbox.b[59] + 2*gbox.b[60] + gbox.b[61] + 2*gbox.b[62] + gbox.b[63];
	    a[47] = -12*gbox.b[0] + 12*gbox.b[1] + 12*gbox.b[2] - 12*gbox.b[3] + 12*gbox.b[4] - 12*gbox.b[5] - 12*gbox.b[6] + 12*gbox.b[7]
		- 6*gbox.b[8] - 6*gbox.b[9] + 6*gbox.b[10] + 6*gbox.b[11] + 6*gbox.b[12] + 6*gbox.b[13] - 6*gbox.b[14] - 6*gbox.b[15]
		- 6*gbox.b[16] + 6*gbox.b[17] - 6*gbox.b[18] + 6*gbox.b[19] + 6*gbox.b[20] - 6*gbox.b[21] + 6*gbox.b[22] - 6*gbox.b[23]
		- 8*gbox.b[24] + 8*gbox.b[25] + 8*gbox.b[26] - 8*gbox.b[27] - 4*gbox.b[28] + 4*gbox.b[29] + 4*gbox.b[30] - 4*gbox.b[31]
		- 3*gbox.b[32] - 3*gbox.b[33] - 3*gbox.b[34] - 3*gbox.b[35] + 3*gbox.b[36] + 3*gbox.b[37] + 3*gbox.b[38] + 3*gbox.b[39]
		- 4*gbox.b[40] - 4*gbox.b[41] + 4*gbox.b[42] + 4*gbox.b[43] - 2*gbox.b[44] - 2*gbox.b[45] + 2*gbox.b[46] + 2*gbox.b[47]
		- 4*gbox.b[48] + 4*gbox.b[49] - 4*gbox.b[50] + 4*gbox.b[51] - 2*gbox.b[52] + 2*gbox.b[53] - 2*gbox.b[54] + 2*gbox.b[55]
		- 2*gbox.b[56] - 2*gbox.b[57] - 2*gbox.b[58] - 2*gbox.b[59] - gbox.b[60] - gbox.b[61] - gbox.b[62] - gbox.b[63];
	    a[48] = 2*gbox.b[0] - 2*gbox.b[4] + gbox.b[24] + gbox.b[28];
	    a[49] = 2*gbox.b[8] - 2*gbox.b[12] + gbox.b[40] + gbox.b[44];
	    a[50] = -6*gbox.b[0] + 6*gbox.b[1] + 6*gbox.b[4] - 6*gbox.b[5] - 4*gbox.b[8] - 2*gbox.b[9] + 4*gbox.b[12] + 2*gbox.b[13]
		- 3*gbox.b[24] + 3*gbox.b[25] - 3*gbox.b[28] + 3*gbox.b[29] - 2*gbox.b[40] - gbox.b[41] - 2*gbox.b[44] - gbox.b[45];
	    a[51] = 4*gbox.b[0] - 4*gbox.b[1] - 4*gbox.b[4] + 4*gbox.b[5] + 2*gbox.b[8] + 2*gbox.b[9] - 2*gbox.b[12] - 2*gbox.b[13]
		+ 2*gbox.b[24] - 2*gbox.b[25] + 2*gbox.b[28] - 2*gbox.b[29] + gbox.b[40] + gbox.b[41] + gbox.b[44] + gbox.b[45];
	    a[52] = 2*gbox.b[16] - 2*gbox.b[20] + gbox.b[48] + gbox.b[52];
	    a[53] = 2*gbox.b[32] - 2*gbox.b[36] + gbox.b[56] + gbox.b[60];
	    a[54] = -6*gbox.b[16] + 6*gbox.b[17] + 6*gbox.b[20] - 6*gbox.b[21] - 4*gbox.b[32] - 2*gbox.b[33] + 4*gbox.b[36] + 2*gbox.b[37]
		- 3*gbox.b[48] + 3*gbox.b[49] - 3*gbox.b[52] + 3*gbox.b[53] - 2*gbox.b[56] - gbox.b[57] - 2*gbox.b[60] - gbox.b[61];
	    a[55] = 4*gbox.b[16] - 4*gbox.b[17] - 4*gbox.b[20] + 4*gbox.b[21] + 2*gbox.b[32] + 2*gbox.b[33] - 2*gbox.b[36] - 2*gbox.b[37]
		+ 2*gbox.b[48] - 2*gbox.b[49] + 2*gbox.b[52] - 2*gbox.b[53] + gbox.b[56] + gbox.b[57] + gbox.b[60] + gbox.b[61];
	    a[56] = -6*gbox.b[0] + 6*gbox.b[2] + 6*gbox.b[4] - 6*gbox.b[6] - 4*gbox.b[16] - 2*gbox.b[18] + 4*gbox.b[20] + 2*gbox.b[22]
		- 3*gbox.b[24] + 3*gbox.b[26] - 3*gbox.b[28] + 3*gbox.b[30] - 2*gbox.b[48] - gbox.b[50] - 2*gbox.b[52] - gbox.b[54];
	    a[57] = -6*gbox.b[8] + 6*gbox.b[10] + 6*gbox.b[12] - 6*gbox.b[14] - 4*gbox.b[32] - 2*gbox.b[34] + 4*gbox.b[36] + 2*gbox.b[38]
		- 3*gbox.b[40] + 3*gbox.b[42] - 3*gbox.b[44] + 3*gbox.b[46] - 2*gbox.b[56] - gbox.b[58] - 2*gbox.b[60] - gbox.b[62];
	    a[58] = 18*gbox.b[0] - 18*gbox.b[1] - 18*gbox.b[2] + 18*gbox.b[3] - 18*gbox.b[4] + 18*gbox.b[5] + 18*gbox.b[6] - 18*gbox.b[7]
		+ 12*gbox.b[8] + 6*gbox.b[9] - 12*gbox.b[10] - 6*gbox.b[11] - 12*gbox.b[12] - 6*gbox.b[13] + 12*gbox.b[14] + 6*gbox.b[15]
		+ 12*gbox.b[16] - 12*gbox.b[17] + 6*gbox.b[18] - 6*gbox.b[19] - 12*gbox.b[20] + 12*gbox.b[21] - 6*gbox.b[22] + 6*gbox.b[23]
		+ 9*gbox.b[24] - 9*gbox.b[25] - 9*gbox.b[26] + 9*gbox.b[27] + 9*gbox.b[28] - 9*gbox.b[29] - 9*gbox.b[30] + 9*gbox.b[31]
		+ 8*gbox.b[32] + 4*gbox.b[33] + 4*gbox.b[34] + 2*gbox.b[35] - 8*gbox.b[36] - 4*gbox.b[37] - 4*gbox.b[38] - 2*gbox.b[39]
		+ 6*gbox.b[40] + 3*gbox.b[41] - 6*gbox.b[42] - 3*gbox.b[43] + 6*gbox.b[44] + 3*gbox.b[45] - 6*gbox.b[46] - 3*gbox.b[47]
		+ 6*gbox.b[48] - 6*gbox.b[49] + 3*gbox.b[50] - 3*gbox.b[51] + 6*gbox.b[52] - 6*gbox.b[53] + 3*gbox.b[54] - 3*gbox.b[55]
		+ 4*gbox.b[56] + 2*gbox.b[57] + 2*gbox.b[58] + gbox.b[59] + 4*gbox.b[60] + 2*gbox.b[61] + 2*gbox.b[62] + gbox.b[63];
	    a[59] = -12*gbox.b[0] + 12*gbox.b[1] + 12*gbox.b[2] - 12*gbox.b[3] + 12*gbox.b[4] - 12*gbox.b[5] - 12*gbox.b[6] + 12*gbox.b[7]
		- 6*gbox.b[8] - 6*gbox.b[9] + 6*gbox.b[10] + 6*gbox.b[11] + 6*gbox.b[12] + 6*gbox.b[13] - 6*gbox.b[14] - 6*gbox.b[15]
		- 8*gbox.b[16] + 8*gbox.b[17] - 4*gbox.b[18] + 4*gbox.b[19] + 8*gbox.b[20] - 8*gbox.b[21] + 4*gbox.b[22] - 4*gbox.b[23]
		- 6*gbox.b[24] + 6*gbox.b[25] + 6*gbox.b[26] - 6*gbox.b[27] - 6*gbox.b[28] + 6*gbox.b[29] + 6*gbox.b[30] - 6*gbox.b[31]
		- 4*gbox.b[32] - 4*gbox.b[33] - 2*gbox.b[34] - 2*gbox.b[35] + 4*gbox.b[36] + 4*gbox.b[37] + 2*gbox.b[38] + 2*gbox.b[39]
		- 3*gbox.b[40] - 3*gbox.b[41] + 3*gbox.b[42] + 3*gbox.b[43] - 3*gbox.b[44] - 3*gbox.b[45] + 3*gbox.b[46] + 3*gbox.b[47]
		- 4*gbox.b[48] + 4*gbox.b[49] - 2*gbox.b[50] + 2*gbox.b[51] - 4*gbox.b[52] + 4*gbox.b[53] - 2*gbox.b[54] + 2*gbox.b[55]
		- 2*gbox.b[56] - 2*gbox.b[57] - gbox.b[58] - gbox.b[59] - 2*gbox.b[60] - 2*gbox.b[61] - gbox.b[62] - gbox.b[63];
	    a[60] = 4*gbox.b[0] - 4*gbox.b[2] - 4*gbox.b[4] + 4*gbox.b[6] + 2*gbox.b[16] + 2*gbox.b[18] - 2*gbox.b[20] - 2*gbox.b[22]
		+ 2*gbox.b[24] - 2*gbox.b[26] + 2*gbox.b[28] - 2*gbox.b[30] + gbox.b[48] + gbox.b[50] + gbox.b[52] + gbox.b[54];
	    a[61] = 4*gbox.b[8] - 4*gbox.b[10] - 4*gbox.b[12] + 4*gbox.b[14] + 2*gbox.b[32] + 2*gbox.b[34] - 2*gbox.b[36] - 2*gbox.b[38]
		+ 2*gbox.b[40] - 2*gbox.b[42] + 2*gbox.b[44] - 2*gbox.b[46] + gbox.b[56] + gbox.b[58] + gbox.b[60] + gbox.b[62];
	    a[62] = -12*gbox.b[0] + 12*gbox.b[1] + 12*gbox.b[2] - 12*gbox.b[3] + 12*gbox.b[4] - 12*gbox.b[5] - 12*gbox.b[6] + 12*gbox.b[7]
		- 8*gbox.b[8] - 4*gbox.b[9] + 8*gbox.b[10] + 4*gbox.b[11] + 8*gbox.b[12] + 4*gbox.b[13] - 8*gbox.b[14] - 4*gbox.b[15]
		- 6*gbox.b[16] + 6*gbox.b[17] - 6*gbox.b[18] + 6*gbox.b[19] + 6*gbox.b[20] - 6*gbox.b[21] + 6*gbox.b[22] - 6*gbox.b[23]
		- 6*gbox.b[24] + 6*gbox.b[25] + 6*gbox.b[26] - 6*gbox.b[27] - 6*gbox.b[28] + 6*gbox.b[29] + 6*gbox.b[30] - 6*gbox.b[31]
		- 4*gbox.b[32] - 2*gbox.b[33] - 4*gbox.b[34] - 2*gbox.b[35] + 4*gbox.b[36] + 2*gbox.b[37] + 4*gbox.b[38] + 2*gbox.b[39]
		- 4*gbox.b[40] - 2*gbox.b[41] + 4*gbox.b[42] + 2*gbox.b[43] - 4*gbox.b[44] - 2*gbox.b[45] + 4*gbox.b[46] + 2*gbox.b[47]
		- 3*gbox.b[48] + 3*gbox.b[49] - 3*gbox.b[50] + 3*gbox.b[51] - 3*gbox.b[52] + 3*gbox.b[53] - 3*gbox.b[54] + 3*gbox.b[55]
		- 2*gbox.b[56] - gbox.b[57] - 2*gbox.b[58] - gbox.b[59] - 2*gbox.b[60] - gbox.b[61] - 2*gbox.b[62] - gbox.b[63];
	    a[63] = 8*gbox.b[0] - 8*gbox.b[1] - 8*gbox.b[2] + 8*gbox.b[3] - 8*gbox.b[4] + 8*gbox.b[5] + 8*gbox.b[6] - 8*gbox.b[7]
		+ 4*gbox.b[8] + 4*gbox.b[9] - 4*gbox.b[10] - 4*gbox.b[11] - 4*gbox.b[12] - 4*gbox.b[13] + 4*gbox.b[14] + 4*gbox.b[15]
		+ 4*gbox.b[16] - 4*gbox.b[17] + 4*gbox.b[18] - 4*gbox.b[19] - 4*gbox.b[20] + 4*gbox.b[21] - 4*gbox.b[22] + 4*gbox.b[23]
		+ 4*gbox.b[24] - 4*gbox.b[25] - 4*gbox.b[26] + 4*gbox.b[27] + 4*gbox.b[28] - 4*gbox.b[29] - 4*gbox.b[30] + 4*gbox.b[31]
		+ 2*gbox.b[32] + 2*gbox.b[33] + 2*gbox.b[34] + 2*gbox.b[35] - 2*gbox.b[36] - 2*gbox.b[37] - 2*gbox.b[38] - 2*gbox.b[39]
		+ 2*gbox.b[40] + 2*gbox.b[41] - 2*gbox.b[42] - 2*gbox.b[43] + 2*gbox.b[44] + 2*gbox.b[45] - 2*gbox.b[46] - 2*gbox.b[47]
		+ 2*gbox.b[48] - 2*gbox.b[49] + 2*gbox.b[50] - 2*gbox.b[51] + 2*gbox.b[52] - 2*gbox.b[53] + 2*gbox.b[54] - 2*gbox.b[55]
		+ gbox.b[56] + gbox.b[57] + gbox.b[58] + gbox.b[59] + gbox.b[60] + gbox.b[61] + gbox.b[62] + gbox.b[63];
	    
	    //for (int j = 0; j < 64; j++) DebugM(4, "a[" << j << "] = " << a[j] << "\n" << endi);
	    //for (int j = 0; j < 64; j++) DebugM(4, "b[" << j << "] = " << gbox.b[j] << "\n" << endi);
	    
	    // Calculate powers of x, y, z for later use
	    // e.g. x[2] = x^2
	    float x[4], y[4], z[4];
	    x[0] = 1; y[0] = 1; z[0] = 1;
	    for (int j = 1; j < 4; j++) {
		x[j] = x[j-1] * loc.x;
		y[j] = y[j-1] * loc.y;
		z[j] = z[j-1] * loc.z;
	    }
	    
	    int ind = 0;
	    f = 0;
	    v = 0;
	    for (int l = 0; l < 4; l++) {
		for (int k = 0; k < 4; k++) {
		    for (int j = 0; j < 4; j++) {
			v += a[ind] * x[j] * y[k] * z[l];
			if (j > 0) f.x -= a[ind] * j * x[j-1] * y[k]   * z[l];
			if (k > 0) f.y -= a[ind] * k * x[j]   * y[k-1] * z[l];
			if (l > 0) f.z -= a[ind] * l * x[j]   * y[k]   * z[l-1];
			ind++;
		    }
		}
	    }
	    
// 	    Force force = scale * Tensor::diagonal(simParams->gridforceScale) * charge * (inv * f);
	    Force force = scale * Tensor::diagonal(simParams->gridforceScale) * charge * (f * inv);	// Must multiply ON THE RIGHT by inv tensor
	    
	    DebugM(4, "f4 = " << f << "\n" << endi);
	    DebugM(4, "force4 = " << force << "\n" << endi);
	    DebugM(4, "v4 = " << v << "\n" << endi);
	    
	    forces[i] += force;
	    extForce += force;
	    Position vpos = homePatch->lattice.reverse_transform(p[i].position, p[i].transform);
	    
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
