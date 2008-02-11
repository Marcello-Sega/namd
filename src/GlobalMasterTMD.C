/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "Molecule.h"
#include "NamdTypes.h"

#include "SimParameters.h"
#include "GlobalMasterTMD.h"
#include "PDB.h"
#include "PDBData.h"
#include "fitrms.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"


class Matrix4 {
  BigReal mat[16];                               ///< the matrix itself
public:
  Matrix4(void) { identity(); }
  Matrix4(const BigReal *m)  { memcpy(mat, m, 16*sizeof(BigReal)); }
  void multpoint(BigReal point[3]) const {
    BigReal tmp[3];
    BigReal itmp3 = 1.0f / (point[0]*mat[3] + point[1]*mat[7] +
                            point[2]*mat[11] + mat[15]);
    tmp[0] = itmp3*point[0];
    tmp[1] = itmp3*point[1];
    tmp[2] = itmp3*point[2];
    point[0]=tmp[0]*mat[0] + tmp[1]*mat[4] + tmp[2]*mat[ 8] + itmp3*mat[12];
    point[1]=tmp[0]*mat[1] + tmp[1]*mat[5] + tmp[2]*mat[ 9] + itmp3*mat[13];
    point[2]=tmp[0]*mat[2] + tmp[1]*mat[6] + tmp[2]*mat[10] + itmp3*mat[14];
  }

  void identity() {
    memset(mat, 0, 16*sizeof(BigReal));
    mat[0]=1.0f;
    mat[5]=1.0f;
    mat[10]=1.0f;
    mat[15]=1.0f;
  }
  void transpose() {
    BigReal tmp[16];
    int i,j;
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) {
        tmp[4*i+j] = mat[i+4*j];
      }
    }
    for(i=0;i<16;i++) mat[i] = tmp[i];
  }
  /// premultiply the matrix by the given matrix, this->other * this
  void multmatrix(const Matrix4 &m) {
    BigReal tmp[4];
    for (int j=0; j<4; j++) {
      tmp[0] = mat[j];
      tmp[1] = mat[4+j];
      tmp[2] = mat[8+j]; 
      tmp[3] = mat[12+j];
      for (int i=0; i<4; i++) {
        mat[4*i+j] = m.mat[4*i]*tmp[0] + m.mat[4*i+1]*tmp[1] +
          m.mat[4*i+2]*tmp[2] + m.mat[4*i+3]*tmp[3]; 
      }
    } 
  }
  void translate(BigReal x, BigReal y, BigReal z) {
    Matrix4 m;		
    m.mat[12] = x;
    m.mat[13] = y;
    m.mat[14] = z;
    multmatrix(m);
  }
  void translate(BigReal d[3]) { translate(d[0], d[1], d[2]); }
};

GlobalMasterTMD::GlobalMasterTMD() {
  DebugM(3,"initialize called\n");
  SimParameters *params = Node::Object()->simParameters;
  outputFreq = params->TMDOutputFreq;
  k = params->TMDk;
  if (!params->TMDInitialRMSD)
    initialRMS = -1; // get from first coordinates
  else
    initialRMS = params->TMDInitialRMSD;
  finalRMS = params->TMDFinalRMSD;

  currentStep = params->firstTimestep;
  firstStep = params->TMDFirstStep;
  lastStep = params->TMDLastStep;

  parseAtoms(params->TMDFile,Node::Object()->molecule->numAtoms);
  k /= numTMDatoms;
  iout << iINFO << numTMDatoms << " TMD ATOMS\n" << endi;
  DebugM(1,"done with initialize\n");
}

void GlobalMasterTMD::parseAtoms(const char *file, int numTotalAtoms) {
  DebugM(3,"parseAtoms called\n");
  PDB tmdpdb(file);
  int numatoms = tmdpdb.num_atoms();
  if (numatoms < 1) 
    NAMD_die("No atoms found in TMDFile\n");
  if (numatoms != numTotalAtoms)
    NAMD_die("The number of atoms in TMDFile must be equal to the total number of atoms in the structure!");
  if ( modifyRequestedAtoms().size() )
    NAMD_bug("GlobalMasterTMD::parseAtoms() modifyRequestedAtoms() not empty");

  numTMDatoms = 0;
  target = new BigReal[3*numatoms];
  aidmap = new int[numatoms];
  Vector *atompos = new Vector[numatoms];
  tmdpdb.get_all_positions(atompos);
  int i;
  for (i=0; i<numatoms; i++) {
#ifdef MEM_OPT_VERSION
    PDBCoreData *atom = tmdpdb.atom(i);
#else
    PDBAtom *atom = tmdpdb.atom(i); // get an atom from the file
#endif
    if (atom->occupancy()) { // if occupancy is not 0, then add it!
      target[3*numTMDatoms  ] = atompos[i].x;
      target[3*numTMDatoms+1] = atompos[i].y;
      target[3*numTMDatoms+2] = atompos[i].z;
      aidmap[i] = numTMDatoms++;
      // add the atom to the list
      modifyRequestedAtoms().add(i);
    }
  }
  delete [] atompos;
  target_aid = new int[numTMDatoms];
  for (i=0; i<numTMDatoms; ++i) {
    target_aid[i] = modifyRequestedAtoms()[i];
  }
  DebugM(1,"done with parseAtoms\n");
}

GlobalMasterTMD::~GlobalMasterTMD() { 
  delete [] target;
  delete [] target_aid;
  delete [] aidmap;
}

void GlobalMasterTMD::calculate() {
  // have to reset my forces every time.  
  modifyAppliedForces().resize(0);
  modifyForcedAtoms().resize(0);

  // see if TMD should be active
  if (currentStep < firstStep || currentStep >= lastStep) {
    currentStep++;
    return;
  }

  // fetch the current coordinates
  AtomIDList::const_iterator a_i = getAtomIdBegin();
  AtomIDList::const_iterator a_e = getAtomIdEnd();
  PositionList::const_iterator p_i = getAtomPositionBegin();
  BigReal *curpos = new BigReal[3*numTMDatoms];
  for ( ; a_i != a_e; ++a_i, ++p_i ) {
    int ind = 3*aidmap[*a_i];
    curpos[ind  ] = p_i->x;
    curpos[ind+1] = p_i->y;
    curpos[ind+2] = p_i->z;
  }

  // align target with current coordinates.   Uses same weight for all
  // atoms.  Maybe instead use weight from occupancy?
  BigReal ttt[16], pre[3], post[3];
  BigReal curRMS = MatrixFitRMS(numTMDatoms, target, curpos, NULL, ttt);

  // Compute targetRMS.
  if (initialRMS < 0) {
    initialRMS = curRMS;
  }
  BigReal frac = (BigReal(currentStep-firstStep)) /
                            (lastStep-firstStep);
  BigReal targetRMS = initialRMS * (1-frac) + frac * finalRMS;

  BigReal maxforce2 = 0.;

  if ((finalRMS < initialRMS && targetRMS <= curRMS) ||
      (finalRMS >= initialRMS && targetRMS > curRMS)) {
    BigReal prefac = k * (targetRMS / curRMS - 1); 
    // compute transformation to align target structure with current structure
    // Might be more stable to instead align current positions with target,
    // although then we have to back-transform the forces.
    int j;
    for (j=0; j<3; j++) {
      post[j] = ttt[4*j+3];
      ttt[4*j+3] = 0;
      pre[j] = ttt[12+j];
      ttt[12+j] = 0;
    }
    Matrix4 result;
    result.translate(pre);
    result.multmatrix(Matrix4(ttt));
    result.translate(post);
  
    // compute forces on each atom
    BigReal myrms = 0;
    for (int i=0; i<numTMDatoms; i++) {
      result.multpoint(target+3*i);
      BigReal dx = curpos[3*i  ] - target[3*i  ];
      BigReal dy = curpos[3*i+1] - target[3*i+1];
      BigReal dz = curpos[3*i+2] - target[3*i+2];
  
      BigReal fvec[3] = { dx, dy, dz };
      Vector force(fvec[0]*prefac, fvec[1]*prefac, fvec[2]*prefac);
      modifyForcedAtoms().add(target_aid[i]);
      modifyAppliedForces().add(force);
      BigReal force2 = force.length2();
      if ( force2 > maxforce2 ) maxforce2 = force2;
    }
  }
  // write output if needed
  if (currentStep % outputFreq == 0) {
    iout << "TMD  " << currentStep << " " << targetRMS << ' ' << curRMS << '\n' << endi;
    // iout << "TMD  " << currentStep << " " << targetRMS << ' ' << curRMS << ' ' << sqrt(maxforce2) << '\n' << endi;
  }
  currentStep++;
  delete [] curpos;
}

