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
#include "GlobalMasterSymmetry.h"
#include "PDB.h"
#include "PDBData.h"
#include "fitrms.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"
//using namespace __gnu_cxx;
using namespace std;
//Read in and parse matrix file
void GlobalMasterSymmetry::parseMatrix(char fileName []){
  
  string line;
  ifstream matrixFile (fileName);
  if (matrixFile.is_open()){
    while (!matrixFile.eof()){
      vector <string> tmp;
      for (int i = 0; i < 4; i++){
      getline(matrixFile, line);
      istringstream iss(line);
      copy(istream_iterator<string>(iss),
            istream_iterator<string>(),
            back_inserter<vector <string>  >(tmp));
      }
      getline(matrixFile, line);
      Matrix4Symmetry tmpmatrix;
      for(int j = 0; j < 16; j++){
        tmpmatrix.mat[j] = atof(tmp[j].c_str());
      }
      tmpmatrix.transpose();
      matrices.add(tmpmatrix);
    }
    matrixFile.close();
  } 
}
GlobalMasterSymmetry::GlobalMasterSymmetry() {
  DebugM(3,"initialize called\n");
  SimParameters *params = Node::Object()->simParameters;
  currentStep = params->firstTimestep;
  firstStep = params->symmetryFirstStep;
  lastStep = params->symmetryLastStep;
  firstFullStep = params->symmetryFirstFullStep;
  lastFullStep = params->symmetryLastFullStep;
  K = params->symmetryk;
  scaleForces = params->symmetryScaleForces;
  if (scaleForces && lastStep == -1){
    NAMD_die("symmetryLastStep must be specified if symmetryScaleForces is enabled!");
  }
  parseAtoms(params->symmetryFile, Node::Object()->molecule->numAtoms);
  if (params->symmetryMatrixFile[0] != '\0'){
    parseMatrix(params->symmetryMatrixFile);
  }
  else{
    initialTransform();
  }
  DebugM(1,"done with initialize\n");
}

//Aligns monomers based on transformation matrices
//found in matrix file
void GlobalMasterSymmetry::alignMonomers(){
  //this is assuming the matrices are written
  //in order of monomer id designation (0, 1, 2,..etc)
  map<int, vector<int> >::iterator mit = dmap.begin();
  for (; mit!=dmap.end(); ++mit){
    for (int i = 0; i < mit->second.size(); i++){
      map<int, BigReal *>::iterator it = posmap.find(mit->second[i]);
      matrices[(mit->first)-1].multpoint(it->second);
    }
  }
}

bool GlobalMasterSymmetry::gluInvertMatrix(const BigReal m[16], BigReal invOut[16])
{
  BigReal inv[16], det;
  int i;

  inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
  + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
  inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
  - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
  inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
  + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
  inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
  - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
  inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
  - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
  inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
  + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
  inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
  - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
  inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
  + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
  inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
  + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
  inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
  - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
  inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
  + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
  inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
  - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
  inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
  - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
  inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
  + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
  inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
  - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
  inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
  + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

  det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
  if (det == 0)
            return false;

            det = 1.0 / det;

            for (i = 0; i < 16; i++)
                      invOut[i] = inv[i] * det;

                      return true;
}

void GlobalMasterSymmetry::initialTransform(){
  map<int, vector<int> >::iterator fit = dmap.find(1);
  BigReal * first = new BigReal [3*fit->second.size()];
  for(int i = 0; i < fit->second.size(); i++){
    map<int, BigReal *>::iterator sit = startmap.find(fit->second[i]);
    first[3*i] = sit->second[0];
    first[3*i+1] = sit->second[1];
    first[3*i+2] = sit->second[2];  
  }

  map<int, vector<int> >::iterator it = dmap.begin();
  for(; it != dmap.end(); ++it){
    BigReal * arr = new BigReal [3*it->second.size()];
    for (int i = 0; i<it->second.size(); i++){
      map<int, BigReal *>::iterator sit = startmap.find(it->second[i]);
      arr[3*i] = sit->second[0];
      arr[3*i+1] = sit->second[1];
      arr[3*i+2] = sit->second[2];  
    }  
      BigReal ttt[16], pre[3], post[3];
      BigReal curRMS = MatrixFitRMS(it->second.size(), arr, first, NULL, ttt);
      int j;
      for (j=0; j<3; j++) {
        post[j] = ttt[4*j+3];
        ttt[4*j+3] = 0;
        pre[j] = ttt[12+j];
        ttt[12+j] = 0;
      }
      Matrix4Symmetry result;
      result.translate(pre);
      result.multmatrix(Matrix4Symmetry(ttt));
      result.translate(post);
      matrices.add(result);

      delete [] arr;
  }
  delete [] first;
}
void GlobalMasterSymmetry::parseAtoms(const char *file, int numTotalAtoms) {
  DebugM(3,"parseAtoms called\n");
  PDB tmdpdb(file);
  int numatoms = tmdpdb.num_atoms();
  if (numatoms < 1) 
    NAMD_die("No atoms found in symmetryFile\n");
  if (numatoms != numTotalAtoms)
    NAMD_die("The number of atoms in symmetryFile must be equal to the total number of atoms in the structure!");
  if ( modifyRequestedAtoms().size() )
    NAMD_bug("GlobalMasterSymmetry::parseAtoms() modifyRequestedAtoms() not empty");

  Vector * atompos = new Vector[numatoms];  
  tmdpdb.get_all_positions(atompos);
  int i;
  for (i=0; i<numatoms; i++) {
#ifdef MEM_OPT_VERSION
    PDBCoreData *atom = tmdpdb.atom(i);
#else
    PDBAtom *atom = tmdpdb.atom(i); // get an atom from the file
#endif
    if (atom->occupancy() && atom->temperaturefactor()) { // if occupancy and beta are not 0, then add it!
      // add the atom to the list
      modifyRequestedAtoms().add(i);
      if(!K){ kmap[i] = atom->occupancy();}

        //check to see if monomer id is already in the map
      map <int, vector<int> >::iterator it = dmap.find((int)atom->temperaturefactor());
      if (it != dmap.end()){
        it->second.push_back(i); //add atomid to vector in proper existing monomer id
      }
      else{
         dmap[(int)atom->temperaturefactor()] = vector <int> (1,i); //create new monomer id with atomid
      }  
    }
  }
  map <int, vector<int> >::iterator it = dmap.begin();
  int atoms_in_monomer = it->second.size();
  for (; it!=dmap.end(); ++it){
    if (it->second.size() != atoms_in_monomer){
      NAMD_die("Every monomer must contain the same number of atoms!");
    }
  }
  delete [] atompos;
}

void GlobalMasterSymmetry::determineAverage() {
   for(int i=0; i<averagePos.size(); i++){delete [] averagePos[i];}
   averagePos.erase(averagePos.begin(), averagePos.end());

   map <int, vector<int> >::iterator it = dmap.begin();
   map <int, BigReal *>::iterator posit;
   int numatoms = it->second.size();
   for (int i = 0; i < numatoms; i++){
     BigReal *arr = new BigReal [3];
     arr[0] = 0;
     arr[1] = 0;
     arr[2] = 0;
     for (; it!=dmap.end(); ++it){
       posit = posmap.find(it->second[i]);
       arr[0] += posit->second[0];
       arr[1] += posit->second[1];
       arr[2] += posit->second[2];
     }
     it = dmap.begin();
     BigReal *avg = new BigReal[3];
     avg[0] = arr[0]/(dmap.size());
     avg[1] = arr[1]/(dmap.size());
     avg[2] = arr[2]/(dmap.size());
     averagePos.push_back(avg);
     delete [] arr;
   }
}

void GlobalMasterSymmetry::backTransform(){
  map <int, vector<int> >::iterator it = dmap.begin();
  int numatoms = it->second.size();
  map <int, BigReal *>::iterator bit = backavg.begin();
  for (; bit != backavg.end(); ++bit){delete [] bit->second;}

  backavg.clear();
  for (; it!=dmap.end(); ++it){
    BigReal *avg = new BigReal [3*averagePos.size()];
    for (int i = 0; i < numatoms; i++){
      avg[3*i] = averagePos[i][0];
      avg[3*i+1] = averagePos[i][1];
      avg[3*i+2] = averagePos[i][2];
    }
    BigReal inverse[16];
    matrices[it->first-1].transpose();
    gluInvertMatrix(matrices[it->first-1].mat, inverse);
    Matrix4Symmetry inv(inverse);
    inv.transpose();
    matrices[it->first-1].transpose();
 
    for (int k = 0; k < numatoms; k++){
      inv.multpoint(avg+3*k);
    }
    backavg[it->first] = avg;
   }
}
GlobalMasterSymmetry::~GlobalMasterSymmetry() { 
  map <int, BigReal *>::iterator pit = posmap.begin();
  for (; pit != posmap.end(); ++pit){delete pit->second;}
  for(int i=0; i<averagePos.size(); i++){delete averagePos[i];}
}
void GlobalMasterSymmetry::calculate() {
  // have to reset my forces every time.  
  modifyAppliedForces().resize(0);
  modifyForcedAtoms().resize(0);

  // see if symmetry restraints should be active
  if (currentStep < firstStep || (currentStep >= lastStep && lastStep != -1)) {
    currentStep++;
    return;
  }

  map<int, Position> positions;
  AtomIDList::const_iterator a_i = getAtomIdBegin();
  AtomIDList::const_iterator a_e = getAtomIdEnd();
  PositionList::const_iterator p_i = getAtomPositionBegin();

  //create mapping of positions now to avoid
  //going through these iterators for each domain
  for ( ; a_i != a_e; ++a_i, ++p_i ){
    positions[*a_i] = *p_i;
  }

  map <int, vector<int> >::iterator it;
  map <int, BigReal *>::iterator pit = posmap.begin();
  for (; pit != posmap.end(); ++pit){delete [] pit->second;}

  posmap.clear();
  for (it = dmap.begin(); it != dmap.end(); ++it){
    // fetch the current coordinates
    for (int i = 0; i < it->second.size(); i++){

      BigReal* arr = new BigReal[3];
      arr[0] = positions[it->second[i]].x;
      arr[1] = positions[it->second[i]].y;
      arr[2] = positions[it->second[i]].z;

      posmap[it->second[i]] = arr;
    } 
}

  alignMonomers();
  determineAverage();
  backTransform();

  //iterate through all domains
  for (it = dmap.begin(); it != dmap.end(); ++it){

  // fetch the current coordinates
  BigReal *curpos = new BigReal[3*(it->second.size())];
  for (int i = 0; i < it->second.size(); i++){
    curpos[3*i  ] = positions[it->second[i]].x;
    curpos[3*i+1] = positions[it->second[i]].y;
    curpos[3*i+2] = positions[it->second[i]].z;  
  }


  BigReal *tmpavg = backavg[it->first];

  for (int i=0; i<it->second.size(); i++) {
    BigReal k;  
    if(!K){
     k = kmap[it->second[i]];  
    }
    else{
     k = K/it->second.size();  
    }
    BigReal maxforce2 = 0.;
    BigReal frac = -k;
    if (scaleForces){

    if (currentStep < firstFullStep){
      BigReal linear_evolve = (BigReal(currentStep-firstStep)) /
                            (firstFullStep-firstStep);
      frac = frac*(linear_evolve);
    }
    if (currentStep > lastFullStep)
    {
      BigReal linear_evolve = (BigReal(currentStep-lastFullStep)) /
                            (lastStep-lastFullStep);
      frac = frac*(1-linear_evolve);
    }

    }
    BigReal dx = curpos[3*i  ] - tmpavg[3*i];
    BigReal dy = curpos[3*i+1] - tmpavg[3*i+1];
    BigReal dz = curpos[3*i+2] - tmpavg[3*i+2];
    BigReal fvec[3] = { dx, dy, dz };
    Vector force(fvec[0]*frac, fvec[1]*frac, fvec[2]*frac);
    modifyForcedAtoms().add(it->second[i]);
    modifyAppliedForces().add(force);
    BigReal force2 = force.length2();
    if ( force2 > maxforce2 ) maxforce2 = force2;
  } 
  delete [] curpos;
 }
  currentStep++;
}
