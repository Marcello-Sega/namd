/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeSMD.h"
#include "Node.h"
#include "SimParameters.h"
#include "ComputeMgr.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "PDB.h"
#include "PDBData.h"
#include "Molecule.h"
#include <iostream.h>
#if !defined(WIN32) || defined(__CYGWIN__)
#include <strstream.h>
#else
#include <strstrea.h>
#endif
#include <string.h>

ComputeSMD::ComputeSMD(ComputeGlobal *h)
: ComputeGlobalMaster(h) {
  SimParameters *simParams = Node::Object()->simParameters;

  moveVel = simParams->SMDVel;
  moveDir = simParams->SMDDir;
  outputFreq = simParams->SMDOutputFreq;
  k = simParams->SMDk;
  currentTime = simParams->firstTimestep;

  char *file = simParams->SMDFile;
  configMsg = new ComputeGlobalConfigMsg;
  parse_atoms(file);  
  iout << iINFO << "SMD ATOMS:  ";
  for (int i=0; i<configMsg->gdef.size()-1; i++)
    iout << configMsg->gdef[i]+1 << ' ';
  iout << "\n" << endi;
}

void ComputeSMD::parse_atoms(char *file) {
  PDB smdpdb(file);
  int numatoms = smdpdb.num_atoms();
  if (numatoms < 1) 
    NAMD_die("No atoms found in SMDFile\n");
  cm.x = cm.y = cm.z = 0;
  Molecule *mol = Node::Object()->molecule;
  BigReal imass = 0;
  for (int i=0; i<numatoms; i++) {
    PDBAtom *atom;
    if ((atom = smdpdb.atom(i))->occupancy()) {
      int ind = atom->serialnumber()-1;
      BigReal mass = mol->atommass(ind); 
      configMsg->gdef.add(atom->serialnumber()-1);
      cm.x += atom->xcoor()*mass;
      cm.y += atom->ycoor()*mass;
      cm.z += atom->zcoor()*mass;
      imass += mass; 
    }
  }
  configMsg->gdef.add(-1);
  if (imass == 0) 
    NAMD_die("SMDFile contained no SMD atoms - occupancy must be nonzero\n");
  cm /= imass;
}

ComputeSMD::~ComputeSMD() { }

void ComputeSMD::initialize() {
  // Send the configMsg we created in parse_atoms
  storedefs(configMsg->gdef);
  host->comm->sendComputeGlobalConfig(configMsg);
  configMsg = 0; 
}

void ComputeSMD::calculate() {

  ComputeGlobalResultsMsg *resultsMsg = new ComputeGlobalResultsMsg;
  resultsMsg->gforce.resize(gmass.size());
  resultsMsg->gforce.setall(Vector(0,0,0));

  Position curcm = gcom[0];
  BigReal diff = (curcm - cm)*moveDir;
  Force f = k*(moveVel*currentTime - diff)*moveDir;
  for (int i=0; i<gmass.size(); i++) 
    resultsMsg->gforce[i] = f;

  if (currentTime % outputFreq == 0)
    output(currentTime, curcm, f);

  // Send results to clients
  host->comm->sendComputeGlobalResults(resultsMsg);
  resultsMsg = 0;
  currentTime++;
}

void ComputeSMD::output(int t, Position p, Force f) {
  if (t % (100*outputFreq) == 0) 
    iout << "SMDTITLE: TS   CURRENT_POSITION         FORCE\n" << endi; 
  iout << "SMD  " << t << ' ' << p << ' ' << f*PNPERKCALMOL << '\n' << endi;
}

