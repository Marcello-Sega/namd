/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
#include "ComputeLCPO.h"
#include "ReductionMgr.h"
#include "PatchMap.h"
#include "ComputeMgr.h"
#include "Molecule.h"
#include "Node.h"
#include "SimParameters.h"
#include "Debug.h"
#include "WorkDistrib.decl.h"
#include "Node.h"
#include "ComputeLCPO.h"
#include "Priorities.h"
#include "PatchMap.inl"
#include "Patch.h"
#include "ComputeMap.h"
#include "LdbCoordinator.h"
#include "common.h"
#include "time.h"

//#define MIN_DEBUG_LEVEL 4
// #define DEBUGM

ComputeLCPO::ComputeLCPO(ComputeID c, PatchID p[], int t[], 
		ComputeNonbondedWorkArrays* _workArrays,
		int minPartition, int maxPartition, int numPartitions, int numPatches)
  : Compute(c), workArrays(_workArrays),
    minPart(minPartition), maxPart(maxPartition),
    strideIg(numPartitions), numParts(numPartitions),
    pairlistsMaxAge(10) {

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  setNumPatches(8);
  SimParameters *simParams = Node::Object()->simParameters;
  pairlistsMaxAge = (simParams->stepsPerCycle-1)/simParams->pairlistsPerCycle;
  surfTen = simParams->surface_tension;
  pairlistsAge = pairlistsMaxAge;

  //CkPrintf("ComputeLCPO[%d] uses patches ",cid);
  for (int i=0; i<getNumPatches(); i++) {
    patchID[i] = p[i];
    //CkPrintf("%02d, ",patchID[i]);
    trans[i] = t[i];
    patch[i] = NULL;
    positionBox[i] = NULL;
    forceBox[i] = NULL;
    lcpoTypeBox[i] = NULL;
  } // for all patches
  //CkPrintf("\n");
  
} // constructor

ComputeLCPO::~ComputeLCPO() {
  DebugM(4, "~ComputeLCPO("<<cid<<") numAtoms("<<patchID[0]<<") = " 
    << numAtoms[0] 
    << " numAtoms("<<patchID[1]<<") = " << numAtoms[1] << "\n" );
  DebugM(4, "~ComputeLCPO("<<cid<<") addr("<<patchID[0]<<") = " 
    << PatchMap::Object()->patch(patchID[0]) << " addr("<<patchID[1]<<") = "
    << PatchMap::Object()->patch(patchID[1]) << "\n");

  for (int i=0; i<getNumPatches(); i++) {
    if (positionBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterPositionPickup(cid,
	 &positionBox[i]);
    }
    if (forceBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterForceDeposit(cid,
		&forceBox[i]);
    }
    if (lcpoTypeBox[i] != NULL) {
      PatchMap::Object()->patch(patchID[i])->unregisterLcpoTypePickup(cid,
		&lcpoTypeBox[i]);
    }
  }
  delete reduction;
} // destructor

void ComputeLCPO::initialize() {
  Compute::initialize();
    // How can we tell if BoxOwner has packed up and left?  Need a mechanism
    // to handle this or do we assume the Boxes have been dumped?
    PatchMap *patchMap = PatchMap::Object();

    for (int i=0; i<8; i++) {
	    if (positionBox[i] == NULL) { // We have yet to get boxes
        patch[i] = PatchMap::Object()->patch(patchID[i]);
	      if (!(patch[i] = PatchMap::Object()->patch(patchID[i]))) {
	        DebugM(5,"invalid patch(" << patchID[i] 
		      << ")  pointer!\n");
	      }
	      positionBox[i] = patch[i]->registerPositionPickup(cid);
	      forceBox[i] = patch[i]->registerForceDeposit(cid);
	      lcpoTypeBox[i] = patch[i]->registerLcpoTypePickup(cid);
        // will need to open a box full of lcpo parameters
	    }
	    numAtoms[i] = patch[i]->getNumAtoms();
    } // for all patches

  DebugM(4, "initialize("<<cid<<") numAtoms("<<patchID[0]<<") = " 
    << numAtoms[0] 
    << " numAtoms(" <<patchID[1]<<") = " << numAtoms[1] << "\n" );

  // set priority

  basePriority = PATCH_PRIORITY(patchID[0]) + PROXY_RESULTS_PRIORITY;

/*
  //get offset between patches TODO
  for (int i=0; i<8; i++) {
    const Lattice &lattice = patch[i]->lattice;
    for (int j=0; j<8; j++) {
      offset[i][j] = lattice.offset(trans[i]) - lattice.offset(trans[j]);
    }
  }
*/

  //get bounds of inner rectangular prism in octet
  //bounds[3][2]; x, y, z; min, max
  bounds[0][0] = 0.5*(patchMap->min_a(patchID[0])+patchMap->max_a(patchID[0]));
  bounds[1][0] = 0.5*(patchMap->min_b(patchID[0])+patchMap->max_b(patchID[0]));
  bounds[2][0] = 0.5*(patchMap->min_c(patchID[0])+patchMap->max_c(patchID[0]));
  bounds[0][1] = 0.5*(patchMap->min_a(patchID[7])+patchMap->max_a(patchID[7]));
  bounds[1][1] = 0.5*(patchMap->min_b(patchID[7])+patchMap->max_b(patchID[7]));
  bounds[2][1] = 0.5*(patchMap->min_c(patchID[7])+patchMap->max_c(patchID[7]));

  oob[0] = bounds[0][0] > bounds[0][1] ? 1 : 0;
  oob[1] = bounds[1][0] > bounds[1][1] ? 1 : 0;
  oob[2] = bounds[2][0] > bounds[2][1] ? 1 : 0;

  pairlistsAge = pairlistsMaxAge;

} // initialize

void ComputeLCPO::atomUpdate() {
  for (int i=0; i<8; i++) {
	  numAtoms[i] = patch[i]->getNumAtoms();
  }
  pairlistsAge = pairlistsMaxAge;
}

//---------------------------------------------------------------------
// doWork
//---------------------------------------------------------------------
void ComputeLCPO::doWork() {
//  CkPrintf("PE%d CID%d doWork()\n", CkMyPe(), cid);
  LdbCoordinator::Object()->startWork(ldObjHandle);
  for (int i=0; i<8; i++) {
    pos[i] = positionBox[i]->open();
    force[i] = forceBox[i]->open();
    posExt[i] = patch[i]->getCompAtomExtInfo();
    lcpoType[i] = lcpoTypeBox[i]->open();
    //CkPrintf("%d\t\tPatch[%d,%d] lcpoTypeBox opened %d\n", sequence(), i, patchID[i], lcpoType[i]);
  }

  doForce();

 // Inform load balancer
  LdbCoordinator::Object()->endWork(ldObjHandle);

  // Close up boxes
  for (int i=0; i<getNumPatches(); i++) {
    positionBox[i]->close(&pos[i]);
    forceBox[i]->close(&force[i]);
    lcpoTypeBox[i]->close(&lcpoType[i]);
  }
} // doWork


int ComputeLCPO::noWork() {

  if ( patch[0]->flags.doFullElectrostatics ) {
    return 0;  // work to do, enqueue as usual
  } else {

    // skip all boxes
    for (int i=0; i<8; i++) {
      positionBox[i]->skip();
      forceBox[i]->skip();
      lcpoTypeBox[i]->skip();
    }

    reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
    reduction->submit();
    LdbCoordinator::Object()->skipWork(ldObjHandle);

    return 1;  // no work to do, do not enqueue
  }
  return 0;
} // noWork

// 1 - yes in bounds, 0 - not in bounds
// this does uniquely assign atoms to ComputeLCPO octets
int ComputeLCPO::isInBounds(Real x, Real y, Real z ) {

  //check x dimension
  if ( bounds[0][0] < bounds[0][1] ) { // internal
    if (x < bounds[0][0] || x >= bounds[0][1] )
      return 0;
  } else { // edge
    if (x < bounds[0][0] && x >= bounds[0][1] )
      return 0;
  }

  //check y dimension
  if ( bounds[1][0] < bounds[1][1] ) { // internal 
    if (y < bounds[1][0] || y >= bounds[1][1] )
      return 0;
  } else { // edge
    if (y < bounds[1][0] && y >= bounds[1][1] )
      return 0;
  }

  //check z dimension
  if ( bounds[2][0] < bounds[2][1] ) { // internal
    if (z < bounds[2][0] || z >= bounds[2][1] )
      return 0;
  } else { // edge
    if (z < bounds[2][0] && z >= bounds[2][1] )
      return 0;
  }

  return 1;
} // isInBounds

inline BigReal calcOverlap( Real r, Real ri, Real rj ) {
  return PI*ri*(2*ri-r-(ri*ri-rj*rj)/r);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//// doForce
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void ComputeLCPO::doForce() {
//CkPrintf("ComputeLCPO[%d]::doForce(%d) %d/%d\n",cid, sequence(),minPart,numParts);

  Real probeRadius = 1.4f;

  Position ngir, ngjr, ngkr;
  Real ri, rj, rk;
  BigReal dxij, dyij, dzij, r2ij, rij;
  BigReal dxik, dyik, dzik, r2ik, rik;
  BigReal dxjk, dyjk, dzjk, r2jk, rjk;

  int numPairs1 = 0;
  int numPairs2 = 0;

//////////////////////////////////////////////////
// Build Pairlists
//////////////////////////////////////////////////
if (pairlistsAge >= pairlistsMaxAge) {
  pairlistsAge = 0;
  //double t_start = 1.0*clock()/CLOCKS_PER_SEC;
  //CkPrintf("C%03d Rebuilding Pairlists %d/%d\n",cid,pairlistsAge,pairlistsMaxAge);
  pairlistsAge = 0;
  int pairsChecked = 0;

  inAtomsPl.reset();
  for (int pI = 0; pI < 8; pI++) {
    pairlists[pI].reset();
    //triplets[pI].reset();
  }
  //find in-bounds atoms in each patch
  for (int pI = 0; pI < 8; pI++) {
    //CkPrintf("Patch[%d] %d %d % 7.3f % 7.3f % 7.3f\n",patchID[pI],trans[pI],offset[pI].x,offset[pI].y,offset[pI].z);
    if (numAtoms[pI] == 0) continue;

    int minIg = 0;
    for (int s = 0; s < minPart; s++) {
      minIg += pos[pI][minIg].nonbondedGroupSize;
    }
    strideIg = numParts;//stride through partitions
    plint *inAtoms = inAtomsPl.newlist(numAtoms[pI]);
    int numAtomsInBounds = 0;

    //iterate over heavy atoms only
    for ( int ngi = minIg; ngi < numAtoms[pI]; /* ngi */) {
      ngir = pos[pI][ngi].position;
      if ( isInBounds(ngir.x, ngir.y, ngir.z) && 
          lcpoType[pI][ngi] > 0 ) {
        inAtoms[numAtomsInBounds++] = ngi;
        int idi = posExt[pI][ngi].id;
        //CkPrintf("C%03d P%03d ATOM %05d\n", cid, pI, idi);

        ri = probeRadius+lcpoParams[ lcpoType[pI][ngi] ][0];

        //find pairs in all 8 patches
        for (int pJ = 0; pJ < 8; pJ++) {
          if (numAtoms[pJ] == 0) continue;

          plint *neighbors = pairlists[pJ].newlist(numAtoms[pJ]); // pI
          int numNeighbors = 0;

          // j atom pairs
          for ( int ngj = 0; ngj < numAtoms[pJ]; /* ngj */) {
            pairsChecked++;
            rj = probeRadius+lcpoParams[ lcpoType[pJ][ngj] ][0];
            int idj = posExt[pJ][ngj].id;

            ngjr = pos[pJ][ngj].position;
            dxij = ngir.x - ngjr.x;
            dyij = ngir.y - ngjr.y;
            dzij = ngir.z - ngjr.z;
            r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
            rij = sqrt(r2ij);

            // radii + margin
            if (rij < (ri+rj+2) && rij > 0.01) {
              //found atom pair
              neighbors[numNeighbors] = ngj;
              //trips[jAtoms]   = ngj;
              numNeighbors++;
              //CkPrintf("C%03d P%03d ATOM %05d P%03d PAIR %05d\n",cid,pI,idi,pJ,idj);

              //CkPrintf("AtomPair[%05d,%05d] (% 7.3f,% 7.3f,% 7.3f) (% 7.3f,% 7.3f,% 7.3f) % 7.3f\n",
              //  posExt[pI][ngi].id,
              //  posExt[pJ][ngj].id,
              //  ngir.x, ngir.y, ngir.z,
              //  ngjr.x, ngjr.y, ngjr.z,
              //  sqrt(r2)
              //);
              numPairs1++;
            } // pair within cutoff

            //jump to next nonbonded group
            ngj += pos[pJ][ngj].nonbondedGroupSize;

          } // for pair atoms
          pairlists[pJ].newsize(numNeighbors); // pI
          //CkPrintf("PI%03d PJ%03d numNeighbors = %d\n", pI, pJ, numNeighbors);
          //CkPrintf("C%03d pairlists[%03d].newsize() %d\n", cid, pJ, numNeighbors);
        } // for patches J
      } // in bounds
      //jump to next nonbonded group for round-robin
      for (int s = 0; s < strideIg; s++) {
        ngi += pos[pI][ngi].nonbondedGroupSize;
      }
    } // for atoms
    inAtomsPl.newsize(numAtomsInBounds);
    //CkPrintf("C%03d inAtomsPl[%03d].newsize() %d\n", cid, pI, numAtomsInBounds);
  } // for patches I
  //double t_stop = 1.0*clock()/CLOCKS_PER_SEC;
  //CkPrintf("LCPO_TIME %7.3f pairs/us %d @ %f\n", pairsChecked*1e-6/(t_stop-t_start),pairsChecked,(t_stop-t_start));
} else {
  //CkPrintf("C%03d Reusing Pairlists %d/%d\n",cid,pairlistsAge,pairlistsMaxAge);
}

  inAtomsPl.reset();

  for (int pI = 0; pI < 8; pI++) {
    pairlists[pI].reset();
  }


  int numTrips = 0;
  BigReal totalSurfaceArea = 0;


  //double t_start = 1.0*clock()/CLOCKS_PER_SEC;
  int numTriplets = 0;

//////////////////////////////////////////////////
// Perform LCPO Calculation
//////////////////////////////////////////////////
  for (int pI = 0; pI < 8; pI++) {
    if (numAtoms[pI] == 0) continue;
    plint *inAtoms;
    int numInAtoms;
    inAtomsPl.nextlist( &inAtoms, &numInAtoms );
    //CkPrintf("C%03d S%d inAtomsPl[%03d].nextlist() %d\n", cid, sequence(), pI, numInAtoms);
    for (int i = 0; i < numInAtoms; i++) {
      int iPairs = 0;
      int iTrips = 0;
      int iIndex = inAtoms[i];
      int idi = posExt[pI][iIndex].id;
      Real xi = pos[pI][iIndex].position.x;
      Real yi = pos[pI][iIndex].position.y;
      Real zi = pos[pI][iIndex].position.z;
      const Real *lcpoParamI = lcpoParams[ lcpoType[pI][iIndex] ];
      ri = probeRadius+lcpoParamI[0];
      //CkPrintf("C%03d Atom[%07d](%2d) = % 7.3f, % 7.3f, % 7.3f\n",cid,idi,lcpoType[pI][iIndex],xi,yi,zi);

      Real P1 = lcpoParamI[1];
      Real P2 = lcpoParamI[2];
      Real P3 = lcpoParamI[3];
      Real P4 = lcpoParamI[4];

//////////////////////////////////////////////////
// S1
//////////////////////////////////////////////////
      BigReal S1 = 4.0*PI*ri*ri;

      //for surface area calculation
      BigReal AijSum = 0;
      BigReal AjkSum = 0;
      BigReal AjkjSum = 0;
      BigReal AijAjkSum = 0;

      //force
      BigReal dAijdrijdxiSum    = 0.0;
      BigReal dAijdrijdyiSum    = 0.0;
      BigReal dAijdrijdziSum    = 0.0;
      BigReal dAijdrijdxiAjkSum = 0.0;
      BigReal dAijdrijdyiAjkSum = 0.0;
      BigReal dAijdrijdziAjkSum = 0.0;


//////////////////////////////////////////////////
// Pairs
//////////////////////////////////////////////////
      plint *neighbors[8];
      int numNeighbors[8];
      for (int pJ = 0; pJ < 8; pJ++) {
        if (numAtoms[pJ] == 0) continue;
        pairlists[pJ].nextlist( &neighbors[pJ], &numNeighbors[pJ] );
        //CkPrintf("Pi%03d Pj%03d numNeighbors = %d\n", pI, pJ, numNeighbors[pJ]);
        //CkPrintf("C%03d pairlists[%03d].nextlist() %d\n", cid, pJ, numNeighbors[pJ]);
      }
      //  pairlists[pI].nextlist( &jAtoms, &numJAtoms );
      for (int pJ = 0; pJ < 8; pJ++) {
        if (numAtoms[pJ] == 0) continue;
        //CkPrintf("C%03d 2PATCH %03d:%03d %03d:%03d\n",cid,pI,numInAtoms,pJ,numNeighbors[pJ]);

        numPairs2 += numNeighbors[pJ];

        for (int j = 0; j < numNeighbors[pJ]; j++) {
          int jIndex = neighbors[pJ][j];
          int idj = posExt[pJ][jIndex].id;
          Real xj = pos[pJ][jIndex].position.x;
          Real yj = pos[pJ][jIndex].position.y;
          Real zj = pos[pJ][jIndex].position.z;
          const Real *lcpoParamJ = lcpoParams[ lcpoType[pJ][jIndex] ];
          rj = probeRadius+lcpoParamJ[0];

          dxij = xj-xi;
          dyij = yj-yi;
          dzij = zj-zi;
          r2ij = dxij*dxij + dyij*dyij + dzij*dzij;
          rij = sqrt(r2ij);
          BigReal rij_1 = 1.f / rij;

          if (rij >= (ri+rj) ||  rij < 0.01) {
            continue;
          }

          //CkPrintf("C%03d AtomPair[%07d,%07d] = % 7.3f\n",cid,idi,idj,rij);
          
//////////////////////////////////////////////////
// S2
//////////////////////////////////////////////////
          BigReal Aij = calcOverlap(rij, ri, rj);
          AijSum += Aij;
          iPairs++;

          //force
          BigReal dAijdrij = PI*ri*(rij_1*rij_1*(ri*ri-rj*rj)-1);
          BigReal dAijdrijdxj = dAijdrij*dxij*rij_1;
          BigReal dAijdrijdyj = dAijdrij*dyij*rij_1;
          BigReal dAijdrijdzj = dAijdrij*dzij*rij_1;


//////////////////////////////////////////////////
// Triplets
//////////////////////////////////////////////////
          BigReal AjkjSum = 0;
          BigReal AjkSum2 = 0;

          //force
          BigReal dAjkdrjkdxjSum = 0.0;
          BigReal dAjkdrjkdyjSum = 0.0;
          BigReal dAjkdrjkdzjSum = 0.0;

          for (int pK = 0; pK < 8; pK++) {
            if (numAtoms[pK] == 0) continue;
            //CkPrintf("C%03d 3PATCH %03d:%03d %03d:%03d %03d:%03d\n",cid,pI,numInAtoms,pJ,numNeighbors[pJ],pK,numNeighbors[pK]);

            for (int k = 0; k < numNeighbors[pK]; k++) {
              int kIndex = neighbors[pK][k];
              //CkPrintf("C%03d kIndex = %d (%d/%d)\n", cid, kIndex, k, numNeighbors[pK]);
              int idk = posExt[pK][kIndex].id;
              Real xk = pos[pK][kIndex].position.x;
              Real yk = pos[pK][kIndex].position.y;
              Real zk = pos[pK][kIndex].position.z;
              const Real *lcpoParamK = lcpoParams[ lcpoType[pK][kIndex] ];
              rk = probeRadius+lcpoParamK[0];

              //CkPrintf("C%03d P%03d ATOM%05d P%03d PAIR%05d P%03d TRIP%05d\n",cid, pI,idi,pJ,idj,pK,idk);

              dxik = xk-xi;
              dyik = yk-yi;
              dzik = zk-zi;
              r2ik = dxik*dxik + dyik*dyik + dzik*dzik;
              rik  = sqrt(r2ik);

              if (rik >= (ri+rk) ||  rik < 0.01) {
                continue;
              }

              dxjk = xk-xj;
              dyjk = yk-yj;
              dzjk = zk-zj;
              r2jk = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
              rjk  = sqrt(r2jk);

              if (rjk >= (rj+rk) ||  rjk < 0.01) {
                continue;
              }

              numTriplets++;

              BigReal rjk_1 = 1.0/rjk;
              BigReal Ajk = calcOverlap(rjk, rj, rk);
//////////////////////////////////////////////////
// S3
//////////////////////////////////////////////////
              AjkSum  += Ajk;
              AjkjSum += Ajk;
              AjkSum2 += Ajk;
              numTrips++;
              iTrips++;

              //force
              BigReal dAjkdrjk = PI*rj*rjk_1*(rjk_1*rjk_1*(rj*rj-ri*ri) - 1.f);

              BigReal dAjkdrjkdxj = dAjkdrjk*(xj-xk);
              BigReal dAjkdrjkdyj = dAjkdrjk*(yj-yk);
              BigReal dAjkdrjkdzj = dAjkdrjk*(zj-zk);

              //f(3*k-2) = f(3*k-2) - dajkddjkdxj*p3p4aij
              //f(3*k-1) = f(3*k-1) - dajkddjkdyj*p3p4aij
              //f(3*k  ) = f(3*k  ) - dajkddjkdzj*p3p4aij

              //p3p4aij = -surften*(p3(i) + p4(i)*aij)*frespa

              force[pK]->f[Results::nbond][kIndex].x -= -dAjkdrjkdxj*(P3+P4*Aij)*surfTen;
              force[pK]->f[Results::nbond][kIndex].y -= -dAjkdrjkdyj*(P3+P4*Aij)*surfTen;
              force[pK]->f[Results::nbond][kIndex].z -= -dAjkdrjkdzj*(P3+P4*Aij)*surfTen;

              /*CkPrintf("S%03d FK[%05d] = % 7.3f, % 7.3f, % 7.3f\n",
                sequence(),
                idk,
                dAjkdrjkdxj*(P3+P4*Aij)*surfTen,
                dAjkdrjkdyj*(P3+P4*Aij)*surfTen,
                dAjkdrjkdzj*(P3+P4*Aij)*surfTen
                ); */

              dAjkdrjkdxjSum += dAjkdrjkdxj;
              dAjkdrjkdyjSum += dAjkdrjkdyj;
              dAjkdrjkdzjSum += dAjkdrjkdzj;

            } // k atoms
//////////////////////////////////////////////////
// S4
//////////////////////////////////////////////////
          } // for patches K
          AijAjkSum += Aij*AjkjSum;

          //force
          dAijdrijdxiSum -= dAijdrijdxj;
          dAijdrijdyiSum -= dAijdrijdyj;
          dAijdrijdziSum -= dAijdrijdzj;

          dAijdrijdxiAjkSum -= dAijdrijdxj*AjkjSum;
          dAijdrijdyiAjkSum -= dAijdrijdyj*AjkjSum;
          dAijdrijdziAjkSum -= dAijdrijdzj*AjkjSum;

          BigReal lastxj = dAijdrijdxj*AjkjSum + Aij*dAjkdrjkdxjSum;
          BigReal lastyj = dAijdrijdyj*AjkjSum + Aij*dAjkdrjkdyjSum;
          BigReal lastzj = dAijdrijdzj*AjkjSum + Aij*dAjkdrjkdzjSum;

          BigReal dAidxj = (P2*dAijdrijdxj + P3*dAjkdrjkdxjSum + P4*lastxj);
          BigReal dAidyj = (P2*dAijdrijdyj + P3*dAjkdrjkdyjSum + P4*lastyj);
          BigReal dAidzj = (P2*dAijdrijdzj + P3*dAjkdrjkdzjSum + P4*lastzj);

          //f(3*j-2) = f(3*j-2) - daidxj
          //f(3*j-1) = f(3*j-1) - daidyj
          //f(3*j  ) = f(3*j  ) - daidzj
          force[pJ]->f[Results::nbond][jIndex].x -= dAidxj*surfTen;
          force[pJ]->f[Results::nbond][jIndex].y -= dAidyj*surfTen;
          force[pJ]->f[Results::nbond][jIndex].z -= dAidzj*surfTen;

          /*CkPrintf("S%03d FJ[%05d] = % 7.3f, % 7.3f, % 7.3f\n",
            sequence(),
            idj,
            -dAidxj*surfTen,
            -dAidyj*surfTen,
            -dAidzj*surfTen
            ); */

        } // for atom j
      } // for patches J

      // final atomic surface area calculation

      BigReal SAi = P1*S1 + P2*AijSum + P3*AjkSum + P4*AijAjkSum;
      SAi = (SAi > 0) ? SAi : 0.0;
      //CkPrintf("SurfaceArea[%05d] = % 7.3f\n", idi, SAi);

//CkPrintf("SA[%05i]=%7.3f (% 7.3e*% 7.3e + % 7.3e*% 7.3e + % 7.3e*% 7.3e + % 7.3e*% 7.3e) P%04i T%05i\n", idi, SAi, P1, S1, P2, AijSum, P3, AjkSum, P4, AijAjkSum, iPairs, iTrips);


      //force
      totalSurfaceArea += SAi;

      BigReal dAidxi = (P2*dAijdrijdxiSum + P4*dAijdrijdxiAjkSum);
      BigReal dAidyi = (P2*dAijdrijdyiSum + P4*dAijdrijdyiAjkSum);
      BigReal dAidzi = (P2*dAijdrijdziSum + P4*dAijdrijdziAjkSum);

      force[pI]->f[Results::nbond][iIndex].x -= dAidxi*surfTen;
      force[pI]->f[Results::nbond][iIndex].y -= dAidyi*surfTen;
      force[pI]->f[Results::nbond][iIndex].z -= dAidzi*surfTen;
      /*CkPrintf("S%03d FI[%05d] = % 7.3f, % 7.3f, % 7.3f\n",
        sequence(),
        idi,
        -dAidxi*surfTen,
        -dAidyi*surfTen,
        -dAidzi*surfTen
        ); */

      //params.ff[0] = r[a]->f[Results::nbond];
      //params.ff[1] = r[b]->f[Results::nbond];


    } // for inAtoms
  } // for patches I
  //double t_stop = 1.0*clock()/CLOCKS_PER_SEC;
  //CkPrintf("LCPO_TIME %7.3f trips/us %d @ %f\n", numTriplets*1e-6/(t_stop-t_start),numTriplets, (t_stop-t_start));

//  CkPrintf("C%03d NUMPAIRS1 %d\n", cid, numPairs1);
//  CkPrintf("C%03d NUMPAIRS2 %d\n", cid, numPairs2);
//  CkPrintf("C%03d NUMTRIPS %d\n", cid, numTrips);
//  CkPrintf("C%03d TOTSURF %f\n",  cid, totalSurfaceArea);



//////////////////////////////////////////////////
//  end calculation by submitting reduction
//////////////////////////////////////////////////
  for ( int i = 0; i < reductionDataSize; ++i )
    reductionData[i] = 0;
  reduction->item(REDUCTION_BC_ENERGY) += totalSurfaceArea * surfTen;
  submitReductionData(reductionData,reduction);
  reduction->submit();

  pairlistsAge++;

}//end do Force


// Lookup table for lcpo paramters
//indices 0 -> 22 are determined in Molecule.C
const Real ComputeLCPO::lcpoParams[23][5] = { //                      neigh
    { 0.00, 0.0000e+00,  0.0000e+00,  0.0000e+00, 0.0000e+00 }, //  0 H
    { 1.70, 7.7887e-01, -2.8063e-01, -1.2968e-03, 3.9328e-04 }, //  1 C sp3 1
    { 1.70, 5.6482e-01, -1.9608e-01, -1.0219e-03, 2.6580e-04 }, //  2 C sp3 2
    { 1.70, 2.3348e-01, -7.2627e-02, -2.0079e-04, 7.9670e-05 }, //  3 C sp3 3
    { 1.70, 0.0000e+00,  0.0000e+00,  0.0000e+00, 0.0000e+00 }, //  4 C sp3 4
    { 1.70, 5.1245e-01, -1.5966e-01, -1.9781e-04, 1.6392e-04 }, //  5 C sp2 2
    { 1.70, 7.0344e-02, -1.9015e-02, -2.2009e-05, 1.6875e-05 }, //  6 C sp2 3
    { 1.60, 7.7914e-01, -2.5262e-01, -1.6056e-03, 3.5071e-04 }, //  7 O sp3 1
    { 1.60, 4.9392e-01, -1.6038e-01, -1.5512e-04, 1.6453e-04 }, //  8 O sp3 2
    { 1.60, 6.8563e-01, -1.8680e-01, -1.3557e-03, 2.3743e-04 }, //  9 O sp2 1
    { 1.60, 8.8857e-01, -3.3421e-01, -1.8683e-03, 4.9372e-04 }, // 10 O=C-O
    { 1.65, 7.8602e-02, -2.9198e-01, -6.5370e-04, 3.6247e-04 }, // 11 N sp3 1
    { 1.65, 2.2599e-01, -3.6648e-02, -1.2297e-03, 8.0038e-05 }, // 12 N sp3 2
    { 1.65, 5.1481e-02, -1.2603e-02, -3.2006e-04, 2.4774e-05 }, // 13 N sp3 3
    { 1.65, 7.3511e-01, -2.2116e-01, -8.9148e-04, 2.5230e-04 }, // 14 N sp2 1
    { 1.65, 4.1102e-01, -1.2254e-01, -7.5448e-05, 1.1804e-04 }, // 15 N sp2 2
    { 1.65, 6.2577e-02, -1.7874e-02, -8.3120e-05, 1.9849e-05 }, // 16 N sp2 3
    { 1.90, 7.7220e-01, -2.6393e-01,  1.0629e-03, 2.1790e-04 }, // 17 S     1
    { 1.90, 5.4581e-01, -1.9477e-01, -1.2873e-03, 2.9247e-04 }, // 18 S     2
    { 1.90, 3.8650e-01, -1.8249e-01, -3.6598e-03, 4.2640e-04 }, // 19 P     3
    { 1.90, 3.8730e-02, -8.9339e-03,  8.3582e-06, 3.0381e-06 }, // 20 P     4
    { 1.80, 9.8318e-01, -4.0437e-01,  1.1249e-04, 4.9901e-04 }, // 21 Cl
    { 0.00, 0.0000e+00,  0.0000e+00,  0.0000e+00, 0.0000e+00 }  // 22 Mg
};

// SASA VMD LCPO gp 50831
// SASA VMD grid gp 41969
// SASA NAMDLCPO gp 48798 GREAT!
