/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifdef DPMTA
#include "Namd.h"
#include "Node.h"
#include "SimParameters.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMap.h"
#include "AtomMap.h"
#include "ComputeDPMTA.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "Communicate.h"
#include "pvmc.h"

#define MIN_DEBUG_LEVEL 1
// #define DEBUGM
#include "Debug.h"

extern Communicate *comm;

void ComputeDPMTA::get_FMA_cube(BigReal *boxsize, Vector *boxcenter)
{
  int max_dim;
  int dim_x,dim_y,dim_z;

  PatchMap *patchMap = PatchMap::Object();
  SimParameters *simParams = Node::Object()->simParameters;

  //  From these extremes, figure out how many patches we will
  //  have to have in each direction
  //  And add the skirt of empty patches by adding 2 patches
  //  in every direction
  dim_x = patchMap->xDimension() + 2;
  dim_y = patchMap->yDimension() + 2;
  dim_z = patchMap->zDimension() + 2;

  max_dim = dim_x;
  if (dim_y>max_dim)	max_dim = dim_y;
  if (dim_z>max_dim)	max_dim = dim_z;

  *boxsize = max_dim*simParams->patchDimension;
  *boxcenter = patchMap->Origin();
  boxcenter->x += *boxsize/2.0 - simParams->patchDimension;
  boxcenter->y += *boxsize/2.0 - simParams->patchDimension;
  boxcenter->z += *boxsize/2.0 - simParams->patchDimension;

  DebugM(2,"cube center: " << (*boxcenter) << " size=" << (*boxsize) << "\n");
}

ComputeDPMTA::ComputeDPMTA(ComputeID c) : ComputeHomePatches(c)
{
  // comm should always be initialized by this point...
  // In the (bug) case that it isn't, then initialize it.
  if (comm == NULL)
  {
    NAMD_die("Communication protocol (Converse, PVM, etc.) not initialized.");
  }

  //  Set everything to 0
  patchData = NULL;
  patchTail = NULL;
  numPatches = 0;
  numDistributed = 0;
  totalAtoms = 0;
  fmaResults = NULL;
  ljResults = NULL;
  local_timestep = 0;

  reduction->Register(REDUCTION_ELECT_ENERGY);

  //  NOTE that the theta value is hardwired to the value of 0.715
  //  as per the recommendation of the Duke developers

  //  NOTE 2: Theta is now an optional config parameter,
  //  but it defaults to 0.715

  int numProcs = CNumPes();
  PmtaInitData pmta_data;
  BigReal boxsize;	// Dimension of FMA cube
  Vector boxcenter;	// Center for FMA cube

  if (CMyPe() != 0)
  {
    slavetids=NULL;
    if (PMTAregister() < 0)
    {
	NAMD_die("PMTARegister failed!!");
    }
    DebugM(1,"DPMTA done PMTAinit.\n");
    return;
  }
  DebugM(1,"DPMTA configuring\n");

  // *****************************************
  // ONLY THE MASTER (NODE 0) NEEDS TO DO THIS:

  slavetids = new int[numProcs];
  if (slavetids == NULL)
  {
    NAMD_die("Memory allocation failed in FMAInterface::FMAInterface");
  }

  // pvm_spawn is a dummy function under Converse.  Just the array is required.
  pvm_spawn(NULL,NULL,0,NULL,numProcs,slavetids);
  DebugM(1,"DPMTA slavetids allocated\n");

  //  Get the size of the FMA cube
  DebugM(1,"DPMTA getting FMA cube\n");
  get_FMA_cube(&boxsize, &boxcenter);
  DebugM(1,"DPMTA got FMA cube\n");

  // reduce function calling time
  SimParameters *simParams = Node::Object()->simParameters;

  //  initialize DPMTA
  pmta_data.nprocs = numProcs;
  pmta_data.nlevels = simParams->FMALevels;
  pmta_data.mp = simParams->FMAMp;
  pmta_data.mp_lj = 4;
  pmta_data.fft = simParams->FMAFFTOn;
  pmta_data.fftblock = simParams->FMAFFTBlock;
  pmta_data.pbc = 0;
  pmta_data.kterm = 0;
  pmta_data.theta = simParams->fmaTheta;
  //  2.5 will allow non-cubical box
  pmta_data.cubelen.x = boxsize;
  pmta_data.cubelen.y = boxsize;
  pmta_data.cubelen.z = boxsize;
  pmta_data.cubectr.x = boxcenter.x;
  pmta_data.cubectr.y = boxcenter.y;
  pmta_data.cubectr.z = boxcenter.z;
  pmta_data.calling_num = pmta_data.nprocs;
  pmta_data.calling_tids = slavetids;

  DebugM(1,"DPMTA calling PMTAinit.\n");
  if (PMTAinit(&pmta_data,slavetids) >= 0)
  {
	iout << iINFO << "SUCCESSFULLY STARTED DPMTA\n" << endi;
  }
  else
  {
	iout << "Unable to start DPMTA!\n" << endi;
  }

  //  Register this master with the other DPMTA processes
  if (PMTAregister() < 0)
  {
	NAMD_die("PMTARegister failed!!");
  }
  DebugM(1,"DPMTA done PMTAinit.\n");

  DebugM(1,"DPMTA configured\n");
}

ComputeDPMTA::~ComputeDPMTA()
{
  DebugM(1,"DPMTA exiting\n");
  //  If this is the master node, then call PMTAexit()
  if (CMyPe() == 0)	PMTAexit();

  //  If there is a list of patch data hanging around, walk down the
  //  list and delete each node
  PatchInfo *ptr, *next;
  if (patchData != NULL)
  {
      ptr = patchData;
      while (ptr != NULL)
      {
         next = ptr->next;
         delete ptr;
         ptr = next;
      }
  }

  if (fmaResults)
	{
	free(fmaResults);
	fmaResults = NULL;
	}
  delete [] ljResults;
  delete [] slavetids;
  DebugM(1,"DPMTA exited\n");

  reduction->unRegister(REDUCTION_ELECT_ENERGY);
}


void ComputeDPMTA::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);
  PmtaParticle *particle_list = NULL;
  SimParameters *simParameters = Node::Object()->simParameters;

  // 0. only run when necessary
  // Skip computations if nothing to do.
  local_timestep++;
  DebugM(2,"fake_seq=" << fake_seq
	<< " timestep=" << local_timestep
	<< " fmaFrequency=" << simParameters->fmaFrequency
	<< "\n");
  if ((!patchList[0].p->flags.doFullElectrostatics)
	|| (local_timestep % simParameters->fmaFrequency != 1))
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      Position *x = (*ap).positionBox->open();
      AtomProperties *a = (*ap).atomBox->open();
      Force *f = (*ap).forceBox->open();
      reduction->submit(fake_seq, REDUCTION_ELECT_ENERGY, 0.0);
      ++fake_seq;
      (*ap).positionBox->close(&x);
      (*ap).atomBox->close(&a);
      (*ap).forceBox->close(&f);
    }
    return;
  }

  DebugM(1,"DPMTA doWork() started at timestep " << fake_seq << "\n");

  // setup
  // 1. get totalAtoms
  totalAtoms = Node::Object()->molecule->numAtoms;
  //  NOTE:  THIS IS LARGER THAN NEEDED, IN FUTURE COUNT ATOMS FIRST!  -JCP

  // 2. setup atom list
  int i,j;
  particle_list = (PmtaParticle *) calloc(totalAtoms, sizeof(PmtaParticle));
  fmaResults = (PmtaPartInfo *) calloc(totalAtoms, sizeof(PmtaPartInfo));
  if (!particle_list || !fmaResults)
	{
	NAMD_die("DPMTA Failed to allocate memory.");
	}

  BigReal unitFactor = sqrt(COLOUMB * ComputeNonbondedUtil::dielectric_1);
  for (i=0, ap = ap.begin(); ap != ap.end(); ap++)
  {
    (*ap).x = (*ap).positionBox->open();
    (*ap).a = (*ap).atomBox->open();

    // store each atom in the particle_list
     for(j=0; j<(*ap).p->getNumAtoms(); j++)
     {
      // explicitly copy -- two different data structures
      particle_list[i].p.x = (*ap).x[j].x;
      particle_list[i].p.y = (*ap).x[j].y;
      particle_list[i].p.z = (*ap).x[j].z;
      particle_list[i].q = (*ap).a[j].charge * unitFactor;
      i++;
      if (i > totalAtoms)
	{
	iout << iERRORF << iPE << " totalAtoms=" << totalAtoms
	     << " but " << i << " atoms are seen!\n" << endi;
	NAMD_die("FMA: atom counts unequal!");
	}
     }

    (*ap).positionBox->close(&(*ap).x);
    (*ap).atomBox->close(&(*ap).a);
  } 

  DebugM(1,"DPMTA doWork() there are " << i << " atoms in this node.\n");

  // 3. (run DPMTA) compute the forces
  if ( PMTAforce(i, particle_list, fmaResults, NULL) < 0 )
    {
      NAMD_die("PMTAforce failed!!");
    }

  // 4. deposit
  i=0;
  BigReal potential=0;
  for (ap = ap.begin(); ap != ap.end(); ap++)
  {
    (*ap).f = (*ap).forceBox->open();

    // deposit here
     for(j=0; j<(*ap).p->getNumAtoms(); j++)
     {
      (*ap).f[j].x += fmaResults[i].f.x;
      (*ap).f[j].y += fmaResults[i].f.y;
      (*ap).f[j].z += fmaResults[i].f.z;
      potential += fmaResults[i].v;
      i++;
     }

    (*ap).forceBox->close(&(*ap).f);
  }

  potential *= 0.5;
  DebugM(4,"Full-electrostatics energy: " << potential << "\n");
  reduction->submit(fake_seq, REDUCTION_ELECT_ENERGY, potential);
  ++fake_seq;

  // 5. clean-up
  if (totalAtoms > 0)
  {
    free(particle_list);
    free(fmaResults);
  }

  DebugM(1,"DPMTA doWork() done\n");
}

#endif

