/***************************************************************************/
/*         (C) Copyright 1996,1997 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#include "common.h"
#include "Namd.h"
#include "Node.h"
#include "SimParameters.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeDPMTA.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "Communicate.h"
#include "pvmc.h"
#include "InfoStream.h"
#include "ProcessorPrivate.h"

#define MIN_DEBUG_LEVEL 2
// #define DEBUGM
#include "Debug.h"

#ifdef DPMTA

// #define DUMP_DPMTA

void ComputeDPMTA::get_FMA_cube(int resize)
{
  Vector boxSize,boxCenter;	// used to see if things change
  PatchMap *patchMap = PatchMap::Object();

  if (usePBC == FALSE)
  {
    //  From these extremes, figure out how many patches we will
    //  have to have in each direction
    SimParameters *simParams = Node::Object()->simParameters;
    int dim_x = patchMap->xDimension();
    int dim_y = patchMap->yDimension();
    int dim_z = patchMap->zDimension();

    boxSize.x = dim_x*simParams->patchDimension;
    boxSize.y = dim_y*simParams->patchDimension;
    boxSize.z = dim_z*simParams->patchDimension;
    BigReal skirt = 2*simParams->patchDimension;

    boxCenter = patchMap->Origin();
    boxCenter.x += boxSize.x/2.0;
    boxCenter.y += boxSize.y/2.0;
    boxCenter.z += boxSize.z/2.0;

    //  add the skirt of empty patches by adding 2 patches in every direction
    boxSize.x += skirt;
    boxSize.y += skirt;
    boxSize.z += skirt;
  }
  else
  {
    DebugM(2,"getting patch info for FMA box\n");

    // determine boxSize from the PBC lattice
    // lattice is the same on all patches, so choose first patch
    ResizeArrayIter<PatchElem> ap(patchList);
    DebugM(2,"getting first patch info for FMA box\n");
    ap = ap.begin();
    DebugM(2,"getting lattice from patch for FMA box\n");
    Lattice lattice = (*ap).p->lattice;
    DebugM(2,"getting patch dimension for FMA box\n");
    boxSize = lattice.dimension();
    DebugM(2,"boxSize is " << boxSize << "\n");
    boxCenter = lattice.origin();
  }

  // don't bother checking if the center has moved since it depends on the size.
  if (boxsize != boxSize)
  {
	DebugM(2,"resetting FMA box\n");
	// reset the size and center
	boxsize = boxSize;
	boxcenter = boxCenter;

	// reset DPMTA (only reset it after it has been initialized!)
	if (resize && usePBC)
	{
	  PmtaVector center,v1,v2,v3;
	  center.x = boxcenter.x;
	  center.y = boxcenter.y;
	  center.z = boxcenter.z;
	  v1.x = boxsize.x;
	  v2.y = boxsize.y;
	  v3.z = boxsize.z;
	  iout << iINFO << "DPMTA box resized:\n";
	  iout << iINFO << "BOX DIMENSIONS = (" << v1.x << ","
		<< v2.y << "," << v3.z << ")\n";
	  iout << iINFO << "BOX CENTER = (" << center.x << ","
		<< center.y << "," << center.z << ")\n";
	  iout << endi;
	  DebugM(2,"calling PMTAresize()\n");
	  PMTAresize(&v1,&v2,&v3,&center);
	  DebugM(2,"called PMTAresize()\n");
	}
  }
  DebugM(2,"cube center: " << boxcenter << " size=" << boxsize << "\n");
}

ComputeDPMTA::ComputeDPMTA(ComputeID c) : ComputeHomePatches(c)
{
  ;
}

void ComputeDPMTA::initialize()
{
  ComputeHomePatches::initialize();

  DebugM(2,"ComputeDPMTA creating\n");
  // comm should always be initialized by this point...
  // In the (bug) case that it isn't, then initialize it.
  if (CpvAccess(comm) == NULL)
  {
    NAMD_die("Communication protocol (Converse, PVM, etc.) not initialized.");
  }

  // **** NOTE: node 0 must initialized before any other nodes register.

  //  Set everything to 0
  totalAtoms = 0;
  fmaResults = NULL;
  ljResults = NULL;
  boxcenter = 1;	// reset the array (no divide by zero)
  boxsize = 1;	// reset the array (no divide by zero)
  usePBC = FALSE;	// assume not...

  // all nodes should init
  reduction->Register(REDUCTION_ELECT_ENERGY);
  reduction->Register(REDUCTION_VIRIAL_SLOW_X);
  reduction->Register(REDUCTION_VIRIAL_SLOW_Y);
  reduction->Register(REDUCTION_VIRIAL_SLOW_Z);

  // Don't need any more initialization  -JCP
  ResizeArrayIter<PatchElem> ap(patchList);
  DebugM(2,"init() getting first patch info for FMA box\n");
  ap = ap.begin();
  DebugM(2,"init() getting lattice from patch for FMA box\n");
  initLattice = (*ap).p->lattice.dimension();
  DebugM(2,"init() initLattice is " << initLattice << "\n");

  //  NOTE that the theta value is hardwired to the value of 0.715
  //  as per the recommendation of the Duke developers

  //  NOTE 2: Theta is now an optional config parameter,
  //  but it defaults to 0.715

  // check for PBC
  usePBC = patchMap->xIsPeriodic()
	 + patchMap->yIsPeriodic()
	 + patchMap->zIsPeriodic();
  if ((usePBC != 0) && (usePBC != 3))
  {
    iout << iERROR << "DPMTA (FMA) does not support " << usePBC
	 << "-dimension PBC.\n" << endi;
  }
  DebugM(2,"Use PBC = " << usePBC << "\n");
  usePBC = (usePBC == 3);	// either PBC "3D" or no PBC
  if ( usePBC ) {
    iout << iWARN << "DPMTA (FMA) pressure tensor is incorrect.\n"
	<< iWARN << "Do not use DPMTA with anisotropic pressure control.\n"
	<< endi;
  }

  //  Get the size of the FMA cube
  DebugM(2,"DPMTA getting FMA cube\n");
  get_FMA_cube(FALSE);
  DebugM(2,"DPMTA got FMA cube\n");

  if (CkMyPe() != 0)
  {
    DebugM(2,"waiting for Init go-ahead\n");
    MIStream *msg2 = CpvAccess(comm)->newInputStream(ANY, DPMTATAG);
    int dummy;
    msg2->get(dummy);
    delete msg2;
    slavetids=NULL;
    if (PMTAregister() < 0)
    {
	NAMD_die("PMTARegister failed!!");
    }
    DebugM(2,"DPMTA done PMTAinit.\n");
    return;
  }
  DebugM(2,"DPMTA configuring\n");

  // *****************************************
  // ONLY THE MASTER (NODE 0) NEEDS TO DO THIS:

  int numProcs = CkNumPes();
  {
    int npatches=(PatchMap::Object())->numPatches();
    if ( numProcs > npatches ) numProcs = npatches;
  }

  slavetids = new int[numProcs];
  if (slavetids == NULL)
  {
    NAMD_die("Memory allocation failed in FMAInterface::FMAInterface");
  }

  // pvm_spawn is a dummy function under Converse.  Just the array is required.
  pvm_spawn(NULL,NULL,0,NULL,numProcs,slavetids);
  DebugM(2,"DPMTA slavetids allocated\n");

  // reduce function calling time
  SimParameters *simParams = Node::Object()->simParameters;

  //  initialize DPMTA
  PmtaInitData pmta_data;
  memset(&pmta_data,0,sizeof(pmta_data));
  pmta_data.nprocs = numProcs;
  pmta_data.nlevels = simParams->FMALevels;
  pmta_data.mp = simParams->FMAMp;
  pmta_data.mp_lj = 4;
  pmta_data.fft = simParams->FMAFFTOn;
  pmta_data.fftblock = simParams->FMAFFTBlock;
  pmta_data.pbc = usePBC;	// use Periodic boundary condition
  pmta_data.kterm = 0;
  pmta_data.theta = simParams->fmaTheta;
  pmta_data.v1.x = boxsize.x;
  pmta_data.v1.y = 0.;
  pmta_data.v1.z = 0.;
  pmta_data.v2.x = 0.;
  pmta_data.v2.y = boxsize.y;
  pmta_data.v2.z = 0.;
  pmta_data.v3.x = 0.;
  pmta_data.v3.y = 0.;
  pmta_data.v3.z = boxsize.z;
  pmta_data.cellctr.x = boxcenter.x;
  pmta_data.cellctr.y = boxcenter.y;
  pmta_data.cellctr.z = boxcenter.z;
  pmta_data.calling_num = pmta_data.nprocs;
  pmta_data.calling_tids = slavetids;

  iout << iINFO << "DPMTA parameters are:\n";
  iout << iINFO << "  LEVELS = " << pmta_data.nlevels << "\n";
  iout << iINFO << "  NUMBER OF MULTIPOLE TERMS = " << pmta_data.mp << "\n";
  iout << iINFO << "  FFT FLAG = " << pmta_data.fft << "\n";
  iout << iINFO << "  FFT BLOCKING FACTOR = " << pmta_data.fftblock << "\n";
  if ( usePBC ) iout << iINFO << "  SYSTEM IS PERIODIC\n" << endi;
  iout << iINFO << "  BOX DIMENSIONS = (" << pmta_data.v1.x << ","
	<< pmta_data.v2.y << "," << pmta_data.v3.z << ")\n";
  iout << iINFO << "  BOX CENTER = (" << pmta_data.cellctr.x << ","
	<< pmta_data.cellctr.y << "," << pmta_data.cellctr.z << ")\n";
  iout << endi;

  if ( usePBC )
  {
    pmta_data.cellctr.x = 0.;
    pmta_data.cellctr.y = 0.;
    pmta_data.cellctr.z = 0.;
  }

#ifdef DUMP_DPMTA
  FILE *fp;
  fp = fopen("DUMP_DPMTA.init","w");
  fwrite(&pmta_data,sizeof(PmtaInitData),1,fp);
  fclose(fp);
#endif

  DebugM(2,"DPMTA calling PMTAinit.\n");
  if (PMTAinit(&pmta_data,slavetids) >= 0)
  {
	iout << iINFO << "SUCCESSFULLY STARTED DPMTA\n" << endi;
  }
  else
  {
	NAMD_die("Unable to start DPMTA!");
  }

  // tell all nodes that it is OK to register
  MOStream *msg = CpvAccess(comm)->newOutputStream(ALL, DPMTATAG, BUFSIZE);
  // don't actually put in data...  Nodes just need it as a flag.
  msg->put(TRUE);
  msg->end();
  delete msg;
  DebugM(2,"Init go-ahead\n");
  MIStream *msg1 = CpvAccess(comm)->newInputStream(ANY, DPMTATAG);
  int dummy1;
  msg1->get(dummy1);
  delete msg1;
  DebugM(2,"got Init go-ahead\n");

  //  Register this master with the other DPMTA processes
  if (PMTAregister() < 0)
  {
	NAMD_die("PMTARegister failed!!");
  }
  DebugM(2,"DPMTA done PMTAinit.\n");
  DebugM(2,"DPMTA configured\n");
}

ComputeDPMTA::~ComputeDPMTA()
{
  DebugM(2,"DPMTA exiting\n");
  //  If this is the master node, then call PMTAexit()
  if (CkMyPe() == 0)	PMTAexit();

  if (fmaResults)
	{
	free(fmaResults);
	fmaResults = NULL;
	}
  delete [] ljResults;
  delete [] slavetids;
  DebugM(2,"DPMTA exited\n");

  reduction->unRegister(REDUCTION_ELECT_ENERGY);
  reduction->unRegister(REDUCTION_VIRIAL_SLOW_X);
  reduction->unRegister(REDUCTION_VIRIAL_SLOW_Y);
  reduction->unRegister(REDUCTION_VIRIAL_SLOW_Z);
}


void ComputeDPMTA::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);
  PmtaParticle *particle_list = NULL;

  // 0. only run when necessary
  // Skip computations if nothing to do.
  if (!patchList[0].p->flags.doFullElectrostatics)
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      Position *x = (*ap).positionBox->open();
      AtomProperties *a = (*ap).atomBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).forceBox->close(&r);
      (*ap).atomBox->close(&a);
      (*ap).positionBox->close(&x);
    }
    reduction->submit(patchList[0].p->flags.seq, REDUCTION_ELECT_ENERGY, 0.0);
    reduction->submit(patchList[0].p->flags.seq, REDUCTION_VIRIAL_SLOW_X, 0.0);
    reduction->submit(patchList[0].p->flags.seq, REDUCTION_VIRIAL_SLOW_Y, 0.0);
    reduction->submit(patchList[0].p->flags.seq, REDUCTION_VIRIAL_SLOW_Z, 0.0);
    return;
  }

  DebugM(2,"DPMTA doWork() started at timestep " << patchList[0].p->flags.seq << "\n");

  // setup
  // 1. get totalAtoms
  for (totalAtoms=0, ap = ap.begin(); ap != ap.end(); ap++)
     totalAtoms += (*ap).p->getNumAtoms();

  Vector newLattice;
  Vector rescaleFactor;
  if (usePBC)
    {
    ap = ap.begin();
    Lattice lattice = (*ap).p->lattice;
    newLattice = lattice.dimension();
    rescaleFactor.x = initLattice.x / newLattice.x;
    rescaleFactor.y = initLattice.y / newLattice.y;
    rescaleFactor.z = initLattice.z / newLattice.z;
    DebugM(2,"Rescale factor = " << initLattice << "/" << newLattice
		<< " = " << rescaleFactor << "\n");
    DebugM(2,"boxcenter = " << boxcenter << "\n");
    }
  else
    {
    rescaleFactor.x = 1;
    rescaleFactor.y = 1;
    rescaleFactor.z = 1;
    }

  // 2. setup atom list
  int i,j;
  particle_list = (PmtaParticle *)calloc(totalAtoms,sizeof(PmtaParticle));
  fmaResults =    (PmtaPartInfo *)calloc(totalAtoms,sizeof(PmtaPartInfo));
  if (!particle_list || !fmaResults)
	{
	NAMD_die("DPMTA Failed to allocate memory.");
	}

  BigReal unitFactor = sqrt(COLOUMB * ComputeNonbondedUtil::dielectric_1);
  DebugM(2,"Charge unit factor = " << unitFactor << "\n");
  for (i=0, ap = ap.begin(); ap != ap.end(); ap++)
  {
    Vector *x = (*ap).positionBox->open();
    AtomProperties *a = (*ap).atomBox->open();

    // store each atom in the particle_list
    Vector pos;
    for(j=0; j<(*ap).p->getNumAtoms(); j++)
    {
      // explicitly copy -- two different data structures
      if (usePBC)
	{
	particle_list[i].p.x = rescaleFactor.x * (x[j].x-boxcenter.x);
	particle_list[i].p.y = rescaleFactor.y * (x[j].y-boxcenter.y);
	particle_list[i].p.z = rescaleFactor.z * (x[j].z-boxcenter.z);
	}
      else
	{
	particle_list[i].p.x = x[j].x;
	particle_list[i].p.y = x[j].y;
	particle_list[i].p.z = x[j].z;
	}
      particle_list[i].q = a[j].charge * unitFactor;
      DebugM(1,"atom[" << i << "]=" << x[j] << " "
	      << a[j].charge*unitFactor << "\n");
      i++;
      if (i > totalAtoms)
	{
	iout << iERRORF << iPE << " totalAtoms=" << totalAtoms
	     << " but " << i << " atoms are seen!\n" << endi;
	NAMD_die("FMA: atom counts unequal!");
	}
    }

    (*ap).atomBox->close(&a);
    (*ap).positionBox->close(&x);
  } 

  if (i != totalAtoms)
  {
    iout << iERRORF << iPE << " totalAtoms=" << totalAtoms
         << " but " << i << " atoms are seen!\n" << endi;
    NAMD_die("FMA: atom counts unequal!");
  }

  DebugM(2,"DPMTA doWork() there are " << totalAtoms << " atoms in this node.\n");

#ifdef DUMP_DPMTA
  FILE *fp;
  char dump_file[32];
  sprintf(dump_file,"DUMP_DPMTA.%d",(int)CkMyPe());
  fp = fopen(dump_file,"w");
  int32 n32 = i;
  fwrite(&n32,sizeof(int32),1,fp);
  fwrite(particle_list,sizeof(PmtaParticle),i,fp);
  fclose(fp);
#endif

  // 3. (run DPMTA) compute the forces
  if ( PMTAforce(i, particle_list, fmaResults, NULL) < 0 )
    {
      NAMD_die("PMTAforce failed!!");
    }
  DebugM(2,"DPMTA forces done.  Now depositing.\n");

  // 4. deposit
  BigReal potential=0;
  for (i=0, ap = ap.begin(); ap != ap.end(); ap++)
  {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];

    // deposit here
    for(j=0; j<(*ap).p->getNumAtoms(); j++)
    {
      if (usePBC)
	{
	f[j].x += fmaResults[i].f.x * rescaleFactor.x * rescaleFactor.x;
	f[j].y += fmaResults[i].f.y * rescaleFactor.y * rescaleFactor.y;
	f[j].z += fmaResults[i].f.z * rescaleFactor.z * rescaleFactor.z;
	potential += fmaResults[i].v * rescaleFactor.x;
	}
      else
	{
	f[j].x += fmaResults[i].f.x;
	f[j].y += fmaResults[i].f.y;
	f[j].z += fmaResults[i].f.z;
	potential += fmaResults[i].v;
	}
      i++;
    }

    (*ap).forceBox->close(&r);
  }

  potential *= 0.5;
  DebugM(4,"Full-electrostatics energy: " << potential << "\n");
  reduction->submit(patchList[0].p->flags.seq, REDUCTION_ELECT_ENERGY, potential);
  // DPMTA won't work correctly if scaled anisotropically anyway.  -JCP
  reduction->submit(patchList[0].p->flags.seq, REDUCTION_VIRIAL_SLOW_X, potential / 3.);
  reduction->submit(patchList[0].p->flags.seq, REDUCTION_VIRIAL_SLOW_Y, potential / 3.);
  reduction->submit(patchList[0].p->flags.seq, REDUCTION_VIRIAL_SLOW_Z, potential / 3.);

  // 5. clean-up
  if (totalAtoms > 0)
  {
    free(particle_list);
    free(fmaResults);
  }

  DebugM(2,"DPMTA doWork() done\n");
}

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1055 $	$Date: 1999/05/11 23:56:19 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDPMTA.C,v $
 * Revision 1.1055  1999/05/11 23:56:19  brunner
 * Changes for new charm version
 *
 * Revision 1.1054  1999/04/23 20:13:52  jim
 * Dies on startup failure.
 *
 * Revision 1.1053  1999/02/17 04:09:55  jim
 * Fixes to make optional force modules work with more nodes than patches.
 *
 * Revision 1.1052  1999/01/06 00:56:20  jim
 * All compute objects except DPMTA now return diagonal of virial tensor.
 *
 * Revision 1.1051  1998/10/24 19:57:22  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1050  1998/06/18 14:48:00  jim
 * Split virial into NORMAL, NBOND, and SLOW parts to match force classes.
 *
 * Revision 1.1049  1998/04/15 22:13:50  jim
 * Make depends returns same results regardless of DPME, DPMTA, TCL or MDCOMM.
 *
 * Revision 1.1048  1997/11/10 16:45:54  milind
 * Made comm a Cpv Variable.
 *
 * Revision 1.1047  1997/10/06 00:12:28  jim
 * Added PatchMap.inl, sped up cycle-boundary tuple code.
 *
 * Revision 1.1046  1997/10/01 16:46:47  milind
 * Removed old NAMD1 messaging and replaced it with new Message Streams library.
 *
 * Revision 1.1045  1997/09/18 22:49:41  jim
 * Added code to dump DPMTA input to files.  Define DUMP_DPMTA to enable.
 *
 * Revision 1.1044  1997/09/12 23:05:51  jim
 * Changes to work with DPMTA-2.6
 *
 * Revision 1.1043  1997/04/04 23:34:15  milind
 * Got NAMD2 to run on Origin2000.
 * Included definitions of class static variables in C files.
 * Fixed alignment bugs by using memcpy instead of assignment in
 * pack and unpack.
 *
 * Revision 1.1042  1997/03/27 17:08:30  nealk
 * Added hydrogen groupings.  Now configuration parameter "splitPatch" determines
 * atom-into-patch distribution.
 *
 * Revision 1.1041  1997/03/27 16:04:49  nealk
 * Removed init() -- no longer necessary.  Thanks Jim!
 * Turned off debugging.
 *
 * Revision 1.1040  1997/03/27 08:04:14  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1039  1997/03/27 03:16:50  jim
 * Added code to check virial calculation, fixed problems with DPMTA and PBC's.
 *
 * Revision 1.1038  1997/03/25 16:57:47  nealk
 * Added PBC scaling to DPMTA.
 * Turned off debugging code in Controller.C.
 *
 * Revision 1.1037  1997/03/20 23:53:33  ari
 * Some changes for comments. Copyright date additions.
 * Hooks for base level update of Compute objects from ComputeMap
 * by ComputeMgr.  Useful for new compute migration functionality.
 *
 * Revision 1.1036  1997/03/20 19:45:08  nealk
 * Moved syncronize point in doWork to fix timing problem.  (Very weird.)
 *
 * Revision 1.1035  1997/03/19 18:47:29  nealk
 * Added log info to Hydrogen.h
 * Fixed ComputeDPMTA.C so node 0 initializes before any other nodes register
 * with the DPMTA library.
 *
 * Revision 1.1034  1997/03/19 18:10:10  nealk
 * Added sorted hydrogen group list to molecule.
 *
 * Revision 1.1033  1997/03/19 11:54:05  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
