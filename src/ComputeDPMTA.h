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

#ifndef COMPUTEDPMTA_H
#define COMPUTEDPMTA_H

#ifdef DPMTA

extern "C"
  {
  #include "dpmta.h"
  }
#include "ComputeHomePatches.h"

class ComputeDPMTA : public ComputeHomePatches {
private:
  int *slavetids;	//  PID for slave processes
  int totalAtoms;	//  Total number of atoms being dealt with
  PmtaPartInfo *fmaResults;	//  Results from the PMTA code
  PmtaPartInfo *ljResults;	//  Results from the PMTA code
  Vector boxsize;	// FMA box size, set by get_FMA_cube()
  Vector boxcenter;	// FMA box center, set by get_FMA_cube()
  int usePBC;		// flag for PBC

  void get_FMA_cube(int resize);

public:
  ComputeDPMTA(ComputeID c);
  virtual ~ComputeDPMTA();
  void doWork();
};

#endif
#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeDPMTA.h,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1013 $	$Date: 1997/03/11 16:37:04 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDPMTA.h,v $
 * Revision 1.1013  1997/03/11 16:37:04  nealk
 * Difference in ElectForce with FMA and non-FMA not due to FMA.
 * (I'm finally convinced.)
 *
 * Revision 1.1012  1997/03/10 17:40:02  ari
 * UniqueSet changes - some more commenting and cleanup
 *
 * Revision 1.1011  1997/03/10 17:35:57  nealk
 * More PBC
 *
 * Revision 1.1010  1997/03/10 15:42:16  nealk
 * Updating once more.
 *
 * Revision 1.1009  1997/03/07 19:20:24  nealk
 * Modified for PBC
 *
 * Revision 1.1008  1997/03/04 15:52:37  nealk
 * Modified get_FMA_cube to allow for a rectanglar region.
 *
 * Revision 1.1007  1997/02/28 20:39:45  nealk
 * Removed local_timestep.  Jim says multiple timestepping is not implemented
 * yet.
 *
 * Revision 1.1006  1997/02/28 17:47:08  nealk
 * More debugging code.
 * Re-added "local_timestep" since "fake_seq" is some weird, nearly random
 * number.
 * Move some uninitialized variables to where they become initialized.
 *
 * Revision 1.1005  1997/02/28 06:57:46  jim
 * DPMTA is now working, except for one little thing.
 * On multiple nodes, the reported energy is wrong, but the
 * trajectory appears to be dead on.  This is an odd one.
 *
 * Revision 1.1004  1997/02/27 20:01:43  nealk
 * DPMTA no longer runs every timestep.
 *
 * Revision 1.1003  1997/02/21 20:45:11  jim
 * Eliminated multiple function for switching and modified 1-4 interactions.
 * Now assumes a switching function, but parameters are such that nothing
 * happens, same for modified 1-4.  Slight penalty for rare simulations
 * in which these features are not used, but otherwise no loss and
 * simplifies code.
 *
 * Revision 1.1002  1997/02/11 18:51:41  ari
 * Modified with #ifdef DPMTA to safely eliminate DPMTA codes
 * fixed non-buffering of migration msgs
 * Migration works on multiple processors
 *
 * Revision 1.1001  1997/02/10 19:36:37  nealk
 * Added DPMTA stuff.
 *
 ***************************************************************************/

