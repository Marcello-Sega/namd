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

extern "C"
  {
  #include "dpmta.h"
  }
#include "ComputeHomePatches.h"

typedef struct patch_info
{
   int pid;                     //  Patch ID number
   int num;                     //  Number of atoms in this patch
   int *indexes;                //  Global indexes for these atoms
   Vector *pos;                 //  Positions for this patch
   struct patch_info *next;     //  Pointer to next link in list
} PatchInfo;

class ComputeDPMTA : public ComputeHomePatches {
private:
  int *slavetids;	//  PID for slave processes
  PatchInfo *patchData;	//  List containing data from patches
  PatchInfo *patchTail;	//  Tail of patch data list
  int numPatches;	//  Number of patches being dealt with
  int numDistributed;	//  Number of patches that we have
			//  distributed forces back to
  int totalAtoms;	//  Total number of atoms being dealt with
  PmtaPartInfo *fmaResults;	//  Results from the PMTA code
  PmtaPartInfo *ljResults;	//  Results from the PMTA code

  void ComputeDPMTA::get_FMA_cube(BigReal *boxsize, Vector *boxcenter);

public:
  ComputeDPMTA(ComputeID c);
  virtual ~ComputeDPMTA();
  void doWork();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeDPMTA.h,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/02/10 19:36:37 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeDPMTA.h,v $
 * Revision 1.1001  1997/02/10 19:36:37  nealk
 * Added DPMTA stuff.
 *
 * Revision 1.778  1997/01/28 01:00:33  ari
 * uplevel
 *
 * Revision 1.2  1997/01/28 01:00:11  ari
 * Adding This again
 *
 * Revision 1.1.2.1  1997/01/27 21:11:36  jim
 * test
 *
 * Revision 1.777  1997/01/17 19:35:45  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.10  1997/01/16 00:55:56  jim
 * Added reduction of energies from ComputeHomeTuples objects, except
 * for ComputeNonbondedExcl which only reports 0 energy.
 * Some problems with ReductionMgr are apparent, but it still runs.
 *
 * Revision 1.9  1997/01/14 17:59:48  jim
 * fixed multiple force additions (no adding to proxies now)
 *
 * Revision 1.8  1996/12/04 18:03:12  jim
 * added AtomProperties checkout
 *
 * Revision 1.7  1996/11/19 06:58:37  jim
 * first compiling templated version, needed ugly void* hack
 *
 * Revision 1.6  1996/11/19 04:24:24  jim
 * first templated version as ComputeHomeTuples<T>
 *
 * Revision 1.5  1996/11/18 21:28:48  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/11/04 20:06:17  nealk
 * Now it compiles :-)
 *
 * Revision 1.3  1996/11/04 19:29:02  nealk
 * Added angleForce() to system, but it is untested.
 *
 * Revision 1.2  1996/11/04 16:55:46  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/11/01 21:20:45  ari
 * Initial revision
 *
 * Revision 1.3  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/19 22:07:49  ari
 * Initial revision
 *
 * Revision 1.4  1996/07/16 01:54:12  ari
 * *** empty log message ***
 *
 * Revision 1.3  96/07/16  01:10:26  01:10:26  ari (Aritomo Shinozaki)
 * Fixed comments, added methods
 * 
 * Revision 1.2  1996/06/25 21:10:48  gursoy
 * *** empty log message ***
 *
 * Revision 1.1  1996/06/24 14:12:26  gursoy
 * Initial revision
 *
 ***************************************************************************/

