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

#ifndef COMPUTE_H
#define COMPUTE_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"

class Patch;

// Message to charmm runtime to for <linkto class=Compute>Compute</linkto>
// object to schedule work for itself.
class LocalWorkMsg : public comm_object
{
private:
      void *Calcfptr;

public:
      LocalWorkMsg();
      ~LocalWorkMsg();
};



// Base class for various forms of Compute objects
// including: <linkto class=ComputeAtoms>ComputeAtoms</linkto> 
// and <linkto class=ComputePatches>ComputePatches</linkto>
class Compute {
   private:
      int numPatches;
      int counter;

   protected :
	  // For use by derived class objects to queue up itself for work
      void enqueueWork() {}
	  // Utility which sets up list of home patches to current Node
	  // and register to the patches.
      void setupHomePatches() {};
	  // Utility which sets up list of given patches and 
	  // registers to the patches.
      void setupPatches(PatchIDList &) {};

   public:
      Compute() {};
      ~Compute() {};

	  // Signal from patch or proxy that data is ready.
	  // When all Patches and Proxies needed by this Compute object
	  // have checked-in, we are ready to enqueueWork()
      void patchReady(PatchID) {};
      void patchReady(void) {};

	  // Actual computation, called by work queue, triggered
	  // by scheduled work msg.
      void virtual doWork(LocalWorkMsg *msg) {};
};


/*
// Compute object for Atom-wise type forces i.e. Bonded Forces in general.
// 
class ComputeAtoms : public Compute
{
	// Array of atoms for which we must report forces.
    AtomArray atom;
	// Array of patches for which we must report all appropo
	// bonded-atom type forces.
    PatchArray patch;
    void (* mapMethod)();
    void mapByAtoms();
    void mapByPatches();
    void mapByHomePatches();
    
    protected:
	    // Take list of atoms and registers to appropriate patches
	    // which are not necessarily the same as atoms
	    // for which forces are computed (generally additional atoms
	    // are needed for computation of forces)
	  registerDependency(AtomArray &);

    public:
	  // void argument constructor implies this object will
	  // do Atom type computation appropriate atoms in Home Patches.
       ComputeAtoms();
	  // Initializes this object to compute appropriate forces
	  // for the atoms in AtomArray.
       ComputeAtoms(PatchArray &);
       ComputeAtoms(AtomArray &);
       ~ComputeAtoms(); 

	  // Signal to take known atoms and register 
	  // dependencies to appropriate patches.
       void mapAtoms();
	  // work method
       void virtual doWork(LocalWorkMsg *msg);
};



class ComputePatches : public Compute
{
    public:
	 // Default constructor imples do computation for all atoms in
	 // a given patch.
       ComputePatches();
	// do computation over all atoms in 
       ComputePatches(PatchArray &);

       ~ComputePatches();

	  // work method
       void virtual doWork(LocalWorkMsg *msg);
};
*/

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Compute.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/19 22:07:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Compute.h,v $
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

