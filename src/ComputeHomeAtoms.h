//-*-c++-*-
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

#ifndef COMPUTEATOM_H
#define COMPUTEATOM_H

#include "main.h"
#include "charm++.h"

#include "NamdTypes.h"

#include "Box.h"
#include "OwnerBox.h"

class Patch;

// Compute object for Atom-wise type
class ComputeHomeAtoms : public Compute
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

#endif

