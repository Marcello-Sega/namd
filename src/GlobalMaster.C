/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "Node.h"
#include "Molecule.h"
#include "NamdTypes.h"
#include "GlobalMaster.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

void GlobalMaster::processData(AtomIDList::iterator a_i,
			       AtomIDList::iterator a_e,
			       PositionList::iterator p_i,
			       PositionList::iterator g_i,
			       PositionList::iterator g_e,
			       AtomIDList::iterator last_atoms_forced_i,
			       AtomIDList::iterator last_atoms_forced_e,
			       ForceList::iterator last_forces_i) {
  atomIdBegin = a_i;
  atomIdEnd = a_e;
  atomPositionBegin = p_i;
  groupPositionBegin = g_i;
  groupPositionEnd = g_e;
  lastAtomsForcedBegin = last_atoms_forced_i;
  lastAtomsForcedEnd = last_atoms_forced_e;
  lastForcesBegin = last_forces_i;

  calculate();

  /* check to make sure the force arrays still match */
  if(appForcesChanged) {
    if(fAtoms.size() != appForces.size())
      NAMD_die("# of atoms forced != # of forces given");
  }
  if(appForcesChanged) {
    if(reqGroups.size() != grpForces.size())
      NAMD_die("# of groups forced != # of groups requested");
  }

  /* check if the groups have changed, and update the masses */
  if(reqGroupsChanged) {
    groupMasses.resize(0);

    // is there a way to do this non-globally?
    Molecule *mol = Node::Object()->molecule;

    // update each group mass
    int g;
    for(g=0;g<reqGroups.size();g++) {
      AtomIDList &atom_list = reqGroups[g];
      BigReal mass = 0; // the total mass of the group
      int a;
      for(a=0;a<atom_list.size();a++) { // get the total
	mass += mol->atommass(a);
      }
      groupMasses.add(mass); // add the mass to the group
    }
  }
}

void GlobalMaster::check() const {
  /* check to make sure the force arrays still match */
  if(fAtoms.size() != appForces.size())
    NAMD_die("# of atoms forced != # of forces given");
  if(reqGroups.size() != grpForces.size())
    NAMD_die("# of groups forced != # of groups requested");
}

void GlobalMaster::clearChanged() {
  reqAtomsChanged = false;
  appForcesChanged = false;
  reqGroupsChanged = false;
}

void GlobalMaster::calculate() {
  NAMD_die("Internal error: pure virtual function called");
}

GlobalMaster::GlobalMaster() {
  clearChanged();
}

bool GlobalMaster::changedAtoms() {
  return reqAtomsChanged;
}

bool GlobalMaster::changedForces() {
  return appForcesChanged;
}

bool GlobalMaster::changedGroups() {
  return reqGroupsChanged;
}

const AtomIDList &GlobalMaster::requestedAtoms() {
  return reqAtoms;
}

AtomIDList &GlobalMaster::modifyRequestedAtoms() {
  reqAtomsChanged = true;
  return reqAtoms;
}

const AtomIDList &GlobalMaster::forcedAtoms() {
  return fAtoms;
}

const ForceList &GlobalMaster::appliedForces() {
  return appForces;
}

const ForceList &GlobalMaster::groupForces() {
  return grpForces;
}

const ResizeArray<AtomIDList> &GlobalMaster::requestedGroups() {
  return reqGroups;
}

AtomIDList &GlobalMaster::modifyForcedAtoms() {
  appForcesChanged = true;
  return fAtoms;
}

ForceList &GlobalMaster::modifyAppliedForces() {
  appForcesChanged = true;
  return appForces;
}

ForceList &GlobalMaster::modifyGroupForces() {
  // XXX should we mark something else here?
  appForcesChanged = true;
  return grpForces;
}

ResizeArray<AtomIDList> &GlobalMaster::modifyRequestedGroups() {
  reqGroupsChanged = true;
  DebugM(1,"Groups have changed.\n");
  return reqGroups;
}

const AtomIDList::iterator GlobalMaster::getAtomIdBegin() {
  return atomIdBegin;
}

const AtomIDList::iterator GlobalMaster::getAtomIdEnd() {
  return atomIdEnd;
}

const PositionList::iterator GlobalMaster::getAtomPositionBegin() {
  return atomPositionBegin;
}

const PositionList::iterator GlobalMaster::getGroupPositionBegin() {
  return groupPositionBegin;
}

const PositionList::iterator GlobalMaster::getGroupPositionEnd() {
  return groupPositionEnd;
}

const ResizeArray<BigReal>::iterator GlobalMaster::getGroupMassBegin()
{
  return groupMasses.begin();
}

const ResizeArray<BigReal>::iterator GlobalMaster::getGroupMassEnd() {
  return groupMasses.end();
}

const AtomIDList::iterator GlobalMaster::getLastAtomsForcedBegin() {
  return lastAtomsForcedBegin;
}

const AtomIDList::iterator GlobalMaster::getLastAtomsForcedEnd() {
  return lastAtomsForcedEnd;
}

const ForceList::iterator GlobalMaster::getLastForcesBegin() {
  return lastForcesBegin;
}
