#include "ComputeConsForce.h"
#include "Molecule.h"
#include "Node.h"


ComputeConsForce::ComputeConsForce(ComputeID c, PatchID pid)
  : ComputePatch(c,pid)  {}


void ComputeConsForce::doForce(CompAtom* p, Results* r)
{ int localID,forceID;
  Molecule *molecule = Node::Object()->molecule;
  int32 *index = molecule->consForceIndexes;  // Indexes into the force array
  Vector *cf = molecule->consForce;  // Force array
  Vector *forces = r->f[Results::normal];

  for (localID=0; localID<numAtoms; ++localID)
    // When the index is -1, it means there's no constant force on this atom
    if ((forceID=index[p[localID].id]) != -1)
      forces[localID] += cf[forceID];
}
