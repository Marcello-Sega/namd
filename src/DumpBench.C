
#include <stdio.h>
#include "DumpBench.h"
#include "SimParameters.h"
#include "ComputeNonbondedUtil.h"
#include "LJTable.h"
#include "Molecule.h"
#include "Node.h"
#include "PatchMap.h"
#include "HomePatch.h"
#include "NamdState.h"
#include "ComputeMap.h"

inline void dump_param(FILE *file, const char *name, int value) {
  fprintf(file,"%s %d\n",name,value);
}

inline void dump_param(FILE *file, const char *name, double value) {
  fprintf(file,"%s %g\n",name,value);
}

inline void dump_param(FILE *file, const char *name, Vector value) {
  fprintf(file,"%s %f %f %f\n",name,value.x,value.y,value.z);
}

int dumpbench(FILE *file) {

  Node *node = Node::Object();

  fprintf(file,"SIMPARAMETERS_BEGIN\n");

  SimParameters *simParams = node->simParameters;

#define SIMPARAM(T,N,V) dump_param(file,#N,simParams->N)
#include "DumpBenchParams.h"
#undef SIMPARAM

  fprintf(file,"SIMPARAMETERS_END\n");

  fprintf(file,"LJTABLE_BEGIN\n");

  const LJTable *ljTable = ComputeNonbondedUtil::ljTable;

  int table_dim = ljTable->get_table_dim();
  fprintf(file,"%d\n",table_dim);

  const LJTable::TableEntry *table = ljTable->get_table();
  int i,j;
  for ( i=0; i < table_dim; ++i) {
    for ( j=i; j < table_dim; ++j)
    {
      const LJTable::TableEntry *curij = &(table[2*(i*table_dim+j)]);
      fprintf(file,"%g %g %g %g\n",curij->A,curij->B,
				(curij+1)->A,(curij+1)->B);
    }
  }

  fprintf(file,"LJTABLE_END\n");

  fprintf(file,"MOLECULE_BEGIN\n");

  const Molecule *mol = node->molecule;

  fprintf(file,"%d %d\n",mol->numAtoms,mol->numCalcExclusions);

  for ( i=0; i<mol->numAtoms; ++i) {
    int vdw = mol->atomvdwtype(i);
    const ExclusionCheck *excl = mol->get_excl_check_for_atom(i);
    int min = excl->min;
    int max = excl->max;
    fprintf(file,"%d %d %d",vdw,min,max);
    if ( min <= max ) {
      int s = max - min + 1;
      const char *f = excl->flags;
      for ( int k=0; k<s; ++k ) {
        int fk = f[k];
        fprintf(file," %d",fk);
      }
    }
    fprintf(file,"\n");
  }

  fprintf(file,"MOLECULE_END\n");

  fprintf(file,"PATCHLIST_BEGIN\n");

  PatchMap *patchMap = PatchMap::Object();
  int numPatches = patchMap->numPatches();
  fprintf(file,"%d\n",numPatches);

  for ( i=0; i<numPatches; ++i) {
    HomePatch *patch = patchMap->homePatch(i);
    fprintf(file,"PATCH_BEGIN\n");
    int numAtoms = patch->getNumAtoms();
    fprintf(file,"%d\n",numAtoms);
    FullAtomList &atoms = patch->getAtomList();
    for ( j=0; j<numAtoms; ++j) {
      CompAtom &a = atoms[j];
      double x,y,z,q;
      int id,hgs,ngia,af,gf,part;
      x = a.position.x;
      y = a.position.y;
      z = a.position.z;
      q = a.charge;
      id = a.id;
      hgs = a.hydrogenGroupSize;
      ngia = a.nonbondedGroupIsAtom;
      af = a.atomFixed;
      gf = a.groupFixed;
      part = a.partition;
      fprintf(file,"%f %f %f %f %d %d %d %d %d %d\n",
        x,y,z,q,id,hgs,ngia,af,gf,part);
    }
    fprintf(file,"PATCH_END\n");
  }

  fprintf(file,"PATCHLIST_END\n");

  fprintf(file,"COMPUTEPAIR_BEGIN\n");

  ComputeMap *computeMap = ComputeMap::Object();
  int numComputes = computeMap->numComputes();
  int numPairComputes = 0;
  for ( i=0; i<numComputes; ++i) {
    if ( computeMap->type(i) == computeNonbondedPairType
         && computeMap->partition(i) == 0 ) ++numPairComputes;
  }
  fprintf(file,"%d\n",numPairComputes);
  for ( i=0; i<numComputes; ++i) {
    if ( computeMap->type(i) == computeNonbondedPairType
         && computeMap->partition(i) == 0 ) {
      int pid1 = computeMap->pid(i,0);
      int trans1 = computeMap->trans(i,0);
      int pid2 = computeMap->pid(i,1);
      int trans2 = computeMap->trans(i,1);
      fprintf(file,"%d %d %d %d\n",pid1,trans1,pid2,trans2);
    }
  }

  fprintf(file,"COMPUTEPAIR_END\n");

  return 0;
}

