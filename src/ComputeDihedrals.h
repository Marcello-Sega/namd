/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEDIHEDRALS_H
#define COMPUTEDIHEDRALS_H

#include "ComputeHomeTuples.h"
#include "ReductionMgr.h"

class Molecule;

class DihedralElem {
public:
    // ComputeHomeTuples interface
    enum { size = 4 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    void computeForce(BigReal*);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void loadTuplesForAtom(void*, AtomID, Molecule*);
    static void getMoleculePointers(Molecule*, int*, int***, Dihedral**);

    // Internal data
    Index dihedralType;

  int hash() const { 
    return 0x7FFFFFFF &((atomID[0]<<24) + (atomID[1]<<16) + (atomID[2]<<8) + atomID[3]);
  }

  enum { dihedralEnergyIndex, TENSOR(virialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_DIHEDRAL_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  inline DihedralElem();
  inline DihedralElem(const Dihedral *a);
  inline DihedralElem(AtomID atom0, AtomID atom1, AtomID atom2, AtomID atom3);
  ~DihedralElem() {};

  inline int operator==(const DihedralElem &a) const;
  inline int operator<(const DihedralElem &a) const;
};

class ComputeDihedrals : public ComputeHomeTuples<DihedralElem,Dihedral>
{
public:

  ComputeDihedrals(ComputeID c) : ComputeHomeTuples<DihedralElem,Dihedral>(c) { ; }

};

#include "ComputeDihedrals.inl"

#endif

