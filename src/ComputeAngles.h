/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEANGLE_H
#define COMPUTEANGLE_H

#include "ComputeHomeTuples.h"
#include "ReductionMgr.h"

class Molecule;

class AngleElem {
public:
    // ComputeHomeTuples interface
    enum { size = 3 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    void computeForce(BigReal*);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void loadTuplesForAtom(void*, AtomID, Molecule*);
    static void getMoleculePointers(Molecule*, int*, int***, Angle**);

    // Internal data
    Index angleType;

  int hash() const { 
    return 0x7FFFFFFF & ((atomID[0]<<22) + (atomID[1]<<11) + (atomID[2])); 
  }

  enum { angleEnergyIndex, TENSOR(virialIndex), reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_ANGLE_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  inline AngleElem();
  inline AngleElem(const Angle *a);
  inline AngleElem(AtomID atom0, AtomID atom1, AtomID atom2);
  ~AngleElem() { };

  inline int operator==(const AngleElem &a) const;
  inline int operator<(const AngleElem &a) const;
};

class ComputeAngles : public ComputeHomeTuples<AngleElem,Angle>
{
public:

  ComputeAngles(ComputeID c) : ComputeHomeTuples<AngleElem,Angle>(c) { ; }

};

#include "ComputeAngles.inl"

#endif

