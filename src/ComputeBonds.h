//-*-c++-*-
/***************************************************************************/
/*       (C) Copyright 1996,1997 The Board of Trustees of the              */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef COMPUTEBOND_H
#define COMPUTEBOND_H

#include "common.h"
#include "NamdTypes.h"
#include "Molecule.h"
#include "ReductionMgr.h"

class TuplePatchElem;

class BondElem {
public:
    // ComputeHomeTuples interface
    enum { size = 2 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    void computeForce(BigReal*);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void loadTuplesForAtom(void*, AtomID, Molecule*);
    static void getMoleculePointers(Molecule*, int*, int***, Bond**);

    int hash() const { return 0x7FFFFFFF & ( (atomID[0]<<16) + (atomID[1])); }

    // Internal data
    Index bondType;

  enum { bondEnergyIndex, virialXIndex, virialYIndex, virialZIndex, reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_BOND_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  inline BondElem();
  inline BondElem(const Bond *a);
  inline BondElem(AtomID atom0, AtomID atom1);
  ~BondElem() {};

  inline int operator==(const BondElem &a) const;
  inline int operator<(const BondElem &a) const;
};

#include "ComputeHomeTuples.h"

class ComputeBonds : public ComputeHomeTuples<BondElem,Bond>
{
public:

  ComputeBonds(ComputeID c) : ComputeHomeTuples<BondElem,Bond>(c) { ; }

};

#include "ComputeBonds.inl"

#endif

