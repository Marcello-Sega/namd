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

#ifndef COMPUTEIMPROPERS_H
#define COMPUTEIMPROPERS_H

#include "ComputeHomeTuples.h"
#include "ReductionMgr.h"

class Molecule;

class ImproperElem {
public:
    // ComputeHomeTuples interface
    enum { size = 4 };
    AtomID atomID[size];
    int    localIndex[size];
    TuplePatchElem *p[size];
    void computeForce(BigReal*);
    // The following is evil, but the compiler chokes otherwise. (JCP)
    static void loadTuplesForAtom(void*, AtomID, Molecule*);
    static void getMoleculePointers(Molecule*, int*, int***, Improper**);

    // Internal data
    Index improperType;

  int hash() const {
    return 0x7FFFFFFF &((atomID[0]<<24) + (atomID[1]<<16) + (atomID[2]<<8) + atomID[3]);
  }
  enum { improperEnergyIndex, virialXIndex, virialYIndex, virialZIndex, reductionDataSize };
  enum { reductionChecksumLabel = REDUCTION_IMPROPER_CHECKSUM };
  static void submitReductionData(BigReal*,SubmitReduction*);

  ImproperElem() {
	atomID[0] = -1;
	atomID[1] = -1;
	atomID[2] = -1;
	atomID[3] = -1;
	p[0] = NULL;
	p[1] = NULL;
	p[2] = NULL;
	p[3] = NULL;
  }
  ImproperElem(const Improper *a) {
    atomID[0] = a->atom1;
    atomID[1] = a->atom2;
    atomID[2] = a->atom3;
    atomID[3] = a->atom4;
    improperType = a->improper_type;
  }

  ImproperElem(AtomID atom0, AtomID atom1, AtomID atom2, AtomID atom3) {
    if (atom0 > atom3) {  // Swap end atoms so lowest is first!
      AtomID tmp = atom3; atom3 = atom0; atom0 = tmp; 
      tmp = atom1; atom1 = atom2; atom2 = tmp;
    }
    atomID[0] = atom0;
    atomID[1] = atom1;
    atomID[2] = atom2;
    atomID[3] = atom3;
  }
  ~ImproperElem() {};

  int operator==(const ImproperElem &a) const {
    return (a.atomID[0] == atomID[0] && a.atomID[1] == atomID[1] &&
        a.atomID[2] == atomID[2] && a.atomID[3] == atomID[3]);
  }

  int operator<(const ImproperElem &a) const {
    return  (atomID[0] < a.atomID[0] ||
            (atomID[0] == a.atomID[0] &&
            (atomID[1] < a.atomID[1] ||
            (atomID[1] == a.atomID[1] &&
            (atomID[2] < a.atomID[2] ||
            (atomID[2] == a.atomID[2] &&
             atomID[3] < a.atomID[3] 
	     ))))));
  }
};

class ComputeImpropers : public ComputeHomeTuples<ImproperElem,Improper>
{
public:

  ComputeImpropers(ComputeID c) : ComputeHomeTuples<ImproperElem,Improper>(c) { ; }

};

#endif

