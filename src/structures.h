//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *     structures.h contains useful structures
 *
 ***************************************************************************/

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "common.h"

// forward references
class IntList;


// status elements, used for Atom status 
#define UnknownAtom      0x00
#define HydrogenAtom     0x01
#define OxygenAtom       0x02
#define HBDonorAtom      0x04
#define HBAcceptorAtom   0x08
#define HBAntecedentAtom 0x10
#define HBHydrogenAtom   0x20


typedef unsigned short Index;		//  Used for index into arrays
					//  or parameters

typedef struct atom_name_info
{
	char *resname;
	char *atomname;
	char *atomtype;
} AtomNameInfo;

typedef struct atom_constants
{
	Real mass;
	Real charge;
	Index vdw_type;
	int32 status;	         // flags telling about this atom
	int32 partner;             // connecting atom, for hydrogens
	IntList *donorList;      // donor-hydrogen pairs this is part of
	IntList *acceptorList;   // acceptor-anteced pairs this is part of
	int32 hydrogenList;	// index of atom in hydrogenGroup list
} Atom;

typedef struct bond
{
	int32 atom1;
	int32 atom2;
	Index bond_type;
} Bond;

typedef struct angle
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	Index angle_type;
} Angle;

typedef struct dihedral
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	int32 atom4;
	Index dihedral_type;
} Dihedral;

typedef struct improper
{
	int32 atom1;
	int32 atom2;
	int32 atom3;
	int32 atom4;
	Index improper_type;
} Improper;

class Exclusion
{
public:
	Exclusion(void) : modified(0) {;}
	Exclusion(int a1, int a2, int mod = 0) :
		atom1(a1), atom2(a2), modified(mod) {;}
	int32 atom1;
	int32 atom2;
	Index modified;
	int hash(void) const
	{
		return atom1 + atom2;
	}
	int operator==(const Exclusion &o) const
	{
		return atom1 == o.atom1 && atom2 == o.atom2;
	}
	int operator<(const Exclusion &o) const
	{
		return
		(
		  ( atom1 < o.atom1 ) ||
		  ( atom1 == o.atom1 && atom2 < o.atom2 ) ||
		  ( atom1 == o.atom1 && atom2 == o.atom2 )
		);
	}
};

#endif

