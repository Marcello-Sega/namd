/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   PDB Class
     Given a PDB file name, read in all the data.
   As of now, you can only search the file for ATOM and 
   HETATM information.  The return value is an IntList (of new'ed
   memory so you have to delete it!) containing the list of all
   fields that match that criterion, indexed by position in the file.
   (Hence, 0 is the 1st ATOM or HETATM record, 10 is the eleventh,
   and so on...).  Note that with these searches there is no
   way to choose ATOM or HETATM; you have to make that distinguishment
   yourself.
*/

#ifndef PDB_H
#define PDB_H

// These are added to the global namespace:
//   whatever is in IntList
//   whatever is in PDBData.h
//   the typedef PDBAtomList (a singly linked list of PDBAtom *
//   the class PDB

#include "IntList.h"
#include "PDBData.h"
#include "Vector.h"
#include "Lattice.h"

typedef PDBAtom *PDBAtomPtr ;
typedef struct PAL {
  PDBAtom *data;
  struct PAL *next;
} PDBAtomList;
  
class PDB {
  private:
    PDBAtomList *atomListHead, *atomListTail;
    PDBAtom **atomArray;
      // this doesn't create a copy 
    void add_atom_element(PDBAtom *newAtom); 
    int atomCount;
    
  public:
    PDB( const char *pdbfilename);   // read in PDB from a file
    ~PDB( void);               // clear everything
    void write(const char *outfilename, const char *commentline=NULL); // write the coordinates to a file
       // the following deals only with ATOMs and HETATMs
    int num_atoms( void);

              // Ways to find an atom based on its characteristics
    IntList *find_atom_serialnumber(int );
    IntList *find_atom_name(const char *);
    IntList *find_atom_alternatelocation(const char *);
    IntList *find_atom_residuename(const char *);
    IntList *find_atom_chain(const char *); 
    IntList *find_atom_residueseq( int);
    IntList *find_atom_insertioncode( const char *);
    IntList *find_atom_segmentname( const char *);

// search on the basis of multiple characteristics; ie, "CZ","PHE",1,"PTI1"
    IntList *find_atom(const char *name=NULL, const char *residue = NULL, 
                  int residueseq = -1,  const char *segment = NULL);

    PDBAtom *atom(int place); // get the nth atom in the PDB file
         // return linked list containing all atoms
    PDBAtomList *atoms(void ) { return atomListHead; }  
         
        // find all atoms within a specific volume of space
        //  in general, an atom at coordinate x is in the region 
        //  [x1, x2) iff x1<=x<x2
    IntList *find_atoms_in_region( Real x1, Real y1, Real z1,
                                  Real x2, Real y2, Real z2, Lattice lat );

	// Find the extreme edges of the molecule
    void find_extremes(BigReal *min, BigReal *max, Vector rec,
                                                  BigReal frac=1.0) const;

    void set_all_positions(Vector *);	//  Reset all the positions in PDB

    void get_all_positions(Vector *);	//  Get all positions in PDB
};

#endif // PDB_H

