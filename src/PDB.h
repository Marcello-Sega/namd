/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   PDB Class
     Given a PDB file name, read in all the data.
*/

#ifndef PDB_H
#define PDB_H

// These are added to the global namespace:
//   whatever is in PDBData.h
//   the typedef PDBAtomList (a singly linked list of PDBAtom *
//   the class PDB

#include "parm.h"

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
    PDB(const char *, Ambertoppar *);  // read AMBER coordinate file
    ~PDB( void);               // clear everything
    void write(const char *outfilename, const char *commentline=NULL); // write the coordinates to a file
       // the following deals only with ATOMs and HETATMs
    int num_atoms( void);

    PDBAtom *atom(int place); // get the nth atom in the PDB file
         // return linked list containing all atoms
    PDBAtomList *atoms(void ) { return atomListHead; }  
         
	// Find the extreme edges of the molecule
    void find_extremes(BigReal *min, BigReal *max, Vector rec,
                                                  BigReal frac=1.0) const;

    void set_all_positions(Vector *);	//  Reset all the positions in PDB

    void get_all_positions(Vector *);	//  Get all positions in PDB
};

#endif // PDB_H

