//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PDB.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.778 $	$Date: 1997/01/28 00:31:02 $
 *
 ***************************************************************************
 * DESCRIPTION:
 * PDB Class
 *   Given a PDB file name, read in all the data.
 * As of now, you can only search the file for ATOM and 
 * HETATM information.  The return value is an IntList (of new'ed
 * memory so you have to delete it!) containing the list of all
 * fields that match that criterion, indexed by position in the file.
 * (Hence, 0 is the 1st ATOM or HETATM record, 10 is the eleventh,
 * and so on...).  Note that with these searches there is no
 * way to choose ATOM or HETATM; you have to make that distinguishment
 * yourself.
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PDB.h,v $
 * Revision 1.778  1997/01/28 00:31:02  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:29  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:37  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.8  1995/07/07 20:06:30  nelson
 * Removed annoying error message in PDBAtomList declaration
 *
 * Revision 1.7  95/03/08  14:47:50  14:47:50  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.6  94/10/12  15:25:40  15:25:40  nelson (Mark T. Nelson)
 * Added get_all_positions, set_all_positions and added comment line
 * to write()
 * 
 * Revision 1.5  94/10/05  17:05:03  17:05:03  dalke (Andrew Dalke)
 * Added a 'write' function to write the PDB coordinates
 * 
 * Revision 1.4  94/09/07  17:16:14  17:16:14  dalke (Andrew Dalke)
 * Converted from linked list to array after reading in PDB
 * 
 * Revision 1.3  94/08/08  15:32:39  15:32:39  nelson (Mark T. Nelson)
 * added find_extremes routine
 * 
 * Revision 1.2  94/07/06  01:31:55  01:31:55  dalke (Andrew Dalke)
 * Added some more "find_atom_*" functions
 * 
 * Revision 1.1  94/07/05  12:52:16  12:52:16  dalke (Andrew Dalke)
 * Initial revision
 * 
 ***************************************************************************/


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
                                  Real x2, Real y2, Real z2 );

	// Find the extreme edges of the molecule
    void find_extremes(Vector *, Vector *);

    void set_all_positions(Vector *);	//  Reset all the positions in PDB

    void get_all_positions(Vector *);	//  Get all positions in PDB
};

#endif // PDB_H
