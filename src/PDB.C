/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   PDB Class
     Code to implement the PDB class.  This reads in a bunch of
   PDBData records from a file, given the filename.  See PDB.h
   for a bit more information.
*/

#include <stdio.h>
#include <strings.h>
#include "common.h"
#include "PDB.h"
#include "SortableResizeArray.h"


// read in a file and stick all the elements on the appropriate list
PDB::PDB( const char *pdbfilename) {
  FILE *infile;
  char buf[160];

  atomCount = 0;
  atomListHead = atomListTail = NULL;
  infile = Fopen(pdbfilename, "r");
  if (! infile) {
     char s[500];
     sprintf(s, "Cannot open file '%s' for input in PDB::PDB.", pdbfilename);
     NAMD_die(s);
  }
  
    // loop through each line in the file and get a record
  while ( fgets(buf, 150, infile) ) {
   PDBData *newelement;
   char *s;
   for (s=buf; *s && *s!='\n'; s++)  // chop off the '\n'
    ;
   *s = 0;
   if ( s == (buf + 149) ) {
     char s[500];
     sprintf( s, "Line too long in pdbfile %s:\n%s\n", pdbfilename, buf);
     NAMD_die(s);
   }
   *(s+1) = 0;  // just to be on the safe side

     // I now have a string; make a PDBData element out of it
   newelement = new_PDBData(buf);
   if (!newelement) {
      NAMD_die("Could not allocate PDBData.\n");
   }
     // I only know how to deal with ATOM and HETATM types; and
     //  I want to throw away the unknown data types; so
   if (newelement -> type() != PDBData::ATOM && 
           newelement -> type() != PDBData::HETATM) {
       delete newelement;
   } else {
       add_atom_element( (PDBAtom *) newelement);
   }
  }  // input while loop
 Fclose(infile);
 
 // now I have a linked list, and I know the size.  However,
 // that's got pretty slow access time for ramdom fetches, so
 // I'll turn it into an array
 {
  atomArray = new PDBAtomPtr[atomCount];
  if ( atomArray == NULL )
  {
    NAMD_die("memory allocation failed in PDB::PDB");
  }
  PDBAtomList *tmp = atomListHead;
  int i=0;                              // just need to copy the pointers
  for (i=0, tmp = atomListHead; tmp != NULL; tmp = tmp -> next, i++) {
    atomArray[i] = tmp -> data;
  }
     // now delete the linked list (w/o deleting the data)
  PDBAtomList *tmp2;
  for (tmp2 = tmp = atomListHead; tmp != NULL; tmp = tmp2) {
    tmp2 = tmp->next;
    delete tmp;
  }
  atomListHead = atomListTail = NULL;
 }  // everything converted
 
}

//  Destructor - delete all the data pointed to by the array
//   and then delete the array
PDB::~PDB( void )
{
	int i;
	for (i=atomCount-1; i>=0; i--)
	   delete atomArray[i];
	delete [] atomArray;
	atomArray = NULL;
	atomCount = 0;
}

// print the PDB file out to a given file name
void PDB::write(const char *outfilename, const char *commentline)
{
	int i;
	char s[200];
	FILE *outfile;
	if ((outfile = fopen(outfilename, "w")) == NULL) {
	   sprintf(s, "Cannot open file '%s' in PDB::write.", outfilename);
	   NAMD_die(s);
	}

	if (commentline != NULL)
	{
		sprintf(s, "REMARK  %s\n", commentline);
		if (fputs(s, outfile) == EOF)
		{
			NAMD_die("EOF in PDB::write writing the comment line - file system full?");
		}
	}

	for (i=0; i<atomCount; i++){ // I only contain ATOM/HETATM records
	  atomArray[i]->sprint(s, PDBData::COLUMNS);
	  if ( (fputs(s, outfile)    == EOF) || 
	       (fputc('\n', outfile) == EOF)    ) {
	    sprintf(s, "EOF in PDB::write line %d - file system full?", i);
	    NAMD_die(s);
	  }
	}
	if (fputs("END\n", outfile) == EOF) {
	   NAMD_die("EOF in PDB::write while printing 'END' -- file system full?");
	}
	if (fclose(outfile) == EOF) {
	   NAMD_die("EOF in PDB::write while closing -- file system full?");
	}
	  
}

// store the info on the linked list
void PDB::add_atom_element( PDBAtom *newAtom)
{
  PDBAtomList *tmp = new PDBAtomList;
  if ( tmp == NULL )
  {
    NAMD_die("memory allocation failed in PDB::add_atom_element");
  }
  tmp -> data = newAtom;
  tmp -> next = NULL;
  
  if (atomListHead == NULL) {        // make the list
    atomListHead = atomListTail = tmp;
  } else {
    atomListTail -> next = tmp;       // add to the tail
    atomListTail = tmp;
  }
  atomCount++;
}


// return the number of atoms found
int PDB::num_atoms( void)
{
  return atomCount;
}

// A macro definition to search through the array for
//   a bit o' data (this macro not for public consumption)
#define find_atom_MACRO(fctnname,vartype,var,comparison)                 \
IntList *PDB::fctnname( vartype var)                                     \
{                                                                        \
  IntList *ilist = new IntList;  /* Basically, I'm going to go through */\
  if ( ilist == NULL )                                                   \
  {                                                                      \
    NAMD_die("memory allocation failed in PDB::fctnname");               \
  }                                                                      \
  PDBAtomPtr *atomPtr;           /* the array looking for              */\
  PDBAtom *atom;                                                         \
  int i;                         /* when "comparison" is true.         */\
  for (i=0, atomPtr = atomArray;  atom = *atomPtr, i<atomCount;          \
                                  atomPtr ++, i++ ) {                    \
     if (comparison)                                                     \
       ilist -> add(i);          /* and keep track of where that was   */\
  }                                                                      \
  return ilist;                                                          \
}                                                                        \

// find an atom on the basis of it serial number
find_atom_MACRO(find_atom_serialnumber, int, serial,
   atom->serialnumber() == serial)
   
// find an atom on the basis of it name
find_atom_MACRO(find_atom_name, const char *, name,
  !strcmp(atom->name(), name) )

// lookup atom on alt. location
find_atom_MACRO(find_atom_alternatelocation, const char *, alt,
  !strcmp(atom->alternatelocation(), alt))
  
// lookup atom on residuename
find_atom_MACRO( find_atom_residuename, const char *, resname,
  !strcmp(atom->residuename(), resname) )

// lookup atom on chain
find_atom_MACRO( find_atom_chain, const char *, chain,
  !strcmp(atom->chain(), chain) )

// lookup atom on residue seq.
find_atom_MACRO( find_atom_residueseq, int, reseq,
  atom->residueseq() == reseq )

// lookup atom on insertion code
find_atom_MACRO( find_atom_insertioncode, const char *, code,
  atom->insertioncode() == code )


// lookup atom on segment name
find_atom_MACRO( find_atom_segmentname, const char *, segname,
  !strcmp(atom->segmentname(), segname) )

// search looking for a match like "CA PRO 34 HRV2"
IntList *PDB::find_atom( const char *name, const char *residue, int reseq, 
                    const char *segment)
{
  IntList *ilist = new IntList;
  if ( ilist == NULL )
  {
    NAMD_die("memory allocation failed in PDB::find_atom");
  }
  PDBAtomPtr *atomptr;
  int i;
  for (i=0, atomptr=atomArray;  i<atomCount; atomptr ++, i++ ) {
     if (
 (!name       || !strcmp(name, (*atomptr)->name())  )      &&
 (!residue    || !strcmp(residue,(*atomptr)->residuename()))&&
 (reseq != -1 ||  (*atomptr)->residueseq() == reseq )     &&
 (!segment    || !strcmp(segment,(*atomptr)->segmentname()))
        )
    ilist -> add(i);
  }
  return ilist;
}

// Find all the atoms in a certain region of space
//  in general, an atom at coordinate x is in the region 
//  [x1, x2) iff x1<=x<x2
IntList *PDB::find_atoms_in_region( Real x1, Real y1, Real z1,
                              Real x2, Real y2, Real z2, Lattice lat )
{
  IntList *ilist = new IntList;
  if ( ilist == NULL )
  {
    NAMD_die("memory allocation failed in PDB::find_atoms_in_region");
  }
  PDBAtomPtr *atomptr;
  Real swaptmp;
  if (x1>x2) { swaptmp=x2; x2=x1; x1=swaptmp; }  // insure x1<x2
  if (y1>y2) { swaptmp=y2; y2=y1; y1=swaptmp; }  // make y1<y2
  if (z1>z2) { swaptmp=z2; z2=z1; z1=swaptmp; }  // ditto for z2 and z2
  int i;
  Vector center(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2));
  for (i=0, atomptr=atomArray; i<atomCount; atomptr++, i++) {
    Vector atom((*atomptr)->xcoor(),(*atomptr)->ycoor(),(*atomptr)->zcoor());
    atom = lat.nearest(atom,center);
    if ( x1 <= atom.x && atom.x < x2 &&
         y1 <= atom.y && atom.y < y2 &&
         z1 <= atom.z && atom.z < z2     ) {
      ilist -> add(i);
    }
  }
  return ilist;
}

// Reset all the atom positions.  This is used in preparation for
// output in cases like the restart files, etc.
void PDB::set_all_positions(Vector *pos)
{
	int i;
	PDBAtomPtr *atomptr;

	for (i=0, atomptr=atomArray; i<atomCount; atomptr++, i++)
	{
		(*atomptr)->xcoor(pos[i].x);
		(*atomptr)->ycoor(pos[i].y);
		(*atomptr)->zcoor(pos[i].z);
	}
}

//  Get all the atom positions into a list of Vectors
void PDB::get_all_positions(Vector *pos)
{
	int i;
	PDBAtomPtr *atomptr;

	for (i=0, atomptr=atomArray; i<atomCount; atomptr++, i++)
	{
		pos[i].x = (*atomptr)->xcoor();
		pos[i].y = (*atomptr)->ycoor();
		pos[i].z = (*atomptr)->zcoor();
	}
}

//  given an index, return that atom
PDBAtom *PDB::atom(int place)
{
  if (place <0 || place >= atomCount)
    return NULL;
  return atomArray[place];
}


// find the lowest and highest bounds for a fraction of the atoms
void PDB::find_extremes(BigReal *min, BigReal *max, Vector rec, BigReal frac) const
{
    SortableResizeArray<Real> coor;
    coor.resize(atomCount);
    SortableResizeArray<Real>::iterator c_i = coor.begin();
    PDBAtomPtr *atomptr = atomArray;
    for (int i=0; i<atomCount; ++i, ++atomptr) {
      PDBAtom *atom = *atomptr;
      Vector pos(atom->xcoor(),atom->ycoor(),atom->zcoor());
      c_i[i] = rec*pos;
    }
    coor.sort();
    int ilow = (1.0 - frac) * atomCount;
    if ( ilow < 0 ) ilow = 0;
    if ( ilow > atomCount/2 ) ilow = atomCount/2;
    *min = coor[ilow];
    int ihigh = atomCount - ilow - 1;
    *max = coor[ihigh];
}

//#define TEST_PDB_CLASS
#ifdef TEST_PDB_CLASS

main()
{
 IntList *ilist;
 PDB *pdb = new PDB("pti.pdb");
 if ( atomArray == NULL )
 {
   NAMD_die("memory allocation failed in main of test PDB class");
 }
 ilist = pdb->find_atom_name("CA");
 if (!ilist)
   printf("None found.\n");
 else {
   int i;
   char s[200];
   for (i=0; i<ilist->num(); i++) {
     pdb->atom((*ilist)[i]) -> sprint(s);
     printf("%s\n", s);
   }
   delete ilist;
 } // if ilist
 
 printf("Now, search through space.\n");
 
 ilist = pdb->find_atoms_in_region(4.38, 19.5, 3.0,
                                   4.40, 20.0, 3.2 );
 if (!ilist)
   printf("None found.\n");
 else {
   int i;
   char s[200];
   printf("%d found\n", ilist -> num());
   for (i=0; i<ilist->num(); i++) {
     pdb->atom((*ilist)[i]) -> sprint(s);
     printf("%s\n", s);
   }
   delete ilist;
 }
}
#endif // TEST_PDB_CLASS

