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
#include <stdlib.h>
#ifndef WIN32
#include <strings.h>
#endif
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


// find the lowest and highest bounds based on a fraction of the atoms
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
    int ihigh = atomCount - ilow - 1;
    BigReal span = coor[ihigh] - coor[ilow];
    BigReal extension = (1.0 - frac) * span / (2.0 * frac - 1.0);
    *max = coor[ihigh] + extension;
    *min = coor[ilow] - extension;
}

//#define TEST_PDB_CLASS
#ifdef TEST_PDB_CLASS

main()
{
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



// This function was borrowed from VMD code in "ReadPARM.C".
// It skips to a new line.
static int readtoeoln(FILE *f) {
  char c;

  /* skip to eoln */
  while((c = getc(f)) != '\n') {
    if (c == EOF) 
      return -1;
  }

  return 0;
}  



// read in an AMBER coordinate file and populate the PDB structure
PDB::PDB( const char *filename, Ambertoppar *amber_data)
{ int i,j,k;
  Real coor[3];
  char buf[13],resname[5],atomname[5];
  FILE *infile;
  PDBAtom *pdb;

  if ((infile=Fopen(filename, "r")) == NULL)
    NAMD_die("Can't open AMBER coordinate file!");

  readtoeoln(infile);  // Skip the first line (title)

  fscanf(infile,"%d",&atomCount);  // Read in num of atoms
  if (atomCount != amber_data->Natom)
    NAMD_die("Num of atoms in coordinate file is different from that in parm file!");
  readtoeoln(infile);

  atomArray = new PDBAtomPtr[atomCount];
  if ( atomArray == NULL )
  {
    NAMD_die("memory allocation failed in PDB::PDB");
  }
  
  // Read in the coordinates, which are in the format of 6F12.7
  // Other fields are copied from "amber_data"
  for (i=0; i<atomCount; ++i)
  { // Read x,y,z coordinates
    for (j=0; j<3; ++j)
    { for (k=0; k<12; ++k)
      { buf[k]=getc(infile);
        if (buf[k]=='\n' || buf[k]=='\0' || buf[k]==EOF)
          NAMD_die("Error reading AMBER coordinate file!");
      }
      buf[12] = '\0';
      coor[j] = atof(buf);
    }
    if (i%2 == 1)
      readtoeoln(infile);
    // Copy name, resname and resid from "amber_data"
    for (j=0; j<4; ++j)
    { resname[j] = amber_data->ResNames[amber_data->AtomRes[i]*4+j];
      atomname[j] = amber_data->AtomNames[i*4+j];
    }
    resname[4] = atomname[4] = '\0';
    // Create a new PDB record, and fill in its entries
    pdb = new PDBAtomRecord("");
    pdb->name(atomname);
    pdb->residuename(resname);
    pdb->serialnumber(i+1);
    pdb->residueseq(amber_data->AtomRes[i]+1);
    pdb->coordinates(coor);
    atomArray[i] = pdb;  // Include the new record into the array
  }
}

#define LINESIZE 100

/* This constructor initializes the PDB data using a Gromacs
   coordinate file, generating an error message if the file
   can't be parsed or if its contents don't jive with what is in
   the topo file <topology>. */
PDB::PDB(const char *filename, const GromacsTopFile *topology) {
  int i;
  char buf[LINESIZE];
  FILE *infile;
  
  /* open up the coordinate file */
  infile=Fopen(filename, "r");
  if (infile == NULL)
    NAMD_die("Can't open GROMACS coordinate file!");

  fgets(buf,LINESIZE-1,infile); // get the title
  if(strcmp(buf,topology->getSystemName()) != 0)
    NAMD_die("System names in topology and coordinate files differ.");

  fgets(buf,LINESIZE-1,infile); // get the number of atoms
  sscanf(buf,"%d",&atomCount);
  if (atomCount != topology->getNumAtoms())
    NAMD_die("Num of atoms in coordinate file is different from that in topology file!");

  /* read in the atoms */
  atomArray = new PDBAtomPtr[atomCount];
  if ( atomArray == NULL )
    NAMD_die("memory allocation failed in PDB::PDB");

  for (i=0;i<atomCount;i++) {
    char *buf2, resname[11], atomname[11], atmtype[11];
    int resnum, typenum;
    Real charge,mass,coor[3];
    PDBAtom *pdb = new PDBAtomRecord("");  
    
    fgets(buf,LINESIZE-1,infile); // get a line
    buf2 = buf+20; // skip three fields to get to the coordinates
    if(3 != sscanf(buf2,"%f%f%f",
		   &coor[0],&coor[1],&coor[2]))
      NAMD_die("Couldn't get three coordinates from file.");
    topology->getAtom(i,&resnum,resname,
		      atomname,atmtype,&typenum,&charge,&mass);
    coor[0] *= 10; // convert to angstroms from nanometers
    coor[1] *= 10;
    coor[2] *= 10;
    
    pdb->name(atomname);
    pdb->residuename(resname);
    pdb->serialnumber(i+1);
    pdb->residueseq(resnum+1);
    pdb->coordinates(coor);
    
    atomArray[i] = pdb;  // Include the new record into the array
  }
}

