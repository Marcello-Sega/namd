/***************************************************************************/
/*     (C) Copyright 1995,1996,1997 The Board of Trustees of the           */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *   The code to implement the various PDBData classes.  These are:
 *  PDBData (base class), PDBUnknown (for unimplemented records),
 *  and PDBAtom, PDBAtomRecord, and PDBHetAtm, for the various types
 *  of atom records.
 *
 ***************************************************************************/

// Here are the routines to manupulate a PDB ATOM record
// It can be created by hand or with a string in the PDB format

#include <stdio.h>    // sprintf and sscanf
#include <strings.h>  // strncpy
#include <stdlib.h>   // atoi and atof
#include "PDBData.h"

//  Define some constants for the class
const int PDBAtom::default_serial = -1;
const int PDBAtom::default_residueseq = -1;
const Real PDBAtom::default_coor = 9999.000; 
const Real PDBAtom::default_occupancy = 1.00;
const Real PDBAtom::default_temperaturefactor = 0.00;
const int PDBAtom::no_footnote = 0;

// write down the names so I won't have to do so again
const char *PDBData::PDBNames[PDBData::UNKNOWN+1] = {
       "HEADER", "OBSLTE", "COMPND",
       "SOURCE", "EXPDTA", "AUTHOR", "REVDAT", "SPRSDE", "JRNL",
       "REMARK", "SEQRES", "FTNOTE", "HET", "FORMUL", "HELIX",
       "SHEET", "TURN", "SSBOND", "SITE", "CRYST1", "ORIGX",
       "SCALE", "MTRIX", "TVECT", "MODEL", "ATOM", "HETATM",
       "SIGATM", "ANISOU", "SIGUIJ", "TER", "ENDMDL", "CONECT",
       "MASTER", "END", "UNKNOWN"};

// Parse the input, a char *, for an integer.
//   The input string is "data" and is of length 'length'.
//   The integer starts at position "start" (ASSUMES 1st character
// is at location 1 !*!*!) and is at most "size" characters long.
//   If the string is not long enough, return the default integer value.
//Example: scan("test 12345", 10, 5, 2, &tempint, 99);
//  would set tempint to 12.
//NOTE:  a blank string ("     ") will return the value 0.  If
//   there is more than one integer in the string, it will return
//   the first one seen and not complain.
void PDBData::scan( const char *data, int length, int start, int size, 
                       int *ans, int defalt)
{
  char tempbuffer[200];       // temp. string buffer
  if (length < start) {       // check if the string is long enough
    *ans = defalt;                   // return the default
    return;
  }
  if (size>199)                      // make sure I won't overflow my array
    size=199;
  strncpy(tempbuffer, data + start-1, size);  // convert the string to an int
  tempbuffer[size]= 0;
  int flg=0;
  for (int i=strlen(tempbuffer)-1; i>=0; i--) { // see if this is a blank string
    if (tempbuffer[i]>' ') {
       flg = 1;
       break;
    }
  }
  if (flg != 1) {  // then it was a blank string
    *ans = defalt;
  } else {
  *ans = atoi(tempbuffer);
  }
}

//  Parse the input for a string; as above, but there is no default string.
// I assume that "*ans" has at least 'size+1' characters.
void PDBData::scan( const char *data, int length, int start, int size, char *ans)
{
  if (length < start) {             // check if the string is long enough
   ans[0] = 0;
   return;
  }
  strncpy(ans, data + start - 1, size);
  ans[size]=0;
                      //  fix up leading and trailing spaces
  int i,j;
  for (i=0,j=0; ans[i]; i++)  //   all this does is strip out _all_ the
    if (ans[i]!=' ' && ans[i]!='\t')   //   spaces -- this is important because I
      ans[j++]=ans[i];  // check that a string is empty by looking to see if 
  ans[j]=0;   // [0] == 0 instead of checking to see if all the elements are spaces
}
//  Parse the input for a Real
void PDBData::scan( const char *data, int length, int start, int size,
                        Real *ans, Real defalt)
{
  char tempbuffer[200];
  if (length < start) {               // check if the string is long enough
   *ans = defalt;                     // return the default
   return;
  }
  if (size>199)                       // make sure I won't overflow my array
    size=199;
  strncpy(tempbuffer, data + start - 1, size);// convert the string to a Real
  tempbuffer[size]= 0;
  int flg=0;
  for (int i=strlen(tempbuffer)-1; i>=0; i--) { // see if this is a blank string
    if (tempbuffer[i]>' ') {
       flg = 1;
       break;
    }
  }
  if (flg != 1) {  // then it was a blank string
    *ans = defalt;
  } else {
    *ans = atof(tempbuffer);  // WARNING : ASSUMES Real <= double!!!
  }
}

/////****/
// Parse the input data looking for field number 'N'
// where a field is defined as a collection of
// non-whitespace (ws == {space, tab}) characters
// and each field is seperated by one or more whitespace characters
// The result is either
//  1)  field 'N' if it exists
//  2)  some string starting with '#' if it does not
// Note: this means that you can use # to mark fields which don't exist
// Also note that this ASSUMES the first field is field 1 !*!*!*
void PDBData::field(const char *data, int fld, char *result)
{
  int i;

  Bool onword = FALSE;
  if (fld<=0)  {     // ask a stupid question, get a stupid answer
    result[0]='#';
    result[1]=0;
    return;
  }
  for (i=0; data[i]; i++)
    if (!onword && data[i] != ' ' && data[i] != '\t')  {  // if I found a field
       onword = TRUE;     // mark that I'm on it
       if (--fld <= 0)    // am I done?
         break;
    } else {
      if (onword && (data[i] == ' ' || data[i] == '\t')) { // left a field
        onword = FALSE;   // mark that I left
      }
    }
  if (fld>0) {  // oh no, didn't find the field!
     result[0] = '#';
     result[1] = 0;
     return;
  }
  
  int cpy=0;  // copy the field to the output
  while (data[i] != ' ' && data[i] != '\t' && data[i])
     result[cpy++] = data[i++];
  result[cpy] = 0;  // terminate and I'm done
}


///*** provide a simple way to print to a specific column
// Note that this ASSUMES the first character is column 1 !*!*!*
// print an integer
void PDBData::sprintcol( char *s, int start, int len, int val)
{
 char temps[100];
 sprintf(temps, "%*d", len, val);  // convert the int to a string
 sprintcol( s, start, len, temps); // copy to the output string
}
// print a Real
void PDBData::sprintcol( char *s, int start, int len, int prec, Real val)
{
 char temps[100];
 sprintf(temps, "%*.*f", len, prec, val);
 sprintcol( s, start, len, temps); // copy to the output string
}
// print a string
void PDBData::sprintcol( char *s, int start, int len, const char *val)
{
 s+=start-1;
 while (len-- >0 && *val)  // copy string up to end of string or len
  *s++ = *val++;
}

/*********************************************************/
//   base class for both PDB ATOM and HETATM      //
/*********************************************************/
void PDBAtom::parse(const char *data)
{
 char tempstr[100];
 field(data, 1, tempstr);  // get info about field #1 (the first one)
 if (tempstr[0] == '#')  {
   parse_field_data(data);
 } else {
   parse_column_data(data);
 } 
}
PDBAtom::PDBAtom( const char *data, PDBPossibleAtoms atomclass)
	: PDBData( atomclass==USE_HETATM ? PDBData::HETATM : PDBData::ATOM)
{
 parse(data);
};

//  This constructor does nothing except default the record
PDBAtom::PDBAtom( void) : PDBData( PDBData:: ATOM)
{
  serialnumber(default_serial);
  name("");
  alternatelocation("");
  residuename("");
  chain("");
  residueseq(default_residueseq);
  insertioncode("");
  xcoor(default_coor);
  ycoor(default_coor);
  zcoor(default_coor);
  occupancy(default_occupancy);
  temperaturefactor(default_temperaturefactor);
  footnote(no_footnote);
  segmentname("");
}


PDBAtom::~PDBAtom( void) {
}


// Create an atom or hetatm record given that it is column based data
void PDBAtom::parse_column_data( const char *data)
{
 int len = strlen(data);  // to check that there is info
 char tempstr[100];   
 int tempint;
 Real tempReal;

     // set the serial number
 scan(data, len, SSERIAL, LSERIAL, &tempint, default_serial);
 serialnumber( tempint );
 
     // set the name
 scan(data, len, SNAME, LNAME, tempstr);
 name( tempstr);
 
     // set the alternate location
 scan(data, len, SALT, LALT, tempstr);
 alternatelocation( tempstr);

     // set the residue name
 scan(data, len, SRESNAME, LRESNAME, tempstr);
 residuename( tempstr);
 
     // set the chain 
 scan(data, len, SCHAIN, LCHAIN, tempstr);
 chain( tempstr);
 
     // set the residue sequence
 {
    // If there are more than 9999 residues, X-Plor uses A000, A001,
    // ..., A999, B000.  Since we try to be X-Plor compatible,
    // we doo the following
    char s[10];
    scan(data, len, SRESSEQ, LRESSEQ, s);
    if (s[0] < '0' || s[0] > '9') {
       static int elvis_count = 0;
       int num = (s[0] - 'A') * 1000 + 10000;
       num += atoi(s+1);
       if (s[0] < 'A' || s[0] > 'Z') {
	  if (elvis_count == 0) {
     iout << iWARN << "Man, tiny Elvis, that number is huge!\n"
     << iWARN << "We don't know how X-Plor represents over Z999 residues\n"
     << iWARN << "And you just tried " << s << " - so we'll fake it as " << num << "\n"
     << iWARN << "This is reversible, but only inside this program.\n" << endi;
               elvis_count = 1;
	   }
       } else {
	  elvis_count = 0;
       }
       residueseq( num);
    } else {
       residueseq(atoi(s));
    }
    // This used to be simply
//    scan(data, len, SRESSEQ, LRESSEQ, &tempint, default_residueseq);
//    residueseq( tempint);
 }
 
     // set the insertion code 
 scan(data, len, SINSERT, LINSERT, tempstr);
 insertioncode( tempstr);
 
     // set the X, Y, and Z coordinates
 scan(data, len, SX, LCOOR, &tempReal, default_coor);
 xcoor( tempReal);
 scan(data, len, SY, LCOOR, &tempReal, default_coor);
 ycoor( tempReal);
 scan(data, len, SZ, LCOOR, &tempReal, default_coor);
 zcoor( tempReal);

     // set the occupancy 
 scan(data, len, SOCC, LOCC, &tempReal, default_occupancy);
 occupancy( tempReal);
 
     // set the temperature factor
 scan(data, len, STEMPF, LTEMPF, &tempReal, default_temperaturefactor);
 temperaturefactor( tempReal);
 
     // set the footnote
 scan(data, len, SFOOT, LFOOT, &tempint, no_footnote);
 footnote( tempint);

   // this is for XPLOR style PDBs which have a segment name
 scan(data, len, SSEGNAME, LSEGNAME, tempstr);
 segmentname( tempstr);
}

void PDBAtom::parse_field_data( const char *data)
{
  char tempstr[100];
  // I already know that the first field starts with a '#' and that
  // the second is either ATOM or HETATM, so I'll start with the third
  field(data, 3, tempstr);
  serialnumber( tempstr[0] != '#' ? atoi(tempstr) : default_serial );

  field(data, 4, tempstr);
  name( tempstr[0] != '#' ? tempstr : "" );

  field(data, 5, tempstr);
  alternatelocation( tempstr[0] != '#' ? tempstr : "" );

  field(data, 6, tempstr);
  residuename( tempstr[0] != '#' ? tempstr : "" );
  		
  field(data, 7, tempstr);
  chain( tempstr[0] != '#' ? tempstr : "" );

  field(data, 8, tempstr);
  residueseq( tempstr[0] != '#' ? atoi(tempstr) : default_residueseq );

  field(data, 9, tempstr);
  insertioncode( tempstr[0] != '#' ? tempstr : "" );

  field(data, 10, tempstr);
  xcoor( tempstr[0] != '#' ?
  	atof( tempstr) : default_coor);  // WARNING: assumes Real <= double
  field(data, 11, tempstr);
  ycoor( tempstr[0] != '#' ?
  	atof( tempstr) : default_coor);  // WARNING: assumes Real <= double
  field(data, 12, tempstr);
  zcoor( tempstr[0] != '#' ?
  	atof( tempstr) : default_coor);  // WARNING: assumes Real <= double

  field(data, 13, tempstr);
  occupancy( tempstr[0] != '#' ? 
	atof( tempstr) : default_occupancy );// WARNING: assumes Real <= double

  field(data, 14, tempstr);
  temperaturefactor( tempstr[0] != '#' ?
  	atof( tempstr) : default_temperaturefactor ); // WARNING: ditto

  field(data, 15, tempstr);
  footnote( tempstr[0] != '#' ? atoi(tempstr) : no_footnote);
  
  field(data, 16, tempstr);
  segmentname( tempstr[0] != '#' ? tempstr : "");
}

  // get/ set the serial number
int PDBAtom:: serialnumber( void) 
{ return myserialnumber; }
void PDBAtom:: serialnumber( int newserialnumber)
{ myserialnumber = newserialnumber; }

  // get/ set the serial number
const char* PDBAtom:: name( void)
{ return myname; }
void PDBAtom:: name( const char *newname)
{ strncpy(myname, newname, LNAME); myname[LNAME]=0; }

  // get/ set the alternate location
const char* PDBAtom:: alternatelocation( void)
{ return myalternatelocation; }
void PDBAtom:: alternatelocation( const char *newalternatelocation)
{ strncpy(myalternatelocation, newalternatelocation, LALT);
  myalternatelocation[LALT]=0;}

  // get/ set the residue name
const char* PDBAtom:: residuename( void)
{ return myresiduename; }
void PDBAtom:: residuename( const char *newresiduename)
{ strncpy(myresiduename, newresiduename, LRESNAME); myresiduename[LRESNAME]=0;}

  // get/ set the chain indentifier
const char* PDBAtom:: chain( void)
{ return mychain; }
void PDBAtom:: chain( const char *newchain)
{ strncpy(mychain, newchain, LCHAIN); mychain[LCHAIN]=0;}

  // get/ set the residue sequence number
int PDBAtom:: residueseq( void)
{ return myresidueseq; }
void PDBAtom:: residueseq( int newresidueseq)
{ myresidueseq = newresidueseq; }

  // get/ set the insertion code
const char* PDBAtom:: insertioncode( void)
{ return myinsertioncode; }
void PDBAtom:: insertioncode( const char *newinsertioncode)
{ strncpy(myinsertioncode, newinsertioncode, LINSERT); 
  myinsertioncode[LINSERT]=0;}

  // get/ set the different coordinates
  // either 1 by 1 ...
Real PDBAtom:: xcoor( void)
{ return mycoor[0]; }
void PDBAtom:: xcoor( Real newxcoor)
{ mycoor[0] = newxcoor; }
Real PDBAtom:: ycoor( void)
{ return mycoor[1]; }
void PDBAtom:: ycoor( Real newycoor)
{ mycoor[1] = newycoor; }
Real PDBAtom:: zcoor( void)
{ return mycoor[2]; }
void PDBAtom:: zcoor( Real newzcoor)
{ mycoor[2] = newzcoor; }
   // ...or all three at once
const Real *PDBAtom:: coordinates( void)
{ return mycoor; }
void PDBAtom:: coordinates(const Real *newcoordinates)
{ for (int i=0; i<3; i++) mycoor[i] = newcoordinates[i]; }

  // get/ set the occupancy
Real PDBAtom:: occupancy( void)
{ return myoccupancy ;}
void PDBAtom:: occupancy( Real newoccupancy)
{ myoccupancy = newoccupancy; }

  // get/ set the temperature factor
Real PDBAtom:: temperaturefactor( void)
{ return mytemperaturefactor; }
void PDBAtom:: temperaturefactor( Real newtemperaturefactor)
{ mytemperaturefactor = newtemperaturefactor; }

  // get/ set the footnote
int PDBAtom:: footnote( void)
{ return myfootnote; }
void PDBAtom:: footnote( int newfootnote)
{ myfootnote = newfootnote; }

  // get/ set the segment name
  // this is not part of the PDB format but is used by XPLOR instead of
  // the chain identifier (see XPLOR 3.1 manual, p 104)
const char* PDBAtom:: segmentname( void)
{ return mysegmentname; }
void PDBAtom:: segmentname( const char *newsegmentname)
{ strncpy(mysegmentname, newsegmentname, LSEGNAME); mysegmentname[LSEGNAME]=0;}

 
// the function to print out an ATOM or HETATM
// size or outstr must be >= 80!
void PDBAtom::sprint_columns( char *outstr)
{ 
  int i;
  for (i=0; i<79; i++)   // dump spaces in outstr -- must be length > 80!!
    outstr[i] = 32;
  outstr[i] = 0;             // and terminate

  sprintcol(outstr, STYPE, LTYPE, PDBNames[type()] );
  sprintcol(outstr, SSERIAL, LSERIAL, serialnumber());
  {
     // For X-Plor compatability, if the name is 1, 2, or 3
     // characters, start it in the 2nd column of the field.
     // 4 letter names use the first column
     if (strlen(name()) == 4) {
	sprintcol(outstr, SNAME, LNAME, name());
     } else {
	sprintcol(outstr, SNAME+1, LNAME-1, name());
     }
  }
  sprintcol(outstr, SALT, LALT, alternatelocation());
  sprintcol(outstr, SRESNAME, LRESNAME, residuename());
  sprintcol(outstr, SCHAIN, LCHAIN, chain());
  {
     // Again, I may have to convert from a number > 9999 to
     // A000 or whatever (see the comments for the residueseq input)
     if (residueseq() <= 9999) {
	sprintcol(outstr, SRESSEQ, LRESSEQ, residueseq());
     } else {
	int val = residueseq() / 1000 - 10;  // integer arithmetic
	int modulo = residueseq() % 1000;
	char s[10];
	sprintf(s, "%c%03d", 'A' + val, modulo);
	sprintcol(outstr, SRESSEQ, LRESSEQ, s);
     }
     // This used to be just ...
//  sprintcol(outstr, SRESSEQ, LRESSEQ, residueseq());
  }
  sprintcol(outstr, SINSERT, LINSERT, insertioncode());
  sprintcol(outstr, SX, LCOOR, LCOORPREC, xcoor());
  sprintcol(outstr, SY, LCOOR, LCOORPREC, ycoor());
  sprintcol(outstr, SZ, LCOOR, LCOORPREC, zcoor());
  sprintcol(outstr, SOCC, LOCC, LOCCPREC, occupancy());
  sprintcol(outstr, STEMPF, LTEMPF, LTEMPFPREC, temperaturefactor());
  if (footnote() == no_footnote)     // special case when no footnote
    sprintcol(outstr, SFOOT, LFOOT, "");
   else
    sprintcol(outstr, SFOOT, LFOOT, footnote() );         
  sprintcol(outstr, SSEGNAME, LSEGNAME, segmentname());
}

/// Print the output by fields; if the values are not known/ are defaults - put a #
void PDBAtom::sprint_fields( char *outstr)
{
  char tmpstr[50];
  sprintf(outstr, "# %s", PDBNames[type()]);
  if (serialnumber() == default_serial)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %i", serialnumber());
  strcat(outstr, tmpstr);
  if (name()[0] == 0)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %s", name());
  strcat(outstr, tmpstr);
  if (alternatelocation()[0] == 0)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %s", alternatelocation());
  strcat(outstr, tmpstr);
  if (residuename()[0] == 0)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %s", residuename());
  strcat(outstr, tmpstr);
  if (chain()[0] == 0)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %s", chain());
  strcat(outstr, tmpstr);
  if (residueseq() == default_residueseq)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %d", residueseq());
  strcat(outstr, tmpstr);
  if (insertioncode()[0] == 0)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %s", insertioncode());
  strcat(outstr, tmpstr);
  if (xcoor() == default_coor)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %*.*f", LCOOR, LCOORPREC, xcoor());
  strcat(outstr, tmpstr);
  if (ycoor() == default_coor)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %*.*f", LCOOR, LCOORPREC, ycoor());
  strcat(outstr, tmpstr);
  if (zcoor() == default_coor)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %*.*f", LCOOR, LCOORPREC, zcoor());
  strcat(outstr, tmpstr);
//  if (occupancy() == default_occupancy)  // no way to tell if the occ. is the default
//      sprintf(tmpstr, " #");
//     else
      sprintf(tmpstr, " %*.*f",  LOCC, LOCCPREC, occupancy());
  strcat(outstr, tmpstr);
//  if (temperaturefactor() == default_temperaturefactor) // ditto previous
//      sprintf(tmpstr, " #");
//     else
      sprintf(tmpstr, " %*.*f", LTEMPF, LTEMPFPREC, temperaturefactor());
  strcat(outstr, tmpstr);
  if (footnote() == no_footnote)     // special case when no footnote
    sprintf(tmpstr, " #");
   else
    sprintf(tmpstr, " %d",  footnote() );         
  strcat(outstr, tmpstr);
  if (segmentname()[0] == 0)
      sprintf(tmpstr, " #");
     else
      sprintf(tmpstr, " %s", segmentname());
  strcat(outstr, tmpstr);

}

void PDBAtom::sprint( char *outstr, PDBFormatStyle usestyle)
{
 if (usestyle == PDBData::COLUMNS)
   sprint_columns( outstr);
 else
   sprint_fields( outstr);
}

//****************** The wrapper for all of the functions ************///
PDBData *new_PDBData(const char *data)  // nasty
{
  char temps1[100];
  char temps2[100];
  char *temps;
  sscanf(data, "%s %s ", temps1, temps2);
  if (temps1[0] == '#')
     temps = temps2;
    else
     temps = temps1;

     // go through the list of possible PDB data types
     //this _should_ be the same as: for(PDBTypes i=HEADER; i<UNKNOWN; i++)
  for (int i=0; i< (int)(sizeof(PDBData::PDBNames) /
                  sizeof(PDBData::PDBNames[0])); i++)
    if (!strcmp(temps, PDBData::PDBNames[i]))
       switch(i) {
         case PDBData::ATOM: return new PDBAtomRecord(data);
         case PDBData::HETATM: return new PDBHetatm(data);
         default: return new PDBUnknown(data);
       }
  // Now, if HETATM is right next to an aton number (like HETATM12345) then the above
  // test will fail, so I have to special case it:
  if (!strncmp(temps, PDBData::PDBNames[PDBData::HETATM], sizeof(PDBData::PDBNames[PDBData::HETATM]))) {
    return new PDBHetatm(data);
  }
   //  Hmm, looks like it isn't any data type, so I'll fake it
   return new PDBUnknown(data);
}

/////********** some test routines left in to show you examples 
//#define TEST_PDBREADER
#ifdef TEST_PDBREADER
main()
{
  char tempstr[100];
  PDBAtomRecord atom("ATOM   6312  CB TALA 3 235I     24.681  54.463 137.827  1.00 51.30      VP3");
  atom.sprint(tempstr, PDBData::COLUMNS);
  cout << tempstr << '\n';
  atom.sprint(tempstr, PDBData::FIELDS);
  cout << tempstr << '\n';
  cout << "Serial number : "  << atom.serialnumber()         <<  "\n";
  cout << "name          : '" << atom.name()                 <<  "'\n";
  cout << "alt. location : '" << atom.alternatelocation()    <<  "'\n";
  cout << "residue name  : '" << atom.residuename()          <<  "'\n";
  cout << "chain         : '" << atom.chain()                <<  "'\n";
  cout << "residue seq   : "  << atom.residueseq()           <<  "\n";
  cout << "insertion code: '" << atom.insertioncode()        << "'\n";
  cout << "X coordinate  : "  << atom.xcoor()                <<  "\n";
  cout << "Y coordinate  : "  << atom.ycoor()                << "\n";
  cout << "Z coordinate  : "  << atom.zcoor()                << "\n";
  cout << "occupancy     : "  << atom.occupancy()            << "\n";
  cout << "temperature factor: " << atom.temperaturefactor() << "\n";
  cout << "footnote      : " << atom.footnote()              << "\n";
  cout << "segment name  : '" << atom.segmentname()          << "'\n";
  cout << '\n';
  
  PDBAtomRecord atom2("# ATOM 6312 CB T ALA 3 235 I 24.681 54.463 137.827 1.00 51.30 # VP3");
  atom2.sprint(tempstr, PDBData::COLUMNS);
  cout << tempstr << '\n';
  atom2.sprint(tempstr, PDBData::FIELDS);
  cout << tempstr << '\n';
  cout << "Serial number : "  << atom2.serialnumber()         <<  "\n";
  cout << "name          : '" << atom2.name()                 <<  "'\n";
  cout << "alt. location : '" << atom2.alternatelocation()    <<  "'\n";
  cout << "residue name  : '" << atom2.residuename()          <<  "'\n";
  cout << "chain         : '" << atom2.chain()                <<  "'\n";
  cout << "residue seq   : "  << atom2.residueseq()           <<  "\n";
  cout << "insertion code: '" << atom2.insertioncode()        << "'\n";
  cout << "X coordinate  : "  << atom2.xcoor()                <<  "\n";
  cout << "Y coordinate  : "  << atom2.ycoor()                << "\n";
  cout << "Z coordinate  : "  << atom2.zcoor()                << "\n";
  cout << "occupancy     : "  << atom2.occupancy()            << "\n";
  cout << "temperature factor: " << atom2.temperaturefactor() << "\n";
  cout << "footnote      : " << atom2.footnote()              << "\n";
  cout << "segment name  : '" << atom2.segmentname()          << "'\n";
  cout << '\n';
  

  PDBAtomRecord atom3("# ATOM # # Q WER # # # # 123.456 # # # 9 LAST anything?");
  atom3.sprint(tempstr, PDBData::COLUMNS);
  cout << tempstr << '\n';
  atom3.sprint(tempstr, PDBData::FIELDS);
  cout << tempstr << '\n';
  cout << "Serial number : "  << atom3.serialnumber()         <<  "\n";
  cout << "name          : '" << atom3.name()                 <<  "'\n";
  cout << "alt. location : '" << atom3.alternatelocation()    <<  "'\n";
  cout << "residue name  : '" << atom3.residuename()          <<  "'\n";
  cout << "chain         : '" << atom3.chain()                <<  "'\n";
  cout << "residue seq   : "  << atom3.residueseq()           <<  "\n";
  cout << "insertion code: '" << atom3.insertioncode()        << "'\n";
  cout << "X coordinate  : "  << atom3.xcoor()                <<  "\n";
  cout << "Y coordinate  : "  << atom3.ycoor()                << "\n";
  cout << "Z coordinate  : "  << atom3.zcoor()                << "\n";
  cout << "occupancy     : "  << atom3.occupancy()            << "\n";
  cout << "temperature factor: " << atom3.temperaturefactor() << "\n";
  cout << "footnote      : " << atom3.footnote()              << "\n";
  cout << "segment name  : '" << atom3.segmentname()          << "'\n";
  cout << '\n';
  
}
#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: PDBData.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1006 $	$Date: 1999/07/09 21:05:49 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: PDBData.C,v $
 * Revision 1.1006  1999/07/09 21:05:49  jim
 * Eliminated compiler warnings.
 *
 * Revision 1.1005  1999/03/10 23:03:14  jim
 * Eliminated warnings sent to stderr.
 *
 * Revision 1.1004  1998/10/24 19:57:47  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1003  1998/03/26 23:28:31  jim
 * Small changes for KCC port.  Altered use of strstream in ComputeFreeEnergy.
 *
 * Revision 1.1002  1997/08/12 22:18:43  milind
 * Made NAMD2 to link on Solaris machines.
 *
 * Revision 1.1001  1997/03/19 11:54:40  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:58:55  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:03  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:36:38  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.7  1995/07/07 19:54:15  nelson
 * Added int to elvis_count declaration
 *
 * Revision 1.6  95/05/19  11:45:27  11:45:27  dalke (Andrew Dalke)
 * made atom name output compatible with X-Plor
 * ie. 1, 2, and 3 letter names start in the 2nd
 * column of the 4 column field.  4 letter namds
 * start in the first
 * 
 * Revision 1.5  95/05/19  10:47:37  10:47:37  dalke (Andrew Dalke)
 * Made residue sequence output compatible with X-Plor when >9999
 * (ie, 10000 == A000, etc.)
 * 
 * Revision 1.4  95/03/08  14:27:42  14:27:42  nelson (Mark T. Nelson)
 * Added copyright message
 * 
 * Revision 1.3  94/10/23  21:57:59  21:57:59  dalke (Andrew Dalke)
 * Fixed non-bug with
 * 1) "    " a blank field in integer and Real input got converted
 *   to 0 or 0.0 while a non-string (ie., past input) got converted
 *   to the default value.  I made the blank string get converted to
 *   the default as well
 * 2) During field output, if the temperature factor was 0, a # sign
 *   was printed.  Since the T.F. field can be used for things other
 *   than temperature factor, I made it print the value, even if it
 *   is the default value.  This is same manner by which I treat
 *   the occupancy field.
 * 
 * Revision 1.2  94/10/19  02:15:03  02:15:03  dalke (Andrew Dalke)
 * Fixed bug with records of the type "HETATM12345".
 * 
 * Revision 1.1  94/07/05  12:52:16  12:52:16  dalke (Andrew Dalke)
 * Initial revision
 * 
 ***************************************************************************/
