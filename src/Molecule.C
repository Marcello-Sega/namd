/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// Molecule.C is compiled twice!
// MOLECULE2_C undefined only compiles first half of file
// MOLECULE2_C defined only compiles second half of file
// This is shameful but it works.  Molecule needs refactoring badly.

/*
   The class Molecule is used to hold all of the structural information
   for a simulation.  This information is read in from a .psf file and
   cross checked with the Parameters object passed in.  All of the structural
   information is then stored in arrays for use.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "InfoStream.h"
#include "Molecule.h"
#include "strlib.h"
#include "MStream.h"
#include "Communicate.h"
// #include "Node.h"
#include "ObjectArena.h"
#include "Parameters.h"
#include "PDB.h"
#include "SimParameters.h"
#include "Hydrogen.h"
#include "UniqueSetIter.h"
#include "charm++.h"
/* BEGIN gf */
#include "ComputeGridForce.h"
#include "GridForceGrid.h"

#include "MGridforceParams.h"
/* END gf */

#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#include "CompressPsf.h"
#include <deque>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

class ResidueLookupElem
{
public:
  char mySegid[11];
  ResidueLookupElem *next;	// stored as a linked list
  int firstResid;		// valid resid is >= firstResid
  int lastResid;		// valid resid is <= lastResid
  ResizeArray<int> atomIndex;	// 0-based index for first atom in residue

  ResidueLookupElem(void) { next = 0; firstResid = -1; lastResid = -1; }
  ~ResidueLookupElem(void) { delete next; }
  int lookup(const char *segid, int resid, int *begin, int *end) const;
  ResidueLookupElem* append(const char *segid, int resid, int aid);
};

#ifndef MOLECULE2_C  // first object file

#ifdef MEM_OPT_VERSION
template int lookupCstPool<AtomSignature>(const vector<AtomSignature>&, const AtomSignature&);
template int lookupCstPool<ExclusionSignature>(const vector<ExclusionSignature>&, const ExclusionSignature&);
#endif

int ResidueLookupElem::lookup(
	const char *segid, int resid, int *begin, int *end) const {
    const ResidueLookupElem *elem = this;
    int rval = -1;  // error
    while ( elem && strcasecmp(elem->mySegid,segid) ) elem = elem->next;
    if ( elem && (resid >= elem->firstResid) && (resid <= elem->lastResid) ) {
      *begin = elem->atomIndex[resid - elem->firstResid];
      *end = elem->atomIndex[resid - elem->firstResid + 1];
      rval = 0;  // no error
    }
    return rval;
}

ResidueLookupElem* ResidueLookupElem::append(
	const char *segid, int resid, int aid) {
    ResidueLookupElem *rval = this;
    if ( firstResid == -1 ) {  // nothing added yet
      strcpy(mySegid,segid);
      firstResid = resid;
      lastResid = resid;
      atomIndex.add(aid);
      atomIndex.add(aid+1);
    } else if ( ! strcasecmp(mySegid,segid) ) {  // same segid
      if ( resid == lastResid ) {  // same resid
        atomIndex[lastResid - firstResid + 1] = aid + 1;
      } else if ( resid < lastResid ) {  // error
	// We can work around this by creating a new segment.
	iout << iWARN << "Residue " << resid <<
	  " out of order in segment " << segid <<
	  ", lookup for additional residues in this segment disabled.\n" << endi;
	rval = next = new ResidueLookupElem;
	next->append(segid,resid,aid);
      } else {  // new resid
        for ( ; lastResid < resid; ++lastResid ) atomIndex.add(aid);
        atomIndex[lastResid - firstResid + 1] = aid + 1;
      }
    } else {  // new segid
      rval = next = new ResidueLookupElem;
      next->append(segid,resid,aid);
    }
    return rval;
}


//  Lookup atom id from segment, residue, and name
int Molecule::get_atom_from_name(
	const char *segid, int resid, const char *aname) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }

  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return -1;
  for ( ; i < end; ++i ) {
    #ifdef MEM_OPT_VERSION    
    Index idx = atomNames[i].atomnameIdx;
    if(!strcasecmp(aname, atomNamePool[idx])) return i;
    #else
    if ( ! strcasecmp(aname,atomNames[i].atomname) ) return i;
    #endif
  }
  return -1;
}

//  Lookup number of atoms in residue from segment and residue
int Molecule::get_residue_size(
	const char *segid, int resid) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }
  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return 0;
  return ( end - i );
}

//  Lookup atom id from segment, residue, and index in residue
int Molecule::get_atom_from_index_in_residue(
	const char *segid, int resid, int index) const {

  if (atomNames == NULL || resLookup == NULL)
  {
    NAMD_die("Tried to find atom from name on node other than node 0");
  }
  int i = 0;
  int end = 0;
  if ( resLookup->lookup(segid,resid,&i,&end) ) return -1;
  if ( index >= 0 && index < ( end - i ) ) return ( index + i );
  return -1;
}

/************************************************************************/
/*                  */
/*      FUNCTION initialize  */
/*                  */
/*  This is the initializer for the Molecule class.  It simply sets */
/*  the counts for all the various parameters to 0 and sets the pointers*/
/*  to the arrays that will store these parameters to NULL, since they  */
/*  have not been allocated yet.          */
/*                  */
/************************************************************************/

void Molecule::initialize(SimParameters *simParams, Parameters *param)
{
  if ( sizeof(int32) != 4 ) { NAMD_bug("sizeof(int32) != 4"); }
  this->simParams = simParams;
  this->params = param;

  /*  Initialize array pointers to NULL  */
  atoms=NULL;
  atomNames=NULL;
  resLookup=NULL;

  // DRUDE
  is_drude_psf = 0;  // assume not Drude model
  drudeConsts=NULL;
  lphosts=NULL;
  anisos=NULL;
  tholes=NULL;
  lphostIndexes=NULL;
  // DRUDE

  //for compressing molecule info
  atomSegResids=NULL;

  if ( simParams->globalForcesOn ) {
    resLookup = new ResidueLookupElem;
  }

  #ifdef MEM_OPT_VERSION
  eachAtomSig = NULL;
  atomSigPoolSize = 0;
  atomSigPool = NULL;
  massPoolSize = 0;
  atomMassPool = NULL;
  eachAtomMass = NULL;
  chargePoolSize = 0;
  atomChargePool = NULL;
  eachAtomCharge = NULL;
  #else
  bonds=NULL;
  angles=NULL;
  dihedrals=NULL;
  impropers=NULL;
  crossterms=NULL;
  #endif

  donors=NULL;
  acceptors=NULL;
  

  #ifndef MEM_OPT_VERSION      
  tmpArena=NULL;
  exclusions=NULL;
  bondsWithAtom=NULL;
  bondsByAtom=NULL;
  anglesByAtom=NULL;
  dihedralsByAtom=NULL;
  impropersByAtom=NULL;
  crosstermsByAtom=NULL;
  // DRUDE
  tholesByAtom=NULL;
  anisosByAtom=NULL;
  // DRUDE
  #endif

  #ifdef MEM_OPT_VERSION
  exclSigPool = NULL;
  exclChkSigPool = NULL;
  exclSigPoolSize = 0;
  eachAtomExclSig = NULL;
  #else
  exclusionsByAtom=NULL;
  fullExclusionsByAtom=NULL;
  modExclusionsByAtom=NULL;
  all_exclusions=NULL;
  #endif

  langevinParams=NULL;
  fixedAtomFlags=NULL;

  #ifdef MEM_OPT_VERSION  
  clusterSigs=NULL;
  #else
  cluster=NULL;  
  #endif
  clusterSize=NULL;

  exPressureAtomFlags=NULL;
  rigidBondLengths=NULL;
  consIndexes=NULL;
  consParams=NULL;
  /* BEGIN gf */
  gridfrcIndexes=NULL;
  gridfrcParams=NULL;
  gridfrcGrid=NULL;
  numGridforces=NULL;
  /* END gf */
  stirIndexes=NULL;
  stirParams=NULL;
  movDragIndexes=NULL;
  movDragParams=NULL;
  rotDragIndexes=NULL;
  rotDragParams=NULL;
  consTorqueIndexes=NULL;
  consTorqueParams=NULL;
  consForceIndexes=NULL;
  consForce=NULL;
//fepb
  fepAtomFlags=NULL;
//fepe

  nameArena = new ObjectArena<char>;
  // nameArena->setAlignment(8);
  // arena->setAlignment(32);
  #ifndef MEM_OPT_VERSION
  arena = new ObjectArena<int32>;
  exclArena = new ObjectArena<char>;
  #endif
  // exclArena->setAlignment(32);

  /*  Initialize counts to 0 */
  numAtoms=0;
  numRealBonds=0;
  numBonds=0;
  numAngles=0;
  numDihedrals=0;
  numImpropers=0;
  numTholes=0;
  numAnisos=0;
  numCrossterms=0;
  numDonors=0;
  numAcceptors=0;
  numExclusions=0;

  // DRUDE
  numLonepairs=0;
  numDrudeAtoms=0;
  numLphosts=0;
  numAnisos=0;
  // DRUDE

  numConstraints=0;
  numStirredAtoms=0;
  numMovDrag=0;
  numRotDrag=0;
  numConsTorque=0;
  numConsForce=0;
  numFixedAtoms=0;
  numFixedGroups=0;
  numExPressureAtoms=0;
  numRigidBonds=0;
  numFixedRigidBonds=0;
  numMultipleDihedrals=0;
  numMultipleImpropers=0;
  numCalcBonds=0;
  numCalcAngles=0;
  numCalcDihedrals=0;
  numCalcImpropers=0;
  numCalcTholes=0;
  numCalcAnisos=0;
  numCalcCrossterms=0;
  numCalcExclusions=0;

//fepb
  numFepInitial = 0;
  numFepFinal = 0;
//fepe

  //fields related with pluginIO-based loading molecule structure
  occupancy = NULL;
  bfactor = NULL;
}

/*      END OF FUNCTION initialize */

/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for the Molecule class. */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param)
{
  initialize(simParams,param);
}

/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for the Molecule class from CHARMM/XPLOR files. */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param, char *filename, ConfigList *cfgList)
{
  initialize(simParams,param);

  if(simParams->useCompressedPsf)
      read_compressed_psf_file(filename, param, cfgList);
  /*else if(simParams->genCompressedPsf){      
      compress_psf_file(this, filename, param, simParams, cfgList);
  }*/      
  else
      read_psf_file(filename, param);
}

/************************************************************************/
/*                                                                      */
/*      FUNCTION Molecule                                               */
/*                                                                      */
/*  This is the constructor for the Molecule class from plugin IO.      */
/*                                                                      */
/************************************************************************/
Molecule::Molecule(SimParameters *simParams, Parameters *param, molfile_plugin_t *pIOHdl, void *pIOFileHdl, int natoms)
{
#ifdef MEM_OPT_VERSION
  NAMD_die("Sorry, plugin IO is not supported in the memory optimized version.");
#endif
    initialize(simParams, param);
    numAtoms = natoms;
    int optflags = MOLFILE_BADOPTIONS;
    molfile_atom_t *atomarray = (molfile_atom_t *) malloc(natoms*sizeof(molfile_atom_t));
    memset(atomarray, 0, natoms*sizeof(molfile_atom_t));

    //1a. read basic atoms information
    int rc = pIOHdl->read_structure(pIOFileHdl, &optflags, atomarray);
    if (rc != MOLFILE_SUCCESS && rc != MOLFILE_NOSTRUCTUREDATA) {
        free(atomarray);
        NAMD_die("ERROR: plugin failed reading structure data");
    }
    if(optflags == MOLFILE_BADOPTIONS) {
        free(atomarray);
        NAMD_die("ERROR: plugin didn't initialize optional data flags");
    }
    if(optflags & MOLFILE_OCCUPANCY) {
        setOccupancyData(atomarray);
    }
    if(optflags & MOLFILE_BFACTOR) {
        setBFactorData(atomarray);
    }
    //1b. load basic atoms information to the molecule object
    plgLoadAtomBasics(atomarray);    
    free(atomarray);

    //2a. read bonds
    //indices are one-based in read_bonds
    int *from, *to;
    float *bondorder;
    if(pIOHdl->read_bonds!=NULL) {
        if(pIOHdl->read_bonds(pIOFileHdl, &numBonds, &from, &to, &bondorder)){
            NAMD_die("ERROR: failed reading bond information.");
        }
    }    
    //2b. load bonds information to the molecule object
    if(numBonds!=0) {
        plgLoadBonds(from,to);
    }

    //3a. read other bonded structures
    int *plgAngles, *plgDihedrals, *plgImpropers, *plgCterms;
    int ctermcols, ctermrows;
    double *angleforces,  *dihedralforces, *improperforces, *ctermforces;

    plgAngles=plgDihedrals=plgImpropers=plgCterms=NULL;
    if(pIOHdl->read_angles!=NULL) {
        if(pIOHdl->read_angles(pIOFileHdl,
                               &numAngles, &plgAngles, &angleforces,
                               &numDihedrals, &plgDihedrals, &dihedralforces,
                               &numImpropers, &plgImpropers, &improperforces,
                               &numCrossterms, &plgCterms, &ctermcols, &ctermrows,
                               &ctermforces)) {
            NAMD_die("ERROR: failed reading angle information.");
        }
    }
    //3b. load other bonded structures to the molecule object
    if(numAngles!=0) plgLoadAngles(plgAngles);
    if(numDihedrals!=0) plgLoadDihedrals(plgDihedrals);
    if(numImpropers!=0) plgLoadImpropers(plgImpropers);
    if(numCrossterms!=0) plgLoadCrossterms(plgCterms);

  numRealBonds = numBonds;
  build_atom_status();

}

/*      END OF FUNCTION Molecule      */

/************************************************************************/
/*                  */
/*        FUNCTION Molecule      */
/*                  */
/*  This is the destructor for the class Molecule.  It simply frees */
/*  the memory allocated for each of the arrays used to store the       */
/*  structure information.            */
/*                  */
/************************************************************************/

Molecule::~Molecule()
{
  /*  Check to see if each array was ever allocated.  If it was   */
  /*  then free it            */
  if (atoms != NULL)
    delete [] atoms;

  if (atomNames != NULL)
  {
    // subarrarys allocated from arena - automatically deleted
    delete [] atomNames;
  }
  delete nameArena;

  if (resLookup != NULL)
    delete resLookup;

  // DRUDE: free arrays read from PSF
  if (drudeConsts != NULL) delete [] drudeConsts;
  if (lphosts != NULL) delete [] lphosts;
  if (anisos != NULL) delete [] anisos;
  if (tholes != NULL) delete [] tholes;
  if (lphostIndexes != NULL) delete [] lphostIndexes;
  // DRUDE

  #ifdef MEM_OPT_VERSION
  if(eachAtomSig) delete [] eachAtomSig;
  if(atomSigPool) delete [] atomSigPool;
  #else
  if (bonds != NULL)
    delete [] bonds;

  if (angles != NULL)
    delete [] angles;

  if (dihedrals != NULL)
    delete [] dihedrals;

  if (impropers != NULL)
    delete [] impropers;

  if (crossterms != NULL)
    delete [] crossterms;

  if (exclusions != NULL)
    delete [] exclusions;
  #endif

  if (donors != NULL)
    delete [] donors;

  if (acceptors != NULL)
    delete [] acceptors;  

  #ifdef MEM_OPT_VERSION
  if(exclSigPool) delete [] exclSigPool;
  if(exclChkSigPool) delete [] exclChkSigPool;
  if(eachAtomExclSig) delete [] eachAtomExclSig;
  #else
  if (bondsByAtom != NULL)
       delete [] bondsByAtom;
  
  if (anglesByAtom != NULL)
       delete [] anglesByAtom;
  
  if (dihedralsByAtom != NULL)
       delete [] dihedralsByAtom;
  
  if (impropersByAtom != NULL)
       delete [] impropersByAtom;
  
  if (crosstermsByAtom != NULL)
       delete [] crosstermsByAtom;  

  if (exclusionsByAtom != NULL)
       delete [] exclusionsByAtom;
  
  if (fullExclusionsByAtom != NULL)
       delete [] fullExclusionsByAtom;
  
  if (modExclusionsByAtom != NULL)
       delete [] modExclusionsByAtom;
  
  if (all_exclusions != NULL)
       delete [] all_exclusions;

  // DRUDE
  if (tholesByAtom != NULL)
       delete [] tholesByAtom;
  if (anisosByAtom != NULL)
       delete [] anisosByAtom;
  // DRUDE
  #endif


  if (fixedAtomFlags != NULL)
       delete [] fixedAtomFlags;

  if (stirIndexes != NULL)
    delete [] stirIndexes;


  #ifdef MEM_OPT_VERSION
  if(clusterSigs != NULL){      
      delete [] clusterSigs;
  }  
  #else
  if (cluster != NULL)
       delete [] cluster;  
  #endif
  if (clusterSize != NULL)
       delete [] clusterSize;

  if (exPressureAtomFlags != NULL)
       delete [] exPressureAtomFlags;

  if (rigidBondLengths != NULL)
       delete [] rigidBondLengths;

//fepb
  if (fepAtomFlags != NULL)
       delete [] fepAtomFlags;
//fepe


  #ifndef MEM_OPT_VERSION
  delete arena;
  delete exclArena;
  #endif
}
/*      END OF FUNCTION Molecule      */

/************************************************************************/
/*                  */
/*        FUNCTION read_psf_file      */
/*                  */
/*   INPUTS:                */
/*  fname - Name of the .psf file to read        */
/*  params - pointer to Parameters object to use to obtain          */
/*     parameters for vdWs, bonds, etc.      */
/*                  */
/*  This function reads a .psf file in.  This is where just about   */
/*   all of the structural information for this class comes from.  The  */
/*   .psf file contains descriptions of the atom, bonds, angles,        */
/*   dihedrals, impropers, and exclusions.  The parameter object is     */
/*   used to look up parameters for each of these entities.    */
/*                  */
/************************************************************************/

void Molecule::read_psf_file(char *fname, Parameters *params)

{
#ifdef MEM_OPT_VERSION
    return;
#else
  char err_msg[512];  //  Error message for NAMD_die
  char buffer[512];  //  Buffer for file reading
  int i;      //  Loop counter
  int NumTitle;    //  Number of Title lines in .psf file
  FILE *psf_file;    //  pointer to .psf file
  int ret_code;    //  ret_code from NAMD_read_line calls

  /* Try and open the .psf file           */
  if ( (psf_file = Fopen(fname, "r")) == NULL)
  {
    sprintf(err_msg, "UNABLE TO OPEN .psf FILE %s", fname);
    NAMD_die(err_msg);
  }

  /*  Read till we have the first non-blank line of file    */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to see if we dropped out of the loop because of a     */
  /*  read error.  This shouldn't happen unless the file is empty */
  if (ret_code!=0)
  {
    sprintf(err_msg, "EMPTY .psf FILE %s", fname);
    NAMD_die(err_msg);
  }

  /*  The first non-blank line should contain the word "psf".    */
  /*  If we can't find it, die.               */
  if (!NAMD_find_word(buffer, "psf"))
  {
    sprintf(err_msg, "UNABLE TO FIND \"PSF\" STRING IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  // DRUDE: set flag if we discover Drude PSF
  if (NAMD_find_word(buffer, "drude"))
  {
    if ( ! simParams->drudeOn ) {
      iout << iINFO << "Warning: Reading PSF supporting DRUDE without "
        "enabling the Drude model in the simulation config file\n";
    }
    is_drude_psf = 1;
  }
  // DRUDE

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to see if we dropped out of the loop because of a     */
  /*  read error.  This shouldn't happen unless there is nothing  */
  /*  but the PSF line in the file        */
  if (ret_code!=0)
  {
    sprintf(err_msg, "MISSING EVERYTHING BUT PSF FROM %s", fname);
    NAMD_die(err_msg);
  }

  /*  This line should have the word "NTITLE" in it specifying    */
  /*  how many title lines there are        */
  if (!NAMD_find_word(buffer, "NTITLE"))
  {
    sprintf(err_msg,"CAN NOT FIND \"NTITLE\" STRING IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  sscanf(buffer, "%d", &NumTitle);

  /*  Now skip the next NTITLE non-blank lines and then read in the*/
  /*  line which should contain NATOM        */
  i=0;

  while ( ((ret_code=NAMD_read_line(psf_file, buffer)) == 0) && 
    (i<NumTitle) )
  {
    if (!NAMD_blank_string(buffer))
      i++;
  }

  /*  Make sure we didn't exit because of a read error    */
  if (ret_code!=0)
  {
    sprintf(err_msg, "FOUND EOF INSTEAD OF NATOM IN PSF FILE %s", 
       fname);
    NAMD_die(err_msg);
  }

  while (NAMD_blank_string(buffer))
  {
    NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we have the line we want      */
  if (!NAMD_find_word(buffer, "NATOM"))
  {
    sprintf(err_msg, "DIDN'T FIND \"NATOM\" IN PSF FILE %s",
       fname);
    NAMD_die(err_msg);
  }

  /*  Read in the number of atoms, and then the atoms themselves  */
  sscanf(buffer, "%d", &numAtoms);

  read_atoms(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NBONDS IN PSF");
  }

  /*  Look for the string "NBOND"          */
  if (!NAMD_find_word(buffer, "NBOND"))
  {
    NAMD_die("DID NOT FIND NBOND AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of bonds and then the bonds themselves  */
  sscanf(buffer, "%d", &numBonds);

  if (numBonds)
    read_bonds(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NTHETA IN PSF");
  }

  /*  Look for the string "NTHETA"        */
  if (!NAMD_find_word(buffer, "NTHETA"))
  {
    NAMD_die("DID NOT FIND NTHETA AFTER BOND LIST IN PSF");
  }

  /*  Read in the number of angles and then the angles themselves */
  sscanf(buffer, "%d", &numAngles);

  if (numAngles)
    read_angles(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NPHI IN PSF");
  }

  /*  Look for the string "NPHI"          */
  if (!NAMD_find_word(buffer, "NPHI"))
  {
    NAMD_die("DID NOT FIND NPHI AFTER ANGLE LIST IN PSF");
  }

  /*  Read in the number of dihedrals and then the dihedrals      */
  sscanf(buffer, "%d", &numDihedrals);

  if (numDihedrals)
    read_dihedrals(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NIMPHI IN PSF");
  }

  /*  Look for the string "NIMPHI"        */
  if (!NAMD_find_word(buffer, "NIMPHI"))
  {
    NAMD_die("DID NOT FIND NIMPHI AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of Impropers and then the impropers  */
  sscanf(buffer, "%d", &numImpropers);

  if (numImpropers)
    read_impropers(psf_file, params);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NDON IN PSF");
  }

  /*  Look for the string "NDON"        */
  if (!NAMD_find_word(buffer, "NDON"))
  {
    NAMD_die("DID NOT FIND NDON AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of hydrogen bond donors and then the donors */
  sscanf(buffer, "%d", &numDonors);

  if (numDonors)
    read_donors(psf_file);

  /*  Read until we find the next non-blank line      */
  ret_code = NAMD_read_line(psf_file, buffer);

  while ( (ret_code==0) && (NAMD_blank_string(buffer)) )
  {
    ret_code = NAMD_read_line(psf_file, buffer);
  }

  /*  Check to make sure we didn't hit the EOF      */
  if (ret_code != 0)
  {
    NAMD_die("EOF ENCOUNTERED LOOKING FOR NACC IN PSF");
  }

  /*  Look for the string "NACC"        */
  if (!NAMD_find_word(buffer, "NACC"))
  {
    NAMD_die("DID NOT FIND NACC AFTER ATOM LIST IN PSF");
  }

  /*  Read in the number of hydrogen bond donors and then the donors */
  sscanf(buffer, "%d", &numAcceptors);

  if (numAcceptors)
    read_acceptors(psf_file);

  /*  look for the explicit non-bonded exclusion section.     */
  while (!NAMD_find_word(buffer, "NNB"))
  {
    ret_code = NAMD_read_line(psf_file, buffer);

    if (ret_code != 0)
    {
      NAMD_die("EOF ENCOUNTERED LOOKING FOR NNB IN PSF FILE");
    }
  }

  /*  Read in the number of exclusions and then the exclusions    */
  sscanf(buffer, "%d", &numExclusions);

  if (numExclusions)
    read_exclusions(psf_file);

  // DRUDE: read lone pair hosts and anisotropic terms from PSF
  if (is_drude_psf)
  {
    while (!NAMD_find_word(buffer, "NUMLP"))
    {
      ret_code = NAMD_read_line(psf_file, buffer);
      if (ret_code != 0)
      {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NUMLP IN DRUDE PSF FILE");
      }
    }
    sscanf(buffer, "%d", &numLphosts);
    if (numLphosts) read_lphosts(psf_file);

    while (!NAMD_find_word(buffer, "NUMANISO"))
    {
      ret_code = NAMD_read_line(psf_file, buffer);
      if (ret_code != 0)
      {
        NAMD_die("EOF ENCOUNTERED LOOKING FOR NUMANISO IN DRUDE PSF FILE");
      }
    }
    sscanf(buffer, "%d", &numAnisos);
    if (numAnisos) read_anisos(psf_file);

  }
  // DRUDE

  /*  look for the cross-term section.     */
  int crossterms_present = 1;
  while (!NAMD_find_word(buffer, "NCRTERM"))
  {
    ret_code = NAMD_read_line(psf_file, buffer);

    if (ret_code != 0)
    {
      // hit EOF before finding cross-term section
      crossterms_present = 0;
      break;
    }
  }

  if ( crossterms_present) {

    /*  Read in the number of cross-terms and then the cross-terms*/
    sscanf(buffer, "%d", &numCrossterms);

    if (numCrossterms)
      read_crossterms(psf_file, params);

  }

  /*  Close the .psf file.  */
  Fclose(psf_file);

  //  analyze the data and find the status of each atom
  numRealBonds = numBonds;
  build_atom_status();

  return;
#endif
}
/*      END OF FUNCTION read_psf_file      */

/************************************************************************/
/*                  */
/*        FUNCTION read_compressed_psf_file      */
/*                  */
/*   INPUTS:                */
/*  fname - Name of the compressed .psf file to read        */
/*  params - pointer to Parameters object to use to obtain          */
/*     parameters for vdWs, bonds, etc.      */
/*                  */
/*  This function reads a compressed .psf file in.  This is where just about   */
/*   all of the structural information for this class comes from.  The  */
/*   .psf file contains descriptions of the atom, bonds, angles,        */
/*   dihedrals, impropers, and exclusions.  The parameter object is     */
/*   used to look up parameters for each of these entities.    */
/*                  */
/************************************************************************/
void Molecule::read_compressed_psf_file(char *fname, Parameters *params, ConfigList *cfgList){
#ifndef MEM_OPT_VERSION
    return;
#else
    FILE *psf_file;    //  pointer to .psf file
    int ret_code;    //  ret_code from NAMD_read_line calls
    char buffer[512];

    /* Try and open the .psf file           */
    if ( (psf_file = Fopen(fname, "r")) == NULL)
    {
        char err_msg[512];
        sprintf(err_msg, "UNABLE TO OPEN THE COMPRESSED .psf FILE %s", fname);
        NAMD_die(err_msg);
    }     

    char strBuf[12];

    //first check compressed psf file format version
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "FORMAT VERSION")) {
        NAMD_die("The compressed psf file format is incorrect, please re-generate!\n");
    }
    float psfVer = 0.0f;
    sscanf(buffer, "FORMAT VERSION: %f\n", &psfVer);
    if(fabs(psfVer - COMPRESSED_PSF_VER)>1e-6) {
        NAMD_die("The compressed psf file format is incorrect, please re-generate!\n");
    }

    //Begin reading constant pools for atoms' basic information
    //1. segment names
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NSEGMENTNAMES"))
        NAMD_die("UNABLE TO FIND NSEGMENTNAMES");
    sscanf(buffer, "%d", &segNamePoolSize);
    if(segNamePoolSize!=0)
        segNamePool = new char *[segNamePoolSize];
    for(int i=0; i<segNamePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%s", strBuf);
        segNamePool[i] = nameArena->getNewArray(strlen(strBuf)+1);
        strcpy(segNamePool[i], strBuf);        
    }    

    //2. residue names
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NRESIDUENAMES"))
        NAMD_die("UNABLE TO FIND NRESIDUENAMES");
    sscanf(buffer, "%d", &resNamePoolSize);
    if(resNamePoolSize!=0)
        resNamePool = new char *[resNamePoolSize];
    for(int i=0; i<resNamePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%s", strBuf);
        resNamePool[i] = nameArena->getNewArray(strlen(strBuf)+1);
        strcpy(resNamePool[i], strBuf);
    }

    //3. atom names
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NATOMNAMES"))
        NAMD_die("UNABLE TO FIND NATOMNAMES");
    sscanf(buffer, "%d", &atomNamePoolSize);
    if(atomNamePoolSize!=0)
        atomNamePool = new char *[atomNamePoolSize];
    for(int i=0; i<atomNamePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%s", strBuf);
        atomNamePool[i] = nameArena->getNewArray(strlen(strBuf)+1);
        strcpy(atomNamePool[i], strBuf);
    }

    //4. atom types
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NATOMTYPES"))
        NAMD_die("UNABLE TO FIND NATOMTYPES");
    sscanf(buffer, "%d", &atomTypePoolSize);
    if(atomTypePoolSize!=0)
        atomTypePool = new char *[atomTypePoolSize];
    for(int i=0; i<atomTypePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%s", strBuf);
        atomTypePool[i] = nameArena->getNewArray(strlen(strBuf)+1);
        strcpy(atomTypePool[i], strBuf);
    }

    //5. charges
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NCHARGES"))
        NAMD_die("UNABLE TO FIND NCHARGES");
    sscanf(buffer, "%d", &chargePoolSize);
    if(chargePoolSize!=0)
        atomChargePool = new Real[chargePoolSize];
    for(int i=0; i<chargePoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%f", atomChargePool+i);        
    }

    //6. masses
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NMASSES"))
        NAMD_die("UNABLE TO FIND NMASSES");
    sscanf(buffer, "%d", &massPoolSize);
    if(massPoolSize!=0)
        atomMassPool = new Real[massPoolSize];
    for(int i=0; i<massPoolSize; i++){
        NAMD_read_line(psf_file, buffer);
        sscanf(buffer, "%f", atomMassPool+i);        
    }

    //Begin reading of atom signatures
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "ATOMSIGS"))
        NAMD_die("UNABLE TO FIND ATOMSIGS");
    sscanf(buffer, "%d", &atomSigPoolSize);
    atomSigPool = new AtomSignature[atomSigPoolSize];
    int typeCnt;
    int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    int tisReal;
    int ttype;
    for(int i=0; i<atomSigPoolSize; i++){
        //BOND SIGS
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NBONDSIGS"))
            NAMD_die("UNABLE TO FIND NBONDSIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].bondCnt = typeCnt;
            atomSigPool[i].bondSigs = new TupleSignature[typeCnt];
        }            
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);            
            sscanf(buffer, "%d | %d | %d", &tmp1, &ttype, &tisReal);
            TupleSignature oneSig(1, BOND, (Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;            
            atomSigPool[i].bondSigs[j]=oneSig;
            if(tisReal) numRealBonds++;
        }

        //ANGLE SIGS
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NTHETASIGS"))
            NAMD_die("UNABLE TO FIND NTHETASIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].angleCnt = typeCnt;
            atomSigPool[i].angleSigs = new TupleSignature[typeCnt];
        }            
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);            
            sscanf(buffer, "%d %d | %d | %d", &tmp1, &tmp2, &ttype, &tisReal);
            TupleSignature oneSig(2,ANGLE,(Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            oneSig.offset[1] = tmp2;            
            atomSigPool[i].angleSigs[j] = oneSig;
        }
        //DIHEDRAL SIGS
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NPHISIGS"))
            NAMD_die("UNABLE TO FIND NPHISIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].dihedralCnt = typeCnt;
            atomSigPool[i].dihedralSigs = new TupleSignature[typeCnt];
        }            
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);            
            sscanf(buffer, "%d %d %d | %d | %d", &tmp1, &tmp2, &tmp3, &ttype, &tisReal);
            TupleSignature oneSig(3,DIHEDRAL,(Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            oneSig.offset[1] = tmp2;
            oneSig.offset[2] = tmp3;
            atomSigPool[i].dihedralSigs[j] = oneSig;
        }
        //IMPROPER SIGS
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NIMPHISIGS"))
            NAMD_die("UNABLE TO FIND NIMPHISIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].improperCnt = typeCnt;
            atomSigPool[i].improperSigs = new TupleSignature[typeCnt];
        }            
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);            
            sscanf(buffer, "%d %d %d | %d | %d", &tmp1, &tmp2, &tmp3, &ttype, &tisReal);
            TupleSignature oneSig(3,IMPROPER,(Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            oneSig.offset[1] = tmp2;
            oneSig.offset[2] = tmp3;
            atomSigPool[i].improperSigs[j] = oneSig;
        }        
        //CROSSTERM SIGS
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NCRTERMSIGS"))
            NAMD_die("UNABLE TO FIND NCRTERMSIGS");
        sscanf(buffer, "%d", &typeCnt);
        if(typeCnt!=0){
            atomSigPool[i].crosstermCnt = typeCnt;
            atomSigPool[i].crosstermSigs = new TupleSignature[typeCnt];
        }            
        for(int j=0; j<typeCnt; j++){
            NAMD_read_line(psf_file, buffer);            
            sscanf(buffer, "%d %d %d %d %d %d %d | %d | %d", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &tmp7, &ttype, &tisReal);
            TupleSignature oneSig(7,CROSSTERM,(Index)ttype, (char)tisReal);
            oneSig.offset[0] = tmp1;
            oneSig.offset[1] = tmp2;
            oneSig.offset[2] = tmp3;
            oneSig.offset[3] = tmp4;
            oneSig.offset[4] = tmp5;
            oneSig.offset[5] = tmp6;
            oneSig.offset[6] = tmp7;
            atomSigPool[i].crosstermSigs[j] = oneSig;
        }        
    }

    //Reading exclusions' signatures    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NEXCLSIGS")){
        NAMD_die("UNABLE TO FIND NEXCLSIGS");
    }
    sscanf(buffer, "%d", &exclSigPoolSize);
    if(exclSigPoolSize>0) exclSigPool = new ExclusionSignature[exclSigPoolSize];
    vector<int> fullExcls;
    vector<int> modExcls;
    for(int i=0; i<exclSigPoolSize; i++){
        int fullExclCnt = NAMD_read_int(psf_file, buffer);
        for(int j=0; j<fullExclCnt; j++)
            fullExcls.push_back(NAMD_read_int(psf_file, buffer));
        int modExclCnt = NAMD_read_int(psf_file, buffer);
        for(int j=0; j<modExclCnt; j++)
            modExcls.push_back(NAMD_read_int(psf_file, buffer));

        //The integers in both vectors should be in the increasing the order
        //since they have been sorted during the compression of psf file
        exclSigPool[i].setOffsets(fullExcls, modExcls);

        fullExcls.clear();
        modExcls.clear();
    }

    //Now read the cluster information
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NCLUSTERS")) {
        NAMD_die("UNABLE TO FIND NCLUSTERS");
    }
    sscanf(buffer, "%d", &numClusters);

    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "CLUSTERCONTIGUITY")) {
        NAMD_die("UNABLE TO FIND CLUSTERCONTIGUITY");
    }
    sscanf(buffer, "%d", &isClusterContiguous);
    if(!isClusterContiguous) {
        clusterSize = new int32[numClusters];
        memset(clusterSize, 0, sizeof(int32)*numClusters);
    }

    //Now begin to read atoms
    //This part could be read in the batch mode. All the above information
    //can be stored in memory which only occupies several KBs.
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NATOM"))
        NAMD_die("UNABLE TO FIND NATOM");
    sscanf(buffer, "%d", &numAtoms);

    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "NHYDROGENGROUP"))
        NAMD_die("UNABLE TO FIND NHYDROGENGROUP");
    sscanf(buffer, "%d", &numHydrogenGroups);

    int isOccupancyValid, isBFactorValid;
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "OCCUPANCYVALID"))
        NAMD_die("UNABLE TO FIND OCCUPANCYVALID");
    sscanf(buffer, "%d", &isOccupancyValid);
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "TEMPFACTORVALID"))
        NAMD_die("UNABLE TO FIND TEMPFACTORVALID");
    sscanf(buffer, "%d", &isBFactorValid);

    if(numAtoms>0){
        atoms = new AtomCstInfo[numAtoms];
        atomNames = new AtomNameIdx[numAtoms];
        eachAtomMass = new Index[numAtoms];
        eachAtomCharge = new Index[numAtoms];
        eachAtomSig = new Index[numAtoms];
        eachAtomExclSig = new Index[numAtoms];

        clusterSigs = new int[numAtoms];

	int *tmpHgSize = new int[numAtoms];
	int *tmpHgGP = new int[numAtoms];
	int *tmpHgSV = new int[numAtoms];
        hydrogenGroup.resize(numAtoms);
        HydrogenGroupID *hg = hydrogenGroup.begin();

        ResidueLookupElem *tmpResLookup = resLookup;

        if(isOccupancyValid) {
            occupancy = new float[numAtoms];
        }
        if(isBFactorValid) {
            bfactor = new float[numAtoms];
        }

        int residue_number; //for residue number
        char *segment_name;
        int read_count;        
        float tmpf[2];
        
#if BINARY_PERATOM_OUTPUT
        //1. open the binary per-atom info file (fname.bin)
        char *binFName = new char[strlen(fname)+10];
        sprintf(binFName, "%s.bin", fname);
        FILE *perAtomFile = Fopen(binFName, "rb");
        if(perAtomFile==NULL) {
            char err_msg[512];
            sprintf(err_msg, "UNABLE TO OPEN THE ASSOCIATED PER-ATOM FILE FOR THE COMPRESSED .psf FILE %s", binFName);            
            NAMD_die(err_msg);
        }                
        //2. read the magic number to determine whether the endian is same or not
        int needFlip = 0;
        int magicNum = COMPRESSED_PSF_MAGICNUM;
        int rMagicNum = COMPRESSED_PSF_MAGICNUM;
        flipNum((char *)&rMagicNum, sizeof(int), 1);
        int fMagicNum;
        fread(&fMagicNum, sizeof(int), 1, perAtomFile);
        if(fMagicNum==magicNum) {
            needFlip = 0;
        }else if(fMagicNum==rMagicNum){
            needFlip = 1;
        }else{
            char err_msg[512];
            sprintf(err_msg, "THE ASSOCIATED PER-ATOM FILE FOR THE COMPRESSED .psf FILE %s IS CORRUPTED", binFName);            
            NAMD_die(err_msg);
        }
        //3. read the file version number
        float verNum =  0.0f;
        fread(&verNum, sizeof(float), 1, perAtomFile);
        if(needFlip) flipNum((char *)&verNum, sizeof(float), 1);
        if(fabs(verNum - COMPRESSED_PSF_VER)>1e-6) {
            char err_msg[512];
            sprintf(err_msg, "THE ASSOCIATED PER-ATOM FILE FOR THE COMPRESSED .psf FILE %s IS INCORRECT, PLEASE RE-GENERATE!\n", binFName);
            NAMD_die(err_msg);
        }
        delete[] binFName;
        //4. read per-atom info
        Index sIdx[8];
        int iIdx[8];
        for(int i=0; i<numAtoms; i++){   
            fread(sIdx, sizeof(Index), 8, perAtomFile);
            fread(iIdx, sizeof(int), 7, perAtomFile);
            fread(tmpf, sizeof(float), 2, perAtomFile);
            if(needFlip) {
                flipNum((char *)sIdx, sizeof(Index), 8);
                flipNum((char *)iIdx, sizeof(int), 7);
                flipNum((char *)tmpf, sizeof(float), 2);
            }
            segment_name = segNamePool[sIdx[0]];                   
            atomNames[i].resnameIdx = sIdx[1];
            atomNames[i].atomnameIdx = sIdx[2];
            atomNames[i].atomtypeIdx = sIdx[3];
            eachAtomCharge[i] = sIdx[4];
            eachAtomMass[i] = sIdx[5];
            eachAtomSig[i] = sIdx[6];
            eachAtomExclSig[i] = sIdx[7];
            residue_number = iIdx[0];
            clusterSigs[i] = iIdx[1];
            if(!isClusterContiguous)
                clusterSize[iIdx[1]]++;

            atoms[i].partner = iIdx[2];
            atoms[i].hydrogenList= iIdx[3];

            tmpHgSize[i] = iIdx[4];
            tmpHgGP[i] = iIdx[5];
            tmpHgSV[i] = iIdx[6];

            //debugExclNum += (exclSigPool[sIdx[7]].fullExclCnt+exclSigPool[sIdx[7]].modExclCnt);
#else
        int idx[15];
        for(int i=0; i<numAtoms; i++){
            NAMD_read_line(psf_file, buffer);
            read_count = sscanf(buffer, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f",
                                idx, idx+1, idx+2, idx+3, idx+4,
                                idx+5, idx+6, idx+7, idx+8, idx+9, idx+10, idx+11, idx+12,
                                idx+13, idx+14, tmpf, tmpf+1);
            if(read_count!=17){
                char err_msg[128];
                sprintf(err_msg, "BAD ATOM LINE FORMAT IN COMPRESSED PSF FILE IN ATOM LINE %d\nLINE=%s",i+1, buffer);
                NAMD_die(err_msg);
            }
            segment_name = segNamePool[idx[0]];                   
            atomNames[i].resnameIdx = (Index)idx[1];
            atomNames[i].atomnameIdx = (Index)idx[2];
            atomNames[i].atomtypeIdx = (Index)idx[3];
            eachAtomCharge[i] = (Index)idx[4];
            eachAtomMass[i] = (Index)idx[5];              
            eachAtomSig[i] = (Index)idx[6];
            eachAtomExclSig[i] = (Index)idx[7];
            residue_number = idx[8];
            clusterSigs[i] = idx[9];
            if(!isClusterContiguous)
                clusterSize[idx[9]]++;

            atoms[i].partner = idx[10];
            atoms[i].hydrogenList= idx[11];

            tmpHgSize[i] = idx[12];
            tmpHgGP[i] = idx[13];
            tmpHgSV[i] = idx[14];
            //debugExclNum += (exclSigPool[idx[7]].fullExclCnt+exclSigPool[idx[7]].modExclCnt);
#endif	    

            if(isOccupancyValid)
                occupancy[i] = tmpf[0];
            if(isBFactorValid)
                bfactor[i] = tmpf[1];

            //Add this atom to residue lookup table
            if(tmpResLookup) tmpResLookup =
                tmpResLookup->append(segment_name, residue_number, i);
            //Determine the type of the atom (H or O)
            Real thisAtomMass = atomMassPool[eachAtomMass[i]];
            atoms[i].status = UnknownAtom;
            if (thisAtomMass <= 0.05) {
              atoms[i].status |= LonepairAtom;
            } else if (thisAtomMass < 1.0) {
              atoms[i].status |= DrudeAtom;
            } else if(thisAtomMass <= 3.5){
                atoms[i].status |= HydrogenAtom;
            }else if(atomNamePool[atomNames[i].atomnameIdx][0]=='O' &&
                     (thisAtomMass >= 14.0) && (thisAtomMass <= 18.0)){
                atoms[i].status |= OxygenAtom;
            }

            //Look up the vdw constants for this atom
            params->assign_vdw_index(atomTypePool[atomNames[i].atomtypeIdx],
                                     &(atoms[i]));
        } //end of reading per-atom information

	//re-sort the hydrogenGroup array based on the recorded "hydrogenList"
	//value
	for(int i=0; i<numAtoms; i++){
	    int hgIdx = atoms[i].hydrogenList;
	    hg[hgIdx].atomID = i;
	    hg[hgIdx].atomsInGroup = tmpHgSize[i];
	    hg[hgIdx].GPID = tmpHgGP[i];
	    hg[hgIdx].waterVal = tmpHgSV[i];
	    hg[hgIdx].isGP = 1;
	    if(tmpHgSize[i]==0) hg[hgIdx].isGP = 0;
	}
	
	delete [] tmpHgSize;
	delete [] tmpHgGP;
	delete [] tmpHgSV;	

	//printf("debugExclNum: %d\n", debugExclNum);
    }
    //build up information for numBonds/numDihedrals.. etc
    numBonds=0;
    numAngles=0;
    numDihedrals=0;
    numImpropers=0;
    numCrossterms=0;
    numTotalExclusions=0;
    for(int i=0; i<numAtoms; i++){
        AtomSignature *thisSig = &atomSigPool[eachAtomSig[i]];
        numBonds += thisSig->bondCnt;
        numAngles += thisSig->angleCnt;
        numDihedrals += thisSig->dihedralCnt;
        numImpropers += thisSig->improperCnt;
        numCrossterms += thisSig->crosstermCnt;

        ExclusionSignature *exclSig = &exclSigPool[eachAtomExclSig[i]];
        numTotalExclusions += (exclSig->fullExclCnt + exclSig->modExclCnt);
    }
    
    numTotalExclusions /= 2;

    //Just reading for the parameters values; extra Bonds, Dihedrals etc.
    //have been taken into account when compressing the molecule object.
    //The actual number of Bonds, Dihedrals etc. will be calculated based
    //on atom signatures.
    if(simParams->extraBondsOn)
	build_extra_bonds(params, cfgList->find("extraBondsFile"));
#if 0
    //This part has been enabled in build_extra_bonds for memory optimized version
    //read extra bond parameters if there is an input of extra bonds (extraBondsOn is true)
    int extraDihedralParamNum = 0;
    int extraImproperParamNum = 0;
    if(simParams->extraBondsOn){
        int numExtraParams=0;

        //read extra bond params
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NEXTRABONDPARAMS")){
            NAMD_die("UNABLE TO FIND NEXTRABONDPARAMS");
        }
        sscanf(buffer, "%d", &numExtraParams);
        if(numExtraParams>0){
            BondValue *newParams = new BondValue[params->NumBondParams+numExtraParams];
            memcpy(newParams, params->bond_array, params->NumBondParams*sizeof(BondValue));
            delete [] params->bond_array;            

            int curNumPs = params->NumBondParams;
            for(int i=0; i<numExtraParams; i++){
                Real k, x0;
                NAMD_read_line(psf_file, buffer);
                sscanf(buffer, "%f %f", &k, &x0);
                newParams[curNumPs+i].k = k;
                newParams[curNumPs+i].x0 = x0;
            }
            params->bond_array = newParams;
            params->NumBondParams += numExtraParams;
        }

        //read extra angle params
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NEXTRAANGLEPARAMS")){
            NAMD_die("UNABLE TO FIND NEXTRAANGLEPARAMS");
        }
        sscanf(buffer, "%d", &numExtraParams);
        if(numExtraParams>0){
            AngleValue *newParams = new AngleValue[params->NumAngleParams+numExtraParams];
            memcpy(newParams, params->angle_array, params->NumAngleParams*sizeof(AngleValue));
            delete [] params->angle_array;

            int curNumPs = params->NumAngleParams;
            for(int i=0; i<numExtraParams; i++) {
                Real k, theta0, k_ub, r_ub;
                NAMD_read_line(psf_file, buffer);
                sscanf(buffer, "%f %f %f %f", &k, &theta0, &k_ub, &r_ub);
                newParams[curNumPs+i].k = k;
                newParams[curNumPs+i].theta0 = theta0;
                newParams[curNumPs+i].k_ub = k_ub;
                newParams[curNumPs+i].r_ub = r_ub;
            }
            params->angle_array = newParams;
            params->NumAngleParams += numExtraParams;
        }
        
        //read extra diheral params
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NEXTRADIHEDRALPARAMS")){
            NAMD_die("UNABLE TO FIND NEXTRADIHEDRALPARAMS");
        }
        sscanf(buffer, "%d", &numExtraParams);
        extraDihedralParamNum = numExtraParams;
        if(numExtraParams>0){
            DihedralValue *newParams = new DihedralValue[params->NumDihedralParams+numExtraParams];
            memcpy(newParams, params->dihedral_array, params->NumDihedralParams*sizeof(DihedralValue));
            delete [] params->dihedral_array;

            int curNumPs = params->NumDihedralParams;
            for(int i=0; i<numExtraParams; i++) {
                int multiplicity;
                Real k, delta;
                int n;
                NAMD_read_line(psf_file, buffer); //read multiplicity
                sscanf(buffer, "%d", &multiplicity);
                newParams[curNumPs+i].multiplicity = multiplicity;
                for(int j=0; j<multiplicity; j++) {
                    NAMD_read_line(psf_file, buffer);
                    sscanf(buffer, "%f %d %f", &k, &n, &delta);
                    newParams[curNumPs+i].values[j].k = k;
                    newParams[curNumPs+i].values[j].n = n;
                    newParams[curNumPs+i].values[j].delta = delta;
                }
            }
            params->dihedral_array = newParams;
        }

        //read extra improper params
        NAMD_read_line(psf_file, buffer);
        if(!NAMD_find_word(buffer, "NEXTRAIMPROPERPARAMS")){
            NAMD_die("UNABLE TO FIND NEXTRAIMPROPERPARAMS");
        }
        sscanf(buffer, "%d", &numExtraParams);
        extraImproperParamNum = numExtraParams;
        if(numExtraParams>0){
            ImproperValue *newParams = new ImproperValue[params->NumImproperParams+numExtraParams];
            memcpy(newParams, params->improper_array, params->NumImproperParams*sizeof(ImproperValue));
            delete [] params->improper_array;

            int curNumPs = params->NumImproperParams;
            for(int i=0; i<numExtraParams; i++) {
                int multiplicity;
                Real k, delta;
                int n;
                NAMD_read_line(psf_file, buffer); //read multiplicity
                sscanf(buffer, "%d", &multiplicity);
                newParams[curNumPs+i].multiplicity = multiplicity;
                for(int j=0; j<multiplicity; j++) {
                    NAMD_read_line(psf_file, buffer);
                    sscanf(buffer, "%f %d %f", &k, &n, &delta);
                    newParams[curNumPs+i].values[j].k = k;
                    newParams[curNumPs+i].values[j].n = n;
                    newParams[curNumPs+i].values[j].delta = delta;
                }
            }
            params->improper_array = newParams;
        }
    }
#endif

    //read DIHEDRALPARAMARRAY and IMPROPERPARAMARRAY    
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "DIHEDRALPARAMARRAY"))
        NAMD_die("UNABLE TO FIND DIHEDRALPARAMARRAY");
    for(int i=0; i<params->NumDihedralParams; i++){
        params->dihedral_array[i].multiplicity = NAMD_read_int(psf_file, buffer);
    }
    //params->NumDihedralParams += extraDihedralParamNum;

    NAMD_read_line(psf_file, buffer); //to read a simple single '\n' line 
    NAMD_read_line(psf_file, buffer);
    if(!NAMD_find_word(buffer, "IMPROPERPARAMARRAY"))
        NAMD_die("UNABLE TO FIND IMPROPERPARAMARRAY");
    for(int i=0; i<params->NumImproperParams; i++){
        params->improper_array[i].multiplicity = NAMD_read_int(psf_file, buffer);
    }
    //params->NumImproperParams += extraImproperParamNum;
    Fclose(psf_file);

    //numRealBonds is set when reading bondSignatures from compressed psf file
    
    build_atom_status();
#endif
}
/*      END OF FUNCTION read_compressed_psf_file      */

/************************************************************************/
/*                  */
/*        FUNCTION read_atoms      */
/*                  */
/*   INPUTS:                */
/*  fd - file pointer to the .psf file        */
/*  params - Parameters object to use for parameters    */
/*                  */
/*  this function reads in the Atoms section of the .psf file.      */
/*   This section consists of numAtoms lines that are of the form:  */
/*     <atom#> <mol> <seg#> <res> <atomname> <atomtype> <charge> <mass> */
/*   Each line is read into the appropriate entry in the atoms array.   */
/*   The parameters object is then used to determine the vdW constants  */
/*   for this atom.              */
/*                  */
/************************************************************************/

void Molecule::read_atoms(FILE *fd, Parameters *params)

{
#ifdef MEM_OPT_VERSION
    return;
#else
  char buffer[512];  // Buffer for reading from file
  int atom_number=0;  // Atom number 
  int last_atom_number=0; // Last atom number, used to assure
        // atoms are in order
  char segment_name[11]; // Segment name
  char residue_number[11]; // Residue number
  char residue_name[11];  // Residue name
  char atom_name[11];  // Atom name
  char atom_type[11];  // Atom type
  Real charge;    // Charge for the current atom
  Real mass;    // Mass for the current atom
  int read_count;    // Number of fields read by sscanf

  /*  Allocate the atom arrays          */
  atoms     = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];
  if(simParams->genCompressedPsf) {
      atomSegResids = new AtomSegResInfo[numAtoms];
  }

  // DRUDE: supplement Atom data
  if (is_drude_psf) {
    drudeConsts = new DrudeConst[numAtoms];
  }
  // DRUDE

  hydrogenGroup.resize(0);

  if (atoms == NULL || atomNames == NULL )
  {
    NAMD_die("memory allocation failed in Molecule::read_atoms");
  }

  ResidueLookupElem *tmpResLookup = resLookup;

  /*  Loop and read in numAtoms atom lines.      */
  while (atom_number < numAtoms)
  {
    /*  Get the line from the file        */
    NAMD_read_line(fd, buffer);

    /*  If its blank or a comment, skip it      */
    if ( (NAMD_blank_string(buffer)) || (buffer[0] == '!') )
      continue;

    /*  Parse up the line          */
    read_count=sscanf(buffer, "%d %s %s %s %s %s %f %f",
       &atom_number, segment_name, residue_number,
       residue_name, atom_name, atom_type, &charge, &mass);

    /*  Check to make sure we found what we were expecting  */
    if (read_count != 8)
    {
      char err_msg[128];

      sprintf(err_msg, "BAD ATOM LINE FORMAT IN PSF FILE IN ATOM LINE %d\nLINE=%s",
         last_atom_number+1, buffer);
      NAMD_die(err_msg);
    }

    // DRUDE: read alpha and thole parameters from atom line
    if (is_drude_psf)
    {
      Real alpha, thole;
      read_count=sscanf(buffer,
//          "%*d %*s %*s %*s %*s %*s %*f %*f %*d %*f %*f %f %f", &alpha, &thole);
                // the two columns preceding alpha and thole will disappear
          "%*d %*s %*s %*s %*s %*s %*f %*f %*d %f %f", &alpha, &thole);
      if (read_count != 2)
      {
        char err_msg[128];

        sprintf(err_msg, "BAD ATOM LINE FORMAT IN PSF FILE "
            "IN ATOM LINE %d\nLINE=%s", last_atom_number+1, buffer);
        NAMD_die(err_msg);
      }
      drudeConsts[atom_number-1].alpha = alpha;
      drudeConsts[atom_number-1].thole = thole;
    }
    // DRUDE

    /*  Check if this is in XPLOR format  */
    int atom_type_num;
    if ( sscanf(atom_type, "%d", &atom_type_num) > 0 )
    {
      NAMD_die("Structure (psf) file is either in CHARMM format (with numbers for atoms types, the X-PLOR format using names is required) or the segment name field is empty.");
    }

    /*  Make sure the atoms were in sequence    */
    if (atom_number != last_atom_number+1)
    {
      char err_msg[128];

      sprintf(err_msg, "ATOM NUMBERS OUT OF ORDER AT ATOM #%d OF PSF FILE",
         last_atom_number+1);
      NAMD_die(err_msg);
    }

    last_atom_number++;

    /*  Dynamically allocate strings for atom name, atom    */
    /*  type, etc so that we only allocate as much space    */
    /*  for these strings as we really need      */
    int reslength = strlen(residue_name)+1;
    int namelength = strlen(atom_name)+1;
    int typelength = strlen(atom_type)+1;

    atomNames[atom_number-1].resname = nameArena->getNewArray(reslength);
    atomNames[atom_number-1].atomname = nameArena->getNewArray(namelength);
    atomNames[atom_number-1].atomtype = nameArena->getNewArray(typelength);
  
    if (atomNames[atom_number-1].resname == NULL)
    {
      NAMD_die("memory allocation failed in Molecule::read_atoms");
    }

    /*  Put the values from this atom into the atoms array  */
    strcpy(atomNames[atom_number-1].resname, residue_name);
    strcpy(atomNames[atom_number-1].atomname, atom_name);
    strcpy(atomNames[atom_number-1].atomtype, atom_type);
    atoms[atom_number-1].mass = mass;
    atoms[atom_number-1].charge = charge;
    atoms[atom_number-1].status = UnknownAtom;

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append(segment_name, atoi(residue_number), atom_number-1);

    if(atomSegResids) { //for compressing molecule information
        AtomSegResInfo *one = atomSegResids + (atom_number - 1);
        memcpy(one->segname, segment_name, strlen(segment_name)+1);
        one->resid = atoi(residue_number);
    }

    /*  Determine the type of the atom (H or O) */
    if (atoms[atom_number-1].mass <= 0.05) {
      atoms[atom_number-1].status |= LonepairAtom;
    } else if (atoms[atom_number-1].mass < 1.0) {
      atoms[atom_number-1].status |= DrudeAtom;
    } else if (atoms[atom_number-1].mass <=3.5) {
      atoms[atom_number-1].status |= HydrogenAtom;
    } else if ((atomNames[atom_number-1].atomname[0] == 'O') && 
         (atoms[atom_number-1].mass >= 14.0) && 
         (atoms[atom_number-1].mass <= 18.0)) {
      atoms[atom_number-1].status |= OxygenAtom;
    }

    /*  Look up the vdw constants for this atom    */
    params->assign_vdw_index(atomNames[atom_number-1].atomtype, 
       &(atoms[atom_number-1]));
        }

  return;
#endif
}
/*      END OF FUNCTION read_atoms      */

/************************************************************************/
/*                  */
/*      FUNCTION read_bonds        */
/*                  */
/*  read_bonds reads in the bond section of the .psf file.  This    */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are bonded together.  Each atom pair is   */
/*  read in.  Then that parameter object is queried to determine the    */
/*  force constant and rest distance for the bond.      */
/*                  */
/************************************************************************/

void Molecule::read_bonds(FILE *fd, Parameters *params)

{
#ifdef MEM_OPT_VERSION
    return;
#else
  int atom_nums[2];  // Atom indexes for the bonded atoms
  char atom1name[11];  // Atom type for atom #1
  char atom2name[11];  // Atom type for atom #2
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
  int origNumBonds = numBonds;   // number of bonds in file header

  /*  Allocate the array to hold the bonds      */
  bonds=new Bond[numBonds];

  if (bonds == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_bonds");
  }

  /*  Loop through and read in all the bonds      */
  while (num_read < numBonds)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "BONDS")-1;

      /*  Check to make sure the index isn't too big  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "BOND INDEX %d GREATER THAN NATOM %d IN BOND # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }
    }

    /*  Get the atom type for the two atoms.  When we query */
    /*  the parameter object, we need to send the atom type */
    /*  that is alphabetically first as atom 1.    */
    if (strcasecmp(atomNames[atom_nums[0]].atomtype, 
         atomNames[atom_nums[1]].atomtype) < 0)
    {
      strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
      strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
    }
    else
    {
      strcpy(atom2name, atomNames[atom_nums[0]].atomtype);
      strcpy(atom1name, atomNames[atom_nums[1]].atomtype);
    }

    /*  Query the parameter object for the constants for    */
    /*  this bond            */
    Bond *b = &(bonds[num_read]);
    params->assign_bond_index(atom1name, atom2name, b);

    /*  Assign the atom indexes to the array element  */
    b->atom1=atom_nums[0];
    b->atom2=atom_nums[1];

    /*  Make sure this isn't a fake bond meant for shake in x-plor.  */
    Real k, x0;
    params->get_bond_params(&k,&x0,b->bond_type);
    if (simParams->watmodel == WAT_SWM4) {
      // need to retain Lonepair bonds for Drude
      if ( k == 0. && !is_lp(b->atom1) && !is_lp(b->atom2)) --numBonds;
      else ++num_read;
    }
    else {
      if ( k == 0. ) --numBonds;  // fake bond
      else ++num_read;  // real bond
    }
  }

  /*  Tell user about our subterfuge  */
  if ( numBonds != origNumBonds ) {
    iout << iWARN << "Ignored " << origNumBonds - numBonds <<
            " bonds with zero force constants.\n" << endi;
    iout << iWARN <<
	"Will get H-H distance in rigid H2O from H-O-H angle.\n" << endi;
  }

  return;
#endif
}
/*      END OF FUNCTION read_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION read_angles        */
/*                  */
/*   INPUTS:                */
/*  fd - File descriptor for .psf file        */
/*  params - Parameters object to query for parameters    */
/*                  */
/*  read_angles reads the angle parameters from the .psf file.      */
/*   This section of the .psf file consists of a list of triplets of    */
/*   atom indexes.  Each triplet represents three atoms connected via   */
/*   an angle bond.  The parameter object is queried to obtain the      */
/*   constants for each bond.            */
/*                  */
/************************************************************************/

void Molecule::read_angles(FILE *fd, Parameters *params)

{
#ifdef MEM_OPT_VERSION
    return;
#else
  int atom_nums[3];  //  Atom numbers for the three atoms
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  register int j;      //  Loop counter
  int num_read=0;    //  Number of angles read so far
  int origNumAngles = numAngles;  // Number of angles in file
  /*  Alloc the array of angles          */
  angles=new Angle[numAngles];

  if (angles == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_angles");
  }

  /*  Loop through and read all the angles      */
  while (num_read < numAngles)
  {
    /*  Loop through the 3 atom indexes in the current angle*/
    for (j=0; j<3; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "ANGLES")-1;

      /*  Check to make sure the atom index doesn't   */
      /*  exceed the Number of Atoms      */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "ANGLES INDEX %d GREATER THAN NATOM %d IN ANGLES # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }
    }

    /*  Place the bond name that is alphabetically first  */
    /*  in the atom1name.  This is OK since the order of    */
    /*  atom1 and atom3 are interchangable.  And to search  */
    /*  the tree of angle parameters, we need the order     */
    /*  to be predictable.          */
    if (strcasecmp(atomNames[atom_nums[0]].atomtype, 
         atomNames[atom_nums[2]].atomtype) < 0)
    {
      strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
      strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
      strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
    }
    else
    {
      strcpy(atom1name, atomNames[atom_nums[2]].atomtype);
      strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
      strcpy(atom3name, atomNames[atom_nums[0]].atomtype);
    }

    /*  Get the constant values for this bond from the  */
    /*  parameter object          */
    params->assign_angle_index(atom1name, atom2name, 
       atom3name, &(angles[num_read]));

    /*  Assign the three atom indices      */
    angles[num_read].atom1=atom_nums[0];
    angles[num_read].atom2=atom_nums[1];
    angles[num_read].atom3=atom_nums[2];

    /*  Make sure this isn't a fake angle meant for shake in x-plor.  */
    Real k, t0, k_ub, r_ub;
    params->get_angle_params(&k,&t0,&k_ub,&r_ub,angles[num_read].angle_type);
    if ( k == 0. && k_ub == 0. ) --numAngles;  // fake angle
    else ++num_read;  // real angle
  }

  /*  Tell user about our subterfuge  */
  if ( numAngles != origNumAngles ) {
    iout << iWARN << "Ignored " << origNumAngles - numAngles <<
            " angles with zero force constants.\n" << endi;
  }

  return;
#endif
}
/*      END OF FUNCTION read_angles      */

/************************************************************************/
/*                  */
/*        FUNCTION read_dihedrals      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for the .psf file        */
/*  params - pointer to parameter object        */
/*                  */
/*  read_dihedreals reads the dihedral section of the .psf file.    */
/*   This section of the file contains a list of quartets of atom       */
/*   numbers.  Each quartet represents a group of atoms that form a     */
/*   dihedral bond.              */
/*                  */
/************************************************************************/

void Molecule::read_dihedrals(FILE *fd, Parameters *params)
{
#ifdef MEM_OPT_VERSION
    return;
#else
  int atom_nums[4];  // The 4 atom indexes
  int last_atom_nums[4];  // Atom numbers from previous bond
  char atom1name[11];  // Atom type for atom 1
  char atom2name[11];  // Atom type for atom 2
  char atom3name[11];  // Atom type for atom 3
  char atom4name[11];  // Atom type for atom 4
  register int j;      // loop counter
  int num_read=0;    // number of dihedrals read so far
  int multiplicity=1;  // multiplicity of the current bond
  Bool duplicate_bond;  // Is this a duplicate of the last bond
  int num_unique=0;   // Number of unique dihedral bonds

  //  Initialize the array used to check for duplicate dihedrals
  for (j=0; j<4; j++)
    last_atom_nums[j] = -1;

  /*  Allocate an array to hold the Dihedrals      */
  dihedrals = new Dihedral[numDihedrals];

  if (dihedrals == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_dihedrals");
  }

  /*  Loop through and read all the dihedrals      */
  while (num_read < numDihedrals)
  {
    duplicate_bond = TRUE;

    /*  Loop through and read the 4 indexes for this bond   */
    for (j=0; j<4; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "DIHEDRALS")-1;

      /*  Check for an atom index that is too large  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "DIHEDRALS INDEX %d GREATER THAN NATOM %d IN DIHEDRALS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      //  Check to see if this atom matches the last bond
      if (atom_nums[j] != last_atom_nums[j])
      {
        duplicate_bond = FALSE;
      }

      last_atom_nums[j] = atom_nums[j];
    }

    /*  Get the atom types for the 4 atoms so we can look  */
    /*  up the constants in the parameter object    */
    strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
    strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
    strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
    strcpy(atom4name, atomNames[atom_nums[3]].atomtype);

    //  Check to see if this is really a new bond or just
    //  a repeat of the last one
    if (duplicate_bond)
    {
      //  This is a duplicate, so increase the multiplicity
      multiplicity++;

      if (multiplicity == 2)
      {
        numMultipleDihedrals++;
      }
    }
    else
    {
      multiplicity=1;
      num_unique++;
    }

    /*  Get the constants for this dihedral bond    */
    params->assign_dihedral_index(atom1name, atom2name, 
       atom3name, atom4name, &(dihedrals[num_unique-1]),
       multiplicity);

    /*  Assign the atom indexes        */
    dihedrals[num_unique-1].atom1=atom_nums[0];
    dihedrals[num_unique-1].atom2=atom_nums[1];
    dihedrals[num_unique-1].atom3=atom_nums[2];
    dihedrals[num_unique-1].atom4=atom_nums[3];

    num_read++;
  }

  numDihedrals = num_unique;

  return;
#endif
}
/*      END OF FUNCTION read_dihedral      */

/************************************************************************/
/*                  */
/*        FUNCTION read_impropers      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*  params - parameter object          */
/*                  */
/*  read_impropers reads the improper section of the .psf file.  */
/*   This section is identical to the dihedral section in that it is    */
/*   made up of a list of quartets of atom indexes that define the      */
/*   atoms that are bonded together.          */
/*                  */
/************************************************************************/

void Molecule::read_impropers(FILE *fd, Parameters *params)
{
#ifdef MEM_OPT_VERSION
    return;
#else
  int atom_nums[4];  //  Atom indexes for the 4 atoms
  int last_atom_nums[4];  //  Atom indexes from previous bond
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  char atom4name[11];  //  Atom type for atom 4
  register int j;      //  Loop counter
  int num_read=0;    //  Number of impropers read so far
  int multiplicity=1;  // multiplicity of the current bond
  Bool duplicate_bond;  // Is this a duplicate of the last bond
  int num_unique=0;   // Number of unique dihedral bonds

  //  Initialize the array used to look for duplicate improper
  //  entries.  Set them all to -1 so we know nothing will match
  for (j=0; j<4; j++)
    last_atom_nums[j] = -1;

  /*  Allocate the array to hold the impropers      */
  impropers=new Improper[numImpropers];

  if (impropers == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_impropers");
  }

  /*  Loop through and read all the impropers      */
  while (num_read < numImpropers)
  {
    duplicate_bond = TRUE;

    /*  Loop through the 4 indexes for this improper  */
    for (j=0; j<4; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "IMPROPERS")-1;

      /*  Check to make sure the index isn't too big  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "IMPROPERS INDEX %d GREATER THAN NATOM %d IN IMPROPERS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      if (atom_nums[j] != last_atom_nums[j])
      {
        duplicate_bond = FALSE;
      }

      last_atom_nums[j] = atom_nums[j];
    }

    /*  Get the atom types so we can look up the parameters */
    strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
    strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
    strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
    strcpy(atom4name, atomNames[atom_nums[3]].atomtype);

    //  Check to see if this is a duplicate improper
    if (duplicate_bond)
    {
      //  This is a duplicate improper.  So we don't
      //  really count this entry, we just update
      //  the parameters object
      multiplicity++;

      if (multiplicity == 2)
      {
        //  Count the number of multiples.
        numMultipleImpropers++;
      }
    }
    else
    {
      //  Not a duplicate
      multiplicity = 1;
      num_unique++;
    }

    /*  Look up the constants for this bond      */
    params->assign_improper_index(atom1name, atom2name, 
       atom3name, atom4name, &(impropers[num_unique-1]),
       multiplicity);

    /*  Assign the atom indexes        */
    impropers[num_unique-1].atom1=atom_nums[0];
    impropers[num_unique-1].atom2=atom_nums[1];
    impropers[num_unique-1].atom3=atom_nums[2];
    impropers[num_unique-1].atom4=atom_nums[3];

    num_read++;
  }

  //  Now reset the numImpropers value to the number of UNIQUE
  //  impropers.  Sure, we waste a few entries in the improper_array
  //  on the master node, but it is very little space . . .
  numImpropers = num_unique;

  return;
#endif
}
/*      END OF FUNCTION read_impropers      */

/************************************************************************/
/*                  */
/*        FUNCTION read_crossterms      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*  params - parameter object          */
/*                  */
/*   This section is identical to the dihedral section in that it is    */
/*   made up of a list of quartets of atom indexes that define the      */
/*   atoms that are bonded together.          */
/*                  */
/************************************************************************/

void Molecule::read_crossterms(FILE *fd, Parameters *params)

{
#ifdef MEM_OPT_VERSION
    return;
#else
  int atom_nums[8];  //  Atom indexes for the 4 atoms
  int last_atom_nums[8];  //  Atom indexes from previous bond
  char atom1name[11];  //  Atom type for atom 1
  char atom2name[11];  //  Atom type for atom 2
  char atom3name[11];  //  Atom type for atom 3
  char atom4name[11];  //  Atom type for atom 4
  char atom5name[11];  //  Atom type for atom 5
  char atom6name[11];  //  Atom type for atom 6
  char atom7name[11];  //  Atom type for atom 7
  char atom8name[11];  //  Atom type for atom 8
  register int j;      //  Loop counter
  int num_read=0;    //  Number of items read so far
  Bool duplicate_bond;  // Is this a duplicate of the last bond

  //  Initialize the array used to look for duplicate crossterm
  //  entries.  Set them all to -1 so we know nothing will match
  for (j=0; j<8; j++)
    last_atom_nums[j] = -1;

  /*  Allocate the array to hold the cross-terms */
  crossterms=new Crossterm[numCrossterms];

  if (crossterms == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_crossterms");
  }

  /*  Loop through and read all the cross-terms      */
  while (num_read < numCrossterms)
  {
    duplicate_bond = TRUE;

    /*  Loop through the 8 indexes for this cross-term */
    for (j=0; j<8; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      atom_nums[j]=NAMD_read_int(fd, "CROSS-TERMS")-1;

      /*  Check to make sure the index isn't too big  */
      if (atom_nums[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "CROSS-TERM INDEX %d GREATER THAN NATOM %d IN CROSS-TERMS # %d IN PSF FILE", atom_nums[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      if (atom_nums[j] != last_atom_nums[j])
      {
        duplicate_bond = FALSE;
      }

      last_atom_nums[j] = atom_nums[j];
    }

    /*  Get the atom types so we can look up the parameters */
    strcpy(atom1name, atomNames[atom_nums[0]].atomtype);
    strcpy(atom2name, atomNames[atom_nums[1]].atomtype);
    strcpy(atom3name, atomNames[atom_nums[2]].atomtype);
    strcpy(atom4name, atomNames[atom_nums[3]].atomtype);
    strcpy(atom5name, atomNames[atom_nums[4]].atomtype);
    strcpy(atom6name, atomNames[atom_nums[5]].atomtype);
    strcpy(atom7name, atomNames[atom_nums[6]].atomtype);
    strcpy(atom8name, atomNames[atom_nums[7]].atomtype);

    //  Check to see if this is a duplicate term
    if (duplicate_bond)
    {
      iout << iWARN << "Duplicate cross-term detected.\n" << endi;
    }

    /*  Look up the constants for this bond      */
    params->assign_crossterm_index(atom1name, atom2name, 
       atom3name, atom4name, atom5name, atom6name,
       atom7name, atom8name, &(crossterms[num_read]));

    /*  Assign the atom indexes        */
    crossterms[num_read].atom1=atom_nums[0];
    crossterms[num_read].atom2=atom_nums[1];
    crossterms[num_read].atom3=atom_nums[2];
    crossterms[num_read].atom4=atom_nums[3];
    crossterms[num_read].atom5=atom_nums[4];
    crossterms[num_read].atom6=atom_nums[5];
    crossterms[num_read].atom7=atom_nums[6];
    crossterms[num_read].atom8=atom_nums[7];

    if(!duplicate_bond) num_read++;
  }

  numCrossterms = num_read;

  return;
#endif
}
/*      END OF FUNCTION read_impropers      */

/************************************************************************/
/*                  */
/*      FUNCTION read_donors        */
/*                  */
/*  read_donors reads in the bond section of the .psf file.  This   */
/*  section contains a list of pairs of numbers where each pair is      */
/*  represents two atoms that are part of an H-bond.  Each atom pair is */
/*  read in.                                                            */
/*                  */
/*  Donor atoms are the heavy atoms to which hydrogens are bonded.      */
/*  There will always be a donor atom for each donor pair.  However,    */
/*  for a united-atom model there may not be an explicit hydrogen       */
/*  present, in which case the second atom index in the pair will be    */
/*  given as 0 in the PSF (and stored as -1 in this program's internal  */
/*  storage).                                                           */
/************************************************************************/

void Molecule::read_donors(FILE *fd)
{
#ifdef MEM_OPT_VERSION
    return;
#else
  int d[2];               // temporary storage of donor atom index
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
  int num_no_hydr=0;      // Number of bonds with no hydrogen given

  /*  Allocate the array to hold the bonds      */
  donors=new Bond[numDonors];

  if (donors == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_donors");
  }

  /*  Loop through and read in all the donors      */
  while (num_read < numDonors)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      d[j]=NAMD_read_int(fd, "DONORS")-1;

      /*  Check to make sure the index isn't too big  */
      if (d[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg,
    "DONOR INDEX %d GREATER THAN NATOM %d IN DONOR # %d IN PSF FILE",
          d[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      /*  Check if there is a hydrogen given */
      if (d[j] < 0)
                          num_no_hydr++;
    }

    /*  Assign the atom indexes to the array element  */
    Bond *b = &(donors[num_read]);
    b->atom1=d[0];
    b->atom2=d[1];

    num_read++;
  }

  return;
#endif
}
/*      END OF FUNCTION read_donors      */


/************************************************************************/
/*                  */
/*      FUNCTION read_acceptors        */
/*                  */
/*  read_acceptors reads in the bond section of the .psf file.      */
/*  This section contains a list of pairs of numbers where each pair is */
/*  represents two atoms that are part of an H-bond.  Each atom pair is */
/*  read in.                                                            */
/*                  */
/*  Acceptor atoms are the heavy atoms to which hydrogens directly      */
/*  orient in a hydrogen bond interaction.  There will always be an     */
/*  acceptor atom for each acceptor pair.  The antecedent atom, to      */
/*  which the acceptor is bound, may not be given in the structure,     */
/*  however, in which case the second atom index in the pair will be    */
/*  given as 0 in the PSF (and stored as -1 in this program's internal  */
/*  storage).                                                           */
/************************************************************************/

void Molecule::read_acceptors(FILE *fd)
{
#ifdef MEM_OPT_VERSION
    return;
#else
  int d[2];               // temporary storage of atom index
  register int j;      // Loop counter
  int num_read=0;    // Number of bonds read so far
        int num_no_ante=0;      // number of pairs with no antecedent

  /*  Allocate the array to hold the bonds      */
  acceptors=new Bond[numAcceptors];

  if (acceptors == NULL)
  {
    NAMD_die("memory allocations failed in Molecule::read_acceptors");
  }

  /*  Loop through and read in all the acceptors      */
  while (num_read < numAcceptors)
  {
    /*  Loop and read in the two atom indexes    */
    for (j=0; j<2; j++)
    {
      /*  Read the atom number from the file.         */
      /*  Subtract 1 to convert the index from the    */
      /*  1 to NumAtoms used in the file to the       */
      /*  0 to NumAtoms-1 that we need    */
      d[j]=NAMD_read_int(fd, "ACCEPTORS")-1;

      /*  Check to make sure the index isn't too big  */
      if (d[j] >= numAtoms)
      {
        char err_msg[128];

        sprintf(err_msg, "ACCEPTOR INDEX %d GREATER THAN NATOM %d IN DONOR # %d IN PSF FILE", d[j]+1, numAtoms, num_read+1);
        NAMD_die(err_msg);
      }

      /*  Check if there is an antecedent given */
      if (d[j] < 0)
                          num_no_ante++;
    }

    /*  Assign the atom indexes to the array element  */
    Bond *b = &(acceptors[num_read]);
    b->atom1=d[0];
    b->atom2=d[1];

    num_read++;
  }

  return;
#endif
}
/*      END OF FUNCTION read_acceptors      */


/************************************************************************/
/*                  */
/*      FUNCTION read_exclusions      */
/*                  */
/*   INPUTS:                */
/*  fd - file descriptor for .psf file        */
/*                  */
/*  read_exclusions reads in the explicit non-bonded exclusions     */
/*  from the .psf file.  This section is a little funky, so hang on.    */
/*  Ok, first there is a list of atom indexes that is NumExclusions     */
/*  long.  These are in some sense the atoms that will be exlcuded.     */
/*  Following this list is a list of NumAtoms length that is a list     */
/*  of indexes into the list of excluded atoms.  So an example.  Suppose*/
/*  we have a 5 atom simulation with 3 explicit exclusions.  The .psf   */
/*  file could look like:            */
/*                  */
/*  3!NNB                */
/*  3 4 5                */
/*  0 1 3 3 3              */
/*                  */
/*  This would mean that atom 1 has no explicit exclusions.  Atom 2     */
/*  has an explicit exclusion with atom 3.  Atom 3 has an explicit      */
/*  exclusion with atoms 4 AND 5.  And atoms 4 and 5 have no explicit   */
/*  exclusions.  Got it!?!  I'm not sure who dreamed this up . . .      */
/*                  */
/************************************************************************/

void Molecule::read_exclusions(FILE *fd)

{
#ifdef MEM_OPT_VERSION
    return;
#else
  int *exclusion_atoms;  //  Array of indexes of excluded atoms
  register int num_read=0;    //  Number fo exclusions read in
  int current_index;  //  Current index value
  int last_index;    //  the previous index value
  register int insert_index=0;  //  index of where we are in exlcusions array

  /*  Allocate the array of exclusion structures and the array of */
  /*  exlcuded atom indexes          */
  exclusions      = new Exclusion[numExclusions];
  exclusion_atoms = new int[numExclusions];

  if ( (exclusions == NULL) || (exclusion_atoms == NULL) )
  {
    NAMD_die("memory allocation failed in Molecule::read_exclusions");
  }

  /*  First, read in the excluded atoms list      */
  for (num_read=0; num_read<numExclusions; num_read++)
  {
    /*  Read the atom number from the file. Subtract 1 to   */
    /*  convert the index from the 1 to NumAtoms used in the*/
    /*  file to the  0 to NumAtoms-1 that we need    */
    exclusion_atoms[num_read]=NAMD_read_int(fd, "IMPROPERS")-1;

    /*  Check for an illegal index        */
    if (exclusion_atoms[num_read] >= numAtoms)
    {
      char err_msg[128];

      sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN PSF FILE", exclusion_atoms[num_read]+1, numAtoms, num_read+1);
      NAMD_die(err_msg);
    }
  }

  /*  Now, go through and read the list of NumAtoms pointers into */
  /*  the array that we just read in        */
  last_index=0;

  for (num_read=0; num_read<numAtoms; num_read++)
  {
    /*  Read in the current index value      */
    current_index=NAMD_read_int(fd, "EXCLUSIONS");

    /*  Check for an illegal pointer      */
    if (current_index>numExclusions)
    {
      char err_msg[128];

      sprintf(err_msg, "EXCLUSION INDEX %d LARGER THAN NUMBER OF EXLCUSIONS %d IN PSF FILE, EXCLUSION #%d\n", 
         current_index+1, numExclusions, num_read);
      NAMD_die(err_msg);
    }

    /*  Check to see if it matches the last index.  If so   */
    /*  than this atom has no exclusions.  If not, then     */
    /*  we have to build some exclusions      */
    if (current_index != last_index)
    {
      /*  This atom has some exlcusions.  Loop from   */
      /*  the last_index to the current index.  This  */
      /*  will include how ever many exclusions this  */
      /*  atom has          */
      for (insert_index=last_index; 
           insert_index<current_index; insert_index++)
      {
        /*  Assign the two atoms involved.      */
        /*  The first one is our position in    */
        /*  the list, the second is based on    */
        /*  the pointer into the index list     */
        int a1 = num_read;
        int a2 = exclusion_atoms[insert_index];
        if ( a1 < a2 ) {
          exclusions[insert_index].atom1 = a1;
          exclusions[insert_index].atom2 = a2;
        } else if ( a2 < a1 ) {
          exclusions[insert_index].atom1 = a2;
          exclusions[insert_index].atom2 = a1;
        } else {
          char err_msg[128];
          sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN PSF FILE\n", a1+1);
          NAMD_die(err_msg);
        }
      }

      last_index=current_index;
    }
  }

  /*  Free our temporary list of indexes        */
  delete [] exclusion_atoms;

  return;
#endif
}
/*      END OF FUNCTION read_exclusions      */

/************************************************************************/
/*                  */
/*        FUNCTION read_lphosts    */
/*                  */
/*   INPUTS:                */
/*  fd - file pointer to the .psf file        */
/*                  */
/*  this function reads in the lone pair host section of the .psf file. */
/*                  */
void Molecule::read_lphosts(FILE *fd)
{
  char buffer[512];  // Buffer for reading from file
  char lptype[8];
  int numhosts, index, i, read_count;
  Real distance, angle, dihedral;

  lphosts = new Lphost[numLphosts];
  if (lphosts == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_lphosts");
  }
  for (i = 0;  i < numLphosts;  i++)
  {
    NAMD_read_line(fd, buffer);
    if ( (NAMD_blank_string(buffer)) || (buffer[0] == '!') ) continue;
    read_count=sscanf(buffer, "%d %d %6s %f %f %f",
        &numhosts, &index, lptype, &distance, &angle, &dihedral);
    if (read_count != 6 || numhosts != 3 || index != 4*i + 1
        || strcmp(lptype,"F") != 0)
    {
      char err_msg[128];
      sprintf(err_msg, "BAD FORMAT FOR LPHOST LINE %d IN PSF FILE LINE\n"
          "LINE=%s\n", i+1, buffer);
      NAMD_die(err_msg);
    }
    lphosts[i].distance = distance;
    lphosts[i].angle = angle * (M_PI/180);        // convert to radians
    lphosts[i].dihedral = dihedral * (M_PI/180);  // convert to radians
  }
  for (i = 0;  i < numLphosts;  i++) {
    lphosts[i].atom1 = NAMD_read_int(fd, "LPHOSTS")-1;
    lphosts[i].atom2 = NAMD_read_int(fd, "LPHOSTS")-1;
    lphosts[i].atom3 = NAMD_read_int(fd, "LPHOSTS")-1;
    lphosts[i].atom4 = NAMD_read_int(fd, "LPHOSTS")-1;
  }
}
/*      END OF FUNCTION read_lphosts    */

/************************************************************************/
/*                  */
/*        FUNCTION read_anisos     */
/*                  */
/*   INPUTS:                */
/*  fd - file pointer to the .psf file        */
/*                  */
/*  this function reads in the anisotropic terms section of .psf file. */
/*                  */
void Molecule::read_anisos(FILE *fd)
{
  char buffer[512];  // Buffer for reading from file
  int numhosts, index, i, read_count;
  Real k11, k22, k33;

  anisos = new Aniso[numAnisos];
  if (anisos == NULL)
  {
    NAMD_die("memory allocation failed in Molecule::read_anisos");
  }
  for (i = 0;  i < numAnisos;  i++)
  {
    NAMD_read_line(fd, buffer);
    if ( (NAMD_blank_string(buffer)) || (buffer[0] == '!') ) continue;
    read_count=sscanf(buffer, "%f %f %f", &k11, &k22, &k33);
    if (read_count != 3)
    {
      char err_msg[128];
      sprintf(err_msg, "BAD FORMAT FOR ANISO LINE %d IN PSF FILE LINE\n"
          "LINE=%s\n", i+1, buffer);
      NAMD_die(err_msg);
    }
    anisos[i].k11 = k11;
    anisos[i].k22 = k22;
    anisos[i].k33 = k33;
  }
  for (i = 0;  i < numAnisos;  i++) {
    anisos[i].atom1 = NAMD_read_int(fd, "ANISOS")-1;
    anisos[i].atom2 = NAMD_read_int(fd, "ANISOS")-1;
    anisos[i].atom3 = NAMD_read_int(fd, "ANISOS")-1;
    anisos[i].atom4 = NAMD_read_int(fd, "ANISOS")-1;
  }
}
/*      END OF FUNCTION read_anisos     */

void Molecule::setOccupancyData(molfile_atom_t *atomarray){
    occupancy = new float[numAtoms];
    for(int i=0; i<numAtoms; i++) {
        occupancy[i] = atomarray[i].occupancy;
    }
}

void Molecule::setBFactorData(molfile_atom_t *atomarray){
    bfactor = new float[numAtoms];
    for(int i=0; i<numAtoms; i++) {
        bfactor[i] = atomarray[i].bfactor;
    }
}

void Molecule::plgLoadAtomBasics(molfile_atom_t *atomarray){
#ifndef MEM_OPT_VERSION
    atoms = new Atom[numAtoms];
    atomNames = new AtomNameInfo[numAtoms];
    if(simParams->genCompressedPsf) {
        atomSegResids = new AtomSegResInfo[numAtoms];
    }    
    hydrogenGroup.resize(0);

    ResidueLookupElem *tmpResLookup = resLookup;

    for(int i=0; i<numAtoms; i++) {
        int reslength = strlen(atomarray[i].resname)+1;
        int namelength = strlen(atomarray[i].name)+1;
        int typelength = strlen(atomarray[i].type)+1;
        atomNames[i].resname = nameArena->getNewArray(reslength);
        atomNames[i].atomname = nameArena->getNewArray(namelength);
        atomNames[i].atomtype = nameArena->getNewArray(typelength);
        strcpy(atomNames[i].resname, atomarray[i].resname);
        strcpy(atomNames[i].atomname, atomarray[i].name);
        strcpy(atomNames[i].atomtype, atomarray[i].type);

        atoms[i].mass = atomarray[i].mass;
        atoms[i].charge = atomarray[i].charge;
        atoms[i].status = UnknownAtom;

        //add this atom to residue lookup table
        if(tmpResLookup) {
            tmpResLookup = tmpResLookup->append(atomarray[i].segid, atomarray[i].resid, i);
        }

        if(atomSegResids) { //for compressing molecule information
            AtomSegResInfo *one = atomSegResids + i;
            memcpy(one->segname, atomarray[i].segid, strlen(atomarray[i].segid)+1);
            one->resid = atomarray[i].resid;
        }
        //Determine the type of the atom
        if(atoms[i].mass <= 0.05) {
            atoms[i].status |= LonepairAtom;
        }else if(atoms[i].mass < 1.0) {
            atoms[i].status |= DrudeAtom;
        }else if(atoms[i].mass <= 3.5) {
            atoms[i].status |= HydrogenAtom;
        }else if((atomNames[i].atomname[0] == 'O') &&
                 (atoms[i].mass>=14.0) && (atoms[i].mass<=18.0)){
            atoms[i].status |= OxygenAtom;
        }
        //Lookk up the vdw constants for this atom
        params->assign_vdw_index(atomNames[i].atomtype, &atoms[i]);
    }
#endif
}

void Molecule::plgLoadBonds(int *from, int *to){
#ifndef MEM_OPT_VERSION
    bonds = new Bond[numBonds];
    char atom1name[11];
    char atom2name[11];
    int realNumBonds = 0;
    for(int i=0; i<numBonds; i++) {
        Bond *thisBond = bonds+realNumBonds;
        thisBond->atom1 = from[i]-1;
        thisBond->atom2 = to[i]-1;
        /* Get the atom type for the two atoms.
         * When we query the parameter object, we
         * need to send the atom type that is alphabetically
         * first as atom 1.
         */
        if(strcasecmp(atomNames[thisBond->atom1].atomtype,
                      atomNames[thisBond->atom2].atomtype)<0) {
            strcpy(atom1name, atomNames[thisBond->atom1].atomtype);
            strcpy(atom2name, atomNames[thisBond->atom2].atomtype);
        }else{
            strcpy(atom2name, atomNames[thisBond->atom1].atomtype);
            strcpy(atom1name, atomNames[thisBond->atom2].atomtype);
        }
        params->assign_bond_index(atom1name, atom2name, thisBond);

        //Make sure this isn't a fake bond meant for shake in x-plor
        Real k, x0;
        params->get_bond_params(&k, &x0, thisBond->bond_type);
        if(simParams->watmodel == WAT_SWM4) {
            //need to retain Lonepair bonds for Drude
            if(k!=0. || is_lp(thisBond->atom1) || 
               is_lp(thisBond->atom2)) {               
                realNumBonds++;
	    }
        }else{
            if(k != 0.) realNumBonds++;
        }
    }

    if(numBonds != realNumBonds) {
        iout << iWARN << "Ignored" << numBonds-realNumBonds <<
            "bonds with zero force constants.\n" <<endi;
        iout << iWARN << "Will get H-H distance in rigid H20 from H-O-H angle.\n" <<endi;
    }
    numBonds = realNumBonds;
#endif
}

void Molecule::plgLoadAngles(int *plgAngles)
{    
#ifndef MEM_OPT_VERSION
    char atom1name[11];
    char atom2name[11];
    char atom3name[11];

    angles=new Angle[numAngles];
    int *atomid = plgAngles;
    int numRealAngles = 0;
    for(int i=0; i<numAngles; i++) {
        Angle *thisAngle = angles+numRealAngles;
        thisAngle->atom1 = atomid[0]-1;
        thisAngle->atom2 = atomid[1]-1;
        thisAngle->atom3 = atomid[2]-1;
        atomid += 3;

        if(strcasecmp(atomNames[thisAngle->atom1].atomtype,
                      atomNames[thisAngle->atom2].atomtype)<0) {
            strcpy(atom1name, atomNames[thisAngle->atom1].atomtype);
            strcpy(atom2name, atomNames[thisAngle->atom2].atomtype);
            strcpy(atom3name, atomNames[thisAngle->atom3].atomtype);
        }else{
            strcpy(atom1name, atomNames[thisAngle->atom3].atomtype);
            strcpy(atom2name, atomNames[thisAngle->atom2].atomtype);
            strcpy(atom3name, atomNames[thisAngle->atom1].atomtype);
        }

        params->assign_angle_index(atom1name, atom2name, atom3name, thisAngle);

        Real k, t0, k_ub, r_ub;
        params->get_angle_params(&k, &t0, &k_ub, &r_ub, thisAngle->angle_type);
        if(k!=0. || k_ub!=0.) numRealAngles++;
    }

    if(numAngles != numRealAngles) {
        iout << iWARN << "Ignored" << numAngles-numRealAngles << 
            " angles with zero force constants.\n" << endi; 
    }
    numAngles = numRealAngles;
#endif
}

void Molecule::plgLoadDihedrals(int *plgDihedrals)
{
#ifndef MEM_OPT_VERSION
    char atom1name[11];
    char atom2name[11];
    char atom3name[11];
    char atom4name[11];
    int lastAtomIds[4];
    int multiplicity = 1; //multiplicity of the current bond

    lastAtomIds[0]=lastAtomIds[1]=lastAtomIds[2]=lastAtomIds[3]=-1;
    dihedrals = new Dihedral[numDihedrals];
    int numRealDihedrals = 0;
    int *atomid = plgDihedrals;
    for(int i=0; i<numDihedrals; i++, atomid+=4) {
        Dihedral *thisDihedral = dihedrals + numRealDihedrals;
        Bool duplicate_bond = TRUE;
        for(int j=0; j<4; j++) {
            if(atomid[j] != lastAtomIds[j]) {
                duplicate_bond = FALSE;
            }
            lastAtomIds[j] = atomid[j];
        }

        strcpy(atom1name, atomNames[atomid[0]-1].atomtype);
        strcpy(atom2name, atomNames[atomid[1]-1].atomtype);
        strcpy(atom3name, atomNames[atomid[2]-1].atomtype);
        strcpy(atom4name, atomNames[atomid[3]-1].atomtype);

        if(duplicate_bond) {
            multiplicity++;
            if(multiplicity==2) {
                numMultipleDihedrals++;
            }
        }else{
            multiplicity=1;
            numRealDihedrals++;
        }

        params->assign_dihedral_index(atom1name, atom2name,
                                      atom3name, atom4name, thisDihedral,
                                      multiplicity);
        thisDihedral->atom1 = atomid[0]-1;
        thisDihedral->atom2 = atomid[1]-1;
        thisDihedral->atom3 = atomid[2]-1;
        thisDihedral->atom4 = atomid[3]-1;
    }

    numDihedrals = numRealDihedrals;
#endif
}

void Molecule::plgLoadImpropers(int *plgImpropers)
{
#ifndef MEM_OPT_VERSION
    char atom1name[11];
    char atom2name[11];
    char atom3name[11];
    char atom4name[11];
    int lastAtomIds[4];
    int multiplicity = 1; //multiplicity of the current bond

    lastAtomIds[0]=lastAtomIds[1]=lastAtomIds[2]=lastAtomIds[3]=-1;
    impropers = new Improper[numImpropers];
    int numRealImpropers = 0;
    int *atomid = plgImpropers;
    for(int i=0; i<numImpropers; i++, atomid+=4) {
        Improper *thisImproper = impropers + numRealImpropers;
        Bool duplicate_bond = TRUE;
        for(int j=0; j<4; j++) {
            if(atomid[j] != lastAtomIds[j]) {
                duplicate_bond = FALSE;
            }
            lastAtomIds[j] = atomid[j];
        }

        strcpy(atom1name, atomNames[atomid[0]-1].atomtype);
        strcpy(atom2name, atomNames[atomid[1]-1].atomtype);
        strcpy(atom3name, atomNames[atomid[2]-1].atomtype);
        strcpy(atom4name, atomNames[atomid[3]-1].atomtype);

        if(duplicate_bond) {
            multiplicity++;
            if(multiplicity==2) {
                numMultipleImpropers++;
            }
        }else{
            multiplicity=1;
            numRealImpropers++;
        }

        params->assign_improper_index(atom1name, atom2name,
                                      atom3name, atom4name, thisImproper,
                                      multiplicity);
        thisImproper->atom1 = atomid[0]-1;
        thisImproper->atom2 = atomid[1]-1;
        thisImproper->atom3 = atomid[2]-1;
        thisImproper->atom4 = atomid[3]-1;
    }

    numImpropers = numRealImpropers;
#endif
}

void Molecule::plgLoadCrossterms(int *plgCterms)
{
#ifndef MEM_OPT_VERSION
    char atom1name[11];
    char atom2name[11];
    char atom3name[11];
    char atom4name[11];
    char atom5name[11];
    char atom6name[11];
    char atom7name[11];
    char atom8name[11];
    int lastAtomIds[8];    

    for(int i=0; i<8; i++)
        lastAtomIds[i]=-1;
    
    crossterms = new Crossterm[numCrossterms];
    int numRealCrossterms = 0;
    int *atomid = plgCterms;
    for(int i=0; i<numCrossterms; i++, atomid+=8) {
        Crossterm *thisCrossterm = crossterms + numRealCrossterms;
        Bool duplicate_bond = TRUE;
        for(int j=0; j<8; j++) {
            if(atomid[j] != lastAtomIds[j]) {
                duplicate_bond = FALSE;
            }
            lastAtomIds[j] = atomid[j];
        }

        strcpy(atom1name, atomNames[atomid[0]-1].atomtype);
        strcpy(atom2name, atomNames[atomid[1]-1].atomtype);
        strcpy(atom3name, atomNames[atomid[2]-1].atomtype);
        strcpy(atom4name, atomNames[atomid[3]-1].atomtype);
        strcpy(atom5name, atomNames[atomid[4]-1].atomtype);
        strcpy(atom6name, atomNames[atomid[5]-1].atomtype);
        strcpy(atom7name, atomNames[atomid[6]-1].atomtype);
        strcpy(atom8name, atomNames[atomid[7]-1].atomtype);

        if(duplicate_bond) {
            iout << iWARN <<"Duplicate cross-term detected.\n" << endi;
        }else
            numRealCrossterms++;

        params->assign_crossterm_index(atom1name, atom2name,
                                       atom3name, atom4name, atom5name,
                                       atom6name, atom7name, atom8name,
                                       thisCrossterm);

        thisCrossterm->atom1 = atomid[0]-1;
        thisCrossterm->atom2 = atomid[1]-1;
        thisCrossterm->atom3 = atomid[2]-1;
        thisCrossterm->atom4 = atomid[3]-1;
        thisCrossterm->atom5 = atomid[4]-1;
        thisCrossterm->atom6 = atomid[5]-1;
        thisCrossterm->atom7 = atomid[6]-1;
        thisCrossterm->atom8 = atomid[7]-1;
    }

    numCrossterms = numRealCrossterms;
#endif
}


/************************************************************************/
/*                  */
/*      FUNCTION print_atoms        */
/*                  */
/*  print_atoms prints out the list of atoms stored in this object. */
/*  It is inteded mainly for debugging purposes.      */
/*                  */
/************************************************************************/

void Molecule::print_atoms(Parameters *params)

{
  register int i;
  Real sigma;
  Real epsilon;
  Real sigma14;
  Real epsilon14;

  DebugM(2,"ATOM LIST\n" \
      << "******************************************\n" \
                  << "NUM  NAME TYPE RES  MASS    CHARGE CHARGE   FEP-CHARGE"  \
      << "SIGMA   EPSILON SIGMA14 EPSILON14\n" \
        << endi);

  for (i=0; i<numAtoms; i++)
  {
    params->get_vdw_params(&sigma, &epsilon, &sigma14, &epsilon14, 
        atoms[i].vdw_type);

#ifdef MEM_OPT_VERSION
     DebugM(2,i+1 << " " << atomNamePool[atomNames[i].atomnameIdx]  \
              << " " << atomTypePool[atomNames[i].atomtypeIdx] << " " \
              << resNamePool[atomNames[i].resnameIdx]  << " "  \
	    << atommass(i)  \
        << " " << atomcharge(i) << " " << sigma \
        << " " << epsilon << " " << sigma14 \
        << " " << epsilon14 << "\n" \
        << endi);
#else
    DebugM(2,i+1 << " " << atomNames[i].atomname  \
              << " " << atomNames[i].atomtype << " " \
              << atomNames[i].resname  << " " << atoms[i].mass  \
        << " " << atoms[i].charge << " " << sigma \
        << " " << epsilon << " " << sigma14 \
        << " " << epsilon14 << "\n" \
        << endi);
#endif
  }
}
/*      END OF FUNCTION print_atoms      */

/************************************************************************/
/*                  */
/*      FUNCTION print_bonds        */
/*                  */
/*  print_bonds prints out the list of bonds stored in this object. */
/*  It is inteded mainly for debugging purposes.      */
/*                  */
/************************************************************************/

void Molecule::print_bonds(Parameters *params)

{
  register int i;
  Real k;
  Real x0;

  DebugM(2,"BOND LIST\n" << "********************************\n" \
      << "ATOM1 ATOM2 TYPE1 TYPE2      k        x0" \
      << endi);

  #ifdef MEM_OPT_VERSION
  for(i=0; i<numAtoms; i++){
      AtomSignature *sig = &atomSigPool[eachAtomSig[i]];
      int bCnt = sig->bondCnt;
      TupleSignature *bSigs = sig->bondSigs;
      for(int j=0; j<bCnt; j++){
          Bond aBond;
          aBond.atom1 = i;
          aBond.atom2 = i+bSigs[j].offset[0];
          aBond.bond_type = bSigs[j].tupleParamType;
          params->get_bond_params(&k, &x0, aBond.bond_type);

          DebugM(2,aBond.atom1+1 << " " \
             << aBond.atom2+1 << " "   \
             << get_atomtype(aBond.atom1) << " "  \
             << get_atomtype(aBond.atom2) << " " << k \
             << " " << x0 << endi);
      }
  }
  #else
  for (i=0; i<numBonds; i++)
  {
    params->get_bond_params(&k, &x0, bonds[i].bond_type);

    DebugM(2,bonds[i].atom1+1 << " " \
       << bonds[i].atom2+1 << " "   \
       << atomNames[bonds[i].atom1].atomtype << " "  \
       << atomNames[bonds[i].atom2].atomtype << " " << k \
       << " " << x0 << endi);
  }
  #endif
}
/*      END OF FUNCTION print_bonds      */

/************************************************************************/
/*                  */
/*      FUNCTION print_exclusions      */
/*                  */
/*  print_exlcusions prints out the list of exlcusions stored in    */
/*  this object.  It is inteded mainly for debugging purposes.    */
/*                  */
/************************************************************************/

void Molecule::print_exclusions()
{
#ifdef MEM_OPT_VERSION
    DebugM(2, "WARNING: this function is not availabe in memory optimized version!\n" << endi);
#else
  register int i;

  DebugM(2,"EXPLICIT EXCLUSION LIST\n" \
      << "********************************\n" \
            << "ATOM1 ATOM2 " \
      << endi);

  for (i=0; i<numExclusions; i++)
  {
    DebugM(2,exclusions[i].atom1+1 << "  " \
       << exclusions[i].atom2+1 << endi);
  }
#endif
}
/*      END OF FUNCTION print_exclusions    */

/************************************************************************/
/*                  */
/*      FUNCTION send_Molecule        */
/*                  */
/*  send_Molecule is used by the Master node to distribute the      */
/*   structural information to all the client nodes.  It is NEVER called*/
/*   by the client nodes.              */
/*                  */
/************************************************************************/


void Molecule::send_Molecule(MOStream *msg)

{
    #ifdef MEM_OPT_VERSION
    //  Now build arrays of indexes into these arrays by atom      
    //generating all the new atom signatures and exclusion signatures only on pe0
    build_lists_by_atom();
    #endif
      
  /*//  Message to send to clients
  int bufSize = BUFSIZE;
  // When the simulation system is very large, then the buffer size should be expanded to reduce the number of one-to-all broadcasts.
  if(numAtoms>=1000000) bufSize=16*BUFSIZE;
  MOStream *msg=com_obj->newOutputStream(ALLBUTME, MOLECULETAG, bufSize);
  if ( msg == NULL )
  {
    NAMD_die("Memory allocation failed in Molecule::send_Molecule");
  }*/

  #ifdef MEM_OPT_VERSION
      msg->put(numAtoms);
      //mass and charge pool needed to be sent to other processors
      //for the sake of function call: build_atom_status
      msg->put(massPoolSize);
      msg->put(massPoolSize, atomMassPool);
      msg->put(numAtoms*sizeof(Index), (char *)eachAtomMass);

      msg->put(chargePoolSize);
      msg->put(chargePoolSize, atomChargePool);
      msg->put(numAtoms*sizeof(Index), (char *)eachAtomCharge);

      //vdw_type, partner etc.
      msg->put(numAtoms*sizeof(AtomCstInfo), (char *)atoms);

      //put atoms' signatures
      msg->put(atomSigPoolSize);
      for(int i=0; i<atomSigPoolSize; i++)
          atomSigPool[i].pack(msg);

      //put atom's exclusion signatures
      msg->put(exclSigPoolSize);
      for(int i=0; i<exclSigPoolSize; i++)
          exclSigPool[i].pack(msg);

      //put eachAtomSig and eachAtomExclSig
      msg->put(numAtoms*sizeof(Index), (char *)eachAtomSig);
      msg->put(numAtoms*sizeof(Index), (char *)eachAtomExclSig);

      msg->put(numAtoms*sizeof(int32), (char *)clusterSigs);
  #else
      msg->put(numAtoms);
      msg->put(numAtoms*sizeof(Atom), (char*)atoms);
  #endif
      //  Send the bond information
      msg->put(numRealBonds);
      msg->put(numBonds);

      #ifndef MEM_OPT_VERSION
      if (numBonds)
      {
        msg->put(numBonds*sizeof(Bond), (char*)bonds);
      }
      #endif

      //  Send the angle information
      msg->put(numAngles);

      #ifndef MEM_OPT_VERSION
      if (numAngles)
      {
        msg->put(numAngles*sizeof(Angle), (char*)angles);
      }
      #endif

      //  Send the dihedral information
      msg->put(numDihedrals);

      #ifndef MEM_OPT_VERSION
      if (numDihedrals)
      {
        msg->put(numDihedrals*sizeof(Dihedral), (char*)dihedrals);
      }
      #endif

      //  Send the improper information
      msg->put(numImpropers);

      #ifndef MEM_OPT_VERSION
      if (numImpropers)
      {
        msg->put(numImpropers*sizeof(Improper), (char*)impropers);
      }
      #endif

      //  Send the crossterm information
      msg->put(numCrossterms);

      #ifndef MEM_OPT_VERSION
      if (numCrossterms)
      {
        msg->put(numCrossterms*sizeof(Crossterm), (char*)crossterms);
      }
      #endif

      // send the hydrogen bond donor information
      msg->put(numDonors);

      if(numDonors)
      {
        msg->put(numDonors*sizeof(Bond), (char*)donors);
      }

      // send the hydrogen bond acceptor information
      msg->put(numAcceptors);

      if(numAcceptors)
      {
        msg->put(numAcceptors*sizeof(Bond), (char*)acceptors);
      }

      //  Send the exclusion information
      #ifndef MEM_OPT_VERSION
      msg->put(numExclusions);

      if (numExclusions)
      {
        msg->put(numExclusions*sizeof(Exclusion), (char*)exclusions);
      }
      #endif

      //hydrogen group info is calculated when generating the compressed molecule
      //information, so this has to be distributed to other nodes instead of
      //being recalculated in build_atom_status in receive_Molecule function
      #ifdef MEM_OPT_VERSION
      msg->put(numHydrogenGroups);      
      msg->put(numAtoms*sizeof(HydrogenGroupID), (char *)hydrogenGroup.begin());      
      #endif
      
      //  Send the constraint information, if used
      if (simParams->constraintsOn)
      {
         msg->put(numConstraints);
         
         msg->put(numAtoms, consIndexes);
         
         if (numConstraints)
         {
           msg->put(numConstraints*sizeof(ConstraintParams), (char*)consParams);
         }
      }
      
      /* BEGIN gf */
      // Send the gridforce information, if used
      if (simParams->mgridforceOn)
      {
	 DebugM(3, "Sending gridforce info\n");
	 msg->put(numGridforceGrids);
	 
	 for (int grid = 0; grid < numGridforceGrids; grid++) {
	     msg->put(numGridforces[grid]);
	     msg->put(numAtoms, gridfrcIndexes[grid]);
	     if (numGridforces[grid])
	     {
		 msg->put(numGridforces[grid]*sizeof(GridforceParams), (char*)gridfrcParams[grid]);
	     }
	     gridfrcGrid[grid]->pack(msg);	// grid object writes its private data to message itself
	 }
      }
      /* END gf */

      //  Send the stirring information, if used
      if (simParams->stirOn)
      {
	       //CkPrintf ("DEBUG: putting numStirredAtoms..\n");
         msg->put(numStirredAtoms);
	       //CkPrintf ("DEBUG: putting numAtoms,stirIndexes.. numAtoms=%d\n",numStirredAtoms);
         msg->put(numAtoms, stirIndexes);
	       //CkPrintf ("DEBUG: if numStirredAtoms..\n");
         if (numStirredAtoms)
         {
           //CkPrintf ("DEBUG: big put, with (char*)stirParams\n");
           msg->put(numStirredAtoms*sizeof(StirParams), (char*)stirParams);
         }
      }
      
      
      //  Send the moving drag information, if used
      if (simParams->movDragOn) {
         msg->put(numMovDrag);
         msg->put(numAtoms, movDragIndexes);
         if (numMovDrag)
         {
           msg->put(numMovDrag*sizeof(MovDragParams), (char*)movDragParams);
         }
      }
      
      //  Send the rotating drag information, if used
      if (simParams->rotDragOn) {
         msg->put(numRotDrag);
         msg->put(numAtoms, rotDragIndexes);
         if (numRotDrag)
         {
           msg->put(numRotDrag*sizeof(RotDragParams), (char*)rotDragParams);
         }
      }
      
      //  Send the "constant" torque information, if used
      if (simParams->consTorqueOn) {
         msg->put(numConsTorque);
         msg->put(numAtoms, consTorqueIndexes);
         if (numConsTorque)
         {
           msg->put(numConsTorque*sizeof(ConsTorqueParams), (char*)consTorqueParams);
         }
      }
      
      // Send the constant force information, if used
      if (simParams->consForceOn)
      { msg->put(numConsForce);
        msg->put(numAtoms, consForceIndexes);
        if (numConsForce)
          msg->put(numConsForce*sizeof(Vector), (char*)consForce);
      }

      //  Send the langevin parameters, if active
      if (simParams->langevinOn || simParams->tCoupleOn)
      {
        msg->put(numAtoms, langevinParams);
      }

      //  Send fixed atoms, if active
      if (simParams->fixedAtomsOn)
      {
        msg->put(numFixedAtoms);
        msg->put(numAtoms, fixedAtomFlags);
	msg->put(numFixedRigidBonds);
      }

      if (simParams->excludeFromPressure) {
        msg->put(numExPressureAtoms);
        msg->put(numAtoms, exPressureAtomFlags);
      }

//fepb
      // send fep atom info
      if (simParams->alchFepOn || simParams->alchThermIntOn || simParams->lesOn || simParams->pairInteractionOn) {
        msg->put(numFepInitial);
        msg->put(numFepFinal);
        msg->put(numAtoms*sizeof(char), (char*)fepAtomFlags);
      }
//fepe

      // DRUDE: send data read from PSF
      msg->put(is_drude_psf);
      if (is_drude_psf) {
        msg->put(numAtoms*sizeof(DrudeConst), (char*)drudeConsts);
        msg->put(numLphosts);
        msg->put(numLphosts*sizeof(Lphost), (char*)lphosts);
        msg->put(numAnisos);
        msg->put(numAnisos*sizeof(Aniso), (char*)anisos);
      }
      // DRUDE

      // Broadcast the message to the other nodes
      msg->end();
      delete msg;

      #ifdef MEM_OPT_VERSION
      build_excl_check_signatures();
      #else      
      //  Now build arrays of indexes into these arrays by atom      
      build_lists_by_atom();
      #endif

}
    /*      END OF FUNCTION send_Molecule      */

    /************************************************************************/
    /*                  */
    /*      FUNCTION receive_Molecule      */
    /*                  */
    /*  receive_Molecule is used by all the clients to receive the  */
    /*   structural data sent out by the master node.  It is NEVER called   */
    /*   by the Master node.            */
    /*                  */
    /************************************************************************/

void Molecule::receive_Molecule(MIStream *msg)
{
      //  Get the atom information
      msg->get(numAtoms);

  #ifdef MEM_OPT_VERSION   
      //are mass and charge pool needed to be sent to other processors???
      msg->get(massPoolSize);
      if(atomMassPool) delete [] atomMassPool;
      atomMassPool = new Real[massPoolSize];
      msg->get(massPoolSize, atomMassPool);
      if(eachAtomMass) delete [] eachAtomMass;
      eachAtomMass = new Index[numAtoms];
      msg->get(numAtoms*sizeof(Index), (char *)eachAtomMass);
      
      msg->get(chargePoolSize);
      if(atomChargePool) delete [] atomChargePool;
      atomChargePool = new Real[chargePoolSize];
      msg->get(chargePoolSize, atomChargePool);
      if(eachAtomCharge) delete [] eachAtomCharge;
      eachAtomCharge = new Index[numAtoms];
      msg->get(numAtoms*sizeof(Index), (char *)eachAtomCharge);

      //vdw_type, partner etc.
      if(atoms) delete [] atoms;
      atoms = new AtomCstInfo[numAtoms];
      msg->get(numAtoms*sizeof(AtomCstInfo), (char *)atoms);

      //get atoms' signatures
      msg->get(atomSigPoolSize);
      if(atomSigPool) delete [] atomSigPool;
      atomSigPool = new AtomSignature[atomSigPoolSize];
      for(int i=0; i<atomSigPoolSize; i++)
          atomSigPool[i].unpack(msg);

      //get exclusions' signatures
      msg->get(exclSigPoolSize);
      if(exclSigPool) delete [] exclSigPool;
      exclSigPool = new ExclusionSignature[exclSigPoolSize];
      for(int i=0; i<exclSigPoolSize; i++)
          exclSigPool[i].unpack(msg);

      //get eachAtomSig and eachAtomExclSig
      if(eachAtomSig) delete [] eachAtomSig;
      eachAtomSig = new Index[numAtoms];
      msg->get(numAtoms*sizeof(Index), (char *)eachAtomSig);
      if(eachAtomExclSig) delete [] eachAtomExclSig;
      eachAtomExclSig = new Index[numAtoms];
      msg->get(numAtoms*sizeof(Index), (char *)eachAtomExclSig);

      if(clusterSigs) delete [] clusterSigs;
      clusterSigs = new int32[numAtoms];
      msg->get(numAtoms*sizeof(int32), (char *)clusterSigs); 
  #else
      delete [] atoms;
      atoms= new Atom[numAtoms];

      msg->get(numAtoms*sizeof(Atom), (char*)atoms);
  #endif

      //  Get the bond information
      msg->get(numRealBonds);
      msg->get(numBonds);

      #ifndef MEM_OPT_VERSION
      if (numBonds)
      {
        delete [] bonds;
        bonds=new Bond[numBonds];

        msg->get(numBonds*sizeof(Bond), (char*)bonds);
      }
      #endif

      //  Get the angle information
      msg->get(numAngles);

      #ifndef MEM_OPT_VERSION
      if (numAngles)
      {
        delete [] angles;
        angles=new Angle[numAngles];

        msg->get(numAngles*sizeof(Angle), (char*)angles);
      }
      #endif

      //  Get the dihedral information
      msg->get(numDihedrals);

      #ifndef MEM_OPT_VERSION
      if (numDihedrals)
      {
        delete [] dihedrals;
        dihedrals=new Dihedral[numDihedrals];

        msg->get(numDihedrals*sizeof(Dihedral), (char*)dihedrals);
      }
      #endif

      //  Get the improper information
      msg->get(numImpropers);

      #ifndef MEM_OPT_VERSION
      if (numImpropers)
      {
        delete [] impropers;
        impropers=new Improper[numImpropers];

        msg->get(numImpropers*sizeof(Improper), (char*)impropers);
      }
      #endif

      //  Get the crossterm information
      msg->get(numCrossterms);

      #ifndef MEM_OPT_VERSION
      if (numCrossterms)
      {
        delete [] crossterms;
        crossterms=new Crossterm[numCrossterms];

        msg->get(numCrossterms*sizeof(Crossterm), (char*)crossterms);
      }
      #endif

      //  Get the hydrogen bond donors
      msg->get(numDonors);

      if (numDonors)
      {
        delete [] donors;
        donors=new Bond[numDonors];

        msg->get(numDonors*sizeof(Bond), (char*)donors);
      }

      //  Get the hydrogen bond acceptors
      msg->get(numAcceptors);

      if (numAcceptors)
      {
        delete [] acceptors;
        acceptors=new Bond[numAcceptors];

        msg->get(numAcceptors*sizeof(Bond), (char*)acceptors);
      }

      //  Get the exclusion information
      #ifndef MEM_OPT_VERSION
      msg->get(numExclusions);

      if (numExclusions)
      {
        delete [] exclusions;
        exclusions=new Exclusion[numExclusions];

        msg->get(numExclusions*sizeof(Exclusion), (char*)exclusions);
      }
      #endif

      #ifdef MEM_OPT_VERSION
      msg->get(numHydrogenGroups);      
      hydrogenGroup.resize(numAtoms);
      msg->get(numAtoms*sizeof(HydrogenGroupID), (char *)hydrogenGroup.begin());     
      #endif
      
      //  Get the constraint information, if they are active
      if (simParams->constraintsOn)
      {
         msg->get(numConstraints);

         delete [] consIndexes;
         consIndexes = new int32[numAtoms];
         
         msg->get(numAtoms, consIndexes);
         
         if (numConstraints)
         {
           delete [] consParams;
           consParams = new ConstraintParams[numConstraints];
      
           msg->get(numConstraints*sizeof(ConstraintParams), (char*)consParams);
         }
      }

      /* BEGIN gf */
      if (simParams->mgridforceOn)
      {
	 DebugM(3, "Receiving gridforce info\n");
	 
	 msg->get(numGridforceGrids);
	 
	 delete [] numGridforces;
	 numGridforces = new int[numGridforceGrids];
	 
	 delete [] gridfrcIndexes;	// Should I be deleting elements of these first?
	 delete [] gridfrcParams;
	 delete [] gridfrcGrid;
	 gridfrcIndexes = new int32*[numGridforceGrids];
	 gridfrcParams = new GridforceParams*[numGridforceGrids];
	 gridfrcGrid = new GridforceGrid*[numGridforceGrids];
	 
	 for (int grid = 0; grid < numGridforceGrids; grid++) {
	     msg->get(numGridforces[grid]);
	     
	     gridfrcIndexes[grid] = new int32[numAtoms];
	     msg->get(numAtoms, gridfrcIndexes[grid]);
	 
	     if (numGridforces[grid])
	     {
		 gridfrcParams[grid] = new GridforceParams[numGridforces[grid]];
		 msg->get(numGridforces[grid]*sizeof(GridforceParams), (char*)gridfrcParams[grid]);
	     }
	     
	     gridfrcGrid[grid] = new GridforceGrid(grid);
	     gridfrcGrid[grid]->unpack(msg);
             CProxy_ComputeGridForceNodeMgr 
               mgr(CkpvAccess(BOCclass_group).computeGridForceNodeMgr);
             GridDepositMsg *outmsg = new GridDepositMsg;
             outmsg->gridnum = grid;
             outmsg->grid = gridfrcGrid[grid];
             outmsg->num_grids = numGridforceGrids;
             mgr[CkMyNode()].depositInitialGrid(outmsg);
	 }
      }
      /* END gf */
      
      //  Get the stirring information, if stirring is  active
      if (simParams->stirOn)
      {
         msg->get(numStirredAtoms);

         delete [] stirIndexes;
         stirIndexes = new int32[numAtoms];
         
         msg->get(numAtoms, stirIndexes);
         
         if (numStirredAtoms)
         {
           delete [] stirParams;
           stirParams = new StirParams[numStirredAtoms];
      
           msg->get(numStirredAtoms*sizeof(StirParams), (char*)stirParams);
         }
      }
      
      //  Get the moving drag information, if it is active
      if (simParams->movDragOn) {
         msg->get(numMovDrag);
         delete [] movDragIndexes;
         movDragIndexes = new int32[numAtoms];
         msg->get(numAtoms, movDragIndexes);
         if (numMovDrag)
         {
           delete [] movDragParams;
           movDragParams = new MovDragParams[numMovDrag];
           msg->get(numMovDrag*sizeof(MovDragParams), (char*)movDragParams);
         }
      }
      
      //  Get the rotating drag information, if it is active
      if (simParams->rotDragOn) {
         msg->get(numRotDrag);
         delete [] rotDragIndexes;
         rotDragIndexes = new int32[numAtoms];
         msg->get(numAtoms, rotDragIndexes);
         if (numRotDrag)
         {
           delete [] rotDragParams;
           rotDragParams = new RotDragParams[numRotDrag];
           msg->get(numRotDrag*sizeof(RotDragParams), (char*)rotDragParams);
         }
      }
      
      //  Get the "constant" torque information, if it is active
      if (simParams->consTorqueOn) {
         msg->get(numConsTorque);
         delete [] consTorqueIndexes;
         consTorqueIndexes = new int32[numAtoms];
         msg->get(numAtoms, consTorqueIndexes);
         if (numConsTorque)
         {
           delete [] consTorqueParams;
           consTorqueParams = new ConsTorqueParams[numConsTorque];
           msg->get(numConsTorque*sizeof(ConsTorqueParams), (char*)consTorqueParams);
         }
      }
      
      // Get the constant force information, if it's active
      if (simParams->consForceOn)
      { msg->get(numConsForce);
        delete [] consForceIndexes;
        consForceIndexes = new int32[numAtoms];
        msg->get(numAtoms, consForceIndexes);
        if (numConsForce)
        { delete [] consForce;
          consForce = new Vector[numConsForce];
          msg->get(numConsForce*sizeof(Vector), (char*)consForce);
        }
      }

      //  Get the langevin parameters, if they are active
      if (simParams->langevinOn || simParams->tCoupleOn)
      {
        delete [] langevinParams;
        langevinParams = new Real[numAtoms];

        msg->get(numAtoms, langevinParams);
      }

      //  Get the fixed atoms, if they are active
      if (simParams->fixedAtomsOn)
      {
        delete [] fixedAtomFlags;
        fixedAtomFlags = new int32[numAtoms];

        msg->get(numFixedAtoms);
        msg->get(numAtoms, fixedAtomFlags);
        msg->get(numFixedRigidBonds);
      }

      if (simParams->excludeFromPressure) {
        exPressureAtomFlags = new int32[numAtoms];
        msg->get(numExPressureAtoms);
        msg->get(numAtoms, exPressureAtomFlags);
      }

//fepb
      //receive fep atom info
      if (simParams->alchFepOn || simParams->lesOn || simParams->alchThermIntOn || simParams->pairInteractionOn) {
        delete [] fepAtomFlags;
        fepAtomFlags = new unsigned char[numAtoms];

        msg->get(numFepInitial);
        msg->get(numFepFinal);
        msg->get(numAtoms*sizeof(unsigned char), (char*)fepAtomFlags);
      }
//fepe

      // DRUDE: receive data read from PSF
      msg->get(is_drude_psf);
      if (is_drude_psf) {
        delete[] drudeConsts;
        drudeConsts = new DrudeConst[numAtoms];
        msg->get(numAtoms*sizeof(DrudeConst), (char*)drudeConsts);
        msg->get(numLphosts);
        delete[] lphosts;
        lphosts = new Lphost[numLphosts];
        msg->get(numLphosts*sizeof(Lphost), (char*)lphosts);
        msg->get(numAnisos);
        delete[] anisos;
        anisos = new Aniso[numAnisos];
        msg->get(numAnisos*sizeof(Aniso), (char*)anisos);
      }
      // DRUDE

      //  Now free the message 
      delete msg;
      
      //  analyze the data and find the status of each atom
      build_atom_status();

      //  Now build arrays of indexes into these arrays by atom
      #ifdef MEM_OPT_VERSION
      build_excl_check_signatures();
      #else
      build_lists_by_atom();      
      #endif

      #ifdef MEM_OPT_VERSION
      delEachAtomSigs();
      delChargeSpace();
      delMassSpace();
      delOtherEachAtomStructs();
      #endif
    }
    /*      END OF FUNCTION receive_Molecule    */


#ifdef MEM_OPT_VERSION
    //Well, the exclusion check signatures could also done on PE0 and
    //sent to other processors through send_Molecule/receive_Molecule 
    //two procedures.
    void Molecule::build_excl_check_signatures(){
       exclChkSigPool = new ExclusionCheck[exclSigPoolSize];
       for(int i=0; i<exclSigPoolSize; i++){
           ExclusionSignature *sig = &exclSigPool[i];
           ExclusionCheck *sigChk = &exclChkSigPool[i];
           if(sig->fullExclCnt){
               if(!sig->modExclCnt){ //only having fullExclusion
                   sigChk->min = sig->fullOffset[0];
                   sigChk->max = sig->fullOffset[sig->fullExclCnt-1];
               }else{ //have both full and modified exclusion
                   int fullMin, fullMax, modMin, modMax;
                   
                   fullMin = sig->fullOffset[0];
                   fullMax = sig->fullOffset[sig->fullExclCnt-1];
               
                   modMin = sig->modOffset[0];
                   modMax = sig->modOffset[sig->modExclCnt-1];
                   
                   if(fullMin < modMin)
                       sigChk->min = fullMin;
                   else
                       sigChk->min = modMin;
                   if(fullMax < modMax)
                       sigChk->max = modMax;
                   else
                       sigChk->max = fullMax;
               }        
           }else{
               if(sig->modExclCnt){
                   sigChk->min = sig->modOffset[0];
                   sigChk->max = sig->modOffset[sig->modExclCnt-1];
               }else{ //both count are 0
                   if(CkMyPe()==0)
                       printf("Warning: an empty exclusion signature with index %d!\n", i);               
                   continue;
               }
           }           

           if((sigChk->max-sigChk->min) > simParams->maxExclusionFlags){
	       printf("The distance of a exclusion check %d exceeds the value %d simulation parameter maxExclusionFlags specifies! Please increase the value\n", sigChk->max-sigChk->min, simParams->maxExclusionFlags);
               NAMD_die("Currently not supporting building exclusion check on the fly for memory optimized version!\n");
           }

           sigChk->flags = new char[sigChk->max-sigChk->min+1];
           memset(sigChk->flags, 0, sizeof(char)*(sigChk->max-sigChk->min+1));
           for(int j=0; j<sig->fullExclCnt; j++){
               int dist = sig->fullOffset[j] - sigChk->min;
               sigChk->flags[dist] = EXCHCK_FULL;
           }
           for(int j=0; j<sig->modExclCnt; j++){
               int dist = sig->modOffset[j] - sigChk->min;
               sigChk->flags[dist] = EXCHCK_MOD;
           }
       }
    }
#endif

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_lists_by_atom      */
    /*                  */
    /*  This function builds O(NumAtoms) arrays that store the bonds,   */
    /*  angles, dihedrals, and impropers, that each atom is involved in.    */
    /*  This is a space hog, but VERY fast.  This will certainly have to    */
    /*  change to make things scalable in memory, but for now, speed is the */
    /*  thing!                */
    /*                  */
    /************************************************************************/
#ifdef MEM_OPT_VERSION

    //In memory optimized version, this function is only performed on PE0 since other nodes
    //may have less memory than the master node.
    //3 tasks is performed in the function: 1) cluster building 2) atom signatures changed due to
    // fixed atoms and alch atoms 3) exclusion signatures (in the format of offset) changed due to
    //hydrogen group/alch/fixed atoms.
    void Molecule::build_lists_by_atom()
       
    {
       register int i;      //  Loop counter

       register int numFixedAtoms = this->numFixedAtoms;
       // if we want forces on fixed atoms then just pretend
       // there are none for the purposes of this routine;
       if ( simParams->fixedAtomsForces ) numFixedAtoms = 0;
                     
       const int pair_self = 
         simParams->pairInteractionOn ? simParams->pairInteractionSelf : 0;

       //Cluster building is done when compressing the psf file
       
       //The safety check for bonds has been already done
       //in the period of generating compressed psf file
       DebugM(3, "Building the atom list for each signature.\n");

       //build up tuplesByAtom structure where each atom's signature may be updated
       //because of fixed atoms, atom signature pool will be changed.       
       //The method to update the current signature pool:

       //For each atom signature, go through those atoms that belong to this signature
       //and split the current signature correspondingly. 
       // Note: 1.The splitted signatures originating from the same atom signature will always
       // be different. But two splitted signatures may be same if they originate from two different
       // atom signatures. Currently, this redundancy is not handled for the sake of convenience
       // and simplicity. I believe practically this will not blow up memory since
       // the space for atom signatures will always remain very low (~10K) 
       //       2. This method is used under the assumption that chance of splitting a current
       // atom signature is small

       // There is another way to do this: (The following code is an implementation of this
       // 1. Finding all the tuples that are not needed to be considered either for 
       // fixed number or for pair_self
       // 2. Build atom signatures for them
       // 3. Add these new signatures to the existing atom signature pool
       // -- Chao Mei       

       numCalcBonds = numBonds;
       numCalcAngles = numAngles;
       numCalcDihedrals = numDihedrals;
       numCalcImpropers = numImpropers;
       numCalcCrossterms = numCrossterms;
       if(numFixedAtoms || pair_self){ //either condition holds                            
           vector<AtomSignature> newAtomSigPool;
           //we assume there are always one empty atom signature which is indexed
           //at 0 in the newAtomSigPool but with atomSigPoolSize in the final combined
           //atom signature pool
           AtomSignature emptyAtomSig;
           newAtomSigPool.push_back(emptyAtomSig);

           //begin to building new atom signatures for those atoms
           //affected by fixed atoms and fepAtomFlags. Besides,
           //the affected atoms' signature index should be updated too
           AtomSignature newAtomSig;           
           for(i=0; i<numAtoms; i++){ //"i" is the id of atom
               AtomSignature *curAtomSig = &atomSigPool[eachAtomSig[i]];

               //deal with fep atoms
               if(pair_self && fepAtomFlags[i]!=1){
                   //all tuples associated with this atom are gone                   
                   //update the number of calculated tuples
                   numCalcBonds -= curAtomSig->bondCnt;
                   numCalcAngles -= curAtomSig->angleCnt;
                   numCalcDihedrals -= curAtomSig->dihedralCnt;
                   numCalcImpropers -= curAtomSig->improperCnt;
                   numCalcCrossterms -= curAtomSig->crosstermCnt;

                   //update this atom's signature
                   eachAtomSig[i] = (Index)atomSigPoolSize;                   
                   
                   continue;
               }

               //deal with fixed atoms
               int sigChanged = 0;
               if(numFixedAtoms && fixedAtomFlags[i]){
                   //examine all the tuples associated with this atom                   
                   newAtomSig = *curAtomSig;

                   //1. bonds
                   int tupleCnt = curAtomSig->bondCnt;
                   TupleSignature *tupleSigs = curAtomSig->bondSigs;
                   for(int j=0; j<tupleCnt; j++){
                       int atom2 = i+tupleSigs[j].offset[0];
                       if(fixedAtomFlags[atom2]){ // find a fixed tuple                           
                           newAtomSig.bondSigs[j].setEmpty();
                           numCalcBonds--;
                           sigChanged = 1;
                       }
                   }

                   //2. angles
                   tupleCnt = curAtomSig->angleCnt;
                   tupleSigs = curAtomSig->angleSigs;
                   for(int j=0; j<tupleCnt; j++){
                       int atom2 = i+tupleSigs[j].offset[0];
                       int atom3 = i+tupleSigs[j].offset[1];
                       if(fixedAtomFlags[atom2] && fixedAtomFlags[atom3]){ // find a fixed tuple                           
                           newAtomSig.angleSigs[j].setEmpty();
                           numCalcAngles--;
                           sigChanged = 1;
                       }
                   }

                   //3. dihedrals
                   tupleCnt = curAtomSig->dihedralCnt;
                   tupleSigs = curAtomSig->dihedralSigs;
                   for(int j=0; j<tupleCnt; j++){
                       int atom2 = i+tupleSigs[j].offset[0];
                       int atom3 = i+tupleSigs[j].offset[1];
                       int atom4 = i+tupleSigs[j].offset[2];
                       if(fixedAtomFlags[atom2] && fixedAtomFlags[atom3]
                          && fixedAtomFlags[atom4]){ // find a fixed tuple                           
                           newAtomSig.dihedralSigs[j].setEmpty();
                           numCalcDihedrals--;
                           sigChanged = 1;
                       }
                   }

                   //4. impropers
                   tupleCnt = curAtomSig->improperCnt;
                   tupleSigs = curAtomSig->improperSigs;
                   for(int j=0; j<tupleCnt; j++){
                       int atom2 = i+tupleSigs[j].offset[0];
                       int atom3 = i+tupleSigs[j].offset[1];
                       int atom4 = i+tupleSigs[j].offset[2];
                       if(fixedAtomFlags[atom2] && fixedAtomFlags[atom3]
                          && fixedAtomFlags[atom4]){ // find a fixed tuple                           
                           newAtomSig.improperSigs[j].setEmpty();
                           numCalcImpropers--;
                           sigChanged = 1;
                       }
                   }

                   //5. crossterms
                   tupleCnt = curAtomSig->crosstermCnt;
                   tupleSigs = curAtomSig->crosstermSigs;
                   for(int j=0; j<tupleCnt; j++){
                       int atom2 = i+tupleSigs[j].offset[0];
                       int atom3 = i+tupleSigs[j].offset[1];
                       int atom4 = i+tupleSigs[j].offset[2];
                       int atom5 = i+tupleSigs[j].offset[3];
                       int atom6 = i+tupleSigs[j].offset[4];
                       int atom7 = i+tupleSigs[j].offset[5];
                       int atom8 = i+tupleSigs[j].offset[6];

                       if(fixedAtomFlags[atom2] && fixedAtomFlags[atom3]
                          && fixedAtomFlags[atom4] && fixedAtomFlags[atom5]
                          && fixedAtomFlags[atom6] && fixedAtomFlags[atom7]
                          && fixedAtomFlags[atom8]){ // find a fixed tuple                           
                           newAtomSig.crosstermSigs[j].setEmpty();
                           numCalcCrossterms--;
                           sigChanged = 1;
                       }
                   }

                   //update the newAtomSig to remove all empty tuple signatures
                   if(sigChanged){
                       newAtomSig.removeEmptyTupleSigs();

                       //add the new atom signature to the newAtomSigPool if it is not in it
                       int poolIndex = lookupCstPool(newAtomSigPool, newAtomSig);
                       if(poolIndex==-1){
                           poolIndex = newAtomSigPool.size();
                           newAtomSigPool.push_back(newAtomSig);
                       }
                       eachAtomSig[i] = (Index)(poolIndex+atomSigPoolSize);
                   }
               }
           }
           //update the atom signature pool from the new one                      
           AtomSignature *tmpAtomSigPool = new AtomSignature[atomSigPoolSize+newAtomSigPool.size()];
           for(i=0; i<atomSigPoolSize; i++)
               tmpAtomSigPool[i] = atomSigPool[i];           
           for(i=0; i<newAtomSigPool.size(); i++)
               tmpAtomSigPool[atomSigPoolSize+i] = newAtomSigPool[i];

           delete [] atomSigPool;
           atomSigPoolSize += newAtomSigPool.size();
           atomSigPool = tmpAtomSigPool;
           newAtomSigPool.clear();

           if(CkMyPe()==0){
               printf("Current number of atom signatures: %d\n", atomSigPoolSize);
           }
       }
          
       DebugM(3,"Building exclusion data.\n");


       //  Build the arrays of exclusions for each atom       
       //inside this function call, unnecessary exclusions are eliminated in the
       //case of fep atoms and fixed atoms
       if(CkMyPe()==0){
           CkPrintf("Current exclusion signatures: %d\n", exclSigPoolSize);
       }
       build_exclusions();       
       
       if(CkMyPe()==0){
           CkPrintf("After stripping, exclusion signatures: %d\n", exclSigPoolSize);
       }
       
       //Finally build the exclChkSigPool of type ExclusionCheck which is separated
       //from this function for saving memory
       //build_excl_check_signature();
    }
#else
    void Molecule::build_lists_by_atom()
       
    {
       register int i;      //  Loop counter

       register int numFixedAtoms = this->numFixedAtoms;
       // if we want forces on fixed atoms then just pretend
       // there are none for the purposes of this routine;
       if ( simParams->fixedAtomsForces ) numFixedAtoms = 0;

//fepb
//     int numFepInitial = this->numFepInitial;
//     int numFepFinal = this->numFepFinal;
//fepe
       tmpArena = new ObjectArena<int32>;
       bondsWithAtom = new int32 *[numAtoms];
       cluster = new int32 [numAtoms];
       clusterSize = new int32 [numAtoms];

       bondsByAtom = new int32 *[numAtoms];
       anglesByAtom = new int32 *[numAtoms];
       dihedralsByAtom = new int32 *[numAtoms];
       impropersByAtom = new int32 *[numAtoms];
       crosstermsByAtom = new int32 *[numAtoms];

       exclusionsByAtom = new int32 *[numAtoms];
       fullExclusionsByAtom = new int32 *[numAtoms];
       modExclusionsByAtom = new int32 *[numAtoms];

       int32 *byAtomSize = new int32[numAtoms];

       const int pair_self = 
         simParams->pairInteractionOn ? simParams->pairInteractionSelf : 0;

       DebugM(3,"Building bond lists.\n");
    
       //  Build the bond lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       for (i=0; i<numRealBonds; i++)
       {
         byAtomSize[bonds[i].atom1]++;
         byAtomSize[bonds[i].atom2]++;
       }
       for (i=0; i<numAtoms; i++)
       {
         bondsWithAtom[i] = tmpArena->getNewArray(byAtomSize[i]+1);
         bondsWithAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numRealBonds; i++)
       {
         int a1 = bonds[i].atom1;
         int a2 = bonds[i].atom2;
         bondsWithAtom[a1][byAtomSize[a1]++] = i;
         bondsWithAtom[a2][byAtomSize[a2]++] = i;
       }

       //  Build cluster information (contiguous clusters)
       for (i=0; i<numAtoms; i++) {
         cluster[i] = i;
       }
       for (i=0; i<numAtoms; i++) {
         int ci = i;
         while ( cluster[ci] != ci ) ci = cluster[ci];
         for ( int32 *b = bondsWithAtom[i]; *b != -1; ++b ) {
           int a = bonds[*b].atom1;
           if ( a == i ) a = bonds[*b].atom2;
           if ( a > i ) {
             int ca = a;
             while ( cluster[ca] != ca ) ca = cluster[ca];
             if ( ca > ci ) cluster[ca] = cluster[ci];
             else cluster[ci] = cluster[ca];
           }
         }
       }
       while ( 1 ) {
         int allok = 1;
         for (i=0; i<numAtoms; i++) {
           int ci = cluster[i];
           if ( cluster[ci] != ci ) {
             allok = 0;
             cluster[i] = cluster[ci];
           }
         }
         if ( allok ) break;
       }
       
       for (i=0; i<numAtoms; i++) {
         clusterSize[i] = 0;
       }       
       for (i=0; i<numAtoms; i++) {           
         clusterSize[cluster[i]] += 1;
       }

/*
       //Getting number of clusters for debugging
       int numClusters=0;
       for(int i=0; i<numAtoms; i++){
           if(clusterSize[i]!=0) numClusters++;
       }
       printf("Num of clusters: %d\n", numClusters);
*/

       //  Build the bond lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcBonds = 0;
       for (i=0; i<numBonds; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[bonds[i].atom1]
                            && fixedAtomFlags[bonds[i].atom2] ) continue;
   
         if ( pair_self && fepAtomFlags[bonds[i].atom1] != 1) continue;
         byAtomSize[bonds[i].atom1]++;
         numCalcBonds++;
       }
       for (i=0; i<numAtoms; i++)
       {
         bondsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         bondsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numBonds; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[bonds[i].atom1]
                            && fixedAtomFlags[bonds[i].atom2] ) continue;
         if ( pair_self && fepAtomFlags[bonds[i].atom1] != 1) continue;
         int a1 = bonds[i].atom1;
         bondsByAtom[a1][byAtomSize[a1]++] = i;
       }
       for (i=0; i<numBonds; ++i) {
         int a1 = bonds[i].atom1;
         int a2 = bonds[i].atom2;
         int j;
         if ( a1 == a2 ) {
           char buff[512];
           sprintf(buff,"Atom %d is bonded to itself", a1+1);
           NAMD_die(buff);
         }
         for ( j = 0; j < byAtomSize[a1]; ++j ) {
           int b = bondsByAtom[a1][j];
           int ba1 = bonds[b].atom1;
           int ba2 = bonds[b].atom2;
           if ( b != i && ( (ba1==a1 && ba2==a2) || (ba1==a2 && ba2==a1) ) ) {
             char buff[512];
             sprintf(buff,"Duplicate bond from atom %d to atom %d", a1+1, a2+1);
             NAMD_die(buff);
           }
         }
         for ( j = 0; j < byAtomSize[a2]; ++j ) {
           int b = bondsByAtom[a2][j];
           int ba1 = bonds[b].atom1;
           int ba2 = bonds[b].atom2;
           if ( b != i && ( (ba1==a1 && ba2==a2) || (ba1==a2 && ba2==a1) ) ) {
             char buff[512];
             sprintf(buff,"Duplicate bond from atom %d to atom %d", a1+1, a2+1);
             NAMD_die(buff);
           }
         }
       }
       
       DebugM(3,"Building angle lists.\n");
    
       //  Build the angle lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcAngles = 0;
       for (i=0; i<numAngles; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[angles[i].atom1]
                            && fixedAtomFlags[angles[i].atom2]
                            && fixedAtomFlags[angles[i].atom3] ) continue;
         if ( pair_self && fepAtomFlags[angles[i].atom1] != 1) continue;
         byAtomSize[angles[i].atom1]++;
         numCalcAngles++;
       }
       for (i=0; i<numAtoms; i++)
       {
         anglesByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         anglesByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numAngles; i++)
       {
         if ( pair_self && fepAtomFlags[angles[i].atom1] != 1) continue;
         if ( numFixedAtoms && fixedAtomFlags[angles[i].atom1]
                            && fixedAtomFlags[angles[i].atom2]
                            && fixedAtomFlags[angles[i].atom3] ) continue;
         int a1 = angles[i].atom1;
         anglesByAtom[a1][byAtomSize[a1]++] = i;
       }
       
       DebugM(3,"Building improper lists.\n");
    
       //  Build the improper lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcImpropers = 0;
       for (i=0; i<numImpropers; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[impropers[i].atom1]
                            && fixedAtomFlags[impropers[i].atom2]
                            && fixedAtomFlags[impropers[i].atom3]
                            && fixedAtomFlags[impropers[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[impropers[i].atom1] != 1) continue;
         byAtomSize[impropers[i].atom1]++;
         numCalcImpropers++;
       }
       for (i=0; i<numAtoms; i++)
       {
         impropersByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         impropersByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numImpropers; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[impropers[i].atom1]
                            && fixedAtomFlags[impropers[i].atom2]
                            && fixedAtomFlags[impropers[i].atom3]
                            && fixedAtomFlags[impropers[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[impropers[i].atom1] != 1) continue;
         int a1 = impropers[i].atom1;
         impropersByAtom[a1][byAtomSize[a1]++] = i;
       }
       
       DebugM(3,"Building dihedral lists.\n");
    
       //  Build the dihedral lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcDihedrals = 0;
       for (i=0; i<numDihedrals; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[dihedrals[i].atom1]
                            && fixedAtomFlags[dihedrals[i].atom2]
                            && fixedAtomFlags[dihedrals[i].atom3]
                            && fixedAtomFlags[dihedrals[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[dihedrals[i].atom1] != 1) continue;
         byAtomSize[dihedrals[i].atom1]++;
         numCalcDihedrals++;
       }
       for (i=0; i<numAtoms; i++)
       {
         dihedralsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         dihedralsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numDihedrals; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[dihedrals[i].atom1]
                            && fixedAtomFlags[dihedrals[i].atom2]
                            && fixedAtomFlags[dihedrals[i].atom3]
                            && fixedAtomFlags[dihedrals[i].atom4] ) continue;
         if ( pair_self && fepAtomFlags[dihedrals[i].atom1] != 1) continue;
         int a1 = dihedrals[i].atom1;
         dihedralsByAtom[a1][byAtomSize[a1]++] = i;
       }
    
       DebugM(3,"Building crossterm lists.\n");
    
       //  Build the crossterm lists
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcCrossterms = 0;
       for (i=0; i<numCrossterms; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[crossterms[i].atom1]
                            && fixedAtomFlags[crossterms[i].atom2]
                            && fixedAtomFlags[crossterms[i].atom3]
                            && fixedAtomFlags[crossterms[i].atom4]
                            && fixedAtomFlags[crossterms[i].atom5]
                            && fixedAtomFlags[crossterms[i].atom6]
                            && fixedAtomFlags[crossterms[i].atom7]
                            && fixedAtomFlags[crossterms[i].atom8] ) continue;
         if ( pair_self && fepAtomFlags[crossterms[i].atom1] != 1) continue;
         byAtomSize[crossterms[i].atom1]++;
         numCalcCrossterms++;
       }
       for (i=0; i<numAtoms; i++)
       {
         crosstermsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         crosstermsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numCrossterms; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[crossterms[i].atom1]
                            && fixedAtomFlags[crossterms[i].atom2]
                            && fixedAtomFlags[crossterms[i].atom3]
                            && fixedAtomFlags[crossterms[i].atom4]
                            && fixedAtomFlags[crossterms[i].atom5]
                            && fixedAtomFlags[crossterms[i].atom6]
                            && fixedAtomFlags[crossterms[i].atom7]
                            && fixedAtomFlags[crossterms[i].atom8] ) continue;
         if ( pair_self && fepAtomFlags[crossterms[i].atom1] != 1) continue;
         int a1 = crossterms[i].atom1;
         crosstermsByAtom[a1][byAtomSize[a1]++] = i;
       }

       // DRUDE: init lphostIndexes array
       if (simParams->drudeOn) {
         // allocate lone pair host index array only if we need it!
         DebugM(3,"Initializing lone pair host index array.\n");
         lphostIndexes = new int32[numAtoms];
         for (i = 0;  i < numAtoms;  i++) {
           lphostIndexes[i] = -1;
         }
         for (i = 0;  i < numLphosts;  i++) {
           int32 index = lphosts[i].atom1;
           lphostIndexes[index] = i;
         }
       }
       // DRUDE
    
       DebugM(3,"Building exclusion data.\n");
    
       //  Build the arrays of exclusions for each atom
       build_exclusions();

       //  Remove temporary structures
       delete [] bondsWithAtom;  bondsWithAtom = 0;
       delete tmpArena;  tmpArena = 0;

       if (exclusions != NULL)
      delete [] exclusions;

       // 1-4 exclusions which are also fully excluded were eliminated by hash table
       numTotalExclusions = exclusionSet.size();
       exclusions = new Exclusion[numTotalExclusions];
       UniqueSetIter<Exclusion> exclIter(exclusionSet);
       for ( exclIter=exclIter.begin(),i=0; exclIter != exclIter.end(); exclIter++,i++ )
       {
         exclusions[i] = *exclIter;
       }
       // Free exclusionSet storage
       // exclusionSet.clear(1);
       exclusionSet.clear();

       DebugM(3,"Building exclusion lists.\n");
    
       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
       }
       numCalcExclusions = 0;
       for (i=0; i<numTotalExclusions; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
         byAtomSize[exclusions[i].atom1]++;
         numCalcExclusions++;
       }
       for (i=0; i<numAtoms; i++)
       {
         exclusionsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         exclusionsByAtom[i][byAtomSize[i]] = -1;
         byAtomSize[i] = 0;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
         int a1 = exclusions[i].atom1;
         exclusionsByAtom[a1][byAtomSize[a1]++] = i;
       }

       int32 *byAtomSize2 = new int32[numAtoms];

       for (i=0; i<numAtoms; i++)
       {
         byAtomSize[i] = 0;
         byAtomSize2[i] = 0;
       }

       for (i=0; i<numTotalExclusions; i++)
       {
         if ( numFixedAtoms && fixedAtomFlags[exclusions[i].atom1]
                            && fixedAtomFlags[exclusions[i].atom2] ) continue;
         if ( exclusions[i].modified ) {
           byAtomSize2[exclusions[i].atom1]++;
           byAtomSize2[exclusions[i].atom2]++;
         } else {
           byAtomSize[exclusions[i].atom1]++;
           byAtomSize[exclusions[i].atom2]++;
         }
       }

       for (i=0; i<numAtoms; i++)
       {
         fullExclusionsByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
         fullExclusionsByAtom[i][0] = 0;
         modExclusionsByAtom[i] = arena->getNewArray(byAtomSize2[i]+1);
         modExclusionsByAtom[i][0] = 0;
       }

       for (i=0; i<numTotalExclusions; i++)
       {
         int a1 = exclusions[i].atom1;
         int a2 = exclusions[i].atom2;
         if ( numFixedAtoms && fixedAtomFlags[a1]
                            && fixedAtomFlags[a2] ) continue;
         int32 *l1, *l2;
         if ( exclusions[i].modified ) {
           l1 = modExclusionsByAtom[a1];
           l2 = modExclusionsByAtom[a2];
         } else {
           l1 = fullExclusionsByAtom[a1];
           l2 = fullExclusionsByAtom[a2];
         }
         l1[++(*l1)] = a2;
         l2[++(*l2)] = a1;
       }

       // DRUDE
       if (is_drude_psf) {

         // build Thole (screened Coulomb) correction terms;
         // they are constructed implicitly from exclusions

         // free the previous Thole array if already allocated
         if (tholes != NULL) delete[] tholes;
         numTholes = 0;

         // count the number of Thole terms
         for (i = 0;  i < numTotalExclusions;  i++) {
           int a1 = exclusions[i].atom1;
           int a2 = exclusions[i].atom2;
           if (a2 < numAtoms-1 && is_drude(a1+1) && is_drude(a2+1)) {
             numTholes++;
           }
         }

         // allocate space for Thole terms
         if (numTholes != 0) tholes = new Thole[numTholes];
         else tholes = NULL;
         int nt = 0;

         Real c = COULOMB*simParams->nonbondedScaling/simParams->dielectric;

         // store Thole terms
         for (i = 0;  i < numTotalExclusions;  i++) {
           int a1 = exclusions[i].atom1;
           int a2 = exclusions[i].atom2;
           // exclusions are stored with a1 < a2
           if (a2 < numAtoms-1 && is_drude(a1+1) && is_drude(a2+1)) {
             Real thsum = drudeConsts[a1].thole + drudeConsts[a2].thole;
             Real aprod = drudeConsts[a1].alpha * drudeConsts[a2].alpha;
             // guard against having alpha==0
             Real apower = (aprod <= 0 ? 0 : powf(aprod, -1.f/6));
             tholes[nt].atom1 = a1;
             tholes[nt].atom2 = a1+1;
             tholes[nt].atom3 = a2;
             tholes[nt].atom4 = a2+1;
             tholes[nt].aa = apower * thsum;
             tholes[nt].qq = c * atoms[a1+1].charge * atoms[a2+1].charge;
             nt++;
           }
         }

         // build Thole lists by atom
         DebugM(3, "Building Thole correction term lists.\n");
         tholesByAtom = new int32 *[numAtoms];

         for (i = 0;  i < numAtoms;  i++) {
           byAtomSize[i] = 0;
         }
         numCalcTholes = 0;
         for (i = 0;  i < numTholes;  i++) {
           if ( numFixedAtoms && fixedAtomFlags[tholes[i].atom1]
                              && fixedAtomFlags[tholes[i].atom2]
                              && fixedAtomFlags[tholes[i].atom3]
                              && fixedAtomFlags[tholes[i].atom4] ) continue;
           if ( pair_self && fepAtomFlags[tholes[i].atom1] != 1) continue;
           byAtomSize[tholes[i].atom1]++;
           numCalcTholes++;
         }
         for (i = 0;  i < numAtoms;  i++) {
           tholesByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
           tholesByAtom[i][byAtomSize[i]] = -1;
           byAtomSize[i] = 0;
         }
         for (i = 0;  i < numTholes;  i++) {
           if ( numFixedAtoms && fixedAtomFlags[tholes[i].atom1]
                              && fixedAtomFlags[tholes[i].atom2]
                              && fixedAtomFlags[tholes[i].atom3]
                              && fixedAtomFlags[tholes[i].atom4] ) continue;
           if ( pair_self && fepAtomFlags[tholes[i].atom1] != 1) continue;
           int a1 = tholes[i].atom1;
           tholesByAtom[a1][byAtomSize[a1]++] = i;
         }

         // build anisotropic lists by atom
         DebugM(3, "Building anisotropic term lists.\n");
         anisosByAtom = new int32 *[numAtoms];

         for (i = 0;  i < numAtoms;  i++) {
           byAtomSize[i] = 0;
         }
         numCalcAnisos = 0;
         for (i = 0;  i < numAnisos;  i++) {
           if ( numFixedAtoms && fixedAtomFlags[anisos[i].atom1]
                              && fixedAtomFlags[anisos[i].atom2]
                              && fixedAtomFlags[anisos[i].atom3]
                              && fixedAtomFlags[anisos[i].atom4] ) continue;
           if ( pair_self && fepAtomFlags[anisos[i].atom1] != 1) continue;
           byAtomSize[anisos[i].atom1]++;
           numCalcAnisos++;
         }
         for (i = 0;  i < numAtoms;  i++) {
           anisosByAtom[i] = arena->getNewArray(byAtomSize[i]+1);
           anisosByAtom[i][byAtomSize[i]] = -1;
           byAtomSize[i] = 0;
         }
         for (i = 0;  i < numAnisos;  i++) {
           if ( numFixedAtoms && fixedAtomFlags[anisos[i].atom1]
                              && fixedAtomFlags[anisos[i].atom2]
                              && fixedAtomFlags[anisos[i].atom3]
                              && fixedAtomFlags[anisos[i].atom4] ) continue;
           if ( pair_self && fepAtomFlags[anisos[i].atom1] != 1) continue;
           int a1 = anisos[i].atom1;
           anisosByAtom[a1][byAtomSize[a1]++] = i;
         }

       }
       // DRUDE

       delete [] byAtomSize;  byAtomSize = 0;
       delete [] byAtomSize2;  byAtomSize2 = 0;


       //  Allocate an array to hold the exclusions for each atom
       all_exclusions = new ExclusionCheck[numAtoms];

       for (i=0; i<numAtoms; i++)
       {
         all_exclusions[i].min = numAtoms;
         all_exclusions[i].max = -1;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         // first atom should alway have lower number!
         int a1 = exclusions[i].atom1;
         int a2 = exclusions[i].atom2;
         if ( numFixedAtoms && fixedAtomFlags[a1]
                            && fixedAtomFlags[a2] ) continue;
         if ( all_exclusions[a1].min > a2 ) all_exclusions[a1].min = a2;
         if ( all_exclusions[a2].min > a1 ) all_exclusions[a2].min = a1;
         if ( a2 > all_exclusions[a1].max ) all_exclusions[a1].max = a2;
         if ( a1 > all_exclusions[a2].max ) all_exclusions[a2].max = a1;
       }
       int exclmem = 0;
       int maxExclusionFlags = simParams->maxExclusionFlags;
       for (i=0; i<numAtoms; i++)
       {
         int s = all_exclusions[i].max - all_exclusions[i].min + 1;
         if ( all_exclusions[i].max != -1 ) {
           if ( s < maxExclusionFlags ) {
             char *f = all_exclusions[i].flags = exclArena->getNewArray(s);
             for ( int k=0; k<s; ++k ) f[k] = 0;
             exclmem += s;
           } else {
             all_exclusions[i].flags = 0;  // need to build on the fly
           }
         } else {
           all_exclusions[i].flags = (char*)-1; // should never dereference
         }
       }
       if ( 0 ) {
         iout << iINFO << numTotalExclusions << " exclusions consume "
            << exclmem << " bytes.\n" << endi;
       }
       for (i=0; i<numTotalExclusions; i++)
       {
         int a1 = exclusions[i].atom1;
         int a2 = exclusions[i].atom2;
         if ( numFixedAtoms && fixedAtomFlags[a1]
                            && fixedAtomFlags[a2] ) continue;
         if ( exclusions[i].modified ) {
           if ( all_exclusions[a1].flags )
             all_exclusions[a1].flags[a2-all_exclusions[a1].min] = EXCHCK_MOD;
           if ( all_exclusions[a2].flags )
             all_exclusions[a2].flags[a1-all_exclusions[a2].min] = EXCHCK_MOD;
         } else {
           if ( all_exclusions[a1].flags )
             all_exclusions[a1].flags[a2-all_exclusions[a1].min] = EXCHCK_FULL;
           if ( all_exclusions[a2].flags )
             all_exclusions[a2].flags[a1-all_exclusions[a2].min] = EXCHCK_FULL;
         }
       }

    }
#endif

    /*    END OF FUNCTION build_lists_by_atom    */

    /****************************************************************/
    /*                */
    /*      FUNCTION build_exclusions    */
    /*                */
    /*  This function builds a list of all the exlcusions       */
    /*  atoms.  These lists include explicit exclusions as well as  */
    /*  exclusions that are calculated based on the bonded structure*/
    /*  and the exclusion flag.  For each pair of atoms that are    */
    /*  excluded, the larger of the 2 atom indexes is stored in the */
    /*  array of the smaller index.  All the arrays are not sorted. */
    /*  Then to determine if two atoms have an exclusion, a linear  */
    /*  search is done on the array of the atom with the smaller    */
    /*  index for the larger index.          */
    /*  If the exclusion policy is set to scaled1-4, there are  */
    /*  actually two lists built.  One contains the pairs of atoms  */
    /*  that are to be exlcuded (i.e., explicit exclusions, 1-2,    */
    /*  and 1-3 interactions) and the other contains just the 1-4   */
    /*  interactions, since they will need to be determined   */
    /*  independantly of the other exclusions.      */
    /*                */
    /****************************************************************/

    void Molecule::build_exclusions()
    {
#ifdef MEM_OPT_VERSION
        numCalcExclusions = numTotalExclusions*2;        

        //stripHGroupExcl() is not needed because the function purpose 
        //is to reduce the memory usage of exclusions. Since now the exclusion
        //signature is used, such function is not critical now
        //What about the following function? Need to be confirmed --Chao Mei
        
        ExclusionSettings exclude_flag;    //  Exclusion policy
		exclude_flag = simParams->exclude;
		int stripHGroupExclFlag = (simParams->splitPatch == SPLIT_PATCH_HYDROGEN);
        if (!simParams->amberOn || !simParams->readExclusions){         	
	        switch (exclude_flag)
	        {
	         case NONE:	           
	         case ONETWO:           
	           break;
	          case ONETHREE:            		    	
	          case ONEFOUR:		    	
	          case SCALED14:            
		    	if ( stripHGroupExclFlag ) stripHGroupExcl();
	            break;
	        }
		}
		else if (stripHGroupExclFlag && exclude_flag!=NONE && exclude_flag!=ONETWO)
			stripHGroupExcl();               

        stripFepFixedExcl();
        return;
#else
      register int i;          //  Loop counter
      ExclusionSettings exclude_flag;    //  Exclusion policy

      exclude_flag = simParams->exclude;
      int stripHGroupExclFlag = (simParams->splitPatch == SPLIT_PATCH_HYDROGEN);

      //  Go through the explicit exclusions and add them to the arrays
      for (i=0; i<numExclusions; i++)
      {
        exclusionSet.add(exclusions[i]);
      }

      // If this is AMBER force field, and readExclusions is TRUE,
      // then all the exclusions were read from parm file, and we
      // shouldn't generate any of them.
      if (!simParams->amberOn || !simParams->readExclusions)
      { //  Now calculate the bonded exlcusions based on the exclusion policy
        switch (exclude_flag)
        {
         case NONE:
           break;
         case ONETWO:
           build12excl();
           break;
          case ONETHREE:
            build12excl();
            build13excl();
	    if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
          case ONEFOUR:
            build12excl();
            build13excl();
            build14excl(0);
	    if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
          case SCALED14:
            build12excl();
            build13excl();
            build14excl(1);
	    if ( stripHGroupExclFlag ) stripHGroupExcl();
            break;
        }
      }
      else if (stripHGroupExclFlag && exclude_flag!=NONE && exclude_flag!=ONETWO)
        stripHGroupExcl();

      stripFepExcl();

      // DRUDE
      if (is_drude_psf) build_inherited_excl();
#endif
    }
    /*      END OF FUNCTION build_exclusions    */


    // Extend exclusions for the Drude model.  The Drude model is generally
    // used with the 1-3 exclusion policy, although the code below also
    // supports the 1-2 exclusion policy.  The use of light (or massless)
    // pseudo-atoms requires the introduction of extra exclusions.
    //
    // Here is the algorithm for determining Drude model exclusions:
    // (1)  Each Drude particle and each lone pair has a single parent atom.
    //      The parent atom must be a heavy atom.
    // (2)  Each Drude particle and lone pair inherit the exclusions of its
    //      parent atom.
    // (3)  If two heavy atoms are excluded and they both have either a
    //      Drude particle or a lone pair, the these light (or massless)
    //      particles are also excluded from interacting with each other.
    void Molecule::build_inherited_excl(void) {
#ifdef MEM_OPT_VERSION
      NAMD_die("Drude and LP particles not supported in memopt version.");
#else
      ExclusionSettings exclude_flag = simParams->exclude;
      int32 *bond1, *bond2, *bond3, *bond4;
      int32 i, j, mid1, mid2, mid3;

      if (exclude_flag == ONEFOUR || exclude_flag == SCALED14) {
        NAMD_die("DRUDE MODEL SUPPORTS ONLY UP TO 1-3 EXCLUSION POLICY");
      }

      if (exclude_flag == NONE) return;

      // validate that each Drude or lone pair particle
      // has a unique parent that is a heavy atom
      for (i = 0;  i < numAtoms;  i++) {

        if (is_drude(i) || is_lp(i)) {
          // find parent (heavy) atom of particle i
          bond1 = bondsWithAtom[i];

          if (-1 == *bond1) {  // i must have one bond
            char err_msg[512];
            const char *idescrip = (is_drude(i) ? "DRUDE" : "LONE PAIR");
            sprintf(err_msg, "FOUND ISOLATED %s PARTICLE %d", idescrip, i+1);
            NAMD_die(err_msg);
          }
          if (-1 != *(bond1+1)) {  // and only one bond
            char err_msg[512];
            const char *idescrip = (is_drude(i) ? "DRUDE" : "LONE PAIR");
            sprintf(err_msg, "FOUND MULTIPLY LINKED %s PARTICLE %d",
                idescrip, i+1);
            NAMD_die(err_msg);
          }

          // mid1 is parent of particle i
          mid1 = bonds[*bond1].atom1;
          if (mid1 == i) mid1 = bonds[*bond1].atom2;

          // make sure that mid1 is a heavy atom
          if (is_drude(mid1) || is_lp(mid1) || is_hydrogen(mid1)) {
            char err_msg[512];
            const char *idescrip = (is_drude(i) ? "DRUDE" : "LONE PAIR");
            sprintf(err_msg, "PARENT ATOM %d of %s PARTICLE %d "
                "IS NOT HEAVY ATOM", mid1+1, idescrip, i+1);
            NAMD_die(err_msg);
          }

          // follow build14excl() code
          bond2 = bondsWithAtom[mid1];

          // loop through all the bonds connected to atom mid1
          while (*bond2 != -1) {
            if (bonds[*bond2].atom1 == mid1) {
              mid2 = bonds[*bond2].atom2;
            }
            else {
              mid2 = bonds[*bond2].atom1;
            }

            // Make sure that we don't double back to where we started from.
            // Doing so causes strange behavior.
            if (mid2 == i) {
              bond2++;
              continue;
            }

            if (exclude_flag == ONETWO) {
              // add (i,mid2) as an exclusion
              if (i < mid2) {
                exclusionSet.add(Exclusion(i, mid2));
              }
              else {
                exclusionSet.add(Exclusion(mid2, i));
              }

              // also exclude any Drude particles or LPs bonded to mid2
              bond3 = bondsWithAtom[mid2];
              while (*bond3 != -1) {
                j = bonds[*bond3].atom1;
                if (is_drude(j) || is_lp(j)) {
                  if      (i < j) exclusionSet.add(Exclusion(i, j));
                  else if (j < i) exclusionSet.add(Exclusion(j, i));
                }
                j = bonds[*bond3].atom2;
                if (is_drude(j) || is_lp(j)) {
                  if      (i < j) exclusionSet.add(Exclusion(i, j));
                  else if (j < i) exclusionSet.add(Exclusion(j, i));
                }
                bond3++;
              }
            }
            else {  // exclude_flag == ONETHREE

              bond3 = bondsWithAtom[mid2];

              // loop through all the bonds connected to mid2
              while (*bond3 != -1) {

                if (bonds[*bond3].atom1 == mid2) {
                  mid3 = bonds[*bond3].atom2;
                }
                else {
                  mid3 = bonds[*bond3].atom1;
                }

                // Make sure we don't double back to where we started.
                // Doing so causes strange behavior.
                if (mid3 == mid1) {
                  bond3++;
                  continue;
                }

                // add (i,mid3) as an exclusion
                if (i < mid3) {
                  exclusionSet.add(Exclusion(i, mid3));
                }
                else if (mid3 < i) {
                  exclusionSet.add(Exclusion(mid3, i));
                }

                // also exclude any Drude particles or LPs bonded to mid3
                bond4 = bondsWithAtom[mid3];
                while (*bond4 != -1) {
                  j = bonds[*bond4].atom1;
                  if (is_drude(j) || is_lp(j)) {
                    if      (i < j) exclusionSet.add(Exclusion(i, j));
                    else if (j < i) exclusionSet.add(Exclusion(j, i));
                  }
                  j = bonds[*bond4].atom2;
                  if (is_drude(j) || is_lp(j)) {
                    if      (i < j) exclusionSet.add(Exclusion(i, j));
                    else if (j < i) exclusionSet.add(Exclusion(j, i));
                  }
                  bond4++;
                }

                ++bond3;
              } // while bond3

            } // else ONETHREE
           
            ++bond2;
          } // while bond2

        } // if i is Drude or LP

      } // for i

#endif
    } 
    // DRUDE


    /************************************************************************/
    /*                  */
    /*      FUNCTION build12excl        */
    /*                  */
    /************************************************************************/

    void Molecule::build12excl(void)
       
    {
#ifndef MEM_OPT_VERSION
       int32 *current_val;  //  Current value to check
       register int i;    //  Loop counter to loop through all atoms
       
       //  Loop through all the atoms marking the bonded interactions for each one
       for (i=0; i<numAtoms; i++)
       {
      current_val = bondsWithAtom[i];
       
      //  Loop through all the bonds for this atom
      while (*current_val != -1)
      {
         if (bonds[*current_val].atom1 == i)
         {
      if (i<bonds[*current_val].atom2)
      {
         exclusionSet.add(Exclusion(i,bonds[*current_val].atom2));
      }
         }
         else
         {
      if (i<bonds[*current_val].atom1)
      {
         exclusionSet.add(Exclusion(i,bonds[*current_val].atom1));
      }
         }
    
         ++current_val;
      }
       }
#endif
    }
    /*      END OF FUNCTION build12excl      */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build13excl        */
    /*                  */
    /************************************************************************/

    void Molecule::build13excl(void)
       
    {
#ifndef MEM_OPT_VERSION
       int32 *bond1, *bond2;  //  The two bonds being checked
       int middle_atom;  //  Common third atom
       register int i;    //  Loop counter to loop through all atoms
       
       //  Loop through all the atoms looking at the bonded connections
       //  for each one
       for (i=0; i<numAtoms; i++)
       {
       bond1 = bondsWithAtom[i];
       
       //  Loop through all the bonds directly connect to atom i
       while (*bond1 != -1)
       {
        if (bonds[*bond1].atom1 == i)
        {
          middle_atom=bonds[*bond1].atom2;
        }
        else
        {
          middle_atom=bonds[*bond1].atom1;
        }

        bond2 = bondsWithAtom[middle_atom];

        //  Now loop through all the bonds connect to the
        //  middle atom
        while (*bond2 != -1)
        {
          if (bonds[*bond2].atom1 == middle_atom)
          {
            if (i < bonds[*bond2].atom2)
            {
              exclusionSet.add(Exclusion(i,bonds[*bond2].atom2));
            }
          }
          else
          {
            if (i < bonds[*bond2].atom1)
            {
              exclusionSet.add(Exclusion(i,bonds[*bond2].atom1));
            }
          }

          ++bond2;
        }

        ++bond1;
      }
       }
#endif
    }
    /*      END OF FUNCTION build13excl      */

    /************************************************************************/
    /*                  */
    /*        FUNCTION build14excl      */
    /*                  */
    /************************************************************************/


    void Molecule::build14excl(int modified)
       
    {
#ifndef MEM_OPT_VERSION
       int32 *bond1, *bond2, *bond3;  //  The two bonds being checked
       int mid1, mid2;    //  Middle atoms
       register int i;      //  Counter to loop through all atoms
       
       //  Loop through all the atoms
       for (i=0; i<numAtoms; i++)
       {  
      // Get all the bonds connect directly to atom i
      bond1 = bondsWithAtom[i];
       
      while (*bond1 != -1)
      {
        if (bonds[*bond1].atom1 == i)
        {
          mid1=bonds[*bond1].atom2;
        }
        else
        {
          mid1=bonds[*bond1].atom1;
        }

        bond2 = bondsWithAtom[mid1];

        //  Loop through all the bonds connected to atom mid1
        while (*bond2 != -1)
        {
          if (bonds[*bond2].atom1 == mid1)
          {
            mid2 = bonds[*bond2].atom2;
          }
          else
          {
            mid2 = bonds[*bond2].atom1;
          }

          //  Make sure that we don't double back to where
          //  we started from.  This causes strange behavior.
          //  Trust me, I've been there . . .
          if (mid2 == i)
          {
            ++bond2;
            continue;
          }

          bond3=bondsWithAtom[mid2];

          //  Loop through all the bonds connected to mid2
          while (*bond3 != -1)
          {
            if (bonds[*bond3].atom1 == mid2)
            {
              //  Make sure that we don't double back to where
              //  we started from.  This causes strange behavior.
              //  Trust me, I've been there . . .
              //  I added this!!!  Why wasn't it there before?  -JCP
              if (bonds[*bond3].atom2 != mid1)
              if (i<bonds[*bond3].atom2)
              {
                 exclusionSet.add(Exclusion(i,bonds[*bond3].atom2,modified));
              }
            }
            else
            {
              //  Make sure that we don't double back to where
              //  we started from.  This causes strange behavior.
              //  Trust me, I've been there . . .
              //  I added this!!!  Why wasn't it there before?  -JCP
              if (bonds[*bond3].atom1 != mid1)
              if (i<bonds[*bond3].atom1)
              {
                 exclusionSet.add(Exclusion(i,bonds[*bond3].atom1,modified));
              }
            }

            ++bond3;
          }

          ++bond2;
        }
    
        ++bond1;
      }
       }
#endif
    }
    /*      END OF FUNCTION build14excl      */


    /************************************************************************/
    /*                                                                      */
    /*        FUNCTION stripHGroupExcl                                      */
    /*                                                                      */
    /*  This function removes all exclusions which are entirely             */
    /*  within a single hydrogen group.  This assumes that they all         */
    /*  exist, which should be true for exclusion policy 1-3 or higher.     */
    /*                                                                      */
    /************************************************************************/

  void Molecule::stripHGroupExcl(void)
  {
#ifdef NAMD_CUDA
    return;
#endif

#ifdef MEM_OPT_VERSION
   //don't eliminate the exclusion signatures but only adjusting numCalcExclusions
   //which is related to the checksum in reductions. HGroupExcl are not calculated
   //at all
    int delExclCnt=0;
    UniqueSet<Exclusion> strippedExcls;
    HydrogenGroup::iterator h_i, h_e, h_j;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      for ( h_j = h_i + 1; h_j != h_e && ! h_j->isGP; ++h_j ) {
		if ( h_i->atomID < h_j->atomID )
	  		strippedExcls.add(Exclusion(h_i->atomID,h_j->atomID));
		else
	  		strippedExcls.add(Exclusion(h_j->atomID,h_i->atomID));
      }
    }
    
    UniqueSetIter<Exclusion> iter(strippedExcls);
    for(iter=iter.begin(); iter!=iter.end(); iter++){
    	int atom1 = iter->atom1;
    	int atom2 = iter->atom2;
    	ExclusionSignature *a1Sig = &exclSigPool[eachAtomExclSig[atom1]];
    	int fullOrMod, idx;
    	idx = a1Sig->findOffset(atom2-atom1, &fullOrMod);
    	if(idx==-1) continue;
        delExclCnt++;
        if(!exclStrippedByFepOrFixedAtoms(atom1, atom2))
            numCalcExclusions -= 2;
    }

    numTotalExclusions -= delExclCnt;    

    return;

#if 0
	UniqueSet<Exclusion> strippedExcls;
	HydrogenGroup::iterator h_i, h_e, h_j;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      for ( h_j = h_i + 1; h_j != h_e && ! h_j->isGP; ++h_j ) {
		if ( h_i->atomID < h_j->atomID )
	  		strippedExcls.add(Exclusion(h_i->atomID,h_j->atomID));
		else
	  		strippedExcls.add(Exclusion(h_j->atomID,h_i->atomID));
      }
    }
    
    vector<int> *strippedList = new vector<int>[numAtoms];
    UniqueSetIter<Exclusion> iter(strippedExcls);
    for(iter=iter.begin(); iter!=iter.end(); iter++){
    	int atom1 = iter->atom1;
    	int atom2 = iter->atom2;
    	strippedList[atom1].push_back(atom2);
    	strippedList[atom2].push_back(atom1);
    }
    
    vector<ExclusionSignature> newExclSigPool;
    for(int i=0; i<numAtoms; i++){
    	ExclusionSignature a1Sig = exclSigPool[eachAtomExclSig[i]];
    	sort(strippedList[i].begin(), strippedList[i].end());    	
    	//eliminate the atoms in the strippedList from the atom signature    	    	
    	int fullOrMod = 0; //0: full; 1: mod    	
        vector<int> whichOffs;
        vector<int> changedOffs;
    	for(int j=0; j<strippedList[i].size(); j++){
    		int atom2 = strippedList[i][j];
    		int offset = atom2-i;
            //fullOffset and modOffset should be in increasing order
            //for findOffset to be correct. So that assignment has
            //to be done later
    		int idx = a1Sig.findOffset(offset, &fullOrMod);
    		if(idx==-1) continue;
    			//NAMD_die("Something wrong in Molecule::stripHGroupExcl!\n");			
			whichOffs.push_back(fullOrMod);
            changedOffs.push_back(idx);
    	}

        for(int j=0; j<whichOffs.size(); j++){
            if(whichOffs[j]==0){ //full exclusion             
				a1Sig.fullOffset[changedOffs[j]] = 0;
			}else{
				a1Sig.modOffset[changedOffs[j]] = 0;
			}
        }
        whichOffs.clear();
        changedOffs.clear();

		if(strippedList[i].size()){
			//indicating a1Sig has changed
			a1Sig.removeEmptyOffset();
			int poolIndex = lookupCstPool(newExclSigPool, a1Sig);
			if(poolIndex==-1){
				poolIndex = newExclSigPool.size();
				newExclSigPool.push_back(a1Sig);
			}
			eachAtomExclSig[i] = poolIndex + exclSigPoolSize;			
		}		   		    		    		    
    	strippedList[i].clear();
    }
    
    addNewExclSigPool(newExclSigPool);

    newExclSigPool.clear();
           
    delete [] strippedList;     
#endif   
#else
    HydrogenGroup::iterator h_i, h_e, h_j;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      for ( h_j = h_i + 1; h_j != h_e && ! h_j->isGP; ++h_j ) {
	if ( h_i->atomID < h_j->atomID )
	  exclusionSet.del(Exclusion(h_i->atomID,h_j->atomID));
	else
	  exclusionSet.del(Exclusion(h_j->atomID,h_i->atomID));
      }
    }
#endif
  }
    /*      END OF FUNCTION stripHGroupExcl      */

    /************************************************************************/
    /*                                                                      */
    /*        FUNCTION stripFepExcl                                         */
    /*                                                                      */
    /************************************************************************/

#ifdef MEM_OPT_VERSION
  void Molecule::stripFepFixedExcl(void)
  {
      int delExclCnt=0;

      vector<ExclusionSignature> newExclSigPool;
      if(simParams->alchFepOn || simParams->alchThermIntOn || simParams->lesOn){
          for(int i=0; i<numAtoms; i++){
              int atom1 = i;
              int t1 = get_fep_type(atom1);
              if(t1==0) continue;
              ExclusionSignature a1Sig = exclSigPool[eachAtomExclSig[atom1]];
              int sigChanged = 0;
              for(int j=0; j<a1Sig.fullExclCnt; j++){
                  int atom2 = atom1 + a1Sig.fullOffset[j];
                  int t2 = get_fep_type(atom2);
                  if(t2 && t1!=t2){
                      a1Sig.fullOffset[j] = 0;
                      sigChanged = 1;
                      delExclCnt++;
                  }                  
              }
              for(int j=0; j<a1Sig.modExclCnt; j++){
                  int atom2 = atom1 + a1Sig.modOffset[j];
                  int t2 = get_fep_type(atom2);
                  if(t2 && t1!=t2){
                      a1Sig.modOffset[j] = 0;
                      sigChanged = 1;
                      delExclCnt++;
                  }
              }

              //deal with the new signature
              if(sigChanged){
                  a1Sig.removeEmptyOffset();
                  int poolIndex = lookupCstPool(newExclSigPool, a1Sig);
                  if(poolIndex==-1){
                      poolIndex = newExclSigPool.size();
                      newExclSigPool.push_back(a1Sig);
                  }
                  eachAtomExclSig[atom1] = poolIndex + exclSigPoolSize;
              }                            
          }
          //integerate new exclusion signatures into the old one
          addNewExclSigPool(newExclSigPool);
          newExclSigPool.clear();
      }else if(simParams->pairInteractionOn){
          for(int i=0; i<numAtoms; i++){
              int atom1 = i;
              int t1 = get_fep_type(atom1);              
              ExclusionSignature a1Sig = exclSigPool[eachAtomExclSig[atom1]];
              int sigChanged = 0;
              for(int j=0; j<a1Sig.fullExclCnt; j++){
                  int atom2 = atom1 + a1Sig.fullOffset[j];
                  int t2 = get_fep_type(atom2);

                  int cond1 = simParams->pairInteractionSelf && (t1!=1 || t2!=1);
                  int cond2 = (t1!=1 || t2!=2) && (t1!=2 || t2!=1);
                  if(cond1 || cond2){                                   
                      a1Sig.fullOffset[j] = 0;
                      sigChanged = 1;
                      delExclCnt++;
                  }                  
              }
              for(int j=0; j<a1Sig.modExclCnt; j++){
                  int atom2 = atom1 + a1Sig.modOffset[j];
                  int t2 = get_fep_type(atom2);
                  int cond1 = simParams->pairInteractionSelf && (t1!=1 || t2!=1);
                  int cond2 = (t1!=1 || t2!=2) && (t1!=2 || t2!=1);
                  if(cond1 || cond2){                  
                      a1Sig.modOffset[j] = 0;
                      sigChanged = 1;
                      delExclCnt++;
                  }
              }

              //deal with the new signature
              if(sigChanged){
                  a1Sig.removeEmptyOffset();
                  int poolIndex = lookupCstPool(newExclSigPool, a1Sig);
                  if(poolIndex==-1){
                      poolIndex = newExclSigPool.size();
                      newExclSigPool.push_back(a1Sig);
                  }
                  eachAtomExclSig[atom1] = poolIndex + exclSigPoolSize;
              }                            
          }
          //integerate new exclusion signatures into the old one
          addNewExclSigPool(newExclSigPool);
          newExclSigPool.clear();
      }

      numTotalExclusions -= (delExclCnt/2); 

      //deal with fixed atoms
    //numCalcExclusions has to be calculated in stripHGroupExcl
      //numCalcExclusions = 0;
      if(!numFixedAtoms){
    	for(int i=0; i<numAtoms; i++){
    	    ExclusionSignature *a1Sig = &exclSigPool[eachAtomExclSig[i]];
    	    //numCalcExclusions += (a1Sig->fullExclCnt + a1Sig->modExclCnt);
    	}
      }else{
          for(int i=0; i<numAtoms; i++){
              int atom1 = i;
              int t1 = fixedAtomFlags[atom1];
              ExclusionSignature a1Sig = exclSigPool[eachAtomExclSig[atom1]];
              if(t1==0){
                  //numCalcExclusions += (a1Sig.fullExclCnt + a1Sig.modExclCnt);
                  continue;
              }              
              int sigChanged = 0;
              for(int j=0; j<a1Sig.fullExclCnt; j++){
                  int atom2 = atom1 + a1Sig.fullOffset[j];
                  int t2 = fixedAtomFlags[atom2];
                  if(t2){
                      a1Sig.fullOffset[j] = 0;
                      sigChanged = 1;
                      delExclCnt++;
                  }                  
              }
              for(int j=0; j<a1Sig.modExclCnt; j++){
                  int atom2 = atom1 + a1Sig.modOffset[j];
                  int t2 = fixedAtomFlags[atom2];
                  if(t2){
                      a1Sig.modOffset[j] = 0;
                      sigChanged = 1;
                      delExclCnt++;
                  }
              }

              //deal with the new signature
              if(sigChanged){
                  a1Sig.removeEmptyOffset();
                  int poolIndex = lookupCstPool(newExclSigPool, a1Sig);
                  if(poolIndex==-1){
                      poolIndex = newExclSigPool.size();
                      newExclSigPool.push_back(a1Sig);
                  }
                  eachAtomExclSig[atom1] = poolIndex + exclSigPoolSize;
              }              
              //numCalcExclusions += (a1Sig.fullExclCnt + a1Sig.modExclCnt);                            
          }
          //integerate new exclusion signatures into the old one
          addNewExclSigPool(newExclSigPool);
          newExclSigPool.clear();
      }      

      numCalcExclusions = (numCalcExclusions - delExclCnt)/2;      
  }
#else
  void Molecule::stripFepExcl(void)
  {      
    UniqueSet<Exclusion> fepExclusionSet;
    UniqueSetIter<Exclusion> exclIter(exclusionSet);

    if ( simParams->alchFepOn || simParams->alchThermIntOn || simParams->lesOn ) {
       for ( exclIter=exclIter.begin(); exclIter != exclIter.end(); exclIter++ )
       {
         int t1 = get_fep_type(exclIter->atom1);
         int t2 = get_fep_type(exclIter->atom2);
         if ( t1 && t2 && t1 != t2 ) {
           fepExclusionSet.add(*exclIter);
         }
       }
    } else if ( simParams->pairInteractionOn ) {
      for ( exclIter=exclIter.begin(); exclIter != exclIter.end(); exclIter++ )
      {
        int ifep_type = get_fep_type(exclIter->atom1);
        int jfep_type = get_fep_type(exclIter->atom2);
        if ( simParams->pairInteractionSelf ) {
          // for pair-self, both atoms must be in group 1.
          if (ifep_type != 1 || jfep_type != 1) {
            fepExclusionSet.add(*exclIter);
          }
        } else {

          // for pair, must have one from each group.
          if (!(ifep_type == 1 && jfep_type == 2) &&
              !(ifep_type == 2 && jfep_type == 1)) {
            fepExclusionSet.add(*exclIter);
          }
        }
       }
    }

    UniqueSetIter<Exclusion> fepIter(fepExclusionSet);
    for ( fepIter=fepIter.begin(); fepIter != fepIter.end(); fepIter++ )
    {
      exclusionSet.del(*fepIter);
    }
  }
#endif
    /*      END OF FUNCTION stripFepExcl      */


/* BEGIN gf */
    /************************************************************************/
    /*                                                                      */
    /*      FUNCTION build_gridforce_params                                 */
    /*                                                                      */
    /*   INPUTS:                                                            */
    /*  gridfrcfile - Value of gridforcefile from config file               */
    /*  gridfrccol - Value of gridforcecol from config file                 */
    /*  gridfrcchrgcol - Value of gridforcechargecol from config file	    */
    /*  potfile - Value of gridforcepotfile from config file                  */
    /*  initial_pdb - PDB object that contains initial positions            */
    /*  cwd - Current working directory                                     */
    /*                                                                      */
    // This function builds all the parameters that are necessary to
    // do gridforcing. This involves looking through a PDB object to
    // determine which atoms are to be gridforced, and what the force
    // multiplier is for each atom.  This information is then stored
    // in the arrays gridfrcIndexes and gridfrcParams.
    /************************************************************************/

void Molecule::build_gridforce_params(StringList *gridfrcfile,
				      StringList *gridfrccol,
				      StringList *gridfrcchrgcol,
				      StringList *potfile,
				      PDB *initial_pdb,
				      char *cwd)
{
    PDB *kPDB;
    register int i;		//  Loop counters
    register int j;
    register int k;

    DebugM(3,  "Entered build_gridforce_params multi...\n");
//     DebugM(3, "\tgridfrcfile = " << gridfrcfile->data << endi);
//     DebugM(3, "\tgridfrccol = " << gridfrccol->data << endi);
    
    MGridforceParams* mgridParams = simParams->mgridforcelist.get_first();
    numGridforceGrids = 0;
    while (mgridParams != NULL) {
	numGridforceGrids++;
	mgridParams = mgridParams->next;
    }
    
    gridfrcIndexes = new int32*[numGridforceGrids];
    gridfrcParams = new GridforceParams*[numGridforceGrids];
    gridfrcGrid = new GridforceGrid*[numGridforceGrids];
    numGridforces = new int[numGridforceGrids];
    
    mgridParams = simParams->mgridforcelist.get_first();
    for (int gridnum = 0; gridnum < numGridforceGrids; gridnum++) {
	int current_index=0;	//  Index into values used
	int kcol = 5;		//  Column to look for force constant in
	int qcol = 0;		//  Column for charge (default 0: use electric charge)
	Real kval = 0;		//  Force constant value retreived
	char filename[129];	//  PDB filename
	char potfilename[129];	//  Potential file name
	
	if (mgridParams == NULL) {
	    NAMD_die("Problem with mgridParams!");
	}
	
	// Now load values from mgridforcelist object
	if (mgridParams->gridforceFile == NULL)
	{
	    kPDB = initial_pdb;
	}
	else
	{
	    DebugM(4, "mgridParams->gridforceFile = " << mgridParams->gridforceFile << "\n" << endi);
	    
	    if ( (cwd == NULL) || (mgridParams->gridforceFile[0] == '/') )
	    {
		strcpy(filename, mgridParams->gridforceFile);
	    }
	    else
	    {
		strcpy(filename, cwd);
		strcat(filename, mgridParams->gridforceFile);
	    }
	
	    kPDB = new PDB(filename);
	    if ( kPDB == NULL )
	    {
		NAMD_die("Memory allocation failed in Molecule::build_gridforce_params");
	    }
	   
	    if (kPDB->num_atoms() != numAtoms)
	    {
		NAMD_die("Number of atoms in grid force PDB doesn't match coordinate PDB");
	    }
	}

	//  Get the column that the force constant is going to be in.  It
	//  can be in any of the 5 floating point fields in the PDB, according
	//  to what the user wants.  The allowable fields are X, Y, Z, O, or
	//  B which correspond to the 1st, 2nd, ... 5th floating point fields.
	//  The default is the 5th field, which is beta (temperature factor)
	if (mgridParams->gridforceCol == NULL)
	{
	    kcol = 5;
	}
	else
	{
	    if (strcasecmp(mgridParams->gridforceCol, "X") == 0)
	    {
		kcol=1;
	    }
	    else if (strcasecmp(mgridParams->gridforceCol, "Y") == 0)
	    {
		kcol=2;
	    }
	    else if (strcasecmp(mgridParams->gridforceCol, "Z") == 0)
	    {
		kcol=3;
	    }
	    else if (strcasecmp(mgridParams->gridforceCol, "O") == 0)
	    {
		kcol=4;
	    }
	    else if (strcasecmp(mgridParams->gridforceCol, "B") == 0)
	    {
		kcol=5;
	    }
	    else
	    {
		NAMD_die("gridforcecol must have value of X, Y, Z, O, or B");
	    }
	}
    
	//  Get the column that the charge is going to be in.
        if (mgridParams->gridforceQcol == NULL)
	{
	    qcol = 0;	// Default: don't read charge from file, use electric charge
	}
	else
	{
	    if (strcasecmp(mgridParams->gridforceQcol, "X") == 0)
	    {
		qcol=1;
	    }
	    else if (strcasecmp(mgridParams->gridforceQcol, "Y") == 0)
	    {
		qcol=2;
	    }
	    else if (strcasecmp(mgridParams->gridforceQcol, "Z") == 0)
	    {
		qcol=3;
	    }
	    else if (strcasecmp(mgridParams->gridforceQcol, "O") == 0)
	    {
		qcol=4;
	    }
	    else if (strcasecmp(mgridParams->gridforceQcol, "B") == 0)
	    {
		qcol=5;
	    }
	    else
	    {
		NAMD_die("gridforcechargecol must have value of X, Y, Z, O, or B");
	    }
	}
    
	if (kcol == qcol) {
	    NAMD_die("gridforcecol and gridforcechargecol cannot have same value");
	}

    
	//  Allocate an array that will store an index into the constraint
	//  parameters for each atom.  If the atom is not constrained, its
	//  value will be set to -1 in this array.
	DebugM(4, "Here!\n" << endi);
	gridfrcIndexes[gridnum] = new int32[numAtoms];
	DebugM(4, "There!\n" << endi);
       
	if (gridfrcIndexes[gridnum] == NULL)
	{
	    NAMD_die("memory allocation failed in Molecule::build_gridforce_params()");
	}

	//  Loop through all the atoms and find out which ones are constrained
	for (i=0; i<numAtoms; i++)
	{
	    //  Get the k value based on where we were told to find it
	    switch (kcol)
	    {
	    case 1:
		kval = (kPDB->atom(i))->xcoor();
		break;
	    case 2:
		kval = (kPDB->atom(i))->ycoor();
		break;
	    case 3:
		kval = (kPDB->atom(i))->zcoor();
		break;
	    case 4:
		kval = (kPDB->atom(i))->occupancy();
		break;
	    case 5:
		kval = (kPDB->atom(i))->temperaturefactor();
		break;
	    }
	   
	    if (kval > 0.0)
	    {
		//  This atom is constrained
		gridfrcIndexes[gridnum][i] = current_index;
		current_index++;
	    }
	    else
	    {
		//  This atom is not constrained
		gridfrcIndexes[gridnum][i] = -1;
	    }
	}
    
	if (current_index == 0)
	{
	    //  Constraints were turned on, but there weren't really any constrained
	    iout << iWARN << "NO GRIDFORCE ATOMS WERE FOUND, BUT GRIDFORCE IS ON . . .\n" << endi;
	}
	else
	{
	    //  Allocate an array to hold the constraint parameters
	    gridfrcParams[gridnum] = new GridforceParams[current_index];
	    if (gridfrcParams[gridnum] == NULL)
	    {
		NAMD_die("memory allocation failed in Molecule::build_gridforce_params");
	    }
	}
    
	numGridforces[gridnum] = current_index;

	//  Loop through all the atoms and assign the parameters for those
	//  that are constrained
	for (i=0; i<numAtoms; i++)
	{
	    if (gridfrcIndexes[gridnum][i] != -1)
	    {
		//  This atom has grid force, so get the k value again
		switch (kcol)
		{
		case 1:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->xcoor();
		    break;
		case 2:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->ycoor();
		    break;
		case 3:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->zcoor();
		    break;
		case 4:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->occupancy();
		    break;
		case 5:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].k = (kPDB->atom(i))->temperaturefactor();
		    break;
		}
	    
		//  Also get charge column
		switch (qcol)
		{
		case 0:
#ifdef MEM_OPT_VERSION
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = atomChargePool[eachAtomCharge[i]];
#else
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = atoms[i].charge;
#endif
		    break;
		case 1:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->xcoor();
		    break;
		case 2:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->ycoor();
		    break;
		case 3:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->zcoor();
		    break;
		case 4:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->occupancy();
		    break;
		case 5:
		    gridfrcParams[gridnum][gridfrcIndexes[gridnum][i]].q = (kPDB->atom(i))->temperaturefactor();
		    break;
		}
	    }
	}
       
	//  If we had to create new PDB objects, delete them now
	if (mgridParams->gridforceFile != NULL)
	{
	    delete kPDB;
	}
    
	//  Now we fill in our grid information
    
	// Open potential file
	if ( (cwd == NULL) || (mgridParams->gridforceVfile[0] == '/') )
	{
	    strcpy(potfilename, mgridParams->gridforceVfile);
	}
	else
	{
	    strcpy(potfilename, cwd);
	    strcat(potfilename, mgridParams->gridforceVfile);
	}
    
//        iout << iINFO << "Allocating grid " << gridnum
//             << "\n" << endi;
             
	gridfrcGrid[gridnum] = new GridforceGrid(gridnum);
	//	gridfrcGrid[gridnum]->initialize(potfilename, simParams, mgridParams);
	gridfrcGrid[gridnum]->init1(potfilename, simParams, mgridParams);
//        ComputeGridForceNodeMgr* mgr = CProxy_ComputeGridForceNodeMgr::
//          ckLocalBranch(CkpvAccess(BOCclass_group).
//          computeGridForceNodeMgr);
//      mgr->depositInitialGrid(gridnum,gridfrcGrid[gridnum],numGridforceGrids);

        CProxy_ComputeGridForceNodeMgr 
          mgr(CkpvAccess(BOCclass_group).computeGridForceNodeMgr);
          
        GridDepositMsg *outmsg = new GridDepositMsg;
        outmsg->gridnum = gridnum;
        outmsg->grid = gridfrcGrid[gridnum];
        outmsg->num_grids = numGridforceGrids;
        mgr[CkMyNode()].depositInitialGrid(outmsg);        
 
	// Finally, get next mgridParams pointer
	mgridParams = mgridParams->next;
    }
}
/* END gf */


#endif  // MOLECULE2_C undefined = first object file
#ifdef MOLECULE2_C  // second object file


    /************************************************************************/
    /*                  */
    /*      FUNCTION build_constraint_params    */
    /*                  */
    /*   INPUTS:                */
    /*  consref - Value of consref parameter from config file    */
    /*  conskfile - Value of conskfile from config file      */
    /*  conskcol - Value of conskcol from config file      */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*  cwd - Current working directory          */
    /*                  */
    /*  This function builds all the parameters that are necessary  */
    /*   to do harmonic constraints.  This involves looking through    */
    /*   one or more PDB objects to determine which atoms are constrained,  */
    /*   and what the force constant and reference position is force each   */
    /*   atom that is constrained.  This information is then stored    */
    /*   in the arrays consIndexes and consParams.        */
    /*                  */
    /************************************************************************/

    void Molecule::build_constraint_params(StringList *consref, 
             StringList *conskfile, 
             StringList *conskcol, 
             PDB *initial_pdb,
             char *cwd)
       
    {
       PDB *refPDB, *kPDB;    //  Pointer to other PDB's if used
       register int i;      //  Loop counter
       int current_index=0;    //  Index into values used
       int kcol = 4;      //  Column to look for force constant in
       Real kval = 0;      //  Force constant value retreived
       char filename[129];    //  PDB filename
       
       //  Get the PDB object that contains the reference positions.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates.  i.e., constrain
       //  the atoms around their initial position.  This is the most likely
       //  case anyway
       if (consref == NULL)
       {
    refPDB = initial_pdb;
       }
       else
       {
    if (consref->next != NULL)
    {
       NAMD_die("Multiple definitions of constraint reference file in configruation file");
    }

    if ( (cwd == NULL) || (consref->data[0] == '/') )
    {
         strcpy(filename, consref->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, consref->data);
    }
    
    refPDB = new PDB(filename);
    if ( refPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_constraint_params");
    }
    
    if (refPDB->num_atoms() != numAtoms)
    {
       NAMD_die("Number of atoms in constraint reference PDB doesn't match coordinate PDB");
    }
       }
       
       //  Get the PDB to read the force constants from.  Again, if the user
       //  gave us another file name, open that one.  Otherwise, just use
       //  the PDB with the initial coordinates
       if (conskfile == NULL)
       {
    kPDB = initial_pdb;
       }
       else
       {
    if (conskfile->next != NULL)
    {
       NAMD_die("Multiple definitions of constraint constant file in configuration file");
    }

    if ( (consref != NULL) && (strcasecmp(consref->data, conskfile->data) == 0) )
    {
       //  Same PDB used for reference positions and force constants
       kPDB = refPDB; 
    }
    else
    {
      if ( (cwd == NULL) || (conskfile->data[0] == '/') )
      {
        strcpy(filename, conskfile->data);
      }
      else
      {
        strcpy(filename, cwd);
        strcat(filename, conskfile->data);
      }

      kPDB = new PDB(filename);
      if ( kPDB == NULL )
      {
        NAMD_die("Memory allocation failed in Molecule::build_constraint_params");
      }
    
      if (kPDB->num_atoms() != numAtoms)
      {
         NAMD_die("Number of atoms in constraint constant PDB doesn't match coordinate PDB");
      }
    }
       }
       
       //  Get the column that the force constant is going to be in.  It
       //  can be in any of the 5 floating point fields in the PDB, according
       //  to what the user wants.  The allowable fields are X, Y, Z, O, or
       //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
       //  The default is the 4th field, which is the occupancy
       if (conskcol == NULL)
       {
    kcol = 4;
       }
       else
       {
    if (conskcol->next != NULL)
    {
       NAMD_die("Multiple definitions of harmonic constraint column in config file");
    }
    
    if (strcasecmp(conskcol->data, "X") == 0)
    {
       kcol=1;
    }
    else if (strcasecmp(conskcol->data, "Y") == 0)
    {
       kcol=2;
    }
    else if (strcasecmp(conskcol->data, "Z") == 0)
    {
       kcol=3;
    }
    else if (strcasecmp(conskcol->data, "O") == 0)
    {
       kcol=4;
    }
    else if (strcasecmp(conskcol->data, "B") == 0)
    {
       kcol=5;
    }
    else
    {
       NAMD_die("conskcol must have value of X, Y, Z, O, or B");
    }
       }
       
       //  Allocate an array that will store an index into the constraint
       //  parameters for each atom.  If the atom is not constrained, its
       //  value will be set to -1 in this array.
       consIndexes = new int32[numAtoms];
       
       if (consIndexes == NULL)
       {
    NAMD_die("memory allocation failed in Molecule::build_constraint_params()");
       }
       
       //  Loop through all the atoms and find out which ones are constrained
       for (i=0; i<numAtoms; i++)
       {
    //  Get the k value based on where we were told to find it
    switch (kcol)
    {
       case 1:
    kval = (kPDB->atom(i))->xcoor();
    break;
       case 2:
    kval = (kPDB->atom(i))->ycoor();
    break;
       case 3:
    kval = (kPDB->atom(i))->zcoor();
    break;
       case 4:
    kval = (kPDB->atom(i))->occupancy();
    break;
       case 5:
    kval = (kPDB->atom(i))->temperaturefactor();
    break;
    }
    
    if (kval > 0.0)
    {
       //  This atom is constrained
       consIndexes[i] = current_index;
       current_index++;
    }
    else
    {
       //  This atom is not constrained
       consIndexes[i] = -1;
    }
       }
       
       if (current_index == 0)
       {
    //  Constraints were turned on, but there weren't really any constrained
    iout << iWARN << "NO CONSTRAINED ATOMS WERE FOUND, BUT CONSTRAINTS ARE ON . . .\n" << endi;
       }
       else
       {
    //  Allocate an array to hold the constraint parameters
    consParams = new ConstraintParams[current_index];
    
    if (consParams == NULL)
    {
       NAMD_die("memory allocation failed in Molecule::build_constraint_params");
    }
       }
       
       numConstraints = current_index;
       
       //  Loop through all the atoms and assign the parameters for those
       //  that are constrained
       for (i=0; i<numAtoms; i++)
       {
    if (consIndexes[i] != -1)
    {
       //  This atom is constrained, so get the k value again
       switch (kcol)
       {
          case 1:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->xcoor();
       break;
          case 2:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->ycoor();
       break;
          case 3:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->zcoor();
       break;
          case 4:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->occupancy();
       break;
          case 5:
       consParams[consIndexes[i]].k = (kPDB->atom(i))->temperaturefactor();
       break;
       }
       
       //  Get the reference position
       consParams[consIndexes[i]].refPos.x = (refPDB->atom(i))->xcoor();
       consParams[consIndexes[i]].refPos.y = (refPDB->atom(i))->ycoor();
       consParams[consIndexes[i]].refPos.z = (refPDB->atom(i))->zcoor();
    }
       }
       
       //  If we had to create new PDB objects, delete them now
       if (consref != NULL)
       {
    delete refPDB;
       }
       
       if ((conskfile != NULL) &&
     !((consref != NULL) && 
       (strcasecmp(consref->data, conskfile->data) == 0)
      )
    )
       {
    delete kPDB;
       }

    }
    /*      END OF FUNCTION build_constraint_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_movdrag_params  */
/*                  */
/*   INPUTS:        */
/*  movDragFile - value of movDragFile from the config file */
/*  movDragCol - value of movDragCol from the config file */
/*  movDragVelFile - value of movDragVelFile from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do moving drag. This involves looking through one or more    */
/*  PDB objects to determine which atoms are dragged,  and what the */
/*  drag parameters for each atom are. This information is then stored */
/*  in the arrays movDragIndexes and movDragParams. */
/*                  */
/************************************************************************/

void Molecule::build_movdrag_params(StringList *movDragFile, 
				    StringList *movDragCol, 
				    StringList *movDragVelFile, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *vPDB;        //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for drag tag in
  Real dtval = 0;          //  Drag tag value retreived
  char mainfilename[129];  //  main moving drag PDB filename
  char velfilename[129];   //  moving drag velocity PDB filename
  
  //  Get the PDB to read the moving drag tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (movDragFile == NULL) {
    tPDB = initial_pdb;
    
  } else {

    if (movDragFile->next != NULL) {
      NAMD_die("Multiple definitions of moving drag tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (movDragFile->data[0] == '/') ) {
      strcpy(mainfilename, movDragFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, movDragFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_movdrag_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in moving drag tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom velocities. If no name given, use
  // movDragFile if it is defined. Can NOT use the PDB coordinate
  // file!
  
  if (movDragVelFile == NULL) {
    if (movDragFile == NULL) {
      NAMD_die("Moving drag velocity file can not be same as coordinate PDB file");
    } else {
      if (movDragVelFile->next != NULL) {
	NAMD_die("Multiple definitions of moving drag velocity file in configuration file");
      };
      vPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (movDragVelFile->data[0] == '/') ) {
      strcpy(velfilename, movDragVelFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, movDragVelFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_movdrag_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in moving drag velocity PDB doesn't match coordinate PDB");
    }
  };
  
  
  //  Get the column that the drag tag is going to be in. If
  //  movDragFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If movDragFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (movDragCol == NULL) {
    dtcol = 4;
  } else {
    if (movDragCol->next != NULL) {
      NAMD_die("Multiple definitions of drag column in config file");
    };
    
    if (movDragFile == NULL
	&& strcasecmp(movDragCol->data, "B")
	&& strcasecmp(movDragCol->data, "O")) {
      NAMD_die("Can not read moving drag tags from X, Y, or Z column of the coordinate or velocity file");
    };
    if (!strcasecmp(movDragCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(movDragCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(movDragCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(movDragCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(movDragCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("movDragCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the drag
  //  parameters for each atom.  If the atom is not dragged, its
  //  value will be set to -1 in this array.
  movDragIndexes = new int32[numAtoms];
    if (movDragIndexes == NULL) {
    NAMD_die("memory allocation failed in Molecule::build_movdrag_params()");
  };
  
  //  Loop through all the atoms and find out which ones are dragged
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is dragged
      movDragIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not dragged
      movDragIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    //  Drag was turned on, but there weren't really any dragged
    iout << iWARN << "NO DRAGGED ATOMS WERE FOUND, BUT MOVING DRAG IS ON . . . " << endi;
  } else {
    //  Allocate an array to hold the drag parameters
    movDragParams = new MovDragParams[current_index];
    if (movDragParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_movdrag_params");
    }
  };
  
  numMovDrag = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are dragged
  for (i=0; i<numAtoms; i++) {
    if (movDragIndexes[i] != -1) {
      movDragParams[movDragIndexes[i]].v[0] = (vPDB->atom(i))->xcoor();
      movDragParams[movDragIndexes[i]].v[1] = (vPDB->atom(i))->ycoor();
      movDragParams[movDragIndexes[i]].v[2] = (vPDB->atom(i))->zcoor();
    };
  };
      
  if (movDragFile != NULL) delete tPDB;
  if (movDragVelFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_movdrag_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_rotdrag_params  */
/*                  */
/*   INPUTS:        */
/*  rotDragFile - value of rotDragFile from the config file */
/*  rotDragCol - value of rotDragCol from the config file */
/*  rotDragAxisFile - value of rotDragAxisFile from the config file */
/*  rotDragPivotFile - value of rotDragPivotFile from the config file */
/*  rotDragVelFile - value of rotDragVelFile from the config file */
/*  rotDragVelCol - value of rotDragVelCol from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do moving drag. This involves looking through one or more    */
/*  PDB objects to determine which atoms are dragged,  and what the */
/*  drag parameters for each atom are. This information is then stored */
/*  in the arrays rotDragIndexes and rotDragParams. */
/*                  */
/************************************************************************/

void Molecule::build_rotdrag_params(StringList *rotDragFile, 
				    StringList *rotDragCol, 
				    StringList *rotDragAxisFile, 
				    StringList *rotDragPivotFile, 
				    StringList *rotDragVelFile, 
				    StringList *rotDragVelCol, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *aPDB, *pPDB, *vPDB; //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for drag tag in
  Real dtval = 0;          //  Drag tag value retreived
  int dvcol = 4;           //  Column to look for angular velocity in
  Real dvval = 0;          //  Angular velocity value retreived
  char mainfilename[129];  //  main rotating drag PDB filename
  char axisfilename[129];  //  rotating drag axis PDB filename
  char pivotfilename[129]; //  rotating drag pivot point PDB filename
  char velfilename[129];   //  rotating drag angular velocity PDB filename
  
  //  Get the PDB to read the rotating drag tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (rotDragFile == NULL) {
    tPDB = initial_pdb;
    
  } else {

    if (rotDragFile->next != NULL) {
      NAMD_die("Multiple definitions of rotating drag tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (rotDragFile->data[0] == '/') ) {
      strcpy(mainfilename, rotDragFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, rotDragFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom rotation axes. If no name given, use
  // rotDragFile if both it AND rotDragPivotFile are defined. Can NOT
  // use the PDB coordinate file, nor rotDragPivotFile!

  if (rotDragAxisFile == NULL) {
    if (rotDragFile == NULL) {
      NAMD_die("Rotating drag axis file can not be same as coordinate PDB file");
    } else {
      if (rotDragAxisFile->next != NULL) {
	NAMD_die("Multiple definitions of rotating drag axis file in configuration file");
      };
      if (rotDragPivotFile == NULL) {
	NAMD_die("Need to specify at least one of rotDragAxisFile and rotDragPivotFile; they can not be same");
      };
      aPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (rotDragAxisFile->data[0] == '/') ) {
      strcpy(axisfilename, rotDragAxisFile->data);
    } else {
      strcpy(axisfilename, cwd);
      strcat(axisfilename, rotDragAxisFile->data);
    }
    
    aPDB = new PDB(axisfilename);
    if ( aPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (aPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag axis PDB doesn't match coordinate PDB");
    }
  };
  
  // Get the PDB to read atom rotation pivot points. If no name given,
  // use rotDragFile if both it AND rotDragAxisFile are defined. Can
  // NOT use the PDB coordinate file, nor rotDragAxisFile!

  if (rotDragPivotFile == NULL) {
    if (rotDragFile == NULL) {
      NAMD_die("Rotating drag pivot point file can not be same as coordinate PDB file");
    } else {
      if (rotDragPivotFile->next != NULL) {
	NAMD_die("Multiple definitions of rotating drag pivot point file in configuration file");
      };
      if (rotDragAxisFile == NULL) {
	NAMD_die("Need to specify at least one of rotDragAxisFile and rotDragPivotFile; they can not be same");
      };
      pPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (rotDragPivotFile->data[0] == '/') ) {
      strcpy(pivotfilename, rotDragPivotFile->data);
    } else {
      strcpy(pivotfilename, cwd);
      strcat(pivotfilename, rotDragPivotFile->data);
    }
    
    pPDB = new PDB(pivotfilename);
    if ( pPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (pPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag pivot point PDB doesn't match coordinate PDB");
    }
  };
  
  
  // Get the PDB to read atom angular velocities. If no name given,
  // use rotDragFile (or the coordinate PDB file if rotDragFile is not
  // defined).

  if (rotDragVelFile == NULL) {
    vPDB = tPDB;
  } else {
    if (rotDragVelFile->next != NULL) {
      NAMD_die("Multiple definitions of rotating drag velocity file in configuration file");
    };
    
    if ( (cwd == NULL) || (rotDragVelFile->data[0] == '/') ) {
      strcpy(velfilename, rotDragVelFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, rotDragVelFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_rotdrag_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in rotating drag velocity PDB doesn't match coordinate PDB");
    }
  };
  
  //  Get the column that the drag tag is going to be in. If
  //  rotDragFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If rotDragFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (rotDragCol == NULL) {
    dtcol = 4;
  } else {
    if (rotDragCol->next != NULL) {
      NAMD_die("Multiple definitions of drag tag column in config file");
    };
    
    if ( rotDragFile == NULL
	 && (!strcasecmp(rotDragCol->data, "X")
	     || !strcasecmp(rotDragCol->data, "Y")
	     || !strcasecmp(rotDragCol->data, "Z"))) {
      NAMD_die("Can not read rotating drag tags from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(rotDragCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(rotDragCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(rotDragCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(rotDragCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(rotDragCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("rotDragCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Get the column that the drag angular velocity is going to be
  //  in. If rotDragVelFile is defined, it can be in any of the 5
  //  floating point fields in the PDB (X, Y, Z, O, or B) which
  //  correspond to the 1st, 2nd, ... 5th floating point fields. If
  //  NEITHER of rotDragVelFile OR rotDragFile is defined, it can
  //  only be O or B fileds. The default is the O (4th) field, which
  //  is the occupancy.

  if (rotDragVelCol == NULL) {
    dvcol = 4;
  } else {
    if (rotDragVelCol->next != NULL) {
      NAMD_die("Multiple definitions of drag angular velocity column in config file");
    };
    
    if (rotDragVelFile == NULL
	&& rotDragFile == NULL
	&& strcasecmp(rotDragCol->data, "B")
	&& strcasecmp(rotDragCol->data, "O")) {
      NAMD_die("Can not read rotating drag angular velocities from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(rotDragVelCol->data, "X")) {
      dvcol=1;
    } else if (!strcasecmp(rotDragVelCol->data, "Y")) {
      dvcol=2;
    } else if (!strcasecmp(rotDragVelCol->data, "Z")) {
      dvcol=3;
    } else if (!strcasecmp(rotDragVelCol->data, "O")) {
      dvcol=4;
    } else if (!strcasecmp(rotDragVelCol->data, "B")) {
      dvcol=5;
    }
    else {
      NAMD_die("rotDragVelCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the drag
  //  parameters for each atom.  If the atom is not dragged, its
  //  value will be set to -1 in this array.
  rotDragIndexes = new int32[numAtoms];
  if (rotDragIndexes == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_rotdrag_params()");
  };
  
  //  Loop through all the atoms and find out which ones are dragged
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is dragged
      rotDragIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not dragged
      rotDragIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    iout << iWARN << "NO DRAGGED ATOMS WERE FOUND, BUT ROTATING DRAG IS ON . . . " << endi;
  } else {
    rotDragParams = new RotDragParams[current_index];
    if (rotDragParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_rotdrag_params");
    }
  };
  
  numRotDrag = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are dragged
  for (i=0; i<numAtoms; i++) {
    if (rotDragIndexes[i] != -1) {
      rotDragParams[rotDragIndexes[i]].a[0] = (aPDB->atom(i))->xcoor();
      rotDragParams[rotDragIndexes[i]].a[1] = (aPDB->atom(i))->ycoor();
      rotDragParams[rotDragIndexes[i]].a[2] = (aPDB->atom(i))->zcoor();
      rotDragParams[rotDragIndexes[i]].p[0] = (pPDB->atom(i))->xcoor();
      rotDragParams[rotDragIndexes[i]].p[1] = (pPDB->atom(i))->ycoor();
      rotDragParams[rotDragIndexes[i]].p[2] = (pPDB->atom(i))->zcoor();
      switch (dvcol) {
      case 1:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->xcoor();
	break;
      case 2:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->ycoor();
	break;
      case 3:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->zcoor();
	break;
      case 4:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->occupancy();
	break;
      case 5:
	rotDragParams[rotDragIndexes[i]].v = (vPDB->atom(i))->temperaturefactor();
	break;
      };
    };
  };
      
  if (rotDragFile != NULL) delete tPDB;
  if (rotDragAxisFile != NULL) delete aPDB;
  if (rotDragPivotFile != NULL) delete pPDB;
  if (rotDragVelFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_rotdrag_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_constorque_params  */
/*                  */
/*   INPUTS:        */
/*  consTorqueFile - value of consTorqueFile from the config file */
/*  consTorqueCol - value of consTorqueCol from the config file */
/*  consTorqueAxisFile - value of consTorqueAxisFile from the config file */
/*  consTorquePivotFile - value of consTorquePivotFile from the config file */
/*  consTorqueValFile - value of consTorqueValFile from the config file */
/*  consTorqueValCol - value of consTorqueValCol from the config file */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory          */
/*                  */
/*  This function builds all the parameters that are necessary  */
/*  to do "constant" torque. This involves looking through one or more    */
/*  PDB objects to determine which atoms are torqued,  and what the */
/*  torque parameters for each atom are. This information is then stored */
/*  in the arrays consTorqueIndexes and consTorqueParams. */
/*                  */
/************************************************************************/

void Molecule::build_constorque_params(StringList *consTorqueFile, 
				    StringList *consTorqueCol, 
				    StringList *consTorqueAxisFile, 
				    StringList *consTorquePivotFile, 
				    StringList *consTorqueValFile, 
				    StringList *consTorqueValCol, 
				    PDB *initial_pdb,
				    char *cwd)
  
{
  PDB *tPDB, *aPDB, *pPDB, *vPDB; //  Pointers to other PDB file(s)
  register int i;          //  Loop counter
  int current_index=0;     //  Index into values used
  int dtcol = 4;           //  Column to look for torque tag in
  Real dtval = 0;          //  Torque tag value retreived
  int dvcol = 4;           //  Column to look for angular velocity in
  Real dvval = 0;          //  Angular velocity value retreived
  char mainfilename[129];  //  main "constant" torque PDB filename
  char axisfilename[129];  //  "constant" torque axis PDB filename
  char pivotfilename[129]; //  "constant" torque pivot point PDB filename
  char velfilename[129];   //  "constant" torque angular velocity PDB filename
  
  //  Get the PDB to read the "constant" torque tags from. Again, if the
  //  user gave us another file name, open that one.  Otherwise, just
  //  use the PDB with the initial coordinates
  if (consTorqueFile == NULL) {
    tPDB = initial_pdb;
    
  } else {

    if (consTorqueFile->next != NULL) {
      NAMD_die("Multiple definitions of \"constant\" torque tag file in configuration file");
    }
    
    if ( (cwd == NULL) || (consTorqueFile->data[0] == '/') ) {
      strcpy(mainfilename, consTorqueFile->data);
    } else {
      strcpy(mainfilename, cwd);
      strcat(mainfilename, consTorqueFile->data);
      }
    
    tPDB = new PDB(mainfilename);
    if ( tPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (tPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque tag PDB doesn't match coordinate PDB");
    }
  }
  
  // Get the PDB to read atom rotation axes. If no name given, use
  // consTorqueFile if both it AND consTorquePivotFile are defined. Can NOT
  // use the PDB coordinate file, nor consTorquePivotFile!

  if (consTorqueAxisFile == NULL) {
    if (consTorqueFile == NULL) {
      NAMD_die("\"Constant\" torque axis file can not be same as coordinate PDB file");
    } else {
      if (consTorqueAxisFile->next != NULL) {
	NAMD_die("Multiple definitions of \"constant\" torque axis file in configuration file");
      };
      if (consTorquePivotFile == NULL) {
	NAMD_die("Need to specify at least one of consTorqueAxisFile and consTorquePivotFile; they can not be same");
      };
      aPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (consTorqueAxisFile->data[0] == '/') ) {
      strcpy(axisfilename, consTorqueAxisFile->data);
    } else {
      strcpy(axisfilename, cwd);
      strcat(axisfilename, consTorqueAxisFile->data);
    }
    
    aPDB = new PDB(axisfilename);
    if ( aPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (aPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque axis PDB doesn't match coordinate PDB");
    }
  };
  
  // Get the PDB to read atom rotation pivot points. If no name given,
  // use consTorqueFile if both it AND consTorqueAxisFile are defined. Can
  // NOT use the PDB coordinate file, nor consTorqueAxisFile!

  if (consTorquePivotFile == NULL) {
    if (consTorqueFile == NULL) {
      NAMD_die("\"Constant\" torque pivot point file can not be same as coordinate PDB file");
    } else {
      if (consTorquePivotFile->next != NULL) {
	NAMD_die("Multiple definitions of \"constant\" torque pivot point file in configuration file");
      };
      if (consTorqueAxisFile == NULL) {
	NAMD_die("Need to specify at least one of consTorqueAxisFile and consTorquePivotFile; they can not be same");
      };
      pPDB = tPDB;
    };

  } else {

    if ( (cwd == NULL) || (consTorquePivotFile->data[0] == '/') ) {
      strcpy(pivotfilename, consTorquePivotFile->data);
    } else {
      strcpy(pivotfilename, cwd);
      strcat(pivotfilename, consTorquePivotFile->data);
    }
    
    pPDB = new PDB(pivotfilename);
    if ( pPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (pPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque pivot point PDB doesn't match coordinate PDB");
    }
  };
  
  
  // Get the PDB to read atom angular velocities. If no name given,
  // use consTorqueFile (or the coordinate PDB file if consTorqueFile is not
  // defined).

  if (consTorqueValFile == NULL) {
    vPDB = tPDB;
  } else {
    if (consTorqueValFile->next != NULL) {
      NAMD_die("Multiple definitions of \"constant\" torque velocity file in configuration file");
    };
    
    if ( (cwd == NULL) || (consTorqueValFile->data[0] == '/') ) {
      strcpy(velfilename, consTorqueValFile->data);
    } else {
      strcpy(velfilename, cwd);
      strcat(velfilename, consTorqueValFile->data);
    }
    
    vPDB = new PDB(velfilename);
    if ( vPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_constorque_params");
    }
    
    if (vPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in \"constant\" torque velocity PDB doesn't match coordinate PDB");
    }
  };
  
  //  Get the column that the torque tag is going to be in. If
  //  consTorqueFile is defined, it can be in any of the 5 floating point
  //  fields in the PDB (X, Y, Z, O, or B) which correspond to the
  //  1st, 2nd, ... 5th floating point fields. If consTorqueFile is NOT
  //  defined, it can only be O or B fileds. The default is the O
  //  (4th) field, which is the occupancy.

  if (consTorqueCol == NULL) {
    dtcol = 4;
  } else {
    if (consTorqueCol->next != NULL) {
      NAMD_die("Multiple definitions of torque tag column in config file");
    };
    
    if ( consTorqueFile == NULL
	 && (!strcasecmp(consTorqueCol->data, "X")
	     || !strcasecmp(consTorqueCol->data, "Y")
	     || !strcasecmp(consTorqueCol->data, "Z"))) {
      NAMD_die("Can not read \"constant\" torque tags from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(consTorqueCol->data, "X")) {
      dtcol=1;
    } else if (!strcasecmp(consTorqueCol->data, "Y")) {
      dtcol=2;
    } else if (!strcasecmp(consTorqueCol->data, "Z")) {
      dtcol=3;
    } else if (!strcasecmp(consTorqueCol->data, "O")) {
      dtcol=4;
    } else if (!strcasecmp(consTorqueCol->data, "B")) {
      dtcol=5;
    }
    else {
      NAMD_die("consTorqueCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Get the column that the torque value is going to be
  //  in. If consTorqueValFile is defined, it can be in any of the 5
  //  floating point fields in the PDB (X, Y, Z, O, or B) which
  //  correspond to the 1st, 2nd, ... 5th floating point fields. If
  //  NEITHER of consTorqueValFile OR consTorqueFile is defined, it can
  //  only be O or B fileds. The default is the O (4th) field, which
  //  is the occupancy.

  if (consTorqueValCol == NULL) {
    dvcol = 4;
  } else {
    if (consTorqueValCol->next != NULL) {
      NAMD_die("Multiple definitions of torque value column in config file");
    };
    
    if (consTorqueValFile == NULL
	&& consTorqueFile == NULL
	&& strcasecmp(consTorqueCol->data, "B")
	&& strcasecmp(consTorqueCol->data, "O")) {
      NAMD_die("Can not read \"constant\" torque values from X, Y, or Z column of the PDB coordinate file");
    };
    if (!strcasecmp(consTorqueValCol->data, "X")) {
      dvcol=1;
    } else if (!strcasecmp(consTorqueValCol->data, "Y")) {
      dvcol=2;
    } else if (!strcasecmp(consTorqueValCol->data, "Z")) {
      dvcol=3;
    } else if (!strcasecmp(consTorqueValCol->data, "O")) {
      dvcol=4;
    } else if (!strcasecmp(consTorqueValCol->data, "B")) {
      dvcol=5;
    }
    else {
      NAMD_die("consTorqueValCol must have value of X, Y, Z, O, or B");
    };
  };
  
  //  Allocate an array that will store an index into the torque
  //  parameters for each atom.  If the atom is not torqued, its
  //  value will be set to -1 in this array.
  consTorqueIndexes = new int32[numAtoms];
  if (consTorqueIndexes == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_constorque_params()");
  };
  
  //  Loop through all the atoms and find out which ones are torqued
  for (i=0; i<numAtoms; i++) {
    switch (dtcol) {
    case 1:
      dtval = (tPDB->atom(i))->xcoor();
      break;
    case 2:
      dtval = (tPDB->atom(i))->ycoor();
      break;
    case 3:
      dtval = (tPDB->atom(i))->zcoor();
      break;
    case 4:
      dtval = (tPDB->atom(i))->occupancy();
      break;
    case 5:
      dtval = (tPDB->atom(i))->temperaturefactor();
      break;
    }
    
    if (dtval != 0.0) {
      //  This atom is torqued
      consTorqueIndexes[i] = current_index;
      current_index++;
    } else {
      //  This atom is not torqued
      consTorqueIndexes[i] = -1;
    }
  }
  
  if (current_index == 0) {
    iout << iWARN << "NO TORQUED ATOMS WERE FOUND, BUT \"CONSTANT\" TORQUE IS ON . . . " << endi;
  } else {
    consTorqueParams = new ConsTorqueParams[current_index];
    if (consTorqueParams == NULL) {
      NAMD_die("memory allocation failed in Molecule::build_constorque_params");
    }
  };
  
  numConsTorque = current_index;
  
  //  Loop through all the atoms and assign the parameters for those
  //  that are torqued
  for (i=0; i<numAtoms; i++) {
    if (consTorqueIndexes[i] != -1) {
      consTorqueParams[consTorqueIndexes[i]].a[0] = (aPDB->atom(i))->xcoor();
      consTorqueParams[consTorqueIndexes[i]].a[1] = (aPDB->atom(i))->ycoor();
      consTorqueParams[consTorqueIndexes[i]].a[2] = (aPDB->atom(i))->zcoor();
      consTorqueParams[consTorqueIndexes[i]].p[0] = (pPDB->atom(i))->xcoor();
      consTorqueParams[consTorqueIndexes[i]].p[1] = (pPDB->atom(i))->ycoor();
      consTorqueParams[consTorqueIndexes[i]].p[2] = (pPDB->atom(i))->zcoor();
      switch (dvcol) {
      case 1:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->xcoor();
	break;
      case 2:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->ycoor();
	break;
      case 3:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->zcoor();
	break;
      case 4:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->occupancy();
	break;
      case 5:
	consTorqueParams[consTorqueIndexes[i]].v = (vPDB->atom(i))->temperaturefactor();
	break;
      };
    };
  };
      
  if (consTorqueFile != NULL) delete tPDB;
  if (consTorqueAxisFile != NULL) delete aPDB;
  if (consTorquePivotFile != NULL) delete pPDB;
  if (consTorqueValFile != NULL) delete vPDB;
}
/*      END OF FUNCTION build_constorque_params    */


/************************************************************************/
/*                  */
/*      FUNCTION build_constant_forces    */
/*                  */
/*   INPUTS:                */
/*  filename - PDB file containing the constant forces    */
/*                  */
/*  This function reads the constant forces from the PDB file.  */
/*   The force vector to be applied on each atom is determined by:    */
/*     occupancy*(X,Y,Z)   */
/*   Only non-zero forces are stored   */
/*                  */
/************************************************************************/

void Molecule::build_constant_forces(char *filename)
{ int i, index;
  PDB *forcePDB;
  
  if (!filename) {
    // then all forces are zero to begin with; may be changed by
    // the consforceconfig command.
    iout << iWARN << "NO CONSTANT FORCES SPECIFIED, BUT CONSTANT FORCE IS ON . . .\n" << endi;
    consForceIndexes = new int32[numAtoms];
    for (i=0; i<numAtoms; i++) consForceIndexes[i] = -1;
    return;
  }

  if ((forcePDB=new PDB(filename)) == NULL)
    NAMD_die("Memory allocation failed in Molecule::build_constant_forces");
  if (forcePDB->num_atoms() != numAtoms)
    NAMD_die("Number of atoms in constant force PDB doesn't match coordinate PDB");

  //  Allocate an array that will store an index into the constant force
  //  array for each atom.  If the atom has no constant force applied, its
  //  value will be set to -1 in this array.
  consForceIndexes = new int32[numAtoms];
  if (consForceIndexes == NULL)
    NAMD_die("memory allocation failed in Molecule::build_constant_forces()");

  //  Loop through all the atoms and find out which ones have constant force
  numConsForce = 0;
  for (i=0; i<numAtoms; i++)
    if ((forcePDB->atom(i)->xcoor()==0 && forcePDB->atom(i)->ycoor()==0 &&
         forcePDB->atom(i)->zcoor()==0) || forcePDB->atom(i)->occupancy()==0)
      //  This atom has no constant force
      consForceIndexes[i] = -1;
    else
      //  This atom has constant force
      consForceIndexes[i] = numConsForce++;

  if (numConsForce == 0)
    // Constant force was turned on, but there weren't really any non-zero forces
    iout << iWARN << "NO NON-ZERO FORCES WERE FOUND, BUT CONSTANT FORCE IS ON . . .\n" << endi;
  else
  { // Allocate an array to hold the forces
    consForce = new Vector[numConsForce];
    if (consForce == NULL)
      NAMD_die("memory allocation failed in Molecule::build_constant_forces");
    // Loop through all the atoms and assign the forces
    for (i=0; i<numAtoms; i++)
      if ((index=consForceIndexes[i]) != -1)
      { //  This atom has constant force on it
        consForce[index].x = forcePDB->atom(i)->xcoor() * forcePDB->atom(i)->occupancy();
        consForce[index].y = forcePDB->atom(i)->ycoor() * forcePDB->atom(i)->occupancy();
        consForce[index].z = forcePDB->atom(i)->zcoor() * forcePDB->atom(i)->occupancy();
      }
  }

  delete forcePDB;
}
/*      END OF FUNCTION build_constant_forces    */


void Molecule::build_langevin_params(BigReal coupling,
    BigReal drudeCoupling, Bool doHydrogen) {

  //  Allocate the array to hold all the data
  langevinParams = new Real[numAtoms];

  if ( (langevinParams == NULL) )
  {
    NAMD_die("memory allocation failed in Molecule::build_langevin_params()");
  }

  //  Loop through all the atoms and get the b value
  for (int i=0; i<numAtoms; i++)
  {
    BigReal bval = coupling;

    if ( (! doHydrogen) && is_hydrogen(i) ) bval = 0;
    else if ( is_lp(i) ) bval = 0;
    else if ( is_drude(i) ) bval = drudeCoupling;

    //  Assign the b value
    langevinParams[i] = bval;
  }

}

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_langevin_params      */
    /*                  */
    /*   INPUTS:                */
    /*  langfile - Value of langevinfile from config file    */
    /*  langcol - Value of langevincol from config file      */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*      cwd - Current working directory          */
    /*                  */
    /*  This function builds the array of b values necessary for  */
    /*   Langevin dynamics.  It takes the name of the PDB file and the      */
    /*   column in the PDB file that contains the b values.  It then  */
    /*   builds the array langevinParams for use during the program.  */
    /*                  */
    /************************************************************************/

    void Molecule::build_langevin_params(StringList *langfile, 
           StringList *langcol, 
           PDB *initial_pdb,
           char *cwd)
       
    {
       PDB *bPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       char filename[129];    //  Filename
       
       //  Get the PDB object that contains the b values.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates.  
       if (langfile == NULL)
       {
    bPDB = initial_pdb;
       }
       else
       {
    if (langfile->next != NULL)
    {
       NAMD_die("Multiple definitions of langvein PDB file in configuration file");
    }

    if ( (cwd == NULL) || (langfile->data[0] == '/') )
    {
         strcpy(filename, langfile->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, langfile->data);
    }
    
    bPDB = new PDB(filename);
    if ( bPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_langevin_params");
    }
    
    if (bPDB->num_atoms() != numAtoms)
    {
       NAMD_die("Number of atoms in langevin parameter PDB doesn't match coordinate PDB");
    }
       }
       
       //  Get the column that the b vaules are in.  It
       //  can be in any of the 5 floating point fields in the PDB, according
       //  to what the user wants.  The allowable fields are X, Y, Z, O, or
       //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
       //  The default is the 4th field, which is the occupancy
       if (langcol == NULL)
       {
    bcol = 4;
       }
       else
       {
    if (langcol->next != NULL)
    {
       NAMD_die("Multiple definitions of langevin parameter column in config file");
    }
    
    if (strcasecmp(langcol->data, "X") == 0)
    {
       bcol=1;
    }
    else if (strcasecmp(langcol->data, "Y") == 0)
    {
       bcol=2;
    }
    else if (strcasecmp(langcol->data, "Z") == 0)
    {
       bcol=3;
    }
    else if (strcasecmp(langcol->data, "O") == 0)
    {
       bcol=4;
    }
    else if (strcasecmp(langcol->data, "B") == 0)
    {
       bcol=5;
    }
    else
    {
       NAMD_die("langevincol must have value of X, Y, Z, O, or B");
    }
       }
       
       //  Allocate the array to hold all the data
       langevinParams = new Real[numAtoms];
       
       if ( (langevinParams == NULL) )
       {
    NAMD_die("memory allocation failed in Molecule::build_langevin_params()");
       }

       //  Loop through all the atoms and get the b value
       for (i=0; i<numAtoms; i++)
       {
    //  Get the k value based on where we were told to find it
    switch (bcol)
    {
       case 1:
    bval = (bPDB->atom(i))->xcoor();
    break;
       case 2:
    bval = (bPDB->atom(i))->ycoor();
    break;
       case 3:
    bval = (bPDB->atom(i))->zcoor();
    break;
       case 4:
    bval = (bPDB->atom(i))->occupancy();
    break;
       case 5:
    bval = (bPDB->atom(i))->temperaturefactor();
    break;
    }
    
    //  Assign the b value
    langevinParams[i] = bval;
       }
       
       //  If we had to create a PDB object, delete it now
       if (langfile != NULL)
       {
    delete bPDB;
       }
    }
    /*      END OF FUNCTION build_langevin_params    */

    /************************************************************************/
    /*                  */
    /*      FUNCTION build_fixed_atoms      */
    /*                  */
    /*   INPUTS:              */
    /*  fixedfile - Value of langevinfile from config file    */
    /*  fixedcol - Value of langevincol from config file    */
    /*  initial_pdb - PDB object that contains initial positions  */
    /*      cwd - Current working directory        */
    /*                  */
    /*  This function builds the list of fixed atoms.      */
    /*   It takes the name of the PDB file and the      */
    /*   column in the PDB file that contains the flags.  It then  */
    /*   builds the array fixedAtomFlags for use during the program.  */
    /*                  */
    /************************************************************************/

    void Molecule::build_fixed_atoms(StringList *fixedfile, 
           StringList *fixedcol, 
           PDB *initial_pdb,
           char *cwd)
       
{
       PDB *bPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       char filename[129];    //  Filename
       
       //  Get the PDB object that contains the b values.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates.  
       if (fixedfile == NULL)
       {
    bPDB = initial_pdb;
       }
       else
       {
    if (fixedfile->next != NULL)
    {
       NAMD_die("Multiple definitions of fixed atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (fixedfile->data[0] == '/') )
    {
         strcpy(filename, fixedfile->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, fixedfile->data);
    }
    
    bPDB = new PDB(filename);
    if ( bPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_fixed_atoms");
    }
    
    if (bPDB->num_atoms() != numAtoms)
    {
       NAMD_die("Number of atoms in fixed atoms PDB doesn't match coordinate PDB");
    }
       }
       
       //  Get the column that the b vaules are in.  It
       //  can be in any of the 5 floating point fields in the PDB, according
       //  to what the user wants.  The allowable fields are X, Y, Z, O, or
       //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
       //  The default is the 4th field, which is the occupancy
       if (fixedcol == NULL)
       {
    bcol = 4;
       }
       else
       {
    if (fixedcol->next != NULL)
    {
       NAMD_die("Multiple definitions of fixed atoms column in config file");
    }
    
    if (strcasecmp(fixedcol->data, "X") == 0)
    {
       bcol=1;
    }
    else if (strcasecmp(fixedcol->data, "Y") == 0)
    {
       bcol=2;
    }
    else if (strcasecmp(fixedcol->data, "Z") == 0)
    {
       bcol=3;
    }
    else if (strcasecmp(fixedcol->data, "O") == 0)
    {
       bcol=4;
    }
    else if (strcasecmp(fixedcol->data, "B") == 0)
    {
       bcol=5;
    }
    else
    {
       NAMD_die("fixedatomscol must have value of X, Y, Z, O, or B");
    }
       }
       
       //  Allocate the array to hold all the data
       fixedAtomFlags = new int32[numAtoms];
       
       if (fixedAtomFlags == NULL)
       {
    NAMD_die("memory allocation failed in Molecule::build_fixed_atoms()");
       }
       
  numFixedAtoms = 0;

       //  Loop through all the atoms and get the b value
       for (i=0; i<numAtoms; i++)
       {
    //  Get the k value based on where we were told to find it
    switch (bcol)
    {
       case 1:
    bval = (bPDB->atom(i))->xcoor();
    break;
       case 2:
    bval = (bPDB->atom(i))->ycoor();
    break;
       case 3:
    bval = (bPDB->atom(i))->zcoor();
    break;
       case 4:
    bval = (bPDB->atom(i))->occupancy();
    break;
       case 5:
    bval = (bPDB->atom(i))->temperaturefactor();
    break;
    }
    
    //  Assign the b value
    if ( bval != 0 ) {
      fixedAtomFlags[i] = 1;
      numFixedAtoms++;
    }
    else {
      fixedAtomFlags[i] = 0;
    }
       }
       
       //  If we had to create a PDB object, delete it now
       if (fixedfile != NULL)
       {
    delete bPDB;
       }

  // now figure out how we interact with rigidBonds 
  // this is mainly for degree of freedom counting
  if ( numRigidBonds ) {
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    int parentIsFixed = 0;
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) {
	parentIsFixed = fixedAtomFlags[h_i->atomID];
	if ( (rigidBondLengths[h_i->atomID] > 0.)  // water
		&& fixedAtomFlags[h_i[1].atomID]
		&& fixedAtomFlags[h_i[2].atomID] ) {
	  ++numFixedRigidBonds;
	}
      } else {
	if ( (rigidBondLengths[h_i->atomID] > 0.)
		&& fixedAtomFlags[h_i->atomID]
		&& parentIsFixed ) {
	  ++numFixedRigidBonds;
	}
      }
    }
  }

  // how many hydrogen groups are completely fixed
  {
    numFixedGroups = 0;
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    int allFixed = 0;
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) {
        if ( allFixed ) ++numFixedGroups;
        allFixed = 1;
      }
      allFixed = allFixed && fixedAtomFlags[h_i->atomID];
    }
    if ( allFixed ) ++numFixedGroups;
  }

}
    /*      END OF FUNCTION build_fixed_atoms    */




/************************************************************************/
/*                                                                      */
/*      FUNCTION build_stirred_atoms                                    */
/*                                                                      */
/*   INPUTS:                                                            */
/*  stirredfile - Value of stirFilename from config file    */
/*  stirredcol - Value of stircol from config file (but B, O only */
/*                    since we need X, Y, Z!     */
/*  initial_pdb - PDB object that contains initial positions  */
/*  cwd - Current working directory        */
/*                  */
/*  This function builds the list of fixed atoms.      */
/*   It takes the name of the PDB file and the      */
/*   column in the PDB file that contains the flags.  It then  */
/*   builds the array fixedAtomFlags for use during the program.  */
/*                                                                      */
/************************************************************************/

    void Molecule::build_stirred_atoms(StringList *stirredfile, 
           StringList *stirredcol, 
           PDB *initial_pdb,
           char *cwd)
       
{
       PDB *sPDB;      //  Pointer to PDB object to use
       int bcol = 4;      //  Column that data is in
       Real bval = 0;      //  b value from PDB file
       int i;      //  Loop counter
       int current_index=0;    //  Index into values used
       char filename[129];    //  Filename
       
       //  Get the PDB object that contains the b values.  If
       //  the user gave another file name, use it.  Otherwise, just use
       //  the PDB file that has the initial coordinates. 
       //  use posssible only if this is 'optional' in simulation parameters
       // dangerous, since restarted simulations will be different
       if (stirredfile == NULL)
       {
    sPDB = initial_pdb;
    // dangerous, since restarted simulations will be different, so warn
        iout << iWARN << "STIRRING USING INITIAL POSITION FILE FOR REFERENCE POSITIONS" << endi;
       }
       else
       {
    if (stirredfile->next != NULL)
    {
       NAMD_die("Multiple definitions of stirred atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (stirredfile->data[0] == '/') )
    {
         strcpy(filename, stirredfile->data);
    }
    else
    {
         strcpy(filename, cwd);
         strcat(filename, stirredfile->data);
    }
        
    //CkPrintf ("DEBUG: the stir filename is %s\n",filename);    
    sPDB = new PDB(filename);

    if ( sPDB == NULL )
    {
      NAMD_die("Memory allocation failed in Molecule::build_stirred_atoms");

    }

    if (sPDB->num_atoms() != numAtoms)
      {
	NAMD_die("Number of atoms in stirred atoms PDB doesn't match coordinate PDB");
      }
 
       }

//  Get the column that the b vaules are in.  It
//  can be in any of the 5 floating point fields in the PDB, according
//  to what the user wants.  The allowable fields are X, Y, Z, O, or
//  B which correspond to the 1st, 2nd, ... 5th floating point fields.
//  The default is the 4th field, which is the occupancy

 
      if (stirredcol == NULL)
      {
	 bcol = 4;
      }
      else
    {
      if (stirredcol->next != NULL)
	{
	  NAMD_die("Multiple definitions of stirred atoms column in config file");
	}
      
      if  (strcasecmp(stirredcol->data, "O") == 0)
	{
	  bcol=4;
	}
      else if (strcasecmp(stirredcol->data, "B") == 0)
	{
	  bcol=5;
	}
      else
	{
	  NAMD_die("stirredAtomsCol must have value of O or B");
	}
    }
          
       //  Allocate an array that will store an index into the stir
       //  parameters for each atom.  If the atom is not stirred, its
       //  value will be set to -1 in this array.
       stirIndexes = new int32[numAtoms];
       
       if (stirIndexes == NULL)
       {
    NAMD_die("memory allocation failed in Molecule::build_stirred_params()");
       }
       
       current_index = 0;
       //  Loop through all the atoms and find out which ones are stirred
       for (i=0; i<numAtoms; i++)
       {



	 //  Get the b value based on where we were told to find it
	 switch (bcol)
	   {

	   case 4:
	     bval = (sPDB->atom(i))->occupancy();
	     break;
	   case 5:
	     bval = (sPDB->atom(i))->temperaturefactor();
	     break;
	   }

	// CkPrintf ("DEBUG--for atom i= %d  bcol= %d    bval= %g   occ= %g   realbval= %g  x= %g numAtoms= %d test= %g\n" ,i        ,bcol, bval, (sPDB->atom(i))->occupancy(), (sPDB->atom(i))->temperaturefactor(),(sPDB->atom(i))->xcoor(), sPDB->num_atoms(), 0.123 );
	 //  Assign the b value
	 if ( bval != 0 )
	 {
	   // This atom is stirred 
	   stirIndexes[i] = current_index;
	   current_index++;
	 }
	 else
	 {
           //This atom is not stirred 
	   stirIndexes[i] = -1;
	 }
       }
	 

    

       
       if (current_index == 0)
       {
    //  Stirring was turned on, but there weren't really any stirred atoms found in file
    iout << iWARN << "NO STIRRED ATOMS WERE FOUND, BUT STIRRING TORQUES ARE ON . . .\n" << endi;
       }
       else
       {
    //  Allocate an array to hold the stirring parameters
    stirParams = new StirParams[current_index];
    
    if (stirParams == NULL)
    {
       NAMD_die("memory allocation failed in Molecule::build_stir_params");
    }
       }
       
       numStirredAtoms = current_index;
       
       //  Loop through all the atoms and assign the parameters for those
       //  that are stirred
       for (i=0; i<numAtoms; i++)
       {
    if (stirIndexes[i] != -1)
    {
       
       //  This atom is stirred, so get the reference position
       stirParams[stirIndexes[i]].refPos.x = (sPDB->atom(i))->xcoor();
       stirParams[stirIndexes[i]].refPos.y = (sPDB->atom(i))->ycoor();
       stirParams[stirIndexes[i]].refPos.z = (sPDB->atom(i))->zcoor();
    }
       }
       
       //  If we had to create a PDB object, delete it now
       if (stirredfile != NULL)
       {
	 delete sPDB;
       }
       

    }

    /*      END OF FUNCTION build_stirred_atoms    */



void Molecule::build_extra_bonds(Parameters *parameters, StringList *file) {
//In the memory optimized version, only the parameters of extraBonds are needed
//to load
  char err_msg[512];
  int a1,a2,a3,a4; float k, ref;
  #ifndef MEM_OPT_VERSION
  ResizeArray<Bond> bonds;
  ResizeArray<Angle> angles;
  ResizeArray<Dihedral> dihedrals;
  ResizeArray<Improper> impropers;
  #endif
  ResizeArray<BondValue> bond_params;
  ResizeArray<AngleValue> angle_params;
  ResizeArray<DihedralValue> dihedral_params;
  ResizeArray<ImproperValue> improper_params;

  if ( ! file ) {
    NAMD_die("NO EXTRA BONDS FILES SPECIFIED");
  }

  for ( ; file; file = file->next ) {  // loop over files
    FILE *f = fopen(file->data,"r");
    if ( ! f ) {
      sprintf(err_msg, "UNABLE TO OPEN EXTRA BONDS FILE %s", file->data);
      NAMD_err(err_msg);
    } else {
      iout << iINFO << "READING EXTRA BONDS FILE " << file->data <<"\n"<<endi;
    }
    
    while ( 1 ) {
      char buffer[512];
      int ret_code;
      do {
        ret_code = NAMD_read_line(f, buffer);
      } while ( (ret_code==0) && (NAMD_blank_string(buffer)) );
      if (ret_code!=0) break;

      char type[512];
      sscanf(buffer,"%s",type);

#define CHECKATOMID(ATOMID) if ( ATOMID < 0 || ATOMID >= numAtoms ) badatom = 1;

      int badline = 0;
      int badatom = 0;
      if ( ! strncasecmp(type,"bond",4) ) {
        if ( sscanf(buffer, "%s %d %d %f %f %s",
	    type, &a1, &a2, &k, &ref, err_msg) != 5 ) badline = 1;
        else {
          CHECKATOMID(a1)
          CHECKATOMID(a2)
        }

        #ifndef MEM_OPT_VERSION              
        Bond tmp;
        tmp.bond_type = parameters->NumBondParams + bonds.size();
        tmp.atom1 = a1;  tmp.atom2 = a2;
        bonds.add(tmp);
        #endif

        BondValue tmpv;
        tmpv.k = k;  tmpv.x0 = ref;
        bond_params.add(tmpv);                
      } else if ( ! strncasecmp(type,"angle",4) ) {
        if ( sscanf(buffer, "%s %d %d %d %f %f %s",
	    type, &a1, &a2, &a3, &k, &ref, err_msg) != 6 ) badline = 1;
        else {
          CHECKATOMID(a1)
          CHECKATOMID(a2)
          CHECKATOMID(a3)
        }
        #ifndef MEM_OPT_VERSION
        Angle tmp;
        tmp.atom1 = a1;  tmp.atom2 = a2;  tmp.atom3 = a3;
        tmp.angle_type = parameters->NumAngleParams + angles.size();
        angles.add(tmp);  
        #endif  

        AngleValue tmpv;
        tmpv.k = k;  tmpv.theta0 = ref / 180. * PI;
        tmpv.k_ub = 0;  tmpv.r_ub = 0;
        angle_params.add(tmpv);      
              
      } else if ( ! strncasecmp(type,"dihedral",4) ) {
        if ( sscanf(buffer, "%s %d %d %d %d %f %f %s",
	    type, &a1, &a2, &a3, &a4, &k, &ref, err_msg) != 7 ) badline = 1;
        else {
          CHECKATOMID(a1)
          CHECKATOMID(a2)
          CHECKATOMID(a3)
          CHECKATOMID(a4)
        }
        #ifndef MEM_OPT_VERSION
        Dihedral tmp;
        tmp.atom1 = a1;  tmp.atom2 = a2;  tmp.atom3 = a3;  tmp.atom4 = a4;
        tmp.dihedral_type = parameters->NumDihedralParams + dihedrals.size();
        dihedrals.add(tmp);
        #endif

        DihedralValue tmpv;
        tmpv.multiplicity = 1;  tmpv.values[0].n = 0;
        tmpv.values[0].k = k;  tmpv.values[0].delta = ref / 180. * PI;
        dihedral_params.add(tmpv);
      } else if ( ! strncasecmp(type,"improper",4) ) {
        if ( sscanf(buffer, "%s %d %d %d %d %f %f %s",
	    type, &a1, &a2, &a3, &a4, &k, &ref, err_msg) != 7 ) badline = 1;
        else {
          CHECKATOMID(a1)
          CHECKATOMID(a2)
          CHECKATOMID(a3)
          CHECKATOMID(a4)
        }
        #ifndef MEM_OPT_VERSION
        Improper tmp;
        tmp.atom1 = a1;  tmp.atom2 = a2;  tmp.atom3 = a3;  tmp.atom4 = a4;
        tmp.improper_type = parameters->NumImproperParams + impropers.size();
        impropers.add(tmp);  
        #endif

        ImproperValue tmpv;
        tmpv.multiplicity = 1;  tmpv.values[0].n = 0;
        tmpv.values[0].k = k;  tmpv.values[0].delta = ref / 180. * PI;
        improper_params.add(tmpv);
      } else if ( ! strncasecmp(type,"#",1) ) {
        continue;  // comment
      } else {
        badline = 1;
      }
#undef CHECKATOMID
      if ( badline ) {
        sprintf(err_msg, "BAD LINE IN EXTRA BONDS FILE %s: %s",
						file->data, buffer);
        NAMD_die(err_msg);
      }
      if ( badatom ) {
        sprintf(err_msg, "BAD ATOM ID IN EXTRA BONDS FILE %s: %s",
						file->data, buffer);
        NAMD_die(err_msg);
      }
    }
    fclose(f);
  }  // loop over files

  // append to parameters and molecule data structures
  int extraNumBonds = bond_params.size();
  if ( extraNumBonds ) {
    iout << iINFO << "READ " << extraNumBonds << " EXTRA BONDS\n" << endi;

    #ifndef MEM_OPT_VERSION
    Bond *newbonds = new Bond[numBonds+extraNumBonds];
    memcpy(newbonds, this->bonds, numBonds*sizeof(Bond));
    memcpy(newbonds+numBonds, bonds.begin(), extraNumBonds*sizeof(Bond));
    delete [] this->bonds;
    this->bonds = newbonds;
    numBonds += extraNumBonds;
    #endif

    BondValue *newbondp = new BondValue[
			parameters->NumBondParams + extraNumBonds];
    memcpy(newbondp, parameters->bond_array,
			parameters->NumBondParams * sizeof(BondValue));
    memcpy(newbondp+parameters->NumBondParams, bond_params.begin(),
			extraNumBonds * sizeof(BondValue));
    delete [] parameters->bond_array;
    parameters->bond_array = newbondp;
    parameters->NumBondParams += extraNumBonds;
  }

  int extraNumAngles = angle_params.size();
  if ( extraNumAngles ) {
    iout << iINFO << "READ " << extraNumAngles << " EXTRA ANGLES\n" << endi;
    #ifndef MEM_OPT_VERSION
    Angle *newangles = new Angle[numAngles+extraNumAngles];
    memcpy(newangles, this->angles, numAngles*sizeof(Angle));
    memcpy(newangles+numAngles, angles.begin(), extraNumAngles*sizeof(Angle));
    delete [] this->angles;
    this->angles = newangles;
    numAngles += extraNumAngles;
    #endif

    AngleValue *newanglep = new AngleValue[
			parameters->NumAngleParams + extraNumAngles];
    memcpy(newanglep, parameters->angle_array,
			parameters->NumAngleParams * sizeof(AngleValue));
    memcpy(newanglep+parameters->NumAngleParams, angle_params.begin(),
			extraNumAngles * sizeof(AngleValue));
    delete [] parameters->angle_array;
    parameters->angle_array = newanglep;
    parameters->NumAngleParams += extraNumAngles;
  }

  int extraNumDihedrals = dihedral_params.size();
  if ( extraNumDihedrals ) {
    iout << iINFO << "READ " << extraNumDihedrals << " EXTRA DIHEDRALS\n" << endi;
    #ifndef MEM_OPT_VERSION
    Dihedral *newdihedrals = new Dihedral[numDihedrals+extraNumDihedrals];
    memcpy(newdihedrals, this->dihedrals, numDihedrals*sizeof(Dihedral));
    memcpy(newdihedrals+numDihedrals, dihedrals.begin(), extraNumDihedrals*sizeof(Dihedral));
    delete [] this->dihedrals;
    this->dihedrals = newdihedrals;
    numDihedrals += extraNumDihedrals;
    #endif

    DihedralValue *newdihedralp = new DihedralValue[
			parameters->NumDihedralParams + extraNumDihedrals];
    memcpy(newdihedralp, parameters->dihedral_array,
			parameters->NumDihedralParams * sizeof(DihedralValue));
    memcpy(newdihedralp+parameters->NumDihedralParams, dihedral_params.begin(),
			extraNumDihedrals * sizeof(DihedralValue));
    delete [] parameters->dihedral_array;
    parameters->dihedral_array = newdihedralp;
    parameters->NumDihedralParams += extraNumDihedrals;
  }

  int extraNumImpropers = improper_params.size();
  if ( extraNumImpropers ) {
    iout << iINFO << "READ " << extraNumImpropers << " EXTRA IMPROPERS\n" << endi;
    #ifndef MEM_OPT_VERSION
    Improper *newimpropers = new Improper[numImpropers+extraNumImpropers];
    memcpy(newimpropers, this->impropers, numImpropers*sizeof(Improper));
    memcpy(newimpropers+numImpropers, impropers.begin(), extraNumImpropers*sizeof(Improper));
    delete [] this->impropers;
    this->impropers = newimpropers;
    numImpropers += extraNumImpropers;
    #endif

    ImproperValue *newimproperp = new ImproperValue[
			parameters->NumImproperParams + extraNumImpropers];
    memcpy(newimproperp, parameters->improper_array,
			parameters->NumImproperParams * sizeof(ImproperValue));
    memcpy(newimproperp+parameters->NumImproperParams, improper_params.begin(),
			extraNumImpropers * sizeof(ImproperValue));
    delete [] parameters->improper_array;
    parameters->improper_array = newimproperp;
    parameters->NumImproperParams += extraNumImpropers;
  }
}// end of Molecule::build_extra_bonds()


//Modifications for alchemical fep
   //
   //  FUNCTION build_fep_flags
   //
   // INPUTS:
   // alchfile - Value of alchfile read from config file
   // alchcol - Value of alch column, read from config file
   // initial_pdb - PDB object that contains the initial positions
   //  cwd - current working directory
   //
   // This function builds the array of state values necessary
   // for FEP or TI. It takes the name of the PDB file and column in
   // the PDB file that contains the alch flag. It then builds
   // the array FepParams for use in the program.
   //
   // function doubles up for TI as well
   
   void Molecule::build_fep_flags(StringList *alchfile,
         StringList *alchcol,
         PDB *initial_pdb,
         char *cwd)
   {
     PDB *bPDB;  //Pointer to PDB object to use
     int bcol = 5;  //Column that the data is in
     Real bval = 0; //flag from PDB file
     int i;         // loop counter
     char filename[129]; // filename

     // get the pdb object that contains the alch flags.
     // if the user gave another filename, use it, else
     // use the pdb file with the initial coordinates
     if (alchfile == NULL) {
       bPDB = initial_pdb;
       strcpy(filename, "coordinate pdb file (default)");
     }
     else {
       if (alchfile->next != NULL) {
        NAMD_die("Multiple definitions of alch PDB file in configuration file");
       }
   
       if ((cwd == NULL) || (alchfile->data[0] == '/')) {
         strcpy(filename, alchfile->data);
       }
       else {
         strcpy(filename, cwd);
         strcat(filename, alchfile->data);
       }

       bPDB = new PDB(filename);
       if (bPDB == NULL) {
         NAMD_die("Memory allocation failed in Molecule:build_fep_flags");
       }

       if (bPDB->num_atoms() != numAtoms) {
         NAMD_die("Number of atoms in alch PDB doesnt match coordinate PDB");
       }
    }
   
    // Get the column that the alch flag is in. It can be in any of the 5 
    // floating point fields in the PDB ie X, Y, Z, O or B.
    // The default is 5th field ie the beta field
    if (alchcol == NULL) {
      bcol = 5;
    }
    else {
      if (alchcol->next != NULL) {
        NAMD_die("Multiple definitions of alch parameter column in config file");
      }

      if (strcasecmp(alchcol->data, "X") == 0) {
       bcol = 1;
      }
      else if (strcasecmp(alchcol->data, "Y") == 0) {
       bcol = 2;
      }
      else if (strcasecmp(alchcol->data, "Z") == 0) {
       bcol = 3;
      }
      else if (strcasecmp(alchcol->data, "O") == 0) {
       bcol = 4;
      }
      else if (strcasecmp(alchcol->data, "B") == 0) {
       bcol = 5;
      }
      else {
       NAMD_die("alchcol must have value of X, Y, Z, O or B");
      }
    }

   iout << iINFO << "To read alch data from file: " << filename << "\n" << endi;
   iout << iINFO << "To read alch flag data from column: " << bcol << "\n" << endi;
 
   //  Allocate the array to hold all the alch data
   fepAtomFlags = new unsigned char[numAtoms];
       
   if (fepAtomFlags == NULL) {
    NAMD_die("Memory allocation failed in Molecule::build_fep_params()");
   }

   double lesMassFactor = 1.0;
   if ( simParams->lesOn && simParams->lesReduceMass ) {
     lesMassFactor = 1.0 / simParams->lesFactor;
   }

   // loop through all the atoms and get the b value
   for (i = 0; i < numAtoms; i++) {
   // Get the alch flag value
    switch (bcol) {
      case 1:
       bval = (bPDB->atom(i))->xcoor();
       break;
      case 2:
       bval = (bPDB->atom(i))->ycoor();
       break;
      case 3:
       bval = (bPDB->atom(i))->zcoor();
       break;
      case 4:
       bval = (bPDB->atom(i))->occupancy();
       break;
      case 5:
       bval = (bPDB->atom(i))->temperaturefactor();
       break;
    }

    // Assign alch flag value
    if (simParams->lesOn) {
      if ( bval == (int) bval && bval > 0 ) {
        if ( bval > simParams->lesFactor ) 
          NAMD_die("LES flag must be less than or equal to lesFactor.");
        fepAtomFlags[i] = (int) bval;
      #ifdef MEM_OPT_VERSION
        Real newMass = atomMassPool[eachAtomMass[i]]*lesMassFactor;
        eachAtomMass[i] = insert_new_mass(newMass);
      #else
        atoms[i].mass *= lesMassFactor;     
      #endif
        numFepFinal++;
        numFepInitial++;
      } else {
        fepAtomFlags[i] = 0;
      }
    } else if (simParams->alchFepOn || simParams->alchThermIntOn) {
      if (bval == 1.0) {
        fepAtomFlags[i] = 1;
        numFepFinal++;
      } else if (bval == -1.0) {
        fepAtomFlags[i] = 2;
        numFepInitial++;
      } else {
        fepAtomFlags[i] = 0;
      }
    } else if (simParams->pairInteractionOn) {
      if (bval == simParams->pairInteractionGroup1) {
        fepAtomFlags[i] = 1;
        numFepInitial++;
      } else if (bval == simParams->pairInteractionGroup2) {
        fepAtomFlags[i] = 2;
        numFepFinal++;
      } else {
        fepAtomFlags[i] = 0;
      }
    } else if (simParams->pressureProfileAtomTypes > 1) {
      fepAtomFlags[i] = (int) bval;
    }
  }

  // if PDB object was created, delete it
  if (alchfile != NULL) {
    delete bPDB;
  }

}
 // End of function build_fep_flags
 
   //
   //
   //  FUNCTION delete_alch_bonded
   //
   // FB - Loop over bonds, angles, dihedrals and impropers, drop any that 
   // contain atoms of both partitions 1 and 2
   // 
   // 

#ifndef MEM_OPT_VERSION
void Molecule::delete_alch_bonded(void)  {

  // Bonds
  suspiciousAlchBonds = 0;  // these really shouldn't exist...?
  for (int i = 0; i < numBonds; i++) {
    int part1 = fepAtomFlags[bonds[i].atom1];
    int part2 = fepAtomFlags[bonds[i].atom2];
    if ((part1 == 1 || part2 == 1 ) &&
      (part1 == 2 || part2 == 2 )) {
      //CkPrintf("-----BOND ATOMS %i %i partitions %i %i \n",bonds[i].atom1, bonds[i].atom2, part1, part2);
      suspiciousAlchBonds++;
    }
  }

  // Angles
  Angle *nonalchAngles;
  nonalchAngles = new Angle[numAngles];
  int nonalchAngleCount = 0;
  alchDroppedAngles = 0;
  for (int i = 0; i < numAngles; i++) {
    int part1 = fepAtomFlags[angles[i].atom1];
    int part2 = fepAtomFlags[angles[i].atom2];
    int part3 = fepAtomFlags[angles[i].atom3];
    if ((part1 == 1 || part2 == 1 || part3 == 1) &&
      (part1 == 2 || part2 == 2 || part3 == 2)) {
      //CkPrintf("-----ANGLE ATOMS %i %i %i partitions %i %i %i\n",angles[i].atom1, angles[i].atom2, angles[i].atom3, part1, part2, part3);
      alchDroppedAngles++;
    }
    else {
      nonalchAngles[nonalchAngleCount++] = angles[i];
    }
  }
  numAngles = nonalchAngleCount;
  delete [] angles;
  angles = new Angle[numAngles];
  for (int i = 0; i < nonalchAngleCount; i++) {
    angles[i]=nonalchAngles[i];
  }
  delete [] nonalchAngles;


  // Dihedrals
  Dihedral *nonalchDihedrals;
  nonalchDihedrals = new Dihedral[numDihedrals];
  int nonalchDihedralCount = 0;
  alchDroppedDihedrals = 0;
  for (int i = 0; i < numDihedrals; i++) {
    int part1 = fepAtomFlags[dihedrals[i].atom1];
    int part2 = fepAtomFlags[dihedrals[i].atom2];
    int part3 = fepAtomFlags[dihedrals[i].atom3];
    int part4 = fepAtomFlags[dihedrals[i].atom4];
    if ((part1 == 1 || part2 == 1 || part3 == 1 || part4 == 1) &&
      (part1 == 2 || part2 == 2 || part3 == 2 || part4 == 2)) {
      //CkPrintf("-----i %i DIHEDRAL ATOMS %i %i %i %i partitions %i %i %i %i\n",i,dihedrals[i].atom1, dihedrals[i].atom2, dihedrals[i].atom3, dihedrals[i].atom4, part1, part2, part3,part4);
      alchDroppedDihedrals++;
    }
    else {
      nonalchDihedrals[nonalchDihedralCount++] = dihedrals[i];
    }
  }
  numDihedrals = nonalchDihedralCount;
  delete [] dihedrals;
  dihedrals = new Dihedral[numDihedrals];
  for (int i = 0; i < numDihedrals; i++) {
    dihedrals[i]=nonalchDihedrals[i];
  }
  delete [] nonalchDihedrals;

  // Impropers
  Improper *nonalchImpropers;
  nonalchImpropers = new Improper[numImpropers];
  int nonalchImproperCount = 0;
  alchDroppedImpropers = 0;
  for (int i = 0; i < numImpropers; i++) {
    int part1 = fepAtomFlags[impropers[i].atom1];
    int part2 = fepAtomFlags[impropers[i].atom2];
    int part3 = fepAtomFlags[impropers[i].atom3];
    int part4 = fepAtomFlags[impropers[i].atom4];
    if ((part1 == 1 || part2 == 1 || part3 == 1 || part4 == 1) &&
      (part1 == 2 || part2 == 2 || part3 == 2 || part4 == 2)) {
      //CkPrintf("-----i %i IMPROPER ATOMS %i %i %i %i partitions %i %i %i %i\n",i,impropers[i].atom1, impropers[i].atom2, impropers[i].atom3, impropers[i].atom4, part1, part2, part3,part4);
      alchDroppedImpropers++;
    }
    else {
      nonalchImpropers[nonalchImproperCount++] = impropers[i];
    }
  }
  numImpropers = nonalchImproperCount;
  delete [] impropers;
  impropers = new Improper[numImpropers];
  for (int i = 0; i < numImpropers; i++) {
    impropers[i]=nonalchImpropers[i];
  }
  delete [] nonalchImpropers;
  
} // end delete_alch_bonded
#endif  

//fepe

void Molecule::build_exPressure_atoms(StringList *fixedfile, 
   StringList *fixedcol, PDB *initial_pdb, char *cwd) {
       
  PDB *bPDB;      //  Pointer to PDB object to use
  int bcol = 4;      //  Column that data is in
  Real bval = 0;      //  b value from PDB file
  int i;      //  Loop counter
  char filename[129];    //  Filename

  //  Get the PDB object that contains the b values.  If
  //  the user gave another file name, use it.  Otherwise, just use
  //  the PDB file that has the initial coordinates.
  if (fixedfile == NULL) {
    bPDB = initial_pdb;
  } else {
    if (fixedfile->next != NULL) {
      NAMD_die("Multiple definitions of excluded pressure atoms PDB file in configuration file");
    }

    if ( (cwd == NULL) || (fixedfile->data[0] == '/') ) {
         strcpy(filename, fixedfile->data);
    } else {
         strcpy(filename, cwd);
         strcat(filename, fixedfile->data);
    }
    bPDB = new PDB(filename);
    if ( bPDB == NULL ) {
      NAMD_die("Memory allocation failed in Molecule::build_exPressure_atoms");
    }

    if (bPDB->num_atoms() != numAtoms) {
      NAMD_die("Number of atoms in excludedPressure atoms PDB doesn't match coordinate PDB");
    }
  }

  //  Get the column that the b vaules are in.  It
  //  can be in any of the 5 floating point fields in the PDB, according
  //  to what the user wants.  The allowable fields are X, Y, Z, O, or
  //  B which correspond to the 1st, 2nd, ... 5th floating point fields.
  //  The default is the 4th field, which is the occupancy
  if (fixedcol == NULL) {
    bcol = 4;
  } else {
    if (fixedcol->next != NULL) {
      NAMD_die("Multiple definitions of excludedPressure atoms column in config file");
    }

    if (strcasecmp(fixedcol->data, "X") == 0) {
       bcol=1;
    } else if (strcasecmp(fixedcol->data, "Y") == 0) {
       bcol=2;
    } else if (strcasecmp(fixedcol->data, "Z") == 0) {
       bcol=3;
    } else if (strcasecmp(fixedcol->data, "O") == 0) {
       bcol=4;
    } else if (strcasecmp(fixedcol->data, "B") == 0) {
       bcol=5;
    } else {
       NAMD_die("excludedPressureFileCol must have value of X, Y, Z, O, or B");
    }
  }

  //  Allocate the array to hold all the data
  exPressureAtomFlags = new int32[numAtoms];

  if (exPressureAtomFlags == NULL) {
    NAMD_die("memory allocation failed in Molecule::build_fixed_atoms()");
  }

  numExPressureAtoms = 0;

  //  Loop through all the atoms and get the b value
  for (i=0; i<numAtoms; i++) {
    //  Get the k value based on where we were told to find it
    switch (bcol) {
       case 1: bval = (bPDB->atom(i))->xcoor(); break;
       case 2: bval = (bPDB->atom(i))->ycoor(); break;
       case 3: bval = (bPDB->atom(i))->zcoor(); break;
       case 4: bval = (bPDB->atom(i))->occupancy(); break;
       case 5: bval = (bPDB->atom(i))->temperaturefactor(); break;
    }

    //  Assign the b value
    if ( bval != 0 ) {
      exPressureAtomFlags[i] = 1;
      numExPressureAtoms++;
    } else {
      exPressureAtomFlags[i] = 0;
    }
  }
  if (fixedfile != NULL) 
    delete bPDB;

  iout << iINFO << "Got " << numExPressureAtoms << " excluded pressure atoms." 
       << endi;
}


    Bool Molecule::is_lp(int anum) {
      return ((atoms[anum].status & LonepairAtom) != 0);
    }

    Bool Molecule::is_drude(int anum) {
      return ((atoms[anum].status & DrudeAtom) != 0);
    }

    Bool Molecule::is_hydrogen(int anum)
    {
  return ((atoms[anum].status & HydrogenAtom) != 0);
    }

    Bool Molecule::is_oxygen(int anum)
    {
  return ((atoms[anum].status & OxygenAtom) != 0);
    }

    Bool Molecule::is_hydrogenGroupParent(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].isGP);
    }

    Bool Molecule::is_water(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].waterVal == 2);
    }

    int Molecule::get_groupSize(int anum)
    {
  return (hydrogenGroup[atoms[anum].hydrogenList].atomsInGroup);
    }

    int Molecule::get_mother_atom(int anum)
    {
  // for efficiency reasons, we are not checking if anum is already 
  // hydrogen or not. This function must be called for hydrogens only;
  return atoms[anum].partner;
    }

void Molecule::reloadCharges(float charge[], int n){
  if ( n != numAtoms )
    NAMD_die("Incorrect number of atoms in Molecule::reloadCharges().");

#ifdef MEM_OPT_VERSION
    delete [] atomChargePool;
    vector<Real> tmpCharges;
    for(int i=0; i<numAtoms; i++){
        int foundIdx=-1;
        //naive searching, better to be binary searching but requiring 
        //inserting charges in increasing/decreasing order
        for(int j=0; j<tmpCharges.size();j++){
        	if(tmpCharges[j] == charge[i]){
        	    foundIdx = j;
        	    break;
        	}
        }
        if(foundIdx==-1){
        	tmpCharges.push_back(charge[i]);
        	foundIdx = tmpCharges.size()-1;
        }
        eachAtomCharge[i] = (Index)foundIdx;
    }
    chargePoolSize = tmpCharges.size();
    atomChargePool = new Real[chargePoolSize];
    for(int i=0; i<chargePoolSize; i++)
        atomChargePool[i] = tmpCharges[i];
#else
  for( int i=0; i<n; ++i ) atoms[i].charge = charge[i];
#endif
}

// go through the molecular structure, analyze the status of each atom,
// and save the data in the Atom structures stored for each atom.  This
// could be built up incrementally while the molecule is being read in,
// but doing it all in one shot allows us to just send the basic info
// over the network and have each node calculate the rest of the data on
// it's own.
void Molecule::build_atom_status(void) {
  register int i;
  int a1, a2, a3;

  // if any atoms have a mass of zero set to 0.001 and warn user
  int numZeroMassAtoms = 0;
  for (i=0; i < numAtoms; i++) {
    #ifdef MEM_OPT_VERSION
    if(atomMassPool[eachAtomMass[i]]<=0){
        ++numZeroMassAtoms;
    }
    #else
    if ( atoms[i].mass <= 0. ) {
      if (simParams->watmodel == WAT_TIP4 ||
          simParams->watmodel == WAT_SWM4) {
        ++numLonepairs;
      } else {
        atoms[i].mass = 0.001;
        ++numZeroMassAtoms;
      }
    }
    else if (atoms[i].mass < 1.) {
      ++numDrudeAtoms;
    }
    #endif
  }
  // DRUDE: verify number of LPs
  if (numLonepairs != numLphosts) {
    NAMD_die("must have same number of LP hosts as lone pairs");
  }
  // DRUDE
  #ifdef MEM_OPT_VERSION
  for(i=0; i<massPoolSize; i++){
      if(atomMassPool[i]<=0.)
          atomMassPool[i] = 0.001;
  }
  #endif

  if (simParams->watmodel == WAT_TIP4 ||
      simParams->watmodel == WAT_SWM4) {
    iout << iWARN << "CORRECTION OF ZERO MASS ATOMS TURNED OFF "
      "BECAUSE LONE PAIRS ARE USED\n" << endi;
  } else {
    if ( numZeroMassAtoms && ! CkMyPe() ) {
      iout << iWARN << "FOUND " << numZeroMassAtoms <<
        " ATOMS WITH ZERO OR NEGATIVE MASSES!  CHANGED TO 0.001\n" << endi;
    }
  }

//In the memory optimization, hydrogen group related information is already
//recorded in the per-Atom file, therefore there's no need re-calculating the
//hydrogen group
//#if 1
#ifndef MEM_OPT_VERSION

  // initialize information for each atom (note that the status has
  // already been initialized during the read/receive phase)
  hydrogenGroup.resize(numAtoms);
  HydrogenGroupID *hg = hydrogenGroup.begin();
  for (i=0; i < numAtoms; i++) {
    atoms[i].partner = (-1);
    hg[i].atomID = i;  // currently unsorted
    hg[i].atomsInGroup = 1;  // currently only 1 in group
    hg[i].isGP = 1;  // assume it is a group parent
    hg[i].GPID = i;  // assume it is a group parent
    hg[i].waterVal = 0;  // for group sorting
  }

  // deal with H-H bonds in a sane manner
  // this information will be rewritten later if bonded elsewhere
  int hhbondcount = 0;
#if 0 
//#ifdef MEM_OPT_VERSION
  for(i=0; i<numAtoms; i++){
    AtomSignature *sig = &atomSigPool[eachAtomSig[i]];
    TupleSignature *bSigs = sig->bondSigs;
    for(int j=0; j<sig->bondCnt; j++){
	a1 = i;
	a2 = i+bSigs[j].offset[0];
    if(!bSigs[j].isReal) continue;
	if (is_hydrogen(a1) && is_hydrogen(a2)) {
	    ++hhbondcount;
	    // make H atoms point at each other for now
	    atoms[a1].partner = a2;
	    atoms[a2].partner = a1;
	    hg[a1].atomsInGroup++;
	    hg[a1].GPID = a2;
	    hg[a2].atomsInGroup++;
	    hg[a2].GPID = a1;
	}
    }
  }
#else
  for (i=0; i < numRealBonds; i++) {
    a1 = bonds[i].atom1;
    a2 = bonds[i].atom2;
    if (is_hydrogen(a1) && is_hydrogen(a2)) {
      ++hhbondcount;
      // make H atoms point at each other for now
      atoms[a1].partner = a2;
      atoms[a2].partner = a1;
      hg[a1].atomsInGroup++;
      hg[a1].GPID = a2;
      hg[a2].atomsInGroup++;
      hg[a2].GPID = a1;
    }
  }
#endif

  if ( hhbondcount && ! CkMyPe() ) {
    iout << iWARN << "Found " << hhbondcount << " H-H bonds.\n" << endi;
  }

  // find which atom each hydrogen is bound to
  // also determine number of atoms in each group
#if 0  
//#ifdef MEM_OPT_VERSION
  for(i=0; i<numAtoms; i++){
    AtomSignature *sig = &atomSigPool[eachAtomSig[i]];
    TupleSignature *bSigs = sig->bondSigs;
    for(int j=0; j<sig->bondCnt; j++){
	a1 = i;
	a2 = i+bSigs[j].offset[0];
    if(!bSigs[j].isReal) continue;
	if (is_hydrogen(a1)) {
	    if (is_hydrogen(a2)) continue;
	    atoms[a1].partner = a2;
	    hg[a2].atomsInGroup++;
	    hg[a1].atomsInGroup = 0;
	    hg[a1].GPID = a2;
	    hg[a1].isGP = 0;
	    // check for waters (put them in their own groups: OH or OHH)
	    if (is_oxygen(a2))  hg[a2].waterVal++;
	}
	if (is_hydrogen(a2)) {
	    atoms[a2].partner = a1;
	    hg[a1].atomsInGroup++;
	    hg[a2].atomsInGroup = 0;
	    hg[a2].GPID = a1;
	    hg[a2].isGP = 0;
	    // check for waters (put them in their own groups: OH or OHH)
	    if (is_oxygen(a1))  hg[a1].waterVal++;
        }
    }
  }
#else
  for (i=0; i < numRealBonds; i++) {
    a1 = bonds[i].atom1;
    a2 = bonds[i].atom2;
    if (is_hydrogen(a1)) {
      if (is_hydrogen(a2)) continue;
      atoms[a1].partner = a2;
      hg[a2].atomsInGroup++;
      hg[a1].atomsInGroup = 0;
      hg[a1].GPID = a2;
      hg[a1].isGP = 0;
      // check for waters (put them in their own groups: OH or OHH)
      if (is_oxygen(a2))  hg[a2].waterVal++;
    }
    if (is_hydrogen(a2)) {
      atoms[a2].partner = a1;
      hg[a1].atomsInGroup++;
      hg[a2].atomsInGroup = 0;
      hg[a2].GPID = a1;
      hg[a2].isGP = 0;
      // check for waters (put them in their own groups: OH or OHH)
      if (is_oxygen(a1))  hg[a1].waterVal++;
    }

    // If we have TIP4P water, check for lone pairs
    if (simParams->watmodel == WAT_TIP4) {
      if (is_lp(a1)) {
        atoms[a1].partner = a2;
        hg[a2].atomsInGroup++;
        hg[a1].atomsInGroup = 0;
        hg[a1].GPID = a2;
        hg[a1].isGP = 0;
      }
      if (is_lp(a2)) {
        atoms[a2].partner = a1;
        hg[a1].atomsInGroup++;
        hg[a2].atomsInGroup = 0;
        hg[a2].GPID = a1;
        hg[a2].isGP = 0;
      }
    }
    // SWM4 water has lone pair and Drude particles
    else if (simParams->watmodel == WAT_SWM4) {
      if (is_lp(a1) || is_drude(a1)) {
        if (is_hydrogen(a2) || is_lp(a2) || is_drude(a2)) {
          char msg[256];
          sprintf(msg, "%s particle %d is bonded to non-parent atom %d",
              (is_lp(a1) ? "Lone pair" : "Drude"), a1+1, a2+1);
          NAMD_die(msg);
        }
        atoms[a1].partner = a2;
        hg[a2].atomsInGroup++;
        hg[a1].atomsInGroup = 0;
        hg[a1].GPID = a2;
        hg[a1].isGP = 0;
      }
      else if (is_lp(a2) || is_drude(a2)) {
        if (is_hydrogen(a1) || is_lp(a1) || is_drude(a1)) {
          char msg[256];
          sprintf(msg, "%s particle %d is bonded to non-parent atom %d",
              (is_lp(a2) ? "Lone pair" : "Drude"), a2+1, a1+1);
          NAMD_die(msg);
        }
        atoms[a2].partner = a1;
        hg[a1].atomsInGroup++;
        hg[a2].atomsInGroup = 0;
        hg[a2].GPID = a1;
        hg[a2].isGP = 0;
      }
    }

  }

#endif

  // check up on our H-H bonds and general sanity check
  int hGPcount = 0;
  for(i=0; i<numAtoms; i++) {
    if ( ! hg[hg[i].GPID].isGP ) {
      char msg[256];
      sprintf(msg, "child atom %d bonded only to child H atoms",i+1);
      NAMD_die(msg);
    }
    if ( hg[i].isGP && is_hydrogen(i) ) {
      if ( hg[i].GPID == i ) continue;  // atomic hydrogen ion
      ++hGPcount;  // molecular hydrogen
      if ( is_hydrogen(hg[i].GPID) && hg[hg[i].GPID].GPID != i ) {
        char msg[256];
        sprintf(msg, "H atom %d bonded only to child H atoms",i+1);
        NAMD_die(msg);
      }
      hg[hg[i].GPID].atomsInGroup = 0;
      hg[hg[i].GPID].isGP = 0;
      hg[i].GPID = i;
      if ( hg[i].atomsInGroup != 2 ) {
        char msg[256];
        sprintf(msg, "H atom %d bonded to multiple H atoms",i+1);
        NAMD_die(msg);
      }
    }
  }
  if ( hGPcount && ! CkMyPe() ) {
    iout << iWARN << "Found " << hGPcount << " H-H molecules.\n" << endi;
  }

  // copy hydrogen groups to migration groups
  for (i=0; i<numAtoms; ++i) {
    if ( hg[i].isGP ) hg[i].GPID = i;  // group parent is its own parent
    else hg[i].waterVal = hg[hg[i].GPID].waterVal;  // copy to children
    hg[i].MPID = hg[i].GPID;
  }

  // determine migration groups based on lone pair hosts
  for (i=0; i<numLphosts; ++i) {
    int a1 = lphosts[i].atom1;
    int a2 = lphosts[i].atom2;
    int a3 = lphosts[i].atom3;
    int a4 = lphosts[i].atom4;
    int m1 = hg[a1].MPID;
    while ( hg[m1].MPID != m1 ) m1 = hg[m1].MPID;
    int m2 = hg[a2].MPID;
    while ( hg[m2].MPID != m2 ) m2 = hg[m2].MPID;
    int m3 = hg[a3].MPID;
    while ( hg[m3].MPID != m3 ) m3 = hg[m3].MPID;
    int m4 = hg[a4].MPID;
    while ( hg[m4].MPID != m4 ) m4 = hg[m4].MPID;
    int mp = m1;
    if ( m2 < mp ) mp = m2;
    if ( m3 < mp ) mp = m3;
    if ( m4 < mp ) mp = m4;
    hg[m1].MPID = mp;
    hg[m2].MPID = mp;
    hg[m3].MPID = mp;
    hg[m4].MPID = mp;
  }
  while ( 1 ) {
    int allok = 1;
    for (i=0; i<numAtoms; ++i) {
      int mp = hg[i].MPID;
      if ( hg[mp].MPID != mp ) {
        allok = 0;
        hg[i].MPID = hg[mp].MPID;
      }
    }
    if ( allok ) break;
  }
  for (i=0; i<numAtoms; ++i) {
    hg[i].isMP = ( hg[i].MPID == i );
    hg[i].atomsInMigrationGroup = 0;
  }
  for (i=0; i<numAtoms; ++i) {
    hg[hg[i].MPID].atomsInMigrationGroup++;
  }

  if ( simParams->splitPatch != SPLIT_PATCH_HYDROGEN ) {
    // every atom its own group
    for (i=0; i<numAtoms; i++) {
      hg[i].isGP = 1;
      hg[i].isMP = 1;
      hg[i].atomsInGroup = 1;
      hg[i].atomsInMigrationGroup = 1;
      hg[i].GPID = i;
      hg[i].MPID = i;
    }
  }

  // count number of groups
  numHydrogenGroups = 0;
  maxHydrogenGroupSize = 0;
  numMigrationGroups = 0;
  maxMigrationGroupSize = 0;
  for(i=0; i<numAtoms; i++)
  {
    if (hg[i].isMP) {
      ++numMigrationGroups;
      int mgs = hg[i].atomsInMigrationGroup;
      if ( mgs > maxMigrationGroupSize ) maxMigrationGroupSize = mgs;
    }
    if (hg[i].isGP) {
      ++numHydrogenGroups;
      int hgs = hg[i].atomsInGroup;
      if ( hgs > maxHydrogenGroupSize ) maxHydrogenGroupSize = hgs;
    }
  }

  hydrogenGroup.sort();

  // sanity checking
  int parentid = -1;
  int hgs = 0;
  for(i=0; i<numAtoms; ++i, --hgs) {
    if ( ! hgs ) {  // expect group parent
      if ( hg[i].isGP ) {
        hgs = hg[i].atomsInGroup;
        parentid = hg[i].atomID;
      } else {
        char buff[512];
        sprintf(buff, "Atom %d has bad hydrogen group size.  "
            "Check for duplicate bonds.", parentid+1);
        NAMD_die(buff);
      }
    } else {  // don't expect group parent
      if ( hg[i].isGP ) {
        char buff[512];
        sprintf(buff, "Atom %d has bad hydrogen group size.  "
            "Check for duplicate bonds.", parentid+1);
        NAMD_die(buff);
      }
    }
  }

  parentid = -1;
  int mgs = 0;
  for(i=0; i<numAtoms; ++i, --mgs) {
    if ( ! mgs ) {  // expect group parent
      if ( hg[i].isMP ) {
        mgs = hg[i].atomsInMigrationGroup;
        parentid = hg[i].atomID;
      } else {
        char buff[512];
        sprintf(buff, "Atom %d has bad migration group size.", parentid+1);
        NAMD_die(buff);
      }
    } else {  // don't expect group parent
      if ( hg[i].isMP ) {
        char buff[512];
        sprintf(buff, "Atom %d has bad migration group size.", parentid+1);
        NAMD_die(buff);
      }
    }
  }


  // finally, add the indexing from atoms[] to hydrogenGroup[]
  for(i=0; i<numAtoms; i++) {
    atoms[hydrogenGroup[i].atomID].hydrogenList = i;
  }

  // check ordering of Drude particles and water
  if (simParams->watmodel == WAT_SWM4) {
    for (i = 0;  i < numAtoms;  i++) {
      if (is_water(hg[i].atomID) && hg[i].isGP) {
        if (i > numAtoms-5
            || ! is_drude(hg[i+1].atomID)
            || ! is_lp(hg[i+2].atomID)
            || ! is_hydrogen(hg[i+3].atomID)
            || ! is_hydrogen(hg[i+4].atomID) ) {
          char msg[256];
          sprintf(msg, "Drude water molecule from HydrogenGroup i=%d "
              "starting at atom %d is not sorted\n", i, hg[i].atomID+1);
          NAMD_die(msg);
        }
        i += 4;  // +1 from loop
        continue;
      } // if water
      else if (is_drude(hg[i].atomID)) {
        if (i < 1 || hg[i-1].atomID != hg[i].GPID) {
          char msg[256];
          sprintf(msg, "Drude particle from HydrogenGroup i=%d must "
              "immediately follow its parent atom %d\n", i, hg[i].GPID+1);
          NAMD_die(msg);
        }
      } // else if Drude
#if 0
      else if (is_lp(hg[i].atomID)) {
        char msg[256];
        sprintf(msg, "Drude lonepair from HydrogenGroup i=%d "
            "at particle %d is NOT from water - unsupported\n",
            i, hg[i].atomID+1);
        NAMD_die(msg);
      }
#endif
    } // for numAtoms
  } // if SWM4

  // set up tail corrections if desired
  if (simParams->LJcorrection) {

      // long-range dispersion correction (FB)
      // Allen and Tildesley, Computer Simulation of Liquids, 1991
      // integrated for NAMD's LJ switch function
      
      // Average values for all particles in the system will be
      // used, as per J Phys Chem B. 2007 111:13052

      // first calculate average A and B coefficients

      int LJtypecount = params->get_num_vdw_params();

      Real A; Real B;
      Real A14; Real B14;
      Real sigma_i, sigma_i14, epsilon_i, epsilon_i14;
      Real sigma_j, sigma_j14, epsilon_j, epsilon_j14;
      Real *ATable = new Real[LJtypecount*LJtypecount];
      Real *BTable = new Real[LJtypecount*LJtypecount];
      int useGeom = simParams->vdwGeometricSigma;
      for (i = 0; i < LJtypecount; i++) {
        for (int j = 0; j < LJtypecount; j++) {
          // A and B calculation code coped from LJTable.C
          if (params->get_vdw_pair_params(i,j, &A, &B, &A14, &B14)) {
            ATable[i*LJtypecount + j] = A;
            BTable[i*LJtypecount + j] = B;
          }
          else {
            params->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
                &epsilon_i14,i);
            params->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14, 
                &epsilon_j14,j);
            BigReal sigma_ij =
              useGeom ? sqrt(sigma_i*sigma_j) : 0.5*(sigma_i+sigma_j);
            BigReal epsilon_ij = sqrt(epsilon_i*epsilon_j);
            sigma_ij *= sigma_ij*sigma_ij;
            sigma_ij *= sigma_ij;

            ATable[i*LJtypecount + j] = 4.0 * sigma_ij * epsilon_ij * sigma_ij;
            BTable[i*LJtypecount + j] = 4.0 * sigma_ij * epsilon_ij;
          }
        }
      }

      int *numAtomsByLjType = new int[LJtypecount];
      for (i = 0; i < LJtypecount; i++) {numAtomsByLjType[i]=0;}
      for (i = 0; i < numAtoms;  i++) {numAtomsByLjType[atoms[i].vdw_type] ++;}

      BigReal sumOfAs = 0; BigReal sumOfBs = 0; int count = 0;
      for (i = 0;  i < LJtypecount;  i++) {
        if (numAtomsByLjType[i]) {
          for (int j = 0;  j < LJtypecount;  j++) {
            if (i == j) {
              sumOfAs += (numAtomsByLjType[i] - 1) * numAtomsByLjType[j] * ATable[i*LJtypecount + j];
              sumOfBs += (numAtomsByLjType[i] - 1) * numAtomsByLjType[j] * BTable[i*LJtypecount + j];
              count += (numAtomsByLjType[i] - 1) * numAtomsByLjType[j];
            }
            else {
              sumOfAs += (numAtomsByLjType[i]) * numAtomsByLjType[j] * ATable[i*LJtypecount + j];
              sumOfBs += (numAtomsByLjType[i]) * numAtomsByLjType[j] * BTable[i*LJtypecount + j];
              count += (numAtomsByLjType[i]) * numAtomsByLjType[j];
            }
          }
        }
      }
      delete [] numAtomsByLjType;
      delete [] ATable;
      delete [] BTable;

      // this naive algorithm (time-consuming n*(n-1)/2) gives the same result
      //BigReal sumOfBs = 0; int count = 0;
      //for (i = 0;  i < numAtoms;  i++) {
      //  for (int j = i + 1; j < numAtoms; j++) {
      //    [get A and B];
      //    sumOfAs += A;
      //    sumOfBs += B;
      //    count ++;
      //  }
      //}


      // discard modified / excluded pairings
      // should be negligible but more correct

      for (i=0; i < numExclusions; i++) {
        int a1 = exclusions[i].atom1;
        int a2 = exclusions[i].atom2;
        if (a1 != a2) { 
          // A and B calculation copied from LJTable.C
          if (params->get_vdw_pair_params(atoms[a1].vdw_type,atoms[a2].vdw_type, 
                &A, &B, &A14, &B14)) {
            sumOfAs -= A;
            sumOfBs -= B;
          }
          else {
            params->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
                &epsilon_i14,atoms[a1].vdw_type);
            params->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14, 
                &epsilon_j14,atoms[a2].vdw_type);
            BigReal sigma_ij =
              useGeom ? sqrt(sigma_i*sigma_j) : 0.5*(sigma_i+sigma_j);
            BigReal epsilon_ij = sqrt(epsilon_i*epsilon_j);

            sigma_ij *= sigma_ij*sigma_ij;
            sigma_ij *= sigma_ij;
            sumOfAs -= 4.0 * sigma_ij * epsilon_ij * sigma_ij;
            sumOfBs -= 4.0 * sigma_ij * epsilon_ij;
          }
          count -= 1;
        }
      }

      BigReal LJAvgA = sumOfAs / count;
      BigReal LJAvgB = sumOfBs / count;
      if ( ! CkMyPe() ) {
        iout << iINFO << "LONG-RANGE LJ: APPLYING ANALYTICAL CORRECTIONS TO "
          << "ENERGY AND PRESSURE\n" << endi; 
        iout << iINFO << "LONG-RANGE LJ: AVERAGE A AND B COEFFICIENTS " 
          << LJAvgA << " AND " << LJAvgB << "\n" << endi; 
      }

      BigReal rcut = simParams->cutoff;
      BigReal rcut2= rcut * rcut;
      BigReal rcut3= rcut2 * rcut;
      BigReal rcut4= rcut2 * rcut2;
      BigReal rcut5= rcut2 * rcut3;
      BigReal rswitch = simParams->switchingDist;
      BigReal rswitch2 = rswitch * rswitch;
      BigReal rswitch3 = rswitch2 * rswitch;
      BigReal rswitch4 = rswitch2 * rswitch2;
      BigReal rswitch5 = rswitch2 * rswitch3;

      if (simParams->switchingActive) {
        tail_corr_ener = tail_corr_virial = (16*numAtoms*numAtoms*PI*(-105*LJAvgB*rcut5*rswitch5 + LJAvgA*(3*rcut4 + 9*rcut3*rswitch + 11*rcut2*rswitch2 + 9*rcut*rswitch3 + 3*rswitch4)))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));
      }
      else {
        tail_corr_virial = -4 * numAtoms * numAtoms * PI * (3 * LJAvgB * rcut3 * rcut3 - 2 * LJAvgA) / (9 * rcut3 * rcut3 * rcut3);
        tail_corr_ener = 2 * numAtoms * numAtoms * PI * (LJAvgA - 3 * LJAvgB * rcut3 * rcut3) / (9 * rcut3 * rcut3 * rcut3); 
      }

    } // LJcorrection 
#endif

  #if 0 
  // debugging code for showing sorted atoms
  if(CkMyPe()==1) {  
  for(i=0; i<numAtoms; i++)
    iout << i << " atomID=" << hydrogenGroup[i].atomID
   << " isGP=" << hydrogenGroup[i].isGP
   << " parent=" << hydrogenGroup[i].GPID
   << " #" << hydrogenGroup[i].atomsInGroup
   << " waterVal=" << hydrogenGroup[i].waterVal
   << " partner=" << atoms[i].partner
   << " hydrogenList=" << atoms[i].hydrogenList
   << "\n" << endi;
  }
  #endif

  // now deal with rigidBonds
  if ( simParams->rigidBonds != RIGID_NONE || simParams->mollyOn ) {
    // temporary variables for use by 4+ site water models
    Real r_oh = -1.0;
    Real r_hh = -1.0;

    delete [] rigidBondLengths;
    rigidBondLengths = new Real[numAtoms];
    if ( ! rigidBondLengths ) {
      NAMD_die("Memory allocation failed in Molecule::build_atom_status()\n");
    }
    for (i=0; i<numAtoms; ++i) rigidBondLengths[i] = 0;
    int mode = simParams->rigidBonds;
    if ( simParams->mollyOn ) mode = RIGID_ALL;

    // add H-mother lengths or 0 if not constrained
#ifdef MEM_OPT_VERSION
  for(i=0; i<numAtoms; i++){
    AtomSignature *sig = &atomSigPool[eachAtomSig[i]];
    TupleSignature *bSigs = sig->bondSigs;
    for(int j=0; j<sig->bondCnt; j++){
	a1 = i;
	a2 = i+bSigs[j].offset[0];
    if(!bSigs[j].isReal) continue;
#else
    for (i=0; i < numRealBonds; i++) {
      a1 = bonds[i].atom1;
      a2 = bonds[i].atom2;
#endif
      Real dum, x0;
    #ifdef MEM_OPT_VERSION
      params->get_bond_params(&dum,&x0,bSigs[j].tupleParamType);
    #else
      params->get_bond_params(&dum,&x0,bonds[i].bond_type);
    #endif
      if (is_hydrogen(a2)) { int tmp = a1;  a1 = a2;  a2 = tmp; } // swap
      if (is_hydrogen(a1)) {
        if ( is_hydrogen(a2) ) {  // H-H
          if ( ! is_water(a2) ) {  // H-H but not water
	    rigidBondLengths[a1] = ( mode == RIGID_ALL ? x0 : 0. );
	    rigidBondLengths[a2] = ( mode == RIGID_ALL ? x0 : 0. );
          }
        } else if ( is_water(a2) || mode == RIGID_ALL ) {
	  rigidBondLengths[a1] = x0;
    if (is_water(a2)) r_oh = rigidBondLengths[a1];
	} else {
	  rigidBondLengths[a1] = 0.;
        }
      }
      // Handle lone pairs if they're allowed
      if (simParams->watmodel == WAT_TIP4) {
        if (is_lp(a2)) { int tmp = a1;  a1 = a2;  a2 = tmp; } // swap
        if (is_lp(a1)) {
          if (! is_water(a2) ) {
            // Currently, lonepairs are only allowed on waters,
            // although this may change in the future
            char err_msg[128];
            sprintf(err_msg, "ILLEGAL LONE PAIR AT INDEX %i\n"
                "LONE PAIRS ARE CURRENTLY ALLOWED ONLY ON WATER MOLECULES\n",
                a1);
            NAMD_die(err_msg);
          } else {
            rigidBondLengths[a1] = x0;
            r_om = x0;
          }
        }
      }
      // Handle SWM4 lone pairs
      // (Drude bonds remain flexible)
      if (simParams->watmodel == WAT_SWM4) {
        if (is_lp(a2)) {
          int tmp = a1;  a1 = a2;  a2 = tmp;  // swap
        }
        if (is_lp(a1)) {
          if (is_water(a2)) {
            // do not count bonds to LPs as rigid, do not set rigidBondLengths[]
            r_om = x0;  // for faster position update routine for LP on water
          }
          else if ( ! simParams->drudeOn) {
            // if not using Drude model, lone pairs allowed only on water
            char msg[128];
            sprintf(msg, "ILLEGAL LONE PAIR AT INDEX %d\n"
                "LONE PAIRS ARE CURRENTLY ALLOWED ONLY ON WATER MOLECULES\n",
                a1+1);
            NAMD_die(msg);
          }
        }
      }

#ifdef MEM_OPT_VERSION
    }
#endif
    }

    // zero out H-H lengths - water handled below
    HydrogenGroup::iterator h_i, h_e;
    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP ) rigidBondLengths[h_i->atomID] = 0.;
    }

    // fill in H-H lengths for water by searching angles - yuck
#ifdef MEM_OPT_VERSION
  for(i=0; i<numAtoms; i++){
    AtomSignature *sig = &atomSigPool[eachAtomSig[i]];
    TupleSignature *aSigs = sig->angleSigs;
    for(int j=0; j<sig->angleCnt; j++){
	a1 = i;
	a2 = i+aSigs[j].offset[0];
	if(!is_water(a2)) continue;
	a3 = i+aSigs[j].offset[1];
	if ( rigidBondLengths[a1] != rigidBondLengths[a3] ) {
	    NAMD_die("Asymmetric water molecule found???  This can't be right.\n");
	}
	Real dum, t0;
	params->get_angle_params(&dum,&t0,&dum,&dum,aSigs[j].tupleParamType);
	rigidBondLengths[a2] = 2. * rigidBondLengths[a1] * sin(0.5*t0);
      r_hh = rigidBondLengths[a2];
    }
}
#else
    for (i=0; i < numAngles; i++) {
      a2 = angles[i].atom2;
      if ( ! is_water(a2) ) continue;
      a1 = angles[i].atom1;
      a3 = angles[i].atom3;
      if (is_lp(a2) || is_lp(a1) || is_lp(a3) ||
          is_drude(a2) || is_drude(a1) || is_drude(a3)) continue;
      if ( rigidBondLengths[a1] != rigidBondLengths[a3] ) {
        if (rigidBondLengths[a1] >0.3 && rigidBondLengths[a3] >0.3) {
          printf("length1: %f length2: %f\n", rigidBondLengths[a1], rigidBondLengths[a3]);

          NAMD_die("Asymmetric water molecule found???  This can't be right.\n");
        }
      }
      Real dum, t0;
      params->get_angle_params(&dum,&t0,&dum,&dum,angles[i].angle_type);
      rigidBondLengths[a2] = 2. * rigidBondLengths[a1] * sin(0.5*t0);
      r_hh = rigidBondLengths[a2];
    }
#endif

    // fill in H-H lengths for waters that are missing angles
    int numBondWaters = 0;
    int numFailedWaters = 0;

#ifdef MEM_OPT_VERSION
  for(i=0; i<numAtoms; i++){
    AtomSignature *sig = &atomSigPool[eachAtomSig[i]];
    TupleSignature *bSigs = sig->bondSigs;
    for(int j=0; j<sig->bondCnt; j++){
	a1 = i;
	a2 = i+bSigs[j].offset[0];
    if(!bSigs[j].isReal) continue;
#else
    for (i=0; i < numRealBonds; i++) {
      a1 = bonds[i].atom1;
      a2 = bonds[i].atom2;
#endif
      if ( ! is_hydrogen(a1) ) continue;
      if ( ! is_hydrogen(a2) ) continue;
      int ma1 = get_mother_atom(a1);
      int ma2 = get_mother_atom(a2);
      if ( ma1 != ma2 ) continue;
      if ( ! is_water(ma1) ) continue;
      if ( rigidBondLengths[ma1] != 0. ) continue;
      Real dum, x0;
    #ifdef MEM_OPT_VERSION
      params->get_bond_params(&dum,&x0,bSigs[j].tupleParamType);
    #else
      params->get_bond_params(&dum,&x0,bonds[i].bond_type);
    #endif
      rigidBondLengths[ma1] = x0;
    }
#ifdef MEM_OPT_VERSION
    }
#endif
    
    // We now should be able to set the parameters needed for water lonepairs
    if (simParams->watmodel == WAT_TIP4) {
      if (r_oh < 0.0 || r_hh < 0.0) {
        printf("ERROR: r_oh %f / r_hh %f\n", r_oh, r_hh);
        NAMD_die("Failed to find water bond lengths\n");
      } 
      r_ohc = sqrt(r_oh * r_oh - 0.25 * r_hh * r_hh);
      printf("final r_om and r_ohc are %f and %f\n", r_om, r_ohc);
    }

    h_i = hydrogenGroup.begin();  h_e = hydrogenGroup.end();
    for( ; h_i != h_e; ++h_i ) {
      if ( h_i->isGP && is_water(h_i->atomID) &&
                     rigidBondLengths[h_i->atomID] == 0. ) {
        if ( h_i + 1 == h_e || h_i + 2 == h_e ||
             h_i[1].isGP || h_i[2].isGP || h_i->atomsInGroup != 3 ) {
          NAMD_die("Abnormal water detected.");
        }
        Bond btmp;
        btmp.atom1 = h_i[1].atomID;
        char atom1name[11];
	strcpy(atom1name,get_atomtype(btmp.atom1));
        btmp.atom2 = h_i[2].atomID;
        char atom2name[11];
	strcpy(atom2name,get_atomtype(btmp.atom2));
        params->assign_bond_index(atom1name,atom2name,&btmp);
        Real k, x0;
	x0 = 0.;
        params->get_bond_params(&k,&x0,btmp.bond_type);
	if ( x0 > 0. ) {
          rigidBondLengths[h_i->atomID] = x0;
	  numBondWaters++;
        } else {
	  numFailedWaters++;
        }
      }
    }
    if ( numBondWaters + numFailedWaters ) {
      iout << iWARN << "Missing angles for " <<
	      ( numBondWaters + numFailedWaters ) << " waters.\n" << endi;
    }
    if ( numBondWaters ) {
      iout << iWARN << "Obtained H-H distance from bond parameters for " <<
	      numBondWaters << " waters.\n" << endi;
    }
    if ( numFailedWaters ) {
      iout << iERROR << "Failed to obtain H-H distance from angles or bonds for " <<
	      numFailedWaters << " waters.\n" << endi;
    }

    // in case both molly and rigidBonds are in use make lengths which
    // are molly-only negative and leave lengths which are both alone
    if ( simParams->mollyOn ) {
      mode = simParams->rigidBonds;
      if ( mode == RIGID_NONE ) {
	for (i=0; i<numAtoms; ++i) rigidBondLengths[i] *= -1;
      } else if ( mode == RIGID_WATER ) {
	for (i=0; i<numAtoms; ++i) {
	  if ( ! is_water(i) ) rigidBondLengths[i] *= -1;
	}
      }
    }

    numRigidBonds = 0;
    for (i=0; i<numAtoms; ++i) {
      if ( rigidBondLengths[i] > 0. ) ++numRigidBonds;
    }

  }

}

#ifdef MEM_OPT_VERSION
//idx1: atom1's exclusion check signature
//to check whether atom1 and atom2 are excluded from each other
int Molecule::checkExclByIdx(int idx1, int atom1, int atom2) const {

  int amin = exclChkSigPool[idx1].min;
  int amax = exclChkSigPool[idx1].max;
  int dist21 = atom2 - atom1;
  if ( dist21 < amin || dist21 > amax ) return 0;
  else return exclChkSigPool[idx1].flags[dist21-amin];

}
#else
int Molecule::checkexcl(int atom1, int atom2) const {

  int amin = all_exclusions[atom1].min;
  int amax = all_exclusions[atom1].max;
  if ( atom2 < amin || atom2 > amax ) return 0;
  else return all_exclusions[atom1].flags[atom2-amin];

}
#endif


/************************************************************************/
/*                  */
/*      FUNCTION Molecule        */
/*                  */
/*  This is the constructor for reading AMBER topology data    */
/*                  */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param, Ambertoppar *amber_data)
{
  initialize(simParams,param);

  read_parm(amber_data);
}
/*      END OF FUNCTION Molecule      */


/************************************************************************/
/*                  */
/*      FUNCTION read_parm   */
/*                  */
/*   INPUTS:                */
/*  amber_data - AMBER data structure    */
/*                  */
/*  This function copys AMBER topology data to the corresponding data  */
/*   structures      */
/*                  */
/************************************************************************/

void Molecule::read_parm(Ambertoppar *amber_data)
{
#ifdef MEM_OPT_VERSION
    NAMD_die("When reading a compressed file or using the memory-optimized version, amber data is not supported!");    
#else
  int i,j,ntheth,nphih,current_index,a1,a2,
      max,min,index,found;

  if (!amber_data->data_read)
    NAMD_die("No data read from parm file yet!");

  // Copy atom informations
  numAtoms = amber_data->Natom;
  atoms = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];

  if(simParams->genCompressedPsf) {
      atomSegResids = new AtomSegResInfo[numAtoms];
  }

  if (atoms == NULL || atomNames == NULL )
    NAMD_die("memory allocation failed when reading atom information");
  ResidueLookupElem *tmpResLookup = resLookup;
  for (i=0; i<numAtoms; ++i)
  { atomNames[i].resname = nameArena->getNewArray(5);
    atomNames[i].atomname = nameArena->getNewArray(5);
    atomNames[i].atomtype = nameArena->getNewArray(5);
    if (atomNames[i].resname == NULL || atomNames[i].atomname == NULL || atomNames[i].atomtype == NULL)
      NAMD_die("memory allocation failed when reading atom information");
    for (j=0; j<4; ++j)
    { atomNames[i].resname[j] = amber_data->ResNames[amber_data->AtomRes[i]*4+j];
      atomNames[i].atomname[j] = amber_data->AtomNames[i*4+j];
      atomNames[i].atomtype[j] = amber_data->AtomSym[i*4+j];
    }
    atomNames[i].resname[4] = atomNames[i].atomname[4] = atomNames[i].atomtype[4] = '\0';
    strtok(atomNames[i].resname," ");
    strtok(atomNames[i].atomname," ");
    strtok(atomNames[i].atomtype," ");
    atoms[i].mass = amber_data->Masses[i];
    // Divide by 18.2223 to convert to charge in units of the electron charge
    atoms[i].charge = amber_data->Charges[i] / 18.2223;
    atoms[i].vdw_type = amber_data->Iac[i] - 1;

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append("MAIN", amber_data->AtomRes[i]+1, i);

    if(atomSegResids) { //for compressing molecule information
        AtomSegResInfo *one = atomSegResids + i;
        memcpy(one->segname, "MAIN", strlen("MAIN")+1);
        one->resid = amber_data->AtomRes[i]+1;
    }
    

    /*  Determine the type of the atom (H or O) */
    atoms[i].status = UnknownAtom; // the default
    if (atoms[i].mass <= 0.05) {
      atoms[i].status |= LonepairAtom;
    } else if (atoms[i].mass < 1.0) {
      atoms[i].status |= DrudeAtom;
    } else if (atoms[i].mass <=3.5) {
      atoms[i].status |= HydrogenAtom;
    } else if ((atomNames[i].atomname[0] == 'O') && 
         (atoms[i].mass >= 14.0) && 
         (atoms[i].mass <= 18.0)) {
      atoms[i].status |= OxygenAtom;
    }
  }

// Note: In AMBER, the atom numbers in bond, angle and dihedral arrays are in fact
// (3*(atnum-1)). So we divide them by 3 to get the real indices of atom array. Also
// note that NAMD indexes arrays from 0 to NumAtoms-1.

  // Copy bond information
  // Fake bonds (bonds with 0 force constant) are ignored
  Real k, x0;
  numBonds = 0;
  if (amber_data->Nbonh + amber_data->Nbona > 0)
  { bonds = new Bond[amber_data->Nbonh + amber_data->Nbona];
    if (bonds == NULL || amber_data->Nbona < 0)
      NAMD_die("memory allocation failed when reading bond information");
    // Bonds WITH hydrogen
    for (i=0; i<amber_data->Nbonh; ++i)
    { bonds[numBonds].atom1 = amber_data->BondHAt1[i] / 3;
      bonds[numBonds].atom2 = amber_data->BondHAt2[i] / 3;
      bonds[numBonds].bond_type = amber_data->BondHNum[i] - 1;
      if (bonds[numBonds].atom1>=numAtoms || bonds[numBonds].atom2>=numAtoms ||
          bonds[numBonds].bond_type>=amber_data->Numbnd)
      { char err_msg[128];
        sprintf(err_msg, "BOND (WITH H) # %d OVERFLOW IN PARM FILE", i+1);
        NAMD_die(err_msg);
      }
      params->get_bond_params(&k,&x0,bonds[numBonds].bond_type);
      // if ( k != 0. ) ++numBonds;  // real bond
      ++numBonds;  // keep all bonds in case needed for rigid water
    }
    // Bonds WITHOUT hydrogen
    for (i=amber_data->Nbonh; i<amber_data->Nbonh+amber_data->Nbona; ++i)
    { bonds[numBonds].atom1 = amber_data->BondAt1[i-amber_data->Nbonh] / 3;
      bonds[numBonds].atom2 = amber_data->BondAt2[i-amber_data->Nbonh] / 3;
      bonds[numBonds].bond_type = amber_data->BondNum[i-amber_data->Nbonh] - 1;
      if (bonds[i].atom1>=numAtoms || bonds[i].atom2>=numAtoms ||
          bonds[i].bond_type>=amber_data->Numbnd)
      { char err_msg[128];
        sprintf(err_msg, "BOND (WITHOUT H) # %d OVERFLOW IN PARM FILE", i+1-amber_data->Nbonh);
        NAMD_die(err_msg);
      }
      params->get_bond_params(&k,&x0,bonds[numBonds].bond_type);
      // if ( k != 0. ) ++numBonds;  // real bond
      ++numBonds;  // keep all bonds in case needed for rigid water
    }
  }
  /*  Tell user about our subterfuge  */
  if ( numBonds !=  amber_data->Nbonh + amber_data->Nbona) {
    iout << iWARN << "Ignored " << amber_data->Nbonh + amber_data->Nbona - numBonds <<
            " bonds with zero force constants.\n" << endi;
    iout << iWARN <<
	"Will get H-H distance in rigid H2O from H-O-H angle.\n" << endi;
  }

  // Copy angle information
  numAngles = amber_data->Ntheth + amber_data->Ntheta;
  if (numAngles > 0)
  { ntheth = amber_data->Ntheth;
    angles = new Angle[numAngles];
    if (angles == NULL || numAngles < ntheth)
      NAMD_die("memory allocation failed when reading angle information");
    // Angles WITH hydrogen
    for (i=0; i<ntheth; ++i)
    { angles[i].atom1 = amber_data->AngleHAt1[i] / 3;
      angles[i].atom2 = amber_data->AngleHAt2[i] / 3;
      angles[i].atom3 = amber_data->AngleHAt3[i] / 3;
      angles[i].angle_type = amber_data->AngleHNum[i] - 1;
      if (angles[i].atom1>=numAtoms || angles[i].atom2>=numAtoms ||
          angles[i].atom3>=numAtoms || angles[i].angle_type>=amber_data->Numang)
      { char err_msg[128];
        sprintf(err_msg, "ANGLE (WITH H) # %d OVERFLOW IN PARM FILE", i+1);
        NAMD_die(err_msg);
      }
    }
    // Angles WITHOUT hydrogen
    for (i=ntheth; i<numAngles; ++i)
    { angles[i].atom1 = amber_data->AngleAt1[i-ntheth] / 3;
      angles[i].atom2 = amber_data->AngleAt2[i-ntheth] / 3;
      angles[i].atom3 = amber_data->AngleAt3[i-ntheth] / 3;
      angles[i].angle_type = amber_data->AngleNum[i-ntheth] - 1;
      if (angles[i].atom1>=numAtoms || angles[i].atom2>=numAtoms ||
          angles[i].atom3>=numAtoms || angles[i].angle_type>=amber_data->Numang)
      { char err_msg[128];
        sprintf(err_msg, "ANGLE (WITHOUT H) # %d OVERFLOW IN PARM FILE", i+1-ntheth);
        NAMD_die(err_msg);
      }
    }
  }

  numExclusions =  0;
  // If readExclusions is TRUE, then we copy exclusions from parm
  // file; otherwise we skip the exclusions here and generate
  // them later in build_exclusions()
  if (simParams->readExclusions)
  { // Copy exclusion information
    // In Amber data structure, Iblo[] is the number of exclusions
    // for each atom; ExclAt[] is the atom index for the excluded atoms.
    exclusions = new Exclusion[amber_data->Nnb];
    if (exclusions == NULL &&  amber_data->Nnb > 0)
      NAMD_die("memory allocation failed when reading exclusion information");
    current_index = 0;
    for (i=0; i<numAtoms; ++i)
      for (j=0; j<amber_data->Iblo[i]; ++j)
      { if (current_index >= amber_data->Nnb)
        { char err_msg[128];
          sprintf(err_msg, "EXCLUSION INDEX EXCEEDS NUMBER OF EXLCUSIONS %d IN AMBER FILE, AT ATOM #%d\n",
            amber_data->Nnb, i+1);
          NAMD_die(err_msg);
        }
        // There's some 0 in the ExclAt[] list, which is strange
        // and redundant. In this case, I simply ignore such entries.
        if (amber_data->ExclAt[current_index] != 0)
        { // Subtract 1 to convert the index from the 1 to NumAtoms
          // used in the file to the 0 to NumAtoms-1 that we need
          a2 = amber_data->ExclAt[current_index] - 1;
          if (a2 < i)
          { // I assume the latter index be larger than the former
            // one, so that the same exclusion won't be double-counted;
            // if not, give error
            char err_msg[128];
            sprintf(err_msg, "Atom #%d has exclusion with atom #%d, in reverse order.", i+1, a2+1);
            NAMD_die(err_msg);
          }
          else if (a2 == i)
          { char err_msg[128];
            sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN AMBER FILE\n", i+1);
              NAMD_die(err_msg);
          }
          else if (a2 >= numAtoms)
          {  char err_msg[128];
             sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN AMBER FILE",
               a2+1, numAtoms, current_index+1);
             NAMD_die(err_msg);
          }
          exclusions[numExclusions].atom1 = i;
          exclusions[numExclusions].atom2 = a2;
          ++numExclusions;
        }
        ++current_index;
      }
    if (current_index < amber_data->Nnb)
    { char err_msg[128];
      sprintf(err_msg, "Num of exclusions recorded (%d) is smaller than what it's supposed to be (%d)",
        current_index,amber_data->Nnb);
      NAMD_die(err_msg);
    }
  }

  // Copy dihedral information
  numDihedrals = amber_data->Nphih + amber_data->Nphia;
  if (numDihedrals > 0)
  { nphih = amber_data->Nphih;
    dihedrals = new Dihedral[numDihedrals];
    if (dihedrals == NULL || numDihedrals < nphih)
      NAMD_die("memory allocation failed when reading dihedral information");
    // Dihedral WITH hydrogen
    for (i=0; i<nphih; ++i)
    { dihedrals[i].atom1 = amber_data->DihHAt1[i] / 3;
      dihedrals[i].atom2 = amber_data->DihHAt2[i] / 3;
      dihedrals[i].atom3 = amber_data->DihHAt3[i] / 3;
      dihedrals[i].atom4 = amber_data->DihHAt4[i] / 3;
      dihedrals[i].dihedral_type = amber_data->DihHNum[i] - 1;
    }
    // Dihedral WITHOUT hydrogen
    for (i=nphih; i<numDihedrals; ++i)
    { dihedrals[i].atom1 = amber_data->DihAt1[i-nphih] / 3;
      dihedrals[i].atom2 = amber_data->DihAt2[i-nphih] / 3;
      dihedrals[i].atom3 = amber_data->DihAt3[i-nphih] / 3;
      dihedrals[i].atom4 = amber_data->DihAt4[i-nphih] / 3;
      dihedrals[i].dihedral_type = amber_data->DihNum[i-nphih] - 1;
    }
  }
  // In AMBER parm file, dihedrals contain 1-4 exclusion infomation:
  // the 1st and 4th atoms have 1-4 nonbond interation. So we should
  // find them in the exclusion array and change their exclusion to
  // 1-4 type. However, there're two exceptions --
  // 1.If the third atom is negative, it means the end group
  //   interactions are to be ignored;
  // 2.If the fourth atom is negative, it means this is an improper.
  // For the above two cases, the actual atom index is the absolute
  // value of the atom number read; and there's no 1-4 interation
  // for these dihedrals.
  // If readExclusions is not TRUE, then we don't worry about
  // exclusions here.
  for (i=0; i<numDihedrals; ++i)
  { if (dihedrals[i].atom3 < 0 || dihedrals[i].atom4 < 0)
    { dihedrals[i].atom3 = abs(dihedrals[i].atom3);
      dihedrals[i].atom4 = abs(dihedrals[i].atom4);
    }
    else if (simParams->readExclusions)
    { if (dihedrals[i].atom1 < dihedrals[i].atom4)
        a1=dihedrals[i].atom1, a2=dihedrals[i].atom4;
      else
        a1=dihedrals[i].atom4, a2=dihedrals[i].atom1;
      // Since in the exclusion array, atom1 is guaranteed to be
      // ordered, we can do a binary serch to find it first.
      found = 0;
      min=0, max=numExclusions-1;
      while (!found && min<=max)
      { index = (min+max)/2;
        if (exclusions[index].atom1 == a1)
          found = 1;
        else if (exclusions[index].atom1 < a1)
          min = index+1;
        else
          max = index-1;
      }
      if (!found)
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
      // After finding atom1, we do a linear serch to find atom2,
      // in both directions.
      for (j=index-1; j>=0 && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; --j);
      if (j<0 || exclusions[j].atom1!=a1)
        for (j=index; j<numExclusions && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; ++j);
      if (j<numExclusions && exclusions[j].atom1==a1)
        exclusions[j].modified = 1;  // Change the exclusion type to 1-4
      else
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
    }
    if (dihedrals[i].atom1>=numAtoms || dihedrals[i].atom2>=numAtoms ||
        dihedrals[i].atom3>=numAtoms || dihedrals[i].atom4>=numAtoms ||
        dihedrals[i].dihedral_type>=amber_data->Nptra)
    { char err_msg[128];
      sprintf(err_msg, "DIHEDRAL # %d OVERFLOW IN PARM FILE", i+1);
      NAMD_die(err_msg);
    }
  }
  
  //  analyze the data and find the status of each atom
  numRealBonds = numBonds;
  build_atom_status();
#endif
}
/*      END OF FUNCTION read_parm    */


/************************************************************************/
/*                                                                      */
/*      FUNCTION Molecule                                               */
/*                                                                      */
/*  This is the constructor for reading GROMACS topology data           */
/*                                                                      */
/************************************************************************/

Molecule::Molecule(SimParameters *simParams, Parameters *param,
		   const GromacsTopFile *gromacsTopFile)
{
  initialize(simParams,param);

  read_parm(gromacsTopFile);
}
/*      END OF FUNCTION Molecule      */

/************************************************************************/
/*                                                                      */
/*      FUNCTION read_parm                                              */
/*                                                                      */
/*   INPUTS:                                                            */
/*  amber_data - AMBER data structure                                   */
/*                                                                      */
/*  This function copys AMBER topology data to the corresponding data   */
/*   structures                                                         */
/*                                                                      */
/************************************************************************/

void Molecule::read_parm(const GromacsTopFile *gf) {
#ifdef MEM_OPT_VERSION
    NAMD_die("When reading a compressed file or using the memory-optimized version, amber data is not supported!");    
#else
  /*  int i,j,ntheth,nphih,current_index,a1,a2,
      max,min,index,found;*/
  int i;
  
  // Initializes the atom array
  numAtoms = gf->getNumAtoms();
  atoms = new Atom[numAtoms];
  atomNames = new AtomNameInfo[numAtoms];

  if(simParams->genCompressedPsf) {
      atomSegResids = new AtomSegResInfo[numAtoms];
  }

  if (atoms == NULL || atomNames == NULL )
    NAMD_die("memory allocation failed when reading atom information");
  ResidueLookupElem *tmpResLookup = resLookup;

  // Copy the individual atoms over
  for (i=0; i<numAtoms; ++i) {
    char *resname = nameArena->getNewArray(11);
    char *atomname = nameArena->getNewArray(11);
    char *atomtype = nameArena->getNewArray(11);
    int resnum,typenum;
    Real charge,mass;

    if (resname == NULL || atomname == NULL || atomtype == NULL)
      NAMD_die("memory allocation failed when reading atom information");

    // get the data out of the GROMACS file
    gf->getAtom(i,&resnum,resname,atomname,atomtype,&typenum,
		&charge,&mass);

    atomNames[i].resname = resname;
    atomNames[i].atomname = atomname;
    atomNames[i].atomtype = atomtype;
    atoms[i].mass = mass;
    atoms[i].charge = charge;
    atoms[i].vdw_type = typenum;

    /*  Add this atom to residue lookup table */
    if ( tmpResLookup ) tmpResLookup =
	tmpResLookup->append("MAIN", resnum+1, i);

    if(atomSegResids) { //for compressing molecule information
        AtomSegResInfo *one = atomSegResids + i;
        memcpy(one->segname, "MAIN", strlen("MAIN")+1);
        one->resid = resnum+1;
    }

    /*  Determine the type of the atom (H or O) */
    // XXX this cannot be done this way in GROMACS
    // For example, in dppc LO2 appears to be an oxygen.
    // And how do the hydrogens in CH3 etc factor in to this?
    atoms[i].status = UnknownAtom; // the default
    if (atoms[i].mass <= 0.05) {
      atoms[i].status |= LonepairAtom;
    } else if (atoms[i].mass < 1.0) {
      atoms[i].status |= DrudeAtom;
    } else if (atoms[i].mass <=3.5) {
      atoms[i].status |= HydrogenAtom;
    } else if ((atomNames[i].atomname[0] == 'O') && 
         (atoms[i].mass >= 14.0) && 
         (atoms[i].mass <= 18.0)) {
      atoms[i].status |= OxygenAtom;
    }
  }

  // Copy bond information
  numBonds = gf->getNumBonds();
  bonds = new Bond[numBonds];
  if (bonds == NULL)
    NAMD_die("memory allocation failed when reading bond information");
  for(i=0;i<numBonds;i++) {
    int type; // to convert the type correctly
    int atom1,atom2;
    gf->getBond(i,&atom1,&atom2,&type);
    bonds[i].atom1 = atom1;
    bonds[i].atom2 = atom2;
    bonds[i].bond_type = (Index)type;
  }

  // Copy angle information
  numAngles = gf->getNumAngles();
  angles = new Angle[numAngles];
  if (angles == NULL)
    NAMD_die("memory allocation failed when reading angle information");
  for(i=0;i<numAngles;i++) {
    int type; // to convert the type correctly
    int atom1,atom2,atom3;
    gf->getAngle(i,&atom1,&atom2,&atom3,&type);

    angles[i].atom1 = atom1;
    angles[i].atom2 = atom2;
    angles[i].atom3 = atom3;

    angles[i].angle_type=type;
  }

  numExclusions =  0;
  exclusions = new Exclusion[numExclusions];

  /*
  // If readExclusions is TRUE, then we copy exclusions from parm
  // file; otherwise we skip the exclusions here and generate
  // them later in build_exclusions()
  if (simParams->readExclusions)
  { // Copy exclusion information
    // In Amber data structure, Iblo[] is the number of exclusions
    // for each atom; ExclAt[] is the atom index for the excluded atoms.
    exclusions = new Exclusion[amber_data->Nnb];
    if (exclusions == NULL &&  amber_data->Nnb > 0)
      NAMD_die("memory allocation failed when reading exclusion information");
    current_index = 0;
    for (i=0; i<numAtoms; ++i)
      for (j=0; j<amber_data->Iblo[i]; ++j)
      { if (current_index >= amber_data->Nnb)
        { char err_msg[128];
          sprintf(err_msg, "EXCLUSION INDEX EXCEEDS NUMBER OF EXLCUSIONS %d IN AMBER FILE, AT ATOM #%d\n",
            amber_data->Nnb, i+1);
          NAMD_die(err_msg);
        }
        // There's some 0 in the ExclAt[] list, which is strange
        // and redundant. In this case, I simply ignore such entries.
        if (amber_data->ExclAt[current_index] != 0)
        { // Subtract 1 to convert the index from the 1 to NumAtoms
          // used in the file to the 0 to NumAtoms-1 that we need
          a2 = amber_data->ExclAt[current_index] - 1;
          if (a2 < i)
          { // I assume the latter index be larger than the former
            // one, so that the same exclusion won't be double-counted;
            // if not, give error
            char err_msg[128];
            sprintf(err_msg, "Atom #%d has exclusion with atom #%d, in reverse order.", i+1, a2+1);
            NAMD_die(err_msg);
          }
          else if (a2 == i)
          { char err_msg[128];
            sprintf(err_msg, "ATOM %d EXCLUDED FROM ITSELF IN AMBER FILE\n", i+1);
              NAMD_die(err_msg);
          }
          else if (a2 >= numAtoms)
          {  char err_msg[128];
             sprintf(err_msg, "EXCLUSION INDEX %d GREATER THAN NATOM %d IN EXCLUSION # %d IN AMBER FILE",
               a2+1, numAtoms, current_index+1);
             NAMD_die(err_msg);
          }
          exclusions[numExclusions].atom1 = i;
          exclusions[numExclusions].atom2 = a2;
          ++numExclusions;
        }
        ++current_index;
      }
    if (current_index < amber_data->Nnb)
    { char err_msg[128];
      sprintf(err_msg, "Num of exclusions recorded (%d) is smaller than what it's supposed to be (%d)",
        current_index,amber_data->Nnb);
      NAMD_die(err_msg);
    }
  }
  */

  // Copy dihedral information
  numDihedrals = gf->getNumDihedrals();
  dihedrals = new Dihedral[numDihedrals];
  if (dihedrals == NULL)
    NAMD_die("memory allocation failed when reading dihedral information");
  for(i=0;i<numDihedrals;i++) {
    int type; // to convert the type correctly
    int atom1,atom2,atom3,atom4;
    gf->getDihedral(i,&atom1,&atom2,&atom3,&atom4,&type);
    dihedrals[i].atom1 = atom1;
    dihedrals[i].atom2 = atom2;
    dihedrals[i].atom3 = atom3;
    dihedrals[i].atom4 = atom4;
    dihedrals[i].dihedral_type = type;
  }

  /*
  // In AMBER parm file, dihedrals contain 1-4 exclusion infomation:
  // the 1st and 4th atoms have 1-4 nonbond interation. So we should
  // find them in the exclusion array and change their exclusion to
  // 1-4 type. However, there're two exceptions --
  // 1.If the third atom is negative, it means the end group
  //   interactions are to be ignored;
  // 2.If the fourth atom is negative, it means this is an improper.
  // For the above two cases, the actual atom index is the absolute
  // value of the atom number read; and there's no 1-4 interation
  // for these dihedrals.
  // If readExclusions is not TRUE, then we don't worry about
  // exclusions here.
  for (i=0; i<numDihedrals; ++i)
  { if (dihedrals[i].atom3 < 0 || dihedrals[i].atom4 < 0)
    { dihedrals[i].atom3 = abs(dihedrals[i].atom3);
      dihedrals[i].atom4 = abs(dihedrals[i].atom4);
    }
    else if (simParams->readExclusions)
    { if (dihedrals[i].atom1 < dihedrals[i].atom4)
        a1=dihedrals[i].atom1, a2=dihedrals[i].atom4;
      else
        a1=dihedrals[i].atom4, a2=dihedrals[i].atom1;
      // Since in the exclusion array, atom1 is guaranteed to be
      // ordered, we can do a binary serch to find it first.
      found = 0;
      min=0, max=numExclusions-1;
      while (!found && min<=max)
      { index = (min+max)/2;
        if (exclusions[index].atom1 == a1)
          found = 1;
        else if (exclusions[index].atom1 < a1)
          min = index+1;
        else
          max = index-1;
      }
      if (!found)
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
      // After finding atom1, we do a linear serch to find atom2,
      // in both directions.
      for (j=index-1; j>=0 && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; --j);
      if (j<0 || exclusions[j].atom1!=a1)
        for (j=index; j<numExclusions && exclusions[j].atom2!=a2 && exclusions[j].atom1==a1; ++j);
      if (j<numExclusions && exclusions[j].atom1==a1)
        exclusions[j].modified = 1;  // Change the exclusion type to 1-4
      else
        NAMD_die("1-4 interaction in dihedral not found in exclusion list!");
    }
    if (dihedrals[i].atom1>=numAtoms || dihedrals[i].atom2>=numAtoms ||
        dihedrals[i].atom3>=numAtoms || dihedrals[i].atom4>=numAtoms ||
        dihedrals[i].dihedral_type>=amber_data->Nptra)
    { char err_msg[128];
      sprintf(err_msg, "DIHEDRAL # %d OVERFLOW IN PARM FILE", i+1);
      NAMD_die(err_msg);
    }
  }
  */
  //  analyze the data and find the status of each atom

  build_atom_status();
#endif
}
/*      END OF FUNCTION read_parm    */

#ifndef MEM_OPT_VERSION
/*
int32 *Molecule::get_bonds_for_atom(int anum){
    NAMD_die("In bonds for atom!");
    return bondsByAtom[anum];
}

Bond *Molecule::get_bond(int bnum){
    NAMD_die("In get_bond!");
    return &bonds[bnum];
}
*/
#endif

#ifdef MEM_OPT_VERSION

void Molecule::delOtherEachAtomStructs(){
    if(CmiMyRank()) return;

    delete [] fixedAtomFlags;
    fixedAtomFlags = NULL;

    //decide whether to free space for hydrogenGroup and atoms fields
    //the condition can be referred to the comment before wrap_coor_int
    //in Output.C
    int peOfCollectionMaster = 0;
    if(CkNumPes()>1 && simParams->shiftIOToOne) peOfCollectionMaster = 1;
    if(CkNumPes()==1 && simParams->wrapWater && !simParams->wrapAll){
	//considering namd runs on a single processor
	return;
    }else if(CkMyPe()==peOfCollectionMaster){
	if(simParams->wrapAll || !simParams->wrapWater){
	    hydrogenGroup.resize(0);
	    delete [] atoms;
	    atoms = NULL;
	} 
    }else{	      
	hydrogenGroup.resize(0);
	delete [] atoms;
	atoms = NULL;
    }
}

//return the index of the new mass in the mass pool
Index Molecule::insert_new_mass(Real newMass){
    //first search
    for(int i=massPoolSize-1; i>=0; i--){
        if(fabs(atomMassPool[i]-newMass)<=1e-6)
            return i;
    }
    //otherwise increase one entry for the new mass
    Real *tmp = new Real[massPoolSize+1];
    tmp[massPoolSize] = newMass;
    memcpy((void *)tmp, (const void *)atomMassPool, sizeof(Real)*massPoolSize);
    delete [] atomMassPool;
    atomMassPool = tmp;
    massPoolSize++;
    return (Index)(massPoolSize-1);
}

/*void Molecule::markClusterIdx(int curClusterIdx, int startAtomID){
    if(cluster[startAtomID] == curClusterIdx) return;
    AtomSignature *sig = &atomSigPool[startAtomID];
    for(int i=0; i<sig->bondCnt; i++){
        TupleSignature *tsig = &(sig->bondSigs[i]);
        int newStartAtomID = startAtomID+tsig->offset[0];
        markClusterIdx(curClusterIdx, newStartAtomID);
    }
}*/

//if the exclusion indicated by atom1 and atom2 will be stripped
//by stripFepExcl function or fixedAtoms, then return 1. Otherwise return 0
int Molecule::exclStrippedByFepOrFixedAtoms(int atom1, int atom2){
    //first test for stripFepExcl
    if(simParams->alchFepOn || simParams->alchThermIntOn || simParams->lesOn){
        int t1 = get_fep_type(atom1);
        int t2 = get_fep_type(atom2);
        if(t1 && t2 && t1!=t2) return 1;
    }else if(simParams->pairInteractionOn){
        int t1 = get_fep_type(atom1);
        int t2 = get_fep_type(atom2);
        if(simParams->pairInteractionSelf){
            if(t1!=1 || t2!=1) return 1;
        }else{
            if((t1!=1 || t2!=2) && (t1!=2 || t2!=1))
                return 1;
        }
    }

    if(numFixedAtoms && fixedAtomFlags[atom1] && fixedAtomFlags[atom2])
        return 1;

    return 0;
}

void Molecule::addNewExclSigPool(const vector<ExclusionSignature>& newExclSigPool){
    ExclusionSignature *tmpExclSigPool = new ExclusionSignature[exclSigPoolSize+newExclSigPool.size()];
    for(int i=0; i<exclSigPoolSize; i++)
        tmpExclSigPool[i] = exclSigPool[i];
    for(int i=0; i<newExclSigPool.size(); i++)
        tmpExclSigPool[i+exclSigPoolSize] = newExclSigPool[i];

    exclSigPoolSize += newExclSigPool.size();
    exclSigPool = tmpExclSigPool;
}
#endif

void TupleSignature::pack(MOStream *msg){
    msg->put((short)tupleType);
    msg->put(numOffset);
    msg->put(numOffset, offset);    
    msg->put(tupleParamType);
    msg->put(isReal);
}

void TupleSignature::unpack(MIStream *msg){
    short ttype;
    msg->get(ttype);
    tupleType = (TupleSigType)ttype;

    msg->get(numOffset);
    delete [] offset;
    offset = new int[numOffset];
    msg->get(numOffset*sizeof(int), (char *)offset);
    
    msg->get(tupleParamType);
    msg->get(isReal);
}

void AtomSignature::pack(MOStream *msg){
    msg->put(bondCnt);
    for(int i=0; i<bondCnt; i++)
        bondSigs[i].pack(msg);

    msg->put(angleCnt);
    for(int i=0; i<angleCnt; i++)
        angleSigs[i].pack(msg);

    msg->put(dihedralCnt);
    for(int i=0; i<dihedralCnt; i++)
        dihedralSigs[i].pack(msg);

    msg->put(improperCnt);
    for(int i=0; i<improperCnt; i++)
        improperSigs[i].pack(msg);

    msg->put(crosstermCnt);
    for(int i=0; i<crosstermCnt; i++)
        crosstermSigs[i].pack(msg);
}

void AtomSignature::unpack(MIStream *msg){
    msg->get(bondCnt);
    delete [] bondSigs;
    if(bondCnt>0){
        bondSigs = new TupleSignature[bondCnt];
        for(int i=0; i<bondCnt; i++)
            bondSigs[i].unpack(msg);
    } else bondSigs = NULL;

    msg->get(angleCnt);
    delete [] angleSigs;
    if(angleCnt>0){
        angleSigs = new TupleSignature[angleCnt];
        for(int i=0; i<angleCnt; i++)
            angleSigs[i].unpack(msg);
    } else angleSigs = NULL;

    msg->get(dihedralCnt);
    delete [] dihedralSigs;
    if(dihedralCnt>0){
        dihedralSigs = new TupleSignature[dihedralCnt];
        for(int i=0; i<dihedralCnt; i++)
            dihedralSigs[i].unpack(msg);
    } else dihedralSigs = NULL;

    msg->get(improperCnt);
    delete [] improperSigs;
    if(improperCnt>0){
        improperSigs = new TupleSignature[improperCnt];
        for(int i=0; i<improperCnt; i++)
            improperSigs[i].unpack(msg);
    } else improperSigs = NULL;

    msg->get(crosstermCnt);
    delete [] crosstermSigs;
    if(crosstermCnt>0){
        crosstermSigs = new TupleSignature[crosstermCnt];
        for(int i=0; i<crosstermCnt; i++)
            crosstermSigs[i].unpack(msg);
    } else crosstermSigs = NULL;
}

void AtomSignature::removeEmptyTupleSigs(){
    int origTupleCnt;
    int idx;
    TupleSignature *tupleSigs;
    TupleSignature *newTupleSigs;

    //bonds
    origTupleCnt = bondCnt;
    tupleSigs= bondSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            bondCnt--;
    }
    if(bondCnt==0){
        delete [] tupleSigs;
        bondSigs = NULL;
    }else if(bondCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[bondCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        bondSigs = newTupleSigs;
    }

    //angles
    origTupleCnt = angleCnt;
    tupleSigs = angleSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            angleCnt--;
    }
    if(angleCnt==0){
        delete [] tupleSigs;
        angleSigs = NULL;
    }else if(angleCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[angleCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        angleSigs = newTupleSigs;
    }

    //dihedrals
    origTupleCnt = dihedralCnt;
    tupleSigs = dihedralSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            dihedralCnt--;
    }
    if(dihedralCnt==0){
        delete [] tupleSigs;
        dihedralSigs = NULL;
    }else if(dihedralCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[dihedralCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        dihedralSigs = newTupleSigs;        
    }


    //impropers
    origTupleCnt = improperCnt;
    tupleSigs = improperSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            improperCnt--;
    }
    if(improperCnt==0){
        delete [] tupleSigs;
        improperSigs = NULL;
    }else if(improperCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[improperCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        improperSigs = newTupleSigs;
    }    

    //crossterms
    origTupleCnt = crosstermCnt;
    tupleSigs = crosstermSigs;
    for(int i=0; i<origTupleCnt; i++){
        if(tupleSigs[i].isEmpty())
            crosstermCnt--;
    }
    if(crosstermCnt==0){
        delete [] tupleSigs;
        crosstermSigs = NULL;
    }else if(crosstermCnt!=origTupleCnt){
        newTupleSigs = new TupleSignature[crosstermCnt];
        idx=0;
        for(int i=0; i<origTupleCnt; i++){
            if(!tupleSigs[i].isEmpty()){
                newTupleSigs[idx] = tupleSigs[i];
                idx++;
            }
        }
        delete [] tupleSigs;
        crosstermSigs = newTupleSigs;
    }    
}

void ExclusionSignature::removeEmptyOffset(){
	int newCnt=0;
	for(int i=0; i<fullExclCnt; i++){
	    if(fullOffset[i]==0) continue;
	    newCnt++;
	}
    if(newCnt==0){
        fullExclCnt = 0;
        delete [] fullOffset;
        fullOffset = NULL;
    }else if(newCnt!=fullExclCnt){
        int *tmpOffset = new int[newCnt];
    	newCnt=0;
    	for(int i=0; i<fullExclCnt; i++){
    	    if(fullOffset[i]==0) continue;
    	    tmpOffset[newCnt] = fullOffset[i];
    	    newCnt++;
    	}
    	delete [] fullOffset;
    	fullOffset = tmpOffset;
        fullExclCnt = newCnt;
    }
	
	
	newCnt=0;
	for(int i=0; i<modExclCnt; i++){
	    if(modOffset[i]==0) continue;
	    newCnt++;
	}
    if(newCnt==0){
        modExclCnt = 0;
        delete [] modOffset;
        modOffset = NULL;
    }else if(newCnt!=modExclCnt){
        int *tmpOffset = new int[newCnt];
        newCnt=0;
        for(int i=0; i<modExclCnt; i++){
            if(modOffset[i]==0) continue;
            tmpOffset[newCnt] = modOffset[i];
            newCnt++;
        }
        delete [] modOffset;
        modOffset = tmpOffset;
        modExclCnt = newCnt;
    }	
}

//returns the index of the offset. If not found, -1 is returned
//fullOrMod indicates where is the offset found. 0 indicates in
//the full exclusion lists, 1 indicates in the modified exclusion
//lists
int ExclusionSignature::findOffset(int offset, int *fullOrMod){
	//assuming all offsets have been sorted increasingly
	//so that binary search could be used	
	int retidx = -1;
	
	*fullOrMod = 0;	
	int low = 0;
	int high = fullExclCnt-1;
	int mid = (low+high)/2;
	while(low<=high){
		if(offset<fullOffset[mid]){
			high = mid-1;
			mid = (high+low)/2;						
		}else if(offset>fullOffset[mid]){
			low = mid+1;
			mid = (high+low)/2;
		}else{
			retidx = mid;
			break;
		}		
	}
	if(retidx!=-1) return retidx;
	
	*fullOrMod = 1;	
	low = 0;
	high = modExclCnt-1;
	mid = (low+high)/2;
	while(low<=high){
		if(offset<modOffset[mid]){
			high = mid-1;
			mid = (high+low)/2;						
		}else if(offset>modOffset[mid]){
			low = mid+1;
			mid = (high+low)/2;
		}else{
			retidx = mid;
			break;
		}		
	}
	return retidx;	
}

void ExclusionSignature::pack(MOStream *msg){
    msg->put(fullExclCnt);    
    msg->put(fullExclCnt, fullOffset);
    msg->put(modExclCnt);    
    msg->put(modExclCnt, modOffset);
}

void ExclusionSignature::unpack(MIStream *msg){    
    msg->get(fullExclCnt);
    delete [] fullOffset;
    fullOffset = new int[fullExclCnt];
    msg->get(fullExclCnt*sizeof(int), (char *)fullOffset);    
    msg->get(modExclCnt);
    delete [] modOffset;
    modOffset = new int[modExclCnt];
    msg->get(modExclCnt*sizeof(int), (char *)modOffset);    
}

#endif  // MOLECULE2_C defined = second object file

