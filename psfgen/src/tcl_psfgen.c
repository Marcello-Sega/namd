
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "psfgen.h"
#include "charmm_parse_topo_defs.h"
#include "topo_mol_output.h"
#include "pdb_file_extract.h"
#include "psf_file_extract.h"
#include "extract_alias.h"

#if defined(NAMD_TCL) || ! defined(NAMD_VERSION)

#include <tcl.h>

/* 
 * Provide user feedback and warnings beyond result values.
 * If we are running interactively, Tcl_Main will take care of echoing results
 * to the console.  If we run a script, we need to output the results
 * ourselves.
 */
void newhandle_msg(void *v, const char *msg) {
  Tcl_Interp *interp = (Tcl_Interp *)v;
  const char *words[2] = {"puts"};
  char *script;
  words[1] = msg;
  script = Tcl_Merge(2,words);
  Tcl_Eval(interp,script);
  Tcl_Free(script);
}

/*
 * Kills molecule to prevent user from saving bogus output.
 */
void psfgen_kill_mol(Tcl_Interp *interp, psfgen_data *data) {
  if (data->mol) {
    Tcl_AppendResult(interp,
	"\nMOLECULE DESTROYED BY FATAL ERROR!  Use resetpsf to start over.",
	NULL);
  }
  topo_mol_destroy(data->mol);
  data->mol = 0;
}

/* This function gets called if/when the Tcl interpreter is deleted. */
static void psfgen_deleteproc(ClientData cd, Tcl_Interp *interp) {
  psfgen_data *data = (psfgen_data *)cd;
  topo_mol_destroy(data->mol);
  topo_defs_destroy(data->defs);
  stringhash_destroy(data->aliases);
  free(data);
}

void psfgen_data_delete_pointer(ClientData cd, Tcl_Interp *interp) {
  psfgen_data **dataptr = (psfgen_data **)cd;
  free(dataptr);
}

static void count_delete_proc(ClientData data, Tcl_Interp *interp) {
  free(data);
}

psfgen_data* psfgen_data_create(Tcl_Interp *interp) {
  char namebuf[128];
  int *countptr;
  int id;
  psfgen_data *data;
  countptr = Tcl_GetAssocData(interp, "Psfgen_count", 0);
  if (!countptr) {
    countptr = (int *)malloc(sizeof(int));
    Tcl_SetAssocData(interp, "Psfgen_count", count_delete_proc, 
      (ClientData)countptr);
    *countptr = 0;
  } 
  id = *countptr;
  data = (psfgen_data *)malloc(sizeof(psfgen_data));
  data->defs = topo_defs_create();
  topo_defs_error_handler(data->defs,interp,newhandle_msg);
  data->aliases = stringhash_create();
  data->mol = topo_mol_create(data->defs);
  topo_mol_error_handler(data->mol,interp,newhandle_msg);
  data->id = id;
  *countptr = id+1;
  sprintf(namebuf,"Psfgen_%d",id);
  Tcl_SetAssocData(interp,namebuf,psfgen_deleteproc,(ClientData)data);
  return data;
}

int tcl_psfcontext(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_topology(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_segment(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_residue(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_mutate(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_multiply(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_coord(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_auto(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_alias(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_pdb(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_coordpdb(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_guesscoord(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_readpsf(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_writepsf(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_writepdb(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_first(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_last(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_patch(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_resetpsf(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_delatom(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);

#if defined(PSFGENTCLDLL_EXPORTS) && defined(_WIN32)
#  undef TCL_STORAGE_CLASS
#  define TCL_STORAGE_CLASS DLLEXPORT

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Window s headers
#include <windows.h>

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
                                         )
{
    return TRUE;
}

EXTERN int Psfgen_Init(Tcl_Interp *interp) {

#else

int Psfgen_Init(Tcl_Interp *interp) {

#endif

  /* Create psfgen data structures; keep in interp so that other libraries
   * can access them.
   */
  psfgen_data **data, *data_0;
  Tcl_SetAssocData(interp, (char *)"Psfgen_count",0,(ClientData)0);
  data = (psfgen_data **)malloc(sizeof(psfgen_data *));
  Tcl_SetAssocData(interp, (char *)"Psfgen_pointer",
		psfgen_data_delete_pointer,(ClientData)data);
  *data = data_0 = psfgen_data_create(interp);

  Tcl_CreateCommand(interp,"psfcontext",tcl_psfcontext,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"topology",tcl_topology,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"readpsf",tcl_readpsf,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"segment",tcl_segment,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"residue",tcl_residue,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"mutate",tcl_mutate,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"multiply",tcl_multiply,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"coord",tcl_coord,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"auto",tcl_auto,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"alias",tcl_alias,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"pdb",tcl_pdb,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"coordpdb",tcl_coordpdb,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"guesscoord",tcl_guesscoord,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"writepsf",tcl_writepsf,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"writepdb",tcl_writepdb,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"first",tcl_first,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"last",tcl_last,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"patch",tcl_patch,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"resetpsf", tcl_resetpsf,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"delatom", tcl_delatom,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
 
  Tcl_PkgProvide(interp, "psfgen", "1.2");

  return TCL_OK;
}

void strtoupper(char *s) {
  while ( *s ) { *s = toupper(*s); ++s; }
}

char* splitcolon(char *s) {
  if ( s ) {
    while ( *s && *s != ':' ) { ++s; }
    if ( *s ) *(s++) = 0; else s = 0;
  }
  return s;
}

int tcl_psfcontext(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {

  int oldid, newid;
  int delold = 0;
  psfgen_data **cur = (psfgen_data **)data;
  char oldidstr[128];
  oldid = (*cur)->id;
  sprintf(oldidstr,"%d",oldid);

  if ( argc == 1 ) {
    Tcl_SetResult(interp,oldidstr,TCL_VOLATILE);
    return TCL_OK;
  }

  if ( argc == 3 ) {
    if ( strcmp(argv[2],"delete") == 0 ) {
      delold = 1;
    } else {
      Tcl_SetResult(interp,"second argument must be delete",TCL_VOLATILE);
      psfgen_kill_mol(interp,*cur);
      return TCL_ERROR;
    }
  }

  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,*cur);
    return TCL_ERROR;
  }

  if (strcmp(argv[1],"new") == 0) {
    psfgen_data *newdata = psfgen_data_create(interp);
    *cur = newdata;
  } else if (Tcl_GetInt(interp,argv[1],&newid) == TCL_OK) {
    psfgen_data *newdata;
    char newkey[128];
    if ( newid == oldid ) {
      if ( delold ) {
        Tcl_SetResult(interp,"specified context in use",TCL_VOLATILE);
        psfgen_kill_mol(interp,*cur);
        return TCL_ERROR;
      } else {
        Tcl_SetResult(interp,oldidstr,TCL_VOLATILE);
        return TCL_OK;
      }
    }
    sprintf(newkey,"Psfgen_%d",newid);
    if ( newdata = Tcl_GetAssocData(interp,newkey,0) ) {
      *cur = newdata;
    } else {
      Tcl_SetResult(interp,"specified context does not exist",TCL_VOLATILE);
      psfgen_kill_mol(interp,*cur);
      return TCL_ERROR;
    }
  } else {
    Tcl_SetResult(interp,"first argument must be existing context or new",TCL_VOLATILE);
    psfgen_kill_mol(interp,*cur);
    return TCL_ERROR;
  }

  if ( delold ) {
    char oldkey[128];
    sprintf(oldkey,"Psfgen_%d",oldid);
    Tcl_DeleteAssocData(interp,oldkey);
    sprintf(oldkey,"deleted %d",oldid);
    Tcl_SetResult(interp,oldkey,TCL_VOLATILE);
    return TCL_OK;
  } else {
    Tcl_SetResult(interp,oldidstr,TCL_VOLATILE);
    return TCL_OK;
  }

}

int tcl_topology(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  FILE *defs_file;
  char *filename;
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc == 1 ) {
    Tcl_SetResult(interp,"no topology file specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( argc > 2 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  filename = argv[1];
  if ( ! ( defs_file = fopen(filename,"r") ) ) {
    sprintf(msg,"ERROR: Unable to open topology file %s\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  } else {
    sprintf(msg,"reading topology file %s\n",filename);
    newhandle_msg(interp,msg);
    charmm_parse_topo_defs(psf->defs,defs_file,interp,newhandle_msg);
    fclose(defs_file);
  }
  return TCL_OK;
}

int tcl_readpsf(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  FILE *psf_file;
  int retval;
  char *filename;
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc == 1 ) {
    Tcl_SetResult(interp,"no psf file specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( argc > 2 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  filename = argv[1];
  if ( ! ( psf_file = fopen(filename,"r") ) ) {
    sprintf(msg,"ERROR: Unable to open psf file %s\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  } else {
    sprintf(msg,"reading structure from psf file %s\n",filename);
    newhandle_msg(interp,msg);
    retval = psf_file_extract(psf->mol, psf_file, interp, newhandle_msg);
    fclose(psf_file);
  }
  if (retval) {
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  return TCL_OK;
}

int tcl_segment(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 3 ) {
    Tcl_SetResult(interp,"arguments: segname { commmands }",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);

  sprintf(msg,"building segment %s",argv[1]);
  newhandle_msg(interp,msg);
  if ( topo_mol_segment(psf->mol,argv[1]) ) {
    Tcl_AppendResult(interp,"ERROR: failed on segment",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  if ( Tcl_Eval(interp,argv[2]) != TCL_OK ) {
    Tcl_AppendResult(interp,"\nERROR: failed while building segment",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  newhandle_msg(interp,"generating structure at end of segment");
  if ( topo_mol_end(psf->mol) ) {
    Tcl_AppendResult(interp,"ERROR: failed on end of segment",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_residue(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 3 ) {
    Tcl_SetResult(interp,"arguments: resid resname",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);
  strtoupper(argv[2]);

  if ( topo_mol_residue(psf->mol,argv[1],argv[2]) ) {
    Tcl_AppendResult(interp,"ERROR: failed on residue",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_mutate(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 3 ) {
    Tcl_SetResult(interp,"arguments: resid resname",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);
  strtoupper(argv[2]);

  if ( topo_mol_mutate(psf->mol,argv[1],argv[2]) ) {
    Tcl_AppendResult(interp,"ERROR: failed on mutate",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_multiply(ClientData data, Tcl_Interp *interp,
                                        int argc, char *argv[]) {
  int i, ncopies, ierr;
  topo_mol_ident_t *targets;
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc<3 || Tcl_GetInt(interp,argv[1],&ncopies) != TCL_OK || ncopies<2 ) {
    Tcl_SetResult(interp,"arguments: ncopies segid?:resid?:atomname? ...",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  targets = (topo_mol_ident_t *) Tcl_Alloc(argc*sizeof(topo_mol_ident_t));
  if ( ! targets ) {
    Tcl_SetResult(interp,"memory allocation failed",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  sprintf(msg,"generating %d copies of selected atoms",ncopies);
  newhandle_msg(interp,msg);
  for ( i=2; i<argc; ++i ) {
    char *ctmp;
    strtoupper(argv[i]);
    targets[i-2].segid = ctmp = argv[i];
    targets[i-2].resid = ctmp = splitcolon(ctmp);
    targets[i-2].aname = splitcolon(ctmp);
  }
  if ( ierr = topo_mol_multiply_atoms(psf->mol,targets,(argc-2),ncopies) ) {
    sprintf(msg,"ERROR: failed to multiply atoms (error=%d)",ierr);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    /* Tcl_AppendResult(interp,"ERROR: failed to multiply atoms",NULL); */
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_coord(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  double x,y,z;
  topo_mol_ident_t target;
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 5 ) {
    Tcl_SetResult(interp,"arguments: segid resid atomname { x y z }",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( argc > 5 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( sscanf(argv[4],"%lf %lf %lf",&x,&y,&z) != 3 ) {
    Tcl_SetResult(interp,"arguments: segid resid atomname { x y z }",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);
  strtoupper(argv[2]);
  strtoupper(argv[3]);

  target.segid = argv[1];
  target.resid = argv[2];
  target.aname = argv[3];
  if ( topo_mol_set_xyz(psf->mol,&target,x,y,z) ) {
    Tcl_AppendResult(interp,"ERROR: failed on coord",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  return TCL_OK;
}


int tcl_auto(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  int i, angles, dihedrals;
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 2 ) {
    Tcl_SetResult(interp,"arguments: ?angles? ?dihedrals? ?none?",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  angles = 0;  dihedrals = 0;
  for ( i = 1; i < argc; ++i ) {
    if ( ! strcmp(argv[i],"angles") ) angles = 1;
    else if ( ! strcmp(argv[i],"dihedrals") ) dihedrals = 1;
    else if ( strcmp(argv[i],"none") ) {
      Tcl_SetResult(interp,"arguments: ?angles? ?dihedrals? ?none?",TCL_VOLATILE);
      psfgen_kill_mol(interp,psf);
      return TCL_ERROR;
    }
  }

  if ( angles ) newhandle_msg(interp,"enabling angle autogeneration");
  else newhandle_msg(interp,"disabling angle autogeneration");
  if ( topo_mol_segment_auto_angles(psf->mol,angles) ) {
    Tcl_AppendResult(interp,"ERROR: failed setting angle autogen",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  if ( dihedrals ) newhandle_msg(interp,"enabling dihedral autogeneration");
  else newhandle_msg(interp,"disabling dihedral autogeneration");
  if ( topo_mol_segment_auto_dihedrals(psf->mol,dihedrals) ) {
    Tcl_AppendResult(interp,"ERROR: failed setting dihedral autogen",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_alias(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 2 ) {
    Tcl_SetResult(interp,"arguments: atom | residue ...",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  if ( ! strcmp(argv[1],"residue") ) {
    if ( argc < 4 ) {
      Tcl_SetResult(interp,"arguments: residue altres realres",TCL_VOLATILE);
      psfgen_kill_mol(interp,psf);
      return TCL_ERROR;
    }
    strtoupper(argv[2]);
    strtoupper(argv[3]);
    sprintf(msg,"aliasing residue %s to %s",argv[2],argv[3]);
    newhandle_msg(interp,msg);
    if ( extract_alias_residue_define(psf->aliases,argv[2],argv[3]) ) {
      Tcl_AppendResult(interp,"ERROR: failed on residue alias",NULL);
      psfgen_kill_mol(interp,psf);
      return TCL_ERROR;
    }
  } else if ( ! strcmp(argv[1],"atom") ) {
    if ( argc < 5 ) {
      Tcl_SetResult(interp,"arguments: atom resname altatom realatom",TCL_VOLATILE);
      psfgen_kill_mol(interp,psf);
      return TCL_ERROR;
    }
    strtoupper(argv[2]);
    strtoupper(argv[3]);
    strtoupper(argv[4]);
    sprintf(msg,"aliasing residue %s atom %s to %s",argv[2],argv[3],argv[4]);
    newhandle_msg(interp,msg);
    if ( extract_alias_atom_define(psf->aliases,argv[2],argv[3],argv[4]) ) {
      Tcl_AppendResult(interp,"ERROR: failed on atom alias",NULL);
      psfgen_kill_mol(interp,psf);
      return TCL_ERROR;
    }
  }

  return TCL_OK;
}

int tcl_pdb(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  FILE *res_file;
  char *filename;
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc == 1 ) {
    Tcl_SetResult(interp,"no pdb file specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( argc > 2 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  filename = argv[1];
  if ( ! ( res_file = fopen(filename,"r") ) ) {
    sprintf(msg,"ERROR: Unable to open pdb file %s to read residues\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  } else {
    sprintf(msg,"reading residues from pdb file %s",filename);
    newhandle_msg(interp,msg);
    if ( pdb_file_extract_residues(psf->mol,res_file,psf->aliases,interp,newhandle_msg) ) {
      Tcl_AppendResult(interp,"ERROR: failed on reading residues from pdb file",NULL);
      fclose(res_file);
      psfgen_kill_mol(interp,psf);
      return TCL_ERROR;
    }
    fclose(res_file);
  }

  return TCL_OK;
}

int tcl_coordpdb(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  FILE *res_file;
  char *filename;
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 2 ) {
    Tcl_SetResult(interp,"arguments: pdbfile ?segid?",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }
  filename = argv[1];
  if ( ! ( res_file = fopen(filename,"r") ) ) {
    sprintf(msg,"ERROR: Unable to open pdb file %s to read coordinates\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  } else {
    const char *segid;
    if (argc == 3) {
      /* Read only coordinates for given segid */
      strtoupper(argv[2]);
      sprintf(msg,"reading coordinates from pdb file %s for segment %s",filename,argv[2]);
      newhandle_msg(interp,msg);
      segid = argv[2];
    } else {
      /* Read all segid's in pdb file */
      sprintf(msg,"reading coordinates from pdb file %s",filename);
      newhandle_msg(interp,msg);
      segid = NULL;
    } 
    if ( pdb_file_extract_coordinates(psf->mol,res_file,segid,psf->aliases,interp,newhandle_msg) ) {
      Tcl_AppendResult(interp,"ERROR: failed on reading coordinates from pdb file",NULL);
      fclose(res_file);
      psfgen_kill_mol(interp,psf);
      return TCL_ERROR;
    }
    fclose(res_file);
  }

  return TCL_OK;

}

int tcl_guesscoord(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  psfgen_data *psf = *(psfgen_data **)data;
  if ( argc > 1 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  newhandle_msg(interp,"guessing coordinates based on topology file");
  if ( topo_mol_guess_xyz(psf->mol) ) {
    Tcl_AppendResult(interp,"ERROR: failed on guessing coordinates",NULL);
    psfgen_kill_mol(interp,psf);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_writepsf(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  FILE *res_file;
  char *filename;
  int charmmfmt;
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc == 1 ) {
    Tcl_SetResult(interp,"no psf file specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  charmmfmt = 0;
  if ( argc == 3 ) {
    if ( strcmp(argv[1],"charmm") == 0 ) charmmfmt = 1;
    else if ( strcmp(argv[1],"x-plor") == 0 ) charmmfmt = 0;
    else {
      sprintf(msg,"ERROR: Unknown psf file format %s (not charmm or x-plor).\n",argv[1]);
      Tcl_SetResult(interp,msg,TCL_VOLATILE);
      return TCL_ERROR;
    }
  }
  filename = argv[argc-1];

  if ( ! ( res_file = fopen(filename,"w") ) ) {
    sprintf(msg,"ERROR: Unable to open psf file %s to write structure\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    return TCL_ERROR;
  }
  sprintf(msg,"writing psf file %s",filename);
  newhandle_msg(interp,msg);
  if ( topo_mol_write_psf(psf->mol,res_file,charmmfmt,interp,newhandle_msg) ) {
    Tcl_AppendResult(interp,"ERROR: failed on writing structure to psf file",NULL);
    fclose(res_file);
    return TCL_ERROR;
  }
  fclose(res_file);

  return TCL_OK;
}

int tcl_writepdb(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  FILE *res_file;
  char *filename;
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc == 1 ) {
    Tcl_SetResult(interp,"no pdb file specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 2 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  filename = argv[1];

  if ( ! ( res_file = fopen(filename,"w") ) ) {
    sprintf(msg,"ERROR: Unable to open pdb file %s to write coordinates\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    return TCL_ERROR;
  }
  sprintf(msg,"writing pdb file %s",filename);
  newhandle_msg(interp,msg);
  if ( topo_mol_write_pdb(psf->mol,res_file,interp,newhandle_msg) ) {
    Tcl_AppendResult(interp,"ERROR: failed on writing coordinates to pdb file",NULL);
    fclose(res_file);
    return TCL_ERROR;
  }
  fclose(res_file);

  return TCL_OK;
}

int tcl_first(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc != 2 ) {
    Tcl_SetResult(interp,"argument: presname",TCL_VOLATILE);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);

  sprintf(msg,"setting patch for first residue to %s",argv[1]);
  newhandle_msg(interp,msg);
  if ( topo_mol_segment_first(psf->mol,argv[1]) ) {
    Tcl_AppendResult(interp,"ERROR: failed to set patch for first residue",NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_last(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc != 2 ) {
    Tcl_SetResult(interp,"argument: presname",TCL_VOLATILE);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);

  sprintf(msg,"setting patch for last residue to %s",argv[1]);
  newhandle_msg(interp,msg);
  if ( topo_mol_segment_last(psf->mol,argv[1]) ) {
    Tcl_AppendResult(interp,"ERROR: failed to set patch for last residue",NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_patch(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  int i;
  topo_mol_ident_t targets[10];
  char msg[128];
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 3 ) {
    Tcl_SetResult(interp,"arguments: presname segid:resid ...",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 10 ) {
    Tcl_SetResult(interp,"too many targets for patch",TCL_VOLATILE);
    return TCL_ERROR;
  }

  strtoupper(argv[1]);
  sprintf(msg,"applying patch %s to %d residues",argv[1],(argc-2));
  newhandle_msg(interp,msg);
  for ( i=2; i<argc; ++i ) {
    strtoupper(argv[i]);
    targets[i-2].segid = argv[i];
    targets[i-2].resid = splitcolon(argv[i]);
    targets[i-2].aname = 0;
    if ( ! targets[i-2].resid ) {
      sprintf(msg,"ERROR: resid missing from patch target %s",argv[i]);
      Tcl_SetResult(interp,msg,TCL_VOLATILE);
      return TCL_ERROR;
    }
  }
  if ( topo_mol_patch(psf->mol,targets,(argc-2),argv[1],0,0,0) ) {
    Tcl_AppendResult(interp,"ERROR: failed to apply patch",NULL);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_resetpsf(ClientData data, Tcl_Interp *interp, int argc, char *argv[]) {
  psfgen_data *psf = *(psfgen_data **)data;

  topo_mol_destroy(psf->mol);
  psf->mol = topo_mol_create(psf->defs);
  topo_mol_error_handler(psf->mol,interp,newhandle_msg);

  return TCL_OK;
}

int tcl_delatom(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  topo_mol_ident_t target;
  psfgen_data *psf = *(psfgen_data **)data;

  if ( argc < 2 ) {
    Tcl_SetResult(interp,"arguments: segid [ resid? [ aname? ]]", TCL_VOLATILE);
    return TCL_ERROR;
  }

  target.segid = argv[1];
  target.resid = argc > 2 ? argv[2] : 0;
  target.aname = argc > 3 ? argv[3] : 0;

  topo_mol_delete_atom(psf->mol, &target);
 
  return TCL_OK;
}
 
#endif

