
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "psfgen.h"
#include "charmm_parse_topo_defs.h"
#include "topo_mol_output.h"
#include "pdb_file_extract.h"
#include "psf_file_extract.h"

#ifndef NAMD_TCL

int main(int argc, char **argv) {
  fprintf(stderr,"%s unavailable on this platform (no Tcl)\n",argv[0]);
  exit(-1);
}

#else

#include <tcl.h>
int Tcl_AppInit(Tcl_Interp *interp);

/* This function gets called if/when the Tcl interpreter is deleted. */
static void psfgen_deleteproc(ClientData cd, Tcl_Interp *interp) {
  psfgen_data *data = (psfgen_data *)cd;
  topo_mol_destroy(data->mol);
  topo_defs_destroy(data->defs);
  stringhash_destroy(data->aliases);
  free(data);
}

void handle_msg(const char *msg) {
  printf("%s\n",msg);
}

int tcl_topology(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_segment(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
int tcl_residue(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
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

int Tcl_AppInit(Tcl_Interp *interp) {

  /* Create psfgen data structures; keep in interp so that other libraries
   * can access them.
   */
  psfgen_data *data = (psfgen_data *)malloc(sizeof(psfgen_data));
  data->defs = topo_defs_create();
  topo_defs_error_handler(data->defs,handle_msg);
  
  data->aliases = stringhash_create();

  data->mol = topo_mol_create(data->defs);
  topo_mol_error_handler(data->mol,handle_msg);
 
  Tcl_SetAssocData(interp, (char *)"Psfgen", psfgen_deleteproc, data);

  Tcl_CreateCommand(interp,"topology",tcl_topology,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"readpsf",tcl_readpsf,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"segment",tcl_segment,
	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateCommand(interp,"residue",tcl_residue,
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

  return TCL_OK;
}

int main(int argc, char **argv) {
  Tcl_Main(argc,argv,Tcl_AppInit);
  /* Never gets here; Tcl_Main evaluates 'exit' in the interpreter */
  exit(0);
}

void strtoupper(char *s) {
  while ( *s ) { *s = toupper(*s); ++s; }
}

char* splitcolon(char *s) {
  while ( *s && *s != ':' ) { ++s; }
  if ( *s ) *(s++) = 0; else s = 0;
  return s;
}

int tcl_topology(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  FILE *defs_file;
  char *filename;
  char msg[128];
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc == 1 ) {
    Tcl_SetResult(interp,"no topology file specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 2 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  filename = argv[1];
  if ( ! ( defs_file = fopen(filename,"r") ) ) {
    sprintf(msg,"ERROR: Unable to open topology file %s\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    return TCL_ERROR;
  } else {
    sprintf(msg,"reading topology file %s\n",filename);
    handle_msg(msg);
    charmm_parse_topo_defs(psf->defs,defs_file,handle_msg);
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
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc == 1 ) {
    Tcl_SetResult(interp,"no psf file specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 2 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  filename = argv[1];
  if ( ! ( psf_file = fopen(filename,"r") ) ) {
    sprintf(msg,"ERROR: Unable to open psf file %s\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    return TCL_ERROR;
  } else {
    sprintf(msg,"reading structure from psf file %s\n",filename);
    handle_msg(msg);
    retval = psf_file_extract(psf->mol, psf_file, handle_msg);
    fclose(psf_file);
  }
  if (retval) 
    return TCL_ERROR;
  return TCL_OK;
}

int tcl_segment(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc < 3 ) {
    Tcl_SetResult(interp,"arguments: segname { commmands }",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);

  sprintf(msg,"building segment %s",argv[1]);
  handle_msg(msg);
  if ( topo_mol_segment(psf->mol,argv[1]) ) {
    Tcl_SetResult(interp,"ERROR: failed on segment",TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( Tcl_Eval(interp,argv[2]) != TCL_OK ) return TCL_ERROR;

  handle_msg("generating structure at end of segment");
  if ( topo_mol_end(psf->mol) ) {
    Tcl_SetResult(interp,"ERROR: failed on end of segment",TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_residue(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc < 3 ) {
    Tcl_SetResult(interp,"arguments: resid resname",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);
  strtoupper(argv[2]);

  sprintf(msg,"adding residue %s:%s",argv[2],argv[1]);
  handle_msg(msg);
  if ( topo_mol_residue(psf->mol,argv[1],argv[2]) ) {
    Tcl_SetResult(interp,"ERROR: failed on residue",TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_coord(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  double x,y,z;
  topo_mol_ident_t target;
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc < 5 ) {
    Tcl_SetResult(interp,"arguments: segid resid atomname { x y z }",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 5 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( sscanf(argv[4],"%lf %lf %lf",&x,&y,&z) != 3 ) {
    Tcl_SetResult(interp,"arguments: segid resid atomname { x y z }",TCL_VOLATILE);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);
  strtoupper(argv[2]);
  strtoupper(argv[3]);

  target.segid = argv[1];
  target.resid = argv[2];
  target.aname = argv[3];
  if ( topo_mol_set_xyz(psf->mol,&target,x,y,z) ) {
    Tcl_SetResult(interp,"ERROR: failed on coord",TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}


int tcl_auto(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  int i, angles, dihedrals;
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc < 2 ) {
    Tcl_SetResult(interp,"arguments: ?angles? ?dihedrals? ?none?",TCL_VOLATILE);
    return TCL_ERROR;
  }

  angles = 0;  dihedrals = 0;
  for ( i = 1; i < argc; ++i ) {
    if ( ! strcmp(argv[i],"angles") ) angles = 1;
    else if ( ! strcmp(argv[i],"dihedrals") ) dihedrals = 1;
    else if ( strcmp(argv[i],"none") ) {
      Tcl_SetResult(interp,"arguments: ?angles? ?dihedrals? ?none?",TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  if ( angles ) handle_msg("enabling angle autogeneration");
  else handle_msg("disabling angle autogeneration");
  if ( topo_mol_segment_auto_angles(psf->mol,angles) ) {
    Tcl_SetResult(interp,"ERROR: failed setting angle autogen",TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( dihedrals ) handle_msg("enabling dihedral autogeneration");
  else handle_msg("disabling dihedral autogeneration");
  if ( topo_mol_segment_auto_dihedrals(psf->mol,dihedrals) ) {
    Tcl_SetResult(interp,"ERROR: failed setting dihedral autogen",TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_alias(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc < 2 ) {
    Tcl_SetResult(interp,"arguments: atom | residue ...",TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( ! strcmp(argv[1],"residue") ) {
    if ( argc < 4 ) {
      Tcl_SetResult(interp,"arguments: residue altres realres",TCL_VOLATILE);
      return TCL_ERROR;
    }
    strtoupper(argv[2]);
    strtoupper(argv[3]);
    sprintf(msg,"aliasing residue %s to %s",argv[2],argv[3]);
    handle_msg(msg);
    if ( extract_alias_residue_define(psf->aliases,argv[2],argv[3]) ) {
      Tcl_SetResult(interp,"ERROR: failed on residue alias",TCL_VOLATILE);
      return TCL_ERROR;
    }
  } else if ( ! strcmp(argv[1],"atom") ) {
    if ( argc < 5 ) {
      Tcl_SetResult(interp,"arguments: atom resname altatom realatom",TCL_VOLATILE);
      return TCL_ERROR;
    }
    strtoupper(argv[2]);
    strtoupper(argv[3]);
    strtoupper(argv[4]);
    sprintf(msg,"aliasing residue %s atom %s to %s",argv[2],argv[3],argv[4]);
    handle_msg(msg);
    if ( extract_alias_atom_define(psf->aliases,argv[2],argv[3],argv[4]) ) {
      Tcl_SetResult(interp,"ERROR: failed on atom alias",TCL_VOLATILE);
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
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc == 1 ) {
    Tcl_SetResult(interp,"no pdb file specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 2 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  filename = argv[1];
  if ( ! ( res_file = fopen(filename,"r") ) ) {
    sprintf(msg,"ERROR: Unable to open pdb file %s to read residues\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    return TCL_ERROR;
  } else {
    sprintf(msg,"reading residues from pdb file %s",filename);
    handle_msg(msg);
    if ( pdb_file_extract_residues(psf->mol,res_file,psf->aliases,handle_msg) ) {
      Tcl_SetResult(interp,"ERROR: failed on reading residues from pdb file",TCL_VOLATILE);
      fclose(res_file);
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
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc < 2 ) {
    Tcl_SetResult(interp,"arguments: pdbfile ?segid?",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 3 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }
  filename = argv[1];
  if ( ! ( res_file = fopen(filename,"r") ) ) {
    sprintf(msg,"ERROR: Unable to open pdb file %s to read coordinates\n",filename);
    Tcl_SetResult(interp,msg,TCL_VOLATILE);
    return TCL_ERROR;
  } else {
    const char *segid;
    if (argc == 3) {
      /* Read only coordinates for given segid */
      strtoupper(argv[2]);
      sprintf(msg,"reading coordinates from pdb file %s for segment %s",filename,argv[2]);
      handle_msg(msg);
      segid = argv[2];
    } else {
      /* Read all segid's in pdb file */
      sprintf(msg,"reading coordinates from pdb file %s",filename);
      handle_msg(msg);
      segid = NULL;
    } 
    if ( pdb_file_extract_coordinates(psf->mol,res_file,segid,psf->aliases,handle_msg) ) {
      Tcl_SetResult(interp,"ERROR: failed on reading coordinates from pdb file",TCL_VOLATILE);
      fclose(res_file);
      return TCL_ERROR;
    }
    fclose(res_file);
  }

  return TCL_OK;

}

int tcl_guesscoord(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  psfgen_data *psf = (psfgen_data *)data;
  if ( argc > 1 ) {
    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
    return TCL_ERROR;
  }

  handle_msg("guessing coordinates based on topology file");
  if ( topo_mol_guess_xyz(psf->mol) ) {
    Tcl_SetResult(interp,"ERROR: failed on guessing coordinates",TCL_VOLATILE);
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
  psfgen_data *psf = (psfgen_data *)data;

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
  handle_msg(msg);
  if ( topo_mol_write_psf(psf->mol,res_file,charmmfmt,handle_msg) ) {
    Tcl_SetResult(interp,"ERROR: failed on writing structure to psf file",TCL_VOLATILE);
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
  psfgen_data *psf = (psfgen_data *)data;

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
  handle_msg(msg);
  if ( topo_mol_write_pdb(psf->mol,res_file,handle_msg) ) {
    Tcl_SetResult(interp,"ERROR: failed on writing coordinates to pdb file",TCL_VOLATILE);
    fclose(res_file);
    return TCL_ERROR;
  }
  fclose(res_file);

  return TCL_OK;
}

int tcl_first(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc != 2 ) {
    Tcl_SetResult(interp,"argument: presname",TCL_VOLATILE);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);

  sprintf(msg,"setting patch for first residue to %s",argv[1]);
  handle_msg(msg);
  if ( topo_mol_segment_first(psf->mol,argv[1]) ) {
    Tcl_SetResult(interp,"ERROR: failed to set patch for first residue",TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_last(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  char msg[128];
  psfgen_data *psf = (psfgen_data *)data;

  if ( argc != 2 ) {
    Tcl_SetResult(interp,"argument: presname",TCL_VOLATILE);
    return TCL_ERROR;
  }
  strtoupper(argv[1]);

  sprintf(msg,"setting patch for last residue to %s",argv[1]);
  handle_msg(msg);
  if ( topo_mol_segment_last(psf->mol,argv[1]) ) {
    Tcl_SetResult(interp,"ERROR: failed to set patch for last residue",TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}

int tcl_patch(ClientData data, Tcl_Interp *interp,
					int argc, char *argv[]) {
  int i;
  topo_mol_ident_t targets[10];
  char msg[128];
  psfgen_data *psf = (psfgen_data *)data;

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
  handle_msg(msg);
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
  if ( topo_mol_patch(psf->mol,targets,(argc-2),argv[1],0) ) {
    Tcl_SetResult(interp,"ERROR: failed to apply patch",TCL_VOLATILE);
    return TCL_ERROR;
  }

  return TCL_OK;
}

#endif  /* NAMD_TCL */

