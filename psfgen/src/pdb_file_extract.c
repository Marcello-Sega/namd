
#include <string.h>
#include <ctype.h>
#include "pdb_file_extract.h"
#include "pdb_file.h"
#include "extract_alias.h"

static void strtoupper(char *s) {
  while ( *s ) { *s = toupper(*s); ++s; }
}

int pdb_file_extract_residues(topo_mol *mol, FILE *file, stringhash *h,
                                void *v,void (*print_msg)(void *,const char *)) {

  char record[PDB_RECORD_LENGTH+2];
  int indx;
  float x,y,z,o,b;
  char name[8], resname[8], chain[8];
  char segname[8], resid[8], insertion[8];
  char oldresid[8];
  const char *realres;
  char msg[128];
  int rcount;

  rcount = 0;
  oldresid[0] = '\0';

  do {
    if((indx = read_pdb_record(file, record)) == PDB_ATOM) {
      get_pdb_fields(record, name, resname, chain,
                   segname, resid, insertion, &x, &y, &z, &o, &b);
      if ( strcmp(oldresid,resid) ) {
        strcpy(oldresid,resid);
        ++rcount;
        strtoupper(resname);
        realres = extract_alias_residue_check(h,resname);
        if ( topo_mol_residue(mol,resid,realres) ) {
          sprintf(msg,"ERROR: failed on residue %s from pdb file",resname);
          print_msg(v,msg);
        }
      }
    }
  } while (indx != PDB_END && indx != PDB_EOF);

  sprintf(msg,"extracted %d residues from pdb file",rcount);
  print_msg(v,msg);
  return 0;
}

int pdb_file_extract_coordinates(topo_mol *mol, FILE *file,
                                const char *segid, stringhash *h,
                                void *v,void (*print_msg)(void *,const char *)) {

  char record[PDB_RECORD_LENGTH+2];
  int indx;
  float x,y,z,o,b;
  char name[8], resname[8], chain[8];
  char segname[8], resid[8], insertion[8];
  topo_mol_ident_t target;
  char msg[128];

  target.segid = segid;

  do {
    if((indx = read_pdb_record(file, record)) == PDB_ATOM) {
      get_pdb_fields(record, name, resname, chain,
                   segname, resid, insertion, &x, &y, &z, &o, &b);
      target.resid = resid;
      strtoupper(resname);
      strtoupper(name);
      target.aname = extract_alias_atom_check(h,resname,name);
      /* Use PDB segid if no segid given */
      if (!segid) {
        target.segid = segname;
      }
      if ( topo_mol_set_xyz(mol,&target,x,y,z) ) {
        sprintf(msg,"Warning: failed to set coordinate for atom %s\t %s:%s\t  %s",name,resname,resid,segid ? segid : segname);
        print_msg(v,msg);
      }
    }
  } while (indx != PDB_END && indx != PDB_EOF);

  return 0;

}


