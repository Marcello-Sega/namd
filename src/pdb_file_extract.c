
#include "pdb_file_extract.h"
#include "pdb_file.h"
#include "extract_alias.h"

int pdb_file_extract_residues(topo_mol *mol, FILE *file, stringhash *h,
                                void (*print_msg)(const char *)) {

  char record[PDB_RECORD_LENGTH+2];
  int indx;
  int atomno;
  float x,y,z,o,b;
  char name[8], resname[8], chain[8];
  char segname[8], resid[8], insertion[8];
  char oldresid[8];
  const char *realres;
  char msg[128];
  int rcount;

  rcount = 0;
  sprintf(oldresid,"");

  do {
    if((indx = read_pdb_record(file, record)) == PDB_ATOM) {
      get_pdb_fields(record, name, resname, chain,
                   segname, resid, insertion, &x, &y, &z, &o, &b);
      if ( strcmp(oldresid,resid) ) {
        strcpy(oldresid,resid);
        ++rcount;
        realres = extract_alias_residue_check(h,resname);
        if ( topo_mol_residue(mol,resid,realres) ) {
          sprintf(msg,"ERROR: failed on residue %s from pdb file",resname);
          print_msg(msg);
        }
      }
    }
  } while (indx != PDB_END && indx != PDB_EOF);

  sprintf(msg,"extracted %d residues from pdb file",rcount);
  print_msg(msg);
  return 0;

}

int pdb_file_extract_coordinates(topo_mol *mol, FILE *file,
                                const char *segid, stringhash *h,
                                void (*print_msg)(const char *)) {

  char record[PDB_RECORD_LENGTH+2];
  int indx;
  int atomno;
  float x,y,z,o,b;
  char name[8], resname[8], chain[8];
  char segname[8], resid[8], insertion[8];
  topo_mol_ident_t target;
  const char *realres, *realatom;
  char msg[128];

  target.segid = segid;

  do {
    if((indx = read_pdb_record(file, record)) == PDB_ATOM) {
      get_pdb_fields(record, name, resname, chain,
                   segname, resid, insertion, &x, &y, &z, &o, &b);
      target.resid = resid;
      target.aname = extract_alias_atom_check(h,resname,name);
      /* Use PDB segid if no segid given */
      if (!segid) {
        target.segid = segname;
      }
      if ( topo_mol_set_xyz(mol,&target,x,y,z) ) {
        sprintf(msg,"ERROR: failed to set coordinate for atom %s in residue %s:%s of segment %s",name,resname,resid,segid ? segid : segname);
        print_msg(msg);
      }
    }
  } while (indx != PDB_END && indx != PDB_EOF);

  return 0;

}


