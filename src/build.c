
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "charmm_parse_topo_defs.h"
#include "topo_defs.h"
#include "topo_mol.h"
#include "pdb_file_extract.h"


void handle_msg(const char *msg) {
  printf("%s\n",msg);
}

void strtoupper(char *s) {
  while ( *s ) { *s = toupper(*s); ++s; }
}

void scanupper(char *s) {
  scanf("%s",s);
  while ( *s ) { *s = toupper(*s); ++s; }
}

void scanupper2(char *s1, char *s2) {
  scanf("%s %s",s1,s2);
  while ( *s1 ) { *s1 = toupper(*s1); ++s1; }
  while ( *s2 ) { *s2 = toupper(*s2); ++s2; }
}

void scanupper3(char *s1, char *s2, char *s3) {
  scanf("%s %s %s",s1,s2,s3);
  while ( *s1 ) { *s1 = toupper(*s1); ++s1; }
  while ( *s2 ) { *s2 = toupper(*s2); ++s2; }
  while ( *s3 ) { *s3 = toupper(*s3); ++s3; }
}

char* splitcolon(char *s) {
  while ( *s && *s != ':' ) { ++s; }
  if ( *s ) *(s++) = 0; else s = 0;
  return s;
}


int main(int argc, char **argv) {

  topo_defs *defs;
  topo_mol *mol;
  stringhash *aliases;
  FILE *defs_file;
  FILE *res_file;
  topo_mol_ident_t target;
  topo_mol_ident_t targets[10];
  char s1[100], s2[100], s3[100], s4[100];
  char msg[128];
  double x,y,z;

  defs = topo_defs_create();
  topo_defs_error_handler(defs,handle_msg);

  aliases = stringhash_create();

  mol = topo_mol_create(defs);
  topo_mol_error_handler(mol,handle_msg);

  while ( scanf("%s",s1) != EOF ) {
    if ( ! strcmp(s1,"topology") ) {  /*  topology filename  */
      scanf("%s",s2);
      if ( ! ( defs_file = fopen(s2,"r") ) ) {
        sprintf(msg,"ERROR: Unable to open topology file %s\n",s2);
        handle_msg(msg);
      } else {
        sprintf(msg,"reading topology file %s\n",s2);
        handle_msg(msg);
        charmm_parse_topo_defs(defs,defs_file,handle_msg);
        fclose(defs_file);
      }
    } else if ( ! strcmp(s1,"segment") ) {  /*  segment foobar  */
      scanupper(s2);
      sprintf(msg,"building segment %s",s2);
      handle_msg(msg);
      if ( topo_mol_segment(mol,s2) ) {
        handle_msg("ERROR: failed on segment");
      }
    } else if ( ! strcmp(s1,"coord") ) {  /*  coord segid resid aname x y z  */
      if ( scanf("%s %s %s %lf %lf %lf",s2,s3,s4,&x,&y,&z) == 6 ) {
        strtoupper(s2);
        strtoupper(s3);
        strtoupper(s4);
        handle_msg("recognized coord");
        target.segid = s2;
        target.resid = s3;
        target.aname = s4;
        if ( topo_mol_set_xyz(mol,&target,x,y,z) ) {
          handle_msg("ERROR: failed on coord");
        }
      } else {
        handle_msg("ERROR: bad format for coord");
      }
    } else if ( ! strcmp(s1,"angles") ) {  /*  angles [no]auto  */
      scanf("%s",s2);
      if ( ! strcmp(s2,"auto") ) {
        handle_msg("enabling angle autogeneration");
        if ( topo_mol_segment_auto_angles(mol,1) ) {
          handle_msg("ERROR: failed on angles auto");
        }
      } else if ( ! strcmp(s2,"noauto") ) {
        handle_msg("disabling angle autogeneration");
        if ( topo_mol_segment_auto_angles(mol,0) ) {
          handle_msg("ERROR: failed on angles noauto");
        }
      }
    } else if ( ! strcmp(s1,"dihedrals") ) {  /*  dihedrals [no]auto  */
      scanf("%s",s2);
      if ( ! strcmp(s2,"auto") ) {
        handle_msg("enabling dihedral autogeneration");
        if ( topo_mol_segment_auto_dihedrals(mol,1) ) {
          handle_msg("ERROR: failed on dihedrals auto");
        }
      } else if ( ! strcmp(s2,"noauto") ) {
        handle_msg("disabling dihedral autogeneration");
        if ( topo_mol_segment_auto_dihedrals(mol,0) ) {
          handle_msg("ERROR: failed on dihedrals noauto");
        }
      }
    } else if ( ! strcmp(s1,"pdb") ) {  /*  pdb  */
      scanf("%s",s2);
      sprintf(msg,"reading residues from pdb file %s",s2);
      handle_msg(msg);
      if ( ! ( res_file = fopen(s2,"r") ) ) {
        handle_msg("ERROR: failed to open pdb file to read residues");
      } else {
        if ( pdb_file_extract_residues(mol,res_file,aliases,handle_msg) ) {
          handle_msg("ERROR: failed on reading residues from pdb file");
        }
      }
      fclose(res_file);
    } else if ( ! strcmp(s1,"coordpdb") ) {  /*  coordinate pdb  */
      scanf("%s %s",s2,s3);
      sprintf(msg,"reading coordinates from pdb file %s for segment %s",s2,s3);
      handle_msg(msg);
      if ( ! ( res_file = fopen(s2,"r") ) ) {
        handle_msg("ERROR: failed to open pdb file to read coordinates");
      } else {
        if ( pdb_file_extract_coordinates(mol,res_file,s3,aliases,handle_msg) ) {
          handle_msg("ERROR: failed on reading coordinates from pdb file");
        }
      }
      fclose(res_file);
    } else if ( ! strcmp(s1,"first") ) {  /*  first  */
      scanupper(s2);
      sprintf(msg,"setting patch for first residue to %s",s2);
      handle_msg(msg);
      if ( topo_mol_segment_first(mol,s2) ) {
        handle_msg("ERROR: to set patch for first residue");
      }
    } else if ( ! strcmp(s1,"last") ) {  /*  last  */
      scanupper(s2);
      sprintf(msg,"setting patch for last residue to %s",s2);
      handle_msg(msg);
      if ( topo_mol_segment_last(mol,s2) ) {
        handle_msg("ERROR: to set patch for last residue");
      }
    } else if ( ! strcmp(s1,"end") ) {  /*  end  */
      handle_msg("generating structure at end of segment");
      if ( topo_mol_end(mol) ) {
        handle_msg("ERROR: failed on end of segment");
      }
    } else if ( ! strcmp(s1,"alias") ) {  /*  aliases  */
      scanf("%s",s2);
      if ( ! strcmp(s2,"residue") ) {
        scanupper2(s1,s2);  /*  altres realres  */
        sprintf(msg,"aliasing residue %s to %s",s1,s2);
        handle_msg(msg);
        if ( extract_alias_residue_define(aliases,s1,s2) ) {
          handle_msg("ERROR: failed on residue alias");
        }
      } else if ( ! strcmp(s2,"atom") ) {
        scanupper3(s1,s2,s3);  /*  resname altatom realatom  */
        sprintf(msg,"aliasing residue %s atom %s to %s",s1,s2,s3);
        handle_msg(msg);
        if ( extract_alias_atom_define(aliases,s1,s2,s3) ) {
          handle_msg("ERROR: failed on atom alias");
        }
      }
    } else if ( ! strcmp(s1,"guesscoord") ) {  /*  guess coordinates  */
      handle_msg("guessing coordinates based on topology file");
      if ( topo_mol_guess_xyz(mol) ) {
        handle_msg("ERROR: failed on guessing coordinates");
      }
    } else if ( ! strcmp(s1,"writepdb") ) {  /*  coordinate pdb  */
      scanf("%s",s2);
      sprintf(msg,"writing pdb file %s",s2);
      handle_msg(msg);
      if ( ! ( res_file = fopen(s2,"w") ) ) {
        handle_msg("ERROR: failed to open pdb file to write coordinates");
      } else {
        if ( topo_mol_write_pdb(mol,res_file,s3,handle_msg) ) {
          handle_msg("ERROR: failed on writing coordinates to pdb file");
        }
      }
      fclose(res_file);
    } else if ( ! strcmp(s1,"writepsf") ) {  /*  structure psf  */
      scanf("%s",s2);
      sprintf(msg,"writing psf file %s",s2);
      handle_msg(msg);
      if ( ! ( res_file = fopen(s2,"w") ) ) {
        handle_msg("ERROR: failed to open psf file to write structure");
      } else {
        if ( topo_mol_write_psf(mol,res_file,s3,handle_msg) ) {
          handle_msg("ERROR: failed on writing structure to psf file");
        }
      }
      fclose(res_file);
    } else if ( ! strcmp(s1,"residue") ) {  /*  residue 1 ALA  */
      scanupper2(s2,s3);
      sprintf(msg,"adding residue %s:%s",s3,s2);
      handle_msg(msg);
      if ( topo_mol_residue(mol,s2,s3) ) {
        handle_msg("ERROR: failed on residue");
      }
    } else if ( ! strcmp(s1,"patch1") ) {  /*  patch1 ASPP BPTI:3  */
      scanupper2(s2,s3);
      sprintf(msg,"applying patch %s to residue %s",s2,s3);
      handle_msg(msg);
      targets[0].segid = s3;
      targets[0].resid = splitcolon(s3);
      targets[0].aname = 0;
      if ( targets[0].resid ) {
        if ( topo_mol_patch(mol,targets,1,s2,0) ) {
          handle_msg("ERROR: failed to apply patch");
        }
      } else {
        handle_msg("ERROR: resid missing from target");
      }
    } else if ( ! strcmp(s1,"patch2") ) {  /*  patch2 DISU BPTI:5 BPTI:55  */
      scanupper3(s2,s3,s4);
      sprintf(msg,"applying patch %s to residues %s and %s",s2,s3,s4);
      handle_msg(msg);
      targets[0].segid = s3;
      targets[0].resid = splitcolon(s3);
      targets[0].aname = 0;
      targets[1].segid = s4;
      targets[1].resid = splitcolon(s4);
      targets[1].aname = 0;
      if ( targets[0].resid && targets[1].resid ) {
        if ( topo_mol_patch(mol,targets,2,s2,0) ) {
          handle_msg("ERROR: failed to apply patch");
        }
      } else {
        handle_msg("ERROR: resid missing from target");
      }
    } else {  /*  ???  */
      sprintf(msg,"ERROR: failed to recognize %s",s1);
      handle_msg(msg);
    }
  }

  topo_mol_destroy(mol);
  topo_defs_destroy(defs);
  stringhash_destroy(aliases);

}

