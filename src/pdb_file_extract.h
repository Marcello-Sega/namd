
#ifndef PDB_FILE_EXTRACT_H
#define PDB_FILE_EXTRACT_H

#include <stdio.h>
#include "stringhash.h"
#include "topo_mol.h"

int pdb_file_extract_residues(topo_mol *mol, FILE *file, stringhash *h,
                                void (*print_msg)(const char *));

int pdb_file_extract_coordinates(topo_mol *mol, FILE *file,
                                const char *segid, stringhash *h,
                                void (*print_msg)(const char *));

#endif

