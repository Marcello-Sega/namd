
#ifndef TOPO_MOL_OUTPUT_H
#define TOPO_MOL_OUTPUT_H

#include <stdio.h>
#include "topo_mol.h"

int topo_mol_write_pdb(topo_mol *mol, FILE *file,
                                void (*print_msg)(const char *));

int topo_mol_write_psf(topo_mol *mol, FILE *file, int charmmfmt,
                                void (*print_msg)(const char *));

#endif

