#include <stdlib.h>

#include "psf_file.h"
#include "psf_file_extract.h"
#include "topo_mol_struct.h"

/* General note: in a few places I read various arrays in reverse order. 
   That's because I want psf files emitted by psfgen to have all the atoms,
   bonds, etc. in the same order as in the original psf file.  We have to 
   reverse it because adding to a linked list reverses the order.  Actually,
   if the original psf file comes from some other program, then we might 
   change the order of the bonds, angles, etc. but we can at least guarantee
   that if we read a psf file written by psfgen, then write it out again,
   the output will match the input exactly.
*/

/* Read in all psf atom information using this struct */
struct psfatom {
  char name[8];
  char atype[8];
  char resname[8];
  char segname[8];
  char resid[8];
  double charge, mass;
};
typedef struct psfatom psfatom;

  
static int extract_bonds(FILE *file, topo_mol *mol, int natoms, 
                         topo_mol_atom_t **molatomlist) {

  int *bonds;
  int i, nbonds;

  /* Build bonds */
  nbonds = psf_start_bonds(file);
  if (nbonds < 0) {
    printf("Couldn't find start of bonds section");
    return -1; 
  }
  bonds = (int *)malloc(2*nbonds*sizeof(int));

  if (psf_get_bonds(file, nbonds, bonds)) {
    printf("Error reading bonds!\n");
    free(bonds);
    return -1;
  }
 
  for (i=nbonds-1; i >= 0; i--) {
    topo_mol_atom_t *atom1, *atom2;
    topo_mol_bond_t *tuple;
    int ind1, ind2; 
  
    ind1 = bonds[2*i]-1; 
    ind2 = bonds[2*i+1]-1;
    if (ind1 < 0 || ind2 < 0 || ind1 >= natoms || ind2 >= natoms) {
      printf("Bad bond indices: %d %d\n", ind1, ind2);
      /* Bad indices, abort now */
      free(bonds);
      return -1;
    }
    
    atom1 = molatomlist[ind1];
    atom2 = molatomlist[ind2];
   
    tuple = memarena_alloc(mol->arena,sizeof(topo_mol_bond_t));
    tuple->next[0] = atom1->bonds;
    tuple->atom[0] = atom1;
    tuple->next[1] = atom2->bonds;
    tuple->atom[1] = atom2;
    tuple->del = 0;
 
    atom1->bonds = tuple; 
    atom2->bonds = tuple;
  }
  free(bonds);
  return 0;
}

static int extract_angles(FILE *file, topo_mol *mol, int natoms, 
                         topo_mol_atom_t **molatomlist) {

  int i, nangles;
  int *angles;
  
  nangles = psf_start_angles(file);
  if (nangles < 0) return -1; 
  angles = (int *)malloc(3*nangles*sizeof(int));

  if (psf_get_angles(file, nangles, angles)) {
    free(angles); 
    return -1; 
  } 
  
  for (i=nangles-1; i >= 0; i--) {
    topo_mol_atom_t *atom1, *atom2, *atom3;
    topo_mol_angle_t *tuple;

    atom1 = molatomlist[angles[3*i]-1];
    atom2 = molatomlist[angles[3*i+1]-1];
    atom3 = molatomlist[angles[3*i+2]-1];
   
    tuple = memarena_alloc(mol->arena,sizeof(topo_mol_angle_t));
    tuple->next[0] = atom1->angles;
    tuple->atom[0] = atom1;
    tuple->next[1] = atom2->angles;
    tuple->atom[1] = atom2;
    tuple->next[2] = atom3->angles;
    tuple->atom[2] = atom3;
    tuple->del = 0;
 
    atom1->angles = tuple; 
    atom2->angles = tuple;
    atom3->angles = tuple;
  }
  free(angles);
  return 0;
}

static int extract_dihedrals(FILE *file, topo_mol *mol, int natoms,
                         topo_mol_atom_t **molatomlist) {

  int i, ndihedrals;
  int *dihedrals;

  ndihedrals = psf_start_dihedrals(file);
  if (ndihedrals < 0) return -1; 
  dihedrals = (int *)malloc(4*ndihedrals*sizeof(int));

  if (psf_get_dihedrals(file, ndihedrals, dihedrals)) {
    free(dihedrals); 
    return -1;
  }
   
  for (i=ndihedrals-1; i >= 0; i--) {
    topo_mol_atom_t *atom1, *atom2, *atom3, *atom4;
    topo_mol_dihedral_t *tuple;

    atom1 = molatomlist[dihedrals[4*i]-1];
    atom2 = molatomlist[dihedrals[4*i+1]-1];
    atom3 = molatomlist[dihedrals[4*i+2]-1];
    atom4 = molatomlist[dihedrals[4*i+3]-1];

    tuple = memarena_alloc(mol->arena,sizeof(topo_mol_dihedral_t));
    tuple->next[0] = atom1->dihedrals;
    tuple->atom[0] = atom1;
    tuple->next[1] = atom2->dihedrals;
    tuple->atom[1] = atom2;
    tuple->next[2] = atom3->dihedrals;
    tuple->atom[2] = atom3;
    tuple->next[3] = atom4->dihedrals;
    tuple->atom[3] = atom4;
    tuple->del = 0;

    atom1->dihedrals = tuple;
    atom2->dihedrals = tuple;
    atom3->dihedrals = tuple;
    atom4->dihedrals = tuple;
  }
  free(dihedrals);
  return 0;
}

static int extract_impropers(FILE *file, topo_mol *mol, int natoms,
                         topo_mol_atom_t **molatomlist) {

  int i, nimpropers;
  int *impropers;
  
  nimpropers = psf_start_impropers(file);
  if (nimpropers < 0) return -1; 
  impropers = (int *)malloc(4*nimpropers*sizeof(int));

  if (psf_get_impropers(file, nimpropers, impropers)) {
    free(impropers); 
    return -1;
  } 
    
  for (i=nimpropers-1; i >= 0; i--) {
    topo_mol_atom_t *atom1, *atom2, *atom3, *atom4;
    topo_mol_improper_t *tuple;

    atom1 = molatomlist[impropers[4*i]-1];
    atom2 = molatomlist[impropers[4*i+1]-1];
    atom3 = molatomlist[impropers[4*i+2]-1];
    atom4 = molatomlist[impropers[4*i+3]-1];
   
    tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
    tuple->next[0] = atom1->impropers;
    tuple->atom[0] = atom1;
    tuple->next[1] = atom2->impropers;
    tuple->atom[1] = atom2;
    tuple->next[2] = atom3->impropers;
    tuple->atom[2] = atom3;
    tuple->next[3] = atom4->impropers;
    tuple->atom[3] = atom4;
    tuple->del = 0;
 
    atom1->impropers = tuple; 
    atom2->impropers = tuple;
    atom3->impropers = tuple;
    atom4->impropers = tuple;
  }
  free(impropers);
  return 0;
}

/* Return the segment corresponding to the given segname.  If the segname
   doesn't exist, add it.  Return NULL on error.
*/
static topo_mol_segment_t *get_segment(topo_mol *mol, const char *segname) {
  int id;
  topo_mol_segment_t *seg = NULL;
  
  if ( (id = hasharray_index(mol->segment_hash, segname)) != HASHARRAY_FAIL) {
    /* Then the segment exists.  Look it up and return it. */
    seg = mol->segment_array[id];
  } else {
    /* Must create new segment */
    id = hasharray_insert(mol->segment_hash, segname);
    if (id != HASHARRAY_FAIL) {
      seg = mol->segment_array[id] =
            (topo_mol_segment_t *) malloc(sizeof(topo_mol_segment_t)); 
      strcpy(seg->segid, segname);
      seg->residue_hash = hasharray_create(
        (void**) &(seg->residue_array), sizeof(topo_mol_residue_t));
      strcpy(seg->pfirst,"");
      strcpy(seg->plast,"");
      seg->auto_angles = 0; 
      seg->auto_dihedrals = 0; 
    }
  }
  return seg;
}

/* Return a new residue with the given resid.  Add it to the given segment.
   If the resid already exists, return NULL.  Return NULL if there's a problem.
*/

static topo_mol_residue_t *get_residue(topo_mol_segment_t *seg, 
        const char *resid) {
  
  int id;
  topo_mol_residue_t *res;
  
  /* Check that the residue doesn't already exist */
  if ( hasharray_index(seg->residue_hash,resid) != HASHARRAY_FAIL ) {
    return NULL; 
  }
  id = hasharray_insert(seg->residue_hash, resid);
  if (id == HASHARRAY_FAIL) {
    return NULL;
  }
  res = &(seg->residue_array[id]);
  strcpy(res->resid, resid);
  
  return res;
}


int psf_file_extract(topo_mol *mol, FILE *file, 
                                void (*print_msg)(const char *)) {

  int i, natoms;
  psfatom *atomlist;
  topo_mol_atom_t **molatomlist;

  natoms = psf_start_atoms(file);
  if (natoms < 0) return -1;
 
  atomlist = (psfatom *)malloc(natoms * sizeof(psfatom));
  molatomlist = (topo_mol_atom_t **)malloc(natoms * sizeof(topo_mol_atom_t *));

  /* Read in all atoms */
  for (i=0; i<natoms; i++) {
    psfatom *atom = atomlist + i;
    if (psf_get_atom(file, atom->name,atom->atype,atom->resname, atom->segname,
                     atom->resid, &atom->charge, &atom->mass)
        < 0) {
      print_msg("error reading atoms");
      return -1;
    }
  }
 
  i=0; 
  while (i < natoms) {
    topo_mol_segment_t *seg;
    topo_mol_residue_t *res;
    topo_mol_atom_t *atomtmp;
    int firstatom, j;
    const char *resid;
    seg = get_segment(mol, atomlist[i].segname);
    if (!seg) { 
      print_msg("ERROR: unable to get segment!");
      break;
    }
    res = get_residue(seg, atomlist[i].resid);
    if (!res) {
      print_msg("Unable to add residue!");
      break;
    }
    resid = atomlist[i].resid;
    strcpy(res->name, atomlist[i].resname);
    res->atoms = 0;
    firstatom = i;
    while (i<natoms && !strcmp(resid, atomlist[i].resid)) {
      /* Add atoms to residue */
      atomtmp = memarena_alloc(mol->arena, sizeof(topo_mol_atom_t));
      atomtmp->bonds = 0;
      atomtmp->angles = 0;
      atomtmp->dihedrals = 0;
      atomtmp->impropers = 0;
      atomtmp->conformations = 0;
      strcpy(atomtmp->name, atomlist[i].name);
      strcpy(atomtmp->type, atomlist[i].atype);
      atomtmp->mass = atomlist[i].mass; 
      atomtmp->charge = atomlist[i].charge;
      atomtmp->x = 0;       
      atomtmp->y = 0;       
      atomtmp->z = 0;       
      atomtmp->xyz_state = TOPO_MOL_XYZ_VOID;
      atomtmp->typeid = 0;
      atomtmp->atomid = 0;

      /* Save pointer to atom in my table so I can put in the bond 
         information without having find the atom.
      */
      molatomlist[i] = atomtmp;
      i++;
    }
    for (j=i-1; j >= firstatom; j--) {
      /* Add new atoms to head of linked list in reverse order, so that
         the linked list is in the order they appear in the psf file. 
      */
      atomtmp = molatomlist[j];
      atomtmp->next = res->atoms;
      res->atoms = atomtmp;
    }  
  }  
  
  /* Check to see if we broke out of the loop prematurely */
  if (i != natoms) {
    free(atomlist);
    free(molatomlist);
    return -1;
  }
    
  if (extract_bonds(file, mol, natoms, molatomlist)) {
    print_msg("Error processing bonds");
    free(atomlist);
    free(molatomlist);
    return -1;
  }
 
  if (extract_angles(file, mol, natoms, molatomlist)) {
    print_msg("Error processing angles");
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  if (extract_dihedrals(file, mol, natoms, molatomlist)) {
    print_msg("Error processing dihedrals");
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  if (extract_impropers(file, mol, natoms, molatomlist)) {
    print_msg("Error processing impropers");
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  free(atomlist);
  free(molatomlist);
  return 0;
}


