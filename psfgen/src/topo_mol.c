
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "topo_defs_struct.h"
#include "topo_mol_struct.h"


topo_mol * topo_mol_create(topo_defs *defs) {
  topo_mol *mol;
  if ( ! defs ) return 0;
  if ( (mol = (topo_mol*) malloc(sizeof(topo_mol))) ) {
    mol->newerror_handler_data = 0;
    mol->newerror_handler = 0;
    mol->defs = defs;
    mol->segment_hash = hasharray_create(
	(void**) &(mol->segment_array), sizeof(topo_mol_segment_t*));
    mol->buildseg = 0;
    mol->arena = memarena_create();
    if ( ! mol->segment_hash || ! mol->arena ) {
      topo_mol_destroy(mol);
      return 0;
    }
  }
  return mol;
}

void topo_mol_destroy(topo_mol *mol) {
  int i,n;
  topo_mol_segment_t *s;
  
  if ( ! mol ) return;
  n = hasharray_count(mol->segment_hash);
  for ( i=0; i<n; ++i ) {
    s = mol->segment_array[i];
    hasharray_destroy(s->residue_hash);
    free(s);
  }
  hasharray_destroy(mol->segment_hash);
  memarena_destroy(mol->arena);
  free((void*)mol);
}

void topo_mol_error_handler(topo_mol *mol, void *v, void (*print_msg)(void *,const char *)) {
  if ( mol ) {
    mol->newerror_handler = print_msg;
    mol->newerror_handler_data = v;
  }
}

/* internal method */
void topo_mol_log_error(topo_mol *mol, const char *msg) {
  if (mol && msg && mol->newerror_handler)
    mol->newerror_handler(mol->newerror_handler_data, msg);
}

int topo_mol_auto_angles(topo_mol *mol, topo_mol_segment_t *seg);
int topo_mol_auto_dihedrals(topo_mol *mol, topo_mol_segment_t *seg);

topo_mol_segment_t * topo_mol_get_seg(topo_mol *mol,
			const topo_mol_ident_t *target) {
  int iseg;
  char errmsg[64 + 3*NAMEMAXLEN];

  if ( ! mol ) return 0;
  iseg = hasharray_index(mol->segment_hash,target->segid);
  if ( iseg == HASHARRAY_FAIL ) {
    sprintf(errmsg,"no segment %s",target->segid);
    topo_mol_log_error(mol,errmsg);
    return 0;
  }
  return mol->segment_array[iseg];
}

topo_mol_residue_t * topo_mol_get_res(topo_mol *mol,
			const topo_mol_ident_t *target, int irel) {
  int nres, ires;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  char errmsg[64 + 3*NAMEMAXLEN];
  seg = topo_mol_get_seg(mol,target);
  if ( ! seg ) return 0;
  nres = hasharray_count(seg->residue_hash);
  ires = hasharray_index(seg->residue_hash,target->resid);
  if ( ires == HASHARRAY_FAIL ) {
    sprintf(errmsg,"no residue %s of segment %s",
					target->resid,target->segid);
    topo_mol_log_error(mol,errmsg);
    return 0;
  }
  if ( (ires+irel) < 0 || (ires+irel) >= nres ) {
    res = seg->residue_array + ires;
    if ( irel < 0 )
      sprintf(errmsg,"no residue %d before %s:%s of segment %s",
		-1*irel,res->name,res->resid,target->segid);
    if ( irel > 0 )
      sprintf(errmsg,"no residue %d past %s:%s of segment %s",
		irel,res->name,res->resid,target->segid);
    topo_mol_log_error(mol,errmsg);
    return 0;
  }

  return (seg->residue_array + ires + irel);
}

topo_mol_atom_t * topo_mol_get_atom(topo_mol *mol,
			const topo_mol_ident_t *target, int irel) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  char errmsg[64 + 3*NAMEMAXLEN];
  res = topo_mol_get_res(mol,target,irel);
  if ( ! res ) return 0;
  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }
  if ( ! atom ) {
    sprintf(errmsg,"no atom %s in residue %s:%s of segment %s",
		target->aname,res->name,res->resid,target->segid);
    topo_mol_log_error(mol,errmsg);
  }
  return atom;
}

int topo_mol_segment(topo_mol *mol, const char *segid) {
  int i;
  topo_mol_segment_t *newitem;
  char errmsg[32 + NAMEMAXLEN];
  if ( ! mol ) return -1;
  mol->buildseg = 0;
  if ( NAMETOOLONG(segid) ) return -2;
  if ( ( i = hasharray_index(mol->segment_hash,segid) ) != HASHARRAY_FAIL ) {
    sprintf(errmsg,"duplicate segment key %s",segid);
    topo_mol_log_error(mol,errmsg);
    return -3;
  } else {
    i = hasharray_insert(mol->segment_hash,segid);
    if ( i == HASHARRAY_FAIL ) return -4;
    newitem = mol->segment_array[i] =
		(topo_mol_segment_t*) malloc(sizeof(topo_mol_segment_t));
    if ( ! newitem ) return -5;
  }
  strcpy(newitem->segid,segid);
  newitem->residue_hash = hasharray_create(
	(void**) &(newitem->residue_array), sizeof(topo_mol_residue_t));
  strcpy(newitem->pfirst,"");
  strcpy(newitem->plast,"");
  newitem->auto_angles = mol->defs->auto_angles;
  newitem->auto_dihedrals = mol->defs->auto_dihedrals;
  mol->buildseg = newitem;
  return 0;
}

int topo_mol_segment_first(topo_mol *mol, const char *rname) {
  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for first patch");
    return -1;
  }
  if ( NAMETOOLONG(rname) ) return -2;
  strcpy(mol->buildseg->pfirst,rname);
  return 0;
}

int topo_mol_segment_last(topo_mol *mol, const char *rname) {
  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for last patch");
    return -1;
  }
  if ( NAMETOOLONG(rname) ) return -2;
  strcpy(mol->buildseg->plast,rname);
  return 0;
}

int topo_mol_segment_auto_angles(topo_mol *mol, int autogen) {
  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for auto angles");
    return -1;
  }
  mol->buildseg->auto_angles = autogen;
  return 0;
}

int topo_mol_segment_auto_dihedrals(topo_mol *mol, int autogen) {
  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for auto dihedrals");
    return -1;
  }
  mol->buildseg->auto_dihedrals = autogen;
  return 0;
}

int topo_mol_residue(topo_mol *mol, const char *resid, const char *rname) {
  int i;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *newitem;
  char errmsg[32 + NAMEMAXLEN];

  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for residue");
    return -1;
  }
  seg = mol->buildseg;
  if ( NAMETOOLONG(resid) ) return -2;
  if ( NAMETOOLONG(rname) ) return -3;
  if ( hasharray_index(seg->residue_hash,resid) != HASHARRAY_FAIL ) {
    sprintf(errmsg,"duplicate residue key %s",resid);
    topo_mol_log_error(mol,errmsg);
    return -3;
  }

  if ( hasharray_index(mol->defs->residue_hash,rname) == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown residue type %s",rname);
    topo_mol_log_error(mol,errmsg);
  }

  i = hasharray_insert(seg->residue_hash,resid);
  if ( i == HASHARRAY_FAIL ) return -4;
  newitem = &(seg->residue_array[i]);
  strcpy(newitem->resid,resid);
  strcpy(newitem->name,rname);
  newitem->atoms = 0;

  return 0;
}

int topo_mol_mutate(topo_mol *mol, const char *resid, const char *rname) {
  int ires;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  char errmsg[32 + 3*NAMEMAXLEN];

  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for mutate");
    return -1;
  }
  seg = mol->buildseg;

  if ( NAMETOOLONG(resid) ) return -2;
  if ( NAMETOOLONG(rname) ) return -3;
  ires = hasharray_index(seg->residue_hash,resid);
  if ( ires == HASHARRAY_FAIL ) {
    sprintf(errmsg,"residue %s does not exist",resid);
    topo_mol_log_error(mol,errmsg);
    return -1;
  }
  res = seg->residue_array + ires;
  sprintf(errmsg,"mutating residue %s from %s to %s",resid,res->name,rname);
  topo_mol_log_error(mol,errmsg);

  if ( hasharray_index(mol->defs->residue_hash,rname) == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown residue type %s",rname);
    topo_mol_log_error(mol,errmsg);
  }

  strcpy(res->name,rname);

  return 0;
}

topo_mol_atom_t * topo_mol_unlink_atom(
		topo_mol_atom_t **atoms, const char *aname) {
  topo_mol_atom_t **atom;
  topo_mol_atom_t *oldatom;
  if ( ! atoms ) return 0;
  for ( atom = atoms ; *atom; atom = &((*atom)->next) ) {
    if ( ! strcmp(aname,(*atom)->name) ) break;
  }
  oldatom = *atom;
  if ( *atom ) *atom = ((*atom)->next);
  return oldatom;
}

int topo_mol_add_atom(topo_mol *mol, topo_mol_atom_t **atoms,
		topo_mol_atom_t **oldatoms, topo_defs_atom_t *atomdef) {
  int idef;
  topo_mol_atom_t *atomtmp;
  topo_defs_type_t *atype;
  char errmsg[128];
  if ( ! mol || ! atoms ) return -1;
  atomtmp = 0;
  if ( oldatoms ) atomtmp = topo_mol_unlink_atom(oldatoms,atomdef->name);
  if ( ! atomtmp ) {
    atomtmp = memarena_alloc(mol->arena,sizeof(topo_mol_atom_t));
    if ( ! atomtmp ) return -2;
    strcpy(atomtmp->name,atomdef->name);
    atomtmp->bonds = 0;
    atomtmp->angles = 0;
    atomtmp->dihedrals = 0;
    atomtmp->impropers = 0;
    atomtmp->conformations = 0;
    atomtmp->x = 0;
    atomtmp->y = 0;
    atomtmp->z = 0;
    atomtmp->xyz_state = TOPO_MOL_XYZ_VOID;
    atomtmp->partition = 0;
    atomtmp->atomid = 0;
  }
  atomtmp->copy = 0;
  atomtmp->next = *atoms;
  atomtmp->charge = atomdef->charge;
  strcpy(atomtmp->type,atomdef->type);
  atomtmp->typeid = 0;
  atomtmp->mass = 0;
  idef = hasharray_index(mol->defs->type_hash,atomtmp->type);
  if ( idef == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown atom type %s",atomtmp->type);
    topo_mol_log_error(mol,errmsg);
    return -3;
  } else {
    atype = &(mol->defs->type_array[idef]);
    atomtmp->typeid = atype->id;
    atomtmp->mass = atype->mass;
  }
  *atoms = atomtmp;
  return 0;
}

topo_mol_bond_t * topo_mol_bond_next(
		topo_mol_bond_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  return 0;
}

topo_mol_angle_t * topo_mol_angle_next(
		topo_mol_angle_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  return 0;
}

topo_mol_dihedral_t * topo_mol_dihedral_next(
		topo_mol_dihedral_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  if ( tuple->atom[3] == atom ) return tuple->next[3];
  return 0;
}

topo_mol_improper_t * topo_mol_improper_next(
		topo_mol_improper_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  if ( tuple->atom[3] == atom ) return tuple->next[3];
  return 0;
}

topo_mol_conformation_t * topo_mol_conformation_next(
		topo_mol_conformation_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  if ( tuple->atom[3] == atom ) return tuple->next[3];
  return 0;
}

void topo_mol_destroy_atom(topo_mol_atom_t *atom) {
  topo_mol_bond_t *bondtmp;
  topo_mol_angle_t *angletmp;
  topo_mol_dihedral_t *dihetmp;
  topo_mol_improper_t *imprtmp;
  topo_mol_conformation_t *conftmp;
  if ( ! atom ) return;
  for ( bondtmp = atom->bonds; bondtmp;
		bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
    bondtmp->del = 1;
  }
  for ( angletmp = atom->angles; angletmp;
		angletmp = topo_mol_angle_next(angletmp,atom) ) {
    angletmp->del = 1;
  }
  for ( dihetmp = atom->dihedrals; dihetmp;
		dihetmp = topo_mol_dihedral_next(dihetmp,atom) ) {
    dihetmp->del = 1;
  }
  for ( imprtmp = atom->impropers; imprtmp;
		imprtmp = topo_mol_improper_next(imprtmp,atom) ) {
    imprtmp->del = 1;
  }
  for ( conftmp = atom->conformations; conftmp;
		conftmp = topo_mol_conformation_next(conftmp,atom) ) {
    conftmp->del = 1;
  }
}

void topo_mol_del_atom(topo_mol_residue_t *res, const char *aname) {
  if ( ! res ) return;
  topo_mol_destroy_atom(topo_mol_unlink_atom(&(res->atoms),aname));
}


int topo_mol_add_bond(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_bond_t *def) {
  topo_mol_bond_t *tuple;
  topo_mol_atom_t *a1, *a2;
  topo_mol_ident_t t1, t2;
  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_bond_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->bonds;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->bonds;
  tuple->atom[1] = a2;
  tuple->del = 0;
  a1->bonds = tuple;
  a2->bonds = tuple;
  return 0;
}

void topo_mol_del_bond(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_bond_t *def) {
  topo_mol_bond_t *tuple;
  topo_mol_atom_t *a1, *a2;
  topo_mol_ident_t t1, t2;
  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  for ( tuple = a1->bonds; tuple;
		tuple = topo_mol_bond_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2 ) tuple->del = 1;
    if ( tuple->atom[0] == a2 && tuple->atom[1] == a1 ) tuple->del = 1;
  }
}


int topo_mol_add_angle(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_angle_t *def) {
  topo_mol_angle_t *tuple;
  topo_mol_atom_t *a1, *a2, *a3;
  topo_mol_ident_t t1, t2, t3;
  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( ! a3 ) return -7;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_angle_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->angles;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->angles;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->angles;
  tuple->atom[2] = a3;
  tuple->del = 0;
  a1->angles = tuple;
  a2->angles = tuple;
  a3->angles = tuple;
  return 0;
}

void topo_mol_del_angle(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_angle_t *def) {
  topo_mol_angle_t *tuple;
  topo_mol_atom_t *a1, *a2, *a3;
  topo_mol_ident_t t1, t2, t3;
  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  for ( tuple = a1->angles; tuple;
		tuple = topo_mol_angle_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2
	&& tuple->atom[2] == a3 ) tuple->del = 1;
    if ( tuple->atom[0] == a3 && tuple->atom[1] == a2
	&& tuple->atom[2] == a1 ) tuple->del = 1;
  }
}


int topo_mol_add_dihedral(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_dihedral_t *def) {
  topo_mol_dihedral_t *tuple;
  topo_mol_atom_t *a1, *a2, *a3, *a4;
  topo_mol_ident_t t1, t2, t3, t4;
  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( ! a3 ) return -7;
  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  if ( ! a4 ) return -9;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_dihedral_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->dihedrals;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->dihedrals;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->dihedrals;
  tuple->atom[2] = a3;
  tuple->next[3] = a4->dihedrals;
  tuple->atom[3] = a4;
  tuple->del = 0;
  a1->dihedrals = tuple;
  a2->dihedrals = tuple;
  a3->dihedrals = tuple;
  a4->dihedrals = tuple;
  return 0;
}

void topo_mol_del_dihedral(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_dihedral_t *def) {
  topo_mol_dihedral_t *tuple;
  topo_mol_atom_t *a1, *a2, *a3, *a4;
  topo_mol_ident_t t1, t2, t3, t4;
  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  for ( tuple = a1->dihedrals; tuple;
		tuple = topo_mol_dihedral_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2
	&& tuple->atom[2] == a3 && tuple->atom[3] == a4 ) tuple->del = 1;
    if ( tuple->atom[0] == a4 && tuple->atom[1] == a3
	&& tuple->atom[2] == a2 && tuple->atom[3] == a1 ) tuple->del = 1;
  }
}


int topo_mol_add_improper(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_improper_t *def) {
  topo_mol_improper_t *tuple;
  topo_mol_atom_t *a1, *a2, *a3, *a4;
  topo_mol_ident_t t1, t2, t3, t4;
  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( ! a3 ) return -7;
  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  if ( ! a4 ) return -9;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->impropers;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->impropers;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->impropers;
  tuple->atom[2] = a3;
  tuple->next[3] = a4->impropers;
  tuple->atom[3] = a4;
  tuple->del = 0;
  a1->impropers = tuple;
  a2->impropers = tuple;
  a3->impropers = tuple;
  a4->impropers = tuple;
  return 0;
}

void topo_mol_del_improper(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_improper_t *def) {
  topo_mol_improper_t *tuple;
  topo_mol_atom_t *a1, *a2, *a3, *a4;
  topo_mol_ident_t t1, t2, t3, t4;
  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  for ( tuple = a1->impropers; tuple;
		tuple = topo_mol_improper_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2
	&& tuple->atom[2] == a3 && tuple->atom[3] == a4 ) tuple->del = 1;
    if ( tuple->atom[0] == a4 && tuple->atom[1] == a3
	&& tuple->atom[2] == a2 && tuple->atom[3] == a1 ) tuple->del = 1;
  }
}


int topo_mol_add_conformation(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_conformation_t *def) {
  topo_mol_conformation_t *tuple;
  topo_mol_atom_t *a1, *a2, *a3, *a4;
  topo_mol_ident_t t1, t2, t3, t4;
  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( ! a3 ) return -7;
  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  if ( ! a4 ) return -9;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_conformation_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->conformations;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->conformations;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->conformations;
  tuple->atom[2] = a3;
  tuple->next[3] = a4->conformations;
  tuple->atom[3] = a4;
  tuple->del = 0;
  tuple->improper = def->improper;
  tuple->dist12 = def->dist12;
  tuple->angle123 = def->angle123;
  tuple->dihedral = def->dihedral;
  tuple->angle234 = def->angle234;
  tuple->dist34 = def->dist34;
  a1->conformations = tuple;
  a2->conformations = tuple;
  a3->conformations = tuple;
  a4->conformations = tuple;
  return 0;
}

void topo_mol_del_conformation(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_conformation_t *def) {
  topo_mol_conformation_t *tuple;
  topo_mol_atom_t *a1, *a2, *a3, *a4;
  topo_mol_ident_t t1, t2, t3, t4;
  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  for ( tuple = a1->conformations; tuple;
		tuple = topo_mol_conformation_next(tuple,a1) ) {
    if ( tuple->improper == def->improper
	&&  tuple->atom[0] == a1 && tuple->atom[1] == a2
	&& tuple->atom[2] == a3 && tuple->atom[3] == a4 ) tuple->del = 1;
    if ( tuple->improper == def->improper
	&& tuple->atom[0] == a4 && tuple->atom[1] == a3
	&& tuple->atom[2] == a2 && tuple->atom[3] == a1 ) tuple->del = 1;
  }
}


int topo_mol_end(topo_mol *mol) {
  int i,n;
  int idef;
  topo_defs *defs;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_defs_residue_t *resdef;
  topo_defs_atom_t *atomdef;
  topo_defs_bond_t *bonddef;
  topo_defs_angle_t *angldef;
  topo_defs_dihedral_t *dihedef;
  topo_defs_improper_t *imprdef;
  topo_defs_conformation_t *confdef;
  topo_mol_ident_t target;
  char errmsg[128];

  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for end");
    return -1;
  }
  seg = mol->buildseg;
  mol->buildseg = 0;
  defs = mol->defs;

  /* add atoms */
  n = hasharray_count(seg->residue_hash);
  for ( i=0; i<n; ++i ) {
    res = &(seg->residue_array[i]);
    idef = hasharray_index(defs->residue_hash,res->name);
    if ( idef == HASHARRAY_FAIL ) {
      sprintf(errmsg,"unknown residue type %s",res->name);
      topo_mol_log_error(mol,errmsg);
      return -1;
    }
    resdef = &(mol->defs->residue_array[idef]);
    if ( resdef->patch ) {
      sprintf(errmsg,"unknown residue type %s",res->name);
      topo_mol_log_error(mol,errmsg);
      return -1;
    }

    /* patches */
    if ( i==0 && ! strlen(seg->pfirst) ) strcpy(seg->pfirst,resdef->pfirst);
    if ( i==(n-1) && ! strlen(seg->plast) ) strcpy(seg->plast,resdef->plast);

    for ( atomdef = resdef->atoms; atomdef; atomdef = atomdef->next ) {
      if ( topo_mol_add_atom(mol,&(res->atoms),0,atomdef) ) { 
        sprintf(errmsg,"add atom failed in residue %s:%s",res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
        return -8;
      }
    }
  }

  for ( i=0; i<n; ++i ) {
    res = &(seg->residue_array[i]);
    idef = hasharray_index(defs->residue_hash,res->name);
    if ( idef == HASHARRAY_FAIL ) {
      sprintf(errmsg,"unknown residue type %s",res->name);
      topo_mol_log_error(mol,errmsg);
      return -1;
    }
    resdef = &(mol->defs->residue_array[idef]);
    target.segid = seg->segid;
    target.resid = res->resid;
    for ( bonddef = resdef->bonds; bonddef; bonddef = bonddef->next ) {
      if ( topo_mol_add_bond(mol,&target,1,bonddef) ) {
        sprintf(errmsg,"Warning: add bond %s(%d) %s(%d) failed in residue %s:%s",
		bonddef->atom1,bonddef->rel1,bonddef->atom2,bonddef->rel2,
		res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
      }
    }
    if ( seg->auto_angles && resdef->angles ) {
      sprintf(errmsg,"Warning: explicit angles in residue %s:%s may conflict with autogeneration",res->name,res->resid);
      topo_mol_log_error(mol,errmsg);
    }
    for ( angldef = resdef->angles; angldef; angldef = angldef->next ) {
      if ( topo_mol_add_angle(mol,&target,1,angldef) ) {
        sprintf(errmsg,"Warning: add angle failed in residue %s:%s",res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
      }
    }
    if ( seg->auto_dihedrals && resdef->dihedrals) {
      sprintf(errmsg,"Warning: explicit dihedrals in residue %s:%s may conflict with autogeneration",res->name,res->resid);
      topo_mol_log_error(mol,errmsg);
    }
    for ( dihedef = resdef->dihedrals; dihedef; dihedef = dihedef->next ) {
      if ( topo_mol_add_dihedral(mol,&target,1,dihedef) ) {
        sprintf(errmsg,"Warning: add dihedral failed in residue %s:%s",res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
      }
    }
    for ( imprdef = resdef->impropers; imprdef; imprdef = imprdef->next ) {
      if ( topo_mol_add_improper(mol,&target,1,imprdef) ) {
        sprintf(errmsg,"Warning: add improper failed in residue %s:%s",res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
      }
    }
    for ( confdef = resdef->conformations; confdef; confdef = confdef->next ) {
      if ( topo_mol_add_conformation(mol,&target,1,confdef) ) {
        sprintf(errmsg,"Warning: add conformation %s-%s-%s-%s failed in residue %s:%s",
		confdef->atom1,confdef->atom2,confdef->atom3,confdef->atom4,
		res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
      }
    }
  }

  /* apply patches */

  res = &(seg->residue_array[0]);
  if ( ! strlen(seg->pfirst) ) strcpy(seg->pfirst,"NONE");

  target.segid = seg->segid;
  target.resid = res->resid;
  if ( topo_mol_patch(mol, &target, 1, seg->pfirst, 1,
	seg->auto_angles, seg->auto_dihedrals) ) return -10;

  res = &(seg->residue_array[n-1]);
  if ( ! strlen(seg->plast) ) strcpy(seg->plast,"NONE");

  target.segid = seg->segid;
  target.resid = res->resid;
  if ( topo_mol_patch(mol, &target, 1, seg->plast, 0,
	seg->auto_angles, seg->auto_dihedrals) ) return -11;

  if (topo_mol_auto_angles(mol, seg)) return -12;
  if (topo_mol_auto_dihedrals(mol, seg)) return -13;

  return 0;
}


int topo_mol_auto_angles(topo_mol *mol, topo_mol_segment_t *seg) {
  int ires, nres;
  topo_mol_residue_t *res;
  topo_mol_bond_t *b1, *b2;
  topo_mol_angle_t *tuple;
  topo_mol_atom_t *atom, *a1, *a2, *a3;

  if (! mol) return -1;
  if (! seg) return -2;
  if ( ! seg->auto_angles ) return 0;

  nres = hasharray_count(seg->residue_hash);
  for ( ires=0; ires<nres; ++ires ) {
    res = &(seg->residue_array[ires]);
    for ( atom = res->atoms; atom; atom = atom->next ) {
      a2 = atom;
      for ( b1 = atom->bonds; b1; b1 = topo_mol_bond_next(b1,atom) ) {
        if ( b1->del ) continue;
        if ( b1->atom[0] == atom ) a1 = b1->atom[1];
        else if ( b1->atom[1] == atom ) a1 = b1->atom[0];
        else return -5;
        b2 = b1;  while ( (b2 = topo_mol_bond_next(b2,atom)) ) {
          if ( b2->del ) continue;
          if ( b2->atom[0] == atom ) a3 = b2->atom[1];
          else if ( b2->atom[1] == atom ) a3 = b2->atom[0];
          else return -6;
          tuple = memarena_alloc(mol->arena,sizeof(topo_mol_angle_t));
          if ( ! tuple ) return -10;
          tuple->next[0] = a1->angles;
          tuple->atom[0] = a1;
          tuple->next[1] = a2->angles;
          tuple->atom[1] = a2;
          tuple->next[2] = a3->angles;
          tuple->atom[2] = a3;
          tuple->del = 0;
          a1->angles = tuple;
          a2->angles = tuple;
          a3->angles = tuple;
        }
      }
    }
  }

  return 0;
}


int topo_mol_auto_dihedrals(topo_mol *mol, topo_mol_segment_t *seg) {
  int ires, nres, found, atomid;
  topo_mol_residue_t *res;
  topo_mol_angle_t *g1, *g2;
  topo_mol_dihedral_t *tuple;
  topo_mol_atom_t *atom, *a1, *a2, *a3, *a4;

  if (! mol) return -1;
  if (! seg) return -2;
  if ( ! seg->auto_dihedrals ) return 0;

  /*  number atoms, needed to avoid duplicate dihedrals below  */
  atomid = 0;
  a1 = a2 = a3 = a4 = 0; /* Make compiler happy */
  nres = hasharray_count(seg->residue_hash);
  for ( ires=0; ires<nres; ++ires ) {
    res = &(seg->residue_array[ires]);
    for ( atom = res->atoms; atom; atom = atom->next ) {
      atom->atomid = ++atomid;
    }
  }

  nres = hasharray_count(seg->residue_hash);
  for ( ires=0; ires<nres; ++ires ) {
    res = &(seg->residue_array[ires]);
    for ( atom = res->atoms; atom; atom = atom->next ) {
      for ( g1 = atom->angles; g1; g1 = topo_mol_angle_next(g1,atom) ) {
        if ( g1->del ) continue;
        if ( g1->atom[1] != atom ) continue;
        for ( g2 = atom->angles; g2; g2 = topo_mol_angle_next(g2,atom) ) {
          if ( g2->del ) continue;
          if ( g2->atom[1] == atom ) continue;
          found = 0;
          if ( g2->atom[0] == atom ) {  /*  XBX BXX  */
            if ( g2->atom[1] == g1->atom[0] ) {  /*  CBA BCD  */
              a1 = g1->atom[2];
              a2 = g1->atom[1];  /* == g2->atom[0] */
              a3 = g1->atom[0];  /* == g2->atom[1] */
              a4 = g2->atom[2];
              found = ( a1->atomid < a4->atomid );
            } else if ( g2->atom[1] == g1->atom[2] ) {  /*  ABC BCD  */
              a1 = g1->atom[0];
              a2 = g1->atom[1];  /* == g2->atom[0] */
              a3 = g1->atom[2];  /* == g2->atom[1] */
              a4 = g2->atom[2];
              found = ( a1->atomid < a4->atomid );
            }
          } else if ( g2->atom[2] == atom ) {  /*  XBX XXB  */
            if ( g2->atom[1] == g1->atom[0] ) {  /*  CBA DCB  */
              a1 = g1->atom[2];
              a2 = g1->atom[1];  /* == g2->atom[2] */
              a3 = g1->atom[0];  /* == g2->atom[1] */
              a4 = g2->atom[0];
              found = ( a1->atomid < a4->atomid );
            } else if ( g2->atom[1] == g1->atom[2] ) {  /*  ABC DCB  */
              a1 = g1->atom[0];
              a2 = g1->atom[1];  /* == g2->atom[2] */
              a3 = g1->atom[2];  /* == g2->atom[1] */
              a4 = g2->atom[0];
              found = ( a1->atomid < a4->atomid );
            }
          } else return -6;
          if ( ! found ) continue;
          tuple = memarena_alloc(mol->arena,sizeof(topo_mol_dihedral_t));
          if ( ! tuple ) return -10;
          tuple->next[0] = a1->dihedrals;
          tuple->atom[0] = a1;
          tuple->next[1] = a2->dihedrals;
          tuple->atom[1] = a2;
          tuple->next[2] = a3->dihedrals;
          tuple->atom[2] = a3;
          tuple->next[3] = a4->dihedrals;
          tuple->atom[3] = a4;
          tuple->del = 0;
          a1->dihedrals = tuple;
          a2->dihedrals = tuple;
          a3->dihedrals = tuple;
          a4->dihedrals = tuple;
        }
      }
    }
  }

  return 0;
}

int topo_mol_patch(topo_mol *mol, const topo_mol_ident_t *targets,
                        int ntargets, const char *rname, int prepend,
			int warn_angles, int warn_dihedrals) {

  int idef;
  topo_defs_residue_t *resdef;
  topo_defs_atom_t *atomdef;
  topo_defs_bond_t *bonddef;
  topo_defs_angle_t *angldef;
  topo_defs_dihedral_t *dihedef;
  topo_defs_improper_t *imprdef;
  topo_defs_conformation_t *confdef;
  topo_mol_residue_t *res;
  char errmsg[128];

  if ( ! mol ) return -1;
  if ( mol->buildseg ) return -2;
  if ( ! mol->defs ) return -3;

  idef = hasharray_index(mol->defs->residue_hash,rname);
  if ( idef == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown patch type %s",rname);
    topo_mol_log_error(mol,errmsg);
    return -4;
  }
  resdef = &(mol->defs->residue_array[idef]);
  if ( ! resdef->patch ) {
    sprintf(errmsg,"unknown patch type %s",rname);
    topo_mol_log_error(mol,errmsg);
    return -5;
  }

  for ( atomdef = resdef->atoms; atomdef; atomdef = atomdef->next ) {
    if ( atomdef->res < 0 || atomdef->res >= ntargets ) return -6;
    res = topo_mol_get_res(mol,&targets[atomdef->res],atomdef->rel);
    if ( ! res ) return -7;
    if ( atomdef->del ) topo_mol_del_atom(res,atomdef->name);
    else if ( topo_mol_add_atom(mol,&(res->atoms),&(res->atoms),atomdef) ) {
      sprintf(errmsg,"add atom failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
      return -8;
    }
  }

  for ( bonddef = resdef->bonds; bonddef; bonddef = bonddef->next ) {
    if ( bonddef->del ) topo_mol_del_bond(mol,targets,ntargets,bonddef);
    else if ( topo_mol_add_bond(mol,targets,ntargets,bonddef) ) {
      sprintf(errmsg,"Warning: add bond failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }
  if ( warn_angles && resdef->angles ) {
    sprintf(errmsg,"Warning: explicit angles in patch %s may conflict with autogeneration",rname);
    topo_mol_log_error(mol,errmsg);
  }
  for ( angldef = resdef->angles; angldef; angldef = angldef->next ) {
    if ( angldef->del ) topo_mol_del_angle(mol,targets,ntargets,angldef);
    else if ( topo_mol_add_angle(mol,targets,ntargets,angldef) ) {
      sprintf(errmsg,"Warning: add angle failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }
  if ( warn_dihedrals && resdef->dihedrals ) {
    sprintf(errmsg,"Warning: explicit dihedrals in patch %s may conflict with autogeneration",rname);
    topo_mol_log_error(mol,errmsg);
  }
  for ( dihedef = resdef->dihedrals; dihedef; dihedef = dihedef->next ) {
    if ( dihedef->del ) topo_mol_del_dihedral(mol,targets,ntargets,dihedef);
    else if ( topo_mol_add_dihedral(mol,targets,ntargets,dihedef) ) {
      sprintf(errmsg,"Warning: add dihedral failed in patch %s",rname);
        topo_mol_log_error(mol,errmsg);
      }
  }
  for ( imprdef = resdef->impropers; imprdef; imprdef = imprdef->next ) {
    if ( imprdef->del ) topo_mol_del_improper(mol,targets,ntargets,imprdef);
    else if ( topo_mol_add_improper(mol,targets,ntargets,imprdef) ) {
      sprintf(errmsg,"Warning: add improper failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }
  for ( confdef = resdef->conformations; confdef; confdef = confdef->next ) {
    if ( confdef->del ) topo_mol_del_conformation(mol,targets,ntargets,confdef);
    else if ( topo_mol_add_conformation(mol,targets,ntargets,confdef) ) {
      sprintf(errmsg,"Warning: add conformation failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }
  return 0;
}

int topo_mol_multiply_atoms(topo_mol *mol, const topo_mol_ident_t *targets,
						int ntargets, int ncopies) {
  int ipass, natoms, iatom, icopy;
  const topo_mol_ident_t *target;
  int itarget;
  topo_mol_atom_t *atom, **atoms;
  topo_mol_residue_t *res;
  topo_mol_segment_t *seg;
  int nres, ires;

  if (!mol) return -1;

  /* two passes needed to find atoms */
  for (ipass=0; ipass<2; ++ipass) {
    if ( ipass ) atoms = memarena_alloc(mol->arena,
				natoms*sizeof(topo_mol_atom_t*));
    natoms = 0;
    /* walk all targets */
    for (itarget=0; itarget<ntargets; ++itarget) {
      target = targets + itarget;
      
      if (!target->resid) { /* whole segment */
        seg = topo_mol_get_seg(mol,target);
        if ( ! seg ) return -2;
        nres = hasharray_count(seg->residue_hash);
        for ( ires=0; ires<nres; ++ires ) {
          res = &(seg->residue_array[ires]);
          for ( atom = res->atoms; atom; atom = atom->next ) {
            if ( ipass ) atoms[natoms] = atom;
            ++natoms;
          }
        }
        continue;
      }

      if (!target->aname) { /* whole residue */
        res = topo_mol_get_res(mol,target,0);
        if ( ! res ) return -3;
        for ( atom = res->atoms; atom; atom = atom->next ) {
          if ( ipass ) atoms[natoms] = atom;
          ++natoms;
        }
        continue;
      }

      /* one atom */
      atom = topo_mol_get_atom(mol,target,0);
      if ( ! atom ) return -4;
      if ( ipass ) atoms[natoms] = atom;
      ++natoms;
    }
  }
  
  /* make one copy on each pass through loop */
  for (icopy=1; icopy<ncopies; ++icopy) {

  /* copy the actual atoms */
  for (iatom=0; iatom<natoms; ++iatom) {
    topo_mol_atom_t *newatom;
    atom = atoms[iatom];
    if ( atom->copy ) {
      topo_mol_log_error(mol,"an atom occurs twice in the selection");
      return -20;
    }
    newatom = memarena_alloc(mol->arena,sizeof(topo_mol_atom_t));
    if ( ! newatom ) return -5;
    memcpy(newatom,atom,sizeof(topo_mol_atom_t));
    atom->next = newatom;
    atom->copy = newatom;
    newatom->bonds = 0;
    newatom->angles = 0;
    newatom->dihedrals = 0;
    newatom->impropers = 0;
    newatom->conformations = 0;
  }

  /* copy associated bonds, etc. */
  for (iatom=0; iatom<natoms; ++iatom) {
    topo_mol_atom_t *a1, *a2, *a3, *a4;
    topo_mol_bond_t *bondtmp;
    topo_mol_angle_t *angletmp;
    topo_mol_dihedral_t *dihetmp;
    topo_mol_improper_t *imprtmp;
    topo_mol_conformation_t *conftmp;
    atom = atoms[iatom];
    for ( bondtmp = atom->bonds; bondtmp;
		bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
      topo_mol_bond_t *tuple;
      if ( bondtmp->del ) continue;
      if ( bondtmp->atom[0] == atom || ( ! bondtmp->atom[0]->copy ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_bond_t));
      if ( ! tuple ) return -6;
      a1 = bondtmp->atom[0]->copy; if ( ! a1 ) a1 = bondtmp->atom[0];
      a2 = bondtmp->atom[1]->copy; if ( ! a2 ) a2 = bondtmp->atom[1];
      tuple->next[0] = a1->bonds;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->bonds;
      tuple->atom[1] = a2;
      tuple->del = 0;
      a1->bonds = tuple;
      a2->bonds = tuple;
    }
    for ( angletmp = atom->angles; angletmp;
		angletmp = topo_mol_angle_next(angletmp,atom) ) {
      topo_mol_angle_t *tuple;
      if ( angletmp->del ) continue;
      if ( angletmp->atom[0] == atom || ( ! angletmp->atom[0]->copy
      && ( angletmp->atom[1] == atom || ( ! angletmp->atom[1]->copy ) ) ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_angle_t));
      if ( ! tuple ) return -7;
      a1 = angletmp->atom[0]->copy; if ( ! a1 ) a1 = angletmp->atom[0];
      a2 = angletmp->atom[1]->copy; if ( ! a2 ) a2 = angletmp->atom[1];
      a3 = angletmp->atom[2]->copy; if ( ! a3 ) a3 = angletmp->atom[2];
      tuple->next[0] = a1->angles;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->angles;
      tuple->atom[1] = a2;
      tuple->next[2] = a3->angles;
      tuple->atom[2] = a3;
      tuple->del = 0;
      a1->angles = tuple;
      a2->angles = tuple;
      a3->angles = tuple;
    }
    for ( dihetmp = atom->dihedrals; dihetmp;
		dihetmp = topo_mol_dihedral_next(dihetmp,atom) ) {
      topo_mol_dihedral_t *tuple;
      if ( dihetmp->del ) continue;
      if ( dihetmp->atom[0] == atom || ( ! dihetmp->atom[0]->copy
      && ( dihetmp->atom[1] == atom || ( ! dihetmp->atom[1]->copy
      && ( dihetmp->atom[2] == atom || ( ! dihetmp->atom[2]->copy ) ) ) ) ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_dihedral_t));
      if ( ! tuple ) return -8;
      a1 = dihetmp->atom[0]->copy; if ( ! a1 ) a1 = dihetmp->atom[0];
      a2 = dihetmp->atom[1]->copy; if ( ! a2 ) a2 = dihetmp->atom[1];
      a3 = dihetmp->atom[2]->copy; if ( ! a3 ) a3 = dihetmp->atom[2];
      a4 = dihetmp->atom[3]->copy; if ( ! a4 ) a4 = dihetmp->atom[3];
      tuple->next[0] = a1->dihedrals;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->dihedrals;
      tuple->atom[1] = a2;
      tuple->next[2] = a3->dihedrals;
      tuple->atom[2] = a3;
      tuple->next[3] = a4->dihedrals;
      tuple->atom[3] = a4;
      tuple->del = 0;
      a1->dihedrals = tuple;
      a2->dihedrals = tuple;
      a3->dihedrals = tuple;
      a4->dihedrals = tuple;
    }
    for ( imprtmp = atom->impropers; imprtmp;
		imprtmp = topo_mol_improper_next(imprtmp,atom) ) {
      topo_mol_improper_t *tuple;
      if ( imprtmp->del ) continue;
      if ( imprtmp->atom[0] == atom || ( ! imprtmp->atom[0]->copy
      && ( imprtmp->atom[1] == atom || ( ! imprtmp->atom[1]->copy
      && ( imprtmp->atom[2] == atom || ( ! imprtmp->atom[2]->copy ) ) ) ) ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
      if ( ! tuple ) return -9;
      a1 = imprtmp->atom[0]->copy; if ( ! a1 ) a1 = imprtmp->atom[0];
      a2 = imprtmp->atom[1]->copy; if ( ! a2 ) a2 = imprtmp->atom[1];
      a3 = imprtmp->atom[2]->copy; if ( ! a3 ) a3 = imprtmp->atom[2];
      a4 = imprtmp->atom[3]->copy; if ( ! a4 ) a4 = imprtmp->atom[3];
      tuple->next[0] = a1->impropers;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->impropers;
      tuple->atom[1] = a2;
      tuple->next[2] = a3->impropers;
      tuple->atom[2] = a3;
      tuple->next[3] = a4->impropers;
      tuple->atom[3] = a4;
      tuple->del = 0;
      a1->impropers = tuple;
      a2->impropers = tuple;
      a3->impropers = tuple;
      a4->impropers = tuple;
    }
    for ( conftmp = atom->conformations; conftmp;
		conftmp = topo_mol_conformation_next(conftmp,atom) ) {
      topo_mol_conformation_t *tuple;
      if ( conftmp->del ) continue;
      if ( conftmp->atom[0] == atom || ( ! conftmp->atom[0]->copy
      && ( conftmp->atom[1] == atom || ( ! conftmp->atom[1]->copy
      && ( conftmp->atom[2] == atom || ( ! conftmp->atom[2]->copy ) ) ) ) ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_conformation_t));
      if ( ! tuple ) return -10;
      a1 = conftmp->atom[0]->copy; if ( ! a1 ) a1 = conftmp->atom[0];
      a2 = conftmp->atom[1]->copy; if ( ! a2 ) a2 = conftmp->atom[1];
      a3 = conftmp->atom[2]->copy; if ( ! a3 ) a3 = conftmp->atom[2];
      a4 = conftmp->atom[3]->copy; if ( ! a4 ) a4 = conftmp->atom[3];
      tuple->next[0] = a1->conformations;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->conformations;
      tuple->atom[1] = a2;
      tuple->next[2] = a3->conformations;
      tuple->atom[2] = a3;
      tuple->next[3] = a4->conformations;
      tuple->atom[3] = a4;
      tuple->del = 0;
      tuple->improper = conftmp->improper;
      tuple->dist12 = conftmp->dist12;
      tuple->angle123 = conftmp->angle123;
      tuple->dihedral = conftmp->dihedral;
      tuple->angle234 = conftmp->angle234;
      tuple->dist34 = conftmp->dist34;
      a1->conformations = tuple;
      a2->conformations = tuple;
      a3->conformations = tuple;
      a4->conformations = tuple;
    }
  }

  /* clean up copy pointers */
  for (iatom=0; iatom<natoms; ++iatom) {
    atom = atoms[iatom];
    if ( atom->partition == 0 ) atom->partition = 1;
    atom->copy->partition = atom->partition + 1;
    atoms[iatom] = atom->copy;
    atom->copy = 0;
  }

  } /* icopy */

  return 0;  /* success */
}

/* API function */
void topo_mol_delete_atom(topo_mol *mol, const topo_mol_ident_t *target) {
  
  topo_mol_residue_t *res;
  topo_mol_segment_t *seg;
  int ires, iseg;
  if (!mol) return;

  iseg = hasharray_index(mol->segment_hash,target->segid);
  if ( iseg == HASHARRAY_FAIL ) {
    char errmsg[50];
    sprintf(errmsg,"no segment %s",target->segid);
    topo_mol_log_error(mol,errmsg);
    return;
  }
  seg = mol->segment_array[iseg];
  
  if (!target->resid) {
    /* Delete this segment */
    if (hasharray_delete(mol->segment_hash, target->segid) < 0) {
      topo_mol_log_error(mol, "Unable to delete segment");
    }
    return;
  }

  ires = hasharray_index(seg->residue_hash,target->resid);
  if ( ires == HASHARRAY_FAIL ) {
    char errmsg[50];
    sprintf(errmsg,"no residue %s of segment %s",
                                        target->resid,target->segid);
    topo_mol_log_error(mol,errmsg);
    return;
  }
  res = seg->residue_array+ires;  
  
  if (!target->aname) {  
    /* Must destroy all atoms in residue, since there may be bonds between
       this residue and other atoms 
    */
    topo_mol_atom_t *atom = res->atoms;
    while (atom) {
      topo_mol_destroy_atom(atom);
      atom = atom->next;
    }
    res->atoms = 0;
    hasharray_delete(seg->residue_hash, target->resid); 
    return;
  }
  /* Just delete one atom */
  topo_mol_destroy_atom(topo_mol_unlink_atom(&(res->atoms),target->aname));
}

int topo_mol_set_xyz(topo_mol *mol, const topo_mol_ident_t *target,
                                        double x, double y, double z) {
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;

  atom = topo_mol_get_atom(mol,target,0);
  if ( ! atom ) return -3;

  atom->x = x;
  atom->y = y;
  atom->z = z;
  atom->xyz_state = TOPO_MOL_XYZ_SET;

  return 0;
}

int topo_mol_clear_xyz(topo_mol *mol, const topo_mol_ident_t *target) {
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;

  atom = topo_mol_get_atom(mol,target,0);
  if ( ! atom ) return -3;

  atom->x = 0;
  atom->y = 0;
  atom->z = 0;
  atom->xyz_state = TOPO_MOL_XYZ_VOID;

  return 0;
}


int topo_mol_guess_xyz(topo_mol *mol) {
  char msg[128];
  int iseg,nseg,ires,nres,ucount,i,nk,nu,gcount,gwild,okwild,wcount,hcount;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom, *a1, *a2, *a3;
  topo_mol_atom_t *ka[4];
  topo_mol_atom_t *ua[4];
  topo_mol_bond_t *bondtmp;
  double dihedral, angle234, dist34;
  topo_mol_atom_t **uatoms;
  topo_mol_conformation_t *conf;
  double r12x,r12y,r12z,r23x,r23y,r23z,ix,iy,iz,jx,jy,jz,kx,ky,kz;
  double tx,ty,tz,a,b,c;

  if ( ! mol ) return -1;

  ucount = 0;
  hcount = 0;
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (hasharray_index(mol->segment_hash, seg->segid) == HASHARRAY_FAIL) 
      continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      for ( atom = res->atoms; atom; atom = atom->next ) {
        if ( atom->xyz_state != TOPO_MOL_XYZ_SET ) {
          ++ucount;
          if ( atom->mass > 2.5 ) ++hcount;
        }
      }
    }
  }
  sprintf(msg,"Warning: guessing coordinates for %d atoms (%d non-hydrogen)",
						ucount, hcount);
  topo_mol_log_error(mol,msg);

  uatoms = (topo_mol_atom_t**) malloc(ucount*sizeof(topo_mol_atom_t*));
  if ( ! uatoms ) return -2;
  ucount = 0;
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (hasharray_index(mol->segment_hash, seg->segid) == HASHARRAY_FAIL) 
      continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      for ( atom = res->atoms; atom; atom = atom->next ) {
        if ( atom->xyz_state != TOPO_MOL_XYZ_SET ) uatoms[ucount++] = atom;
      }
    }
  }

  for ( i=0; i<ucount; ++i ) uatoms[i]->xyz_state = TOPO_MOL_XYZ_VOID;

  /* everything below based on atom 4 unknown, all others known */

  /* from the CHARMM docs:

    Normal IC table entry:
                I
                 \
                  \
                   J----K
                         \
                          \
                           L
        values (Rij),(Tijk),(Pijkl),(Tjkl),(Rkl)

    Improper type of IC table entry:
                I        L
                 \     /
                  \   /
                   *K
                   |
                   |
                   J
        values (Rik),(Tikj),(Pijkl),T(jkl),(Rkl)

  */

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

  gcount = 1;
  okwild = 0;
  wcount = 0;
  hcount = 0;
  while ( gcount || ! okwild ) {
   if ( gcount == 0 ) { if ( okwild ) break; else okwild = 1; }
   gcount = 0;
   for ( i=0; i<ucount; ++i ) { atom = uatoms[i];
    if ( atom->xyz_state != TOPO_MOL_XYZ_VOID ) continue;
    for ( conf = atom->conformations; conf;
		conf = topo_mol_conformation_next(conf,atom) ) {
      if ( conf->del ) continue;
      else if ( conf->atom[0] == atom &&
		conf->atom[1]->xyz_state != TOPO_MOL_XYZ_VOID &&
		conf->atom[2]->xyz_state != TOPO_MOL_XYZ_VOID &&
		conf->atom[3]->xyz_state != TOPO_MOL_XYZ_VOID ) {
        if ( conf->improper ) {
          a1 = conf->atom[3]; a2 = conf->atom[1]; a3 = conf->atom[2];
          dist34 = conf->dist12;
          angle234 = conf->angle123 * (M_PI/180.0);
          dihedral = -1.0 * conf->dihedral * (M_PI/180.0);
        } else {
          a1 = conf->atom[3]; a2 = conf->atom[2]; a3 = conf->atom[1];
          dist34 = conf->dist12;
          angle234 = conf->angle123 * (M_PI/180.0);
          dihedral = conf->dihedral * (M_PI/180.0);
        } 
      } 
      else if ( conf->atom[3] == atom &&
		conf->atom[2]->xyz_state != TOPO_MOL_XYZ_VOID &&
		conf->atom[1]->xyz_state != TOPO_MOL_XYZ_VOID &&
		conf->atom[0]->xyz_state != TOPO_MOL_XYZ_VOID ) {
        if ( conf->improper ) {
          a1 = conf->atom[0]; a2 = conf->atom[1]; a3 = conf->atom[2];
          dist34 = conf->dist34;
          angle234 = conf->angle234 * (M_PI/180.0);
          dihedral = conf->dihedral * (M_PI/180.0);
        } else {
          a1 = conf->atom[0]; a2 = conf->atom[1]; a3 = conf->atom[2];
          dist34 = conf->dist34;
          angle234 = conf->angle234 * (M_PI/180.0);
          dihedral = conf->dihedral * (M_PI/180.0);
        } 
      } 
      else continue;

      gwild = 0;
      if ( dist34 == 0.0 ) { dist34 = 1.0; gwild = 1; }
      if ( angle234 == 0.0 ) { angle234 = 109.0*M_PI/180.0; gwild = 1; }

      r12x = a2->x - a1->x;
      r12y = a2->y - a1->y;
      r12z = a2->z - a1->z;
      r23x = a3->x - a2->x;
      r23y = a3->y - a2->y;
      r23z = a3->z - a2->z;
      a = sqrt(r23x*r23x + r23y*r23y + r23z*r23z);
      if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
      ix = a * r23x;
      iy = a * r23y;
      iz = a * r23z;
      tx = r12y*r23z - r12z*r23y;
      ty = r12z*r23x - r12x*r23z;
      tz = r12x*r23y - r12y*r23x;
      a = sqrt(tx*tx + ty*ty + tz*tz);
      if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
      kx = a * tx;
      ky = a * ty;
      kz = a * tz;
      tx = ky*iz - kz*iy;
      ty = kz*ix - kx*iz;
      tz = kx*iy - ky*ix;
      a = sqrt(tx*tx + ty*ty + tz*tz);
      if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
      jx = a * tx;
      jy = a * ty;
      jz = a * tz;
      a = -1.0 * dist34 * cos(angle234);
      b = dist34 * sin(angle234) * cos(dihedral);
      c = dist34 * sin(angle234) * sin(dihedral);

      if ( gwild && ! okwild ) continue;
      if ( okwild ) {
        ++wcount;
        if ( atom->mass > 2.5 ) ++hcount;
      }

      atom->x = a3->x + a * ix + b * jx + c * kx;
      atom->y = a3->y + a * iy + b * jy + c * ky;
      atom->z = a3->z + a * iz + b * jz + c * kz;
      atom->xyz_state = TOPO_MOL_XYZ_GUESS; ++gcount;
    }
   }
  }

  for ( i=0; i<ucount; ++i ) { atom = uatoms[i];
    if ( atom->xyz_state != TOPO_MOL_XYZ_VOID ) continue;

    /* pick heaviest atom we are bonded to (to deal with water) */
    a1 = 0;
    for ( bondtmp = atom->bonds; bondtmp;
		bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
      if ( bondtmp->atom[0] == atom ) a2 = bondtmp->atom[1];
      else a2 = bondtmp->atom[0];
      if ( a1 == 0 || a2->mass > a1->mass ) a1 = a2;
    }
    if ( a1 == 0 || a1->xyz_state == TOPO_MOL_XYZ_VOID ) continue;
    atom = a1;

    /* find all bonded atoms known and unknown coordinates */
    nk = 0;  nu = 0;
    for ( bondtmp = atom->bonds; bondtmp;
		bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
      if ( bondtmp->del ) continue;
      if ( bondtmp->atom[0] == atom ) a2 = bondtmp->atom[1];
      else a2 = bondtmp->atom[0];
      if ( a2->xyz_state == TOPO_MOL_XYZ_VOID ) {
        if ( nu < 4 ) ua[nu++] = a2;
      } else {
        if ( nk < 4 ) ka[nk++] = a2;
      }
    }

    if ( nu + nk > 4 ) continue;  /* no intuition beyond this case */

    if ( nk == 0 ) {  /* not bonded to any known atoms */
      a1 = ua[0];
      a1->x = atom->x + 1.0;
      a1->y = atom->y;
      a1->z = atom->z;
      a1->xyz_state = TOPO_MOL_XYZ_GUESS;
      ++gcount;  ++wcount;
      if ( a1->mass > 2.5 ) ++hcount;
      continue;
    }

    if ( nk == 1 ) {  /* bonded to one known atom */
      a1 = ka[0];
      ix = a1->x - atom->x;
      iy = a1->y - atom->y;
      iz = a1->z - atom->z;
      a = sqrt(ix*ix+iy*iy+iz*iz);
      if ( a ) a = 1.0 / a;  else continue;
      ix *= a; iy *= a; iz *= a;
      jx = -1.0 * iy;  jy = ix;  jz = 0;
      if ( jx*jx + jy*jy + jz*jz < 0.1 ) {
        jx = 0;  jy = -1.0 * iz;  jz = iy;
      }
      a = sqrt(jx*jx+jy*jy+jz*jz);
      if ( a ) a = 1.0 / a;  else continue;
      jx *= a; jy *= a; jz *= a;
      angle234 = 109.0*M_PI/180.0;
      a = cos(angle234);
      b = sin(angle234);
      a2 = ua[0];
      a2->x = atom->x + a * ix + b * jx;
      a2->y = atom->y + a * iy + b * jy;
      a2->z = atom->z + a * iz + b * jz;
      a2->xyz_state = TOPO_MOL_XYZ_GUESS;
      ++gcount;  ++wcount;
      if ( a2->mass > 2.5 ) ++hcount;
      continue;
    }

    if ( nk == 2 ) {  /* bonded to two known atoms */
      a1 = ka[0];
      ix = a1->x - atom->x;
      iy = a1->y - atom->y;
      iz = a1->z - atom->z;
      a = sqrt(ix*ix+iy*iy+iz*iz);
      if ( a ) a = 1.0 / a;  else continue;
      ix *= a; iy *= a; iz *= a;
      jx = ix;  jy = iy;  jz = iz;
      a1 = ka[1];
      ix = a1->x - atom->x;
      iy = a1->y - atom->y;
      iz = a1->z - atom->z;
      a = sqrt(ix*ix+iy*iy+iz*iz);
      if ( a ) a = 1.0 / a;  else continue;
      ix *= a; iy *= a; iz *= a;
      kx = jx - ix;  ky = jy - iy;  kz = jz - iz;
      jx += ix;  jy += iy;  jz += iz;
      a = sqrt(jx*jx+jy*jy+jz*jz);
      if ( a ) a = 1.0 / a;  else continue;
      jx *= a; jy *= a; jz *= a;
      if ( nu == 1 ) {  /* one unknown atom */
        a2 = ua[0];
        a2->x = atom->x - jx;
        a2->y = atom->y - jy;
        a2->z = atom->z - jz;
        a2->xyz_state = TOPO_MOL_XYZ_GUESS;
        ++gcount;  ++wcount;
        if ( a2->mass > 2.5 ) ++hcount;
      } else {  /* two unknown atoms */
        a = sqrt(kx*kx+ky*ky+kz*kz);
        if ( a ) a = 1.0 / a;  else continue;
        kx *= a; ky *= a; kz *= a;
        ix = jy*kz - jz*ky;
        iy = jz*kx - jx*kz;
        iz = jx*ky - jy*kx;
        a = sqrt(ix*ix+iy*iy+iz*iz);
        if ( a ) a = 1.0 / a;  else continue;
        ix *= a; iy *= a; iz *= a;
        angle234 = (180.0-0.5*109.0)*M_PI/180.0;
        a = sin(angle234);
        b = cos(angle234);
        a1 = ua[0];
        a2 = ua[1];
        a1->x = atom->x + a * ix + b * jx;
        a1->y = atom->y + a * iy + b * jy;
        a1->z = atom->z + a * iz + b * jz;
        a2->x = atom->x - a * ix + b * jx;
        a2->y = atom->y - a * iy + b * jy;
        a2->z = atom->z - a * iz + b * jz;
        a1->xyz_state = TOPO_MOL_XYZ_GUESS;
        ++gcount;  ++wcount;
        if ( a1->mass > 2.5 ) ++hcount;
        a2->xyz_state = TOPO_MOL_XYZ_GUESS;
        ++gcount;  ++wcount;
        if ( a2->mass > 2.5 ) ++hcount;
      }
      continue;
    }

    if ( nk == 3 ) {  /* bonded to three known atoms */
      a1 = ka[0];
      ix = a1->x - atom->x;
      iy = a1->y - atom->y;
      iz = a1->z - atom->z;
      a = sqrt(ix*ix+iy*iy+iz*iz);
      if ( a ) a = 1.0 / a;  else continue;
      ix *= a; iy *= a; iz *= a;
      jx = ix;  jy = iy;  jz = iz;
      a1 = ka[1];
      ix = a1->x - atom->x;
      iy = a1->y - atom->y;
      iz = a1->z - atom->z;
      a = sqrt(ix*ix+iy*iy+iz*iz);
      if ( a ) a = 1.0 / a;  else continue;
      ix *= a; iy *= a; iz *= a;
      jx += ix;  jy += iy;  jz += iz;
      a1 = ka[2];
      ix = a1->x - atom->x;
      iy = a1->y - atom->y;
      iz = a1->z - atom->z;
      a = sqrt(ix*ix+iy*iy+iz*iz);
      if ( a ) a = 1.0 / a;  else continue;
      ix *= a; iy *= a; iz *= a;
      jx += ix;  jy += iy;  jz += iz;
      a = sqrt(jx*jx+jy*jy+jz*jz);
      if ( a ) a = 1.0 / a;  else continue;
      a2 = ua[0];
      a2->x = atom->x - a * jx;
      a2->y = atom->y - a * jy;
      a2->z = atom->z - a * jz;
      a2->xyz_state = TOPO_MOL_XYZ_GUESS;
      ++gcount;  ++wcount;
      if ( a2->mass > 2.5 ) ++hcount;
      continue;
    }

  }

  gcount = 0;
  for ( i=0; i<ucount; ++i ) {
    if ( uatoms[i]->xyz_state == TOPO_MOL_XYZ_VOID ) ++gcount;
  }
  if ( wcount ) {
    sprintf(msg,"Warning: poorly guessed coordinates for %d atoms (%d non-hydrogen)", wcount, hcount);
    topo_mol_log_error(mol,msg);
  }
  if ( gcount ) {
    sprintf(msg,"Warning: failed to guess coordinates for %d atoms",gcount);
    topo_mol_log_error(mol,msg);
  }

  free((void*)uatoms);

  return 0;
}


