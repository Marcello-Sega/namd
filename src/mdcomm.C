/**********************************************************************
 *         Copyright (C) 1995  The Board of Trustees of               *
 *                             the University of Illinois             *
 *                                                                    *
 *  This file is part of the RAPP software package, a library and     *
 *  associated programs for coordinating client/server applications.  *
 *                                                                    *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <pwd.h>
#include <string.h>

#include "mdcomm.h"
//#include "mdcomm_user_commands.h"

#ifdef MDCOMM
#include <rapp.h>

typedef struct heap_list{
  void *mem;
  struct heap_list *next;
} heap_list;
static struct heap_list *list = NULL, *tail=NULL;
static int add_to_list(void *);
static void delete_list(int);

#ifdef DEBUG
int recv_error(int, int, int);
#endif

/* 
 * RAPP MDComm data transfer functions
 */

int
mdcomm_send_static(rapp_active_socket_t *sock, void *static_ptr)
{
  struct vmd_static_data *sdata = (struct vmd_static_data *) static_ptr;
  int dummy;

  if (rapp_send(sock, XFER_NAMELEN, &sdata->nameLen, 
                sizeof(sdata->nameLen), RAPP_INT) != sizeof(sdata->nameLen))
    return -1;

  if (rapp_send(sock, XFER_NATOMS, &sdata->numAtoms, 
                sizeof(sdata->numAtoms), RAPP_INT) != sizeof(sdata->numAtoms))
    return -1;

  if (rapp_send(sock, XFER_NATOMNAMES, &sdata->numAtomNames, 
                sizeof(sdata->numAtomNames), RAPP_INT)
      != sizeof(sdata->numAtomNames))
    return -1;

  if (rapp_send(sock, XFER_NATOMTYPES, &sdata->numAtomTypes, 
                sizeof(sdata->numAtomTypes), RAPP_INT)
      != sizeof(sdata->numAtomTypes))
    return -1;

  /* Hydrogen bonds not implemented in NAMD - following two calls
     to rapp_send reflect this */

  dummy = 0;
  if (rapp_send(sock, XFER_NHBONDDONORS, &dummy, 
                sizeof(dummy), RAPP_INT) != sizeof(dummy))
    return -1;

  if (rapp_send(sock, XFER_NHBONDACCEPTORS, &dummy, 
                sizeof(dummy), RAPP_INT) != sizeof(dummy))
    return -1;

  if (rapp_send(sock, XFER_NBONDS, &sdata->numBonds, 
                sizeof(sdata->numBonds), RAPP_INT) != sizeof(sdata->numBonds))
    return -1;

  if (rapp_send(sock, XFER_ATOMNAMELIST, sdata->atomNames, 
                sdata->numAtomNames * sdata->nameLen, RAPP_BYTE)
      != sdata->numAtomNames * sdata->nameLen)
    return -1;

  if (rapp_send(sock, XFER_ATOMNAMEIDX, sdata->atomNameIndexes, 
                sdata->numAtoms * sizeof(int), RAPP_INT)
      != sdata->numAtoms * sizeof(int))
    return -1;

  if (rapp_send(sock, XFER_ATOMTYPELIST, sdata->atomTypes, 
                sdata->numAtomTypes * sdata->nameLen, RAPP_BYTE)
      != sdata->numAtomTypes * sdata->nameLen)
    return -1;

  if (rapp_send(sock, XFER_ATOMTYPEIDX, sdata->atomTypeIndexes, 
                sdata->numAtoms * sizeof(int), RAPP_INT)
      != sdata->numAtoms * sizeof(int))
    return -1;


  /* Hydrogen bonds not implemented in NAMD - following four calls
     to rapp_send reflect this */

  if (rapp_send(sock, XFER_DONOR, NULL, 0, RAPP_INT) != 0)
    return -1;

  if (rapp_send(sock, XFER_DONORH, NULL, 0, RAPP_INT) != 0)
    return -1;

  if (rapp_send(sock, XFER_ACCEPTOR, NULL, 0, RAPP_INT) != 0)
    return -1;

  if (rapp_send(sock, XFER_ACCEPTORA, NULL, 0, RAPP_INT) != 0)
    return -1;

  if (rapp_send(sock, XFER_RESIDS, sdata->resIds, 
                sdata->numAtoms * sdata->nameLen, RAPP_BYTE)
      != sdata->numAtoms * sdata->nameLen)
    return -1;

  if (rapp_send(sock, XFER_RESNAMES, sdata->resNames, 
                sdata->numAtoms * sdata->nameLen, RAPP_BYTE)
      != sdata->numAtoms * sdata->nameLen)
    return -1;

  if (rapp_send(sock, XFER_RADII, sdata->radii, 
                sdata->numAtoms * sizeof(float), RAPP_FLOAT)
      != sdata->numAtoms * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_BONDS, sdata->bonds, 
                sdata->numBonds * 2 * sizeof(int), RAPP_INT)
      != sdata->numBonds * 2 * sizeof(int))
    return -1;

  if (rapp_send(sock, XFER_CHARGE, sdata->charge, 
                sdata->numAtoms * sizeof(float), RAPP_FLOAT)
      != sdata->numAtoms * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_MASS, sdata->mass, 
                sdata->numAtoms * sizeof(float), RAPP_FLOAT)
      != sdata->numAtoms * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_OCCUPANCY, sdata->occupancy, 
                sdata->numAtoms * sizeof(float), RAPP_FLOAT)
      != sdata->numAtoms * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_BETA, sdata->beta, 
                sdata->numAtoms * sizeof(float), RAPP_FLOAT)
      != sdata->numAtoms * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_SEGIDS, sdata->segIds, 
                sdata->numAtoms * sdata->nameLen, RAPP_BYTE)
      != sdata->numAtoms * sdata->nameLen)
    return -1;

  if (rapp_send(sock, XFER_MAXPATCHES, &sdata->maxNumPatches, 
                sizeof(sdata->maxNumPatches), RAPP_INT)
      != sizeof(sdata->maxNumPatches))
    return -1;

  return 0;
}


int
mdcomm_recv_static(rapp_active_socket_t *sock, void *static_ptr)
{
  struct vmd_static_data *sdata = (struct vmd_static_data *) static_ptr;
  int tag;

  if ((rapp_recv(sock, &tag, &sdata->nameLen, RAPP_INT)
       != sizeof(sdata->nameLen)) || (tag != XFER_NAMELEN))
    return -1;

  if ((rapp_recv(sock, &tag, &sdata->numAtoms, RAPP_INT)
       != sizeof(sdata->numAtoms)) || (tag != XFER_NATOMS))
    return -1;

  if ((rapp_recv(sock, &tag, &sdata->numAtomNames, RAPP_INT)
       != sizeof(sdata->numAtomNames)) || (tag != XFER_NATOMNAMES))
    return -1;

  if ((rapp_recv(sock, &tag, &sdata->numAtomTypes, RAPP_INT)
       != sizeof(sdata->numAtomTypes)) || (tag != XFER_NATOMTYPES))
    return -1;

  if ((rapp_recv(sock, &tag, &sdata->nHBondDonors, RAPP_INT)
       != sizeof(sdata->nHBondDonors)) || (tag != XFER_NHBONDDONORS))
    return -1;

  if ((rapp_recv(sock, &tag, &sdata->nHBondAcceptors, RAPP_INT) 
       != sizeof(sdata->nHBondAcceptors)) || (tag !=  XFER_NHBONDACCEPTORS))
    return -1;

  if ((rapp_recv(sock, &tag, &sdata->numBonds, RAPP_INT)
       != sizeof(sdata->numBonds)) || (tag != XFER_NBONDS))
    return -1;

  if ((sdata->atomNames =
       (char *) malloc(sdata->numAtomNames * sdata->nameLen)) == NULL)
    return -1;
  add_to_list(sdata->atomNames);

  if ((rapp_recv(sock, &tag, sdata->atomNames, RAPP_BYTE)
       != sdata->numAtomNames * sdata->nameLen) || (tag != XFER_ATOMNAMELIST))
    return -1;

  if ((sdata->atomNameIndexes = (int *) malloc(sdata->numAtoms * sizeof(int)))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->atomNameIndexes);
  if ((rapp_recv(sock, &tag, sdata->atomNameIndexes, RAPP_INT)
       != sdata->numAtoms * sizeof(int)) || (tag != XFER_ATOMNAMEIDX)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->atomTypes =
       (char *) malloc(sdata->numAtomTypes * sdata->nameLen))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->atomTypes);
  if ((rapp_recv(sock, &tag, sdata->atomTypes, RAPP_BYTE)
      != sdata->numAtomTypes * sdata->nameLen) || (tag != XFER_ATOMTYPELIST)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->atomTypeIndexes =
       (int *) malloc(sdata->numAtoms * sizeof(int))) == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->atomTypeIndexes);
  if ((rapp_recv(sock, &tag, sdata->atomTypeIndexes, RAPP_INT)
      != sdata->numAtoms * sizeof(int)) || (tag != XFER_ATOMTYPEIDX)) {
    delete_list(1);
    return -1;
  }

  if (sdata->nHBondDonors> 0) {
    if ((sdata->donor = (int *) malloc(sdata->nHBondDonors * sizeof(int)))
	== NULL) {
      delete_list(1);
      return -1;
    }
    add_to_list(sdata->donor);
  }
  if ((rapp_recv(sock, &tag, sdata->donor, RAPP_INT)
       != sdata->nHBondDonors * sizeof(int)) || (tag != XFER_DONOR)) {
    delete_list(1);
    return -1;
  }

  if (sdata->nHBondDonors > 0) {
    if ((sdata->donorh = (int *) malloc(sdata->nHBondDonors * sizeof(int)))
	== NULL) {
      delete_list(1);
      return -1;
    }
    add_to_list(sdata->donorh);
  }
  if ((rapp_recv(sock, &tag, sdata->donorh, RAPP_INT)
       != sdata->nHBondDonors * sizeof(int)) || (tag != XFER_DONORH)) {
    delete_list(1);
    return -1;
  }

  if (sdata->nHBondAcceptors > 0) {
    if ((sdata->acceptor =
	 (int *) malloc(sdata->nHBondAcceptors * sizeof(int)))
	== NULL) {
      delete_list(1);
      return -1;
    }
    add_to_list(sdata->acceptor);
  }
  if ((rapp_recv(sock, &tag, sdata->acceptor, RAPP_INT)
       != sdata->nHBondAcceptors * sizeof(int)) || (tag != XFER_ACCEPTOR)) {
    delete_list(1);
    return -1;
  }

  if (sdata->nHBondAcceptors > 0) {
    if ((sdata->acceptora =
	 (int *) malloc(sdata->nHBondAcceptors * sizeof(int)))
	== NULL) {
      delete_list(1);
      return -1;
    }
    add_to_list(sdata->acceptora);
  }
  if ((rapp_recv(sock, &tag, sdata->acceptora, RAPP_INT)
       != sdata->nHBondAcceptors * sizeof(int)) || (tag != XFER_ACCEPTORA)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->resIds = (char *) malloc(sdata->numAtoms * sdata->nameLen))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->resIds);
  if ((rapp_recv(sock, &tag, sdata->resIds, RAPP_BYTE)
       != sdata->numAtoms * sdata->nameLen) || (tag != XFER_RESIDS)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->resNames = (char *) malloc(sdata->numAtoms * sdata->nameLen))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->resNames);
  if ((rapp_recv(sock, &tag, sdata->resNames,  RAPP_BYTE)
       != sdata->numAtoms * sdata->nameLen) || (tag != XFER_RESNAMES)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->radii = (float *) malloc(sdata->numAtoms * sizeof(float)))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->radii);
  if ((rapp_recv(sock, &tag, sdata->radii, RAPP_FLOAT)
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_RADII)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->bonds = (int *) malloc(sdata->numBonds * 2 * sizeof(int)))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->bonds);
  if ((rapp_recv(sock, &tag, sdata->bonds, RAPP_INT)
       != sdata->numBonds * 2 * sizeof(int)) || (tag != XFER_BONDS)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->charge = (float *) malloc(sdata->numAtoms * sizeof(float)))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->charge);
  if ((rapp_recv(sock, &tag, sdata->charge, RAPP_FLOAT)
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_CHARGE)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->mass = (float *) malloc(sdata->numAtoms * sizeof(float)))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->mass);
  if ((rapp_recv(sock, &tag, sdata->mass, RAPP_FLOAT)
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_MASS)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->occupancy = (float *) malloc(sdata->numAtoms * sizeof(float)))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->occupancy);
  if ((rapp_recv(sock, &tag, sdata->occupancy, RAPP_FLOAT)
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_OCCUPANCY)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->beta = (float *) malloc(sdata->numAtoms * sizeof(float)))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->beta);
  if ((rapp_recv(sock, &tag, sdata->beta, RAPP_FLOAT)
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_BETA)) {
    delete_list(1);
    return -1;
  }

  if ((sdata->segIds = (char *) malloc(sdata->numAtoms * sdata->nameLen))
      == NULL) {
    delete_list(1);
    return -1;
  }
  add_to_list(sdata->segIds);
  if ((rapp_recv(sock, &tag, sdata->segIds, RAPP_BYTE)
       != sdata->numAtoms * sdata->nameLen) || (tag != XFER_SEGIDS)) {
    delete_list(1);
    return -1;
  }

  if ((rapp_recv(sock, &tag, &sdata->maxNumPatches, RAPP_INT)
       != sizeof(sdata->maxNumPatches)) || (tag != XFER_MAXPATCHES)) {
    delete_list(1);
    return -1;
  }

  delete_list(0);
  return 0;
}


int
mdcomm_send_dynamic(rapp_active_socket_t *sock,
		    void *static_ptr, void *dynamic_ptr)
{
  struct vmd_static_data *sdata = (struct vmd_static_data *) static_ptr;
  struct vmd_dyn_data *ddata = (struct vmd_dyn_data *) dynamic_ptr;
  float dummy = 0.0;

  if (rapp_send(sock, XFER_STEP, ddata->timestep, 
                sizeof(int), RAPP_INT) != sizeof(int)) {
    perror("send_dynamic");
    return -1;
  }

  if (rapp_send(sock, XFER_DT, ddata->elapsed_time, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_CPU, ddata->elapsed_cpu, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_TEMP, ddata->T, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_EPOT, ddata->Epot, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_ETOT, ddata->Etot, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_EVDW, ddata->Evdw, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_EELEC, ddata->Eelec, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_ESUP, &dummy, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_EHBO, &dummy, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_EBOND, ddata->Ebond, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_EANGLE, ddata->Eangle, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_EDIHED, ddata->Edihe, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_EIMPR, ddata->Eimpr, 
                sizeof(float), RAPP_FLOAT) != sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_XCOORDS, ddata->X, 
                sdata->numAtoms * sizeof(float), RAPP_FLOAT)
      != sdata->numAtoms * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_YCOORDS, ddata->Y, 
                sdata->numAtoms * sizeof(float), RAPP_FLOAT)
      != sdata->numAtoms * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_ZCOORDS, ddata->Z, 
                sdata->numAtoms * sizeof(float), RAPP_FLOAT)
      != sdata->numAtoms * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_NPATCHES, ddata->numPatches, 
                sizeof(int), RAPP_INT) != sizeof(int))
    return -1;

  if (rapp_send(sock, XFER_PATCHX, ddata->pXOrigins, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_PATCHY, ddata->pYOrigins, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_PATCHZ, ddata->pZOrigins, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_PATCHLENGTH, ddata->patchLength, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_PATCHWIDTH, ddata->patchWidth, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_PATCHHEIGHT, ddata->patchHeight, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_PATCHLOAD, ddata->patchLoads, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_PATCHATOMS, ddata->patchAtomNums, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  if (rapp_send(sock, XFER_PATCHNODE, ddata->patchNode, 
                *ddata->numPatches * sizeof(float), RAPP_FLOAT)
      != *ddata->numPatches * sizeof(float))
    return -1;

  return 0;
}

#ifdef DEBUG

int
mdcomm_recv_dynamic(rapp_active_socket_t *sock,
		    void *static_ptr, void *dynamic_ptr)
{
  struct vmd_static_data *sdata = (struct vmd_static_data *) static_ptr;
  struct vmd_dyn_data *ddata = (struct vmd_dyn_data *) dynamic_ptr;
  int tag;
  int ret;

  if (((ret = rapp_recv(sock, &tag, ddata->timestep, RAPP_INT))
       != sizeof(int)) || (tag != XFER_STEP)) {
    recv_error(tag, ret, XFER_STEP);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->elapsed_time, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_DT)) {
    recv_error(tag, ret, XFER_DT);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->elapsed_cpu, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_CPU)) {
    recv_error(tag, ret, XFER_CPU);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->T, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_TEMP)) {
    recv_error(tag, ret, XFER_TEMP);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Epot, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_EPOT)) {
    recv_error(tag, ret, XFER_EPOT);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Etot, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_ETOT)) {
    recv_error(tag, ret, XFER_ETOT);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Evdw, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_EVDW)) {
    recv_error(tag, ret, XFER_EVDW);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Eelec, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_EELEC)) {
    recv_error(tag, ret, XFER_EELEC);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Esup, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_ESUP)) {
    recv_error(tag, ret, XFER_ESUP);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Ehbo, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_EHBO)) {
    recv_error(tag, ret, XFER_EHBO);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Ebond, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_EBOND)) {
    recv_error(tag, ret, XFER_EBOND);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Eangle, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_EANGLE)) {
    recv_error(tag, ret, XFER_EANGLE);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Edihe, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_EDIHED)) {
    recv_error(tag, ret, XFER_EDIHED);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Eimpr, RAPP_FLOAT))
       != sizeof(float)) || (tag != XFER_EIMPR)) {
    recv_error(tag, ret, XFER_EIMPR);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->X, RAPP_FLOAT))
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_XCOORDS)) {
    recv_error(tag, ret, XFER_XCOORDS);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Y, RAPP_FLOAT))
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_YCOORDS)) {
    recv_error(tag, ret, XFER_YCOORDS);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->Z, RAPP_FLOAT))
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_ZCOORDS)) {
    recv_error(tag, ret, XFER_ZCOORDS);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->numPatches, RAPP_INT))
       != sizeof(int)) || (tag != XFER_NPATCHES)) {
    recv_error(tag, ret, XFER_NPATCHES);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->pXOrigins, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHX)) {
    recv_error(tag, ret, XFER_PATCHX);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->pYOrigins, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHY)) {
    recv_error(tag, ret, XFER_PATCHY);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->pZOrigins, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHZ)) {
    recv_error(tag, ret, XFER_PATCHZ);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->patchLength, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHLENGTH)) {
    recv_error(tag, ret, XFER_PATCHLENGTH);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->patchWidth, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHWIDTH)) {
    recv_error(tag, ret, XFER_PATCHWIDTH);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->patchHeight, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHHEIGHT)) {
    recv_error(tag, ret, XFER_PATCHHEIGHT);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->patchLoads, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHLOAD)) {
    recv_error(tag, ret, XFER_PATCHLOAD);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->patchAtomNums, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHATOMS)) {
    recv_error(tag, ret, XFER_PATCHATOMS);
    return -1;
  }

  if (((ret = rapp_recv(sock, &tag, ddata->patchNode, RAPP_FLOAT))
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHNODE)) {
    recv_error(tag, ret, XFER_PATCHNODE);
    return -1;
  }

  return 0;
}

int
recv_error(int recd_tag, int recd_bytes, int exp_tag)
{
  fprintf(stderr,"Expected tag %d (received %d), received %d bytes.\n",
	  exp_tag, recd_tag, recd_bytes);
  return 0;
}

#else

int
mdcomm_recv_dynamic(rapp_active_socket_t *sock,
		    void *static_ptr, void *dynamic_ptr)
{
  struct vmd_static_data *sdata = (struct vmd_static_data *) static_ptr;
  struct vmd_dyn_data *ddata = (struct vmd_dyn_data *) dynamic_ptr;
  int tag;

  if ((rapp_recv(sock, &tag, ddata->timestep, RAPP_INT)
       != sizeof(int)) || (tag != XFER_STEP))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->elapsed_time, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_DT))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->elapsed_cpu, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_CPU))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->T, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_TEMP))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Epot, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_EPOT))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Etot, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_ETOT))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Evdw, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_EVDW))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Eelec, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_EELEC))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Esup, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_ESUP))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Ehbo, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_EHBO))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Ebond, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_EBOND))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Eangle, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_EANGLE))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Edihe, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_EDIHED))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Eimpr, RAPP_FLOAT)
       != sizeof(float)) || (tag != XFER_EIMPR))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->X, RAPP_FLOAT)
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_XCOORDS))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Y, RAPP_FLOAT)
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_YCOORDS))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->Z, RAPP_FLOAT)
       != sdata->numAtoms * sizeof(float)) || (tag != XFER_ZCOORDS))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->numPatches, RAPP_INT)
       != sizeof(int)) || (tag != XFER_NPATCHES))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->pXOrigins, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHX))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->pYOrigins, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHY))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->pZOrigins, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHZ))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->patchLength, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHLENGTH))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->patchWidth, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHWIDTH))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->patchHeight, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHHEIGHT))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->patchLoads, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHLOAD))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->patchAtomNums, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHATOMS))
    return -1;

  if ((rapp_recv(sock, &tag, ddata->patchNode, RAPP_FLOAT)
       != *ddata->numPatches * sizeof(float)) || (tag != XFER_PATCHNODE))
    return -1;

  return 0;
}

#endif

/*
 * RAPP MDComm dynamic data memory requirement utility function
 */

int
mdcomm_mem_calc(void *static_ptr)
{
  struct vmd_static_data *sdata = (struct vmd_static_data *) static_ptr;
  int bytes_needed;

  bytes_needed =
    sdata->numAtoms * 3 * sizeof(float) +     /* Atom coordinates       */
    2 * sizeof(int) +                         /* Timestep, numPatches   */
    13 * sizeof(float) +                      /* Elapsed_time,          */
                                              /* elapsed_cpu,           */
                                              /* temperature,           */
                                              /* 10 energies            */
    sdata->maxNumPatches * 9 * sizeof(float); /* Patch X, Y, Z origins, */
                                              /* length, width, height, */
                                              /* load, atomNums, node   */

#ifdef DEBUG
  printf("MDCOMM MEM_CALC: need %d bytes of shared memory\n",bytes_needed);
#endif

  return bytes_needed;
}


/*
 * RAPP MDComm shared memory mapping function
 */

int
mdcomm_map_shmem(void *static_ptr, void *dynamic_ptr, shmaddr_t shmem)
{
  struct vmd_static_data *sdata = (struct vmd_static_data *) static_ptr;
  struct vmd_dyn_data *ddata = (struct vmd_dyn_data *) dynamic_ptr;
  int coordmem, patchmem;

  coordmem = sdata->numAtoms * sizeof(float);
  patchmem = sdata->maxNumPatches * sizeof(float);

  ddata->timestep = (int *) shmem;
  ddata->elapsed_time = (float *) (((char *) ddata->timestep) + sizeof(int));
  ddata->elapsed_cpu = (float *) ((char *) ddata->elapsed_time + sizeof(float));
  ddata->T = (float *) ((char *) ddata->elapsed_cpu + sizeof(float));
  ddata->Epot = (float *) ((char *) ddata->T + sizeof(float));
  ddata->Etot = (float *) ((char *) ddata->Epot + sizeof(float));
  ddata->Evdw = (float *) ((char *) ddata->Etot + sizeof(float));
  ddata->Eelec = (float *) ((char *) ddata->Evdw + sizeof(float));
  ddata->Esup = (float *) ((char *) ddata->Eelec + sizeof(float));
  ddata->Ehbo = (float *) ((char *) ddata->Esup + sizeof(float));
  ddata->Ebond = (float *) ((char *) ddata->Ehbo + sizeof(float));
  ddata->Eangle = (float *) ((char *) ddata->Ebond + sizeof(float));
  ddata->Edihe = (float *) ((char *) ddata->Eangle + sizeof(float));
  ddata->Eimpr = (float *) ((char *) ddata->Edihe + sizeof(float));
  ddata->X = (float *) ((char *) ddata->Eimpr + sizeof(float));
  ddata->Y = (float *) ((char *) ddata->X + coordmem);
  ddata->Z = (float *) ((char *) ddata->Y + coordmem);
  ddata->numPatches = (int *) ((char *) ddata->Z + coordmem);
  ddata->pXOrigins = (float *) ((char *) ddata->numPatches + sizeof(int));
  ddata->pYOrigins = (float *) ((char *) ddata->pXOrigins + patchmem);
  ddata->pZOrigins = (float *) ((char *) ddata->pYOrigins + patchmem);
  ddata->patchLength = (float *) ((char *) ddata->pZOrigins + patchmem);
  ddata->patchWidth = (float *) ((char *) ddata->patchLength + patchmem);
  ddata->patchHeight = (float *) ((char *) ddata->patchWidth + patchmem);
  ddata->patchAtomNums = (float *) ((char *) ddata->patchHeight + patchmem);
  ddata->patchLoads = (float *) ((char *) ddata->patchAtomNums + patchmem);
  ddata->patchNode = (float *) ((char *) ddata->patchLoads + patchmem);

  return 0;
}


/*
 * RAPP MDComm application instantiation function
 */

#define CONFIG_FILE "NAMD.conf"
#define CONFIG_FMT "%-16s = %s\n"

int
mdcomm_exec(rapp_appd_handle_t handle, rapp_paramlist_t *list)
{
  int i, multivalued;
  struct timeval tp;
  struct timezone tzp;
  char *login, *ptr;
  struct passwd *pwent;
  FILE *fd;

  /* This is a fix to force NAMD to work properly on
   * the Convex Exemplar, where the syntax of the command-line
   * arguments that NAMD accepts is different than on other machines
   */
#if !defined(EXEMPLAR) && !defined(CSPP)
  char *namd_argv[3];
#else
  char *namd_argv[4];
#endif

  if ((fd = fopen(CONFIG_FILE, "w")) == NULL) 
    return -1;

  gettimeofday(&tp, &tzp);

  if ((login = getlogin()) == NULL)
    if ((pwent = getpwuid(getuid())) == (struct passwd *) NULL)
      login = "unknown user";
    else
      login = pwent->pw_name;

  fprintf(fd, "# NAMD configuration file -- \n");
  fprintf(fd, "# Generated by MDComm for %s at ", login);
  time_t tmp_t = time(NULL);
  fprintf(fd, "%s", ctime(&tmp_t));
  fprintf(fd, "#\n");

  for (i=0; i<list->count; i++) {
    ptr = list->param[i]->val;
    multivalued = FALSE;
    do {
      if (*ptr == ' ') multivalued = TRUE;
    } while (*++ptr);
    if (multivalued) {
      ptr = list->param[i]->val;
      while ((ptr = strtok(ptr," ")) != NULL) {
	fprintf(fd, CONFIG_FMT, list->param[i]->keyword, ptr);
	ptr = NULL;
      }
    }
    else
      fprintf(fd, CONFIG_FMT, list->param[i]->keyword, list->param[i]->val);
  }
  fprintf(fd, CONFIG_FMT, "MDComm", "on");
  fclose(fd);

  ptr = strrchr(list->exe,'/');
  namd_argv[0] = (ptr == NULL ? list->exe : ptr+1);
#ifndef EXEMPLAR
  namd_argv[1] = CONFIG_FILE;
  namd_argv[2] = NULL;
#else
  namd_argv[1] = CONFIG_FILE;
  namd_argv[2] = "1";
  namd_argv[3] = NULL;
#endif

  execv(list->exe, namd_argv);

  return -1;
}


/* 
 * Some useful utility functions
 */

static int
add_to_list(void *ptr)
{
  if (list == NULL) {
    if ((list = tail = (struct heap_list *) malloc(sizeof(struct heap_list)))
	== NULL) 
      return -1;
  }
  else {
    if ((tail->next = (struct heap_list *) malloc(sizeof(struct heap_list)))
	== NULL) {
      delete_list(1);
      return -1;
    }
    tail = tail->next;
  }
  tail->mem = ptr;
  tail->next = NULL;
    
  return 0;
}

static void
delete_list(int flag)
{
  struct heap_list *tmp;

  while (list != NULL) {
    tmp = list->next;
    if (flag) free(list->mem);
    free(list);
    list = tmp;
  }    
  
  return;
}

#endif /* MDCOMM */

