/**********************************************************************
 *         Copyright (C) 1995  The Board of Trustees of               *
 *                             the University of Illinois             *
 *                                                                    *
 *  This file is part of the RAPP software package, a library and     *
 *  associated programs for coordinating client/server applications.  *
 *                                                                    *
 **********************************************************************/

#ifndef MDCOMM_H
#define MDCOMM_H

#ifdef MDCOMM

#include <rapp.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef int Bool;

#define VMD_NAME_LEN 4                  /* Length of name fields for VMD */

/* This is all the static data that needs to be transferred to VMD */

typedef struct vmd_static_data
{
  int nameLen;                  /*  Length of each atom name             */
  int numAtoms;			/*  Total number of atoms                */
  int numAtomNames;		/*  Number of unique atom names          */
  int numAtomTypes;		/*  Number of unique atom types          */
  int numBonds;			/*  Number of linear bonds               */
  char *atomNames;		/*  List of unique atom names            */
  int *atomNameIndexes;		/*  Indexes for each atom into the       */
				/*  atomNames array                      */
  char *atomTypes;		/*  List of unique atom types            */
  int *atomTypeIndexes;		/*  Indexes for each atom into the       */
                                /*  atomTypes array                      */
  int *bonds;			/*  List of bonded atoms                 */
  char *resIds;			/*  Residue ids for each atom            */
  char *resNames;               /*  Residue names for each atom          */
  float *radii;			/*  Calculated vdw radii                 */
  float *charge;		/*  Charge for each atom                 */
  float *mass;			/*  Mass for each atom                   */
  float *occupancy;		/*  Occupancy for each atom              */
  float *beta;			/*  Beta coupling for each atom          */
  char *segIds;			/*  Segment IDs for each atom            */
  int maxNumPatches;            /*  Max number of patches                */

  /******** The remaining structure members are NOT used by NAMD *********/
  /********    They are included for VMD compatibility reasons   *********/
  int nHBondDonors;             /*  List of hydrogen bond donors         */
  int *donor;                   /*  List of h-donors (index into atoms)  */
  int *donorh;                  /*  List of associated h-donors          */
  int nHBondAcceptors;          /*  List of hydrogen bond acceptors      */
  int *acceptor;                /*  List of h-acceptors (idx into atoms) */
  int *acceptora;               /*  List of associated h-acceptors       */
  /***********************************************************************/
} VmdStaticData;

/* This structure is used to transfer dynamic data to VMD */

typedef struct vmd_dyn_data
{
  Bool energiesArrived;         /*  Flag TRUE->energies have been set    */
  Bool coordsArrived;           /*  Flag TRUE->coords have been set      */

  int *timestep;                /*  Current timestep                     */
  float *elapsed_time;          /*  Total elapsed simulation time (?)    */
  float *elapsed_cpu;	        /*  Total CPU time.  Not yet implemented */
  float *T;                     /*  Temperature                          */
  float *Epot;                  /*  Potential energy                     */
  float *Etot;                  /*  Total energy                         */
  float *Evdw;                  /*  Van der Waals energy                 */
  float *Eelec;                 /*  Electrostatic energy                 */
  float *Esup;                  /*  Not used by NAMD                     */
  float *Ehbo;                  /*  Not used by NAMD                     */
  float *Ebond;                 /*  Linear bond energy                   */
  float *Eangle;                /*  Angle bond energy                    */
  float *Edihe;                 /*  Dihedral bond energy                 */
  float *Eimpr;                 /*  Improper bond energy                 */
  float *X;                     /*  Array of X coordinates               */
  float *Y;                     /*  Array of Y coordinates               */
  float *Z;                     /*  Array of Z coordinates               */
  int *numPatches;              /*  Number of patches is simulation      */
  float *pXOrigins;             /*  Origins of patches                   */
  float *pYOrigins;             /*  Origins of patches                   */
  float *pZOrigins;             /*  Origins of patches                   */
  float *patchLength;           /*  Length of each patch                 */
  float *patchWidth;            /*  Length of each patch                 */
  float *patchHeight;           /*  Length of each patch                 */
  float *patchAtomNums;         /*  Number of atoms on each patch        */
  float *patchLoads;            /*  Number of atoms on each patch        */
  float *patchNode;             /*  Number of atoms on each patch        */
} VmdDynData;


enum mdcomm_tags {
  XFER_NAMELEN,
  XFER_NATOMS,
  XFER_NATOMNAMES,
  XFER_NATOMTYPES,
  XFER_NHBONDDONORS,
  XFER_NHBONDACCEPTORS,
  XFER_NBONDS,
  XFER_ATOMNAMELIST,
  XFER_ATOMNAMEIDX,
  XFER_ATOMTYPELIST,
  XFER_ATOMTYPEIDX,
  XFER_DONOR,
  XFER_DONORH,
  XFER_ACCEPTOR,
  XFER_ACCEPTORA,
  XFER_RESIDS,
  XFER_RESNAMES,
  XFER_RADII,
  XFER_BONDS,
  XFER_CHARGE,
  XFER_MASS,
  XFER_OCCUPANCY,
  XFER_BETA,
  XFER_SEGIDS,
  XFER_MAXPATCHES,
  XFER_STEP,
  XFER_DT,
  XFER_CPU,
  XFER_TEMP,
  XFER_EPOT,
  XFER_ETOT,
  XFER_EVDW,
  XFER_EELEC,
  XFER_ESUP,
  XFER_EHBO,
  XFER_EBOND,
  XFER_EANGLE,
  XFER_EDIHED,
  XFER_EIMPR,
  XFER_XCOORDS,
  XFER_YCOORDS,
  XFER_ZCOORDS,
  XFER_NPATCHES,
  XFER_PATCHX,
  XFER_PATCHY,
  XFER_PATCHZ,
  XFER_PATCHLENGTH,
  XFER_PATCHWIDTH,
  XFER_PATCHHEIGHT,
  XFER_PATCHLOAD,
  XFER_PATCHATOMS,
  XFER_PATCHNODE
};
typedef enum mdcomm_tags mdcomm_tags_t;

/* 
 * RAPP MDComm data transfer functions
 */

int mdcomm_send_static(rapp_active_socket_t *, void *);
int mdcomm_recv_static(rapp_active_socket_t *, void *);
int mdcomm_send_dynamic(rapp_active_socket_t *, void *, void *);
int mdcomm_recv_dynamic(rapp_active_socket_t *, void *, void *);

/*
 * RAPP MDComm application instantiation function
 */
int mdcomm_exec(rapp_appd_handle_t , rapp_paramlist_t *);

/*
 * RAPP MDComm dynamic data memory requirement utility function
 */
int mdcomm_mem_calc(void *);

/*
 * RAPP MDComm shared memory mapping function
 */
int mdcomm_map_shmem(void *, void *, shmaddr_t);


#ifdef __cplusplus
}
#endif

#endif /* MDCOMM */

#endif /* MDCOMM_H */
