/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
/* $Id: dpme2_paralibsW.c,v 1.1 1997/04/23 17:05:17 nealk Exp $
 */

/***********************************************************
* This file contain routines necessary to starting off 
* the pvm process and group. A pvm group of workers (slaves)
* is named (direct_solvers) and it performs the direct_sum
* in parallel.
*/

#include "dpme2.h"

/**********************************************************
* This routine is called by all processes and have sections
* that are executed only by the master process (group id =0)
* this section can be modified for spawning on the cray T3D
*/

void 
dpme_register(char *spawnme, char **args, PeInfo *pe_info)
{
  int i,mytid, cnode, nslaves;
  char ProgName[80];
#if 0
  char charbuff[80];
  FILE *infile;
#endif
 int *mygroupid, nproc,  *tidarray,  *ncube;
/*--------------------------------------------*/
  nproc=  pe_info->myproc.nprocs;
  ncube=  &(pe_info->myproc.ncube);
  mygroupid =   &(pe_info->myproc.node);
  tidarray = pe_info->inst_node;
/*--------------------------------------------*/
  /* setup also in main() ayt 3/97  */
  nslaves = nproc ;

  /* enroll in PVM */
  if ( (mytid= pvm_mytid()) <0)  { 
    printf (" * The instant is %d\n",mytid);
    fprintf (stderr,"error: enrolling PVM, exiting \n");
    exit(1);
  }; 
  
  /* Join the group; all the group of direct solvers */
  *mygroupid= pvm_joingroup(GROUP);  
  
#if VERBOSE
  printf (" * NODE %d :  instant is %d \n", *mygroupid, mytid);  
  sprintf(ProgName, "%s%s", WORKINGDIR , spawnme);
  printf (" %s \n", ProgName);
  fflush(stdout);
#endif
  /********** two ways to spawn, top is general, bottom for COW only ***********/               
  /* spawn is only valid for the  C.O.W. cluster, take out for T3D */
  if ( (*mygroupid==0) && (nslaves!=0) ){

#if 1
    cnode=pvm_spawn(ProgName, args, 0, "", nslaves,tidarray);
    if ((cnode!= nslaves)&&(nslaves!=0)) {
      fprintf(stderr,"Error: spawned %d of %d processes\n",cnode,nslaves);
      error_handler("Check your WORKINGDIR def");
    }
#else
    /* another method to spawn slaves using a file "host.active" that 
     * has a ranked list in ascending order of load  */
    if ( (infile = fopen("host.active","r") ) == NULL) {
      fprintf(stderr,"Unable to open the host.active file: !!! \n");
      exit(1); 
    }
    sprintf(ProgName, "%s%s", WORKINGDIR , spawnme);
    for (i=nslaves-1;i>=0; i--){
      fscanf(infile,"%s",charbuff);
      printf(" i will initiate node %s\n",charbuff);
      cnode=pvm_spawn(ProgName, args,1,charbuff,1,tidarray);  
      printf(" status= %d: The instant for node %d is %d > 0 ?!\n",cnode,i,tidarray[i]);  
    }
#endif
    /************************************************************************/
  } /* end if master */
    
  /* exit if not being used */
  if ( *mygroupid > nslaves ) {
    printf(" Slave Error Exiting now \n");
    fflush(stdout);
    pvm_lvgroup(GROUP); 
    pvm_exit();
    exit(0);
  }
  /* initialize tidarray to -1, it will be filled as needed by each PE 
     it will be indexed by group number. */
  for (i=0;i<=nproc+1;i++)
    tidarray[i]=-1;

  /* all should sync now */
  if(  (pvm_barrier(GROUP,nproc+1)) <0) {
    printf(" Error: at the   barrier \n");
    fflush(stdout);
    pvm_exit();
    exit(0);
  }
  *ncube =0;
  while ( (pow(2.0,*ncube) < nproc )) (*ncube)++;

} /* end dpme_register() */
/*********************************************************************/
/* Merging integer values across the parallel system, efficient for power of 2 config
*/
 
int merge_i(int *ivalue,  MYPROC *myproc, int *tidarray)
{
  int i1, i2;
  int itmp, ipartner, i;
  int ipartner_tid;
  
  i1 = myproc->ncube - 1;
  for (i = 0; i <= i1; ++i) {
    i2= pow(2.0,(double)i);
    ipartner = myproc->node ^ i2 ; /* XOR = ieor */
    /*********** PVM- WRITE ********************/
    pvm_initsend(PvmDataType);
    pvm_pkint(ivalue,1,1);
    ipartner_tid=  tidarray[ipartner];
    if (ipartner_tid <0) {
      ipartner_tid = pvm_gettid(GROUP,ipartner);
      tidarray[ipartner]= ipartner_tid;
    }
    pvm_send(ipartner_tid,MSG_100);
    /*********** PVM - READ ********************/
    pvm_recv(ipartner_tid,MSG_100);
    pvm_upkint(&itmp,1,1);
    *ivalue += itmp;
  }
  return 0;
} /* merge_i */

/********************************************************************/
/* Merging integer values across the parallel system                
 * This merge is intended for non-power of 2 config
 * where we simply make the master gets all the results 
*/
 
int merge_i2(int *ivalue,  MYPROC *myproc, int *tidarray)
{
  int itmp, master, i;
  int master_tid,num_procs;
   

  if (myproc->nprocs == (int)pow(2.0,myproc->ncube)){
    merge_i(ivalue,myproc,tidarray);
  }
  else{
    master = 0; /* ie the master */
    master_tid = tidarray[master];
    if (master_tid <0) {
      master_tid = pvm_gettid(GROUP,master);
      tidarray[master] = master_tid;
    }
    num_procs = myproc->nprocs;

    /*********** PVM- WRITE ********************/
    if (myproc->node!=master){
      pvm_initsend(PvmDataType);
      pvm_pkint(ivalue,1,1);
      pvm_send(master_tid,MSG_100);
      /* fprintf(stderr,"send in merge_i from %i to %i(tid=%i) \n",
	 myproc->node,master,master_tid); fflush(stderr); */
    }
    /*********** PVM - READ ********************/
    else{  
      for (i=1;i<num_procs;i++){
	pvm_recv(-1,MSG_100);
	pvm_upkint(&itmp,1,1);
	*ivalue += itmp;
      }
    }
  } /* 1st else */
  return 0;
} /* merge_i2 */
/*********************************************************************/
/* swap particles that are left my box with  ones entering it        */

int swap(int node, double *sbuf, int islen, int isnode, 
	 double *rbuf, int *irlen, int irnode, int *tidarray)
{
  
  int ipartner_tid;
  int mybuf;
  int numdouble;

  numdouble=0;
  if (isnode != node) {
    /*** PVM SEND list ***/
   
    pvm_initsend(PvmDataType);
    pvm_pkint(&islen,1,1);
    if (islen) {
#if DPME_DEBUG
      printf("\n $ SND %i TO %i Size(%i) \n",node, isnode,islen); fflush(stdout);
#endif
      pvm_pkdouble(sbuf,islen,1);
    }
    
    ipartner_tid=  tidarray[isnode];
    if (ipartner_tid <0) {
      ipartner_tid = pvm_gettid(GROUP,isnode);
      tidarray[isnode] =   ipartner_tid; 
    }
    pvm_send(ipartner_tid,MSG_200);

    /*** PVM RECV list ***/
    /* printf("$ RCV %i From %i \n",node, irnode); */
    ipartner_tid=  tidarray[irnode];
    if (ipartner_tid <0) {
      ipartner_tid = pvm_gettid(GROUP,irnode);
      tidarray[irnode]= ipartner_tid;
    }
    mybuf=pvm_recv(ipartner_tid,MSG_200);
    pvm_upkint(&numdouble,1,1);

#if DPME_DEBUG
    printf("$$ PE%d to rcv From PE%d Size(%d) \n",node, irnode,numdouble);fflush(stdout);
#endif

    if (numdouble) {     
      pvm_upkdouble(rbuf,numdouble,1);
    }
  }
  *irlen = numdouble; /* num of doubles */

#if 0
/*----------------------------------------------------------------*/
  {
  int i;
  printf("%%%%%%%%%%%%%%%%%%%%%Node %d rcvd%%%%%%%%%%%%%%%%%% \n",node);
  --rbuf; 
  for (i=0;i<numdouble/5.;i++)
    printf("(%f, %f, %f, %f, %f)  ",*(++rbuf), *(++rbuf), *(++rbuf), *(++rbuf),
	   *(++rbuf));
  }
/*----------------------------------------------------------------*/
#endif


  return 0;
} /* swap */
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/* Global reduction of a double, efficient for power of 2 configs */
int merge_d(double *dvalue, MYPROC myproc,int *tidarray, int msg_tag)
{
  int i,i1,i2,ipartner, ipartner_tid ;
  double tmp;

  i1 = myproc.ncube - 1;
  for (i = 0; i <= i1; ++i) {
    i2 = pow(2.0,(double)i);
    ipartner = myproc.node ^ i2; /* XOR */
    /*********** PVM- WRITE ********************/
    pvm_initsend(PvmDataType);
    pvm_pkdouble(dvalue,1,1);
    ipartner_tid=  tidarray[ipartner];
    if (ipartner_tid <0) {
      ipartner_tid = pvm_gettid(GROUP,ipartner);
      tidarray[ipartner]= ipartner_tid;
    }
    pvm_send(ipartner_tid,msg_tag);
    /*********** PVM - READ ********************/
    pvm_recv(ipartner_tid,msg_tag);
    pvm_upkdouble(&tmp,1,1);
    *dvalue += tmp;
  }
  return 0;
} /* merge_d*/
/*********************************************************************/
/* Global reduction of a double */
/* targeted for configs of non-power of 2 */
/* basically just let the master do the sum */

int merge_d2(double *dvalue, MYPROC myproc, int *tidarray, int msg_tag)
{
  int  master,i;
  double tmp;
  int master_tid,num_procs;

  /* this will allow for two outstanding merges to be recieved properly ;
   * otherwise, if u have 2 consecutive merge_d, mesgs may get mixed up,
   * eg, in adj_dir and self 
   */

  if (myproc.nprocs == (int)pow(2.0,myproc.ncube)){
#if VEBOSE
    printf("calling merge_d...\n");
#endif
    merge_d(dvalue,myproc,tidarray,msg_tag);
  }
  else{ 
    master = 0; /* ie the master */
    master_tid = tidarray[master];
    if (master_tid <0) { /* it will <0 due to initiaize */
      master_tid = pvm_gettid(GROUP,master);
      tidarray[master]=  master_tid; 
    }
    
    num_procs =myproc.nprocs;
    
    /*********** PVM- WRITE ********************/
    if (myproc.node!=0){
      pvm_initsend(PvmDataType);
      pvm_pkdouble(dvalue,1,1);
      pvm_send(master_tid,msg_tag);
    }
    /*********** PVM - READ ********************/
    else{  
      for (i=1;i<num_procs;i++){
	pvm_recv(-1,msg_tag);
	pvm_upkdouble(&tmp,1,1);
	*dvalue += tmp;
	
#if DPME_DEBUG 
	printf("merge_d2: Msg=%d rcvd %f total %f \n",msg_tag,tmp,*dvalue);
#endif
	
      }
    }
  } /* end 1st else */
  return 0;
} /* merge_d2*/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/* Global reduction of a double vector, mainly used in summing the virial
 * targeted for configs of non-power of 2  basically just let the master do the sum 
 */

int merge_v2(double *dvalue, MYPROC myproc, int *tidarray, int vec_length,int flag)
{
  int  master,i,j;
  double *tmp;
  int master_tid,num_procs;
  int msg_tag;
  
  /* this will allow for two outstanding merges to be recieved properly ;
   * otherwise, if u have 2 consecutive merge_d, mesgs may get mixed up,
   * eg, in adj_dir and self 
   */

  switch (flag){
  case 1: /* special flage for adjusted dir-virial only */
    msg_tag= MSG_500; 
    break;
  case 2: /* special flage for adj-virial only */
    msg_tag= MSG_510;
    break;
 case 3: /* special flage for recip-virial only */
    msg_tag= MSG_520;
    break;
  default: /* any tag >3 */
    msg_tag = MSG_500;
    break;
  }
 
  master = 0; /* ie the master */
  master_tid = tidarray[master];
  if (master_tid <0) { /* it will <0 due to initiaize */
    master_tid = pvm_gettid(GROUP,master);
    tidarray[master]=  master_tid; 
  }
 
  num_procs =myproc.nprocs;

  /* alloc memory for tmp array */
  if (!((tmp) = (double *) malloc( vec_length * sizeof(double))) )
    fprintf(stderr,"Error in allocating space for tmp array !!!\n"); 
  
  /*********** PVM- WRITE ********************/
  if (myproc.node!=0){
    pvm_initsend(PvmDataType);
    pvm_pkdouble(dvalue,vec_length,1);
    pvm_send(master_tid,msg_tag);
  }
  /*********** PVM - READ ********************/
  else{ /* Pe=0  */
  
    if ( flag==1) { /* dir vir */
      for (i=1;i<num_procs;i++){
	pvm_recv(-1,msg_tag);
	pvm_upkdouble(tmp,vec_length,1);
	
	for (j=0;j<vec_length;j++) 
	  dvalue[j] += tmp[j];
      }
    } /* if  flag ==1 */

    else { /* ie flag ==2,3 */
      pvm_recv(-1,msg_tag);
      pvm_upkdouble(dvalue,vec_length,1);
    } /* else */

  } /* 1st else */
  
  free(tmp);

  return 0;
} /* merge_v2*/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/* Workers pack their direct-sum forces to send to Master */
void dpme_snd_dir(AtomInfo atom_info, int *list,PeInfo pe_info,
		    Pme2Particle *Myparticles, PmeVector *mydirectF,
		    double *mytime)
{
  
  int i,ii,rcvr_tid;
  int atompnt,nlocal;
  MYPROC myproc;
  int *tidarray;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
  /*-----------------------------------*/
  atompnt = atom_info.atompnt;
  nlocal=  (atom_info.nlocal);
  tidarray= pe_info.inst_node;
  myproc=  pe_info.myproc;
/*-----------------------------------*/
  pvm_initsend(PvmDataType);
  pvm_pkint(&myproc.node,1,1);
  pvm_pkint(&nlocal, 1,1);
  i = atompnt;
  for (ii = 1; ii <=  nlocal ; ++ii) {
    pvm_pkdouble(&Myparticles[i].id,1,1);
    i=list[i];
  }  
  i= atompnt;
  for (ii = 1; ii <=  nlocal ; ++ii) {
    pvm_pkdouble(&(mydirectF[i].x),3,1);
    i=list[i];
  }
   
  rcvr_tid=tidarray[myproc.nprocs];
  if (rcvr_tid <0) {
    rcvr_tid=pvm_gettid(GROUP,myproc.nprocs);
    tidarray[myproc.nprocs] = rcvr_tid;
  } 
  
  /* printf(" $$$$PE=%d will be sending dir_force to PE=%d/tid=%d of %d atoms \n",
   * myproc.node,myproc.nprocs,tidarray[myproc.nprocs],nlocal);
   */

  pvm_send(rcvr_tid,MSG_600);
#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif
} /* end worker_snd_dir() */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* Master rcv direct-sum force data */
void dpme_rcv_dir(PePnt **pe_atompnt, double **pe_atoms,  PmeVector **directF,
		    MYPROC myproc, int numatoms, double *mytime)
{
  int i,ii;
  int my_atomid,atompnt,local_atoms,pe_id;
  static int firsttime =0; 
  static PmeVector *mydirectF;
  static PePnt *mype_atompnt;
  static double *mype_atoms;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif


  if (!firsttime) {
    if (!(mydirectF = (PmeVector *) malloc((numatoms+1) * sizeof(PmeVector))) )
      fprintf(stderr,"Error in allocating space for Master  directF  Data !!!\n");
    if (!(mype_atompnt = (PePnt *) malloc((myproc.nprocs) * sizeof(PePnt))) )
      fprintf(stderr,"Error in allocating space for Master pe_atompnt Data !!!\n");
    if (!(mype_atoms = (double *) malloc((numatoms+1) * sizeof(double))) )
      fprintf(stderr,"Error in allocating space for Master pe_atoms Data !!!\n");
    firsttime++;
  }

  (*directF)= mydirectF;
  (*pe_atompnt)= mype_atompnt;
  (*pe_atoms)= mype_atoms;


  atompnt=0;
  for (i=1;i<=myproc.nprocs;i++){ /* for (i=1;i<=myproc.nprocs;i++){ */
    pvm_recv(-1,MSG_600);
    pvm_upkint(&pe_id,1,1); /* proc no. */
    pvm_upkint(&local_atoms,1,1); /* how many atoms I am getting */

    mype_atompnt[pe_id].pnt = atompnt; /* set PE to point to start of its atom-list */
    mype_atompnt[pe_id].tot = local_atoms; /* set PE to point to total atoms it has */

    pvm_upkdouble(&mype_atoms[atompnt],local_atoms,1); /* Atom_id's only */
#if DPME_DEBUG
    printf("$$$ Master rcvd dirfrc=%d atoms from PE=%d \n",local_atoms,pe_id);
#endif

    for (ii=1;ii<=local_atoms; ii++){
      my_atomid=(int)mype_atoms[atompnt];
      pvm_upkdouble(&(mydirectF[my_atomid].x),3,1);
      atompnt ++; /* exiting the loop atompnt =+ local_atoms */
    }
  } /* for loop */

#if TIMEME
  gettimeofday(&(time2),&tzp);
  *mytime=swatch(time1,time2);
#endif

} /* master_rcv_dir */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* master sends pack to workers updated particles positions */
void dpme_snd_xyz(PePnt *pe_atompnt, double *pe_atoms, Pme2Particle *Myparticles,
		 PeInfo pe_info, double *mytime) 
{
  int i,ii,atompnt,rcvr_tid;
  MYPROC myproc;
  int *tidarray;

#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
/*----------------------------------------*/
  myproc= pe_info.myproc;
  tidarray= pe_info.inst_node;
/*----------------------------------------*/
  for (i=0;i<myproc.nprocs;i++){
    pvm_initsend(PvmDataType);
    atompnt=pe_atompnt[i].pnt; /* pointer to start of PE's list */
    
    for (ii=1; ii<=pe_atompnt[i].tot; ii++) {
      /* remeber that the master has its Myparticles list indexed
       * by the global atom_id #, so are the dir and rcp forces 
       */
      pvm_pkdouble(&(Myparticles[(int)pe_atoms[atompnt]].x),3,1);
      atompnt++; /* next atom in pe_atoms[] */
    }
#if 0
    rcvr_tid=tidarray[i];
    if (rcvr_tid <0) {
      rcvr_tid = pvm_gettid(GROUP,i);
      tidarray[i]= rcvr_tid;
    } 
#else
    rcvr_tid = pvm_gettid(GROUP,i);
#endif
    
    /* printf("$$$ Master sending XYZ of %d atoms to PE%d(tid=%d) \n",
     *  pe_atompnt[i].tot,i,rcvr_tid);
     */

    pvm_send(rcvr_tid,MSG_610);
  } /* for all Workers */
#if TIMEME
  gettimeofday(&(time2),&tzp);
  *mytime=swatch(time1,time2);
#endif
} /*end master_snd_xyz */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* workers rcv new (updated) xyz coordinates (from master) in the order they 
 * sent their directF list to the master.
 */
void dpme_rcv_xyz(AtomInfo atom_info, int *list, Pme2Particle *Myparticles,
		    double *mytime) 
{
  int i,ii;
  int atompnt,nlocal;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
  /*------------------------------------------------------------------*/
  atompnt = (atom_info.atompnt);
  nlocal = atom_info.nlocal;
  
  pvm_recv(-1,MSG_610);
  i= atompnt;
  for (ii = 1; ii <=  nlocal ; ++ii) {
    pvm_upkdouble(&(Myparticles[i].x),3,1); 
    i=list[i];
  }
#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif
} /* end worker_rcv_xyz()  */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* call this function anytime the application wants to get the energies 
* of the ewald sum, other wise there is no need to do call it 
*/
void collect_ene_results(double *dir_ene, double *adj_dir, double *self_ene,
			 double *rcp_ene, double *adj_rcp, PeInfo pe_info)
{
  MYPROC myproc;
  int *tidarray;
  int i,master_tid;
  double tmp_ene;
/*----------------------------------------*/
  myproc= pe_info.myproc;
  tidarray= pe_info.inst_node;
  master_tid=pvm_gettid(GROUP,0);
/*----------------------------------------*/
#if VERBOSE
  printf("calling collect_ene_results()....\n");
#endif
  if (myproc.node==myproc.nprocs){ /* ie master PE doing recip */
      pvm_initsend(PvmDataType);
      pvm_pkdouble(rcp_ene,1,1);
      pvm_pkdouble(adj_rcp,1,1);
      pvm_send(master_tid,MSG_700);
  }
  else if  (myproc.node==0) { /* ie this PE was the first to run */
    for (i=1;i<myproc.nprocs;i++){
      pvm_recv(-1,MSG_710);
      pvm_upkdouble(&tmp_ene,1,1);
      *dir_ene += tmp_ene;
      pvm_upkdouble(&tmp_ene,1,1);
      *adj_dir += tmp_ene;
      pvm_upkdouble(&tmp_ene,1,1);
      *self_ene += tmp_ene;
    }
      pvm_recv(-1,MSG_700);
      pvm_upkdouble(rcp_ene,1,1);
      pvm_upkdouble(adj_rcp,1,1);
  }
  else { /* any other PE */
      pvm_initsend(PvmDataType);
      pvm_pkdouble(dir_ene,1,1);
      pvm_pkdouble(adj_dir,1,1);
      pvm_pkdouble(self_ene,1,1);
      pvm_send(master_tid,MSG_710);
  } 

} /* collect_ene_results */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* collect all virial from all processors */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
void collect_virials(double *adj_vir, double *dir_vir, double *rec_vir,
		     PeInfo pe_info, double *mytime)
{
  int i;
  MYPROC myproc;
  int *tidarray;

#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
/*------------------------------*/
  myproc= pe_info.myproc;
  tidarray= pe_info.inst_node;
/*-----------------------------*/
#if VERBOSE
  printf(" collecting virials ....\n"); 
#endif
 if (myproc.node==myproc.nprocs) {         /* ie master PE doing recip */
   merge_v2(adj_vir,myproc,tidarray,6,2);
   merge_v2(rec_vir,myproc,tidarray,6,3);
 }
 else if  (myproc.node==0) {              /* ie this PE was the first to run */
   merge_v2(dir_vir,myproc,tidarray,6,1);
   merge_v2(adj_vir,myproc,tidarray,6,2);
   merge_v2(rec_vir,myproc,tidarray,6,3);
 }
 else {                                   /* any other PE */
   merge_v2(dir_vir,myproc,tidarray,6,1);
 }

  if (myproc.node==0)
    for (i=0;i<6;i++){
      printf("  dir_vir[%d]=%f  ",i,dir_vir[i]);
      printf("rcp_vir[%d]=%f  ",i,rec_vir[i]);
      printf("adj_vir[%d]=%f  \n",i,adj_vir[i]);
    }

#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif
} 
/* collect_virials */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
