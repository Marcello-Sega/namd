/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
/* $Id: dpme2_paralibs.c,v 1.1 1997/03/14 15:07:01 nealk Exp $
 */

/***********************************************************
* This file contain routines necessary to starting off 
* the pvm process and group. A pvm group of workers (slaves)
* is named (direct_solvers) and it performs the direct_sum
* in parallel.
*/

#include "dpme2.h"
#include "pvm3.h"

/**********************************************************
* This routine is called by all processes and have sections
* that are executed only by the master process (group id =0)
* this section can be modified for spawning on the cray T3D
*/

void 
dpme_register(char *spawnme, int *mygroupid, int nproc, 
	      int *tidarray, char **args, int *ncube)
{
  
  int i,mytid, cnode, nslaves;
  char ProgName[80];
  
  nslaves = nproc -1 ;
  /* enroll in PVM */
  
  if ( (mytid= pvm_mytid()) <0) 
    { printf (" * The instant is %d\n",mytid);
      fprintf (stderr,"error: enrolling PVM, exiting \n");
      exit(1);
    }; 
  
  /* Join the group; all the group of direct solvers */
  *mygroupid= pvm_joingroup(GROUP);  
  
#if VERBOSE
  printf (" * NODE %d :  instant is %d \n", *mygroupid, mytid);
#endif
                
  /* spawn is only valid for the  N.O.W. cluster, take out for T3D */

  if ( (*mygroupid==0) && (nslaves!=0) ){
#if 1
  
#if VERBOSE
    sprintf(ProgName, "%s%s", WORKINGDIR , spawnme);
    printf (" %s \n", ProgName);
#endif
    cnode=pvm_spawn(ProgName, args, 0, "", nslaves,tidarray);
    if ((cnode!= nslaves)&&(nslaves!=0)) {
      fprintf(stderr,"Error: spawned %d of %d processes\n",cnode,nslaves);
      error_handler("Check your WORKINGDIR def");
    }
#endif
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
  for (i=0;i<nproc;i++)
    tidarray[i]=-1;

  /* all should sync now */
  if(  (pvm_barrier(GROUP,nproc)) <0) {
    printf(" Error: at the   barrier \n");
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
    if (myproc->node!=0){
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
  printf("%%%%%%%%%%%%%%%%%%%%%Node %d rcvd%%%%%%%%%%%%%%%%%% \n",node);
  --rbuf; 
  for (i=0;i<numdouble/5.;i++)
    printf("(%f, %f, %f, %f, %f)  ",*(++rbuf), *(++rbuf), *(++rbuf), *(++rbuf),
	   *(++rbuf));
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

int merge_v2(double *dvalue, MYPROC myproc, int *tidarray, int vec_length, int msg_tag)
{
  int  master,i,j;
  double *tmp;
  int master_tid,num_procs;

  
  /* this will allow for two outstanding merges to be recieved properly ;
   * otherwise, if u have 2 consecutive merge_d, mesgs may get mixed up,
   * eg, in adj_dir and self 
   */

  switch (msg_tag){
  case 1: 
    msg_tag= MSG_500;
    break;
  case 2:
    msg_tag= MSG_510;
    break;
  default: /* any tag >2 */
    msg_tag = MSG_520;
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
  else{  
    for (i=1;i<num_procs;i++){
      pvm_recv(-1,msg_tag);
      pvm_upkdouble(tmp,vec_length,1);
      
      for (j=0;j<vec_length;j++) {
	dvalue[j] += tmp[j];
#if DPME_DEBUG 
      printf("merge_v2: Msg=%d rcvd %f total %f \n",msg_tag,tmp[j],dvalue[j]);
#endif
      }
    }
  }
  
  free(tmp);

  return 0;
} /* merge_v2*/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

