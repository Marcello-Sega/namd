/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
/* $Id: dpme2_utility.c,v 1.1 1997/03/14 15:07:05 nealk Exp $
 */

/*******************************************************************************
 * A. (Nour) Toukmaji - Duke University; ECE Dept. 1996
 * Acknowledgment: the subroutines in this file have been incorporated / modified
 * from original PME source code by Tom Darden NIEHS, N.C.
 **********************************************************
 * This file has utility subroutines that are used for by the driver application
 ********************************************************************************/

#include "dpme2.h"
#include "pvm3.h"
#define FIXIT 1 /* this was added to fix the swap boundaries when NPE=3 */
/* ------------------------------------------------------------------ */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* findes the ewald coeffecient */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
 int find_ewaldcof(double *cutoff, double *dtol, 
	double *ewaldcof)
{
 
    /* Local variables */
    double myerfc, term;
    int i, n;
    double x, y, xhi, xlo;
    
/* first get direct sum tolerance. How big must ewaldcof be to get */
/* terms outside the cutoff below tol */
  
    x = .5;
    i = 0;

    do {
    x *= 2.;
    ++i;
    y = x * *cutoff;
    myerfc=erfc(y);
    term = myerfc / *cutoff;
  }
    while(term >= *dtol);
    
/* binary search tolerance is 2 to the -60th */
    n = i + 60;
    xlo = 0.;
    xhi = x;
    for (i = 1; i <= n; ++i) {
	x = (xlo + xhi) / 2;
	y = x * *cutoff;
	myerfc=erfc(y);
	term = myerfc / *cutoff;
	if (term >= *dtol) {
	    xlo = x;
	} else {
	  xhi = x;
	}
      }
    *ewaldcof = x;
    return 0;
  } /* find_ewaldcof */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* ---------------------------------------------------------- */
/* This function sets up general parameters that has to do with
 * with geometry of the system 
 */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int setup_general(double *cutoff, double *dtol, double *ewaldcof, int *numatoms, 
		  PmeBVector *box, PmeVector *prd, PmeVector *prd2, int tsteps, 
		  Pme2Particle **ParticlePtr, PmeVector *mymc2)
{
  double factor, cgh, cgo;
  char DataFile[80], s1[81];
  FILE  *infile; 
  int i; 
  static Pme2Particle *myparticles;
  
  /* For this value of cutoff and direct space tolerance 
   * find the ewald coefficient ewaldcof.                
   */
  
  find_ewaldcof(cutoff, dtol, ewaldcof);
  

  /* water charges O and H */
  factor = sqrt(332.17752);
  cgh = factor * .417;
  cgo = cgh * -2.;

  /* open particle data  file */
  sprintf(DataFile,"%s%s",DATADIR,"small2.pdb");
  if ( (infile = fopen(DataFile,"r") ) == NULL) {
    fprintf(stderr,"Unable to open the input file: !!! \n");
    exit(1); 
  }
  
  fscanf(infile,"%s %s %s",s1,s1,s1); /* will read the text */
  fscanf(infile,"%lf%lf%lf",&(box->x), &(box->y), &(box->z)); /* ortho-box length */
  fscanf(infile,"%d",numatoms);
  fscanf(infile,"%lf",&(box->origin));
#if VERBOSE
  printf("Particle Sytem size= %d ..... Ewald coefficient = %f\n",
	  *numatoms, *ewaldcof); 
  printf("Box_dim[%f,%f,%f] .... cutoff = %f ..... tol=%f\n",box->x,box->y,box->z,
	 *cutoff,*dtol);
#endif 

  prd->x = box->x;
  prd->y = box->y;
  prd->z = box->z;

  prd2->x = 0.5 * box->x;
  prd2->y = 0.5 * box->y;
  prd2->z = 0.5 * box->z;

  mymc2->x= 2. *  ( prd->x - *cutoff);
  mymc2->y= 2. *  ( prd->y - *cutoff);
  mymc2->z= 2. *  ( prd->z - *cutoff);
/********************************************/
/* allocate main particle info array (x,y,z,cg), using double pointer */
  if (!((myparticles) = (Pme2Particle *) malloc( (*numatoms+1) * sizeof(Pme2Particle))) )
    fprintf(stderr,"Error in allocating space for Particle Data !!!\n");
  (*ParticlePtr)= myparticles;

  
  for (i =1; i<= *numatoms ;i++){
    fscanf(infile,"%lf %lf %lf",&(myparticles[i].x),&(myparticles[i].y),
	   &(myparticles[i].z));

    if( (i-4) % 3 == 0) 
      myparticles[i].cg= cgo;
    else 
      myparticles[i].cg= cgh;
    
   /* printf("atom[%i]=(%lf,%lf,%lf,%lf) \n",i,myparticles[i].x,myparticles[i].y,
	   myparticles[i].z,myparticles[i].cg); */
  } /* for i */
 fclose(infile);
 
  /*  you may need to call barrier here! */ 
/********************************************/ 
  return 0;
} /* setup_general */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* This routine maps the simulation box to the number of processors
 * it assumes that the processors are configured in cube of dimension
 * ncube, and hence it defines geometrical boundaries for each 
 * processor ( basically the borders )
 */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int setup_parallel(int numatoms, int *npdim, double rs, PmeVector prd, PmeVector prd2, 
		   int *need, double *border, MYPROC *myproc, int *mpart, int *nswap, 
		   PmeVector *nbin, PmeVector *binsize, PmeVector *mbin, PmeVector *mbinlo,
		   double *boundlo, double *boundhi, int *spart, int *rpart, int igrid, 
		   int *tidarray)
{
  double smallnum = 1.e-6;
  int ineigh, ibin; /* flags decide alt methods to use */
  double myprd[3]; /* store prd structure for looping */
  int k,kk,i_1,iflag;
  int me[3];
  double blo, bhi,btmp,bmid;
  int ix,iy,iz,nbox;

#if (VERBOSE)
  printf(" calling setup_parallel()... \n");
#endif


  if (numatoms >1000 ) { 
    ineigh= 1; /* use Link_cell */
#if (VERBOSE)
    printf(" You are using Link_cell Method\n");
    printf(" The UPDATE_TIME=%d\n",UPDATE_TIME);
#endif
  }
  else {
    ineigh =0; /* use Verlet-list */
#if (VERBOSE)
    printf(" You are using Verlet-list Method \n");
    printf(" The UPDATE_TIME=%d\n",UPDATE_TIME);
#endif
  }
  
  ibin =0; /* program will figure out number of bins */
 

 if (igrid==0) {  /* program will figure out processor-grid config */
   npdim[0]= pow(2.0,(myproc->ncube/3));
   npdim[1]= pow(2.0,((myproc->ncube+1)/3));
   npdim[2]= pow(2.0,((myproc->ncube+2)/3));
 }
  if (npdim[0]*npdim[1]*npdim[2]!= (myproc->nprocs)){
    error_handler("Bad grid of processors, try power of 2, or specify the proc-grid");
  }

#if VERBOSE
    printf("npdim=(%d,%d,%d) total procs=%d \n",npdim[0],npdim[1],npdim[2],
	    (myproc->nprocs));
#endif

 /* setting up a 3_D mesh with pbc */
  mesh_3d(npdim, &npdim[1], &npdim[2], &(myproc->node), 
	  me, &me[1], &me[2], mpart, &mpart[1], 
	  &mpart[2], &mpart[3], &mpart[4], &mpart[5]);
 
 

  border[0] = (double) me[0] / npdim[0] * prd.x - prd2.x;
  border[1] = (double) (me[0] + 1) / npdim[0] *  prd.x - prd2.x;
  border[2] = (double) me[1] / npdim[1] *  prd.y - prd2.y;
  border[3] = (double) (me[1] + 1) / npdim[1] * prd.y - prd2.y;
  border[4] = (double) me[2] / npdim[2] * prd.z - prd2.z;
  border[5] = (double) (me[2] + 1) / npdim[2] * prd.z - prd2.z;

  need[0] = rs  / (prd.x / npdim[0]) + 1.;
  need[1] = rs / (prd.y / npdim[1]) + 1.;
  need[2] = rs / (prd.z / npdim[2]) + 1.;

/* don't exchange if only 1 box in a dimension */
    if (npdim[0] == 1) need[0] = 0;
   
    if (npdim[1] == 1) need[1] = 0;
   
    if (npdim[2] == 1) 	need[2] = 0;
   
/*don't exchange more than 1/2 way over 
  (e.g. 3 boxes away when npdim = 5 ) */

    if (need[0]*2 > npdim[0]) --need[0];
   
    if (need[1]*2 > npdim[1]) --need[1];
    
    if (need[2]*2 > npdim[2]) --need[2];

/* setup 4 parameters for each exchange: (spart,rpart,boundlo,boundhi) */
/*  (1,2) nodes to swap with (3,4) slab boundaries 
 * (in correct dimension) of atoms that will be sent*/
/* 1st part of if is sending to the west (south,down) */
/* 2nd part of if is sending to the east (north,up) */
/* nbox = box (in this dimension) who originally owned the atoms */
/*   I will be sending in this swap */
/* kk/2+1 vs npdim/2 test is to make sure a node doesn't get more than 1/2 */
/*   of a box from both directions (else would get some atoms twice) */

  myprd[0] = prd.x;
  myprd[1] = prd.y;
  myprd[2] = prd.z;
  *nswap = 0;
  for (k = 1; k <= 3; ++k) {
    i_1 = (need[k - 1] * 2) - 1;
    for (kk = 0; kk <= i_1; ++kk) {
      ++ (*nswap);
      if (*nswap <= nsmax) {
	if (kk % 2 == 0) {
	  spart[*nswap] = mpart[(k << 1) - 2];
	  rpart[*nswap] = mpart[(k << 1) - 1];
	  nbox = me[k - 1] + kk / 2;
	  if (nbox >= npdim[k - 1]) {
	    nbox -= npdim[k - 1];
	  }
	  blo = (double) nbox / npdim[k - 1] * 
	    myprd[k - 1] - myprd[k - 1] / 2.0;
	  if (kk / 2 + 1 < need[k - 1]) {
	    bhi = (double) (nbox + 1) / npdim[k - 1] * 
	      myprd[k - 1] - myprd[k - 1] / 2.0;
	  } else {
	    bhi = border[(k << 1) - 2] + rs ;
	    if (bhi >= myprd[k - 1] / 2.0) {
	      bhi -= myprd[k - 1];
	    }
	  }
	  if (kk / 2 + 1 == npdim[k - 1] / 2) {
	    btmp = (double) (nbox + 1) / npdim[k - 1] * 
	      myprd[k - 1] - myprd[k - 1] / 2.0;
	    
	    bmid = (blo + btmp) / 2.0;
#if !FIXIT
	    bhi = dmin(bhi,bmid);
#endif 

#if FIXIT
	    if ((npdim[k-1] %2 != 0) && (npdim[k-1] %3 != 0) )
	      bhi= dmax(bhi,bmid);
	    else {
	      bhi = dmin(bhi,bmid);
	    }
	    if (npdim[k-1]==3) bhi=btmp;
#endif   
	  }
	} else {
	  spart[*nswap ] = mpart[(k << 1) - 1];
	  rpart[*nswap] = mpart[(k << 1) - 2];
	  nbox = me[k - 1] - kk / 2;
	  if (nbox < 0) {
	    nbox += npdim[k - 1];
	  }
	  bhi = (double) (nbox + 1) / npdim[k - 1] * 
	    myprd[k - 1] - myprd[k - 1] / 2.0;
	  if (kk / 2 + 1 < need[k - 1]) {
	    blo = (double) nbox / npdim[k - 1] * 
	      myprd[k - 1] - myprd[k - 1] / 2.0;
	  } else {
	    blo = border[(k << 1) - 1] - rs ;
	    if (blo < -(double)myprd[k - 1] / 2.0) {
	      blo += myprd[k - 1];
	    }
	  }
	  if (kk / 2 + 1 == npdim[k - 1] / 2) {
	    btmp = (double) nbox / npdim[k - 1] * 
	      myprd[k - 1] - myprd[k - 1] / 2.0;
	    bmid = (btmp + bhi) / 2.0;
#if !FIXIT
	    blo = dmax(blo,bmid);
#endif 

#if FIXIT
	   /* this is the new fix for swap boundaries */
	    if ((npdim[k-1] %2 != 0) && (npdim[k-1] %3 != 0) )
	      blo= dmin(blo,bmid);
	    else {
	      blo = dmax(blo,bmid);
	    }
	    if (npdim[k-1]==3) blo=btmp;
#endif
	  }
	}
	boundlo[*nswap - 1] = blo;
	boundhi[*nswap - 1] = bhi;
      }
    }
  }
 

/* setup neighbor binning parameters in box owned by each processor */
/*  addding small -> bins slightly larger */
/*  prevents round-off error (bin # ix = nbinx) when atoms are binned */
  
  if (ineigh == 1) {
    /* For LC method */
    if (ibin == 0) {
      nbin->x = (int) (prd.x / (rs));
      nbin->y = (int) (prd.y / (rs));
      nbin->z = (int) (prd.z / (rs));
    }
    binsize->x = (prd.x + prd.x * smallnum) /nbin->x;
    binsize->y = (prd.y + prd.y * smallnum) /nbin->y;
    binsize->z = (prd.z + prd.z * smallnum) /nbin->z;
    mbinlo->x = (int) ((border[0] + prd2.x) / binsize->x) - 1;
    mbinlo->y = (int) ((border[2] + prd2.y) / binsize->y) - 1;
    mbinlo->z = (int) ((border[4] + prd2.z) / binsize->z) - 1;
    ix = (int) ((border[1] + prd2.x) / binsize->x) + 1;
    iy = (int) ((border[3] + prd2.y) / binsize->y) + 1;
    iz = (int) ((border[5] + prd2.z) / binsize->z) + 1;
    if (ix > nbin->x) {
      ix = nbin->x;
    }
    if (iy > nbin->y) {
      iy = nbin->y;
    }
    if (iz > nbin->z) {
      iz = nbin->z;
    }
    mbin->x = (int) (ix - mbinlo->x + 1);
    mbin->y = (int) (iy - mbinlo->y + 1);
    mbin->z = (int) (iz - mbinlo->z + 1);
    mbin->x = min(mbin->x,nbin->x);
    mbin->y = min(mbin->y,nbin->y);
    mbin->z = min(mbin->z,nbin->z);
    if (mbinlo->x < 0) {
      mbinlo->x = (int)nbin->x - 1;
    }
    if (mbinlo->y < 0) {
      mbinlo->y = (int)nbin->y - 1;
    }
    if (mbinlo->z < 0) {
      mbinlo->z = (int)nbin->z - 1;
    }
  }
/**************DPME_DEBUG***************/
#if DPME_DEBUG2
  printf("setup_parallel Proc-%i:\n borders %f %f %f %f %f %f \n",myproc->node,
	 border[0],border[1],border[2],
	 border[3],border[4],border[5]);
 for (k=1;k<=*nswap;k++)
   printf("K=%i, spart=%i, rpart=%i \n",k,spart[k],rpart[k]);
fflush(stdout);
#endif
  /************ handle parameter paroblems *****************/
  /* ierrorflag = 0; not used cuz error_handle exits anyway*/
  if (rs  >= prd2.x || rs  >= prd2.y || rs  >=  prd2.z) {
    error_handler("Outer cutoff >= 1/2 box size");
  }
  if (ineigh == 1 && (nbin->x <= 2 || nbin->y <= 2 || nbin->z <= 2)) {
    error_handler("Two or less bins in a dimension");
  }
 
  if (*nswap > nsmax) {
    error_handler("Swap array too small - boost nsmax");
  }
  iflag = 0;
  if (ineigh == 1 && mbin->x * mbin->y * mbin->z > nbmax) {
    iflag = -1;
  }
  
  merge_i2(&iflag,myproc,tidarray);
  if (iflag < 0) {
    error_handler("Too many local bins - boost nbmax");
  }

  return(ineigh);
} /* setup_parallel */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* Each processor read  positions  of all atoms               */
/*  but will only store ones in its own box                   */
/* NO VELOCITIES  in this code, can add independently         */
/* Reads in the input file small.pdb and each node keeps atoms* 
 * that are within its borders 
 */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int setup_atom(int **mylist, Pme2Particle **ParticlePtr, PmeBVector box, int numatoms, 
	       double *border, PmeVector prd2, MYPROC *myproc, int *nlocal, 
	       int *freepnt, int *atompnt, int *tidarray)
{
  int  *list;
  Pme2Particle *TempoParticlePtr;
  Pme2Particle *MyParticlePtr;

  /* Local variables */
  int itmp, max_list_sz;
  double xtmp, ytmp, ztmp;
  int i,n, iflag;
  char mystr[100];

  /* note that list was maxed with npmax however here is a good
     estimate for it--- not used yet! */
  max_list_sz= (int)(1.5* numatoms/(myproc->nprocs)) + 1;
  
  /* allocate List which indexes particles*/
  if (!(list = (int  *) malloc((npmax+1) * sizeof(int))) )
    fprintf(stderr,"Error in allocating space for mylist  Data !!!\n");
  (*mylist)= list; /* return the global index */

  for (i = 1; i <= npmax; ++i) {
    list[i] = i+1 ; 
  }
  
  *freepnt = 1;
  list[npmax] = 0;
  *nlocal = 0;
  *atompnt = npmax+1 ;

#if (VERBOSE) 
  printf("setupatom:borders %f %f %f %f %f %f \n",border[0],border[1],border[2],
	 border[3],border[4],border[5]);
  printf(" NOTE: ALL COORDS ARE SCALED from  -L/2 -> L/2\n");
#endif
  
  TempoParticlePtr = (*ParticlePtr);
 
  for (i =1; i<= numatoms;i++){
    
    xtmp= TempoParticlePtr[i].x - prd2.x - box.origin;
    ytmp= TempoParticlePtr[i].y - prd2.y - box.origin;
    ztmp= TempoParticlePtr[i].z - prd2.z - box.origin;

    if (xtmp >= border[0] && xtmp < border[1] && 
	ytmp >= border[2] && ytmp < border[3] && 
	ztmp >= border[4] && ztmp < border[5]) {
      ++ (*nlocal);
      if (*freepnt != 0) {
        itmp = *atompnt;
        *atompnt = *freepnt;
        *freepnt = list[*freepnt];
        list[*atompnt] = itmp;
#if (DPME_DEBUG3)
	printf("list[%d]=%d \n",*atompnt,   list[*atompnt]);
#endif
	TempoParticlePtr[i].x = xtmp;
	TempoParticlePtr[i].y = ytmp;
	TempoParticlePtr[i].z = ztmp;
      } /* if freepnt */
    } /* if xtmp... */
    else {
      TempoParticlePtr[i].cg=0. ; /* reset to exclude later */ 
    }
#if DPME_DEBUG3   
    printf("COORDS %d %f %f %f %f\n",i,TempoParticlePtr[i].x,TempoParticlePtr[i].y,
	   TempoParticlePtr[i].z, TempoParticlePtr[i].cg ); 
#endif     

  } /* for i: numatoms */

/* copy ALL particle array to new array of size nlocal and free */
/****************************************************************/
  /* allocate main particle info array (x,y,z,cg)*/
  if (!(MyParticlePtr = (Pme2Particle *) malloc((namax+1) * sizeof(Pme2Particle))) )
    fprintf(stderr,"Error in allocating space for My Particle Data !!!\n");

  i=1; /* Myparticleptr index */
  for (n =1; n<= numatoms;n++){
    if (TempoParticlePtr[n].cg !=0.) {
      MyParticlePtr[i].x=TempoParticlePtr[n].x;/* copy contents */
      MyParticlePtr[i].y=TempoParticlePtr[n].y;
      MyParticlePtr[i].z=TempoParticlePtr[n].z;
      MyParticlePtr[i].cg=TempoParticlePtr[n].cg;
      MyParticlePtr[i].id= n;
      i++; /* update to next atom */
    }
 } /* for n: numatoms */
   free(ParticlePtr[0]); 
  (*ParticlePtr) =  MyParticlePtr;  /* return new pointer */
  

/* Handle errors (if any )                  */
/********************************************/ 
  iflag = 0;
  if (numatoms > 1000000 )
    pvm_barrier(GROUP,-1);/* synchro all grouup members */ 
  if (*nlocal > npmax) iflag = -1;
  merge_i2(&iflag,myproc,tidarray);
  if (iflag < 0) {
    sprintf(mystr,"Too many atoms/processor(%d) - boost npmax(%d)",*nlocal,npmax);
    error_handler(mystr);
  }  
#if VERBOSE
  printf("Num of particle I own=%i \n",*nlocal);
#endif
  return 0;
} /* setup_atom */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* send out atoms that have left my box, receive ones entering */
/*  my box since last reneighboring  done in all 3 directions  */
/*  done before every reneighboring.                           */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int exchange(int *nlocal, Pme2Particle **AllParticles, double *border, int *list, 
	     int *npdim, int *freepnt, int *atompnt, int *mpart, MYPROC myproc,
	     int *tidarray,double *mytime)
{
  
  int icnt, itmp;
  int i, j, k, iprev, ii;
  double bhi, blo;
  int ndelete;
  double tmp;
  double  *sndbuf, *rcvbuf;
  double *particle;
  double particle_pos;
  static int nexcmax =0; /* max # of atoms leave'g mybox (not total) */
  Pme2Particle *Myparticles;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
  /***************************************************************/  
#if (VERBOSE)
  printf(" calling exchange()... \n");
#endif 

  Myparticles= (*AllParticles); /* pointer setup */

  /* allocate memory for send and recv buffers */
  /* may have to come up with a better estimate for allocc/realloc */
  if (!((sndbuf) = (double *) malloc( nfmax * sizeof(double))) )
    fprintf(stderr,"Error in allocating space for SNDBUF !!!\n"); 

   if (!((rcvbuf) = (double *) malloc( 2*nfmax * sizeof(double))) )
    fprintf(stderr,"Error in allocating space for RCVBUF !!!\n"); 

  for (k = 1; k <= 3; ++k) {
    if (npdim[k - 1] > 1) {
      blo = border[(k << 1) - 2];
      bhi = border[(k << 1) - 1];
      ndelete = 0;
      iprev = 0;
      j = 0;
      i = *atompnt;
      

      /* fill buffer with atoms leaving my box, update local list */
      for (ii = 1; ii <= *nlocal; ++ii) {

	particle= &(Myparticles[i].x) ; /* refs the whole (xyz cg id ) strc */
	particle_pos= *(particle + k - 1);
	
	/* printf("k=%d, particle_pos=%f blo=%f, bhi=%f \n",
	 *    k,particle_pos,blo,bhi);
	 */

	if ( particle_pos < blo || particle_pos >= bhi) { 
	  sndbuf[j] = *particle; /* x */
	  sndbuf[j + 1] = *(++particle); /* y */
	  sndbuf[j + 2] = *(++particle); /* z */
	  sndbuf[j + 3] =  *(++particle); /* cg */
	  sndbuf[j + 4] = *(++particle); /* new double id */

#if DPME_DEBUG 
	printf("XCHNG: proc(%d) deleting atom(%f)\n",myproc.node,Myparticles[i].id);    
	printf("XCHNG: del atom (%f,%f,%f,%f,%f)\n ",   sndbuf[j],
	     sndbuf[j + 1], sndbuf[j + 2], sndbuf[j + 3],sndbuf[j + 4] ); 
#endif	   

	  j += 5;
	  ++ndelete;
	  itmp = list[i];
	  if (iprev == 0) {
	    *atompnt = itmp;
	  } else {
	    list[iprev] = itmp;
	  }
	  itmp = list[i];
	  list[i] = *freepnt;
	  *freepnt = i;
	  i = itmp;
	} else {
	  iprev = i;
	  i = list[i];
	}
      } /* for ii */
       
      *nlocal -= ndelete;
      /* set the max number of atoms ever leaving my box any tstep */
      nexcmax = max(nexcmax,j/5); /* may be useful */
      
      if (j > nfmax) {
	printf("Sending too many exchange atoms: %d, %d, %d\n ",
	       myproc.node,mpart[(k << 1) - 2],k);
	error_handler("Try increasing nfmax size");
      }
      
      /* send in both directions (if neighboring nodes different) */

      icnt=0; /* initalize */	 
      itmp=0; 
      swap(myproc.node, sndbuf, j, mpart[(k << 1) - 2], rcvbuf, &icnt, 
	   mpart[(k << 1) - 1],tidarray);

      if (npdim[k - 1] > 2) {
	swap(myproc.node, sndbuf, j, mpart[(k << 1) - 1], &rcvbuf[icnt],
	     &itmp, mpart[(k << 1) - 2],tidarray);
	icnt += itmp;
      }

      /* check incoming atoms if they are in my box (could be in node's */
      /*  box on other side of the sender) */

      for (j = 1; j <= icnt; j += 5) {
	/* printf("UNPACK RCV of sz %d \n",icnt); */
	tmp = rcvbuf[j + k - 2]; /* which coordinate to consider */
	
	if (tmp >= blo && tmp < bhi){ 
	  ++(*nlocal);
	  if (*freepnt != 0) {
	    itmp = *atompnt;
	    *atompnt = *freepnt;
	    *freepnt = list[*freepnt];
	    list[*atompnt] = itmp;
	    Myparticles[*atompnt].x=  rcvbuf[j - 1];
	    Myparticles[*atompnt].y=  rcvbuf[j];
	    Myparticles[*atompnt].z=  rcvbuf[j +  1];
	    Myparticles[*atompnt].cg= rcvbuf[j +  2];
	    Myparticles[*atompnt].id= rcvbuf[j +  3];
#if DPME_DEBUG
	    printf("XCHNG: proc(%d) adding atom(%f)\n",myproc.node, 
		   Myparticles[*atompnt].id);

	    printf("XCHNG add (%f,%f,%f,%f,%f)  ", Myparticles[*atompnt].x,
		   Myparticles[*atompnt].y, Myparticles[*atompnt].z,
		   Myparticles[*atompnt].cg,Myparticles[*atompnt].id);
#endif
	  }
	}
      } /* for j */

      if (*nlocal > npmax) {
	printf("Node %d, has now  %d atoms \n", myproc.node, *nlocal);
	error_handler("XCHNG:Received too many exchange atoms, may be boost npmax");
      }
    } /* if npdim */
  } /* for k */

  free(sndbuf);
  free(rcvbuf);

#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif
  return 0;
} /* exchange */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* swap slabs of atoms with other processors in all 3 directions */
/*  done every timestep */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int communicate(int nswap, int *nslist, int *slist, MYPROC myproc, Pme2Particle *Myparticles, 
		int *spart, int *rpart,int *tidarray,double *mytime)
{
  int i_2,i,ii,j,jj,k,icnt,ipnt;
  double  *sndbuf,*rcvbuf;

#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif

  if (!((sndbuf) = (double *) malloc( nfmax * sizeof(double))) )
    fprintf(stderr,"Error in allocating space for SNDBUF !!!\n"); 
  if (!((rcvbuf) = (double *) malloc( 2*nfmax * sizeof(double))) )
    fprintf(stderr,"Error in allocating space for RCVBUF !!!\n");
  
  ipnt = npmax; /* a global constant */
  for (k = 1; k <= nswap; ++k) {
    j = 0;
    i_2 = nslist[k+1] - 1;
    for (ii = nslist[k]; ii <= i_2; ++ii) {
      i = slist[ii];
      sndbuf[j]= Myparticles[i].x;
      sndbuf[j+1]= Myparticles[i].y;
      sndbuf[j+2]= Myparticles[i].z;
      sndbuf[j+3]= Myparticles[i].cg;
      sndbuf[j+4]= Myparticles[i].id;
      
      j += 5;
    }
    icnt = 0; /* total num of doubles recvd */
    swap(myproc.node, sndbuf, j, spart[k], rcvbuf, &icnt, rpart[k], tidarray);
    /* unpack buffer recvd, see notes at borders()  */
    jj=ipnt+1;
    for (j = 1; j <= icnt; j +=5,jj++) {
      Myparticles[jj].x=rcvbuf[j-1];
      Myparticles[jj].y= rcvbuf[j];
      Myparticles[jj].z=rcvbuf[j +  1];
      Myparticles[jj].cg=rcvbuf[j +  2];
      Myparticles[jj].id=rcvbuf[j +  3];
      /* printf("*rcvd[%d] (%f,%f,%f,%f,%i)  ",jj,Myparticles[jj].x, 
	 Myparticles[jj].y ,Myparticles[jj].z, Myparticles[jj].cg,
	 Myparticles[jj].id); */
    }
    ipnt += icnt /5;
  } /* for k */
  free(sndbuf);
  free(rcvbuf);

#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif

  return 0;
} /* communicate */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* make lists of nearby atoms to send to nbr nodes at every timestep*/
/*  one list made for every swap that will be made */
/*  as list is made, actually do swaps */
/*  this does equivalent of a communicate (so don't need to explicitly */
/*   call communicate routine on reneighboring timestep) */
/*  3 atom positions, charge and atom-id  are what is swapped */
/*  this routine called before every reneighboring */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int borders(int *nswap, double *boundlo, double *boundhi, int *atompnt, int *nlocal, 
	    int *nother, Pme2Particle **AllParticles, int *list, MYPROC myproc, int *need, 
	    int *nslist, int *slist, int *spart, int *rpart,int *tidarray,double *mytime)
{
 int npnt,icnt;
 int i_1,i_2,ii,i,j,jj,k,kk;
 static int nswpmax=0; /* max # of atoms ever sent in one swap */   
 Pme2Particle *Myparticles;
 double  *particle,particle_pos,blo,bhi;
 double  *sndbuf,*rcvbuf;
 char str[120];
#if TIMEME
 struct timeval time1,time2;
 struct timezone tzp;
 gettimeofday(&(time1),&tzp);
#endif

#if (VERBOSE)
  printf(" calling borders()... \n");
#endif 

 *nswap = 0;
 *nother = 0;
 npnt = 1;
 Myparticles= (*AllParticles); /* pointer setup */
 
 /* allocate memory for send and recv buffers */
 /* may have to come up with a better estimate for allocc/realloc */

 if (!((sndbuf) = (double *) malloc( nfmax * sizeof(double))) )
   fprintf(stderr,"Error in allocating space for SNDBUF !!!\n"); 
 if (!((rcvbuf) = (double *) malloc( 2*nfmax * sizeof(double))) )
   fprintf(stderr,"Error in allocating space for RCVBUF !!!\n");

/* buffer up all atoms I own (plus those previously received) 
	that are inside slab boundaries. 
	also store pointers to those atoms in slist for communicate 
	routine to use in future timesteps
*/

 for (k = 1; k <= 3; ++k) {
   i_1 = (need[k - 1] << 1) - 1;
   for (kk = 0; kk <= i_1; ++kk) {
     ++(*nswap);
     nslist[*nswap] = npnt;
     blo = boundlo[*nswap - 1];
     bhi = boundhi[*nswap - 1];
#if DPME_DEBUG
     printf(" \n NODE %d with blo=%f and bhi=%f \n",myproc.node,blo,bhi);
#endif
     j = 0;
     i = *atompnt;
     i_2 = *nlocal + *nother;
     for (ii = 1; ii <= i_2; ++ii) {
       particle= &(Myparticles[i].x); /* refs the whole (xyz cg id ) strc */   
       particle_pos= *(particle + k - 1);
       /* particle within my borders */
       if ( particle_pos >= blo && particle_pos < bhi) { 
	 if (npnt <= nemax)  slist[npnt] = i;
	 ++npnt;
	 if (j <= nfmax) {
	   sndbuf[j] = *particle; /* x */
	   sndbuf[j + 1] =  *(++particle); /* y */
	   sndbuf[j + 2] =  *(++particle); /* z  */
	   sndbuf[j + 3] =  *(++particle); /* cg */
	   sndbuf[j + 4] =  *(++particle); /* double id */
	   /* sndbuf[j + 4] = *((int *)(++particle));  int id */
	   /* printf("Packing atom%d.k%d(%f)=%f  ",i,k,sndbuf[j + 4], particle_pos);*/
	 }
	 j += 5;
       }
       if (ii <= *nlocal) {
	 i = list[i];
       } else {
	 ++i;
       }
     } /* for ii */
     /* Computing MAX , may be useful ?*/
     nswpmax = max(nswpmax,j/5);
     
     if (npnt > nemax) {
#if VERBOSE
       printf(" Node:%i Num to swap %i > nemax (%i) \n",
	      myproc.node,npnt,nemax);
#endif
       sprintf(str,"borders():Too many atoms(%d) in border list, boost nemax(%d)?",
	       npnt,nemax);
       error_handler(str);
     }
     if (j > nfmax){   
       sprintf(str,"borders():Sending too many border atoms(%d/5), boost nfmax(%d)?",
	       j,nfmax);
       error_handler(str);
     }
     /* swap atoms, put incoming ones at end of my position array */
     icnt = 0; /* set to # of atoms*5 rcvd in swap */
     
     swap(myproc.node, sndbuf,j, spart[*nswap], rcvbuf,
	  &icnt, rpart[*nswap],tidarray);

     jj= npmax+ *nother + 1;

     /* unpack buffer recvd: originally u could pass &Myparticles[1+nother+npmax]
      * to the swap routine, but the atom id kept getting messed up */
     for (j = 1; j <= icnt; j +=5,jj++) {
       Myparticles[jj].x=rcvbuf[j-1];
       Myparticles[jj].y= rcvbuf[j];
       Myparticles[jj].z=rcvbuf[j +  1];
       Myparticles[jj].cg=rcvbuf[j +  2];
       Myparticles[jj].id=rcvbuf[j +  3];
     
       /* printf("*rcvd[%d] (%f,%f,%f,%f,%f)  ",jj,Myparticles[jj].x, Myparticles[jj].y ,
	* Myparticles[jj].z, Myparticles[jj].cg,Myparticles[jj].id);
	*/
     }
     *nother += icnt / 5;
     if (*nother >= nomax) {
       sprintf(str,"borders():Received too many border atoms(%d), nomax=%d",
	       *nother,nomax);
       error_handler(str);
     }

   } /* for kk */
 } /* for k */
 nslist[*nswap+1] = npnt;
 free(sndbuf);
 free(rcvbuf);
#if TIMEME
 gettimeofday(&(time2),&tzp);
*mytime=swatch(time1,time2);
#endif
 return 0;
} /* borders*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*  no binning, Newton's 3rd law */
/*  N^2 / 2 search for neighbor pairs in my box */
/*  pair stored once in list if atoms i AND j are in my box (and i < j) */
/*  pair stored by me if j is NOT in my box (also stored by node owning j) */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int build_verlet(int *nlocal, int *nother, int *list, Pme2Particle **AllParticles, int *nnlist, int *nlist, PmeVector myprd, double rs, int *atompnt,double *mytime)
{
  Pme2Particle *Myparticles;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double cutsq2,rsq,xms2,yms2,zms2;
  int i,i_2,ii,j,jj,npnt;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
#if (VERBOSE)
  printf(" calling neighbors0()... \n");
#endif   

  xms2= 2. * (myprd.x-rs);
  yms2= 2. * (myprd.y-rs);
  zms2= 2. * (myprd.z-rs);
  cutsq2=rs * rs;

  Myparticles= (*AllParticles); /* pointer setup */
  npnt = 1;
  i = *atompnt;
  for (ii = 1; ii <= *nlocal; ii++) {
    nnlist[ii] = npnt;
    xtmp = Myparticles[i].x;
    ytmp =  Myparticles[i].y;
    ztmp =  Myparticles[i].z;
    j = list[i];
    i_2 = *nlocal + *nother;

    for (jj = ii + 1; jj <= i_2; ++jj) {
      delx = xtmp -  Myparticles[j].x;
      dely = ytmp -  Myparticles[j].y;
      delz = ztmp -  Myparticles[j].z;
            
      delx -= myprd.x * Nint(delx/xms2);
      dely -= myprd.y * Nint(dely/yms2);
      delz -= myprd.z * Nint(delz/zms2);
      
      rsq = delx * delx + dely * dely + delz * delz;


      if (rsq <= cutsq2) {
	nlist[npnt] = j;
	++npnt;
      }
      if (jj <= *nlocal) {
	j = list[j];
      } else {
	++j;
      }
    }
    if (npnt > npmax*nnmax) 
      error_handler("build_verlet():Neighbor list too big");
    i = list[i];
  } /* for ii */
  nnlist[*nlocal+1] = npnt;
  
  /* print the nbr atoms of each atom */
  i= *atompnt;
  for (ii = 1; ii <= *nlocal; ii++) {
    /* printf("\n %% build_verlet():local index atom-%d has nbr: \n",i); */
    for (jj=nnlist[ii]; jj<= nnlist[ii+1]-1; jj++){
      j=nlist[jj];
     /*  printf(" #%d  ",j); */
    }
    i= list[i];
  }
#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif
    return 0;
} /* build_verlet*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*  Link-Cell Binning (Full space 27 cells ) Newton's 3rd law   *
 *  all of mine and nearby atoms binned once                    *
 *  each owned atom i checks 27 surrounding boxes pair stored   *
 *  once in list if atoms i AND j are in my box (and i < j)     * 
 *  pair stored by me if j is NOT in my box                     *
 *  (also stored by node owning j)                              */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int build_linkcell(int *atompnt, int *nlocal, int *nother, PmeVector mbin, PmeVector binsize, 
	      int *list, int *nlist, int *nnlist, Pme2Particle *Myparticles, int **bin, 
	      int **binpnt, PmeVector nbin, PmeVector myprd, PmeVector myprd2, 
	      PmeVector mbinlo, double rs,double *mytime)
{
  static int *mybin,*mybinpnt;
  static int firsttime=0;
  int i_1,i,ii,j,k,npnt,ib;
  int ix,iy,iz,ixx,iyy,izz;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,cutsq2;
  double xms2,yms2,zms2;
  char str[100];
  
#if TIMEME
  struct timeval intime[4];
  struct timezone tzp;
  struct timeval time1,time2;
  gettimeofday(&(time1),&tzp);
#endif

#if (VERBOSE)
  printf(" calling neighbors1()... \n");
#endif

  if (!firsttime) {
    if (!(mybin = (int *) malloc((namax+1) * sizeof(int))) )
      fprintf(stderr,"Error in allocating space for bin Data !!!\n");
    if (!(mybinpnt = (int *) malloc((nbmax+1) * sizeof(int))) )
      fprintf(stderr,"Error in allocating space for binpnt Data !!!\n");
    firsttime ++;
  }
  /* pointer setup */
  (*bin)=mybin;
  (*binpnt)=mybinpnt;

  xms2= 2. * (myprd.x-rs);
  yms2= 2. * (myprd.y-rs);
  zms2= 2. * (myprd.z-rs);
  cutsq2= rs *rs;

  /* DPME_DEBUG */
  /* printf("mbin(%f,%f,%f), mbinlo(%f,%f,%f), binsz(%f,%f,%f),nbin(%f,%f,%f)\n",
     mbin.x,mbin.y,mbin.z,mbinlo.x,mbinlo.y,mbinlo.z, binsize.x,binsize.y,
     binsize.z,nbin.x,nbin.y,nbin.z);
  */

  i_1 = (int)(mbin.x * mbin.y * mbin.z);
  for (i = 1; i <= i_1; ++i) {
    mybinpnt[i] = 0;
  }

#if TIMEME*DPME_DEBUG
  gettimeofday(&intime[0],&tzp);
#endif
 
  i = *atompnt;
  i_1 = *nlocal + *nother;
  for (ii = 1; ii <= i_1; ++ii) {
    ix = (Myparticles[i].x + myprd2.x) / binsize.x;
    iy = (Myparticles[i].y + myprd2.y) / binsize.y;
    iz = (Myparticles[i].z + myprd2.z) / binsize.z;

    ix -= mbinlo.x;
    if (ix < 0)  ix += nbin.x;
    iy -= mbinlo.y;
    if (iy < 0)   iy += nbin.y;
    iz -= mbinlo.z;
    if (iz < 0)  iz += nbin.z;

    if (ix <= mbin.x && iy <= mbin.y && iz <= mbin.z){
      ib = iz * mbin.y * mbin.x + iy * mbin.x + ix + 1;
      mybin[i] = mybinpnt[ib];
      mybinpnt[ib] = i;
    }
    if (ii <= *nlocal) {
      i = list[i];
    } else {
      ++i;
    }
  }

#if TIMEME*DPME_DEBUG
  gettimeofday(&intime[1],&tzp);
#endif


  npnt = 1;
  i = *atompnt;
  for (ii = 1; ii <=  *nlocal; ++ii) {
    nnlist[ii] = npnt;

    xtmp =Myparticles[i].x;
    ytmp =Myparticles[i].y;
    ztmp =Myparticles[i].z;

    ixx = (xtmp + myprd2.x) / binsize.x;
    iyy = (ytmp + myprd2.y) / binsize.y;
    izz = (ztmp + myprd2.z) / binsize.z;
    for (k = 0; k <= 26; ++k) {
      ix = ixx + (k % 3) - 1;
      iy = iyy + ((k / 3) % 3) - 1;
      iz = izz + (k / 9) - 1;
      ix -= mbinlo.x;
      if (ix < 0) {
	ix += nbin.x;
      }
      if (ix == nbin.x) {
	ix = 0;
      }
      iy -= mbinlo.y;
      if (iy < 0) {
	iy += nbin.y;
      }
      if (iy == nbin.y) {
	iy = 0;
      }
      iz -= mbinlo.z;
      if (iz < 0) {
	iz += nbin.z;
      }
      if (iz == nbin.z) {
	iz = 0;
      }
      ib = iz * mbin.y * mbin.x + iy * mbin.x + ix  + 1;
      j = mybinpnt[ib];
      while (j != 0) {
	if (j <= i) j = mybin[j];
	else
	  { 
	    delx = xtmp - Myparticles[j].x;
	    dely = ytmp -  Myparticles[j].y;
	    delz = ztmp -  Myparticles[j].z;

	    delx -= myprd.x * Nint(delx/xms2);
	    dely -= myprd.y * Nint(dely/yms2);
	    delz -= myprd.z * Nint(delz/zms2);
	    rsq = delx * delx + dely * dely + delz * delz;
	    if (rsq <= cutsq2) {
	      nlist[npnt] = j;
	      ++npnt;
	    }
	    j = mybin[j];
	  } /* else j > i */
	
      } /* while j !=0 */
    } /* for k */
    if (npnt > npmax*nnmax) {
      sprintf(str,"Build_Linkcell list too big(%d)...boost nnmax(%d)npmax(%d)!",
	      npnt,nnmax,npmax);
      error_handler(str);
    }
    i = list[i];
  }
  nnlist[*nlocal+1] = npnt;

#if (TIMEME*DPME_DEBUG)
  gettimeofday(&intime[2],&tzp);
  printf("NBR1, time 1st big loop=%f \n",swatch(intime[0],intime[1]));
  printf("NBR1, time 2nd big loop=%f \n",swatch(intime[1],intime[2]));
#endif

#if (TIMEME)
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif
  return 0;
} /* build_linkcell*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* Print error messages  and exits pvm                          */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

int error_handler(char *errstring)
{
  fprintf(stderr,"%s !\n Exiting now ... \n", errstring);
  fflush(stdout);
  pvm_lvgroup(GROUP); 
  pvm_exit();
  exit(0);
  return 0;
 }
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* Config the processors as a 3D mesh (if possible) w/PBC       */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int mesh_3d(int *nx, int *ny, int *nz, int *node, 
	      int *ix, int *iy, int *iz, int *iwest, int *ieast, 
	      int *isouth, int *inorth, int *idown, int *iup)
{
  *ix = *node % *nx;
  *iy = *node / *nx % *ny;
  *iz = *node / *nx / *ny;
  *iwest = *ix - 1;
  if (*iwest == -1) {
    *iwest = *nx - 1;
  }
  *ieast = *ix + 1;
  if (*ieast == *nx) {
    *ieast = 0;
  }
  *isouth = *iy - 1;
  if (*isouth == -1) {
    *isouth = *ny - 1;
  }
  *inorth = *iy + 1;
  if (*inorth == *ny) {
    *inorth = 0;
  }
  *idown = *iz - 1;
  if (*idown == -1) {
    *idown = *nz - 1;
  }
  *iup = *iz + 1;
  if (*iup == *nz) {
    *iup = 0;
  }
  *iwest = *nx * *ny * *iz + *nx * *iy + *iwest;
  *ieast = *nx * *ny * *iz + *nx * *iy + *ieast;
  *isouth = *nx * *ny * *iz + *nx * *isouth + *ix;
  *inorth = *nx * *ny * *iz + *nx * *inorth + *ix;
  *idown = *nx * *ny * *idown + *nx * *iy + *ix;
  *iup = *nx * *ny * *iup + *nx * *iy + *ix;

  return 0;
} /* mesh_3d*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* calculate ewald_direct potential (w/o corrections) ...
   This requires doing a global reduction to collect from nodes 
   ************** DON'T USE IT, ONLY FOR TESTING  **************
   USE DIR_FORCE() to calc force and energy.
 */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double dir_energy(int *atompnt, int *nlocal, int *nnlist, int *nlist, int *list, 
		  MYPROC myproc, PmeVector myprd, double cutoff, double ewaldcof, 
		  Pme2Particle *Myparticles, PmeVector mymc2,int *tidarray)
{
  int i,i_2,ii,j,k; 
  double direct_eng,factor,rsq,delx,dely,delz;
  double cutsq1,myr,myx,ewaldpot;
  extern double sqrt(double), erfc(double);
  
  cutsq1= cutoff * cutoff;
  direct_eng = 0.0;
  i = *atompnt;
  for (ii = 1; ii <=  *nlocal; ++ii) {
    i_2 = nnlist[ii+1] - 1;
    for (j = nnlist[ii]; j <= i_2; ++j) {
      k=nlist[j];
      delx= Myparticles[i].x-Myparticles[k].x;
      dely= Myparticles[i].y-Myparticles[k].y;
      delz= Myparticles[i].z-Myparticles[k].z;
     
      delx -= myprd.x * Nint(delx/mymc2.x);
      dely -= myprd.y * Nint(dely/mymc2.y);
      delz -= myprd.z * Nint(delz/mymc2.z);
      rsq= delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq1) {
	factor = 1.0;
	if (k > npmax)   factor = .5;
	myr=sqrt(rsq);
	myx= myr * ewaldcof;
	ewaldpot = erfc(myx) / myr;
	direct_eng +=  Myparticles[i].cg * Myparticles[k].cg * factor * ewaldpot;
      }
    } /* for j */
    i = list[i];
  } /* for ii */
  merge_d2(&direct_eng,myproc,tidarray,MSG_300);
  return(direct_eng);
} /* dir_energy*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* calculate the direct_sum force and energy on atoms you own, use Newtons
 * 3rd law to avoid double computing. 
 */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double dir_force(int *atompnt, int *nlocal, int *list, int *nlist, int *nnlist, 
		 Pme2Particle *Myparticles, PmeVector myprd, PmeVector mymc2, 
		 double cutoff, double ewaldcof, PmeVector **directF, double *virial, 
		 MYPROC myproc,int time_count, int tsteps, int *tidarray,double *mytime)
                   
{
  int i,i_2,ii,j,k;
  double xtmp,ytmp,ztmp,delx,dely,delz,cutsq1,rsq;
  static int firsttime=0;
  static PmeVector *mydirectF;
  double myr,myx,merfc,term,tempcg,tmp,tempx,tempy,tempz;
  double factor, ewaldpot,direct_eng; /* testing ewaldpot */
  double inv_myr,pi_ewaldcof;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif

#if VIRIAL
  --virial;
  for (i =1; i <=6; ++i) {
    virial[i] = 0.;
  }
#endif

  pi_ewaldcof= TwoBySqrtPi * ewaldcof;
  direct_eng = 0.0; /* testing ewaldpot */

  if (!firsttime) {
#if DPME_DEBUG
    printf(" Allocating Memory to directF \n");
#endif
    /* alloc force array firsttime only, later can realloc on renieghbor */
    if (!(mydirectF = (PmeVector *) malloc((npmax+1) * sizeof(PmeVector))) )
      fprintf(stderr,"Error in allocating space for directF  Data !!!\n");
    firsttime ++;
  }
  (*directF)= mydirectF;
  cutsq1= cutoff * cutoff;
  i = *atompnt;
  for (ii = 1; ii <= *nlocal ; ++ii) {
    mydirectF[i].x= 0.;
    mydirectF[i].y= 0.;  
    mydirectF[i].z= 0.;
    i = list[i];
  }
  i = *atompnt;

  for (ii = 1; ii <=  *nlocal ; ++ii) {
    xtmp = Myparticles[i].x;
    ytmp = Myparticles[i].y;
    ztmp = Myparticles[i].z;
    i_2 = nnlist[ii+1] - 1;
    
    for (k = nnlist[ii]; k <= i_2; ++k) {
      j = nlist[k];
      delx = xtmp - Myparticles[j].x;
      dely = ytmp - Myparticles[j].y;
      delz = ztmp - Myparticles[j].z;
      
      delx -= myprd.x * Nint(delx/mymc2.x);
      dely -= myprd.y * Nint(dely/mymc2.y);
      delz -= myprd.z * Nint(delz/mymc2.z);

      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq <= cutsq1) {
	myr= sqrt(rsq);
	inv_myr= 1.0/myr;
	myx= myr * ewaldcof;
	merfc= erfc(myx);
	term= pi_ewaldcof  * exp(-(myx*myx))*inv_myr + merfc*inv_myr*inv_myr;
	tempcg =  Myparticles[i].cg * Myparticles[j].cg ;
	tmp= tempcg * term * inv_myr;
	tempx= tmp * delx;
	tempy= tmp * dely;
	tempz= tmp * delz;

        mydirectF[i].x += tempx;
	mydirectF[i].y += tempy;
	mydirectF[i].z += tempz;
	factor = .5; /* testing ewaldpot */
	if (j <= npmax) { /* if j is local store it */
	  mydirectF[j].x -= tempx;
	  mydirectF[j].y -= tempy;
	  mydirectF[j].z -= tempz;
	  factor = 1.0; /* testing ewaldpot */
	}
	ewaldpot= merfc * inv_myr; /* testing Pot */
	direct_eng += ewaldpot * factor * tempcg; 
#if VIRIAL
	virial[1] += -tempx * delx * factor;
	virial[2] += -tempy * delx * factor;
	virial[3] += -tempz * delx * factor;
	virial[4] += -tempy * dely * factor;
	virial[5] += -tempz * dely * factor;
	virial[6] += -tempz * delz * factor; 
#endif
      }
    } /* for k */
    i = list[i];
  } /* for ii */
  
  /* collect all contributions */
  /* if (time_count==tsteps) */ /* can do this only at last step or every x-steps */
#if DPME_DEBUG
  printf("$$ node %d had dir_eng(uncorrected)=%f \n",myproc.node,direct_eng);
#endif
  merge_d2(&direct_eng,myproc,tidarray,MSG_300);

#if VERBOSE
  printf("\nAt t=%d..........direct_energy2= %f \n",time_count,direct_eng);
#endif

#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif
  return (direct_eng);
} /* dir_force*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* allocate a double vector w/ range v[nl..nh]  with offset     *
*  free using free(v-1)                                         */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double *dvector( int nl, int nh)
{
  double *v;

#if SUN4SOL2
  v= valloc((nh-nl+2)*sizeof(double)); /* this works but not ANSI*/
#else
  v=(double *)malloc(((nh-nl+2)*sizeof(double)));
#endif

  
  if (!v) { 
    fprintf(stderr,"ERROR in allocating a vector,exiting!!!!\n");
    exit(0);
  }
  return v-nl+1;
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* allocate a double vector (same as dvector() ) without offset */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double *dvector2( int nl, int nh)
{
  double *v;
  v=(double *)malloc(((nh-nl+1)*sizeof(double)));
  if (!v) { 
    fprintf(stderr,"ERROR in allocating a vector,exiting!!!!\n");
    exit(0);
  }
  return v-nl;
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* compute wall clock time when using gettimeofday() */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double swatch(struct timeval st,struct timeval end)
{
  double total;
  
 total= (                 (double)end.tv_sec
                        - (double)st.tv_sec 
                        + ((double)end.tv_usec
                        - (double)st.tv_usec) / 1000000) ;
return total;
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* user clock time */
double uswatch(struct tms startbuf, struct tms endbuf)
{
  
  return  (endbuf.tms_utime - startbuf.tms_utime)/(float)CLK_TCK;
}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* this routines perturb the x,y,z coordinates in some fasion to emulate 
 *  time integration- this eventually sh'd be replaced by an integrator 
 * Note that we allow for periodic wrap-around, ie if particle goes out of
 * the box its wrapped back into it from the other side ( tested works fine)
 * only trouble is when orginal coordinate ==0
 */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void  update_coordinates(int *atompnt, int nlocal, int *list, Pme2Particle *particlelist,
			 PmeVector prd, PmeVector prd2) 
{
  int i,ii;

#if VERBOSE
  printf ("calling update_coordinates()...\n");
#endif

  i = *atompnt;

  for (ii=1; ii<=nlocal; ii++) {
    particlelist[i].x += 0.001;
    if (particlelist[i].x > prd2.x) {
      particlelist[i].x -=prd.x;
      printf("CORRECTING atom-x1 %f out of box \n",particlelist[i].id);
    }
    else 
      if (particlelist[i].x <= -prd2.x) {
	particlelist[i].x +=prd.x; 
	printf("CORRECTING atom-x2 %f out of box \n",particlelist[i].id);
      }
 
     particlelist[i].y  += 0.001;  
     if (particlelist[i].y > prd2.y) { 
       particlelist[i].y -=prd.y; 
       printf("CORRECTING atom-y1 %f out of box \n",particlelist[i].id);
     }
     else 
       if (particlelist[i].y <= -prd2.y) { 
	 particlelist[i].y +=prd.y; 
	 printf("CORRECTING atom-y2 %f out of box \n",particlelist[i].id);
       }
     
     particlelist[i].z += 0.001;
     if (particlelist[i].z > prd2.z) { 
       particlelist[i].z -=prd.z;
       printf("CORRECTING atom-z1 %f out of box \n",particlelist[i].id);
     }
     else 
       if (particlelist[i].z <= -prd2.z) { 
	 particlelist[i].z +=prd.z; 
	 printf("CORRECTING atom-z2 %f out of box \n",particlelist[i].id);
       }
     i = list[i];
   } /* for */
  
} /* update_coordinate */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* Correct the energy/forces assuming a system of water.
 * This is equilvalent to adjust_direct and adjust_recip in version 1
 * and are provided by way of example.
 * Only adjust the energy by 50% if an atom is not owned by a PE.
 */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double adjust_dir_recip(int *atompnt, int *nlocal, int *list, int *nlist, int *nnlist, 
			Pme2Particle *Myparticles, PmeVector myprd, PmeVector mymc2, 
			double cutoff, double ewaldcof, PmeVector *mydirectF, 
			double *recip,	MYPROC myproc,int time_count, int tsteps,
			int *tidarray, 	double *dir_vir,double *adj_vir, 
			double *my_adj_rcp_eng, double *my_adj_dir_eng, 
			PmeVector **adjustF,double *mytime)
                   
{
  static int firsttime=0;
  static PmeVector *myadjustF;
  float full_interact;
  double adj_dir_eng,adj_recip_eng;
  int i,i_2,ii,j,k,need1,need2;
  int atom_type,atom_id,found;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif

  if (!firsttime) {
#if DPME_DEBUG
    printf(" Allocating Memory to adjustF \n");
#endif
    /* alloc force array firsttime only, later can realloc on renieghbor */
    if (!(myadjustF = (PmeVector *) malloc((npmax+1) * sizeof(PmeVector))) )
      fprintf(stderr,"Error in allocating space for adjustF  Data !!!\n");
    firsttime ++;
  }
  (*adjustF)= myadjustF;


  adj_dir_eng = 0.0; 
  adj_recip_eng = 0.0; 

  for (i=0; i<6;i++)
    adj_vir[i]=0.0; /* initialize virial every time step */

  i = *atompnt;
  for (ii = 1; ii <= *nlocal ; ++ii) {
    myadjustF[i].x= 0.;
    myadjustF[i].y= 0.;  
    myadjustF[i].z= 0.;
    i = list[i];
  }

  i = *atompnt;
  for (ii = 1; ii <=  *nlocal ; ++ii) {

    /* Water has 3 atoms, Oxygen, H1, H2 they ordered as O,H,H,O,H,H,O,... 
     *  in small.pdb , so we can use the id tag to figure out which water it is
     * if atom_type=1 then  Oxygen; atom_type=2 then H1; atom_type=0 then H2.
     * Each atom has 2 corrections to make, namely O-H1, O-H2, H1-H2.
     */
    atom_id=  (int)(Myparticles[i].id); 
    atom_type= atom_id%3; 
    switch(atom_type) 
      {
      case 1: {
	need1=  atom_id + 1 ;
	need2=  need1 +1;
	break;
      }
      case 2: {
	need1=  atom_id - 1 ;
	need2=  need1 + 2 ;
	break;
      }
      case 0: {
	need1=  atom_id - 1 ;
	need2=  need1 - 1;
	break;
      }
      }
    found=2; /* count # of bounded interactions */

    i_2 = nnlist[ii+1] - 1;
    for (k = nnlist[ii]; (k <= i_2) && (found!=0); ++k) { 

      /* if found==0 then mv on to next Myparticle[i] */
      j = nlist[k];
      if (j <= npmax ) 
	full_interact=1.0; /* store the force interaction */
      else 
	full_interact=0.5;
     
      if (Myparticles[j].id == need1) {
	correct_water_dir(i,j,Myparticles[i],Myparticles[j],cutoff,ewaldcof,myprd,mymc2,
			  &adj_dir_eng,mydirectF,recip,dir_vir,full_interact);
	correct_water_recip(i,j,Myparticles[i],Myparticles[j],ewaldcof,
			    &adj_recip_eng,myadjustF,adj_vir,full_interact);
	found--;
      }
      else 
	if (Myparticles[j].id == need2) {
	  correct_water_dir(i,j,Myparticles[i],Myparticles[j],cutoff,ewaldcof,myprd,mymc2,
			    &adj_dir_eng,mydirectF,recip,dir_vir,full_interact);
	  correct_water_recip(i,j,Myparticles[i],Myparticles[j],ewaldcof,
			      &adj_recip_eng,myadjustF,adj_vir,full_interact);
	  found--;
	}  
    } /* for k */
    i = list[i];
  } /* for ii */
  
  /* collect all contributions */
  /* if (time_count==tsteps) */ /* can do this only at last step or every x-steps */
  merge_d2(&adj_recip_eng,myproc,tidarray,MSG_320); 
  merge_d2(&adj_dir_eng,myproc,tidarray,MSG_350); 
#if VERBOSE
  printf(" The adjust_dir_energy=%f \n",adj_dir_eng);
  printf(" The adjust_recip_energy =%f \n",adj_recip_eng);
#endif

  *my_adj_rcp_eng = adj_recip_eng;
  *my_adj_dir_eng = adj_dir_eng;
#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif
  return (0);
} /* adj_dir_eng*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* this routine will correct the Direct_sum forces to exclude the bonded
* interactions as was done in adjust_recip() for the recip_sum.
* Note: this works only for water and provided by way of example */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int correct_water_dir( int i, int j,Pme2Particle particle1, Pme2Particle particle2, 
		  double cutoff, double ewaldcof, PmeVector myprd, PmeVector mymc2, double *ene, 
		  PmeVector *mydirectF, double *recip,double *virial,float full_interact)
{
  double grad[3];
  double r,cutoffsqr;
  register double tempo; /* ayt aug 8 */
  double crd[3], xij, yij, zij, pot; 
 
#if VIRIAL
   --virial;
#endif

  /* Function Body */
  
  cutoffsqr=  cutoff *  cutoff;
  
  xij =  particle1.x - particle2.x;
  yij =  particle1.y -  particle2.y;
  zij =  particle1.z -  particle2.z;

  xij -= myprd.x * Nint(xij/mymc2.x);
  yij -= myprd.y * Nint(yij/mymc2.y);
  zij -= myprd.z * Nint(zij/mymc2.z);
  
  r = (xij*xij + yij*yij + zij*zij); /* this is r_sqr */

/* it seems there is no need to check cutoff , we know that if choose the atoms properly
 * r will be always <= cutoffsqr at least for water where O-H ~ 1 A, leave it for genrality
 */

  if (r <= cutoffsqr) { 
    r=sqrt(r); 
    xij -= myprd.x * Nint(xij * recip[0]);
    yij -= myprd.y * Nint(yij * recip[4]);
    zij -= myprd.z * Nint(zij * recip[8]);
    crd[0] = xij;
    crd[1] = yij;
    crd[2] = zij;
  
    ew_direct(&ewaldcof, crd, &pot, grad,&r);
    
    /* if full_interact (ie I own both particles) 
     *  then 100% intercation else 50% 
     */
    tempo=  - particle1.cg * particle2.cg; /* note minus sign , ayt aug 9 */
    *ene -= full_interact * (-tempo) * pot; /* note minus sign , ayt aug 9 */
  
    mydirectF[i].x -=tempo * grad[0];  
    mydirectF[i].y -=tempo * grad[1];
    mydirectF[i].z -=tempo * grad[2];
  

    if (full_interact==1.0) {
      mydirectF[j].x +=tempo * grad[0];
      mydirectF[j].y +=tempo * grad[1];
      mydirectF[j].z +=tempo * grad[2];
    }
#if VIRIAL
   
    virial[1] += tempo * xij * grad[0] * full_interact;
    virial[2] += tempo * xij * grad[1] * full_interact;
    virial[3] += tempo * xij * grad[2] * full_interact;
    virial[4] += tempo * yij * grad[1] * full_interact;
    virial[5] += tempo * yij * grad[2] * full_interact;
    virial[6] += tempo * zij * grad[2] * full_interact;
#endif
  } 

  return 0;
} /* corect_water_dir */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* this routine evaluates a pair of Direct sum interaction */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
int ew_direct(double *ewaldcof, double *crd, 
	double *pot, double *grad, double *r)
{
    /* Local variables */
    double term;
    int i;
    double x, merfc; 
    
    /* Parameter adjustments */
    --crd;
    --grad;

    x = *r * *ewaldcof;
    merfc=erfc(x);
    *pot = merfc / *r;
    term = (TwoBySqrtPi) * *ewaldcof * exp(-(x * x)) / *r + merfc / (*r * *r);
    for (i = 1; i <= 3; ++i) {
	grad[i] = -term * crd[i] / *r;
    }
    return 0;
} /* ew_direct */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* this routine is called by adjust_recip() used to subtract bonded interactions 
 * this was formerly get_adjust_pair()
 */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
int correct_water_recip( int i,  int j, Pme2Particle particle1, Pme2Particle particle2,
			double ewaldcof, double *ene, PmeVector *afparticle, double *vir,
			float full_interact)
{
  double grad[3],temp;
  double crd[3], pot;
  
#if VIRIAL
  --vir;
#endif
  
  /* Function Body */
  crd[0] = particle1.x - particle2.x;
  crd[1] = particle1.y -  particle2.y;
  crd[2] = particle1.z -  particle2.z;
  
  ew_adjust(&ewaldcof, crd, &pot, grad); 
  temp = particle1.cg * particle2.cg;
  *ene +=  full_interact * temp * pot;
  
  afparticle[i].x -=  temp * grad[0];
  afparticle[i].y -=  temp * grad[1];
  afparticle[i].z -=  temp * grad[2];

  /* if full_interact (ie I own both particles) 
   *  then 100% intercation else 50% 
   */
  if (full_interact==1.0) {
    afparticle[j].x +=  temp  * grad[0];
    afparticle[j].y +=  temp * grad[1];
    afparticle[j].z +=  temp  * grad[2];
  }
#if VIRIAL
    vir[1] += particle1.cg * particle2.cg * crd[0] * grad[0] * full_interact;
    vir[2] += particle1.cg * particle2.cg * crd[0] * grad[1] * full_interact;
    vir[3] += particle1.cg * particle2.cg * crd[0] * grad[2] * full_interact;
    vir[4] += particle1.cg * particle2.cg * crd[1] * grad[1] * full_interact;
    vir[5] += particle1.cg * particle2.cg * crd[1] * grad[2] * full_interact;
    vir[6] += particle1.cg * particle2.cg * crd[2] * grad[2] * full_interact;
#endif

    return 0;
} /* correct_water_recip */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* this routine called by get_adjust_pair() used in correcting recip_sum */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
 int ew_adjust(double *ewaldcof, double *crd, double *pot, double *grad)
{
    /* Builtin functions */
    double sqrt(double), exp(double);

    /* Local variables */
    double term;
    int i;
    double r, x, merfc;
    
    /* Parameter adjustments */
    --crd;
    --grad;

/* assume x,y,z are minimum image */
    r = sqrt( crd[1] * crd[1]  +crd[2] * crd[2] + crd[3] * crd[3] );
    x = r * *ewaldcof;
    merfc= erfc(x);
    *pot = -(1. - merfc) / r;
    term = -(TwoBySqrtPi) * *ewaldcof * exp(-(x * x )) / r + (1. - merfc)
	     / (r * r);
    for (i = 1; i <= 3; ++i) {
      grad[i] = term * crd[i] / r;
    }
    return 0;
} /* ew_adjust */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* calculates the self term of Ewald sum which only contributes to the energy 
* and not forces */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
double self(Pme2Particle *particlelist,  int numatoms,int atompnt, int *list,
	 double ewaldcof, MYPROC myproc, int *tidarray,double *mytime)
{
  /* Local variables */
  int i,n;
  double ee;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
  
  ee = 0.0;
  
  n = atompnt;
  
  for (i = 1; i <= ( numatoms); ++i) {
    ee += (particlelist[n].cg * particlelist[n].cg);
    n=list[n];
  }
  
  ee *= -ewaldcof / SQRT_PI;
  
  /* global sum across processors, used tag-2 to avoid confusion with previous merge */
  merge_d2(&ee,myproc,tidarray,MSG_380);
#if VERBOSE
  printf(" The self_energy=%f \n",ee);   
#endif
#if TIMEME
  gettimeofday(&(time2),&tzp);
 *mytime=swatch(time1,time2);
#endif

  return(ee);
} /* self */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* collect all virial from all processors */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
void collect_virials(double *adj_vir, double *dir_vir, double *rec_vir,
		     MYPROC myproc, int *tidarray,double *mytime)
{
  int i;
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
#if VERBOSE
  printf(" collecting virials ....\n"); 
#endif
  merge_v2(adj_vir,myproc,tidarray,6,1);
  merge_v2(dir_vir,myproc,tidarray,6,2);
  /* no need to recip_vir because each proc has all recip interactions 
   * merge_v2(rec_vir,myproc,tidarray,6,3);
   */
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
/* reads in a formatted output file generated by exact Ewald sum and has
 * the dir forces followed by adjust-recip and recip forces 
 * offset can be 1,2,or 3 depending which forces you would like to read in
*/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
int read_file(int numatoms,int offset,PmeVector *myforce)
{
  FILE  *infile; 
  int i,atomid; 
  char DataFile[80];
  double x,y,z;

  /* open data  file */
  sprintf(DataFile,"%s%s",DATADIR,"nourexact.txt");
  if ( (infile = fopen(DataFile,"r") ) == NULL) {
    fprintf(stderr,"Unable to open the input file: !!! \n");
    exit(1); 
  }
  /* seek the correct segment */
  for (i=1; i<=(offset-1)*numatoms;i++)
      fscanf(infile,"%d%lf%lf%lf",&atomid,&x,&y,&z);
  for (i=1+(offset-1)*numatoms;i<=(offset*numatoms);i++){
    fscanf(infile,"%d%lf%lf%lf",&atomid,&x,&y,&z);
    myforce[atomid].x=x;
    myforce[atomid].y=y;
    myforce[atomid].z=z;
   
  }

    fclose(infile);
    return 0;
} /* read_file */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* after u read the exact data file using (read_file) you can compare
* your force results with the exact force results using this routine
* the BaseF vector is used as the denomintor in RMS error formula and thus
* can be either ExactF or the sum of all recip/dir/adjust forces
*/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
int check_forces(int numatoms, int atompnt, int *list, PmeVector *ApprxF, 
		 PmeVector *ExactF, PmeVector *BaseF)
{
    /* System generated locals */
    
    double d_1, d_2, d_3; /* ok */

    /* Local variables */
    int i,n,atomidx,atomidy,atomidz;
    double rms, rms_den, rms_num;
    register double xmax, ymax, zmax;
  
    rms_num = 0.;
    rms_den = 0.;
    xmax=0.;
    ymax=0.;
    zmax=0.;
    
    i=atompnt;

    for (n = 1; n <= (numatoms); ++n) {
      d_1 = ApprxF[i].x - ExactF[i].x;
      d_2 =  ApprxF[i].y - ExactF[i].y;
      d_3 =  ApprxF[i].z - ExactF[i].z;
      if (fabs(d_1) > xmax) {
	xmax=fabs(d_1);
	atomidx=i;
      }
      if (fabs(d_2) > ymax) { 
	ymax=fabs(d_2);
	atomidy=i;
      }
      if (fabs(d_3) > zmax) { 
	zmax=fabs(d_3);
	atomidz=i;
      }
      /* printf(" at i = %d , fx1=  %lf , fx2= %lf \n",i,fx1[i],fx2[i]); */ 
      rms_num = rms_num + d_1 * d_1 + d_2 * d_2 + d_3 * d_3;
      rms_den = rms_den +(BaseF[i].x*BaseF[i].x) +(BaseF[i].y*BaseF[i].y)+
	(BaseF[i].z*BaseF[i].z);
      
      i=list[i];
    } /* for */
    
    rms = sqrt(rms_num / rms_den);
    printf("_________________ERROR INFO____________ \n");
    printf("rms (relative) force err  = %3.10f , avegarged= %3.10f\n",
	   rms,rms/ numatoms);
    printf("rms absolute force err averaged= %3.10f \n",sqrt(rms_num)/ numatoms);
    printf("Max error in X direction = +/-  %3.10f at %d \n",xmax,atomidx); 
    printf("Max error in Y direction = +/-  %3.10f at %d \n",ymax,atomidy); 
    printf("Max error in Z direction = +/-  %3.10f at %d \n",zmax,atomidz); 
    printf("----------------------------------------\n");
    return 0;
} 
/* check_forces */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/* this routine compares the approximate forces computed here versus 
 *"exact" Ewald forces computed with another program. This routine is just
 * a utility to compute the RMS errors.
 */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
int verify_results(int nlocal, int atompnt, int *list,
		   PmeVector *directF, PmeVector *recipF, PmeVector *adjustF)
{
  int i, ii,offset;
  PmeVector *ExactF,*RExactF, *DExactF,*AExactF, *ApprxF;

#if VERBOSE
  printf("calling verify results...\n");
#endif

   /* first read the exact date file, offset =1 if direct_force, 2=adjust_forces
    * or 3=if recip_forces. Don't forget to pass the proper array in check_force
    */
   offset=3; /*recipF*/
   if (!(RExactF = (PmeVector *) malloc((nlocal+1) * sizeof(PmeVector))) )
     fprintf(stderr,"Error in allocating space for ExactF  Data !!!\n");
   read_file(nlocal,offset,RExactF);

   offset=1;/*directF*/
  if (!(DExactF = (PmeVector *) malloc((nlocal+1) * sizeof(PmeVector))) )
     fprintf(stderr,"Error in allocating space for ExactF  Data !!!\n");
   read_file(nlocal,offset,DExactF);
  
   offset=2; /*adjustF */
   if (!(AExactF = (PmeVector *) malloc((nlocal+1) * sizeof(PmeVector))) )
     fprintf(stderr,"Error in allocating space for ExactF  Data !!!\n");
   read_file(nlocal,offset,AExactF);

  if (!(ExactF = (PmeVector *) malloc((nlocal+1) * sizeof(PmeVector))) )
     fprintf(stderr,"Error in allocating space for ExactF  Data !!!\n"); 
  if (!(ApprxF = (PmeVector *) malloc((nlocal+1) * sizeof(PmeVector))) )
     fprintf(stderr,"Error in allocating space for ApprxF  Data !!!\n");

   i=atompnt;
   for (ii=1;ii<=nlocal;ii++){
     ExactF[i].x= DExactF[i].x+ RExactF[i].x+ AExactF[i].x;
     ExactF[i].y= DExactF[i].y+ RExactF[i].y+ AExactF[i].y;  
     ExactF[i].z= DExactF[i].z+ RExactF[i].z+ AExactF[i].z;
     ApprxF[i].x= directF[i].x + recipF[i].x + adjustF[i].x;
     ApprxF[i].y= directF[i].y + recipF[i].y + adjustF[i].y;
     ApprxF[i].z= directF[i].z + recipF[i].z + adjustF[i].z;
     i=list[i];
   }
  printf("Checking recip_force vs exact_recip in terms of total force\n");
  check_forces(nlocal,atompnt,list,recipF,RExactF,ExactF); 
  printf("Checking apprxtotal vs exact_total in terms of total force\n");
  check_forces(nlocal,atompnt,list,ApprxF,ExactF,ExactF); 

  free(RExactF); free(DExactF);  free(AExactF); 
  free(ExactF);free(ApprxF);

  return(0);
}
/* verify_results */
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
