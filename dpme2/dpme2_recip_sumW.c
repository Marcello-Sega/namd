/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
/* $Id: dpme2_recip_sumW.c,v 1.1 1997/04/23 17:05:19 nealk Exp $
 */


/* This is a modified version for dpme2 - by ayt 4/16/96  */
/* modified for Master/Worker on COW 3/21/97 */
/* this should work as long as u pass nlocal = numatoms
 * this code uses a generic 3D fft (serial) 
 */
/************************************************************
* written by A. (Nour) Toukmaji,  1995
*  Duke University, ECE Dept 
************************************************************
This program will perform the Reciprocal_space sum of PME
method(T. Darden NIEHS N.C.). It handles all the calls to recip_sum domain.
Called subroutines are both in utility.c and pme_recip.a
*************************************************************/

#include "dpme2.h"

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* This routine will setup the data structures  sizes for the 
* recip-sum calculations of PME method. More importantly, packup
* fr1,2,3[] arrays (fractional coordinates) and the charge array
* and broadcast to all processors.
* This should be called once.
*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
int recip_sum_setup2(int numatoms,int nlocal, 
		     double *recip, Pme2Particle *Myparticles,
		     Grid nfftgrd, int order, MYPROC myproc, double **fr1, double **fr2, 
		     double **fr3, double *fftable, double *bsp_mod1, double *bsp_mod2, 
		     double *bsp_mod3)
{
  int nfft1,nfft2,nfft3,sizfftab, sizffwrk;
  static int firsttime=0;
  double dumffwork[3];
  static double *myfr1,*myfr2,*myfr3;
  


  /* Grid Dimensions */
  nfft1 = nfftgrd.x; /* assuming an orthogonal system */
  nfft2 = nfftgrd.y;
  nfft3 = nfftgrd.z;
  
  if (!firsttime){ 
#if DPME_DEBUG
    printf("  ** Alloc mem to myfr123...\n");
#endif
    myfr1=dvector2(0,numatoms); /* fractional coordinates of all atoms*/
    myfr2=dvector2(0,numatoms); /* the first nlocal elements */
    myfr3=dvector2(0,numatoms); /* of fr123 are locol to mynode */
  }

  (*fr1)=myfr1;  (*fr2)=myfr2;  (*fr3)=myfr3; 

     
  /* u don't need tmpcg, but u must chagne get_scaled_frcac2 */

  /* initialize the fr[] every time step, set fr of my atoms only */

  get_scaled_fractionals2(&nlocal,&Myparticles[1], recip, &nfft1, &nfft2, &nfft3, 
			  myfr1, myfr2, myfr3);


  /* initialize the bsp coeffs */
  if (!firsttime){ 
    pmesh_kspace_setup(bsp_mod1, bsp_mod2,bsp_mod3,fftable,dumffwork, 
		       &nfft1, &nfft2, &nfft3, &order, &sizfftab, &sizffwrk);
    /* for (i=0;i<(3*(4*nfft+15));i++) fprintf(stderr,"%f\n", fftable[i]); */
    firsttime++; /* only now its updated */
  }

  return 0;  
} /* end recip_sum_setup2 */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* This routine calculates the recip_sum force/energy.
 * Here we recieve the bcast of fractional coordinates and construct
 * a full Q matrix and do as usual the FFT stuf. This is redundant
 * but avoids using a central master to do this and limits the 
 * communication to one multicast/timestep ( no return messages ).
*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double recip_sum_calc2(double ewaldcof, double volume, double *recip, MYPROC myproc, 
		       int nlocal, 
		       double **fr1, double **fr2, double **fr3, int numatoms,
		       int order, Grid nfftgrd, Pme2Particle *Myparticles, double *fftable,
		       double *virial, double *bsp_mod1, double *bsp_mod2,
		       double *bsp_mod3, PmeVector *rfparticle,int last_step)
{
  double *myfr1,*myfr2,*myfr3;
  static double *theta1, *theta2, *theta3,*dtheta1,*dtheta2,*dtheta3;
  static int nfftable,nffwork;
  static int nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3;
  int siz_q;
  double eer; /* recip energy */
  static int firsttime=0;
#if TIMEME
  struct timeval intime[8];
  struct timezone tzp;
  extern double swatch(struct timeval st, struct timeval end); 
#endif
  static double *q; /* the charge grid */
  static double *ffwork; /*scratch space for fft */

  /* fr123 are the fractional coordinates */
  myfr1=(*fr1);   myfr2=(*fr2);   myfr3=(*fr3); 
  
  if(!firsttime){
    /* alloc arrays for spline interpolation */
    theta1=dvector(0,numatoms*order);
    theta2=dvector(0,numatoms*order);
    theta3=dvector(0,numatoms*order);
    dtheta1=dvector(0,numatoms*order);
    dtheta2=dvector(0,numatoms*order);
    dtheta3=dvector(0,numatoms*order);

    /* Grid Dimensions will be used in allocating q */
    nfft1 = nfftgrd.x; /* assuming an orthogonal system */
    nfft2 = nfftgrd.y;
    nfft3 = nfftgrd.z;
    /* to get the size of q,fold in pmesh_kspace_get_sizes &get_fftdims*/
    /* u only need to do this once */
    nffwork=  max(nfft3,(max(nfft1,nfft2))); /* added 12/27 */
    nfftable= 4* nffwork + 15; 

    nfftdim1=nfft1;  
    nfftdim2=nfft2;  
    nfftdim3=nfft3;
    
    if (nfft1 == (2*(nfft1/2))) nfftdim1= nfft1+1;
    if (nfft2 == (2*(nfft2/2))) nfftdim2= nfft2+1;
    if (nfft3 == (2*(nfft3/2))) nfftdim3= nfft3+1;
    siz_q=2*nfftdim1*nfftdim2*nfftdim3;
    q=dvector(0,siz_q); /* alloc q  now */
#if (DEBUG_DPME)
    printf(" alloc  ffwork  now of size %i\n",siz_q);
#endif
    ffwork= dvector(0,siz_q);  /* alloc ffwork  now */

    firsttime++;
  }
 
#if TIMEME
  gettimeofday(&intime[1],&tzp);
#endif 

  /* verified that fr123 theta1 are composed correctly */
  get_bspline_coeffs2(&numatoms, myfr1,myfr2,myfr3,&order,
		      theta1,theta2,theta3,dtheta1,dtheta2,dtheta3);

#if TIMEME
  gettimeofday(&intime[2],&tzp);
#endif

  fill_charge_grid2(&numatoms, nlocal,
		    Myparticles, theta1, theta2, 
		    theta3,  myfr1, myfr2, myfr3, &order, &nfft1, &nfft2, &nfft3, 
		    &nfftdim1, &nfftdim2, &nfftdim3, q);
 

#if TIMEME
  gettimeofday(&intime[3],&tzp);
#endif

  fft_back((doublecomplex *)q, fftable, (doublecomplex *)ffwork, &nfft1, &nfft2, &nfft3, &nfftdim1,
	   &nfftdim2, &nfftdim3, &nfftable, &nffwork);


#if TIMEME
  gettimeofday(&intime[4],&tzp);
#endif

  scalar_sum(q, &ewaldcof, &volume, recip, bsp_mod1, 
	     bsp_mod2, bsp_mod3, &nfft1, &nfft2, &nfft3, &nfftdim1, 
	     &nfftdim2, &nfftdim3, &eer, virial); 
#if TIMEME
  gettimeofday(&intime[5],&tzp);
#endif  
  fft_forward((doublecomplex *)q, fftable, (doublecomplex *)ffwork, &nfft1, &nfft2, &nfft3,
	      &nfftdim1, &nfftdim2, &nfftdim3, &nfftable, &nffwork);
#if TIMEME
  gettimeofday(&intime[6],&tzp);
#endif   
  grad_sum(&nlocal,
	   Myparticles, recip, theta1, theta2, theta3,
	   dtheta1, dtheta2, dtheta3, rfparticle,
	   myfr1, myfr2, myfr3, &order, &nfft1, &nfft2, &nfft3, 
	   &nfftdim1, &nfftdim2, &nfftdim3, q);
#if TIMEME
  gettimeofday(&intime[7],&tzp);
#endif 
  /************************************************/
#if (DPME_DEBUG)
  if (myproc.node==0) { 
    /*
       for  (i=0; i<siz_q;i++)
       fprintf(stderr,"%f\n",q[i]);
       */  
    /* 
     *  ok verified fr1 to be composed right, net necesserly in order 
     *  for (i=0;i<numatoms;i++)
     *  fprintf(stderr,"fr1[%i]=%f  ",i,fr1[i]);
     *  for (i=0;i<numatoms*order;i++)
     *  printf("%f \n",theta1[i]); 
     */
    
  }
#endif
  
  if(last_step){ 
    /* if this is the last time step -> free */
#if DPME_DEBUG
    printf("\n FREE Arrays Q/ffwork/THETA/DTHETA NOW! \n");
#endif
    /* alloc via dvector return ptr+1 */
    free(theta1-1);free(theta2-1);free(theta3-1);
    free(dtheta1-1);free(dtheta2-1);free(dtheta3-1); 
    free(q-1);
    free(ffwork-1);
  }
  /* free every time step */
  
#if (TIMEME*VERBOSE)
  printf("~~~~~~~~~~~~~~~~~RCP SUM wall time[secs]~~~~~~~~~~~~~~~~~~~~~\n");
  printf("bspline-cofs2 =%f  fill_grd =%f, FFT-bk=%f  \n",swatch(intime[1],intime[2]),
	 swatch(intime[2],intime[3]),swatch(intime[3],intime[4]));
  printf("scalar_sum=%f, FFT-fwd=%f , grad_sum=%f\n",
	 swatch(intime[4],intime[5]),swatch(intime[5],intime[6]),
	 swatch(intime[6],intime[7]));
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
#endif
  return (eer); /* recip energy */
} /* end recip_sum_calc 2*/

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* this routine will calculate the recip_sum on each PE  */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

double dpme_eval_recip(AtomInfo atom_info, Pme2Particle *Myparticles, 
		      PmeVector **recipF, double *virial,GridInfo grid_info,
		      BoxInfo box_info, PeInfo pe_info,
		      int time_count, int tsteps, double *mytime)
{
  static int firsttime=0;
  static int last_step=0; /* if this is last time step */
  static double  *fr1, *fr2,  *fr3,*bsp_mod1,*bsp_mod2,*bsp_mod3;
  static double  *fftable;
  double recip_ene; /* recip pot */
  static PmeVector *myrecipF; 
  int  i,order, atompnt,nlocal,numatoms;
  double ewaldcof,*recip,volume;
  MYPROC myproc;
  Grid nfftgrd; 
#if TIMEME
  struct timeval time1,time2;
  struct timezone tzp;
  gettimeofday(&(time1),&tzp);
#endif
  
  /************************************/
  atom_info.nlocal =atom_info.numatoms;  /* i have all the atoms for recip sum */
  grid_info.volume= box_info.box.x * box_info.box.y * box_info.box.z; 
  atompnt = atom_info.atompnt;
  nlocal=  (atom_info.nlocal);
  numatoms = atom_info.numatoms;
  ewaldcof= box_info.ewaldcof;
  recip = box_info.recip;
  myproc= pe_info.myproc;
  volume = grid_info.volume;
  order= grid_info.order;
  nfftgrd= grid_info.nfftgrd;
  /*-----------------------------------------*/
#if VIRIAL
  --virial;
#endif

  if(!firsttime){

#if DPME_DEBUG
    printf(" allocating memory for bsp_mod, fr, fftable, recipF... \n");
#endif
    bsp_mod1=dvector(0,nfftgrd.x);/* keep don't free */
    bsp_mod2=dvector(0,nfftgrd.y);/* keep don't free */
    bsp_mod3=dvector(0,nfftgrd.z);/* keep don't free */ 
    /* keep don't free */
    fftable= dvector(0,3*(4*max(nfftgrd.z,(max(nfftgrd.x,nfftgrd.y)))+15));

    /* alloc force array firsttime only, later can realloc on renieghbor */
    if (!(myrecipF = (PmeVector *) malloc((numatoms+1) * sizeof(PmeVector))) )
      fprintf(stderr,"Error in allocating space for recipF  Data !!!\n");
    firsttime++;
  } /* if firsttime */
  
  (*recipF)=myrecipF;
  
  if (time_count==tsteps) last_step=1;


    recip_sum_setup2(numatoms,nlocal,
		     recip,Myparticles,nfftgrd,order,myproc,&fr1,
		     &fr2,&fr3, fftable,bsp_mod1,bsp_mod2,bsp_mod3);


    for (i=0;i<=numatoms;i++){
      myrecipF[i].x=0.;
      myrecipF[i].y=0.;
      myrecipF[i].z=0.;
    }

    recip_ene=recip_sum_calc2(ewaldcof,volume,recip,myproc, nlocal,
			      &fr1,&fr2,&fr3,
			      numatoms,order,nfftgrd,Myparticles,fftable,virial,
			      bsp_mod1,bsp_mod2,bsp_mod3,myrecipF,last_step);


#if VERBOSE
  printf("\n At t=%d..........recip_energy=%f\n",time_count,recip_ene);
#endif


  /* check if this is the last time-step before freeing */

  if (last_step) {
#if (DPME_DEBUG) 
      printf("freeing fr123, bsp_mod123,fftable...\n");
#endif
      free(fr1);
      free(fr2);
      free(fr3);
      free(bsp_mod1 -1); /* since dvector returns ptr +1 */
      free(bsp_mod2 -1);
      free(bsp_mod3 -1);
      free(fftable -1);
    }
#if TIMEME
  gettimeofday(&(time2),&tzp);
  *mytime=swatch(time1,time2);
#endif

  
  return(recip_ene);
}
/* end dpme_eval_recip() */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/





