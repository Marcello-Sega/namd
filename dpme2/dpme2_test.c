/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
 /* $Id: dpme2_test.c,v 1.1 1997/03/14 15:07:04 nealk Exp $ 
  */
 /*************************************************************************
  * A. (Nour) Toukmaji - Duke University; ECE Dept. 1995
  * Acknowledgment: the subroutines in this file have been incorporated / modified
  * from original PME source code by Tom Darden NIEHS, N.C., 
  * and Steve Plimpton (Sandia National Labs).
  * This release of DPME2 comprises the distributed version of PME method.
  * Currently, the direct_sum is distributed, while the recip_sum
  * is redundanlty performed by each processor due to serial 3D-FFT.

  * This driver program  calls up the major subroutines in DPME2 
  * and evaluates the Ewald sum for a system of water.
  * This version runs the Ewald sum for m-timesteps updating the 
  * neighbor-list data  every (UPDATE_TIME) timesteps.

  * NOTE: There is no integrator in this code and thus particle locations are
  * updated simply by modifying the particle coordinates by a small perturbation.
  *************************************************************************/
 /* Important notes ******************************************************
 *i- All particles locations are scaled from L/2 to -L/2 
 *ii- made the id in pmeparticle structure a double (not int) to fascilitate
     pvm-pk of recip_sum_setup
 *iii- At every timestep you either do 1-exchange()/border()/Nbr() then bcast fr123 
     or 2- bcast fr123 and then communicate(). This way u always bcast an updated
     version of coordinates.
 ************************************************************************/

 #include "dpme2.h"


 void main (int argc, char **argv)
 {

   struct tms startbuf,endbuf;
   struct timeval time1, time2;
   struct timezone tzp;
   double time_used[8];

   /* virial arrays */
   double adj_vir[6], recip_vir[6], dir_vir[6];

   /* Local variables */
   char *args[11];
   int nfft;
   int order;
   Grid nfftgrd;
   int numatoms, i;
   double  self_ene, dtol, ewaldcof, recip[9] /* was [3][3] */;
   double volume, adj_recip_ene, cutoff,min_box;
   PmeBVector box; /* unit cell dimensions in x y z, origin */
  
   /* parallel variables */
   int nslaves; 
   int inst_node[MAXPROC];

   /* mutliple time step variables */
   int time_count; /* time index for mts version */
   static int tsteps; /* desired number of time steps */

   /* vars for interfacing */
   Pme2Particle *ParticlePtr; /* contains (x,y,z,cg,id) per particle */
   PmeVector    *directF; /* contains x,y,z components of Direct-space sum */
   PmeVector    *recipF; /* contains x,y,z components of Recip-space sum */
   PmeVector    *adjustF; /* contains x,y,z components of recip adjust-forces  */

   /* vars added for new method of spatial decompostion */
   int npdim[3]; /* # of processors in each dim */
   int need[3]; /* # boxes away in each dimension for swapping */
   double border[6]; /* boundaries of processors box in each Dim */
   int mpart[6]; /* node# of nbr proc in each dim */
   PmeVector prd,prd2; /* contains xprd,yprd,zprd = box dimensions(half) */
   PmeVector mc2;  /* box min. inner cutoff */
   MYPROC myproc; /* info about parallel config */
   int nswap; /* number of swaps in each direction */
   PmeVector nbin, binsize; /*  bin size in  each dimension */
   PmeVector mbin, mbinlo; /* # bins in my box and 1st bin addrs */
   double boundlo[nsmax],boundhi[nsmax]; /* atom-position boundaries */
   int spart[nsmax], rpart[nsmax]; /* node's to send/rcv data from */
   double rs; /* this is the shell cutoff in multi-time step */
   int *list; /* linked list of owned atoms in a proc */
   int slist[nemax]; /* the swap list of atoms to send out */
   int nslist[nsmax+2]; /* pointer to beginnig of swap list per swap */
   int nlist[npmax*nnmax + nnmax]; /* nbr list of my atoms */
   int nnlist[npmax+2]; /* pointers to start of nbr list of my atoms */
   int nlocal; /* num of atoms in a processor */
   int nother; /* number of nearby atoms that My node stores */
   int freepnt, atompnt; /* index to 1st freespace , atom in list */
   int ineigh; /* specify which nbr method to use Verlet/LC*/
   double direct_energy,recip_energy; /* ewaldpot direct and recip */
   double adj_dir_ene; /* dir_sum energy correction for bonded interactions */
   int *bin, *binpnt; /* bins of atoms in Link_cell method */
   int igrid; /* flag indicates if user will specify processor grid */
 

 /*##########################################################################*/

 #if (DPME_DEBUG)
   pvm_catchout(stdout);  /* if def=1  show all slaves output on this screen */
 #endif

   for (i=0;i<argc-1;i++) /* copy the args to be passed to spawned procs */
     args[i]= argv[i+1]; 
   args[argc-1]=  (char*)0;  /* null char to end */

  if( argc < 6) {
    fprintf(stderr,"*****************************************************\n");
    fprintf(stderr,"usage:  %s  [-c#  -t# -n#[-nx# -ny# -nz#] -o#  -m# -s# -x# -y# -z#]\n",\
	    argv[0]); 
    fprintf(stderr,"  -c  direct space cutoff, -t direct space tolerance\n");
    fprintf(stderr,"  -n No. grid pts (cube) OR -nx -ny -nz (rectangle) (in pwr 2,3,5)\n");
    fprintf(stderr,"  -o order of interpolation\n");
    fprintf(stderr,"  -s No. of slaves to run, excluding the master process\n"); 
    fprintf(stderr,"  -m  number of timesteps \n"); 
    fprintf(stderr,"  -x,y,z the processors grid (if desired) \n");
    pvm_exit();
    exit(0);
  }

  else { /* command line arguments */

    nfft=0; /* defualt system is cube */

    i = 1; 
    while( i < argc ) {
      switch( argv[i][1] )
	{
	case 'c': cutoff= atof( argv[i]+2);
	  break;
	case 't': dtol= atof( argv[i]+2 );
	  break;
	  
	case 'o': order = atoi( argv[i]+2);
	  break;
	case 's': nslaves =  atoi( argv[i]+2 ); 
	  break;
	case 'm':  tsteps = atoi( argv[i]+2 );
	  break;
	case 'x':  npdim[0] = atoi( argv[i]+2 );
	  break;
	case 'y':  npdim[1] = atoi( argv[i]+2 );
	  break;
	case 'z':  npdim[2] = atoi( argv[i]+2 );
	  break;
	case 'n': 
	  switch( argv[i][2] ) 
	    {
	    case 'x': nfftgrd.x= atoi( argv[i]+3);
	      break;
	    case 'y': nfftgrd.y= atoi( argv[i]+3);
	      break;
	    case 'z': nfftgrd.z= atoi( argv[i]+3);
	      break;  
	    default: nfft= atoi( argv[i]+2);
	      break;
	    }
	  break;
	default: 
	  {   
	    fprintf(stderr,"Unknown command line argument: %s\n",argv[i]);
	    exit(-1);
	  }
	}
      i++;
    }
  } /* else */  
   igrid=0;
   if ( (argc >=10) ) igrid=1; /* user specifies prcoessor grid */
#if DPME_DEBUG
   printf("argc=%d, igrid=%d, nfft=%d\n",argc,igrid,nfft);
#endif
 /******************************************************************
 * All tasks should call dpme_register() to register with pvm
 * and with the work group "direct_solver" this includes
 * the master process
 */
   myproc.nprocs = nslaves +1 ;

   dpme_register("dpme2_test", &(myproc.node), myproc.nprocs, inst_node,args,
		 &myproc.ncube);

   /* set grid dimensions assuming orthogonal */
   if (nfft) {
     nfftgrd.x=nfft;
     nfftgrd.y=nfft;  
     nfftgrd.z=nfft;
   }
   
   /* program timer */
   /* timer begins here to be compatible with dpme version 1 */
   gettimeofday(&(time1),&tzp);
   times(&startbuf);

   /* find ewaldcof, read system coordinates, assign charge */
   setup_general(&cutoff, &dtol, &ewaldcof, &numatoms, &box,&prd,&prd2,tsteps,
		 &(ParticlePtr),&mc2);

   rs = cutoff; /* rs must be bounded as : (cutoff < rs < L/2), set to cutoff for now */

   /* find minimum box dim */
   min_box= ((box.x < box.y)? box.x:box.y);
   min_box = (( min_box < box.z) ? min_box: box.z); 
   if (rs>(min_box/2.)) error_handler("rs should be less than L/2");
   /***************************************************************/
   /* some global setup for running the water box example */
   for (i = 0; i <= 8; i++) {
     recip[i] = 0.;
   }
   recip[0]=1.0/box.x;
   recip[4]=1.0/box.y;
   recip[8]=1.0/box.z;

   volume= box.x * box.y * box.z; /* assume an orthogonal box */

   /* configure the geometrical boundaries for each PE, 
    * define its North/East/South/West/Up/Down neighbor PEs
    */
   ineigh=setup_parallel(numatoms, npdim,rs,prd,prd2,need,border,&myproc,mpart,&nswap,
		  &nbin,&binsize, &mbin,&mbinlo,boundlo,boundhi,spart,rpart,igrid,inst_node);
 
   /* read the config file small2.pdb and keep those particles within my boundary */
   setup_atom(&list,&ParticlePtr,box,numatoms,border,prd2,&myproc,&nlocal,&freepnt,
	      &atompnt,inst_node); 
   
   /*********************************************************************/
  /******************** start Multiple time step loop *******************/
  /*********************************************************************/
  for (time_count=1; time_count<=tsteps; time_count++) {  


    if (((time_count-1) % (UPDATE_TIME) )!=0){
      printf (" ****** COMMUNICATE at time-step %d ********\n",time_count);

      /* call the eval_recip with flag=1 so that broadcast is sent        */
      /* don't compute recip-ene now , do it with another call with flag=2 */
      eval_recip_sum(ParticlePtr,recip,&recipF,recip_vir,nfftgrd,numatoms,
		     order,volume,ewaldcof,myproc,nlocal,atompnt,list,1,
		     time_count,tsteps,&time_used[0]);
#if TIMEME
       rcp_time_tot += time_used[0];
#endif
      /* call communicate between reneighobring steps: exchange border atoms */
      communicate(nswap,nslist,slist,myproc,ParticlePtr,spart,rpart,inst_node,&time_used[1]); 
    }


    else {
      
      printf(" ****** RE-Nieghboring at time-step %d  ******\n",time_count);
      
      /* exchange atoms leaving / entering my box */
      exchange(&nlocal,&ParticlePtr,border,list,npdim,&freepnt,&atompnt,mpart,
	       myproc,inst_node,&time_used[1]);

      /* make up list of border atoms and exchange them */
      borders(&nswap,boundlo,boundhi,&atompnt,&nlocal,&nother,&ParticlePtr,list,myproc,
	      need,nslist,slist,spart,rpart,inst_node,&time_used[2]);

      /* make Particle list of my atoms using Verlet-like  or LinkCell-like method  */
      if (ineigh == 0) {
	build_verlet(&nlocal,&nother,list,&ParticlePtr,nnlist,nlist,prd,rs,&atompnt,&time_used[3]);
      }
      else {
	build_linkcell(&atompnt,&nlocal,&nother,mbin,binsize,list,nlist,nnlist,
		  ParticlePtr,&bin,&binpnt,nbin,prd,prd2,mbinlo,rs,&time_used[3]); 
      }

      /* call the eval_recip with flag=1 so that broadcast is sent        */
      /* don't compute recip-ene now , do it with another call with flag=2 */
      eval_recip_sum(ParticlePtr,recip,&recipF,recip_vir,nfftgrd,numatoms,
		     order,volume,ewaldcof,myproc,nlocal,atompnt,list,1,
		     time_count,tsteps,&time_used[0]);
#if TIMEME
      rcp_time_tot += time_used[0];
#endif
    } /* end else */
/****************************************************************************/
    /* calc the direct_sum forces and energy for atoms I own only */
    direct_energy= dir_force(&atompnt,&nlocal,list,nlist,nnlist,ParticlePtr,
			     prd,mc2,cutoff,ewaldcof,&directF,dir_vir,myproc,
			     time_count,tsteps,inst_node,&time_used[4]);
#if TIMEME
    dir_time_tot += time_used[4];
#endif
/***********************************************************************/
/* after you exchanged  the fractional coordinates now do the recip sum on each PE */
/* calc recip_sum; call w/flag=2 to calc,after u have called it earlier w/flag=1 */
    recip_energy=eval_recip_sum(ParticlePtr,recip,&recipF,recip_vir,nfftgrd,
				numatoms,order,volume,ewaldcof,myproc,nlocal,atompnt,
				list,2,time_count,tsteps,&time_used[0]);
#if TIMEME
    rcp_time_tot +=time_used[0];
#endif
/***************************************************************************/
/* adjust the direct and recip sum energy and forces in parallel works for Water ONLY */

    adjust_dir_recip(&atompnt,&nlocal,list,nlist,nnlist,ParticlePtr,
		     prd, mc2,cutoff,ewaldcof,directF,recip,myproc,
		     time_count,tsteps,inst_node,dir_vir,adj_vir,
		     &adj_recip_ene,&adj_dir_ene,&adjustF,&time_used[5]);
/***************************************************************************/   
/* do the self energy term in parallel, report only my PE's share */
    self_ene= self(ParticlePtr,nlocal,atompnt,list,ewaldcof,myproc,inst_node,&time_used[6]);
/***************************************************************************/
/* get the virials : direct, recip, adjust from all PEs */
#if VIRIAL
    collect_virials(adj_vir,dir_vir,recip_vir,myproc,inst_node,&time_used[7]);
#else 
    time_used[7]=0.0; /* if virial not active */
#endif
/***************************************************************************/
  /* update coordinates ONLY your particles for next time-step, fake EOM integration */
    update_coordinates(&atompnt,nlocal,list,ParticlePtr,prd,prd2);
/***************************************************************************/

  }  /* end time count loop */

   times(&endbuf);  /* prog time */
   gettimeofday(&(time2),&tzp);  /* prog time */
   /* ******************************************************************** */
   /********************************* END MTS ******************************/ 
   /* ******************************************************************** */
 
   printf("*********  Results at last time step **********\n");
   printf("* Total Adj_direct_energy =          %f\n",direct_energy+adj_dir_ene);
   printf("* Total Recip_energy =                 %f\n",recip_energy);
   printf("* Total Adjust_energy =                %f\n",adj_recip_ene);
   printf("* Total Self_energy =                %f \n",self_ene);
   printf("***********************************************\n");
  
  /**************************************************************************/  
 /* I added this barrier because to synchro results in the last time step */

   pvm_barrier(GROUP,-1);/* synchro all group members */ 
/**************************************************************************/ 

/* PRINT THE FORCES, COORDINATES: use for verification */
 
#if DPME_DEBUG2 
/* print the atom number and its xyz force components for the direct/recip sums
 * and the adjust force.  */
   i=atompnt;
   for (ii=1;ii<=nlocal;ii++){
     /* fprintf(stdout," Dir%i(%d)(%f,%f,%f)\n", i, (int)ParticlePtr[i].id, 
	directF[i].x, directF[i].y, directF[i].z);
	fprintf(stdout," Rcp%i(%d)(%f,%f,%f)\n", i, (int)ParticlePtr[i].id, 
	recipF[i].x, recipF[i].y, recipF[i].z);
	fprintf(stdout," Adj%i(%d)(%f,%f,%f)\n", i, (int)ParticlePtr[i].id, 
	adjustF[i].x, adjustF[i].y, adjustF[i].z); 
      */
     i=list[i];
   }

/* print the coordinates this PE  has ( for all coordinates let slaves=0) */
   i=atompnt;
   for (ii=1;ii<=nlocal;ii++){
     fprintf(stdout,
  	     "Atom %d = %2.3f\n%2.3f\n%2.3f\n",(int)ParticlePtr[i].id,
	     ParticlePtr[i].x+prd2.x+box.origin, 
	     ParticlePtr[i].y+prd2.y+box.origin, 
	     ParticlePtr[i].z+prd2.z+box.origin);
     i=list[i];
   }
#endif


#if (TIMEME*VERBOSE)
   printf("===================Last step breakdown timings=======================\n");
   printf("Wall time in communicate or exchange= %f \n",time_used[1]);
   printf("Wall time in borders = %f \n",time_used[2]);
   printf("Wall time in build_verlet/lc = %f \n",time_used[3]);
   printf("Wall time in adj/slf/vir = %f \n",time_used[5]+time_used[6]+time_used[7]);
   printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
#if TIMEME
  printf("Total Wall/User  time = %f / %f [secs]\n",
	 swatch(time1,time2),uswatch(startbuf,endbuf));
  printf("Total DIR_SUM walltime on 1 PE=%f [sec]\n",dir_time_tot);
  printf("Total RCP_SUM walltime on 1 PE=%f [sec]\n",rcp_time_tot);
#endif

#if DPME_DEBUG
   printf("Clean up: free Memory \n");
   printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif
/**************************************************************************/ 
  free(directF); /* alloc'ed in dir_force() */
  free(recipF); /* alloc'ed in  eval_recip_sum() */
  free(adjustF);/* alloc'ed in adjust_dir_recip() */
  free(ParticlePtr); /* alloc'ed in setup_atom() */
  free(bin);  /* alloced in build_linkcell */
  free(binpnt);/* alloced in build_linkcell */
  free(list); 

  pvm_exit(); 
  exit(0);

} /* END MAIN */
/***************************************************************/
