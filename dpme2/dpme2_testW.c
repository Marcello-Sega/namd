/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
 /* $Id: dpme2_testW.c,v 1.1 1997/04/23 17:05:20 nealk Exp $ 
  */
/**********************************************************************
*
* New * New * New * New * New * New * New * New * New * New * New * New
* 
* This Code implements the Master/Worker paradigm , where the Master
* Does the recip sum only while the dir_sum is done by workers.
* This approach is targeted for cluster of WorkStations that  don't 
* have a parallel 3D FFT
* Note that the master is the last PE that joins the group and not the 
* first (or PE=0) as is traditionally done.
*/
/* TO DO :
 * 1. make nlist/nnlist and all other big array allocated statically to be 
 * dynamically allocated/free in their respective routines
 */
 /*************************************************************************
  * A. (Nour) Toukmaji - Duke University; ECE Dept. 1995
  * Acknowledgment: some subroutines in this file have been incorporated / modified
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
  
   /* Local variables */
   char *args[11];
   int  i;
   double  self_ene ; /* the self energy */
   double adj_recip_ene; /* the adjusted reciprocal energy */
      
   /* timing vars */
   struct tms startbuf,endbuf;
   struct timeval time1, time2;
   struct timezone tzp;
   double dir_time_tot=0.0;
   double rcp_time_tot=0.0;
   double time_used[10];
  
   /* mutliple time step variables */
   int time_count; /* time index for mts version */
   static int tsteps; /* desired number of time steps */

   /* vars for interfacing */
   Pme2Particle *ParticlePtr; /* contains (x,y,z,cg,id) per particle */
   PmeVector    *directF; /* contains x,y,z components of Direct-space sum */
   PmeVector    *recipF; /* contains x,y,z components of Recip-space sum */
   PmeVector    *adjustF; /* contains x,y,z components of recip adjust-forces  */
   PePnt        *pe_atompnt; /* master-list of PE's atom pointers to pe_atom[] */
   double       *pe_atoms; /* master-list of atoms that arrived from each worker */

   /* vars added for new method of spatial decompostion */
   int *list; /* linked list of owned atoms in a proc */
 
   double direct_energy,recip_energy; /* ewaldpot direct and recip */
   double adj_dir_ene; /* dir_sum energy correction for bonded interactions */
   double adj_vir[6], recip_vir[6], dir_vir[6]; /* virial arrays */
   int max_used; /* max storage to alloc for directF */
   
  
   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
   /* 3/97 new data structures defined in dpme2def.h*/

   AtomInfo atom_info; /* has main system info about atoms */
   BoxInfo box_info; /* has system box info */
   PeInfo pe_info; /* has processing elements  data */
   BndryInfo bndry_info; /* info about the spatial-decomp boxes/borders to PEs */
   SwapInfo swap_info; /* info about which and where atoms to be swapped */
   BinInfo bin_info; /* linkcell bins and their supporting lists */
   GridInfo grid_info; /* the reciprocal space grid info */
 /*##########################################################################*/
   
#if (DPME_DEBUG) 
   printf("pvm_catchout is on\n");
   pvm_catchout(stdout);  /* if def=1  show all slaves output on tasker screen */
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

    grid_info.nfft=0; /* defualt system is cube */

    i = 1; 
    while( i < argc ) {
      switch( argv[i][1] )
	{
	case 'c': box_info.cutoff= atof( argv[i]+2);
	  break;
	case 't': box_info.dtol= atof( argv[i]+2 );
	  break;
	  
	case 'o': grid_info.order = atoi( argv[i]+2);
	  break;
	case 's':  pe_info.myproc.nprocs =  atoi( argv[i]+2 ); 
	  break;
	case 'm':  tsteps = atoi( argv[i]+2 );
	  break;
	case 'x':  bndry_info.npdim[0] = atoi( argv[i]+2 );
	  break;
	case 'y':  bndry_info.npdim[1] = atoi( argv[i]+2 );
	  break;
	case 'z':  bndry_info.npdim[2] = atoi( argv[i]+2 );
	  break;
	case 'n': 
	  switch( argv[i][2] ) 
	    {
	    case 'x': grid_info.nfftgrd.x= atoi( argv[i]+3);
	      break;
	    case 'y':  grid_info.nfftgrd.y= atoi( argv[i]+3);
	      break;
	    case 'z':  grid_info.nfftgrd.z= atoi( argv[i]+3);
	      break;  
	    default:  grid_info.nfft= atoi( argv[i]+2);
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
   pe_info.igrid=0;
   if ( (argc >=10) ) pe_info.igrid=1; /* user specifies prcoessor grid */

   /* set grid dimensions assuming orthogonal */
   if (grid_info.nfft) {
     grid_info.nfftgrd.x=grid_info.nfft;
     grid_info.nfftgrd.y=grid_info.nfft;  
     grid_info.nfftgrd.z=grid_info.nfft;
   }
  
   /******************************************************************
    * All tasks should call dpme_register() to register with pvm
    * and with the work group "direct_solver" this includes
    * the master process
    */
   
   dpme_register("dpme2_test",args,&(pe_info));

   /************************* program timer******************************* */
   /* timer begins here to be compatible with dpme version 1 */
   gettimeofday(&(time1),&tzp);
   times(&startbuf);

   /* reads the data file and setup the parallel domain decomposistion */

   dpme_setup(&atom_info, &box_info, &(ParticlePtr),&swap_info,&bndry_info,
	      &bin_info,&pe_info,&time_used[0]);
 
   
   /******************MASTER MASTER MASTER MASTER MASTER MASTER******************/
   /*********************************************************************/
   /******************MASTER MASTER MASTER MASTER MASTER MASTER******************/
   /*********************************************************************/
   /******************MASTER MASTER MASTER MASTER MASTER MASTER******************/ 
   
   if (pe_info.myproc.node==pe_info.myproc.nprocs) {  /* IMP: Master is the last PE */   
    
     /*********************************************************************/
     /******************** start Multiple time step loop *******************/
     /*********************************************************************/
     for (time_count=1; time_count<=tsteps; time_count++) {  
       printf (" ****** Master at time-step %d ********\n",time_count);
           
       /* evaluate the recip sum (serially) */

       recip_energy=dpme_eval_recip(atom_info, ParticlePtr,&recipF,recip_vir,
				    grid_info,box_info,pe_info,
				    time_count,tsteps,&time_used[1]);
#if TIMEME
       rcp_time_tot += time_used[1];
#endif
       /***********************************************************************/
       /* adjust the recip sum energy and forces works for Water ONLY */
       /* adjust_recip call the one of DPME1 */
       
       dpme_adjust_recip(atom_info,box_info,ParticlePtr,pe_info,time_count,tsteps,
		    adj_vir, &adj_recip_ene,&adjustF,&time_used[2]);
       /***********************************************************************/
       /* now Master recvs the adjusted dir-force from slaves */

       dpme_rcv_dir(&pe_atompnt,&pe_atoms,&directF,pe_info.myproc, 
		      atom_info.numatoms,&time_used[3]); 
       
       /***********************************************************************/
       /* update coordinates ALL  particles for next time-step, 
	*  fakes EOM integration 
	*/
       
       update_coordinates(atom_info,ParticlePtr,box_info); 
       
       /***********************************************************************/
       /* Master packs and sends to slaves their new updated coordinates */
       
       dpme_snd_xyz (pe_atompnt, pe_atoms, ParticlePtr, pe_info, &time_used[4]); 
       /***********************************************************************/
       /* All PE's send thier energies/virials  to PE#0 to print only */
       time_used[7]=0.0; /* if virial not active */
       if ( (((time_count-1) % (UPDATE_TIME) )==0) ||  (time_count==tsteps) ) {

	 collect_ene_results(&direct_energy,&adj_dir_ene, &self_ene, &recip_energy, 
			     &adj_recip_ene,pe_info);
#if VIRIAL
       collect_virials(adj_vir,dir_vir,recip_vir,pe_info,&time_used[7]);
#endif
       }
       /****************************************************************************/
     }  /* end time count loop */
 
#if TIMEME*VERBOSE
     printf("============================================================\n");
     printf("====== Timings break down of the last time step ============\n");
     printf("============================================================\n");
     printf("====== dpme_setup= %f [sec] \n",time_used[0]);
     printf("====== dpme_eval_recip=%f [sec]\n",time_used[1]);
     printf("====== adjust_recip=%f  \n",time_used[2]);
     printf("====== snd_XYZ=%f   rcv_dir=%f   virial=%f   \n",
	    time_used[4],time_used[3],time_used[7]);
     printf("============================================================\n");
#endif
    

#if DPME_DEBUG
     /* reads in exact resutls and compares with this run,
      * attached file is only for 2k system 
      */
     verify_results(atom_info.numatoms,directF,recipF,adjustF);
     
     /* dump all forces , note that directF (direct-sum force is already adjusted) 
      * while the recip sum force (recipF) becomes adjusted when u add the adjustF
      * forces to it. The adjustment comes about from evaluating the bonded interactions
      */
     
     for (i=1;i<=atom_info.numatoms;i++) {
       printf("Atm.id=%d x (%6.3f,%6.3f,%6.3f) x (%6.3f,%6.3f,%6.3f) x (%6.3f,%6.3f,%6.3f) \n",
	      (int)ParticlePtr[i].id,directF[i].x,directF[i].y,directF[i].z,
	      adjustF[i].x,adjustF[i].y, adjustF[i].z,
	      recipF[i].x,recipF[i].y,recipF[i].z);
     }
#endif
   } /*END MASTER REGION */
   /* ******************************************************************** */
   /************************** MASTER END MTS ******************************/ 
   /* ******************************************************************** */
   


   /* ******************************************************************** */
   /****************** WORKER WORKER WORKER  WORKER WORKER WORKER ****************/
   /* ******************************************************************** */
   /****************** WORKER WORKER WORKER  WORKER WORKER WORKER ****************/
   /* ******************************************************************** */
   /****************** WORKER WORKER WORKER  WORKER WORKER WORKER ****************/
   /* ******************************************************************** */

   else { /* SLAVE REGION */
       
     /*********************************************************************/
     /******************** SLAVE start Multiple time step loop *******************/
     /*********************************************************************/
     for (time_count=1; time_count<=tsteps; time_count++) {  

       max_used = dpme_reneighbor(&atom_info,&list,&box_info,&ParticlePtr,&pe_info,
				  &bndry_info,&swap_info,&bin_info,time_count,
				  tsteps,time_used);
       /****************************************************************************/
       /* calc the direct_sum forces and energy for atoms I own only */
       
       direct_energy= dpme_dir_force(&atom_info,list,swap_info,ParticlePtr,
				     box_info,&directF,dir_vir,pe_info,
				     time_count,tsteps,&time_used[4],max_used);
#if TIMEME
       dir_time_tot += time_used[4];
#endif

       /***********************************************************************/
       /* adjust direct sum energy and forces in parallel works for Water ONLY */
       /* can choose to send the directForces w/o adjusting for bonded interactions*/
       
       dpme_adjust_dir(&atom_info,swap_info, ParticlePtr,box_info,directF,pe_info,
		  time_count,tsteps,dir_vir,list, &adj_dir_ene,&time_used[5]);
       
       /***********************************************************************/
       /****** slave should send back adjusted dir-force to master *******/
       
       dpme_snd_dir(atom_info,list,pe_info,ParticlePtr,directF,&time_used[8]);
       
       /***********************************************************************/  
       /* do the self energy term in parallel, report only my PE's share */
       self_ene= 
	 dpme_self(atom_info,ParticlePtr,list,box_info.ewaldcof,pe_info,&time_used[6]);
     
       /***************************************************************************/
       /**** worker rcv new updated set of coordinates from the master ************/
       
       dpme_rcv_xyz(atom_info,list,ParticlePtr,&time_used[9]); 
       /***********************************************************************/
       /* All PE's send thier energies/virial to PE#0 only if the application cares 
	* about these values, currently is called every time the list is redone
	* and at the last step
	*/   
       time_used[7]=0.0; /* if virial not active */
       if (  (((time_count-1) % (UPDATE_TIME) )==0) ||  (time_count==tsteps) ) {
	 collect_ene_results(&direct_energy,&adj_dir_ene, &self_ene, &recip_energy, 
			     &adj_recip_ene,pe_info);
#if VIRIAL
	 collect_virials(adj_vir,dir_vir,recip_vir,pe_info,&time_used[7]);
#endif
       }
       /****************************************************************************/
     }  /* end time count loop */
     
#if TIMEME*VERBOSE
     printf("============================================================\n");
     printf("====== Timings break down of the last time step ============\n");
     printf("============================================================\n");
     printf("====== dpme_setup= %f [sec] \n",time_used[0]);
     printf("====== exchange/communicate %f [sec]\n",time_used[1]);
     printf("====== border= %f      build_LC=%f  \n",time_used[2],time_used[3]);
     printf("====== adjust_dir=%f   self=%f   virial=%f   \n",
	    time_used[5],time_used[6],time_used[7]);
     printf("====== snd_dir=%f   rcv_XYZ=%f      \n",time_used[8],time_used[9]);
     printf("============================================================\n");
#endif
     /* ******************************************************************** */
     /********************************* END MTS ******************************/ 
     /* ******************************************************************** */
   } /* END SLAVE REGION */
   

   times(&endbuf);  /* prog time */
   gettimeofday(&(time2),&tzp);  /* prog time */
 /****************************************************************************/

   printf("*********  Results at last time step **********\n");
   printf("* Total Adj_direct_energy =          %f\n",direct_energy+adj_dir_ene);
   printf("* Total Recip_energy =                 %f\n",recip_energy);
   printf("* Total Adjust_energy =                %f\n",adj_recip_ene);
   printf("* Total Self_energy =                %f \n",self_ene);
   printf("* TOTAL ENE =                       %f\n",
	  direct_energy+adj_dir_ene+recip_energy+adj_recip_ene+self_ene);
   printf("***********************************************\n");
  
   /**************************************************************************/  
   /* I added this barrier to synchro results in the last time step */

   if( pvm_barrier(GROUP,(pe_info.myproc.nprocs+1))  <0) { /* synchro all group members */ 
     printf(" Error: at the last barrier \n");
     pvm_exit();
     exit(0);
   }
   /**************************************************************************/ 
#if TIMEME
  printf("Total Wall/User  time = %f / %f [secs]\n",
	 swatch(time1,time2),uswatch(startbuf,endbuf));
  printf("Total DIR_SUM walltime on 1 PE=%f [sec]\n",dir_time_tot);
  printf("Total RCP_SUM walltime on 1 PE=%f [sec]\n",rcp_time_tot);
#endif
/**************************************************************************/ 

#if DPME_DEBUG
   printf("Clean up: free Memory \n");
   printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
#endif

   free(directF); /* alloc'ed in dir_force()[slave] and master_rcv_dir[mstr] */
   free(ParticlePtr); /* alloc'ed in setup_atom()[slave] setup_general[mstr] */
   
   if (pe_info.myproc.node == pe_info.myproc.nprocs ){ /* master */ 
     free(recipF); /* alloc'ed in  eval_recip_sum() */
     free(adjustF);/* alloc'ed in adjust_dir_recip() */
     free(pe_atompnt); /* alloc'ed in master_rcv_dir() */
     free(pe_atoms); /* alloc'ed in master_rcv_dir() */
   }
   else { /* worker */
     free(bin_info.bin);  /* alloced in build_linkcell */
     free(bin_info.binpnt);/* alloced in build_linkcell */
     free(list); 
     free(swap_info.nlist);   /* alloc'd in dpme_reneighbor */
     free(swap_info.imglist); /* alloc'd in dpme_reneighbor */
   }
   pvm_exit(); 
   exit(0);
   
 } /* END MAIN */
/***************************************************************/
