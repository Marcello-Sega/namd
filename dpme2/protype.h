/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
/* $Id: protype.h,v 1.4 1999/03/18 22:06:32 jim Exp $
 */

/* prototype of functions and subroutines used in DPME2 */

/**************** dpme2_recip_sum.c *****************************/
int recip_sum_setup2(int numatoms,int nlocal,
		     double *recip, Pme2Particle *Myparticles,
		     Grid nfftgrd, int order, MYPROC myproc, double **fr1, double **fr2,
		     double **fr3, double *fftable, double *bsp_mod1, double *bsp_mod2,
		     double *bsp_mod3);

double recip_sum_calc2(double ewaldcof, double volume, double *recip, MYPROC myproc, 
		       int nlocal,
		       double **fr1, double **fr2, double **fr3, int numatoms,
		       int order, Grid nfftgrd, Pme2Particle *Myparticles, double *fftable,
		       double *virial, double *bsp_mod1, double *bsp_mod2, 
		       double *bsp_mod3,PmeVector *rfparticle, int last_step);

double dpme_eval_recip(AtomInfo atom_info, Pme2Particle *Myparticles, 
		      PmeVector **recipF, double *virial,GridInfo grid_info,
		      BoxInfo box_info, PeInfo pe_info,
		      int time_count, int tsteps, double *mytime);

/*...........................below only in dpme2_recip_sum2.c ..............................*/
int get_bspline_coeffs3(int atom_index, int *numatoms, double *fr1, 
			double *fr2, double *fr3, int *order, double *theta1, 
			double *theta2, double *theta3, double *dtheta1, 
			double *dtheta2, double *dtheta3);
int fill_charge_grid3(int atom_index,int *numatoms, int nlocal, Pme2Particle *ParticlePtr,
		      double *theta1, double *theta2, double *theta3, 
		      double *fr1, double *fr2, double *fr3,  int *order, 
		      int *nfft1,  int *nfft2,  int *nfft3,  int *nfftdim1, 
		      int *nfftdim2, int *nfftdim3, double *q, double *tmpcg,int clear_q_grid);

/********************* dpme2_utility.c *************************/
int find_ewaldcof(double *cutoff, double *dtol, 
		  double *ewaldcof);
int setup_general(AtomInfo *atom_info,double *cutoff, double *dtol, double *ewaldcof, 
		  PmeBVector *box, PmeVector *prd, PmeVector *prd2, 
		  Pme2Particle **ParticlePtr, PmeVector *mymc2);

int setup_parallel(AtomInfo *atom_info, BoxInfo *box_info, BndryInfo *bndry_info,
		   PeInfo *pe_info, SwapInfo *swap_info, BinInfo *bin_info);

int setup_atom( AtomInfo *atom_info,int **mylist, Pme2Particle **ParticlePtr, PmeBVector box, 
	       double *border, PmeVector prd2, MYPROC *myproc,
	       int *tidarray);

int exchange(AtomInfo *atom_info, Pme2Particle **AllParticles, double *border, int *list, 
	     int *npdim, int *mpart, MYPROC myproc,
	     int *tidarray,double *mytime);

int communicate(int nswap, int *nslist, int *slist, MYPROC myproc, Pme2Particle *Myparticles, 
		int *spart, int *rpart, int *tidarray,double *mytime);

int borders(AtomInfo *atom_info,int *nswap, double *boundlo, double *boundhi, 
	    Pme2Particle **AllParticles, int *list, MYPROC myproc, int *need, 
	    int *nslist, int *slist, int *spart, int *rpart,int *tidarray,double *mytime);

int build_verlet(AtomInfo *atom_info, int *list, Pme2Particle **AllParticles, int *nnlist, 
		 int *nlist, PmeVector myprd, double rs, double *mytime);

int build_linkcell(AtomInfo *atom_info , PmeVector mbin, PmeVector binsize, 
		   int *list, int *nlist, int *nnlist, short int *imglist, 
		   Pme2Particle *Myparticles, int **bin, 
		   int **binpnt, PmeVector nbin, PmeVector myprd, PmeVector myprd2, 
		   PmeVector mbinlo, double rs,double *mytime);

int error_handler(char *errstring);

int mesh_3d(int *nx, int *ny, int *nz, int *node, 
	      int *ix, int *iy, int *iz, int *iwest, int *ieast, 
	      int *isouth, int *inorth, int *idown, int *iup);

double dir_energy(int *atompnt, int *nlocal, int *nnlist, int *nlist, int *list, MYPROC myproc, 
		  PmeVector myprd, double cutoff, double ewaldcof, Pme2Particle *Myparticles, 
		  PmeVector mymc2, int *tidarray);
double dpme_dir_force(AtomInfo *atom_info, int *list, SwapInfo swap_info,
		 Pme2Particle *Myparticles, BoxInfo box_info, 
		 PmeVector **directF, double *virial, 
		 PeInfo pe_info,int time_count, int tsteps, double *mytime,
		 int max_used);


double *dvector( int nl, int nh);

double *dvector2( int nl, int nh);

void   update_coordinates(AtomInfo atom_info, Pme2Particle *particlelist,BoxInfo box_info);

double adjust_dir_recip(int *atompnt, int *nlocal, int *list, int *nlist, int *nnlist, 
			Pme2Particle *Myparticles, PmeVector myprd, PmeVector mymc2, 
			double cutoff, double ewaldcof, PmeVector *mydirectF, 
			double *recip, MYPROC myproc,int time_count, int tsteps, 
			int *tidarray, double *dir_vir,double *adj_vir, 
			double *my_adj_rcp_eng, double *my_adj_dir_ene, 
			PmeVector **adjustF,double *mytime,int max_used);

double dpme_adjust_recip(AtomInfo atom_info,BoxInfo box_info,Pme2Particle *Myparticles,
		    PeInfo pe_info, int time_count, int tsteps, double *adj_vir, 
		    double *my_adj_rcp_eng, PmeVector **adjustF,double *mytime);


double dpme_adjust_dir(AtomInfo *atom_info, SwapInfo swap_info,
		  Pme2Particle *Myparticles, BoxInfo box_info, 
		  PmeVector *mydirectF, PeInfo pe_info,
		  int time_count, int tsteps,
		  double *dir_vir,  int *list,
		  double *my_adj_dir_eng, double *mytime);

int correct_water_dir( int i, int j,Pme2Particle particle1, Pme2Particle particle2, 
		  double cutoff, double ewaldcof, PmeVector myprd, PmeVector mymc2, double *ene, 
		  PmeVector *mydirectF, double *recip,double *virial,float full_interact);

int correct_water_recip( int i,  int j, Pme2Particle particle1, Pme2Particle particle2,
			double ewaldcof, double *ene, PmeVector *afparticle, double *vir);
		
int ew_direct(double *ewaldcof, double *crd, double *pot, double *grad, double *r);
int ew_adjust(double *ewaldcof, double *crd, double *pot, double *grad);

double dpme_self(AtomInfo atom_info, Pme2Particle *particlelist,  int *list,
	 double ewaldcof, PeInfo pe_info, double *mytime);
int read_file(int numatoms,int offset,PmeVector *myforce);
int check_forces(int numatoms, PmeVector *ApprxF,PmeVector *ExactF, PmeVector *BaseF);
int verify_results(int nlocal,PmeVector *directF, PmeVector *recipF, PmeVector *adjustF);
double swatch(struct timeval st,struct timeval end);
double uswatch(struct tms startbuf, struct tms endbuf);

int dpme_setup(AtomInfo *atom_info, BoxInfo *box_info, Pme2Particle **ParticlePtr,
	       SwapInfo *swap_info, BndryInfo *bndry_info, BinInfo *bin_info, 
	       PeInfo *pe_info, double *mytime);
int dpme_reneighbor(AtomInfo *atom_info, int **list, BoxInfo *box_info, 
		    Pme2Particle **ParticlePtr,PeInfo *pe_info,BndryInfo *bndry_info,
		    SwapInfo *swap_info,BinInfo *bin_info, int time_count,int tsteps,
		    double *time_used);

/************************ dpme2_paralibs.c ************************/

void dpme_register(char *spawnme, char **args, PeInfo *pe_info);
int merge_i(int *ivalue, MYPROC *myproc, int *tidarray); 
int merge_i2(int *ivalue, MYPROC *myproc, int *tidarray);
int swap(int node, double *sbuf, int islen, int isnode, double *rbuf, int *irlen, 
	 int irnode,int *tidarray);
int merge_d(double *dvalue, MYPROC myproc, int *tidarray,int msg_tag);
int merge_d2(double *dvalue, MYPROC myproc, int *tidarray,int msg_tag );
int merge_v2(double *dvalue, MYPROC myproc, int *tidarray, int vec_length, int msg_tag);
/* added for the master /worker  version */

void dpme_snd_dir(AtomInfo atom_info, int *list,PeInfo pe_info,
		    Pme2Particle *Myparticles, PmeVector *mydirectF, double *mytime);
void dpme_rcv_dir(PePnt **pe_atompnt, double **pe_atoms,  PmeVector **directF,
		    MYPROC myproc, int numatoms, double *mytime);
void dpme_snd_xyz(PePnt *pe_atompnt, double *pe_atoms, Pme2Particle *Myparticles,
		    PeInfo pe_info, double *mytime) ;
void dpme_rcv_xyz(AtomInfo atom_info, int *list, Pme2Particle *Myparticles,
		    double *mytime ) ;
void collect_ene_results(double *dir_ene, double *adj_dir, double *self_ene,
			 double *rcp_ene, double *adj_rcp, PeInfo pe_info);
void collect_virials(double *adj_vir, double *dir_vir, double *rec_vir,
		     PeInfo pe_info, double *mytime);
/*****************  erfcfun.c ***********************************/
double  erfcfun(double *x);

/**************** Dpme_recip dir ***************************/
 int get_bspline_coeffs2( int *numatoms, double *fr1, 
	double *fr2, double *fr3,  int *order, double *theta1, 
	double *theta2, double *theta3, double *dtheta1, 
	double *dtheta2, double *dtheta3);

 int fill_bspline(double w,   int order, double 
	*array, double *darray);

 int init(double *c, double *x,  int *order);

 int one_pass(double *c, double *x,  int *k);

 int diff(double *c, double *d,  int *order);

int fill_charge_grid2(int numatoms, int nlocal,
		      Pme2Particle *ParticlePtr,
	 double *theta1, double *theta2, double *theta3, 
	double *fr1, double *fr2, double *fr3,  int order, 
	int nfft1,  int nfft2,  int nfft3,  int nfftdim1, 
	int nfftdim2, int nfftdim3, double *q);

 int clearq(double *q, int *ntot);


 int scalar_sum(double *q, double *ewaldcof, 
	double *volume, double *recip, double *bsp_mod1, 
	double *bsp_mod2, double *bsp_mod3, int *nfft1, 
	int *nfft2, int *nfft3, int *nfftdim1, int *nfftdim2, 
	int *nfftdim3, double *eer, double *vir);

 int grad_sum( int numatoms, 
	      Pme2Particle *ParticlePtr, 
	double *recip, double *theta1, double *theta2, double 
	*theta3, double *dtheta1, double *dtheta2, double *
	dtheta3, PmeVector *rfparticle, double *
	fr1, double *fr2, double *fr3,  int order,  int nfft1,
	  int nfft2,  int nfft3,  int nfftdim1,  int nfftdim2,
	  int nfftdim3, double *q);


 int get_fftdims( int *nfft1,  int *nfft2,  int *
	nfft3,  int *nfftdim1,  int *nfftdim2,  int *nfftdim3, 
	 int *nfftable,  int *nffwork,  int *sizfftab,  int *
	sizffwrk);

 int fft_setup(double *array, double *fftable, 
	double *ffwork,  int *nfft1,  int *nfft2,  int *nfft3, 
	 int *nfftdim1,  int *nfftdim2,  int *nfftdim3,  int *
	nfftable,  int *nffwork);


 int fft_forward(doublecomplex *array, double *fftable, 
	doublecomplex *ffwork,  int *nfft1,  int *nfft2,  int *nfft3, 
	 int *nfftdim1,  int *nfftdim2,  int *nfftdim3,  int *
	nfftable,  int *nffwork);

 int fft_back(doublecomplex *array, double *fftable, 
	 doublecomplex *ffwork,  int *nfft1,  int *nfft2,  int *nfft3, 
	 int *nfftdim1,  int *nfftdim2,  int *nfftdim3,  int *
	nfftable,  int *nffwork);

int pmesh_kspace_get_sizes( int *nfft1,  int *nfft2, 
     int *nfft3,  int *numatoms,  int *order,  int *sizfftab, 
     int *sizffwrk,  int *siztheta,  int *siz_q,  int *
    sizheap,  int *sizstack);

 int pmesh_kspace_setup(double *bsp_mod1, double *
	bsp_mod2, double *bsp_mod3, double *fftable, double *
	ffwork,  int *nfft1,  int *nfft2,  int *nfft3,  int *
	order,  int *sizfftab,  int *sizffwrk);


 int get_scaled_fractionals2( int *numatoms,
	Pme2Particle *ParticlePtr, double *recip,  int *nfft1, 
	 int *nfft2,  int *nfft3, double *fr1, double *fr2, 
	double *fr3);

 int load_bsp_moduli(double *bsp_mod1, double *
	bsp_mod2, double *bsp_mod3,  int *nfft1,  int *nfft2, 
	 int *nfft3,  int *order);

 int dftmod(double *bsp_mod, double *bsp_arr, 
	 int *nfft);
double *dvector_pme(int nh);
/*****************************************/
