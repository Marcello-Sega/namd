/* 
 *  Abdulnour Toukmaji 
 *  Copyright (c) 1996,1997 Duke University
 *  All rights reserved
 */
/* $Id: protype.h,v 1.1 1997/03/14 15:07:20 nealk Exp $
 */

/* prototype of functions and subroutines used in DPME2 */

#include "dpme2def.h"

/**************** dpme2_recip_sum.c *****************************/
int recip_sum_setup2(int numatoms,int nlocal,int atompnt, int *list, 
		     double *recip, Pme2Particle *Myparticles,
		     Grid nfftgrd, int order, MYPROC myproc, double **fr1, double **fr2,
		     double **fr3, double *fftable, double *bsp_mod1, double *bsp_mod2,
		     double *bsp_mod3);

double recip_sum_calc2(double ewaldcof, double volume, double *recip, MYPROC myproc, 
		       int nlocal,int atompnt, int *list, 
		       double **fr1, double **fr2, double **fr3, int numatoms,
		       int order, Grid nfftgrd, Pme2Particle *Myparticles, double *fftable,
		       double *virial, double *bsp_mod1, double *bsp_mod2, 
		       double *bsp_mod3,PmeVector *rfparticle, int last_step);

double eval_recip_sum(Pme2Particle *Myparticles, double *recip, PmeVector **recipF, 
		      double *virial, Grid nfftgrd, int numatoms, int order, double volume, 
		      double ewaldcof, MYPROC myproc, int nlocal, int atompnt, int *list,
		      int op_flag,int time_count, int tsteps,double *mytime);
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

int setup_general(double *cutoff, double *dtol, double *ewaldcof, int *numatoms,
		  PmeBVector *box, PmeVector *prd, PmeVector *prd2, int tsteps, 
		  Pme2Particle **ParticlePtr,PmeVector *mymc2);

int setup_parallel(int numatoms, int *npdim, double rs, PmeVector prd, PmeVector prd2, 
		   int *need, double *border, MYPROC *myproc, int *mpart, int *nswap, 
		   PmeVector *nbin, PmeVector *binsize, PmeVector *mbin, 
		   PmeVector *mbinlo, double *boundlo, double *boundhi, int *spart, 
		   int *rpart,int igrid, int *tidarray);

int setup_atom(int **list, Pme2Particle **ParticlePtr, PmeBVector box, int numatoms, double *border, PmeVector prd2, MYPROC *myproc, int *nlocal, int *freepnt, int *atompnt, int *tidarray);

int exchange(int *nlocal, Pme2Particle **AllParticles, double *border, int *list, int *npdim, 
	     int *freepnt, int *atompnt, int *mpart, MYPROC myproc, int *tidarray,double *mytime);

int communicate(int nswap, int *nslist, int *slist, MYPROC myproc, Pme2Particle *Myparticles, 
		int *spart, int *rpart, int *tidarray,double *mytime);

int borders(int *nswap, double *boundlo, double *boundhi, int *atompnt, int *nlocal, int *nother, 
	    Pme2Particle **AllParticles, int *list, MYPROC myproc, int *need, int *nslist, 
	    int *slist, int *spart, int *rpart, int *tidarray,double *mytime);

int build_verlet(int *nlocal, int *nother, int *list, Pme2Particle **AllParticles, int *nnlist, 
		 int *nlist, PmeVector myprd, double rs, int *atompnt,double *mytime);

int build_linkcell(int *atompnt, int *nlocal, int *nother, PmeVector mbin, PmeVector binsize, 
		   int *list, int *nlist, int *nnlist, Pme2Particle *Myparticles, int **bin, 
		   int **binpnt, PmeVector nbin, PmeVector myprd, PmeVector myprd2, 
		   PmeVector mbinlo, double rs,double *mytime);

int error_handler(char *errstring);

int mesh_3d(int *nx, int *ny, int *nz, int *node, 
	      int *ix, int *iy, int *iz, int *iwest, int *ieast, 
	      int *isouth, int *inorth, int *idown, int *iup);

double dir_energy(int *atompnt, int *nlocal, int *nnlist, int *nlist, int *list, MYPROC myproc, 
		  PmeVector myprd, double cutoff, double ewaldcof, Pme2Particle *Myparticles, 
		  PmeVector mymc2, int *tidarray);


double dir_force(int *atompnt, int *nlocal, int *list, int *nlist, int *nnlist, 
		 Pme2Particle *Myparticles, PmeVector myprd, PmeVector mymc2, 
		 double cutoff, double ewaldcof, PmeVector **directF, double *virial,
		 MYPROC myproc, int time_count, int tsteps, int *tidarray,double *mytime);

double *dvector( int nl, int nh);

double *dvector2( int nl, int nh);

void  update_coordinates(int *atompnt, int nlocal, int *list, Pme2Particle *particlelist,
			 PmeVector prd, PmeVector prd2);

double adjust_dir_recip(int *atompnt, int *nlocal, int *list, int *nlist, int *nnlist, 
			Pme2Particle *Myparticles, PmeVector myprd, PmeVector mymc2, 
			double cutoff, double ewaldcof, PmeVector *mydirectF, 
			double *recip, MYPROC myproc,int time_count, int tsteps, 
			int *tidarray, double *dir_vir,double *adj_vir, 
			double *my_adj_rcp_eng, double *my_adj_dir_ene, 
			PmeVector **adjustF,double *mytime);

int correct_water_dir( int i, int j,Pme2Particle particle1, Pme2Particle particle2, 
		  double cutoff, double ewaldcof, PmeVector myprd, PmeVector mymc2, double *ene, 
		  PmeVector *mydirectF, double *recip,double *virial,float full_interact);
int correct_water_recip( int i,  int j, Pme2Particle particle1, Pme2Particle particle2,
			double ewaldcof, double *ene, PmeVector *afparticle, double *vir,
			float full_interact);
int ew_direct(double *ewaldcof, double *crd, double *pot, double *grad, double *r);
int ew_adjust(double *ewaldcof, double *crd, double *pot, double *grad);
double self(Pme2Particle *particlelist,  int numatoms, int atompnt, int *list,
	    double ewaldcof, MYPROC myproc,int *tidarray,double *mytime);
void collect_virials(double *adj_vir, double *dir_vir, double *rec_vir,
		     MYPROC myproc, int *tidarray,double *mytime);
int read_file(int numatoms,int offset,PmeVector *myforce);
int check_forces(int numatoms, int atompnt, int *list, PmeVector *ApprxF, 
		 PmeVector *ExactF, PmeVector *BaseF);
int verify_results(int nlocal, int atompnt, int *list,
		   PmeVector *directF, PmeVector *recipF, PmeVector *adjustF);
double swatch(struct timeval st,struct timeval end);
double uswatch(struct tms startbuf, struct tms endbuf);
/************************ dpme2_paralibs.c ************************/

void dpme_register(char *spawnme, int *mygroupid, int nproc, int *inst_node, 
		   char **args, int *ncube);
int merge_i(int *ivalue, MYPROC *myproc, int *tidarray); 
int merge_i2(int *ivalue, MYPROC *myproc, int *tidarray);
int swap(int node, double *sbuf, int islen, int isnode, double *rbuf, int *irlen, 
	 int irnode,int *tidarray);
int merge_d(double *dvalue, MYPROC myproc, int *tidarray,int msg_tag);
int merge_d2(double *dvalue, MYPROC myproc, int *tidarray,int msg_tag );
int merge_v2(double *dvalue, MYPROC myproc, int *tidarray, int vec_length, int msg_tag);
/*****************  erfcfun.c ***********************************/
double  erfcfun(double *x);

/**************** Dpme_recip dir ***************************/
 int get_bspline_coeffs2( int *numatoms, double *fr1, 
	double *fr2, double *fr3,  int *order, double *theta1, 
	double *theta2, double *theta3, double *dtheta1, 
	double *dtheta2, double *dtheta3);

 int fill_bspline(double *w,   int *order, double 
	*array, double *darray);

 int init(double *c, double *x,  int *order);

 int one_pass(double *c, double *x,  int *k);

 int diff(double *c, double *d,  int *order);

int fill_charge_grid2(int *numatoms, int nlocal,int atompnt,int *list, 
		      Pme2Particle *ParticlePtr,
	 double *theta1, double *theta2, double *theta3, 
	double *fr1, double *fr2, double *fr3,  int *order, 
	int *nfft1,  int *nfft2,  int *nfft3,  int *nfftdim1, 
	int *nfftdim2, int *nfftdim3, double *q, double *tmpcg);

 int clearq(double *q, int *ntot);


 int scalar_sum(double *q, double *ewaldcof, 
	double *volume, double *recip, double *bsp_mod1, 
	double *bsp_mod2, double *bsp_mod3, int *nfft1, 
	int *nfft2, int *nfft3, int *nfftdim1, int *nfftdim2, 
	int *nfftdim3, double *eer, double *vir);

 int grad_sum( int *numatoms, int atompnt, int *list, 
	      Pme2Particle *ParticlePtr, 
	double *recip, double *theta1, double *theta2, double 
	*theta3, double *dtheta1, double *dtheta2, double *
	dtheta3, PmeVector *rfparticle, double *
	fr1, double *fr2, double *fr3,  int *order,  int *nfft1,
	  int *nfft2,  int *nfft3,  int *nfftdim1,  int *nfftdim2,
	  int *nfftdim3, double *q);


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


 int get_scaled_fractionals2( int *numatoms,int *atompnt, int *list, 
	double *tmpcg, Pme2Particle *ParticlePtr, double *recip,  int *nfft1, 
	 int *nfft2,  int *nfft3, double *fr1, double *fr2, 
	double *fr3);

 int load_bsp_moduli(double *bsp_mod1, double *
	bsp_mod2, double *bsp_mod3,  int *nfft1,  int *nfft2, 
	 int *nfft3,  int *order);

 int dftmod(double *bsp_mod, double *bsp_arr, 
	 int *nfft);
double *dvector_pme(int nh);
/*****************************************/
