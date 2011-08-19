
#include "ComputeNonbondedCUDAKernel.h"
#include <stdio.h>

#ifdef NAMD_CUDA


__constant__ unsigned int const_exclusions[MAX_CONST_EXCLUSIONS];

static unsigned int *overflow_exclusions;

#define SET_EXCL(EXCL,BASE,DIFF) \
         (EXCL)[((BASE)+(DIFF))>>5] |= (1<<(((BASE)+(DIFF))&31))

void cuda_bind_exclusions(const unsigned int *t, int n) {

  cudaMalloc((void**) &overflow_exclusions, n*sizeof(unsigned int));
  cuda_errcheck("malloc overflow_exclusions");
  cudaMemcpy(overflow_exclusions, t,
		n*sizeof(unsigned int), cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to overflow_exclusions");
  int nconst = ( n < MAX_CONST_EXCLUSIONS ? n : MAX_CONST_EXCLUSIONS );
  cudaMemcpyToSymbol(const_exclusions, t, nconst*sizeof(unsigned int), 0);
  cuda_errcheck("memcpy to const_exclusions");
}


texture<float2, 1, cudaReadModeElementType> lj_table;
int lj_table_size;

void cuda_bind_lj_table(const float2 *t, int _lj_table_size) {
    static float2 *ct;
    static int lj_table_alloc;
    lj_table_size = _lj_table_size;
    if ( ct && lj_table_alloc < lj_table_size ) {
      cudaFree(ct);
      cuda_errcheck("freeing lj table");
      ct = 0;
    }
    if ( ! ct ) {
      lj_table_alloc = lj_table_size;
      cudaMalloc((void**) &ct, lj_table_size*lj_table_size*sizeof(float2));
      cuda_errcheck("allocating lj table");
    }
    cudaMemcpy(ct, t, lj_table_size*lj_table_size*sizeof(float2),
                                            cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to lj table");

    lj_table.normalized = false;
    lj_table.addressMode[0] = cudaAddressModeClamp;
    lj_table.filterMode = cudaFilterModePoint;

    cudaBindTexture((size_t*)0, lj_table, ct,
        lj_table_size*lj_table_size*sizeof(float2));
    cuda_errcheck("binding lj table to texture");
}


texture<float4, 1, cudaReadModeElementType> force_table;
texture<float4, 1, cudaReadModeElementType> energy_table;

void cuda_bind_force_table(const float4 *t, const float4 *et) {
    static cudaArray *ct;
    static cudaArray *ect;
    if ( ! ct ) {
      cudaMallocArray(&ct, &force_table.channelDesc, FORCE_TABLE_SIZE, 1);
      cuda_errcheck("allocating force table");
    }
    if ( ! ect ) {
      cudaMallocArray(&ect, &energy_table.channelDesc, FORCE_TABLE_SIZE, 1);
      cuda_errcheck("allocating energy table");
    }
    cudaMemcpyToArray(ct, 0, 0, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    // cudaMemcpy(ct, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to force table");
    cudaMemcpyToArray(ect, 0, 0, et, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to energy table");

    force_table.normalized = true;
    force_table.addressMode[0] = cudaAddressModeClamp;
    force_table.addressMode[1] = cudaAddressModeClamp;
    force_table.filterMode = cudaFilterModeLinear;

    energy_table.normalized = true;
    energy_table.addressMode[0] = cudaAddressModeClamp;
    energy_table.addressMode[1] = cudaAddressModeClamp;
    energy_table.filterMode = cudaFilterModeLinear;

    cudaBindTextureToArray(force_table, ct);
    cuda_errcheck("binding force table to texture");

    cudaBindTextureToArray(energy_table, ect);
    cuda_errcheck("binding energy table to texture");
}

static int patch_pairs_size;
static patch_pair *patch_pairs;
static float *virial_buffers;  // one per patch pair
static float *slow_virial_buffers;  // one per patch pair

static int block_flags_size;
static unsigned int *block_flags;

static int force_lists_size;
static force_list *force_lists;
static unsigned int *force_list_counters;

static int force_buffers_size;
static float4 *force_buffers;
static float4 *slow_force_buffers;

static int atoms_size;
static atom *atoms;
static atom_param *atom_params;
static float4 *forces;
static float4 *slow_forces;
static float *virials;  // one per patch
static float *slow_virials;  // one per patch

static int patch_pairs_alloc;
static int block_flags_alloc;
static int force_buffers_alloc;
static int force_lists_alloc;
static int atoms_alloc;

static int max_atoms_per_patch;

// static cudaStream_t stream;
cudaStream_t stream;
 
void cuda_init() {
  forces = 0;
  slow_forces = 0;
  virials = 0;
  slow_virials = 0;
  atom_params = 0;
  atoms = 0;
  force_buffers = 0;
  slow_force_buffers = 0;
  force_lists = 0;
  force_list_counters = 0;
  patch_pairs = 0;
  virial_buffers = 0;
  slow_virial_buffers = 0;
  block_flags = 0;

  patch_pairs_alloc = 0;
  block_flags_alloc = 0;
  force_buffers_alloc = 0;
  force_lists_alloc = 0;
  atoms_alloc = 0;

  cudaStreamCreate(&stream);
  cuda_errcheck("cudaStreamCreate");
}

void cuda_bind_patch_pairs(const patch_pair *pp, int npp,
                        const force_list *fl, int nfl,
                        int atoms_size_p, int force_buffers_size_p,
                        int block_flags_size_p, int max_atoms_per_patch_p) {

  patch_pairs_size = npp;
  force_buffers_size = force_buffers_size_p;
  force_lists_size = nfl;
  atoms_size = atoms_size_p;
  block_flags_size = block_flags_size_p;
  max_atoms_per_patch = max_atoms_per_patch_p;

#if 0
 printf("%d %d %d %d %d %d %d %d\n",
      patch_pairs_size , patch_pairs_alloc ,
      force_buffers_size , force_buffers_alloc ,
      force_lists_size , force_lists_alloc ,
      atoms_size , atoms_alloc );
#endif

 if ( patch_pairs_size > patch_pairs_alloc ||
      block_flags_size > block_flags_alloc ||
      force_buffers_size > force_buffers_alloc ||
      force_lists_size > force_lists_alloc ||
      atoms_size > atoms_alloc ) {

  block_flags_alloc = (int) (1.2 * block_flags_size);
  patch_pairs_alloc = (int) (1.2 * patch_pairs_size);
  force_buffers_alloc = (int) (1.2 * force_buffers_size);
  force_lists_alloc = (int) (1.2 * force_lists_size);
  atoms_alloc = (int) (1.2 * atoms_size);

  // if ( forces ) cudaFree(forces);
  // if ( slow_forces ) cudaFree(slow_forces);
  forces = slow_forces = 0;
  if ( atom_params ) cudaFree(atom_params);
  if ( atoms ) cudaFree(atoms);
  if ( force_buffers ) cudaFree(force_buffers);
  if ( slow_force_buffers ) cudaFree(slow_force_buffers);
  if ( force_lists ) cudaFree(force_lists);
  if ( force_list_counters ) cudaFree(force_list_counters);
  // if ( virials ) cudaFree(virials);
  virials = slow_virials = 0;
  if ( patch_pairs ) cudaFree(patch_pairs);
  if ( virial_buffers ) cudaFree(virial_buffers);
  if ( slow_virial_buffers ) cudaFree(slow_virial_buffers);
  if ( block_flags ) cudaFree(block_flags);
  cuda_errcheck("free everything");

#if 0
  int totalmem = patch_pairs_alloc * sizeof(patch_pair) +
		force_lists_alloc * sizeof(force_list) +
		2 * force_buffers_alloc * sizeof(float4) +
		atoms_alloc * sizeof(atom) +
		atoms_alloc * sizeof(atom_param) +
		2 * atoms_alloc * sizeof(float4);
  // printf("allocating %d MB of memory on GPU\n", totalmem >> 20);
  printf("allocating %d MB of memory for block flags\n",
				(block_flags_alloc * 4) >> 20);
#endif

  cudaMalloc((void**) &block_flags, block_flags_alloc * 4);
  cudaMalloc((void**) &virial_buffers, patch_pairs_alloc * 16*sizeof(float));
  cudaMalloc((void**) &slow_virial_buffers, patch_pairs_alloc * 16*sizeof(float));
  cudaMalloc((void**) &patch_pairs, patch_pairs_alloc * sizeof(patch_pair));
  // cudaMalloc((void**) &virials, 2 * force_lists_alloc * 16*sizeof(float));
  // slow_virials = virials + force_lists_size * 16;
  cudaMalloc((void**) &force_lists, force_lists_alloc * sizeof(force_list));
  cudaMalloc((void**) &force_list_counters, force_lists_alloc * sizeof(unsigned int));
  cudaMalloc((void**) &force_buffers, force_buffers_alloc * sizeof(float4));
  cudaMalloc((void**) &slow_force_buffers, force_buffers_alloc * sizeof(float4));
  cudaMalloc((void**) &atoms, atoms_alloc * sizeof(atom));
  cudaMalloc((void**) &atom_params, atoms_alloc * sizeof(atom_param));
  // cudaMalloc((void**) &forces, atoms_alloc * sizeof(float4));
  // cudaMalloc((void**) &slow_forces, atoms_alloc * sizeof(float4));
  cuda_errcheck("malloc everything");

 }

  cudaMemcpy(patch_pairs, pp, npp * sizeof(patch_pair),
				cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to patch_pairs");

  cudaMemcpy(force_lists, fl, nfl * sizeof(force_list),
				cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to force_lists");

  cudaMemset(force_list_counters, 0, nfl * sizeof(unsigned int));
  cuda_errcheck("memset force_list_counters");
}

void cuda_bind_atom_params(const atom_param *t) {
  cudaMemcpyAsync(atom_params, t, atoms_size * sizeof(atom_param),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to atom_params");
}

void cuda_bind_atoms(const atom *a) {
  cuda_errcheck("before memcpy to atoms");
  cudaMemcpyAsync(atoms, a, atoms_size * sizeof(atom),
				cudaMemcpyHostToDevice, stream);
  cuda_errcheck("memcpy to atoms");
}

void cuda_bind_forces(float4 *f, float4 *f_slow) {
  cudaHostGetDevicePointer(&forces, f, 0);
  cuda_errcheck("cudaHostGetDevicePointer forces");
  cudaHostGetDevicePointer(&slow_forces, f_slow, 0);
  cuda_errcheck("cudaHostGetDevicePointer slow_forces");
}

void cuda_bind_virials(float *v) {
  cudaHostGetDevicePointer(&virials, v, 0);
  cuda_errcheck("cudaHostGetDevicePointer virials");
  slow_virials = virials + force_lists_size*16;
}

#if 0
void cuda_load_forces(float4 *f, float4 *f_slow, int begin, int count) {
  // printf("load forces %d %d %d\n",begin,count,atoms_size);
  cudaMemcpyAsync(f+begin, forces+begin, count * sizeof(float4),
				cudaMemcpyDeviceToHost, stream);
  if ( f_slow ) {
    cudaMemcpyAsync(f_slow+begin, slow_forces+begin, count * sizeof(float4),
				cudaMemcpyDeviceToHost, stream);
  }
  cuda_errcheck("memcpy from forces");
}

void cuda_load_virials(float *v, int doSlow) {
  int count = force_lists_size;
  if ( doSlow ) count *= 2;
  cudaMemcpyAsync(v, virials, count * 16*sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cuda_errcheck("memcpy from virials");
}
#endif

#if 0
__host__ __device__ static int3 patch_coords_from_id(
        dim3 PATCH_GRID, int id) {

  return make_int3( id % PATCH_GRID.x,
                ( id / PATCH_GRID.x ) % PATCH_GRID.y,
                id / ( PATCH_GRID.x * PATCH_GRID.y ) );
}

__host__ __device__ static int patch_id_from_coords(
        dim3 PATCH_GRID, int3 coords) {

  // handles periodic boundaries
  int x = (coords.x + 4 * PATCH_GRID.x) % PATCH_GRID.x;
  int y = (coords.y + 4 * PATCH_GRID.y) % PATCH_GRID.y;
  int z = (coords.z + 4 * PATCH_GRID.z) % PATCH_GRID.z;

  return ( z * PATCH_GRID.y + y ) * PATCH_GRID.x + x;
}

__host__ __device__ static int3 patch_offset_from_neighbor(int neighbor) {

  // int3 coords = patch_coords_from_id(make_uint3(3,3,3), 13 + neighbor);
  int3 coords = patch_coords_from_id(make_uint3(3,3,3), neighbor);
  return make_int3(coords.x - 1, coords.y - 1, coords.z - 1);

}
#endif
 
#define BLOCK_SIZE 128
#define SHARED_SIZE 32


#define MAKE_PAIRLIST
#define DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef MAKE_PAIRLIST
#define DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_SLOW
#define DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"
#undef DO_ENERGY
#include "ComputeNonbondedCUDAKernelBase.h"


void cuda_nonbonded_forces(float3 lata, float3 latb, float3 latc,
		float cutoff2, float plcutoff2,
		int cbegin, int ccount, int pbegin, int pcount,
		int doSlow, int doEnergy, int usePairlists, int savePairlists) {

 if ( ccount ) {
   if ( usePairlists ) {
     if ( ! savePairlists ) plcutoff2 = 0.;
   } else {
     plcutoff2 = cutoff2;
   }
   int grid_dim = 65535;  // maximum allowed
   for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
     if ( grid_dim > ccount - cstart ) grid_dim = ccount - cstart;
     // printf("%d %d %d\n",cbegin+cstart,grid_dim,patch_pairs_size);

#define CALL(X) X<<< grid_dim, BLOCK_SIZE, 0, stream \
	>>>(patch_pairs+cbegin+cstart,atoms,atom_params,force_buffers, \
	     (doSlow?slow_force_buffers:0), block_flags, \
             virial_buffers, (doSlow?slow_virial_buffers:0), \
             overflow_exclusions, force_list_counters, force_lists, \
             forces, virials, \
             (doSlow?slow_forces:0), (doSlow?slow_virials:0), \
             lj_table_size, \
	     lata, latb, latc, cutoff2, plcutoff2, doSlow)

     if ( doEnergy ) {
       if ( doSlow ) {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_slow_energy_pairlist);
         else CALL(dev_nonbonded_slow_energy);
       } else {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_energy_pairlist);
         else CALL(dev_nonbonded_energy);
       }
     } else {
       if ( doSlow ) {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_slow_pairlist);
         else CALL(dev_nonbonded_slow);
       } else {
         if ( plcutoff2 != 0. ) CALL(dev_nonbonded_pairlist);
         else CALL(dev_nonbonded);
       }
     }

     cuda_errcheck("dev_nonbonded");
   }
 }

#if 0
 if ( pcount ) {
  // printf("%d %d %d\n",pbegin,pcount,force_lists_size);
  dev_sum_forces<<< pcount, BLOCK_SIZE, 0, stream
	>>>(atoms,force_lists+pbegin,force_buffers,
                virial_buffers,forces,virials);
  if ( doSlow ) {
    dev_sum_forces<<< pcount, BLOCK_SIZE, 0, stream
	>>>(atoms,force_lists+pbegin,slow_force_buffers,
                slow_virial_buffers,slow_forces,slow_virials);
  }
  cuda_errcheck("dev_sum_forces");
 }
#endif

}


int cuda_stream_finished() {
  return ( cudaStreamQuery(stream) == cudaSuccess );
}


#else  // NAMD_CUDA

// for make depends
#include "ComputeNonbondedCUDAKernelBase.h"

#endif  // NAMD_CUDA

