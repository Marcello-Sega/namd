
#include "ComputeNonbondedCUDAKernel.h"
#include <stdio.h>

#ifdef NAMD_CUDA


__constant__ unsigned int exclusions[MAX_EXCLUSIONS];

#define SET_EXCL(EXCL,BASE,DIFF) \
         (EXCL)[((BASE)+(DIFF))>>5] |= (1<<(((BASE)+(DIFF))&31))

void cuda_bind_exclusions(const unsigned int *t, int n) {

  cudaMemcpyToSymbol(exclusions, t, n*sizeof(unsigned int), 0);
  cuda_errcheck("memcpy to exclusions");
}


texture<float4, 1, cudaReadModeElementType> force_table;

void cuda_bind_force_table(const float4 *t) {
    static cudaArray *ct;
    if ( ! ct ) {
      cudaMallocArray(&ct, &force_table.channelDesc, FORCE_TABLE_SIZE, 1);
      cuda_errcheck("allocating force table");
    }     cudaMemcpyToArray(ct, 0, 0, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);     // cudaMemcpy(ct, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    cuda_errcheck("memcpy to force table");

    force_table.normalized = true;
    force_table.addressMode[0] = cudaAddressModeClamp;
    force_table.addressMode[1] = cudaAddressModeClamp;
    force_table.filterMode = cudaFilterModeLinear;

    cudaBindTextureToArray(force_table, ct);
    cuda_errcheck("binding force table to texture");
}

static int patch_pairs_size;
static patch_pair *patch_pairs;

static int force_lists_size;
static force_list *force_lists;

static int force_buffers_size;
static float4 *force_buffers;
static float4 *slow_force_buffers;

static int atoms_size;
static atom *atoms;
static atom_param *atom_params;
static float4 *forces;
static float4 *slow_forces;

static int patch_pairs_alloc;
static int force_buffers_alloc;
static int force_lists_alloc;
static int atoms_alloc;

static int max_atoms_per_patch;

// static cudaStream_t stream;
cudaStream_t stream;
 
void cuda_init() {
  forces = 0;
  slow_forces = 0;
  atom_params = 0;
  atoms = 0;
  force_buffers = 0;
  slow_force_buffers = 0;
  force_lists = 0;
  patch_pairs = 0;

  patch_pairs_alloc = 0;
  force_buffers_alloc = 0;
  force_lists_alloc = 0;
  atoms_alloc = 0;

  cudaStreamCreate(&stream);
  cuda_errcheck("cudaStreamCreate");
}

void cuda_bind_patch_pairs(const patch_pair *pp, int npp,
                        const force_list *fl, int nfl,
                        int atoms_size_p, int force_buffers_size_p,
			int max_atoms_per_patch_p) {

  patch_pairs_size = npp;
  force_buffers_size = force_buffers_size_p;
  force_lists_size = nfl;
  atoms_size = atoms_size_p;
  max_atoms_per_patch = max_atoms_per_patch_p;

#if 0
 printf("%d %d %d %d %d %d %d %d\n",
      patch_pairs_size , patch_pairs_alloc ,
      force_buffers_size , force_buffers_alloc ,
      force_lists_size , force_lists_alloc ,
      atoms_size , atoms_alloc );
#endif

 if ( patch_pairs_size > patch_pairs_alloc ||
      force_buffers_size > force_buffers_alloc ||
      force_lists_size > force_lists_alloc ||
      atoms_size > atoms_alloc ) {

  patch_pairs_alloc = (int) (1.2 * patch_pairs_size);
  force_buffers_alloc = (int) (1.2 * force_buffers_size);
  force_lists_alloc = (int) (1.2 * force_lists_size);
  atoms_alloc = (int) (1.2 * atoms_size);

  if ( forces ) cudaFree(forces);
  if ( slow_forces ) cudaFree(slow_forces);
  if ( atom_params ) cudaFree(atom_params);
  if ( atoms ) cudaFree(atoms);
  if ( force_buffers ) cudaFree(force_buffers);
  if ( slow_force_buffers ) cudaFree(slow_force_buffers);
  if ( force_lists ) cudaFree(force_lists);
  if ( patch_pairs ) cudaFree(patch_pairs);
  cuda_errcheck("free everything");

#if 1
  int totalmem = patch_pairs_alloc * sizeof(patch_pair) +
		force_lists_alloc * sizeof(force_list) +
		2 * force_buffers_alloc * sizeof(float4) +
		atoms_alloc * sizeof(atom) +
		atoms_alloc * sizeof(atom_param) +
		2 * atoms_alloc * sizeof(float4);
  // printf("allocating %d MB of memory on GPU\n", totalmem >> 20);
#endif

  cudaMalloc((void**) &patch_pairs, patch_pairs_alloc * sizeof(patch_pair));
  cudaMalloc((void**) &force_lists, force_lists_alloc * sizeof(force_list));
  cudaMalloc((void**) &force_buffers, force_buffers_alloc * sizeof(float4));
  cudaMalloc((void**) &slow_force_buffers, force_buffers_alloc * sizeof(float4));
  cudaMalloc((void**) &atoms, atoms_alloc * sizeof(atom));
  cudaMalloc((void**) &atom_params, atoms_alloc * sizeof(atom_param));
  cudaMalloc((void**) &forces, atoms_alloc * sizeof(float4));
  cudaMalloc((void**) &slow_forces, atoms_alloc * sizeof(float4));
  cuda_errcheck("malloc everything");

 }

  cudaMemcpy(patch_pairs, pp, npp * sizeof(patch_pair),
				cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to patch_pairs");

  cudaMemcpy(force_lists, fl, nfl * sizeof(force_list),
				cudaMemcpyHostToDevice);
  cuda_errcheck("memcpy to force_lists");
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

__global__ static void dev_nonbonded(
	const patch_pair *patch_pairs,
	const atom *atoms,
	const atom_param *atom_params,
	float4 *force_buffers,
	float4 *slow_force_buffers,
        float3 lata, float3 latb, float3 latc,
	float cutoff2) {
// call with two blocks per patch_pair
// call with BLOCK_SIZE threads per block
// call with no shared memory

  #define jpqs jpqu.d
  __shared__ union {
    atom d[SHARED_SIZE];
    unsigned int i[BLOCK_SIZE];
    float f[BLOCK_SIZE];
  } jpqu;

  #define japs japu.d
  __shared__ union {
    atom_param d[SHARED_SIZE];
    unsigned int i[BLOCK_SIZE];
  } japu;

  #define myPatchPair pp.pp
  __shared__ union { patch_pair pp; unsigned int i[PATCH_PAIR_SIZE]; } pp;

  if ( threadIdx.x < PATCH_PAIR_USED ) {
    unsigned int tmp = ((unsigned int*)patch_pairs)[
			(sizeof(patch_pair)>>2)*blockIdx.x+threadIdx.x];
    pp.i[threadIdx.x] = tmp;
  }
  __syncthreads();

  // convert scaled offset with current lattice
  if ( threadIdx.x == 0 ) {
    float offx = myPatchPair.offset.x * lata.x
               + myPatchPair.offset.y * latb.x
               + myPatchPair.offset.z * latc.x;
    float offy = myPatchPair.offset.x * lata.y
               + myPatchPair.offset.y * latb.y
               + myPatchPair.offset.z * latc.y;
    float offz = myPatchPair.offset.x * lata.z
               + myPatchPair.offset.y * latb.z
               + myPatchPair.offset.z * latc.z;
    myPatchPair.offset.x = offx;
    myPatchPair.offset.y = offy;
    myPatchPair.offset.z = offz;
  }

  __syncthreads();

  for ( int blocki = 0;
        blocki < myPatchPair.patch1_force_size;
        blocki += BLOCK_SIZE ) {

  atom ipq;
  struct {
    float sqrt_epsilon;
    float half_sigma;
    int index; } iap;

  // load patch 1
  if ( blocki + threadIdx.x < myPatchPair.patch1_force_size ) {
    int i = myPatchPair.patch1_atom_start + blocki + threadIdx.x;
    float4 tmpa = ((float4*)atoms)[i];

    ipq.position.x = tmpa.x + myPatchPair.offset.x;
    ipq.position.y = tmpa.y + myPatchPair.offset.y;
    ipq.position.z = tmpa.z + myPatchPair.offset.z;
    ipq.charge = tmpa.w;

    uint4 tmpap = ((uint4*)atom_params)[i];

    iap.sqrt_epsilon = __int_as_float(tmpap.x);
    iap.half_sigma = __int_as_float(tmpap.y);
    iap.index = tmpap.z;
  }

  float4 ife, ife_slow;
  ife.x = 0.f;
  ife.y = 0.f;
  ife.z = 0.f;
  ife.w = 0.f;
  ife_slow.x = 0.f;
  ife_slow.y = 0.f;
  ife_slow.z = 0.f;
  ife_slow.w = 0.f;

  for ( int blockj = 0;
        blockj < myPatchPair.patch2_size;
        blockj += SHARED_SIZE ) {

  int shared_size = myPatchPair.patch2_size - blockj;
  if ( shared_size > SHARED_SIZE ) shared_size = SHARED_SIZE;

  // load patch 2
  // sync needed because of loop, could avoid with double-buffering
  __syncthreads();

  if ( threadIdx.x < 4 * shared_size ) {
    int j = myPatchPair.patch2_atom_start + blockj;
    jpqu.i[threadIdx.x] = ((unsigned int *)(atoms + j))[threadIdx.x];
    japu.i[threadIdx.x] = ((unsigned int *)(atom_params + j))[threadIdx.x];
  }
  __syncthreads();

  // calc forces on patch 1
  if ( blocki + threadIdx.x < myPatchPair.patch1_force_size ) {

// be careful not to use // comments inside macros!
#define FORCE_INNER_LOOP(IPQ,IAP,DO_SLOW) \
    for ( int j = 0; j < shared_size; ++j ) { \
      /* actually calculate force */ \
      float tmpx = jpqs[j].position.x - IPQ.position.x; \
      float tmpy = jpqs[j].position.y - IPQ.position.y; \
      float tmpz = jpqs[j].position.z - IPQ.position.z; \
      float r2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz; \
      if ( r2 < cutoff2 ) { \
        float4 fi = tex1D(force_table, rsqrtf(r2)); \
        bool excluded = false; \
        int indexdiff = (int)(IAP.index) - (int)(japs[j].index); \
        if ( abs(indexdiff) <= (int) japs[j].excl_maxdiff ) { \
          indexdiff += japs[j].excl_index; \
          excluded = ((exclusions[indexdiff>>5] & (1<<(indexdiff&31))) != 0); \
        } \
        float e = IAP.half_sigma + japs[j].half_sigma;  /* sigma */ \
        e *= e*e;  /* sigma^3 */ \
        e *= e;  /* sigma^6 */ \
        e *= ( e * fi.z + fi.y );  /* s^12 * fi.z - s^6 * fi.y */ \
        e *= IAP.sqrt_epsilon * japs[j].sqrt_epsilon;  /* full L-J */ \
        float e_slow = IPQ.charge * jpqs[j].charge; \
        e += e_slow * fi.x; \
        if ( DO_SLOW ) e_slow *= fi.w; \
        if ( ! excluded ) { \
          ife.w += r2 * e; \
          ife.x += tmpx * e; \
          ife.y += tmpy * e; \
          ife.z += tmpz * e; \
          if ( DO_SLOW ) { \
          ife_slow.w += r2 * e_slow; \
          ife_slow.x += tmpx * e_slow; \
          ife_slow.y += tmpy * e_slow; \
          ife_slow.z += tmpz * e_slow; \
          } \
        } \
      }  /* cutoff */ \
    } /* end of FORCE_INNER_LOOP macro */

    if ( slow_force_buffers ) {
      FORCE_INNER_LOOP(ipq,iap,1)
    } else {
      FORCE_INNER_LOOP(ipq,iap,0)
    }

  } // if
  } // blockj loop

  if ( blocki + threadIdx.x < myPatchPair.patch1_force_size ) {
    int i_out = myPatchPair.patch1_force_start + blocki + threadIdx.x;
    force_buffers[i_out] = ife;
    if ( slow_force_buffers ) {
      slow_force_buffers[i_out] = ife_slow;
    }
  }

  } // blocki loop

}


__global__ static void dev_sum_forces(
	const force_list *force_lists,
	const float4 *force_buffers,
	float4 *forces) {
// call with one block per patch
// call multiple of 64 threads per block
// call with no shared memory

  __shared__ force_list myForceList;

  if ( threadIdx.x == 0 ) {
    myForceList = force_lists[blockIdx.x];
  }
  __syncthreads();

  for ( int j = threadIdx.x; j < myForceList.patch_size; j += blockDim.x ) {

    const float4 *fbuf = force_buffers + myForceList.force_list_start + j;
    float4 fout;
    fout.x = 0.f;
    fout.y = 0.f;
    fout.z = 0.f;
    fout.w = 0.f;
    for ( int i=0; i < myForceList.force_list_size; ++i ) {
      float4 f = *fbuf;
      fout.x += f.x;
      fout.y += f.y;
      fout.z += f.z;
      fout.w += f.w;
      fbuf += myForceList.patch_size;
    }

    forces[myForceList.force_output_start + j] = fout;

  }
}


void cuda_nonbonded_forces(float3 lata, float3 latb, float3 latc, float cutoff2,
		int cbegin, int ccount, int pbegin, int pcount, int doSlow) {

 if ( ccount ) {
   int grid_dim = 65535;  // maximum allowed
   for ( int cstart = 0; cstart < ccount; cstart += grid_dim ) {
     if ( grid_dim > ccount - cstart ) grid_dim = ccount - cstart;
     // printf("%d %d %d\n",cbegin+cstart,grid_dim,patch_pairs_size);
     dev_nonbonded<<< grid_dim, BLOCK_SIZE, 0, stream
	>>>(patch_pairs+cbegin+cstart,atoms,atom_params,force_buffers,
	     (doSlow?slow_force_buffers:0), lata, latb, latc, cutoff2);
     cuda_errcheck("dev_nonbonded");
   }
 }

 if ( pcount ) {
  // printf("%d %d %d\n",pbegin,pcount,force_lists_size);
  dev_sum_forces<<< pcount, 128, 0, stream
	>>>(force_lists+pbegin,force_buffers,forces);
  if ( doSlow ) {
    dev_sum_forces<<< pcount, 128, 0, stream
	>>>(force_lists+pbegin,slow_force_buffers,slow_forces);
  }
  cuda_errcheck("dev_sum_forces");
 }

}


int cuda_stream_finished() {
  return ( cudaStreamQuery(stream) == cudaSuccess );
}


#endif  // NAMD_CUDA

