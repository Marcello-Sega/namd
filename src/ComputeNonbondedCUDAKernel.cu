
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
    }
    cudaMemcpyToArray(ct, 0, 0, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
    // cudaMemcpy(ct, t, FORCE_TABLE_SIZE*sizeof(float4), cudaMemcpyHostToDevice);
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

  if ( forces ) cudaFree(forces);
  if ( slow_forces ) cudaFree(slow_forces);
  if ( atom_params ) cudaFree(atom_params);
  if ( atoms ) cudaFree(atoms);
  if ( force_buffers ) cudaFree(force_buffers);
  if ( slow_force_buffers ) cudaFree(slow_force_buffers);
  if ( force_lists ) cudaFree(force_lists);
  if ( force_list_counters ) cudaFree(force_list_counters);
  if ( virials ) cudaFree(virials);
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
  cudaMalloc((void**) &virials, 2 * force_lists_alloc * 16*sizeof(float));
  slow_virials = virials + force_lists_size * 16;
  cudaMalloc((void**) &force_lists, force_lists_alloc * sizeof(force_list));
  cudaMalloc((void**) &force_list_counters, force_lists_alloc * sizeof(unsigned int));
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

__device__ __forceinline__ static void dev_sum_forces(
        const int force_list_index,
	const atom *atoms,
	const force_list *force_lists,
	const float4 *force_buffers,
	const float *virial_buffers,
	float4 *forces, float *virials);

__global__ static void dev_nonbonded(
	const patch_pair *patch_pairs,
	const atom *atoms,
	const atom_param *atom_params,
	float4 *force_buffers,
	float4 *slow_force_buffers,
	unsigned int *block_flags,
	float *virial_buffers,
	float *slow_virial_buffers,
        unsigned int *force_list_counters,
        const force_list *force_lists,
        float4 *forces, float *virials,
        float4 *slow_forces, float *slow_virials,
        float3 lata, float3 latb, float3 latc,
	float cutoff2, float plcutoff2, int doSlow) {
// call with one block per patch_pair
// call with BLOCK_SIZE threads per block
// call with no shared memory

#ifdef __DEVICE_EMULATION__
  #define myPatchPair (*(patch_pair*)(&pp.i))
#else
  #define myPatchPair pp.pp
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    patch_pair pp;
#endif
    unsigned int i[PATCH_PAIR_SIZE];
  } pp;

 { // start of nonbonded calc

  #define pl plu.c
  __shared__ union {
    unsigned int i[BLOCK_SIZE];
    char c[4*BLOCK_SIZE];
  } plu;

  volatile __shared__ union {
    float a2d[32][3];
    float a1d[32*3];
  } sumf;

  volatile __shared__ union {
    float a2d[32][3];
    float a1d[32*3];
  } sumf_slow;

#ifdef __DEVICE_EMULATION__
  #define jpqs ((atom*)(jpqu.i))
#else
  #define jpqs jpqu.d
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    atom d[SHARED_SIZE];
#endif
    unsigned int i[4*SHARED_SIZE];
    float f[4*SHARED_SIZE];
  } jpqu;

#ifdef __DEVICE_EMULATION__
  #define japs ((atom_param*)(japu.i))
#else
  #define japs japu.d
#endif
  __shared__ union {
#ifndef __DEVICE_EMULATION__
    atom_param d[SHARED_SIZE];
#endif
    unsigned int i[4*SHARED_SIZE];
  } japu;

  if ( threadIdx.x < PATCH_PAIR_USED ) {
    unsigned int tmp = ((unsigned int*)patch_pairs)[
			PATCH_PAIR_SIZE*blockIdx.x+threadIdx.x];
    pp.i[threadIdx.x] = tmp;
  }

  if ( threadIdx.x < 96 ) { // initialize net force in shared memory
    sumf.a1d[threadIdx.x] = 0.f;
    sumf_slow.a1d[threadIdx.x] = 0.f;
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

  // avoid syncs by having all warps load pairlist
  if ( plcutoff2 == 0 ) {
    int i_pl = (blocki >> 2) + myPatchPair.block_flags_start;
    plu.i[threadIdx.x] = block_flags[i_pl + (threadIdx.x & 31)];
  } else {
    plu.i[threadIdx.x] = 0;
  }
  int pli = 4 * ( threadIdx.x & 96 );

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
        blockj += SHARED_SIZE, ++pli ) {

#ifdef __DEVICE_EMULATION__
  if ( plcutoff2 == 0 && threadIdx.x == 0 ) printf("%d %d %d %d %d %d %d\n", blockIdx.x, blocki, blockj, pli, pl[pli], (pli+128)&255, pl[(pli+128)&255]);
#endif
  if ( plcutoff2 == 0 && pl[pli] == 0 ) continue;

  int shared_size = myPatchPair.patch2_size - blockj;
  if ( shared_size > SHARED_SIZE ) shared_size = SHARED_SIZE;

  // load patch 2
  __syncthreads();

  if ( threadIdx.x < 4 * shared_size ) {
    int j = myPatchPair.patch2_atom_start + blockj;
    jpqu.i[threadIdx.x] = ((unsigned int *)(atoms + j))[threadIdx.x];
    japu.i[threadIdx.x] = ((unsigned int *)(atom_params + j))[threadIdx.x];
  }
  __syncthreads();

  if ( plcutoff2 == 0 && (pl[pli] & (1 << (threadIdx.x >> 5))) == 0 ) continue;

  // calc forces on patch 1
  if ( blocki + threadIdx.x < myPatchPair.patch1_force_size ) {

// be careful not to use // comments inside macros!
#define FORCE_INNER_LOOP(IPQ,IAP,DO_SLOW,DO_PAIRLIST) \
    for ( int j = 0; j < shared_size; ++j ) { \
      /* actually calculate force */ \
      float tmpx = jpqs[j].position.x - IPQ.position.x; \
      float tmpy = jpqs[j].position.y - IPQ.position.y; \
      float tmpz = jpqs[j].position.z - IPQ.position.z; \
      float r2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz; \
      DO_PAIRLIST \
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
      } }  /* cutoff */ \
    } /* end of FORCE_INNER_LOOP macro */

    if ( plcutoff2 == 0 ) {  // use pairlist
      if ( doSlow ) {
        FORCE_INNER_LOOP(ipq,iap,1,{)
      } else {
        FORCE_INNER_LOOP(ipq,iap,0,{)
      }
    } else {  // create pairlist
      bool plpli = 0;
      if ( doSlow ) {
        FORCE_INNER_LOOP(ipq,iap,1,if(r2<plcutoff2){plpli=1;)
      } else {
        FORCE_INNER_LOOP(ipq,iap,0,if(r2<plcutoff2){plpli=1;)
      }
      if ( plpli ) pl[pli] = 1;
    }

  } // if
  } // blockj loop

  if ( blocki + threadIdx.x < myPatchPair.patch1_force_size ) {
    int i_out = myPatchPair.patch1_force_start + blocki + threadIdx.x;
    force_buffers[i_out] = ife;
    if ( doSlow ) {
      slow_force_buffers[i_out] = ife_slow;
    }
    // accumulate net force to shared memory, warp-synchronous
    const int subwarp = threadIdx.x >> 2;  // 32 entries in table
    const int thread = threadIdx.x & 3;  // 4 threads share each entry
    for ( int g = 0; g < 4; ++g ) {
      if ( thread == g ) {
        sumf.a2d[subwarp][0] += ife.x;
        sumf.a2d[subwarp][1] += ife.y;
        sumf.a2d[subwarp][2] += ife.z;
        if ( doSlow ) {
          sumf_slow.a2d[subwarp][0] += ife_slow.x;
          sumf_slow.a2d[subwarp][1] += ife_slow.y;
          sumf_slow.a2d[subwarp][2] += ife_slow.z;
        }
      }
    }
  }
  if ( plcutoff2 != 0 ) {
    __syncthreads();  // all shared pairlist writes complete
    unsigned int pltmp;
    if ( threadIdx.x < 32 ) {
      pltmp = plu.i[threadIdx.x];
      pltmp |= plu.i[threadIdx.x+32] << 1;
      pltmp |= plu.i[threadIdx.x+64] << 2;
      pltmp |= plu.i[threadIdx.x+96] << 3;
    }
    __syncthreads();  // all shared pairlist reads complete
    if ( threadIdx.x < 32 ) {
      int i_pl = (blocki >> 2) + myPatchPair.block_flags_start;
      block_flags[i_pl + threadIdx.x] = pltmp;
    }
  }

  } // blocki loop

  __syncthreads();
  if ( threadIdx.x < 24 ) { // reduce forces, warp-synchronous
                            // 3 components, 8 threads per component
    const int i_out = myPatchPair.virial_start + threadIdx.x;
    {
      float f;
      f = sumf.a1d[threadIdx.x] + sumf.a1d[threadIdx.x + 24] + 
          sumf.a1d[threadIdx.x + 48] + sumf.a1d[threadIdx.x + 72];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 12];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 6];
      sumf.a1d[threadIdx.x] = f;
      f += sumf.a1d[threadIdx.x + 3];
      f *= 0.5f;  // compensate for double-counting
      // calculate virial contribution on first 3 threads
      sumf.a2d[threadIdx.x][0] = f * myPatchPair.offset.x;
      sumf.a2d[threadIdx.x][1] = f * myPatchPair.offset.y;
      sumf.a2d[threadIdx.x][2] = f * myPatchPair.offset.z;
      if ( threadIdx.x < 9 ) {  // write out output buffer
        virial_buffers[i_out] = sumf.a1d[threadIdx.x];
      }
    }
    if ( doSlow ) { // repeat above for slow forces
      float fs;
      fs = sumf_slow.a1d[threadIdx.x] + sumf_slow.a1d[threadIdx.x + 24] + 
           sumf_slow.a1d[threadIdx.x + 48] + sumf_slow.a1d[threadIdx.x + 72];
      sumf_slow.a1d[threadIdx.x] = fs;
      fs += sumf_slow.a1d[threadIdx.x + 12];
      sumf_slow.a1d[threadIdx.x] = fs;
      fs += sumf_slow.a1d[threadIdx.x + 6];
      sumf_slow.a1d[threadIdx.x] = fs;
      fs += sumf_slow.a1d[threadIdx.x + 3];
      fs *= 0.5f;
      sumf_slow.a2d[threadIdx.x][0] = fs * myPatchPair.offset.x;
      sumf_slow.a2d[threadIdx.x][1] = fs * myPatchPair.offset.y;
      sumf_slow.a2d[threadIdx.x][2] = fs * myPatchPair.offset.z;
      if ( threadIdx.x < 9 ) {
        slow_virial_buffers[i_out] = sumf_slow.a1d[threadIdx.x];
      }
    }
  }

 } // end of nonbonded calc

 { // start of force sum

  // make sure forces are visible in global memory
  __threadfence();

  __shared__ bool sumForces;

  if (threadIdx.x == 0) {
    int fli = myPatchPair.patch1_force_list_index;
    int fls = myPatchPair.patch1_force_list_size;
    int old = atomicInc(force_list_counters+fli,fls-1);
    sumForces = ( old == fls - 1 );
  }

  __syncthreads();

  if ( sumForces ) {
    dev_sum_forces(myPatchPair.patch1_force_list_index,
       atoms,force_lists,force_buffers,
       virial_buffers,forces,virials);

    if ( doSlow ) {
      dev_sum_forces(myPatchPair.patch1_force_list_index,
         atoms,force_lists,slow_force_buffers,
         slow_virial_buffers,slow_forces,slow_virials);
    }
  }

 } // end of force sum
}


__device__ __forceinline__ static void dev_sum_forces(
        const int force_list_index,
	const atom *atoms,
	const force_list *force_lists,
	const float4 *force_buffers,
	const float *virial_buffers,
	float4 *forces, float *virials) {
// call with one block per patch
// call BLOCK_SIZE threads per block
// call with no shared memory

  #define myForceList fl.fl
  __shared__ union {
    force_list fl;
    unsigned int i[FORCE_LIST_SIZE];
  } fl;

  if ( threadIdx.x < FORCE_LIST_USED ) {
    unsigned int tmp = ((unsigned int*)force_lists)[
                        FORCE_LIST_SIZE*force_list_index+threadIdx.x];
    fl.i[threadIdx.x] = tmp;
  }

  volatile __shared__ union {
    float a3d[32][3][3];
    float a2d[32][9];
    float a1d[32*9];
  } virial;

  for ( int i = threadIdx.x; i < 32*9; i += BLOCK_SIZE ) {
    virial.a1d[i] = 0.f;
  }

  __syncthreads();

  float vxx = 0.f;
  float vxy = 0.f;
  float vxz = 0.f;
  float vyx = 0.f;
  float vyy = 0.f;
  float vyz = 0.f;
  float vzx = 0.f;
  float vzy = 0.f;
  float vzz = 0.f;

  for ( int j = threadIdx.x; j < myForceList.patch_size; j += BLOCK_SIZE ) {

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
      fbuf += myForceList.patch_stride;
    }

    // compiler will use st.global.f32 instead of st.global.v4.f32
    // if forcedest is directly substituted in the assignment
    const int forcedest = myForceList.force_output_start + j;
    forces[forcedest] = fout;

    float4 pos = ((float4*)atoms)[myForceList.atom_start + j];

    // accumulate per-atom virials to registers
    vxx += fout.x * pos.x;
    vxy += fout.x * pos.y;
    vxz += fout.x * pos.z;
    vyx += fout.y * pos.x;
    vyy += fout.y * pos.y;
    vyz += fout.y * pos.z;
    vzx += fout.z * pos.x;
    vzy += fout.z * pos.y;
    vzz += fout.z * pos.z;

  }

  { // accumulate per-atom virials to shared memory, warp-synchronous
    const int subwarp = threadIdx.x >> 2;  // 32 entries in table
    const int thread = threadIdx.x & 3;  // 4 threads share each entry
    for ( int g = 0; g < 4; ++g ) {
      if ( thread == g ) {
        virial.a3d[subwarp][0][0] += vxx;
        virial.a3d[subwarp][0][1] += vxy;
        virial.a3d[subwarp][0][2] += vxz;
        virial.a3d[subwarp][1][0] += vyx;
        virial.a3d[subwarp][1][1] += vyy;
        virial.a3d[subwarp][1][2] += vyz;
        virial.a3d[subwarp][2][0] += vzx;
        virial.a3d[subwarp][2][1] += vzy;
        virial.a3d[subwarp][2][2] += vzz;
      }
    }
  }
  __syncthreads();
  { // accumulate per-compute virials to shared memory, data-parallel
    const int halfwarp = threadIdx.x >> 4;  // 8 half-warps
    const int thread = threadIdx.x & 15;
    if ( thread < 9 ) {
      for ( int i = halfwarp; i < myForceList.force_list_size; i += 8 ) {
        virial.a2d[halfwarp][thread] +=
          virial_buffers[myForceList.virial_list_start + 16*i + thread];
      }
    }
  }
  __syncthreads();
  { // reduce virials in shared memory, warp-synchronous
    const int subwarp = threadIdx.x >> 3;  // 16 quarter-warps
    const int thread = threadIdx.x & 7;  // 8 threads per component
    if ( subwarp < 9 ) {  // 9 components
      float v;
      v = virial.a2d[thread][subwarp] + virial.a2d[thread+8][subwarp] +
          virial.a2d[thread+16][subwarp] + virial.a2d[thread+24][subwarp];
      virial.a2d[thread][subwarp] = v;
      v += virial.a2d[thread+4][subwarp];
      virial.a2d[thread][subwarp] = v;
      v += virial.a2d[thread+2][subwarp];
      virial.a2d[thread][subwarp] = v;
      v += virial.a2d[thread+1][subwarp];
      virial.a2d[thread][subwarp] = v;
    }
  }
  __syncthreads();
  if ( threadIdx.x < 9 ) {  // 9 components
    virials[myForceList.virial_output_start + threadIdx.x] =
                                              virial.a2d[0][threadIdx.x];
  }

}


void cuda_nonbonded_forces(float3 lata, float3 latb, float3 latc,
		float cutoff2, float plcutoff2,
		int cbegin, int ccount, int pbegin, int pcount,
		int doSlow, int usePairlists, int savePairlists) {

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
     dev_nonbonded<<< grid_dim, BLOCK_SIZE, 0, stream
	>>>(patch_pairs+cbegin+cstart,atoms,atom_params,force_buffers,
	     (doSlow?slow_force_buffers:0), block_flags,
             virial_buffers, (doSlow?slow_virial_buffers:0),
             force_list_counters, force_lists,
             forces, virials,
             (doSlow?slow_forces:0), (doSlow?slow_virials:0),
	     lata, latb, latc, cutoff2, plcutoff2, doSlow);
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


#endif  // NAMD_CUDA

