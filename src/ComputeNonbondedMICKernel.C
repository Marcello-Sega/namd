// NOTE: See ComputeNonbondedMICKernel.h for a general description of NAMD on MIC.

// Only compile the contents of this file if this build supports MIC.
#ifdef NAMD_MIC


#include <stdlib.h>
#include <stdio.h>
#include <offload.h>
#include <string.h>
#include <stdint.h>

#ifndef __MIC__
  #include "charm++.h"
#endif

#include "ComputeNonbondedMICKernel.h"

#if defined(__MIC__) || defined(__MIC2__)
  #include <immintrin.h>
#endif

#if 0 && MIC_HANDCODE_FORCE != 0
  #include "vec_debug.h"
#endif

#include <assert.h>
#include <math.h>
#if MIC_LOCK_MEMORY
  #include <sys/mman.h>
  #include <errno.h>
#endif

#if USE_DMK_COMMON != 0
  #pragma offload_attribute (push, target(mic))
  #include "dmk_common.C"
  #pragma offload_attribute (pop)
#endif

#if (VTUNE_INFO != 0) && (VTUNE_TASKS != 0)
  #include "ittnotify.h"
  __itt_domain *vtune_domain = __itt_domain_create("NAMD");
  __itt_string_handle *vtune_frame = __itt_string_handle_create("frame.timestep");
  __itt_string_handle *vtune_offload_task = __itt_string_handle_create("non-bonded.offload.task");
  __thread int vtune_frameStartedFlag = 0;
  //void endOffloadTask() { __itt_task_end(vtune_domain); }
  void endOffloadTask() { __itt_frame_end_v3(vtune_domain, NULL); }
#endif


// Setup __ASSUME_ALIGNED macro to take the appropriate action based on macro flags
#if USE_ASSUME_ALIGNED && CHECK_ASSUME_ALIGNED
  #define __ASSUME_ALIGNED(v) assert(((unsigned long long int)(v)) % 64 == 0); __assume_aligned((v), MIC_ALIGN);
#elif USE_ASSUME_ALIGNED
  #define __ASSUME_ALIGNED(v) __assume_aligned((v), MIC_ALIGN);
#else
  #define __ASSUME_ALIGNED(v)
#endif

// Setup __ASSUME macro to take the appropriate action based on macro flags
#if USE_ASSUME && CHECK_ASSUME
  #define __ASSUME(s) assert(s); __assume(s);
#elif USE_ASSUME
  #define __ASSUME(s) __assume(s);
#else
  #define __ASSUME(s)
#endif

// Setup RESTRICT macro to take the appropriate action based on macro flags
#if USE_RESTRICT
  #define RESTRICT __restrict
#else
  #define RESTRICT
#endif

#ifdef WIN32
  #define __thread __declspec(thread)
#endif


////////////////////////////////////////////////////////////////////////////////
// Global data used during simulation (tables, constants, etc. that are setup
//   during startup but essentially read-only and constant throughout the
//   steady-state simulation.

__thread __attribute__((target(mic))) double *device__table_four = NULL;
__thread __attribute__((target(mic))) float *device__table_four_float = NULL;
__thread __attribute__((target(mic))) int device__table_four_n_16 = 0;

__thread __attribute__((target(mic))) char *device__lj_table = NULL;
__thread __attribute__((target(mic))) float *device__lj_table_float = NULL;
__thread __attribute__((target(mic))) int device__lj_table_dim = 0;
__thread __attribute__((target(mic))) int device__lj_table_size = 0;

__thread __attribute__((target(mic))) unsigned int *device__exclusion_bits = NULL;
__thread __attribute__((target(mic))) long int device__exclusion_bits_size = 0;

__thread __attribute__((target(mic))) mic_constants *device__constants = NULL;


////////////////////////////////////////////////////////////////////////////////
// Device variables which exist both on the host and the MIC device and/or are
//   used to manage moving data betwen the host and the MIC device.

__thread __attribute__((target(mic))) int device__myPE = -1;
__thread __attribute__((target(mic))) int device__myNode = -1;

__thread __attribute__((target(mic))) const patch_pair *device__patch_pairs;
__thread __attribute__((target(mic))) int device__patch_pairs_size;
__thread __attribute__((target(mic))) int device__patch_pairs_alloc_size;

__thread __attribute__((target(mic))) unsigned long long int device__pairlists;  // NOTE: Don't use a global in case the device is shared between multiple threads, so each thread's patch pair lists has its own pairlist pointer
__thread __attribute__((target(mic))) int device__pairlists_alloc_size;

__thread __attribute__((target(mic))) int** device__pl_array;
__thread __attribute__((target(mic))) int* device__pl_size;
__thread __attribute__((target(mic))) double** device__r2_array;

__thread __attribute__((target(mic))) force_list *device__force_lists;
__thread __attribute__((target(mic))) int device__force_lists_size;
__thread __attribute__((target(mic))) int device__force_lists_alloc_size;

__thread __attribute__((target(mic))) atom *device__atoms;
__thread __attribute__((target(mic))) atom_param *device__atom_params;
__thread __attribute__((target(mic))) double4 *device__forces;
__thread __attribute__((target(mic))) double4 *device__slow_forces;
__thread __attribute__((target(mic))) int device__atoms_size;
__thread __attribute__((target(mic))) int device__atoms_alloc_size;

__thread __attribute__((target(mic))) double4 *device__force_buffers;
__thread __attribute__((target(mic))) double4 *device__slow_force_buffers;
__thread __attribute__((target(mic))) int device__force_buffers_req_size;
__thread __attribute__((target(mic))) int device__force_buffers_alloc_size;

// DMK - TODO | FIXME - See if these variables can be moved to mic_constants
__thread __attribute__((target(mic))) mic_position3_t device__lata;
__thread __attribute__((target(mic))) mic_position3_t device__latb;
__thread __attribute__((target(mic))) mic_position3_t device__latc;

__thread __attribute__((target(mic))) mic_kernel_data *device__remote_kernel_data;
__thread __attribute__((target(mic))) mic_kernel_data *device__local_kernel_data;

__thread __attribute__((target(mic))) double *device__virial_energy;

__thread __attribute__((target(mic))) int device__numOMPThreads;

__thread int tag_atom_params;
__thread int tag_remote_kernel;
__thread int tag_local_kernel;

__thread int patch_pairs_copySize;
__thread int force_lists_copySize;
__thread int atom_params_copySize;


// DMK - DEBUG - Based on the number of kernels invoked, track the timestep number.
__thread __attribute__((target(mic))) int device__timestep;

// DMK - DEBUG - Counters for some statistics when they are collected.
__thread __attribute__((target(mic))) int* device__loop_counts;
__thread __attribute__((target(mic))) int* device__cutoff_clear_counts;
__thread __attribute__((target(mic))) int* device__iTest_counts;

// DMK - TRACING / TIMING
#include <sys/time.h>
__declspec(target(mic))
double getCurrentTime() {
  timeval now;
  gettimeofday(&now, NULL);
  return (double)(now.tv_sec + (now.tv_usec * 1.0e-6));
}

// DMK - TRACING
#if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
  #if MIC_DEVICE_TRACING_DETAILED != 0
    __thread __attribute__((target(mic))) double* device__device_times_computes; // Two values per "entry" (start, stop)
    __thread __attribute__((target(mic))) double* device__device_times_patches;
  #endif
  __thread __attribute__((target(mic))) double* device__device_times_start;
#endif

// DMK - DEBUG
#if MIC_PRAGMA_TIMING_STATS

  __thread double pragma_timing_start = -1.0;
  __thread double pragma_timing_issued = -1.0;
  __thread double pragma_issue_time_total = 0.0;
  __thread double pragma_wait_time_total = 0.0;
  __thread __attribute__((target(mic))) double pragma_deviceExec_time[1] = { 0.0 };
  __thread double pragma_deviceExec_time_total = 0.0;
  __thread int pragma_num_samples = 0;

  #define MIC_PRAGMA_TIMING_DEVICE_CLAUSE   out(pragma_deviceExec_time[0:1] : alloc_if(0) free_if(0))

  void pragma_timing_submitFinish(double finishTime, int pe, int rank) {

    //// DMK - DEBUG
    //printf("submitFinish :: &pragma_deviceExec_time = %p\n", &pragma_deviceExec_time);
    printf("submitFinish :: pragma_deviceExec_time = %le\n", pragma_deviceExec_time[0]);

    pragma_issue_time_total += (pragma_timing_issued - pragma_timing_start);
    pragma_wait_time_total += (finishTime - pragma_timing_issued);
    pragma_deviceExec_time_total += pragma_deviceExec_time[0];
    pragma_num_samples++;
    if (pragma_num_samples >= 200) {
      printf("[MIC-TIMING] :: PE:%d, RANK:%d :: issue: %.3lf ms, wait: %.3lf ms, deviceExec:%.3lf\n",
             pe, rank,
             (pragma_issue_time_total / pragma_num_samples) * 1.0e3,
             (pragma_wait_time_total / pragma_num_samples) * 1.0e3,
             (pragma_deviceExec_time_total / pragma_num_samples) * 1.0e3
            );
      pragma_issue_time_total = 0.0;
      pragma_wait_time_total = 0.0;
      pragma_deviceExec_time_total = 0.0;
      pragma_num_samples = 0;
    }
  }

#else
  #define MIC_PRAGMA_TIMING_DEVICE_CLAUSE
#endif


// DMK - TODO : A device version of the die function, which does not have access
//   to the host info.  Only to be used within the kernel itself.  Perhaps add
//   some startup code to send host/PE# and device# to device and then target
//   mic_die to the device as well, removing the need for this separate function.
__attribute__((target(mic)))
void mic_dev_die(const char * const str) {
  #ifdef __MIC__
    const char * const loc = "on device";
  #else
    const char * const loc = "on host";
  #endif
  if (str != NULL) {
    printf("[MIC_DIE] :: \"%s\" (%s)\n", str, loc);
  } else {
    printf("[MIC_DIE] :: mic_dev_die called (%s)\n", loc);
  }
  fflush(NULL);
  abort();
}


__attribute__((target(mic)))
void mic_print_config() {

  #if MIC_PRINT_CONFIG != 0

    // DMK - TODO | FIXME : Create a mechanism, so that it is only printed once if
    //   there are multiple MIC devices being used.
  
    printf("device :: MULTIPLE_THREADS  (%d)\n", MULTIPLE_THREADS);
    printf("device :: # OpenMP Threads : %d\n", device__numOMPThreads);

    printf("device :: USE_ASSUME          (%d)\n", USE_ASSUME);
    printf("device :: USE_ASSUME_ALIGNED  (%d)\n", USE_ASSUME_ALIGNED);
    printf("device :: USE_RESTRICT        (%d)\n", USE_RESTRICT);

    printf("device :: MIC_HANDCODE_FORCE                     (%d)\n", MIC_HANDCODE_FORCE);
    printf("device ::   MIC_HANDCODE_FORCE_PFDIST            (%d)\n", MIC_HANDCODE_FORCE_PFDIST);
    printf("device ::   MIC_HANDCODE_FORCE_USEGATHER_NBTBL   (%d)\n", MIC_HANDCODE_FORCE_USEGATHER_NBTBL);
    printf("device ::   MIC_HANDCODE_FORCE_USEGATHER         (%d)\n", MIC_HANDCODE_FORCE_USEGATHER);
    printf("device ::   MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT  (%d)\n", MIC_HANDCODE_FORCE_LOAD_VDW_UPFRONT);
    printf("device ::   MIC_HANDCODE_FORCE_COMBINE_FORCES    (%d)\n", MIC_HANDCODE_FORCE_COMBINE_FORCES);
    printf("device ::   MIC_HANDCODE_FORCE_TRANS_AS_MADD     (%d)\n", MIC_HANDCODE_FORCE_TRANS_AS_MADD);

    printf("device :: MIC_HANDCODE_FORCE_SINGLE            (%d)\n", MIC_HANDCODE_FORCE_SINGLE);
    printf("device :: MIC_HANDCODE_FORCE_CALCR2TABLE       (%d)\n", MIC_HANDCODE_FORCE_CALCR2TABLE);
    printf("device :: MIC_HANDCODE_FORCE_SOA_VS_AOS        (%d)\n", MIC_HANDCODE_FORCE_SOA_VS_AOS);

    printf("device :: MIC_SPLIT_WITH_HOST  (%d)\n", MIC_SPLIT_WITH_HOST);

    printf("device :: MIC_SORT_ATOMS  (%d)\n", MIC_SORT_ATOMS);
    printf("device :: MIC_SORT_LISTS  (%d)\n", MIC_SORT_LISTS);

    printf("device :: MIC_HANDCODE_PLGEN    (%d)\n", MIC_HANDCODE_PLGEN);
    printf("device :: MIC_TILE_PLGEN        (%d)\n", MIC_TILE_PLGEN);
    printf("device :: MIC_CONDITION_NORMAL  (%d)\n", MIC_CONDITION_NORMAL);
    printf("device :: MIC_PAD_PLGEN         (%d)\n", MIC_PAD_PLGEN);

    printf("device :: MULTIPLE_THREADS  (%d)\n", MULTIPLE_THREADS);

    printf("device :: MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING  (%d)\n", MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING);
    printf("device ::   MIC_SPECIFIC_COMPUTE_PARTITIONING__SELF_P1  (%d)\n", MIC_SPECIFIC_COMPUTE_PARTITIONING__SELF_P1);
    printf("device ::   MIC_SPECIFIC_COMPUTE_PARTITIONING__PAIR_P1  (%d)\n", MIC_SPECIFIC_COMPUTE_PARTITIONING__PAIR_P1);
    printf("device ::   MIC_SPECIFIC_COMPUTE_PARTITIONING__PAIR_P2  (%d)\n", MIC_SPECIFIC_COMPUTE_PARTITIONING__PAIR_P2);

    printf("device :: MIC_MULTIPLE_KERNELS  (%d)\n", MIC_MULTIPLE_KERNELS);

    printf("device :: MIC_SORT_COMPUTES  (%d)\n", MIC_SORT_COMPUTES);

    printf("device :: MIC_VIRIAL_ENERGY_ALT  (%d)\n", MIC_VIRIAL_ENERGY_ALT);

    printf("device :: REFINE_PAIRLISTS            (%d)\n", REFINE_PAIRLISTS);
    printf("device ::   REFINE_PAIRLIST_HANDCODE  (%d)\n", REFINE_PAIRLIST_HANDCODE);
    printf("device ::   REFINE_PAIRLISTS_XYZ      (%d)\n", REFINE_PAIRLISTS_XYZ);

    printf("device :: MIC_FULL_CHECK     (%d)\n", MIC_FULL_CHECK);
    printf("device :: MIC_EXCL_CHECKSUM  (%d)\n", MIC_EXCL_CHECKSUM);

    printf("device :: MIC_ALIGN  (%d)\n", MIC_ALIGN);

    fflush(NULL);

  #endif // MIC_PRINT_CONFIG
}


// Function to initialize a given MIC device by initializing various device
//   variables, arrays, etc.
// Input:
//  - pe : the PE number for the host core associated with the device
//  - deviceNum : the device number to be initialized
void mic_init_device(const int pe, const int node, const int deviceNum) {

  // Record the PE associated with the device
  device__myPE = pe;
  device__myNode = node;

  // Initialize the various device variables
  device__patch_pairs = NULL;
  device__patch_pairs_size = 0;
  device__patch_pairs_alloc_size = 0;
  device__pairlists = 0; // NULL / Unset
  device__pairlists_alloc_size = 0;
  device__force_lists = NULL;
  device__force_lists_size = 0;
  device__force_lists_alloc_size = 0;
  device__atoms = NULL;
  device__atom_params = NULL;
  device__forces = NULL;
  device__slow_forces = NULL;
  device__atoms_alloc_size = 0;
  device__force_buffers = NULL;
  device__slow_force_buffers = NULL;
  device__force_buffers_req_size = 0;
  device__force_buffers_alloc_size = 0;

  device__remote_kernel_data = new mic_kernel_data[1];
  device__local_kernel_data = new mic_kernel_data[1];

  // DMK - DEBUG - Device timestep counter (only used on device, so initialize to invalid value on host)
  device__timestep = 0;

  // DMK - NOTE: These values are only used on the device, so are just
  //   set to NULL on the host
  device__pl_array = NULL;
  device__pl_size = NULL;
  device__r2_array = NULL;

  // Initialize flags for offload buffers that are not copied every timestep
  patch_pairs_copySize = 0;
  force_lists_copySize = 0;
  atom_params_copySize = 0;

  // DMK - DEBUG
  device__loop_counts = NULL;
  device__cutoff_clear_counts = NULL;
  device__iTest_counts = NULL;

  device__virial_energy = NULL;

  device__numOMPThreads = 0;

  // DMK - TIMING
  #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
    #if MIC_DEVICE_TRACING_DETAILED != 0
      device__device_times_computes = NULL;
      device__device_times_patches = NULL;
    #endif
    device__device_times_start = NULL;
  #endif

  // Initialize the device itself via an offload pragma section
  #pragma offload target(mic:deviceNum) \
    in(device__myNode) in(device__myPE) in(deviceNum) \
    in(device__remote_kernel_data[0:1] : alloc_if(1) free_if(0)) \
    in(device__local_kernel_data[0:1] : alloc_if(1) free_if(0)) \
    nocopy(device__pairlists) nocopy(device__pairlists_alloc_size) \
    nocopy(device__pl_array) nocopy(device__pl_size) nocopy(device__r2_array) \
    nocopy(device__timestep) \
    nocopy(device__loop_counts) nocopy(device__cutoff_clear_counts) nocopy(device__iTest_counts) \
    nocopy(device__force_buffers_req_size) nocopy(device__force_buffers_alloc_size) \
    nocopy(device__virial_energy) \
    nocopy(device__force_buffers) nocopy(device__slow_force_buffers) 
  {
    #if !defined(__MIC__)
      mic_dev_die("mic_init_device :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
    #endif

    // Initialize the force buffer arrays and sizes
    device__force_buffers_alloc_size = 0;
    device__force_buffers_req_size = 0;
    device__force_buffers = NULL;
    device__slow_force_buffers = NULL;

    // Initialize the pairlist pointer array (an array of pointers, one entry per pairlist)
    device__pairlists = NULL;
    device__pairlists_alloc_size = 0;

    // Get the number of threads available to the code
    const int numThreads = omp_get_max_threads();
    #if MULTIPLE_THREADS != 0
      device__numOMPThreads = numThreads;
    #else
      device__numOMPThreads = 1;
    #endif

    // Initialize the r2 arrays (scratch buffers used by computes to refine pairlists
    //   each timestep, when pairlist refinement is enabled)
    #if REFINE_PAIRLISTS != 0
      device__pl_array = new int*[numThreads];     __ASSERT(device__pl_array != NULL);
      device__pl_size = new int[numThreads];       __ASSERT(device__pl_size != NULL);
      device__r2_array = new double*[numThreads];  __ASSERT(device__r2_array != NULL);
      for (int i = 0; i < numThreads; i++) {
        device__pl_array[i] = NULL;
        device__pl_size[i] = 0;
        device__r2_array[i] = NULL;
      }
    #else
      device__pl_array = NULL;
      device__pl_size = 0;
      device__r2_array = NULL;
    #endif

    // DMK - DEBUG - Initialize a counter to track timestep count (for debug output reference)
    device__timestep = 0;

    // DMK - DEBUG
    #if MIC_STATS_LOOP_COUNTS != 0
      int lc_size = 6 * MIC_STATS_LOOP_COUNTS_NUM_BINS;
      device__loop_counts = new int[lc_size];  __ASSERT(device__loop_counts != NULL);
      device__cutoff_clear_counts = new int[lc_size];  __ASSERT(device__cutoff_clear_counts != NULL);
      device__iTest_counts = new int[lc_size];  __ASSERT(device__iTest_counts != NULL);
      memset(device__loop_counts, 0, sizeof(int) * lc_size);
      memset(device__cutoff_clear_counts, 0, sizeof(int) * lc_size);
      memset(device__iTest_counts, 0, sizeof(int) * lc_size);
    #else
      device__loop_counts = NULL;
      device__cutoff_clear_counts = NULL;
      device__iTest_counts = NULL;
    #endif

    #if MIC_VIRIAL_ENERGY_ALT != 0
      device__virial_energy = NULL;
    #endif

    // DMK : NOTE | TODO | FIXME - This will print the configuration of the MICs, but the condition
    //   assumes that there will be at least one MIC on the first node (node 0).  Further, it assumes
    //   that all MICs will be the same.  If there are cases where these assumptions do not hold, we
    //   need to expand this code to cover those cases.  For now, leaving it this way because it
    //   reduces the amount of output during the run (espeically when scaling to multiple nodes and/or
    //   multiple MICs per node).
    if (node == 0 && deviceNum == 0) { mic_print_config(); }

  } // end pragma offload

  #if MIC_PRAGMA_TIMING_STATS != 0
    #pragma offload target(mic:deviceNum) \
      inout(pragma_deviceExec_time[0:1] : alloc_if(1) free_if(0) align(64))
    {
      pragma_deviceExec_time[0] = 0.0;
    }
  #endif
}


// Function that returns 0 if executed on the host, 1 if executed on the MIC device.  Used to
//   indicate if a target is available or not during application startup.
// Input: N/A
// Output:
//  - 0 if executed on host, 1 if executed on device
__attribute__((target(mic))) int mic_check_internal() {
    int retval;
    #ifdef __MIC__
        retval = 1;
    #else
        retval = 0;
    #endif
    return retval;
}


// Function to check that the given device is available for offloading.
// Input:
//  - dev : The device number (0+) to check
// Output:
//  - 1 if device is available (will call mic_dev_die if not available to prevent further application progress)
int mic_check(int dev) {
  int dev_ok = 0;
  #pragma offload target(mic: dev) inout(dev_ok)
  {
    #if !defined(__MIC__)
      mic_dev_die("mic_check :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
    #endif
    dev_ok = mic_check_internal();
  }
  return dev_ok;
}


// Function called during startup to push the table_four data (lookup table used
//   during force computation) that has been calculated on the host down to the
//   given device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - table_four : pointer to the table data itself
//   - table_four_n : dimension of the table_four table (used for indexing and
//       size calculations)
// Output: N/A
void mic_bind_table_four(const int deviceNum,
                         const double *table_four,
                         const int table_four_n,
                         const int table_four_n_16
                        ) {

  // NOTE: This function should only be called once per device at startup,
  //   so check for multiple calls
  if (device__table_four != NULL) {
    printf("[PE:%d,DEV:%d] :: ERROR :: !!! Attempt to write device__table_four multiple times !!!\n", device__myPE, deviceNum);
    fflush(NULL);
    return;
  }

  // Copy the table pointer and dimension information into the device variables.
  //   Note that there are actually several sub-tables within the table_four
  //   table, each a multiple of table_four_n elements in length.
  device__table_four = (double*)table_four;
  device__table_four_float = NULL;
  device__table_four_n_16 = table_four_n_16;
  int numTableFourElems = 61 * table_four_n_16;  // WARNING !!! Must match ComputeNonbondedUtil.C (const 61) !!!

  // Transfer the table_four data and dimension data to the given device
  #pragma offload target(mic:deviceNum) \
    in(device__table_four_n_16) \
    in(device__table_four[0:numTableFourElems] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
    nocopy(device__table_four_float) in(numTableFourElems)
  {
    #if !defined(__MIC__)
      mic_dev_die("mic_bind_table_four :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
    #endif

    //#if (MIC_HANDCODE_FORCE != 0) && (MIC_HANDCODE_FORCE_SINGLE != 0)
    #if MIC_HANDCODE_FORCE_SINGLE != 0
      //device__table_four_float = (float*)(_mm_malloc(sizeof(float) * numTableFourElems, 64));
      device__table_four_float = (float*)(_MM_MALLOC_WRAPPER(sizeof(float) * numTableFourElems, 64, "device__table_four_float"));
      assert(device__table_four_float != NULL);
      for (int i = 0; i < numTableFourElems; i++) {
        device__table_four_float[i] = (float)(device__table_four[i]);
      }
    #else
      device__table_four_float = NULL;
    #endif
  }
}


// Function called during startup to push the lj_table data (lookup table used
//   during force computation) that has been calculated on the host down to the
//   given device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - lj_table : pointer to the table data itself
//   - lj_table_dim : dimension of the lj_table (used for indexing)
//   - lj_table_size : size of the lj_table
// Output: N/A
void mic_bind_lj_table(const int deviceNum,
                       const char *lj_table,
                       const int lj_table_dim,
                       const int lj_table_size
                      ) {

  // NOTE: This function should only be called once per device at startup,
  //   so check for multiple calls
  if (device__lj_table != NULL) {
    printf("[PE:%d,DEV:%d] :: ERROR :: !!! Attempt to write device__lj_table multiple times !!!\n", device__myPE, deviceNum);
    fflush(NULL);
    return;
  }

  // Copy the table pointer, dimension, and size information into the device variables.
  device__lj_table = (char*)lj_table;
  device__lj_table_float = NULL;
  device__lj_table_dim = lj_table_dim;
  device__lj_table_size = lj_table_size;

  // Transfer the table data, dimension, and size information to the given device.
  #pragma offload target(mic:deviceNum) \
    in(device__lj_table_dim) \
    in(device__lj_table_size) \
    in(device__lj_table[0:lj_table_size] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
    nocopy(device__lj_table_float)
  {
    #if !defined(__MIC__)
      mic_dev_die("mic_bind_lj_table :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
    #endif

    // If single precision is being used, the convert the entries of the table
    //   from double to single
    //#if (MIC_HANDCODE_FORCE != 0) && (MIC_HANDCODE_FORCE_SINGLE != 0)
    #if MIC_HANDCODE_FORCE_SINGLE != 0
      int numElements = device__lj_table_size * sizeof(char) / sizeof(double);
      //device__lj_table_float = (float*)(_mm_malloc(sizeof(float) * numElements, 64));
      device__lj_table_float = (float*)(_MM_MALLOC_WRAPPER(sizeof(float) * numElements, 64, "device__lj_table_float"));
      assert(device__lj_table_float != NULL);
      double *srcPtr = (double*)(device__lj_table);
      for (int i = 0; i < numElements; i++) {
        device__lj_table_float[i] = (float)(srcPtr[i]);
      }
    #else
      device__lj_table_float = NULL;
    #endif
  }
}


// Function called during startup to push the exclusion data (lookup table used
//   during pairlist generation) that has been calculated on the host down to the
//   given device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - exclusion_bits : pointer to the exclusion data itself
//   - exclusion_bits_size : size of the exclsion data
// Output: N/A
void mic_bind_exclusions(const int deviceNum,
                         unsigned int *exclusion_bits,
                         const long int exclusion_bits_size
                        ) {

  // NOTE: This function should only be called once per device at startup,
  //   so check for multiple calls
  if (device__exclusion_bits != NULL) {
    printf("[PE:%d,DEV:%d] :: ERROR :: !!! Attempt to write device__exclusion_bits multiple times !!!\n", device__myPE, deviceNum);
    fflush(NULL);
    return;
  }

  // Copy the exclusion pointer and size information in to the device variables.
  device__exclusion_bits = exclusion_bits;
  device__exclusion_bits_size = exclusion_bits_size;

  // Transfer the exclusion data and size information down to the given device.
  #pragma offload target(mic:deviceNum) \
    in(device__exclusion_bits_size) \
    in(device__exclusion_bits[0:exclusion_bits_size] : alloc_if(1) free_if(0) align(MIC_ALIGN))
  {
    #if !defined(__MIC__)
      mic_dev_die("mic_bind_exclusions :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
    #endif
  }
}


// Function called during startup to push (and calculate) several variables that
//   are constant and used throughout the simulation.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   ... : remaining parameters are simply simulation constants
// Output: N/A
void mic_bind_constants(const int deviceNum,
                        const double cutoff2,
                        const double dielectric_1,
                        const double scaling,
                        const double scale14,
                        const double r2_delta,
                        const int r2_delta_exp,
                        const int commOnly
                       ) {
  
  // NOTE: This function should only be called once per device at startup,
  //   so check for multiple calls
  if (device__constants != NULL) {
    printf("[PE:%d,DEV:%d] :: ERROR :: !!! Attempt to write device__constants multiple times !!!\n", device__myPE, deviceNum);
    fflush(NULL);
    return;
  }

  // Create a mic_constants data structure and copy (or calculate) the
  //   various constants that are to be pushed to the device in to this
  //   data structure.
  device__constants = new mic_constants[1];
  device__constants->cutoff2 = cutoff2;
  device__constants->dielectric_1 = dielectric_1;
  device__constants->scaling = scaling;
  device__constants->scale14 = scale14;
  device__constants->modf_mod = 1.0 - scale14;
  device__constants->r2_delta = r2_delta;
  device__constants->r2_delta_exp = r2_delta_exp;
  device__constants->r2_delta_expc = 64 * (r2_delta_exp - 1023);
  device__constants->commOnly = commOnly;

  // Transfer the constants to the given device.
  #pragma offload target(mic:deviceNum) \
    in(device__constants[0:1] : alloc_if(1) free_if(0) align(MIC_ALIGN))
  {
    #if !defined(__MIC__)
      mic_dev_die("mic_bind_constants :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
    #endif

    // Correct r2_delta_expc for use with single precision if using single precision on the device
    //#if (MIC_HANDCODE_FORCE != 0) && (MIC_HANDCODE_FORCE_SINGLE != 0)
    #if MIC_HANDCODE_FORCE_SINGLE != 0
      device__constants->r2_delta_expc = 64 * (device__constants->r2_delta_exp - 127);
    #endif
  }
}


// Function called after patch pairs have been modified on the host to push those
//   modifications to the MIC device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - patch_pairs : array containing individual patch_pair records (data to copy)
//   - patch_pairs_size : the length (in elements) of the patch_pairs array (valid/used elements)
//   - patch_pairs_bufSize : the actual length (in elements) of the patch_pairs array (allocated)
// Output: N/A
void mic_bind_patch_pairs_only(const int deviceNum,
                               const patch_pair *patch_pairs,
                               const int patch_pairs_size,
                               const int patch_pairs_bufSize
			      ) {

  //// DMK - DEBUG
  //printf("DEBUG 0.0 - mic_bind_patch_pairs_only() called...\n");
  //fflush(NULL);

  // Validate the parameters
  #if MIC_FULL_CHECK != 0
    assert(patch_pairs_bufSize >= patch_pairs_size);
  #endif

  // Check if the buffer currently allocated on the device is too small.  If so,
  //   reallocate the buffer to be larger.
  if (patch_pairs_bufSize > device__patch_pairs_alloc_size) {  // If buffer is too small, reallocate

    // Free the old buffer
    if (device__patch_pairs != NULL) {
      #pragma offload target(mic:deviceNum) \
        in(device__patch_pairs[0:0] : alloc_if(0) free_if(1))
      {
        #if !defined(__MIC__)
          mic_dev_die("mic_bind_patch_pairs (0) :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
        #endif
      }
    }

    // Record the new info
    device__patch_pairs = patch_pairs;
    device__patch_pairs_size = patch_pairs_size;
    device__patch_pairs_alloc_size = patch_pairs_bufSize;

    // Allocate the new memory
    // NOTE: This will allocate an amount of memory equal to the amount allocated
    //   on the host (i.e. if the host over allocates expecting future growth, then
    //   the device will match that over allocation).
    #pragma offload target(mic:deviceNum) \
      in(device__patch_pairs_size) \
      in(device__patch_pairs_alloc_size) \
      in(device__patch_pairs[0:patch_pairs_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      nocopy(device__virial_energy)
    {
      #if !defined(__MIC__)
        mic_dev_die("mic_bind_patch_pairs (1) :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
      #endif

      // Make sure there are enough pairlists pointers (one per pairlist type per patch_pair)
      const int numPairlists = NUM_PAIRLIST_TYPES * device__patch_pairs_size;
      if (numPairlists > device__pairlists_alloc_size) {
        int **new_pairlists = new int*[numPairlists];
        int initStart = 0;
        if (device__pairlists != 0) {
          int **old_pairlists = (int**)(device__pairlists);
          memcpy((void*)new_pairlists, (void*)old_pairlists, sizeof(int*) * device__pairlists_alloc_size);
          delete [] old_pairlists;
          initStart = device__pairlists_alloc_size;
        }
        for (int i = initStart; i < numPairlists; i++) { new_pairlists[i] = NULL; }
        device__pairlists = (unsigned long long int)new_pairlists;
        device__pairlists_alloc_size = numPairlists;
      }

      // Grow the virial/energy buffer, if need be
      #if MIC_VIRIAL_ENERGY_ALT != 0
        if (device__virial_energy != NULL) { _mm_free(device__virial_energy); }
        int virialEnergySize = 16 * (device__numOMPThreads + device__patch_pairs_alloc_size);
        //device__virial_energy = (double*)(_mm_malloc(virialEnergySize * sizeof(double), 64));
        device__virial_energy = (double*)(_MM_MALLOC_WRAPPER(virialEnergySize * sizeof(double), 64, "device__virial_energy"));
      #endif
    }

    // DMK - TIMING
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)

      // Free old buffers if need be
      if (device__device_times_start != NULL) {
        #if MIC_DEVICE_TRACING_DETAILED != 0
          #pragma offload target(mic:deviceNum) \
            out(device__device_times_computes[0:0] : alloc_if(0) free_if(1)) \
            out(device__device_times_start[0:0] : alloc_if(0) free_if(1))
          { }
        #else
          #pragma offload target(mic:deviceNum) \
            out(device__device_times_start[0:0] : alloc_if(0) free_if(1))
          { }
        #endif
      }

      //device__device_times_computes = (double*)_mm_malloc(patch_pairs_bufSize * 2 * sizeof(double), 64);
      //device__device_times_start = (double*)_mm_malloc(2 * 2 * sizeof(double), 64);
      #if MIC_DEVICE_TRACING_DETAILED != 0
        device__device_times_computes = (double*)_MM_MALLOC_WRAPPER(patch_pairs_bufSize * 2 * sizeof(double), 64, "device__device_times_computes");
      #endif
      device__device_times_start = (double*)_MM_MALLOC_WRAPPER(10 * 2 * sizeof(double), 64, "device__device_times_start");

      // Allocate new buffers
      #if MIC_DEVICE_TRACING_DETAILED != 0
        #pragma offload target(mic:deviceNum) \
          in(device__device_times_computes[0:2*patch_pairs_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
          in(device__device_times_start[0:10] : alloc_if(1) free_if(0) align(MIC_ALIGN))
        { }
      #else
        #pragma offload target(mic:deviceNum) \
          in(device__device_times_start[0:10] : alloc_if(1) free_if(0) align(MIC_ALIGN))
        { }
      #endif

    #endif
  }

  // Record the current 'used' length of the array and flag the array for data transfer
  device__patch_pairs_size = patch_pairs_size;
  patch_pairs_copySize = patch_pairs_size;
}

// Function called after force lists have been modified on the host to push those
//   modifications to the MIC device
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - force_lists : array containing individual force_list records (data to copy)
//   - force_lists_size : the length (in elements) of the force_lists array (valid/used elements)
//   - force_lists_bufSize : the actual length (in elements) of the force_lists array (allocated)
// Output: N/A
void mic_bind_force_lists_only(const int deviceNum,
                               force_list *force_lists,
                               const int force_lists_size,
                               const int force_lists_bufSize
                              ) {

  // Validate the parameters
  #if MIC_FULL_CHECK != 0
    assert(force_lists_bufSize >= force_lists_size);
  #endif

  // Check if the buffer currently allocated on the device is too small.  If so,
  //   reallocate the buffer to be larger.
  if (force_lists_bufSize > device__force_lists_alloc_size) {

    // Free the old buffer
    if (device__force_lists != NULL) {
      #pragma offload target(mic:deviceNum) in(device__force_lists[0:0] : alloc_if(0) free_if(1))
      {
        #if !defined(__MIC__)
          mic_dev_die("mic_bind_patch_pairs (2) :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
        #endif
      }
    }

    // Record the new info
    device__force_lists = force_lists;
    device__force_lists_size = force_lists_size;
    device__force_lists_alloc_size = force_lists_bufSize; //(int)(1.2 * force_lists_size);

    // Allocate the new memory
    // NOTE: This will allocate an amount of memory equal to the amount allocated
    //   on the host (i.e. if the host over allocates expecting future growth, then
    //   the device will match that over allocation).
    #pragma offload target(mic:deviceNum) \
      in(device__force_lists_size) \
      in(device__force_lists_alloc_size) \
      in(device__force_lists[0:force_lists_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN))
    {
      #if !defined(__MIC__)
        mic_dev_die("mic_bind_patch_pairs (3) :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
      #endif
    }

    // DMK - TIMING
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)

      // Free old buffers if need be
      if (device__device_times_patches != NULL) {
        #pragma offload target(mic:deviceNum) \
          out(device__device_times_patches[0:0] : alloc_if(0) free_if(1))
        { }
      }

      //device__device_times_patches = (double*)_mm_malloc(force_lists_bufSize * 2 * sizeof(double), 64);
      device__device_times_patches = (double*)_MM_MALLOC_WRAPPER(force_lists_bufSize * 2 * sizeof(double), 64, "device__device_times_patches");

      // Allocate new buffers
      #pragma offload target(mic:deviceNum) \
        in(device__device_times_patches[0:2*force_lists_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN))
      { }

    #endif
  }

  // Record the current 'used' length of the array and flag the array for data transfer
  device__force_lists_size = force_lists_size;
  force_lists_copySize = force_lists_size;
}

// Function called after the atom list has been modified on the host to push those
//   modifications to the MIC device.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - atoms : array containing individual atom records (data to copy)
//   - atom_params : array containing individual atom_param records (data to copy)
//   - forces : array that will contain output data from the device (data to copy back)
//   - slow_forces : array that will contain output data from the device (data to copy back)
//   - atoms_size : the length (in elements) of the specified arrays (valid/used elements)
//   - atoms_bufSize : the actual length (in elements) of the specified arrays (allocated)
// Output: N/A
void mic_bind_atoms_only(const int deviceNum,
                         atom *atoms,
                         atom_param *atom_params,
                         double4 *forces,
                         double4 *slow_forces,
                         const int atoms_size,
                         const int atoms_bufSize
			) {

  // Validate the parameters
  #if MIC_FULL_CHECK != 0
    assert(atoms_bufSize >= atoms_size);
  #endif

  // Check if the buffers currently allocated on the device are too small.  If so,
  //   reallocate the buffers to be larger.
  if (atoms_bufSize > device__atoms_alloc_size) {

    // Unlock the old buffers
    #if MIC_LOCK_MEMORY
      if (device__atoms != NULL) {
        if (0 != munlock(device__atoms, device__atoms_alloc_size * sizeof(atom))) { printf("[MIC-WARNING] :: failed to munlock device__atoms\n"); }
        if (0 != munlock(device__atom_params, device__atoms_alloc_size * sizeof(atom_param))) { printf("[MIC-WARNING] :: failed to munlock device__atom_params\n"); }
        if (0 != munlock(device__forces, device__atoms_alloc_size * sizeof(double4))) { printf("[MIC-WARNING] :: failed to munlock device__forces\n"); }
        if (0 != munlock(device__slow_forces, device__atoms_alloc_size * sizeof(double4))) { printf("[MIC-WARNING] :: failed to munlock device__slow_forces\n"); }
      }
    #endif

    // Free the old buffers
    if (device__atoms != NULL) {
      #pragma offload target(mic:deviceNum) \
        in(device__atoms[0:0] : alloc_if(0) free_if(1)) \
        in(device__atom_params[0:0] : alloc_if(0) free_if(1)) \
        in(device__forces[0:0] : alloc_if(0) free_if(1)) \
        in(device__slow_forces[0:0] : alloc_if(0) free_if(1))
      {
        #if !defined(__MIC__)
          mic_dev_die("mic_bind_patch_pairs (4) :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
        #endif
      }
    }

    // Record the new info
    device__atoms = atoms;
    device__atom_params = atom_params;
    device__forces = forces;
    device__slow_forces = slow_forces;
    device__atoms_size = atoms_size;
    device__atoms_alloc_size = atoms_bufSize;

    // Lock the new buffers
    #if MIC_LOCK_MEMORY
      if (0 != mlock(device__atoms, device__atoms_alloc_size * sizeof(atom))) { perror("[MIC-WARNING] :: failed to mlock device__atoms - "); }
      if (0 != mlock(device__atom_params, device__atoms_alloc_size * sizeof(atom_param))) { printf("[MIC-WARNING] :: failed to mlock device__atom_params\n"); }
      if (0 != mlock(device__forces, device__atoms_alloc_size * sizeof(double4))) { printf("[MIC-WARNING] :: failed to mlock device__forces\n"); }
      if (0 != mlock(device__slow_forces, device__atoms_alloc_size * sizeof(double4))) { printf("[MIC-WARNING] :: failed to mlock device__slow_forces\n"); }
    #endif

    // Allocate the new memory on the device
    // NOTE: This will allocate an amount of memory for each buffer equal to the ammount
    //   allocated on the host (i.e. if the host over allocates expecting future growth,
    //   then the device will match that over allocation).
    #pragma offload target(mic:deviceNum) \
      in(device__atoms_size) \
      in(device__atoms_alloc_size) \
      in(device__atoms[0:atoms_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      in(device__atom_params[0:atoms_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      in(device__forces[0:atoms_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      in(device__slow_forces[0:atoms_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN))
    {
      #if !defined(__MIC__)
        mic_dev_die("mic_bind_patch_pairs (5) :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
      #endif
    }
  }

  // Record the current 'used' length of the arrays and flag those arrays for data transfer
  // NOTE: The 'atoms', 'forces', and 'slow_forces' arrays will be copied every timestep, as
  //   the data within them changes every timestep.  The 'atom_params' array is only copied
  //   periodically (when atoms migrate).  As such, only flag the transfer of atom_params here.
  device__atoms_size = atoms_size;
  atom_params_copySize = atoms_size;
}

// Function called to allocate memory on the MIC device related to the individual and private
//   force buffers for each of the compute objects.  Note that this function simply allocates
//   the memory on the device, as there is no host-side equivalent buffer (i.e. device only).
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - force_buffers_size : the amount of memory to be allocated on the device for the buffer
void mic_bind_force_buffers_only(const int deviceNum,
                                 const int force_buffers_size
                                ) {

  if (device__force_buffers_req_size < force_buffers_size) {
    device__force_buffers_req_size = (int)(force_buffers_size * 1.2f);
    device__force_buffers_req_size = (device__force_buffers_req_size + 4095) & (~4095);  // Round up to 4K

    //// DMK - DEBUG
    //printf("[MIC-HOST] :: adjusting device__force_buffers_req_size (%d)\n", device__force_buffers_req_size);
  }

  /*
  // Check if the buffers currently allocated on the device are too small.  If so,
  //   reallocate the buffers to be larger.
  if (force_buffers_size > device__force_buffers_alloc_size) {

    // Free the old buffers
    if (device__force_buffers != NULL) {
      #pragma offload target(mic:deviceNum) \
        in(device__force_buffers[0:0] : alloc_if(0) free_if(1)) \
        in(device__slow_force_buffers[0:0] : alloc_if(0) free_if(1))
      {
        #if !defined(__MIC__)
          mic_dev_die("mic_bind_patch_pairs (6) :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
        #endif
      }
      delete [] device__force_buffers;
      delete [] device__slow_force_buffers;
    }

    // Recod the new info
    int force_buffers_bufSize = (int)(force_buffers_size * 1.1);
    device__force_buffers = new double4[force_buffers_bufSize];
    device__slow_force_buffers = new double4[force_buffers_bufSize];
    device__force_buffers_alloc_size = force_buffers_bufSize;

    // Allocate the new memory
    #pragma offload target(mic:deviceNum) \
      in(device__force_buffers_alloc_size) \
      in(device__force_buffers[0:force_buffers_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN)) \
      in(device__slow_force_buffers[0:force_buffers_bufSize] : alloc_if(1) free_if(0) align(MIC_ALIGN))
    {
      #if !defined(__MIC__)
        mic_dev_die("mic_bind_patch_pairs (7) :: Attempt to execute offload on host (not supported yet)... exiting."); fflush(NULL);
      #endif
    }
  }
  */
}


#if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0

void mic_submit_patch_data(const int deviceNum,
                           void * const hostPtr,
                           void* &hostPtr_prev,
                           const int transferBytes,
                           const int allocBytes_host,
                           int &allocBytes_device,
                           uint64_t &devicePtr
                          ) {

  //// DMK - DEBUG
  //printf("[DEBUG] :: mic_submit_patch_data(dN:%d, hP:%p, hPp:%p, tB:%d, aBh:%d, aBd:%d, dP:%ld) - Called...\n",
  //       deviceNum, hostPtr, hostPtr_prev, transferBytes, allocBytes_host, allocBytes_device, devicePtr
  //      );

  // NOTE: This function is called once per patch per timestep by one of the host threads
  // TODO : Since signals are thread specific, cannot use them here since the master does
  //   the checking for signal completions.  Figure out a way to make this asynchronous,
  //   but just do synchronous transfer for now to see what happens.

  // Check to see if the device's copy of the buffer needs to be reallocated
  // NOTE: If we are in a position where the buffer might grow in size, the contents
  //   need to be rewritten in full (i.e. atom migration has occured), so don't worry
  //   about the previous contents of the buffer
  int allocFlag = 0;
  if (hostPtr != hostPtr_prev || allocBytes_host > allocBytes_device) {

    // Free the old device buffer if it has been allocated previously
    if (devicePtr != 0) {
      char* __hostPtr_prev = (char*)hostPtr_prev;
      #pragma offload target(mic:deviceNum) \
        out(__hostPtr_prev[0:0] : alloc_if(0) free_if(1))
      { }
    }

    // Flag that allocation is needed
    allocFlag = 1;
  }

  // Transfer the data (allocating if needed)
  char* __hostPtr = (char*)hostPtr;
  uint64_t __devicePtr = 0;
  int __length = transferBytes;
  if (allocFlag != 0) {
    __length = allocBytes_host;
    allocBytes_device = allocBytes_host;
  }
  #pragma offload target(mic:deviceNum) \
    in(__hostPtr[0:__length] : alloc_if(allocFlag) free_if(0) align(MIC_ALIGN)) \
    out(__devicePtr)
  { __devicePtr = (uint64_t)__hostPtr; }

  devicePtr = __devicePtr;
  hostPtr_prev = hostPtr;


  /*
  // Verify the parameters
  __ASSERT(deviceNum >= 0 && hostPtr != NULL && pid >= 0 && transferBytes > 0 && allocBytes > 0);

  // Manage the device's array of pointers to patch data (device__atomData), growing it if it isn't large enough
  if (pid > atomData_allocSize[deviceNum]) {
    pthread_mutex_lock(&(atomData_locks[deviceNum]));
    if (pid > atomData_allocSize[deviceNum]) { // Check again in case someone else had the lock and did it already

      // Free the old device memory (pulling the data)
      if (device__atomData != NULL) {
        #pragma offload_transfer target(mic:deviceNum) \
          out(device__atomData[0:atomData_allocSize[deviceNum]] : alloc_if(0) free_if(1))
	{ }
      }

      // Place old pointers in temp variables
      uint64_t *oldAtomData = device__atomData;
      int *oldAtomDataLen = atomDataLen;
      int oldAtomData_allocSize = atomData_allocSize[deviceNum];

      // Calculate a new allocation size for the arrays
      int newAllocSize = (int)(pid * 1.2f);
      newAllocSize = ((newAllocSize + 127) & (~127)); // Round up to a multiple of 128 (minimizes mallocs if process starts with small pids)

      // Allocate new host memory, copy in the old data, and free the old data
      uint64_t *device__atomData = (uint64_t*)(_mm_malloc(newAllocSize * sizeof(uint64_t), MIC_ALIGN));
      int *atomDataLen = (int*)(_mm_malloc(newAllocSize * sizeof(int), MIC_ALIGN));
      __ASSERT(device__atomData != NULL && atomDataLen != NULL);

      // Copy the old data and initialize new entries
      memcpy(device__atomData, oldAtomData, oldAtomData_allocSize * sizeof(uint64_t));
      for (int i = oldAtomData_allocSize; i < newAllocSize; i++) { device__atomData[i] = 0L; }
      memcpy(atomDataLen, oldAtomDataLen, oldAtomData_allocSize * sizeof(int));
      for (int i = oldAtomData_allocSize; i < newAllocSize; i++) { atomDataLen[i] = 0; }

      // Free the old host memory
      _mm_free(oldAtomData);
      _mm_free(oldAtomDataLen);

      // Record the new allocation size
      atomData_allocSize[deviceNum] = newAllocSize;

      // Allocate the new data on the device
      #pragma offload_transfer target(mic:deviceNum) \
        in(device__atomData[0:newAllocSize] : alloc_if(1) free_if(0) align(MIC_ALIGN))
      { }
    }
    pthread_mutex_unlock(&(atomData_locks[deviceNum]));
  }

  // DMK - Continue here... this approach won't work because hostPtr, etc. aren't declared
  //   with (and cannot be declared with) __declspec(target(mic)).  So... try creating a
  //   global array of pointers on the host, one per atom buffer on the device, and using
  //   each individual pointer in an offload in clause... for example...
  //       in( (device_patchBuffers[patchIndex])[0:numAtoms] : alloc_if(?) free_if(?) ... )

  // Determin if free/alloc needs to occur on the device
  int allocFlag = ((device__atomData[pid] == 0L) || (transferBytes > 
  */



  /*
  // Determine if free/alloc needs to occur on device
  int allocFlag = ((devicePtr == 0L) || (transferBytes > allocBytes) || (hostPtr != hostPtr_prev));
  int freeFlag = ((devicePtr != 0L) && (allocFlag));  // NOTE: need to reallocate and existing device data

  // If reallocation needs to occur, start by deallocating the old data
  // NOTE: This should be an increasinly rare event as the buffer only grows in size
  //   and the length is a function of the number of atoms per patch
  if (freeFlag) {
    #pragma offload_transfer target(mic:deviceNum) \
      out(hostPtr_prev[0:0] : alloc_if(0) free_if(1))
    { }
  }

  // Transfer the specified data to the device
  int hostPtrBytes = transferBytes;
  if (allocFlag) { hostPtrBytes = allocBytes; }
  #pragma offload target(mic:deviceNum) \
    in(hostPtr[0:hostPtrBytes] : alloc_if(allocFlag) free_if(0) align(MIC_ALIGN)) \
    inout(devicePtr) \
    signal(hostPtr)
  { devicePtr = (uint64_t)hostPtr; }

  // Update the "previous" hostPtr to point to the new buffer now that the device's
  //   memory has been reallocated
  hostPtr_prev = hostPtr;
  */
}
#endif  // MIC_SUBMIT_ATOMS_ON_ARRIVAL


// DMK - DEBUG - Branching and lane utilization stats
#if MIC_ACTIVE_CUTOFF_STATS != 0
  __declspec(target(mic)) int activeCount = 0;
  __declspec(target(mic)) int activePossibleCount = 0;
  __declspec(target(mic)) int cutoffCount = 0;
  __declspec(target(mic)) int cutoffPossibleCount = 0;
#endif

// DMK - DEBUG - Unique cacheline access for gather/scatter stats
#if MIC_GATHER_SCATTER_STATS != 0

  __declspec(target(mic)) int gather_single[17] = { 0 };
  __declspec(target(mic)) int gather_double[17] = { 0 };
  __declspec(target(mic)) int scatter_single[17] = { 0 };
  __declspec(target(mic)) int scatter_double[17] = { 0 };

  // NOTE: This can be slow, just collects stats, not used in production
  #define __SUBMIT_GATHER_SCATTER_STATS(type, index_vec, mask_mask, scale_shift, numEntries, str) \
  { \
    int numUnique = 0; \
    int index[16] __attribute__((aligned(64))); \
    _mm512_store_epi32(index, _mm512_srli_epi32(index_vec, 6 - scale_shift)); \
    int mask = _mm512_mask2int(mask_mask); \
    for (int i = 0; i < numEntries; i++) { \
      if (((mask >> i) & 0x1) == 0) { continue; } \
      int isUnique = 1; \
      for (int j = i + 1; j < numEntries; j++) { \
        if (((mask >> j) & 0x1) == 0) { continue; } \
        if (index[i] == index[j]) { isUnique = 0; break; } \
      } \
      numUnique += isUnique; \
    } \
    type[numUnique]++; \
  }

  #define SUBMIT_GATHER_SINGLE_STATS(index_vec, mask_mask) __SUBMIT_GATHER_SCATTER_STATS(gather_single, index_vec, mask_mask, 2, 16, "GS")
  #define SUBMIT_GATHER_DOUBLE_STATS(index_vec, mask_mask) __SUBMIT_GATHER_SCATTER_STATS(gather_double, index_vec, mask_mask, 3, 8, "GD")
  #define SUBMIT_SCATTER_SINGLE_STATS(index_vec, mask_mask) __SUBMIT_SCATTER_SCATTER_STATS(scatter_single, index_vec, mask_mask, 2, 16, "SS")
  #define SUBMIT_SCATTER_DOUBLE_STATS(index_vec, mask_mask) __SUBMIT_SCATTER_SCATTER_STATS(scatter_double, index_vec, mask_mask, 3, 8, "SD")

#else

  #define SUBMIT_GATHER_SINGLE_STATS(index_vec, mask_mask)
  #define SUBMIT_GATHER_DOUBLE_STATS(index_vec, mask_mask)
  #define SUBMIT_SCATTER_SINGLE_STATS(index_vec, mask_mask)
  #define SUBMIT_SCATTER_DOUBLE_STATS(index_vec, mask_mask)

#endif


// Kernel bodies : Include CopmuteNonbondedMICKernelBase.h to declare each
//   version of the force computation functions (calc_pair, calc_self, etc.).
//   Each time the header file is included, different macros are set/unset
//   to control the exact contents of 'that version of the kernel.'

#define NBPAIR 1
#define NBSELF 2

#define NBTYPE NBPAIR
  #include "ComputeNonbondedMICKernelBase.h"
  #define CALCENERGY
    #include "ComputeNonbondedMICKernelBase.h"
  #undef CALCENERGY
  #define FULLELECT
    #include "ComputeNonbondedMICKernelBase.h"
    #define CALCENERGY
      #include "ComputeNonbondedMICKernelBase.h"
    #undef CALCENERGY
  #undef FULLELECT
#undef NBTYPE

#define NBTYPE NBSELF
  #include "ComputeNonbondedMICKernelBase.h"
  #define CALCENERGY
    #include "ComputeNonbondedMICKernelBase.h"
  #undef CALCENERGY
  #define FULLELECT
    #include "ComputeNonbondedMICKernelBase.h"
    #define CALCENERGY
      #include "ComputeNonbondedMICKernelBase.h"
    #undef CALCENERGY
  #undef FULLELECT
#undef NBTYPE


// This function is the main "computation" function, called each timestep to
//   initiate computation on the MIC device.  A signal is setup as part of
//   the offload pragma, which is periodically checked by the
//   mic_check_remote_kernel_complete and mic_check_local_kernel_complete
//   functions (which are called periodically by the Charm++ runtime system),
//   so progress can continue once the device has finished its computation.
// Input:
//   - deviceNum : the device to push the data to (number within node)
//   - isRemote : flag indicating a 'remote' (1) or 'local' (0) kernel invocation (currently, must be 'local')
//   - lata : lattice information ('a' along with latb and latc)
//   - latb : lattice information ('b' along with lata and latc)
//   - latc : lattice information ('c' along with lata and latb)
//   - doSlow : flag indicating if slow forces should be calculated in this kernel invocation
//   - doEnergy : flag indiciating if energy information should be calculated in this kernel invocation
//   - usePairlists : flag indicating if pairlists should be used within this kernel invocation (currently, must be)
//   - savePairlists : flag indicating if pairlists should be saved across kernel invocations (currently, must be)
//   - atomsChanged : flag indicating whether or not the atoms (and related) buffers have changed (length, not content)
// Output: N/A
void mic_nonbonded_forces(const int deviceNum,
                          const int isRemote,
                          const int numLocalAtoms,
                          const int numLocalComputes,
                          const int numLocalPatches,
                          const mic_position3_t lata,
                          const mic_position3_t latb,
                          const mic_position3_t latc,
                          const int doSlow,
                          const int doEnergy,
                          const int usePairlists,
                          const int savePairlists,
                          const int atomsChanged
                         ) {

  //// DMK - DEBUG
  //printf("[DEBUG:NBF-%d] :: mic_nonbonded_forces() - Called...\n", isRemote);
  //printf("[DEBUG:NBF-%d] ::   isRemote:%d\n", isRemote, isRemote);
  //printf("[DEBUG:NBF-%d] ::   numLocalAtoms:%d\n", isRemote, numLocalAtoms);
  //printf("[DEBUG:NBF-%d] ::   numLocalComputes:%d\n", isRemote, numLocalComputes);
  //printf("[DEBUG:NBF-%d] ::   numLocalPatches:%d\n", isRemote, numLocalPatches);
  //fflush(NULL);

  // DMK - NOTE : At the moment, information regarding which patch pairs
  //   are local vs remote is not being passed into this function (just
  //   if it should to local or remote patch pairs).  Add that support
  //   as the kernel is written.

  const int numAtoms = device__atoms_size;
  const int numPatchPairs = device__patch_pairs_size;
  const int numForceLists = device__force_lists_size;

  // Setup the lattice information
  // DMK - TODO | FIXME : See if these variables can be shifted to mic_constants
  device__lata = lata;  
  device__latb = latb;
  device__latc = latc;

  // Get and set the tag to use
  int *tag_kernel = ((isRemote) ? (&tag_remote_kernel) : (&tag_local_kernel));
  *tag_kernel = 1;

  // Setup the kernel data structures
  // NOTE: The remote kernel will be called first, so set them both up at the same time during remote issue
  mic_kernel_data *kernel_data = ((isRemote) ? (device__remote_kernel_data) : (device__local_kernel_data));
  kernel_data->isRemote = isRemote;
  kernel_data->numLocalAtoms = numLocalAtoms;
  kernel_data->numLocalComputes = numLocalComputes;
  kernel_data->numLocalPatches = numLocalPatches;
  kernel_data->doSlow = doSlow;
  kernel_data->doEnergy = doEnergy;
  kernel_data->usePairlists = usePairlists;
  kernel_data->savePairlists = savePairlists;
  kernel_data->virial_xx = 0.0;
  kernel_data->virial_xy = 0.0;
  kernel_data->virial_xz = 0.0;
  kernel_data->virial_yy = 0.0;
  kernel_data->virial_yz = 0.0;
  kernel_data->virial_zz = 0.0;
  kernel_data->fullElectVirial_xx = 0.0;
  kernel_data->fullElectVirial_xy = 0.0;
  kernel_data->fullElectVirial_xz = 0.0;
  kernel_data->fullElectVirial_yy = 0.0;
  kernel_data->fullElectVirial_yz = 0.0;
  kernel_data->fullElectVirial_zz = 0.0;
  kernel_data->vdwEnergy = 0.0;
  kernel_data->electEnergy = 0.0;
  kernel_data->fullElectEnergy = 0.0;
  int remoteKernelDataLen = (isRemote) ? (1) : (0);
  int localKernelDataLen = (isRemote) ? (0) : (1);

  //// DMK - DEBUG
  //printf("[DEBUG] :: kernel_data:%p (r:%p, l:%p)\n", kernel_data, device__remote_kernel_data, device__local_kernel_data);
  //printf("[DEBUG] :: kernelDataLen - r:%d, l:%d\n", remoteKernelDataLen, localKernelDataLen);
  //fflush(NULL);

  //// Select the kernel data struct and initialize member variables for this timestep
  ///mic_kernel_data *kernel_data = ((isRemote) ? (device__remote_kernel_data) : (device__local_kernel_data));
  //device__local_kernel_data->patchPairStart = 0; // DMK - NOTE : 0 for now, must change when/if local vs remote support added
  //device__local_kernel_data->patchPairEnd = 0;   // DMK - NOTE : 0 for now, must change when/if local vs remote support added
  //device__local_kernel_data->doSlow = doSlow;
  //device__local_kernel_data->doEnergy = doEnergy;
  //device__local_kernel_data->usePairlists = usePairlists;
  //device__local_kernel_data->savePairlists = savePairlists;
  //device__local_kernel_data->virial_xx = 0.0;
  //device__local_kernel_data->virial_xy = 0.0;
  //device__local_kernel_data->virial_xz = 0.0;
  //device__local_kernel_data->virial_yy = 0.0;
  //device__local_kernel_data->virial_yz = 0.0;
  //device__local_kernel_data->virial_zz = 0.0;
  //device__local_kernel_data->fullElectVirial_xx = 0.0;
  //device__local_kernel_data->fullElectVirial_xy = 0.0;
  //device__local_kernel_data->fullElectVirial_xz = 0.0;
  //device__local_kernel_data->fullElectVirial_yy = 0.0;
  //device__local_kernel_data->fullElectVirial_yz = 0.0;
  //device__local_kernel_data->fullElectVirial_zz = 0.0;
  //device__local_kernel_data->vdwEnergy = 0.0;
  //device__local_kernel_data->electEnergy = 0.0;
  //device__local_kernel_data->fullElectEnergy = 0.0;

  // For buffers that are only periodically copied/updated, get the
  //   flags that indicate whether or not those buffers should be
  //   transfered during this kernel invocation, and then reset them
  const int slowForcesNumAtoms = ((doSlow) ? (numAtoms) : (0));
  const int _patch_pairs_copySize = patch_pairs_copySize;  // from __thread variable to stack variable for use within offload pragma
  const int _force_lists_copySize = force_lists_copySize;  // from __thread variable to stack variable for use within offload pragma
  const int _atom_params_copySize = atom_params_copySize;  // from __thread variable to stack variable for use within offload pragma
  patch_pairs_copySize = 0;
  force_lists_copySize = 0;
  atom_params_copySize = 0;

  // NOTE: A compute is "remote" if either patch that it depends on is "remote," so some of the "remote"
  //   computes will require "local" patch data (i.e. copy all patch data at start of the "remote" kernel).
  #if MIC_MULTIPLE_KERNELS != 0
    int toCopySize_atoms = (isRemote) ? (numAtoms) : (0);  // The remote kernel will copy in all atoms
    int toCopySize_atom_params = (isRemote) ? (_atom_params_copySize) : (0);  // If there are atom_params to copy, the remote kernel will copy in all atom_params
    int toCopySize_patch_pairs = (isRemote) ? (_patch_pairs_copySize) : (0);  // This could potentially be split between kernels (but is relatively rare, so don't worry about it for now)
    int toCopySize_force_lists = (isRemote) ? (_force_lists_copySize) : (0);  // This could potentially be split between kernels (but is relatively rare, so don't worry about it for now)
    int toCopyStart_forces = (isRemote) ? (numLocalAtoms) : (0);
    int toCopySize_forces = (isRemote) ? (numAtoms - numLocalAtoms) : (numLocalAtoms);
    int toCopyStart_slow_forces = (doSlow) ? (toCopyStart_forces) : (0);
    int toCopySize_slow_forces = (doSlow) ? (toCopySize_forces) : (0);
  #else
    int toCopySize_atoms = numAtoms;
    int toCopySize_atom_params = _atom_params_copySize;
    int toCopySize_patch_pairs = _patch_pairs_copySize;
    int toCopySize_force_lists = _force_lists_copySize;
    int toCopyStart_forces = 0;
    int toCopySize_forces = numAtoms;
    int toCopyStart_slow_forces = 0;
    int toCopySize_slow_forces = slowForcesNumAtoms;
  #endif

  #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
    #if MIC_DEVICE_TRACING_DETAILED != 0
      int toCopySize__device_times_computes = ((isRemote) ? (0) : (2 * device__patch_pairs_size));
      int toCopySize__device_times_patches = ((isRemote) ? (0) : (2 * device__force_lists_size));
    #endif
    int toCopySize__device_times_start = ((isRemote) ? (0) : (10));
    #if MIC_DEVICE_TRACING_DETAILED != 0
      #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE \
        out(device__device_times_computes[0:toCopySize__device_times_computes] : alloc_if(0) free_if(0)) \
        out(device__device_times_patches[0:toCopySize__device_times_patches] : alloc_if(0) free_if(0)) \
        out(device__device_times_start[0:toCopySize__device_times_start] : alloc_if(0) free_if(0))
    #else
      #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE \
        out(device__device_times_start[0:toCopySize__device_times_start] : alloc_if(0) free_if(0))
    #endif
  #else
    #define MIC_DEVICE_TIMING_PRAGMA_CLAUSE
  #endif

  //// DMK - DEBUG
  //printf("[DEBUG:NBF-%d] ::   atoms[0:%d]\n", isRemote, toCopySize_atoms);
  //printf("[DEBUG:NBF-%d] ::   atom_params[0:%d]\n", isRemote, toCopySize_atom_params);
  //printf("[DEBUG:NBF-%d] ::   patch_pairs[0:%d]\n", isRemote, toCopySize_patch_pairs);
  //printf("[DEBUG:NBF-%d] ::   force_lists[0:%d]\n", isRemote, toCopySize_force_lists);
  //printf("[DEBUG:NBF-%d] ::   forces[%d:%d]\n", isRemote, toCopyStart_forces, toCopySize_forces);
  //printf("[DEBUG:NBF-%d] ::   slow_forces[%d:%d]\n", isRemote, toCopyStart_slow_forces, toCopySize_slow_forces);

  #if (VTUNE_INFO != 0) && (VTUNE_TASKS != 0)
    __itt_frame_begin_v3(vtune_domain, NULL);
    //__itt_task_begin(vtune_domain, __itt_null, __itt_null, vtune_offload_task);
  #endif

  #if MIC_TRACING != 0
    double pragma_start = CmiWallTimer();
  #endif

  #if MIC_PRAGMA_TIMING_STATS != 0
    #if MIC_TRACING != 0
      pragma_timing_start = pragma_start;
    #else
      pragma_timing_start = CmiWallTimer();
    #endif
  #endif

  #if MIC_KERNEL_DATA_TRANSFER_STATS != 0
    int transfer_in = sizeof(isRemote)
      + (3 * sizeof(device__lata))
      + sizeof(int) // device__force_buffers_req_size
      + (sizeof(atom) * toCopySize_atoms)
      + (sizeof(atom_param) * toCopySize_atom_params)
      + (sizeof(patch_pair) * toCopySize_patch_pairs)
      + (sizeof(force_list) * toCopySize_force_lists);
    int transfer_inout = sizeof(mic_kernel_data);
    int transfer_out = sizeof(double4) * (toCopySize_forces + toCopySize_slow_forces)
      + sizeof(int); // device__timestep (debug output, not actually required)
    printf("[MIC-DEBUG] :: PE %04d-%05d-%1d :: transfering %d bytes (in:%d, inout:%d, out:%d)\n",
           device__myPE, device__timestep, isRemote,
           transfer_in + transfer_inout + transfer_out, transfer_in, transfer_inout, transfer_out
          );
  #endif

  //// DMK - DEBUG
  //printf("[DEBUG] :: << issuing %s kernel >>\n", (isRemote) ? ("remote") : ("local")); fflush(NULL);

  // Trigger computation on the MIC device (asynchronous)
  //  in(device__atoms[0:numAtoms] : alloc_if(0) free_if(0) align(MIC_ALIGN))
  //  in(device__atom_params[0:_atom_params_copySize] : alloc_if(0) free_if(0) align(MIC_ALIGN))
  //  in(device__patch_pairs[0:_patch_pairs_copySize] : alloc_if(0) free_if(0) align(MIC_ALIGN))
  //  in(device__force_lists[0:_force_lists_copySize] : alloc_if(0) free_if(0) align(MIC_ALIGN))
  //  out(device__forces[0:numAtoms] : alloc_if(0) free_if(0) align(MIC_ALIGN))
  //  out(device__slow_forces[0:slowForcesNumAtoms] : alloc_if(0) free_if(0) align(MIC_ALIGN))
  #pragma offload target(mic:deviceNum) \
    in(isRemote) \
    in(device__lata) in(device__latb) in(device__latc) \
    in(device__force_buffers_req_size) \
    in(device__atoms[0:toCopySize_atoms] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
    in(device__atom_params[0:toCopySize_atom_params] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
    in(device__patch_pairs[0:toCopySize_patch_pairs] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
    in(device__force_lists[0:toCopySize_force_lists] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
    inout(device__local_kernel_data[0:localKernelDataLen] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
    inout(device__remote_kernel_data[0:remoteKernelDataLen] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
    out(device__forces[toCopyStart_forces:toCopySize_forces] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
    out(device__slow_forces[toCopyStart_slow_forces:toCopySize_slow_forces] : alloc_if(0) free_if(0) align(MIC_ALIGN)) \
    MIC_PRAGMA_TIMING_DEVICE_CLAUSE \
    MIC_DEVICE_TIMING_PRAGMA_CLAUSE \
    nocopy(device__force_buffers_alloc_size) \
    nocopy(device__force_buffers[0:0] : alloc_if(0) free_if(0)) \
    nocopy(device__slow_force_buffers[0:0] : alloc_if(0) free_if(0)) \
    nocopy(device__patch_pairs_size) \
    nocopy(device__force_lists_size) \
    nocopy(device__atoms_size) \
    nocopy(device__pairlists) nocopy(device__pairlists_alloc_size) \
    nocopy(device__table_four[0:0] : alloc_if(0) free_if(0)) \
    nocopy(device__table_four_n_16) \
    nocopy(device__lj_table_dim) \
    nocopy(device__lj_table[0:0] : alloc_if(0) free_if(0)) \
    nocopy(device__exclusion_bits[0:0] : alloc_if(0) free_if(0)) \
    nocopy(device__constants[0:0] : alloc_if(0) free_if(0)) \
    nocopy(device__pl_array) nocopy(device__pl_size) nocopy(device__r2_array) \
    out(device__timestep) \
    nocopy(device__loop_counts) nocopy(device__cutoff_clear_counts) nocopy(device__iTest_counts) \
    nocopy(device__table_four_float) nocopy(device__lj_table_float) \
    nocopy(device__virial_energy) \
    signal(tag_kernel)
  {

    // DMK - DEBUG
    #if (MIC_PRAGMA_TIMING_STATS != 0)
      double __device_start = getCurrentTime();
    #endif

    // DMK - TRACING
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(0):(1))] = getCurrentTime();
    #endif

    //// DMK - DEBUG
    //printf("[DEBUG:MIC] :: Starting %s kernel for timestep %d...\n",
    //       (isRemote) ? ("remote") : ("local"),
    //       device__timestep
    //      ); fflush(NULL);

    #if !defined(__MIC__)
      mic_dev_die("mic_nonbonded_forces :: Attempt to execute MIC kernels on host (not supported yet)... exiting."); fflush(NULL);
    #endif

    // DMK - DEBUG
    #if MIC_ACTIVE_CUTOFF_STATS != 0
      activeCount = 0; activePossibleCount = 0;
      cutoffCount = 0; cutoffPossibleCount = 0;
    #endif

    // DMK - DEBUG
    #if MIC_GATHER_SCATTER_STATS != 0
      for (int i = 0; i < 17; i++) {
        gather_single[i] = gather_double[i] = scatter_single[i] = scatter_double[i] = 0;
      }
    #endif

    // DMK - NOTE | TODO | FIXME : All of the kernel calls will be "local" for now, so just use
    //   that kernel data.  However, a mechanism for having two kernel data structs and knowing which
    //   one is in use needs to be created.  Perhaps creating a function instead of a code block and
    //   making it a parameter(?).
    #if MIC_MULTIPLE_KERNELS != 0
      mic_kernel_data *kernel_data = (isRemote) ? (device__remote_kernel_data) : (device__local_kernel_data);
    #else
      mic_kernel_data *kernel_data = device__local_kernel_data;
    #endif

    //// DMK - DEBUG
    //printf("[DEBUG:MIC-%d-%d] ::   kernel_data...\n", device__timestep, isRemote);
    //printf("[DEBUG:MIC-%d-%d] ::     isRemote:%d(%d)\n", device__timestep, isRemote, kernel_data->isRemote, isRemote);
    //printf("[DEBUG:MIC-%d-%d] ::     numLocalAtoms:%d\n", device__timestep, isRemote, kernel_data->numLocalAtoms);
    //printf("[DEBUG:MIC-%d-%d] ::     numLocalComputes:%d\n", device__timestep, isRemote, kernel_data->numLocalComputes);
    //printf("[DEBUG:MIC-%d-%d] ::     numLocalPatches:%d\n", device__timestep, isRemote, kernel_data->numLocalPatches);
    //fflush(NULL);

    // Make sure there is enough memory allocated for the force buffers
    #if MIC_FULL_CHECK != 0
      assert(device__force_buffers_req_size > 0);
    #endif
    if (device__force_buffers_req_size > device__force_buffers_alloc_size) {

      //// DMK - DEBUG
      //printf("[MIC] :: reallocating force buffers (%d->%d)\n",
      //       device__force_buffers_alloc_size, device__force_buffers_req_size
      //      );

      if (device__force_buffers != NULL) { _mm_free(device__force_buffers); }
      if (device__slow_force_buffers != NULL) { _mm_free(device__slow_force_buffers); }
      device__force_buffers_alloc_size = device__force_buffers_req_size;
      //device__force_buffers = (double4*)(_mm_malloc(sizeof(double) * 4 * device__force_buffers_alloc_size, 64));
      //device__slow_force_buffers = (double4*)(_mm_malloc(sizeof(double) * 4 * device__force_buffers_alloc_size, 64));
      device__force_buffers = (double4*)(_MM_MALLOC_WRAPPER(sizeof(double) * 4 * device__force_buffers_alloc_size, 64, "device__force_buffers"));
      device__slow_force_buffers = (double4*)(_MM_MALLOC_WRAPPER(sizeof(double) * 4 * device__force_buffers_alloc_size, 64, "device__slow_force_buffers"));
      #if MIC_FULL_CHECK != 0
        assert(device__force_buffers != NULL && device__slow_force_buffers != NULL);
      #endif
    }

    //// Make sure there are enough pairlists pointers (several per patch_pair)
    //const int numPairlists = NUM_PAIRLIST_TYPES * device__patch_pairs_size;
    //if (numPairlists > device__pairlists_alloc_size) {
    //  int **new_pairlists = new int*[numPairlists];
    //  int initStart = 0;
    //  if (device__pairlists != 0) {
    //    int **old_pairlists = (int**)(device__pairlists);
    //    memcpy((void*)new_pairlists, (void*)old_pairlists, sizeof(int*) * device__pairlists_alloc_size);
    //    delete [] old_pairlists;
    //    initStart = device__pairlists_alloc_size;
    //  }
    //  for (int i = initStart; i < numPairlists; i++) { new_pairlists[i] = NULL; }
    //  device__pairlists = (unsigned long long int)new_pairlists;
    //  device__pairlists_alloc_size = numPairlists;
    //}

    // Initialize the overall virial summation/accumulation variables that will
    //   eventually be passed back to the host
    #if MIC_VIRIAL_ENERGY_ALT == 0
      double virial_xx = 0.0;
      double virial_xy = 0.0;
      double virial_xz = 0.0;
      double virial_yy = 0.0;
      double virial_yz = 0.0;
      double virial_zz = 0.0;
      double fullElectVirial_xx = 0.0;
      double fullElectVirial_xy = 0.0;
      double fullElectVirial_xz = 0.0;
      double fullElectVirial_yy = 0.0;
      double fullElectVirial_yz = 0.0;
      double fullElectVirial_zz = 0.0;
      double vdwEnergy = 0.0;
      double electEnergy = 0.0;
      double fullElectEnergy = 0.0;
    #endif

    //// DMK - DEBUG
    //printf("[MIC-DEBUG] :: Starting computes on device (timestep: %d)...\n", device__timestep);
    //fflush(NULL);

    #if MIC_MULTIPLE_KERNELS != 0
      const int ppI_start = (isRemote) ? (kernel_data->numLocalComputes) : (0);
      const int ppI_end = (isRemote) ? (device__patch_pairs_size) : (kernel_data->numLocalComputes);
    #else
      const int ppI_start = 0;
      const int ppI_end = device__patch_pairs_size;
    #endif

    //// DMK - DEBUG
    //printf("[DEBUG:MIC-%d-%d] ::   processing computes[%d->%d] - device__patch_pairs_size:%d\n",
    //       device__timestep, isRemote,
    //       ppI_start, ppI_end, device__patch_pairs_size
    //      ); fflush(NULL);

    // DMK - TRACING
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(2):(3))] = getCurrentTime();
    #endif

    // Compute each compute (patch pair) in parallel with one another
    #pragma novector
    #if MULTIPLE_THREADS != 0
      #if MIC_VIRIAL_ENERGY_ALT != 0
        #pragma omp parallel for schedule(dynamic, 1)
      #else
        #pragma omp parallel for schedule(dynamic, 1) \
          reduction(+ : vdwEnergy, electEnergy, fullElectEnergy) \
          reduction(+ : virial_xx, virial_xy, virial_xz, virial_yy, virial_yz, virial_zz) \
          reduction(+ : fullElectVirial_xx, fullElectVirial_xy, fullElectVirial_xz, fullElectVirial_yy, fullElectVirial_yz, fullElectVirial_zz )
      #endif
    #endif
    for (int ppI = ppI_start; ppI < ppI_end; ppI++) {

      #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
        device__device_times_computes[ppI * 2] = getCurrentTime();
      #endif

      // DMK - DEBUG
      //printf("[MIC-DEBUG] ::   ppI: %d\n", ppI); fflush(NULL);

      // Create and populate the mic_params data structure to be passed into the kernel function
      mic_params params;

      // Initialize the pointer to the patch_pairs structure for this compute
      params.pp = (patch_pair*)(device__patch_pairs + ppI);

      #if MIC_SORT_LISTS != 0
        const int abSwapFlag = ((params.pp->patch1_size < params.pp->patch2_size) ? (1) : (0));
        #define ABSWAP(t, f)  ((abSwapFlag)?(t):(f))
      #else
        #define ABSWAP(t, f)  (f)
      #endif

      // DMK - DEBUG - A few values used for debugging
      params.ppI = ppI;
      params.p1 = ABSWAP(params.pp->p2, params.pp->p1);
      params.p2 = ABSWAP(params.pp->p1, params.pp->p2);

      // In some cases, buffers used by the compute kernel functions (calc_self, etc.) only
      //   require scratch buffers (i.e. not used for input or output).  As such, we can
      //   allocate these buffers on a one-per-thread basis rather than a one-per-compute
      //   basis to save some memory (and get memory location reuse for any given thread).
      //   These fields relate to the pointer and sizes for those buffers.
      const int threadIndex = omp_get_thread_num();
      // If pairlist refinement is enabled, initialize the related scratch buffer pointers
      #if REFINE_PAIRLISTS != 0
        params.plArray_ptr = device__pl_array + threadIndex;
        params.plSize_ptr = device__pl_size + threadIndex;
        params.r2Array_ptr = device__r2_array + threadIndex;
      #endif

      // Initialize pointers to the various lookup tables, constants, etc. that are used
      //   from within the compute functions
      //#if (MIC_HANDCODE_FORCE != 0) && (MIC_HANDCODE_FORCE_SINGLE != 0)
      #if MIC_HANDCODE_FORCE_SINGLE != 0
        params.table_four_base_ptr = (void*)device__table_four_float;
        params.lj_table_base_ptr = (void*)device__lj_table_float;
      #else
        params.table_four_base_ptr = (void*)device__table_four;
        params.lj_table_base_ptr = (void*)device__lj_table;
      #endif
      params.table_four_n_16 = device__table_four_n_16;
      params.lj_table_dim = device__lj_table_dim;
      params.exclusion_bits = device__exclusion_bits;
      params.constants = device__constants;

      // Setup pairlist flags and pairlist buffer (an array of pointers, one for each pairlist,
      //   for each compute)
      params.usePairlists = kernel_data->usePairlists; //device__usePairlists;
      params.savePairlists = kernel_data->savePairlists; //device__savePairlists;
      params.pairlists_ptr = &(((int**)device__pairlists)[NUM_PAIRLIST_TYPES * ppI]);

      // Setup the sizes of the atom lists and force lists
      int n0 = ABSWAP(params.pp->patch2_size, params.pp->patch1_size);
      int n1 = ABSWAP(params.pp->patch1_size, params.pp->patch2_size);
      int n0_16 = (n0 + 15) & (~15); // Round up to nearest multiple of 16
      int n1_16 = (n1 + 15) & (~15);
      params.numAtoms[0] = n0;
      params.numAtoms[1] = n1;
      params.numAtoms_16[0] = n0_16;
      params.numAtoms_16[1] = n1_16;

      // Setup the pointers to the input particle data for this compute
      // NOTE: All of the particle data is sent in two chunks (atoms and atom_params, where
      //   atoms changes every timestep but atom_params only changes periodically).  Each
      //   of these buffers contains several sub-arrays (one for each field) in an
      //   structure-of-arrays (SOA) format.  This code, along with the force code
      //   below it, is simply creating the array pointers for each "field's array."
      mic_position_t *pBase0 = (mic_position_t*)(device__atoms + ABSWAP(params.pp->patch2_atom_start, params.pp->patch1_atom_start));
      mic_position_t *pBase1 = (mic_position_t*)(device__atoms + ABSWAP(params.pp->patch1_atom_start, params.pp->patch2_atom_start));
      int *pExtBase0 = (int*)(device__atom_params + ABSWAP(params.pp->patch2_atom_start, params.pp->patch1_atom_start));
      int *pExtBase1 = (int*)(device__atom_params + ABSWAP(params.pp->patch1_atom_start, params.pp->patch2_atom_start));
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        params.p[0] = (atom*)pBase0;
        params.p[1] = (atom*)pBase1;
        params.pExt[0] = (atom_param*)pExtBase0;
        params.pExt[1] = (atom_param*)pExtBase1;
      #else
        params.p_x[0] = pBase0 + (0 * n0_16); params.p_x[1] = pBase1 + (0 * n1_16);
        params.p_y[0] = pBase0 + (1 * n0_16); params.p_y[1] = pBase1 + (1 * n1_16);
        params.p_z[0] = pBase0 + (2 * n0_16); params.p_z[1] = pBase1 + (2 * n1_16);
        params.p_q[0] = pBase0 + (3 * n0_16); params.p_q[1] = pBase1 + (3 * n1_16);
        params.pExt_vdwType[0] = pExtBase0 + (0 * n0_16);
        params.pExt_index[0] = pExtBase0 + (1 * n0_16);
        params.pExt_exclIndex[0] = pExtBase0 + (2 * n0_16);
        params.pExt_exclMaxDiff[0] = pExtBase0 + (3 * n0_16);
        params.pExt_vdwType[1] = pExtBase1 + (0 * n1_16);
        params.pExt_index[1] = pExtBase1 + (1 * n1_16);
        params.pExt_exclIndex[1] = pExtBase1 + (2 * n1_16);
        params.pExt_exclMaxDiff[1] = pExtBase1 + (3 * n1_16);
      #endif

      // Setup the pointers to the output force data for this compute
      // NOTE: Forces are output every timestep, but slow forces (fullf) are output only
      //   during some timesteps.
      params.doSlow = kernel_data->doSlow;  // Flag indicating if slow forces should be calculate this timestep or not
      params.fl[0] = (force_list*)(device__force_lists + ABSWAP(params.pp->patch2_force_list_index, params.pp->patch1_force_list_index));
      params.fl[1] = (force_list*)(device__force_lists + ABSWAP(params.pp->patch1_force_list_index, params.pp->patch2_force_list_index));
      double *ffBase0 = (double*)(device__force_buffers + ABSWAP(params.pp->patch2_force_start, params.pp->patch1_force_start));
      double *ffBase1 = (double*)(device__force_buffers + ABSWAP(params.pp->patch1_force_start, params.pp->patch2_force_start));
      double *fullfBase0 = (double*)(device__slow_force_buffers + ABSWAP(params.pp->patch2_force_start, params.pp->patch1_force_start));
      double *fullfBase1 = (double*)(device__slow_force_buffers + ABSWAP(params.pp->patch1_force_start, params.pp->patch2_force_start));
      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
        params.ff[0] = (double4*)ffBase0;
        params.ff[1] = (double4*)ffBase1;
        params.fullf[0] = (double4*)fullfBase0;
        params.fullf[1] = (double4*)fullfBase1;
      #else
        params.ff_x[0] = ffBase0 + (0 * n0_16); params.ff_x[1] = ffBase1 + (0 * n1_16);
        params.ff_y[0] = ffBase0 + (1 * n0_16); params.ff_y[1] = ffBase1 + (1 * n1_16);
        params.ff_z[0] = ffBase0 + (2 * n0_16); params.ff_z[1] = ffBase1 + (2 * n1_16);
        params.ff_w[0] = ffBase0 + (3 * n0_16); params.ff_w[1] = ffBase1 + (3 * n1_16);
        params.fullf_x[0] = fullfBase0 + (0 * n0_16); params.fullf_x[1] = fullfBase1 + (0 * n1_16);
        params.fullf_y[0] = fullfBase0 + (1 * n0_16); params.fullf_y[1] = fullfBase1 + (1 * n1_16);
        params.fullf_z[0] = fullfBase0 + (2 * n0_16); params.fullf_z[1] = fullfBase1 + (2 * n1_16);
        params.fullf_w[0] = fullfBase0 + (3 * n0_16); params.fullf_w[1] = fullfBase1 + (3 * n1_16);
      #endif

      // Create the offsets for the first list of particles.
      // NOTE: In this version of the nonbonded code, the positions of the atoms are stored
      //   as offsets within the patch in which they are located.  These offsets represent
      //   the offsets between the two patches being interacted (including periodic boundaries).
      params.offset.x = params.pp->offset.x * device__lata.x
                      + params.pp->offset.y * device__latb.x
                      + params.pp->offset.z * device__latc.x;
      params.offset.y = params.pp->offset.x * device__lata.y
                      + params.pp->offset.y * device__latb.y
                      + params.pp->offset.z * device__latc.y;
      params.offset.z = params.pp->offset.x * device__lata.z
                      + params.pp->offset.y * device__latb.z
                      + params.pp->offset.z * device__latc.z;
      #if MIC_SORT_LISTS != 0
        params.offset.x *= ABSWAP(-1.0f, 1.0f);
        params.offset.y *= ABSWAP(-1.0f, 1.0f);
        params.offset.z *= ABSWAP(-1.0f, 1.0f);
      #endif

      // DMK - DEBUG - Record the center point for each of the input patches
      float patch1_center_x = params.pp->patch1_center_x * device__lata.x
                             + params.pp->patch1_center_y * device__latb.x
                             + params.pp->patch1_center_z * device__latc.x;
      float patch1_center_y = params.pp->patch1_center_x * device__lata.y
                             + params.pp->patch1_center_y * device__latb.y
                             + params.pp->patch1_center_z * device__latc.y;
      float patch1_center_z = params.pp->patch1_center_x * device__lata.z
                             + params.pp->patch1_center_y * device__latb.z
                             + params.pp->patch1_center_z * device__latc.z;
      float patch2_center_x = params.pp->patch2_center_x * device__lata.x
                             + params.pp->patch2_center_y * device__latb.x
                             + params.pp->patch2_center_z * device__latc.x;
      float patch2_center_y = params.pp->patch2_center_x * device__lata.y
                             + params.pp->patch2_center_y * device__latb.y
                             + params.pp->patch2_center_z * device__latc.y;
      float patch2_center_z = params.pp->patch2_center_x * device__lata.z
                             + params.pp->patch2_center_y * device__latb.z
                             + params.pp->patch2_center_z * device__latc.z;
      params.patch1_center_x = ABSWAP(patch2_center_x, patch1_center_x);
      params.patch2_center_x = ABSWAP(patch1_center_x, patch2_center_x);
      params.patch1_center_y = ABSWAP(patch2_center_y, patch1_center_y);
      params.patch2_center_y = ABSWAP(patch1_center_y, patch2_center_y);
      params.patch1_center_z = ABSWAP(patch2_center_z, patch1_center_z);
      params.patch2_center_z = ABSWAP(patch1_center_z, patch2_center_z);

      // Initialize the virial accumulators for this compute
      params.virial_xx = 0;
      params.virial_xy = 0;
      params.virial_xz = 0;
      params.virial_yy = 0;
      params.virial_yz = 0;
      params.virial_zz = 0;
      params.fullElectVirial_xx = 0;
      params.fullElectVirial_xy = 0;
      params.fullElectVirial_xz = 0;
      params.fullElectVirial_yy = 0;
      params.fullElectVirial_yz = 0;
      params.fullElectVirial_zz = 0;
      params.vdwEnergy = 0.0;
      params.electEnergy = 0.0;
      params.fullElectEnergy = 0.0;

      // Select the version of the kernel to call based on the timestep's requirements
      //   and what type of compute this compute is
      int isSelf = (params.pp->patch1_force_list_index == params.pp->patch2_force_list_index);  // NOTE: Many ways to check this (arbitrary test used here)
      int selfBit = ((isSelf) ? (0x01) : (0x00));
      int doSlowBit = ((kernel_data->doSlow) ? (0x02) : (0x00));
      int doEnergyBit = ((kernel_data->doEnergy) ? (0x04) : (0x00));
      int kernelSelect = selfBit | doSlowBit | doEnergyBit;
      switch (kernelSelect) {
        case 0x00: calc_pair(params); break;
        case 0x01: calc_self(params); break;
        case 0x02: calc_pair_fullelect(params); break;
        case 0x03: calc_self_fullelect(params); break;
        case 0x04: calc_pair_energy(params); break;
        case 0x05: calc_self_energy(params); break;
        case 0x06: calc_pair_energy_fullelect(params); break;
        case 0x07: calc_self_energy_fullelect(params); break;
        default:
          mic_dev_die("!!! INVALID KERNEL SELECTION ON MIC DEVICE !!!\n");
          break;
      } // end switch (kernelSelect)

      // Contribute this compute's virial summations into the overall virial summation
      //   that will be passed back to the host
      #if MIC_SORT_LISTS != 0
        if (abSwapFlag) {
          #if MIC_VIRIAL_ENERGY_ALT != 0
            int veI = device__numOMPThreads + ppI;
            device__virial_energy[veI +  0] = -1 * params.virial_xx;
            device__virial_energy[veI +  1] = -1 * params.virial_xy;
            device__virial_energy[veI +  2] = -1 * params.virial_xz;
            device__virial_energy[veI +  3] = -1 * params.virial_yy;
            device__virial_energy[veI +  4] = -1 * params.virial_yz;
            device__virial_energy[veI +  5] = -1 * params.virial_zz;
            device__virial_energy[veI +  6] = -1 * params.fullElectVirial_xx;
            device__virial_energy[veI +  7] = -1 * params.fullElectVirial_xy;
            device__virial_energy[veI +  8] = -1 * params.fullElectVirial_xz;
            device__virial_energy[veI +  9] = -1 * params.fullElectVirial_yy;
            device__virial_energy[veI + 10] = -1 * params.fullElectVirial_yz;
            device__virial_energy[veI + 11] = -1 * params.fullElectVirial_zz;
            device__virial_energy[veI + 12] = -1 * params.vdwEnergy;
            device__virial_energy[veI + 13] = -1 * params.electEnergy;
            device__virial_energy[veI + 14] = -1 * params.fullElectEnergy;
          #else
            virial_xx -= params.virial_xx;
            virial_xy -= params.virial_xy;
            virial_xz -= params.virial_xz;
            virial_yy -= params.virial_yy;
            virial_yz -= params.virial_yz;
            virial_zz -= params.virial_zz;
            fullElectVirial_xx -= params.fullElectVirial_xx;
            fullElectVirial_xy -= params.fullElectVirial_xy;
            fullElectVirial_xz -= params.fullElectVirial_xz;
            fullElectVirial_yy -= params.fullElectVirial_yy;
            fullElectVirial_yz -= params.fullElectVirial_yz;
            fullElectVirial_zz -= params.fullElectVirial_zz;
            vdwEnergy -= params.vdwEnergy;
            electEnergy -= params.electEnergy;
            fullElectEnergy -= params.fullElectEnergy;
          #endif
        } else {
      #endif
          #if MIC_VIRIAL_ENERGY_ALT != 0
            int veI = device__numOMPThreads + ppI;
            device__virial_energy[veI +  0] = params.virial_xx;
            device__virial_energy[veI +  1] = params.virial_xy;
            device__virial_energy[veI +  2] = params.virial_xz;
            device__virial_energy[veI +  3] = params.virial_yy;
            device__virial_energy[veI +  4] = params.virial_yz;
            device__virial_energy[veI +  5] = params.virial_zz;
            device__virial_energy[veI +  6] = params.fullElectVirial_xx;
            device__virial_energy[veI +  7] = params.fullElectVirial_xy;
            device__virial_energy[veI +  8] = params.fullElectVirial_xz;
            device__virial_energy[veI +  9] = params.fullElectVirial_yy;
            device__virial_energy[veI + 10] = params.fullElectVirial_yz;
            device__virial_energy[veI + 11] = params.fullElectVirial_zz;
            device__virial_energy[veI + 12] = params.vdwEnergy;
            device__virial_energy[veI + 13] = params.electEnergy;
            device__virial_energy[veI + 14] = params.fullElectEnergy;
          #else
            virial_xx += params.virial_xx;
            virial_xy += params.virial_xy;
            virial_xz += params.virial_xz;
            virial_yy += params.virial_yy;
            virial_yz += params.virial_yz;
            virial_zz += params.virial_zz;
            fullElectVirial_xx += params.fullElectVirial_xx;
            fullElectVirial_xy += params.fullElectVirial_xy;
            fullElectVirial_xz += params.fullElectVirial_xz;
            fullElectVirial_yy += params.fullElectVirial_yy;
            fullElectVirial_yz += params.fullElectVirial_yz;
            fullElectVirial_zz += params.fullElectVirial_zz;
            vdwEnergy += params.vdwEnergy;
            electEnergy += params.electEnergy;
            fullElectEnergy += params.fullElectEnergy;
          #endif
      #if MIC_SORT_LISTS != 0
	}
      #endif

      #undef ABSWAP

      // DMK - TRACING
      #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
        device__device_times_computes[ppI * 2 + 1] = getCurrentTime();

        //// DMK - DEBUG
        //if (device__timestep == 3) {
        //  printf("[DEV-TIME] pp[%04d] %.6lf -> %.6lf (%.6lf)\n",
        //         ppI,
        //         device__device_times_computes[ppI * 2], device__device_times_computes[ppI * 2 + 1],
        //         device__device_times_computes[ppI * 2 + 1] - device__device_times_computes[ppI * 2]
	//        );
	//}

      #endif

    } // end parallel for (ppI < device__patch_pairs_size)

    // DMK - TRACING
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(4):(5))] = getCurrentTime();

      //// DMK - DEBUG
      //if (device__timestep == 3) {
      //  printf("[DEV-TIME] computes %0.6lf -> %.6lf (%.6lf)\n",
      //         device__device_times_start[((isRemote)?(2):(3))], device__device_times_start[((isRemote)?(4):(5))],
      //         device__device_times_start[((isRemote)?(4):(5))] - device__device_times_start[((isRemote)?(2):(3))]
      //        );
      //}

    #endif

    // DMK - DEBUG
    #if MIC_ACTIVE_CUTOFF_STATS != 0
      printf("Mask Stats: active: %d / %d , cutoff: %d / %d\n", activeCount, activePossibleCount, cutoffCount, cutoffPossibleCount);
    #endif

    // DMK - DEBUG
    #if MIC_GATHER_SCATTER_STATS != 0
      printf("Gather Double Stats: (0) %d _ %d %d %d %d _ %d %d %d %d (8)\n",
             gather_double[0], gather_double[1], gather_double[2], gather_double[3],
             gather_double[4], gather_double[5], gather_double[6], gather_double[7],
             gather_double[8]
            );
      printf("Gather Single Stats: (0) %d _ %d %d %d %d _ %d %d %d %d _ %d %d %d %d _ %d %d %d %d (16)\n",
             gather_single[ 0], gather_single[ 1], gather_single[ 2], gather_single[ 3],
             gather_single[ 4], gather_single[ 5], gather_single[ 6], gather_single[ 7],
             gather_single[ 8], gather_single[ 9], gather_single[10], gather_single[11],
             gather_single[12], gather_single[13], gather_single[14], gather_single[15],
             gather_single[16]
            );
      printf("Scatter Double Stats: (0) %d _ %d %d %d %d _ %d %d %d %d (8)\n",
             scatter_double[0], scatter_double[1], scatter_double[2], scatter_double[3],
             scatter_double[4], scatter_double[5], scatter_double[6], scatter_double[7],
             scatter_double[8]
            );
      printf("Scatter Single Stats: (0) %d _ %d %d %d %d _ %d %d %d %d _ %d %d %d %d _ %d %d %d %d (16)\n",
             scatter_single[ 0], scatter_single[ 1], scatter_single[ 2], scatter_single[ 3],
             scatter_single[ 4], scatter_single[ 5], scatter_single[ 6], scatter_single[ 7],
             scatter_single[ 8], scatter_single[ 9], scatter_single[10], scatter_single[11],
             scatter_single[12], scatter_single[13], scatter_single[14], scatter_single[15],
             scatter_single[16]
            );
    #endif

    // Store the reduced virial values into the kernel_data structure to be passed
    //   back to the host core
    #if MIC_VIRIAL_ENERGY_ALT != 0
      if (isRemote == 0) {
        __ASSUME_ALIGNED(device__virial_energy);
        for (int ppI = 0; ppI < device__patch_pairs_size; ppI++) {
          int veI = device__numOMPThreads + ppI;
          kernel_data->virial_xx += device__virial_energy[veI + 0];
          kernel_data->virial_xy += device__virial_energy[veI + 1];
          kernel_data->virial_xz += device__virial_energy[veI + 2];
          kernel_data->virial_yy += device__virial_energy[veI + 3];
          kernel_data->virial_yz += device__virial_energy[veI + 4];
          kernel_data->virial_zz += device__virial_energy[veI + 5];
          kernel_data->fullElectVirial_xx += device__virial_energy[veI + 6];
          kernel_data->fullElectVirial_xy += device__virial_energy[veI + 7];
          kernel_data->fullElectVirial_xz += device__virial_energy[veI + 8];
          kernel_data->fullElectVirial_yy += device__virial_energy[veI + 9];
          kernel_data->fullElectVirial_yz += device__virial_energy[veI + 10];
          kernel_data->fullElectVirial_zz += device__virial_energy[veI + 11];
          kernel_data->vdwEnergy = device__virial_energy[veI + 12];
          kernel_data->electEnergy = device__virial_energy[veI + 13];
          kernel_data->fullElectEnergy = device__virial_energy[veI + 14];
        }
      }
    #else
      kernel_data->virial_xx = virial_xx;
      kernel_data->virial_xy = virial_xy;
      kernel_data->virial_xz = virial_xz;
      kernel_data->virial_yy = virial_yy;
      kernel_data->virial_yz = virial_yz;
      kernel_data->virial_zz = virial_zz;
      kernel_data->fullElectVirial_xx = fullElectVirial_xx;
      kernel_data->fullElectVirial_xy = fullElectVirial_xy;
      kernel_data->fullElectVirial_xz = fullElectVirial_xz;
      kernel_data->fullElectVirial_yy = fullElectVirial_yy;
      kernel_data->fullElectVirial_yz = fullElectVirial_yz;
      kernel_data->fullElectVirial_zz = fullElectVirial_zz;
      kernel_data->vdwEnergy = vdwEnergy;
      kernel_data->electEnergy = electEnergy;
      kernel_data->fullElectEnergy = fullElectEnergy;
    #endif

    //// DMK - DEBUG
    //printf("[MIC-DEBUG] :: Starting force reductions (timestep: %d)...\n", device__timestep);
    //fflush(NULL);

    #if MIC_MULTIPLE_KERNELS != 0
      int flI_start = (isRemote) ? (kernel_data->numLocalPatches) : (0);
      int flI_end = (isRemote) ? (device__force_lists_size) : (kernel_data->numLocalPatches);
    #else
      int flI_start = 0;
      int flI_end = device__force_lists_size;
    #endif

    //// DMK - DEBUG
    //printf("[DEBUG:MIC-%d-%d] ::   processing patches[%d->%d] - device__force_lists_size:%d\n",
    //       device__timestep, isRemote,
    //       flI_start, flI_end, device__force_lists_size
    //      ); fflush(NULL);

    // Once all of the computes (patch pairs) are complete, reduce the forces for
    //   each patch (force list)
    const int numThreads = device__numOMPThreads; //omp_get_max_threads();
    const int numForceLoopIters = flI_end - flI_start;
    int numForceLoopSplits = 1;
    if ((numForceLoopIters > 0) && ((2 * numThreads) > numForceLoopIters)) {
      numForceLoopSplits = (2 * numThreads) / numForceLoopIters;  // NOTE: 2 * numThreads to break it up more (smaller chunks of work)
      if (numForceLoopSplits < 2) { numForceLoopSplits = 2; }  // At least split in half
      if (numForceLoopSplits > 16) { numForceLoopSplits = 16; }  // Don't split any single patch too fine (threading overhead costs)
    }
    flI_start *= numForceLoopSplits;
    flI_end *= numForceLoopSplits;

    // DMK - TRACING
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(6):(7))] = getCurrentTime();
    #endif

    #pragma novector
    #if MULTIPLE_THREADS != 0
      #pragma omp parallel for schedule(dynamic)
    #endif
    for (int _flI = flI_start; _flI < flI_end; _flI++) {

      // Calculate this "task's" portion of the patch object's work (flIPartILo -> flIPartIHi)
      int flI = _flI / numForceLoopSplits;
      int flIPart = _flI % numForceLoopSplits;
      const force_list &fl = device__force_lists[flI];
      const int f_len = fl.patch_stride * 4; // NOTE : number of individual doubles
      __ASSUME(f_len % 16 == 0);
      int flIPartJLo = (int)(((float)(f_len)) * (((float)(flIPart    )) / ((float)(numForceLoopSplits))));
      int flIPartJHi = (int)(((float)(f_len)) * (((float)(flIPart + 1)) / ((float)(numForceLoopSplits))));
      // NOTE: Force flIPartJLo and flIPartJHi to cacheline boundaries to avoid false sharing
      flIPartJLo = (flIPartJLo + 7) & (~7);
      flIPartJHi = (flIPartJHi + 7) & (~7);
      if (flIPartJHi > f_len) { flIPartJHi = f_len; }
      __ASSUME(flIPartJLo % 8 == 0);
      __ASSUME(flIPartJHi % 8 == 0);  // NOTE: true after capping at f_len since f_len % 16 == 0

      // DMK - TRACING - NOTE: ONLY RECORDING ONE PART OF EACH SPLIT SET, SO PROJECTIONS WILL ONLY
      //   SOME OF THE FORCE "TASKS" !!!!!!!!!!!! (for now, find a solution - problem is numForceLoopSplits
      //   differs when multiple kernels are being used).
      #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
        if (flIPart == flI % numForceLoopSplits) {  // NOTE: Don't select the same part each time (e.g. only flIPart == 0)
          device__device_times_patches[flI * 2] = getCurrentTime();
	}
      #endif

      // Sum the output for each contributing compute
      {
        // Setup the pointer to the final force array that will be passed back up to the host
        double * RESTRICT fDst = (double*)(device__forces + fl.force_output_start);
        __ASSUME_ALIGNED(fDst);

        // Setup the pointer to the list of arrays, each with output from one of
        //   the compute objects that contributed to this patch's force data
        double *fSrcBase = (double*)(device__force_buffers + fl.force_list_start);
        __ASSUME_ALIGNED(fSrcBase);

        // Accumulate the forces from the various computes contributing to this patch
        memset(fDst + flIPartJLo, 0, sizeof(double) * (flIPartJHi - flIPartJLo));
        for (int i = 0; i < fl.force_list_size; i++) {
          #pragma simd
          for (int j = flIPartJLo; j < flIPartJHi; j++) {
            fDst[j] += fSrcBase[i * f_len + j];
          }
        }
      }

      // Sum the output for each contributing compute
      if (kernel_data->doSlow) {

        // Setup the pointer to the final force array that will be passed back up to the host
        double * RESTRICT fDst = (double*)(device__slow_forces + fl.force_output_start);
        __ASSUME_ALIGNED(fDst);

        // Setup the pointer to the list of arrays, each with output from one of
        //   the compute objects that contributed to this patch's force data
        double *fSrcBase = (double*)(device__slow_force_buffers + fl.force_list_start);
        __ASSUME_ALIGNED(fSrcBase);

        // Accumulate the forces from the various computes contributing to this patch
        memset(fDst + flIPartJLo, 0, sizeof(double) * (flIPartJHi - flIPartJLo));
        for (int i = 0; i < fl.force_list_size; i++) {
          #pragma simd
          for (int j = flIPartJLo; j < flIPartJHi; j++) {
            fDst[j] += fSrcBase[i * f_len + j];
          }
        }
      }
     
      // DMK - TRACING
      #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0) && (MIC_DEVICE_TRACING_DETAILED != 0)
        device__device_times_patches[flI * 2 + 1] = getCurrentTime();
      #endif

    } // end parallel for (flI < device__force_lists_size)

    // DMK - TRACING
    #if (MIC_TRACING != 0 && MIC_DEVICE_TRACING != 0)
      device__device_times_start[((isRemote)?(8):(9))] = getCurrentTime();
    #endif

    // DMK - DEBUG
    #if 0
      for (int flI = 0; flI < device__force_lists_size; flI++) {
        printf(":: Force_List %d...\n", flI);
        const force_list *fl = device__force_lists + flI;
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          const double4 *f = device__forces + fl->force_output_start;
          for (int i = 0; i < fl->patch_stride /*&& i < 5*/; i++) {
            printf("::   force[%3d,%3d] = { x:%+.3le, y:%+.3le, z:%+.3le, w:%+.3le }\n", flI, i, f[i].x, f[i].y, f[i].z, f[i].w);
	  }
        #else
          const double *f_x = (double*)(device__forces + fl->force_output_start);
          const double *f_y = f_x + fl->patch_stride;
          const double *f_z = f_y + fl->patch_stride;
          const double *f_w = f_z + fl->patch_stride;
          for (int i = 0; i < fl->patch_stride /*&& i < 5*/; i++) {
            printf("::   force[%3d,%3d] = { x:%+.3le, y:%+.3le, z:%+.3le, w:%+.3le }\n", flI, i, f_x[i], f_y[i], f_z[i], f_w[i]);
	  }
        #endif
      }

    #endif

    //// DMK - DEBUG
    //printf("[DEBUG:MIC] :: Finishing %s kernel for timestep %d...\n",
    //       (isRemote) ? ("remote") : ("local"),
    //       device__timestep
    //      ); fflush(NULL);

    // DMK - DEBUG
    if (kernel_data->isRemote == 0) { device__timestep++; }

    #if MIC_PRAGMA_TIMING_STATS != 0
      pragma_deviceExec_time[0] = getCurrentTime() - __device_start;

      //// DMK - DEBUG
      //printf("pragma_deviceExec_time = %.3lf ms\n", pragma_deviceExec_time * 1.0e3);

    #endif

  } // end pragma offload

  #if MIC_TRACING != 0
    double pragma_end = CmiWallTimer();
    MIC_TRACING_RECORD(MIC_EVENT_OFFLOAD_PRAGMA, pragma_start, pragma_end);
  #endif

  #if MIC_PRAGMA_TIMING_STATS != 0
    #if MIC_TRACING != 0
      pragma_timing_issued = pragma_end;
    #else
      pragma_timing_issued = CmiWallTimer();
    #endif
  #endif

  #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
    #undef MIC_DEVICE_TIMING_PRAGMA_CLAUSE
  #endif

  //// DMK - DEBUG
  //printf("mic_nonbonded :: &pragma_deviceExec_time = %p\n", &pragma_deviceExec_time);

}

int mic_check_remote_kernel_complete(const int deviceNum) {
  // DMK - NOTE : Disable these warnings for now since the device code is only called once,
  //   but the state variable goes from 0 -> 1 -> 2 -> 0, which checks at both 1 and 2
  //if (!tag_remote_kernel) { printf("WARNING :: mic_check_remote_kernel_complete :: called when kernel not active.\n"); }
  if (_Offload_signaled(deviceNum, &tag_remote_kernel)) {
    tag_remote_kernel = 0;
    return 1;
  }
  return 0;
}

int mic_check_local_kernel_complete(const int deviceNum) {
  // DMK - NOTE : Disable these warnings for now since the device code is only called once,
  //   but the state variable goes from 0 -> 1 -> 2 -> 0, which checks at both 1 and 2
  //if (!tag_local_kernel) { printf("WARNING :: mic_check_local_kernel_complete :: called when kernel not active.\n"); }
  if (_Offload_signaled(deviceNum, &tag_local_kernel)) {
    tag_local_kernel = 0;
    return 1;
  }
  return 0;
}

void _mic_print_stats(const int deviceNum) {

  // DMK - DEBUG
  #if MIC_STATS_LOOP_COUNTS != 0
    #pragma offload target(mic:deviceNum) \
      nocopy(device__loop_counts) nocopy(device__cutoff_clear_counts) nocopy(device__iTest_counts)
    {
      int totals[6] = { 0 };

      // Loop Counts
      memset(totals, 0, 6 * sizeof(int));
      printf("[MIC-LOOP-COUNTS] :: Start\n");
      printf("[MIC-LOOP-COUNTS], bin-min, self-normal, self-modified, self-excluded, pair-normal, pair-modified, pair-excluded\n");
      for (int i = 0; i < MIC_STATS_LOOP_COUNTS_NUM_BINS; i++) {
        int baseI = 6 * i;
        printf("[MIC-LOOP-COUNTS], %7d, %11d, %13d, %13d, %11d, %13d, %13d\n",
               i * MIC_STATS_LOOP_COUNTS_BIN_SIZE,
               device__loop_counts[baseI  ], device__loop_counts[baseI+1], device__loop_counts[baseI+2],
               device__loop_counts[baseI+3], device__loop_counts[baseI+4], device__loop_counts[baseI+5]
              );
        for (int j = 0; j < 6; j++) { totals[j] += device__loop_counts[baseI+j]; }
      }
      printf("[MIC-LOOP-COUNTS],  totals, %11d, %13d, %13d, %11d, %13d, %13d\n",
             totals[0], totals[1], totals[2], totals[3], totals[4], totals[5]
            );
      printf("[MIC-LOOP-COUNTS] :: End\n");

      // Cutoff Clear Counts
      memset(totals, 0, 6 * sizeof(int));
      printf("[MIC-CUTOFF-CLEAR-COUNTS] :: Start\n");
      printf("[MIC-CUTOFF-CLEAR-COUNTS], bin-min, self-normal, self-modified, self-excluded, pair-normal, pair-modified, pair-excluded\n");
      for (int i = 0; i < MIC_STATS_LOOP_COUNTS_NUM_BINS; i++) {
        int baseI = 6 * i;
        printf("[MIC-CUTOFF-CLEAR-COUNTS], %7d, %11d, %13d, %13d, %11d, %13d, %13d\n",
               i * MIC_STATS_LOOP_COUNTS_BIN_SIZE,
               device__cutoff_clear_counts[baseI  ], device__cutoff_clear_counts[baseI+1], device__cutoff_clear_counts[baseI+2],
               device__cutoff_clear_counts[baseI+3], device__cutoff_clear_counts[baseI+4], device__cutoff_clear_counts[baseI+5]
              );
        for (int j = 0; j < 6; j++) { totals[j] += device__cutoff_clear_counts[baseI+j]; }
      }
      printf("[MIC-CUTOFF-CLEAR-COUNTS],  totals, %11d, %13d, %13d, %11d, %13d, %13d\n",
             totals[0], totals[1], totals[2], totals[3], totals[4], totals[5]
            );
      printf("[MIC-CUTOFF-CLEAR-COUNTS] :: End\n");

      // iTest counts
      memset(totals, 0, 6 * sizeof(int));
      printf("[MIC-ITEST-COUNTS] :: Start\n");
      printf("[MIC-ITEST-COUNTS], bin-min, self-normal, self-modified, self-excluded, pair-normal, pair-modified, pair-excluded\n");
      for (int i = 0; i < MIC_STATS_LOOP_COUNTS_NUM_BINS; i++) {
        int baseI = 6 * i;
        printf("[MIC-ITEST-COUNTS], %7d, %11d, %13d, %13d, %11d, %13d, %13d\n",
               i * MIC_STATS_LOOP_COUNTS_BIN_SIZE,
               device__iTest_counts[baseI  ], device__iTest_counts[baseI+1], device__iTest_counts[baseI+2],
               device__iTest_counts[baseI+3], device__iTest_counts[baseI+4], device__iTest_counts[baseI+5]
              );
        for (int j = 0; j < 6; j++) { totals[j] += device__iTest_counts[baseI+j]; }
      }
      printf("[MIC-ITEST-COUNTS],  totals, %11d, %13d, %13d, %11d, %13d, %13d\n",
             totals[0], totals[1], totals[2], totals[3], totals[4], totals[5]
            );
      printf("[MIC-ITEST-COUNTS] :: End\n");
    }
  #endif
}

void mic_free_device(const int deviceNum) {

  // Cleanup kernel data (for buffers allocated via offload pragmas, use "free_if(1)")
  #pragma offload target(mic:deviceNum) \
    in(device__remote_kernel_data[0:0] : alloc_if(0) free_if(1)) \
    in(device__local_kernel_data[0:0] : alloc_if(0) free_if(1)) \
    in(device__table_four[0:0] : alloc_if(0) free_if(1)) \
    in(device__lj_table[0:0] : alloc_if(0) free_if(1)) \
    in(device__exclusion_bits[0:0] : alloc_if(0) free_if(1)) \
    in(device__constants[0:0] : alloc_if(0) free_if(1)) \
    in(device__patch_pairs[0:0] : alloc_if(0) free_if(1)) \
    in(device__force_lists[0:0] : alloc_if(0) free_if(1)) \
    in(device__atoms[0:0] : alloc_if(0) free_if(1)) \
    in(device__atom_params[0:0] : alloc_if(0) free_if(1)) \
    in(device__forces[0:0] : alloc_if(0) free_if(1)) \
    in(device__slow_forces[0:0] : alloc_if(0) free_if(1)) \
    in(device__force_buffers[0:0] : alloc_if(0) free_if(1)) \
    in(device__slow_force_buffers[0:0] : alloc_if(0) free_if(1)) \
    nocopy(device__pairlists) nocopy(device__pairlists_alloc_size) \
    nocopy(device__pl_array) nocopy(device__pl_size) nocopy(device__r2_array) \
    nocopy(device__loop_counts) nocopy(device__cutoff_clear_counts) nocopy(device__iTest_counts) \
    nocopy(device__table_four_float) nocopy(device__lj_table_float) \
    nocopy(device__virial_energy)
  {
    const int numThreads = omp_get_max_threads();

    // Cleanup pairlist memory
    if (device__pairlists != NULL) {
      int **pairlists = (int**)(device__pairlists);
      for (int i = 0; i < device__pairlists_alloc_size; i++) {
        if (pairlists[i] != NULL) { delete [] pairlists[i]; }
      }
      delete [] pairlists;
      device__pairlists = 0;
      device__pairlists_alloc_size = 0;
    }

    // Cleanup refined pairlist memory
    #if REFINE_PAIRLISTS != 0
      for (int i = 0; i < numThreads; i++) {
        if (device__pl_array[i] != NULL) { delete [] device__pl_array[i]; }
        if (device__r2_array[i] != NULL) { delete [] device__r2_array[i]; }
      }
      delete [] device__pl_array; device__pl_array = NULL;
      delete [] device__r2_array; device__r2_array = NULL;
      delete [] device__pl_size; device__pl_size = NULL;
    #endif

    //#if (MIC_HANDCODE_FORCE != 0) && (MIC_HANDCODE_FORCE_SINGLE != 0)
    #if MIC_HANDCODE_FORCE_SINGLE != 0
      if (device__table_four_float != NULL) { _mm_free(device__table_four_float); device__table_four_float = NULL; }
      if (device__lj_table_float != NULL) { _mm_free(device__lj_table_float); device__lj_table_float = NULL; }
    #endif

    // DMK - DEBUG
    #if MIC_STATS_LOOP_COUNTS != 0
      if (device__loop_counts != NULL) { delete [] device__loop_counts; device__loop_counts = NULL; }
      if (device__cutoff_clear_counts != NULL) { delete [] device__cutoff_clear_counts; device__cutoff_clear_counts = NULL; }
      if (device__iTest_counts != NULL) { delete [] device__iTest_counts; device__iTest_counts = NULL; }
    #endif

    #if MIC_VIRIAL_ENERGY_ALT != 0
      if (device__virial_energy != NULL) { _mm_free(device__virial_energy); device__virial_energy = NULL; }
    #endif
  }

  if (device__remote_kernel_data != NULL) { delete [] device__remote_kernel_data; device__remote_kernel_data = NULL; }
  if (device__local_kernel_data != NULL) { delete [] device__local_kernel_data; device__local_kernel_data = NULL; }
}

#else  // NAMD_MIC

// for make depends
#include "ComputeNonbondedMICKernel.h"
#include "ComputeNonbondedMICKernelBase.h"

#endif  // NAMD_MIC
