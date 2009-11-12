
#include "common.h"
#include "charm++.h"

#ifdef NAMD_CUDA
#include <cuda_runtime.h>
#include <cuda.h>
#endif

#include "WorkDistrib.h"
#include "ComputeMgr.h"
#include "ComputeNonbondedCUDA.h"
#include "ComputeNonbondedCUDAKernel.h"
#include "ObjectArena.h"

#ifdef NAMD_CUDA

extern cudaStream_t stream;

void cuda_errcheck(const char *msg) {
  cudaError_t err;
  if ((err = cudaGetLastError()) != cudaSuccess) {
    char errmsg[1024];
    sprintf(errmsg,"CUDA error %s: %s", msg, cudaGetErrorString(err));
    NAMD_die(errmsg);
  }
}

char *devicelist;
static int usedevicelist;
static int ignoresharing;

void cuda_getargs(char **argv) {
  devicelist = 0;
  usedevicelist = CmiGetArgStringDesc(argv, "+devices", &devicelist,
	"comma-delimited list of CUDA device numbers such as 0,2,1,2");
  ignoresharing = CmiGetArgFlag(argv, "+ignoresharing");
}

static int shared_gpu;
static int first_pe_sharing_gpu;
static int next_pe_sharing_gpu;

static int gpu_is_mine;

void cuda_initialize() {

  char host[128];
#ifdef NOHOSTNAME
  sprintf(host,"unknown");
#else
  gethostname(host, 128);  host[127] = 0;
#endif

  int myPhysicalNodeID = CmiPhysicalNodeID(CkMyPe());
  int myRankInPhysicalNode;
  int numPesOnPhysicalNode;
  int *pesOnPhysicalNode;
  CmiGetPesOnPhysicalNode(myPhysicalNodeID,
                           &pesOnPhysicalNode,&numPesOnPhysicalNode);
  {
    int i;
    for ( i=0; i < numPesOnPhysicalNode; ++i ) {
      if ( i && (pesOnPhysicalNode[i] <= pesOnPhysicalNode[i-1]) ) {
        i = numPesOnPhysicalNode;
        break;
      }
      if ( pesOnPhysicalNode[i] == CkMyPe() ) break;
    }
    if ( i == numPesOnPhysicalNode || i != CmiPhysicalRank(CkMyPe()) ) {
      CkPrintf("Bad result from CmiGetPesOnPhysicalNode!\n");
      for ( i=0; i < numPesOnPhysicalNode; ++i ) {
        CkPrintf("pe %d physnode rank %d of %d is %d\n", CkMyPe(),
          i, numPesOnPhysicalNode, pesOnPhysicalNode[i]);
      }
      myRankInPhysicalNode = 0;
      numPesOnPhysicalNode = 1;
      pesOnPhysicalNode = new int[1];
      pesOnPhysicalNode[0] = CkMyPe();
    } else {
      myRankInPhysicalNode = i;
    }
  }
  // CkPrintf("Pe %d ranks %d in physical node\n",CkMyPe(),myRankInPhysicalNode);

  int deviceCount = 0;
  cudaGetDeviceCount(&deviceCount);
  if ( deviceCount <= 0 ) {
    NAMD_die("No CUDA devices found.");
  }

  int *devices;
  int ndevices = 0;
  int nexclusive = 0;
  if ( usedevicelist ) {
    devices = new int[strlen(devicelist)];
    int i = 0;
    while ( devicelist[i] ) {
      ndevices += sscanf(devicelist+i,"%d",devices+ndevices);
      while ( devicelist[i] && isdigit(devicelist[i]) ) ++i;
      while ( devicelist[i] && ! isdigit(devicelist[i]) ) ++i;
    }
  } else {
    if ( ! CkMyPe() ) {
      CkPrintf("Did not find +devices i,j,k,... argument, using all\n");
    }
    devices = new int[deviceCount];
    for ( int i=0; i<deviceCount; ++i ) {
      int dev = i % deviceCount;
#if CUDA_VERSION >= 2020
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, dev);
      if ( deviceProp.computeMode != cudaComputeModeProhibited
           && deviceProp.multiProcessorCount > 2 ) {  // exclude weak cards
        devices[ndevices++] = dev;
      }
      if ( deviceProp.computeMode == cudaComputeModeExclusive ) {
        ++nexclusive;
      }
#else
      devices[ndevices++] = dev;
#endif
    }
  }

  if ( ! ndevices ) {
    NAMD_die("All CUDA devices are in prohibited mode.");
  }

  shared_gpu = 0;
  gpu_is_mine = 1;
  first_pe_sharing_gpu = CkMyPe();
  next_pe_sharing_gpu = CkMyPe();

 if ( (ndevices >= numPesOnPhysicalNode) || (nexclusive == 0) ) {

  int dev;
  if ( numPesOnPhysicalNode > 1 ) {
    dev = devices[myRankInPhysicalNode % ndevices];
    if ( ! ignoresharing ) {
     for ( int i = (myRankInPhysicalNode + 1) % numPesOnPhysicalNode;
          i != myRankInPhysicalNode;
          i = (i + 1) % numPesOnPhysicalNode ) {
      if (devices[i % ndevices] == dev) {
        shared_gpu = 1;
        next_pe_sharing_gpu = pesOnPhysicalNode[i];
        break;
      }
     }
    }
    if ( shared_gpu ) {
      for ( int i = 0; i < numPesOnPhysicalNode; ++i ) {
        if (devices[i % ndevices] == dev) {
          first_pe_sharing_gpu = pesOnPhysicalNode[i];
          break;
        }
      }
      CkPrintf("Pe %d sharing CUDA device %d first %d next %d\n",
		CkMyPe(), dev, first_pe_sharing_gpu, next_pe_sharing_gpu);
    }
  } else {  // in case phys node code is lying
    dev = devices[CkMyPe() % ndevices];
  }

  // disable token-passing but don't submit local until remote finished
  // if shared_gpu is true, otherwise submit all work immediately
  first_pe_sharing_gpu = CkMyPe();
  next_pe_sharing_gpu = CkMyPe();

  gpu_is_mine = ( first_pe_sharing_gpu == CkMyPe() ); 

  if ( dev >= deviceCount ) {
    char buf[256];
    sprintf(buf,"Pe %d unable to bind to CUDA device %d on %s because only %d devices are present",
		CkMyPe(), dev, host, deviceCount);
    NAMD_die(buf);
  }

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, dev);
  CkPrintf("Pe %d physical rank %d binding to CUDA device %d on %s: '%s'  Mem: %dMB  Rev: %d.%d\n",
             CkMyPe(), myRankInPhysicalNode, dev, host,
             deviceProp.name, deviceProp.totalGlobalMem / (1024*1024),
             deviceProp.major, deviceProp.minor);

  cudaSetDevice(dev);
  cudaError_t err;
  if ((err = cudaGetLastError()) != cudaSuccess) {
    char errmsg[1024];
    sprintf(errmsg,"CUDA error binding to device %d on pe %d: %s",
			dev, CkMyPe(), cudaGetErrorString(err));
    NAMD_die(errmsg);
  }

 }  // just let CUDA pick a device for us

  if ( sizeof(patch_pair) & 15 ) NAMD_die("sizeof(patch_pair) % 16 != 0");
  if ( sizeof(force_list) & 15 ) NAMD_die("sizeof(force_list) % 16 != 0");
  if ( sizeof(atom) & 15 ) NAMD_die("sizeof(atom) % 16 != 0");
  if ( sizeof(atom_param) & 15 ) NAMD_die("sizeof(atom_param) % 16 != 0");

  cuda_init();

}


void build_cuda_force_table() {
  ComputeNonbondedCUDA::build_force_table();
}

void ComputeNonbondedCUDA::build_force_table() {  // static

  float4 t[FORCE_TABLE_SIZE];

  const BigReal r2_delta = ComputeNonbondedUtil:: r2_delta;
  const int r2_delta_exp = ComputeNonbondedUtil:: r2_delta_exp;
  // const int r2_delta_expc = 64 * (r2_delta_exp - 127);
  const int r2_delta_expc = 64 * (r2_delta_exp - 1023);

  double r2list[FORCE_TABLE_SIZE];  // double to match cpu code
  for ( int i=1; i<FORCE_TABLE_SIZE; ++i ) {
    double r = ((double) FORCE_TABLE_SIZE) / ( (double) i + 0.5 );
    r2list[i] = r*r + r2_delta;
  }

  union { double f; int32 i[2]; } byte_order_test;
  byte_order_test.f = 1.0;  // should occupy high-order bits only
  int32 *r2iilist = (int32*)r2list + ( byte_order_test.i[0] ? 0 : 1 );

  for ( int i=1; i<FORCE_TABLE_SIZE; ++i ) {
    double r = ((double) FORCE_TABLE_SIZE) / ( (double) i + 0.5 );
    int table_i = (r2iilist[2*i] >> 14) + r2_delta_expc;  // table_i >= 0

    if ( r > cutoff ) {
      t[i].x = 0.;
      t[i].y = 0.;
      t[i].z = 0.;
      t[i].w = 0.;
      continue;
    }

    BigReal diffa = r2list[i] - r2_table[table_i];

    // coulomb 1/r or fast force
    // t[i].x = 1. / (r2 * r);  // -1/r * d/dr r^-1
    {
      // BigReal table_a = fast_table[4*table_i];
      BigReal table_b = fast_table[4*table_i+1];
      BigReal table_c = fast_table[4*table_i+2];
      BigReal table_d = fast_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].x = 2. * grad;
    }


    // pme correction for slow force
    // t[i].w = 0.;
    {
      // BigReal table_a = scor_table[4*table_i];
      BigReal table_b = scor_table[4*table_i+1];
      BigReal table_c = scor_table[4*table_i+2];
      BigReal table_d = scor_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].w = 2. * grad;
    }


    // vdw 1/r^6
    // t[i].y = 6. / (r8);  // -1/r * d/dr r^-6
    {
      // BigReal table_a = vdwb_table[4*table_i];
      BigReal table_b = vdwb_table[4*table_i+1];
      BigReal table_c = vdwb_table[4*table_i+2];
      BigReal table_d = vdwb_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].y = 2. * -1. * grad;
    }


    // vdw 1/r^12
    // t[i].z = 12e / (r8 * r4 * r2);  // -1/r * d/dr r^-12
    {
      // BigReal table_a = vdwa_table[4*table_i];
      BigReal table_b = vdwa_table[4*table_i+1];
      BigReal table_c = vdwa_table[4*table_i+2];
      BigReal table_d = vdwa_table[4*table_i+3];
      BigReal grad =
		( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
      t[i].z = 2. * grad;
    }

    // CkPrintf("%d %g %g %g %g %g %g\n", i, r, diffa,
    //   t[i].x, t[i].y, t[i].z, t[i].w);

/*
    double r2 = r * r;
    double r4 = r2 * r2;
    double r8 = r4 * r4;

    t[i].x = 1. / (r2 * r);  // -1/r * d/dr r^-1
    t[i].y = 6. / (r8);  // -1/r * d/dr r^-6
    t[i].z = 12. / (r8 * r4 * r2);  // -1/r * d/dr r^-12
    t[i].w = 0.;
*/
  }

  t[0].x = 0.f;
  t[0].y = 0.f;
  t[0].z = 0.f;
  t[0].w = 0.f;

  cuda_bind_force_table(t);

  if ( ! CkMyPe() ) {
    CkPrintf("Info: Updated CUDA force table with %d elements.\n", FORCE_TABLE_SIZE);
  }
}

static ComputeNonbondedCUDA* cudaCompute = 0;
static ComputeMgr *computeMgr = 0;

static ushort2 *exclusionsByAtom;

void ComputeNonbondedCUDA::build_exclusions() {
  Molecule *mol = Node::Object()->molecule;
  int natoms = mol->numAtoms; 
  exclusionsByAtom = new ushort2[natoms];

  // create unique sorted lists

  ObjectArena<int32> listArena;
  ResizeArray<int32*> unique_lists;
  SortableResizeArray<int32> curList;
  int totalbits = 0;
  for ( int i=0; i<natoms; ++i ) {
    const int32 *mol_list = mol->get_full_exclusions_for_atom(i);
    curList.resize(0);
    int n = mol_list[0] + 1;
    curList.add(0);  // always excluded from self
    for ( int j=1; j<n; ++j ) { curList.add(mol_list[j] - i); }
    curList.sort();

    int j;
    for ( j=0; j<unique_lists.size(); ++j ) {
      if ( n != unique_lists[j][0] ) continue;  // no match
      int k;
      for ( k=0; k<n; ++k ) {
        if ( unique_lists[j][k+3] != curList[k] ) break;
      }
      if ( k == n ) break;  // found match
    }
    if ( j == unique_lists.size() ) {  // no match
      int32 *list = listArena.getNewArray(n+3);
      list[0] = n;
      int maxdiff = 0;
      maxdiff = -1 * curList[0];
      if ( curList[n-1] > maxdiff ) maxdiff = curList[n-1];
      list[1] = maxdiff;
      list[2] = totalbits + maxdiff;
      totalbits += 2*maxdiff + 1;
      for ( int k=0; k<n; ++k ) {
        list[k+3] = curList[k];
      }
      unique_lists.add(list);
    }
    exclusionsByAtom[i].x = unique_lists[j][1];  // maxdiff
    exclusionsByAtom[i].y = unique_lists[j][2];  // start
  }

  if ( totalbits & 31 ) totalbits += ( 32 - ( totalbits & 31 ) );

  if ( ! CkMyPe() ) {
  CkPrintf("Info: Found %d unique exclusion lists needing %d bytes\n",
		unique_lists.size(), totalbits / 8);
  }

#define SET_EXCL(EXCL,BASE,DIFF) \
         (EXCL)[((BASE)+(DIFF))>>5] |= (1<<(((BASE)+(DIFF))&31))

  unsigned int *exclusion_bits = new unsigned int[totalbits/32];
  memset(exclusion_bits, 0, totalbits/8);

  int base = 0;
  for ( int i=0; i<unique_lists.size(); ++i ) {
    base += unique_lists[i][1];
    if ( base != unique_lists[i][2] ) {
      NAMD_bug("ComputeNonbondedCUDA::build_exclusions base != stored");
    }
    int n = unique_lists[i][0];
    for ( int j=0; j<n; ++j ) {
      SET_EXCL(exclusion_bits,base,unique_lists[i][j+3]);
    }
    base += unique_lists[i][1] + 1;
  }

  cuda_bind_exclusions(exclusion_bits, totalbits/32);

  delete [] exclusion_bits;
}


void register_cuda_compute_self(ComputeID c, PatchID pid) {

  if ( ! cudaCompute ) NAMD_die("register_self called early");

  cudaCompute->requirePatch(pid);

  ComputeNonbondedCUDA::compute_record cr;
  cr.c = c;
  cr.pid[0] = pid;  cr.pid[1] = pid;
  cr.offset = 0.;
  if ( cudaCompute->patchMap->node(pid) == CkMyPe() ) {
    cudaCompute->localComputeRecords.add(cr);
  } else {
    cudaCompute->remoteComputeRecords.add(cr);
  }
}

void register_cuda_compute_pair(ComputeID c, PatchID pid[], int t[]) {

  if ( ! cudaCompute ) NAMD_die("register_pair called early");
 
  cudaCompute->requirePatch(pid[0]);
  cudaCompute->requirePatch(pid[1]);

  ComputeNonbondedCUDA::compute_record cr, cr2;
  cr.c = c;  cr2.c = c;
  cr.pid[0] = pid[0];  cr.pid[1] = pid[1];
  cr2.pid[0] = pid[1];  cr2.pid[1] = pid[0];

  int t1 = t[0];
  int t2 = t[1];
  Vector offset = cudaCompute->patchMap->center(pid[0])
                - cudaCompute->patchMap->center(pid[1]);
  offset.x += (t1%3-1) - (t2%3-1);
  offset.y += ((t1/3)%3-1) - ((t2/3)%3-1);
  offset.z += (t1/9-1) - (t2/9-1);
  cr.offset = offset;
  cr2.offset = -1. * offset;
    
  if ( cudaCompute->patchMap->node(pid[0]) == CkMyPe() ) {
    cudaCompute->localComputeRecords.add(cr);
  } else {
    cudaCompute->remoteComputeRecords.add(cr);
  }
  if ( cudaCompute->patchMap->node(pid[1]) == CkMyPe() ) {
    cudaCompute->localComputeRecords.add(cr2);
  } else {
    cudaCompute->remoteComputeRecords.add(cr2);
  }
}

void unregister_cuda_compute(ComputeID c) {  // static

  NAMD_die("unregister_compute unimplemented");

}

static int atomsChanged = 0;
static int computesChanged = 0;

static cudaEvent_t start_upload;
static cudaEvent_t start_calc;
static cudaEvent_t end_remote_calc;
static cudaEvent_t end_remote_download;
static cudaEvent_t end_local_calc;
static cudaEvent_t end_local_download;

ComputeNonbondedCUDA::ComputeNonbondedCUDA(ComputeID c, ComputeMgr *mgr) : Compute(c) {
  // CkPrintf("create ComputeNonbondedCUDA\n");
  cudaCompute = this;
  computeMgr = mgr;
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();
  reduction = 0;
  build_exclusions();

  SimParameters *params = Node::Object()->simParameters;
  if (params->pressureProfileOn) {
    NAMD_die("pressure profile not supported in CUDA");
  }

  atomsChanged = 1;
  computesChanged = 1;
  workStarted = 0;
  basePriority = PROXY_DATA_PRIORITY;

  cudaEventCreate(&start_upload);
  cudaEventCreate(&start_calc);
  cudaEventCreate(&end_remote_calc);
  cudaEventCreate(&end_remote_download);
  cudaEventCreate(&end_local_calc);
  cudaEventCreate(&end_local_download);
}


ComputeNonbondedCUDA::~ComputeNonbondedCUDA() { ; }

void ComputeNonbondedCUDA::requirePatch(int pid) {

  if ( ! reduction ) {
    reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  }
  computesChanged = 1;
  patch_record &pr = patchRecords.item(pid);
  if ( pr.refCount == 0 ) {
    pr.isLocal = ( patchMap->node(pid) == CkMyPe() );
    if ( pr.isLocal ) {
      localActivePatches.add(pid);
    } else {
      remoteActivePatches.add(pid);
    }
    activePatches.add(pid);
    setNumPatches(activePatches.size());
    pr.patchID = pid;
    pr.p = patchMap->patch(pid);
    pr.positionBox = pr.p->registerPositionPickup(cid);
    pr.forceBox = pr.p->registerForceDeposit(cid);
    pr.x = NULL;
    pr.xExt = NULL;
    pr.r = NULL;
    pr.f = NULL;
  }
  pr.refCount += 1;
}

static ResizeArray<patch_pair> patch_pairs;
static ResizeArray<force_list> force_lists;

static int num_atom_records_allocated;
static atom_param* atom_params;
static atom* atoms;
static float4* forces;
static float4* slow_forces;

static int cuda_timer_count;
static double cuda_timer_total;
static double kernel_time;

static int ccd_index_remote_download;
static int ccd_index_local_download;
#define CUDA_CONDITION CcdPERIODIC

void cuda_check_remote_progress(void *arg, double) {
  if ( cudaEventQuery(end_remote_download) == cudaSuccess ) {
    // ((ComputeNonbondedCUDA *) arg)->finishWork();
    WorkDistrib::messageEnqueueWork((ComputeNonbondedCUDA *) arg);
  } else {
    ccd_index_remote_download = CcdCallOnCondition(CUDA_CONDITION, cuda_check_remote_progress, arg);
  }
}

void cuda_check_local_progress(void *arg, double) {
  if ( cudaEventQuery(end_local_download) == cudaSuccess ) {
    kernel_time += CkWallTimer();
    WorkDistrib::messageEnqueueWork((ComputeNonbondedCUDA *) arg);
  } else {
    ccd_index_local_download = CcdCallOnCondition(CUDA_CONDITION, cuda_check_local_progress, arg);
  }
}

#if 0
// don't use this one unless timer is part of stream, above is better
void cuda_check_progress(void *arg, double) {
  if ( cuda_stream_finished() ) {
    kernel_time += CkWallTimer();
    CcdCallOnCondition(CUDA_CONDITION, ccd_index);
    // ((ComputeNonbondedCUDA *) arg)->finishWork();
    WorkDistrib::messageEnqueueWork((ComputeNonbondedCUDA *) arg);
  }
}
#endif

void ComputeNonbondedCUDA::atomUpdate() { atomsChanged = 1; }

static int kernel_launch_state = 0;

void ComputeNonbondedCUDA::doWork() {

// CkPrintf("Pe %d doWork %d\n", CkMyPe(), workStarted);

  if ( workStarted ) {
    if ( finishWork() ) {  // finished
      workStarted = 0;
      basePriority = PROXY_DATA_PRIORITY;  // higher to aid overlap
    } else {  // need to call again
      workStarted = 2;
      basePriority = COMPUTE_HOME_PRIORITY;  // lower for local
      if ( kernel_launch_state > 2 ) ccd_index_local_download = CcdCallOnCondition(CUDA_CONDITION,cuda_check_local_progress,this);
    }
    return;
  }
  workStarted = 1;
  basePriority = COMPUTE_PROXY_PRIORITY;

  Molecule *mol = Node::Object()->molecule;
  Parameters *params = Node::Object()->parameters;
  SimParameters *simParams = Node::Object()->simParameters;

  if ( simParams->useFlexibleCell ) {
    NAMD_die("flexible-cell constant pressure not supported in CUDA");
  }

  for ( int i=0; i<activePatches.size(); ++i ) {
    patch_record &pr = patchRecords[activePatches[i]];
    pr.x = pr.positionBox->open();
    pr.xExt = pr.p->getCompAtomExtInfo();
  }

 if ( atomsChanged || computesChanged ) {
  int npatches = activePatches.size();

  if ( computesChanged ) {
    computesChanged = 0;

    int *ap = activePatches.begin();
    for ( int i=0; i<localActivePatches.size(); ++i ) {
      *(ap++) = localActivePatches[i];
    }
    for ( int i=0; i<remoteActivePatches.size(); ++i ) {
      *(ap++) = remoteActivePatches[i];
    }

    int nlc = localComputeRecords.size();
    int nrc = remoteComputeRecords.size();
    computeRecords.resize(nlc+nrc);
    compute_record *cr = computeRecords.begin();
    for ( int i=0; i<nlc; ++i ) {
      *(cr++) = localComputeRecords[i];
    }
    for ( int i=0; i<nrc; ++i ) {
      *(cr++) = remoteComputeRecords[i];
    }

    force_lists.resize(npatches);
    for ( int i=0; i<npatches; ++i ) {
      patchRecords[activePatches[i]].localIndex = i;
      force_lists[i].force_list_size = 0;
    }

    int ncomputes = computeRecords.size();
    patch_pairs.resize(ncomputes);
    for ( int i=0; i<ncomputes; ++i ) {
      ComputeNonbondedCUDA::compute_record &cr = computeRecords[i];
      int lp1 = patchRecords[cr.pid[0]].localIndex;
      int lp2 = patchRecords[cr.pid[1]].localIndex;
      force_lists[lp1].force_list_size++;
      patch_pair &pp = patch_pairs[i];
      pp.offset.x = cr.offset.x;
      pp.offset.y = cr.offset.y;
      pp.offset.z = cr.offset.z;
    }

   if ( simParams->outputCudaTiming ) {
    CkPrintf("Pe %d has %d local and %d remote patches and %d local and %d remote computes.\n",
	CkMyPe(), localActivePatches.size(), remoteActivePatches.size(),
	localComputeRecords.size(), remoteComputeRecords.size());
   }
  }

  int istart = 0;
  int flstart = 0;
  int max_atoms_per_patch = 0;
  int i;
  for ( i=0; i<npatches; ++i ) {
    if ( i == localActivePatches.size() ) {
       num_local_atom_records = istart;
    }
    force_lists[i].force_list_start = flstart;
    force_lists[i].force_output_start = istart;
    patch_record &pr = patchRecords[activePatches[i]];
    pr.localStart = istart;
    int natoms = pr.p->getNumAtoms();
    int nfreeatoms = natoms;
    if ( fixedAtomsOn ) {
      const CompAtomExt *aExt = pr.xExt;
      for ( int j=0; j<natoms; ++j ) {
        if ( aExt[j].atomFixed ) --nfreeatoms;
      }
    }
    if ( natoms > max_atoms_per_patch ) max_atoms_per_patch = natoms;
    pr.numAtoms = natoms;
    pr.numFreeAtoms = nfreeatoms;
    if ( natoms & 15 ) { natoms += 16 - (natoms & 15); }
    if ( nfreeatoms & 15 ) { nfreeatoms += 16 - (nfreeatoms & 15); }
    force_lists[i].patch_size = nfreeatoms;
    flstart += nfreeatoms * force_lists[i].force_list_size;
    istart += natoms;  // already rounded up
    force_lists[i].force_list_size = 0;  // rebuild below
  }
  if ( i == localActivePatches.size() ) {
     num_local_atom_records = istart;
  }
  num_force_records = flstart;
  num_atom_records = istart;
  num_remote_atom_records = num_atom_records - num_local_atom_records;
  if ( num_atom_records > num_atom_records_allocated ) {
    if ( num_atom_records_allocated ) {
      cudaFreeHost(atom_params);
      cudaFreeHost(atoms);
      cudaFreeHost(forces);
      cudaFreeHost(slow_forces);
    }
    num_atom_records_allocated = 1.1 * num_atom_records + 1;
    cudaMallocHost((void**)&atom_params,sizeof(atom_param)*num_atom_records_allocated);
    cudaMallocHost((void**)&atoms,sizeof(atom)*num_atom_records_allocated);
    cudaMallocHost((void**)&forces,sizeof(float4)*num_atom_records_allocated);
    cudaMallocHost((void**)&slow_forces,sizeof(float4)*num_atom_records_allocated);
  }

#if 0
  if ( max_atoms_per_patch > MAX_ATOMS_PER_PATCH ) {
    char errstr[1024];
    sprintf(errstr,"Found patch with %d atoms; limit is %d for CUDA."
	"  Try enabling twoAwayX, twoAwayY, and/or twoAwayZ to fix this.",
	max_atoms_per_patch, MAX_ATOMS_PER_PATCH);
    NAMD_die(errstr);
  }
#endif

  int ncomputes = computeRecords.size();
  for ( int i=0; i<ncomputes; ++i ) {
    ComputeNonbondedCUDA::compute_record &cr = computeRecords[i];
    int p1 = cr.pid[0];
    int p2 = cr.pid[1];
    int lp1 = patchRecords[p1].localIndex;
    int lp2 = patchRecords[p2].localIndex;
    patch_pair &pp = patch_pairs[i];
    pp.patch1_atom_start = patchRecords[p1].localStart;
    pp.patch1_force_start = force_lists[lp1].force_list_start +
	force_lists[lp1].patch_size * force_lists[lp1].force_list_size;
    pp.patch1_size = patchRecords[p1].numAtoms;
    pp.patch1_force_size = patchRecords[p1].numFreeAtoms;
    // istart += pp.patch1_size;
    // if ( istart & 15 ) { istart += 16 - (istart & 15); }
    pp.patch2_atom_start = patchRecords[p2].localStart;
    // pp.patch2_force_start = force_lists[lp2].force_list_start +
// 	force_lists[lp2].patch_size * force_lists[lp2].force_list_size;
    pp.patch2_size = patchRecords[p2].numAtoms;
    // pp.patch2_force_size = patchRecords[p2].numFreeAtoms;
    // this must happen at end to get self computes right
    force_lists[lp1].force_list_size++;
    // if ( lp1 != lp2 || cr.t[0] != cr.t[1] ) {
    //   force_lists[lp2].force_list_size++;
    // }
    // istart += pp.patch2_size;
    // if ( istart & 15 ) { istart += 16 - (istart & 15); }
  }

#if 0
  CkPrintf("Pe %d cuda_bind_patch_pairs %d %d %d %d %d\n", CkMyPe(),
	patch_pairs.size(), force_lists.size(),
	num_atom_records, num_force_records,
	max_atoms_per_patch);
#endif

  int totalmem = patch_pairs.size() * sizeof(patch_pair) +
                force_lists.size() * sizeof(force_list) +
                num_force_records * sizeof(float4) +
                num_atom_records * sizeof(atom) +
                num_atom_records * sizeof(atom_param) +
                num_atom_records * sizeof(float4);
  int totalcopy = num_atom_records * ( sizeof(atom) + sizeof(float4) );
/*
  CkPrintf("Pe %d allocating %d MB of GPU memory, will copy %d kB per step\n",
			CkMyPe(), totalmem >> 20, totalcopy >> 10);
*/

  cuda_bind_patch_pairs(patch_pairs.begin(), patch_pairs.size(),
			force_lists.begin(), force_lists.size(),
			num_atom_records, num_force_records,
			max_atoms_per_patch);

 }  // atomsChanged || computesChanged

  double charge_scaling = sqrt(COLOUMB * scaling * dielectric_1);


  for ( int i=0; i<activePatches.size(); ++i ) {
    patch_record &pr = patchRecords[activePatches[i]];

    int start = pr.localStart;
    int n = pr.numAtoms;
    const CompAtom *a = pr.x;
    const CompAtomExt *aExt = pr.xExt;
    if ( atomsChanged ) {
      atom_param *ap = atom_params + start;
      int k = 0;
      for ( int j=0; j<n; ++j ) {
        // put free atoms first
        if ( fixedAtomsOn && aExt[j].atomFixed ) continue;
        ap[k].index = aExt[j].id;
        ap[k].excl_index = exclusionsByAtom[aExt[j].id].y;
        ap[k].excl_maxdiff = exclusionsByAtom[aExt[j].id].x;
        int vdwtype = mol->atomvdwtype(aExt[j].id);
        Real sig, eps, sig14, eps14;
        params->get_vdw_params(&sig,&eps,&sig14,&eps14,vdwtype);
        ap[k].sqrt_epsilon = sqrt(4.0 * scaling * eps);
        ap[k].half_sigma = 0.5 * sig;
// CkPrintf("%d %d %f %f\n", k, vdwtype, sig, eps);
        ++k;
      }
      if ( fixedAtomsOn ) for ( int j=0; j<n; ++j ) {
        // put fixed atoms at end
        if ( ! (fixedAtomsOn && aExt[j].atomFixed) ) continue;
        ap[k].index = aExt[j].id;
        ap[k].excl_index = exclusionsByAtom[aExt[j].id].y;
        ap[k].excl_maxdiff = exclusionsByAtom[aExt[j].id].x;
        int vdwtype = mol->atomvdwtype(aExt[j].id);
        Real sig, eps, sig14, eps14;
        params->get_vdw_params(&sig,&eps,&sig14,&eps14,vdwtype);
        ap[k].sqrt_epsilon = sqrt(4.0 * scaling * eps);
        ap[k].half_sigma = 0.5 * sig;
// CkPrintf("%d %d %f %f\n", k, vdwtype, sig, eps);
        ++k;
      }
    }
    {
      Vector center =
        pr.p->flags.lattice.unscale(cudaCompute->patchMap->center(pr.patchID));
      atom *ap = atoms + start;
      int k = 0;
      for ( int j=0; j<n; ++j ) {
        // put free atoms first
        if ( fixedAtomsOn && aExt[j].atomFixed ) continue;
        ap[k].position.x = a[j].position.x - center.x;
        ap[k].position.y = a[j].position.y - center.y;
        ap[k].position.z = a[j].position.z - center.z;
        ap[k].charge = charge_scaling * a[j].charge;
        ++k;
      }
      if ( fixedAtomsOn ) for ( int j=0; j<n; ++j ) {
        // put fixed atoms at end
        if ( ! (fixedAtomsOn && aExt[j].atomFixed) ) continue;
        ap[k].position.x = a[j].position.x - center.x;
        ap[k].position.y = a[j].position.y - center.y;
        ap[k].position.z = a[j].position.z - center.z;
        ap[k].charge = charge_scaling * a[j].charge;
        ++k;
      }
    }
  }

  kernel_time = -1. * CkWallTimer();
#if 0
  kernel_launch_state = 3;

  cudaEventRecord(start_upload, stream);

  if ( atomsChanged ) {
    cuda_bind_atom_params(atom_params);
  }

  XXX - THIS PATH NOT UPDATED FOR SLOW FORCES - DO NOT COMPILE
  cuda_bind_atoms(atoms);
  cudaEventRecord(start_calc, stream);
  cuda_nonbonded_forces(cutoff2,
	localComputeRecords.size(),remoteComputeRecords.size(),
	localActivePatches.size(),remoteActivePatches.size());
  cudaEventRecord(end_remote_calc, stream);
  cuda_load_forces(forces,num_local_atom_records,num_remote_atom_records);
  cudaEventRecord(end_remote_download, stream);
  cuda_nonbonded_forces(cutoff2,
	0,localComputeRecords.size(),
	0,localActivePatches.size());
  cudaEventRecord(end_local_calc, stream);
  cuda_load_forces(forces,0,num_local_atom_records);
  cudaEventRecord(end_local_download, stream);

  if ( cuda_stream_finished() ) {
    CkPrintf("CUDA not overlapping with CPU work.\n");
  }

  // finishWork();

  ccd_index_remote_download = CcdCallOnCondition(CUDA_CONDITION,cuda_check_remote_progress,this);
#else

  kernel_launch_state = 1;
  if ( gpu_is_mine ) recvYieldDevice(-1);

#endif
}

static int ccd_index_remote_calc;
void cuda_check_remote_calc(void *arg, double) {
  // in theory we only need end_remote_calc, but overlap isn't reliable
  // if ( cudaEventQuery(end_remote_calc) == cudaSuccess ) {
  if ( cudaEventQuery(end_remote_download) == cudaSuccess ) {
// CkPrintf("Pe %d yielding to %d after remote calc\n", CkMyPe(), next_pe_sharing_gpu);
    computeMgr->sendYieldDevice(next_pe_sharing_gpu);
// CkPrintf("Pe %d yielded to %d after remote calc\n", CkMyPe(), next_pe_sharing_gpu);
  } else {
    ccd_index_remote_calc = CcdCallOnCondition(CUDA_CONDITION, cuda_check_remote_calc, arg);
  }
}

static int ccd_index_local_calc;
void cuda_check_local_calc(void *arg, double) {
  // in theory we only need end_local_calc, but overlap isn't reliable
  // if ( cudaEventQuery(end_local_calc) == cudaSuccess ) {
  if ( cudaEventQuery(end_local_download) == cudaSuccess ) {
// CkPrintf("Pe %d yielding to %d after local calc\n", CkMyPe(), next_pe_sharing_gpu);
    computeMgr->sendYieldDevice(next_pe_sharing_gpu);
// CkPrintf("Pe %d yielded to %d after local calc\n", CkMyPe(), next_pe_sharing_gpu);
  } else {
    ccd_index_local_calc = CcdCallOnCondition(CUDA_CONDITION, cuda_check_local_calc, arg);
  }
}

// computeMgr->sendYieldDevice(next_pe_sharing_gpu);

void ComputeNonbondedCUDA::recvYieldDevice(int pe) {

// CkPrintf("Pe %d entering state %d via yield from pe %d\n",
//           CkMyPe(), kernel_launch_state, pe);

  Flags &flags = patchRecords[activePatches[0]].p->flags;
  int doSlow = flags.doFullElectrostatics;

  Lattice &lattice = flags.lattice;
  float3 lata, latb, latc;
  lata.x = lattice.a().x;
  lata.y = lattice.a().y;
  lata.z = lattice.a().z;
  latb.x = lattice.b().x;
  latb.y = lattice.b().y;
  latb.z = lattice.b().z;
  latc.x = lattice.c().x;
  latc.y = lattice.c().y;
  latc.z = lattice.c().z;

  switch ( kernel_launch_state ) {
  case 1:
    ++kernel_launch_state;
    gpu_is_mine = 0;
    cudaEventRecord(start_upload, stream);

    if ( atomsChanged ) {
      cuda_bind_atom_params(atom_params);
    }

    cuda_bind_atoms(atoms);
    cudaEventRecord(start_calc, stream);
    cuda_nonbonded_forces(lata, latb, latc, cutoff2,
	localComputeRecords.size(),remoteComputeRecords.size(),
	localActivePatches.size(),remoteActivePatches.size(), doSlow);
    cudaEventRecord(end_remote_calc, stream);
    cuda_load_forces(forces, (doSlow ? slow_forces : 0 ),
        num_local_atom_records,num_remote_atom_records);
    cudaEventRecord(end_remote_download, stream);
    ccd_index_remote_download = CcdCallOnCondition(CUDA_CONDITION,cuda_check_remote_progress,this);
    if ( shared_gpu ) {
      ccd_index_remote_calc = CcdCallOnCondition(CUDA_CONDITION,cuda_check_remote_calc,this);
      break;
    }
 
  case 2:
    ++kernel_launch_state;
    gpu_is_mine = 0;
    cuda_nonbonded_forces(lata, latb, latc, cutoff2,
	0,localComputeRecords.size(),
	0,localActivePatches.size(), doSlow);
    cudaEventRecord(end_local_calc, stream);
    cuda_load_forces(forces, (doSlow ? slow_forces : 0 ),
        0,num_local_atom_records);
    cudaEventRecord(end_local_download, stream);
    if ( workStarted == 2 ) ccd_index_local_download = CcdCallOnCondition(CUDA_CONDITION,cuda_check_local_progress,this);
    if ( shared_gpu ) {
      ccd_index_local_calc = CcdCallOnCondition(CUDA_CONDITION,cuda_check_local_calc,this);
      break;
    }

  default:
    gpu_is_mine = 1;
    break;
  }

}


int ComputeNonbondedCUDA::finishWork() {

  cuda_errcheck("cuda stream completed");

  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;

  Flags &flags = patchRecords[activePatches[0]].p->flags;
  int doSlow = flags.doFullElectrostatics;

  ResizeArray<int> &patches( workStarted == 1 ?
				remoteActivePatches : localActivePatches );

  for ( int i=0; i<patches.size(); ++i ) {
    patch_record &pr = patchRecords[patches[i]];
    pr.r = pr.forceBox->open();
    pr.f = pr.r->f[Results::nbond];
  }

  // long long int wcount = 0;
  double virial = 0.;
  double virial_slow = 0.;

  for ( int i=0; i<patches.size(); ++i ) {
    patch_record &pr = patchRecords[patches[i]];
    int start = pr.localStart;
    int n = pr.numAtoms;
    Force *f = pr.f;
    Force *f_slow = pr.r->f[Results::slow];
    const CompAtom *a = pr.x;
    const CompAtomExt *aExt = pr.xExt;
    float4 *af = forces + start;
    float4 *af_slow = slow_forces + start;
    int k = 0;
    for ( int j=0; j<n; ++j ) {
      // only free atoms return forces
      if ( fixedAtomsOn && aExt[j].atomFixed ) continue;
      f[j].x += af[k].x;
      f[j].y += af[k].y;
      f[j].z += af[k].z;
      // wcount += af[k].w;
      virial += af[k].w;
      if ( doSlow ) {
        f_slow[j].x += af_slow[k].x;
        f_slow[j].y += af_slow[k].y;
        f_slow[j].z += af_slow[k].z;
        virial_slow += af_slow[k].w;
      }
      ++k;
    }

#if 0
    // check exclusions reported as w
    if ( CkNumPes() == 1 ) {
      const CompAtomExt *aExt = pr.xExt;
      int k = 0;
      for ( int j=0; j<n; ++j ) {
        // only free atoms return forces
        if ( fixedAtomsOn && aExt[j].atomFixed ) continue;
        int excl_expected = mol->get_full_exclusions_for_atom(aExt[j].id)[0] + 1;
        if ( af[k].w != excl_expected ) {
          CkPrintf("%d:%d(%d) atom %d found %d exclusions but expected %d\n",
		i, j, k, aExt[j].id, (int)af[k].w, excl_expected );
        }
        ++k;
      }
    }
#endif
    

#if 0
    if ( i % 31 == 0 ) for ( int j=0; j<3; ++j ) {
      CkPrintf("Pe %d patch %d atom %d (%f %f %f) force %f\n", CkMyPe(), i,
	j, pr.x[j].position.x, pr.x[j].position.y, pr.x[j].position.z,
	af[j].w);
    }
#endif
    pr.positionBox->close(&(pr.x));
    pr.forceBox->close(&(pr.r));
  }

  virial *= (-1./6.);
  reduction->item(REDUCTION_VIRIAL_NBOND_XX) += virial;
  reduction->item(REDUCTION_VIRIAL_NBOND_YY) += virial;
  reduction->item(REDUCTION_VIRIAL_NBOND_ZZ) += virial;
  if ( doSlow ) {
    virial_slow *= (-1./6.);
    reduction->item(REDUCTION_VIRIAL_SLOW_XX) += virial_slow;
    reduction->item(REDUCTION_VIRIAL_SLOW_YY) += virial_slow;
    reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += virial_slow;
  }

  if ( workStarted == 1 ) return 0;  // not finished, call again

  atomsChanged = 0;
  reduction->submit();

  // int natoms = mol->numAtoms; 
  // double wpa = wcount;  wpa /= natoms;

  // CkPrintf("Pe %d CUDA kernel %f ms, total %f ms, wpa %f\n", CkMyPe(),
	// 	kernel_time * 1.e3, time * 1.e3, wpa);

  float upload_ms, remote_calc_ms, remote_download_ms;
  float local_calc_ms, local_download_ms, total_ms;
  cuda_errcheck("before event timers");
  cudaEventElapsedTime(&upload_ms, start_upload, start_calc);
  cuda_errcheck("event timer 1");
  cudaEventElapsedTime(&remote_calc_ms, start_calc, end_remote_calc);
  cuda_errcheck("event timer 2");
  cudaEventElapsedTime(&remote_download_ms, end_remote_calc, end_remote_download);
  cuda_errcheck("event timer 3");
  cudaEventElapsedTime(&local_calc_ms, end_remote_download, end_local_calc);
  cuda_errcheck("event timer 4");
  cudaEventElapsedTime(&local_download_ms, end_local_calc, end_local_download);
  cuda_errcheck("event timer 5");
  cudaEventElapsedTime(&total_ms, start_upload, end_local_download);
  cuda_errcheck("event timer 6");
  cuda_errcheck("event timers");

  cuda_timer_total += kernel_time;
  if ( simParams->outputCudaTiming &&
	cuda_timer_count % simParams->outputCudaTiming == 0 ) {
    cuda_timer_total /= (cuda_timer_count + 1);
    CkPrintf("CUDA EVENT TIMING: %d %f %f %f %f %f %f\n",
	CkMyPe(), upload_ms, remote_calc_ms, remote_download_ms,
			local_calc_ms, local_download_ms, total_ms);
    CkPrintf("CUDA TIMING: %f ms/step on node %d\n",
			cuda_timer_total * 1.e3, CkMyPe());
    cuda_timer_count = 0;
    cuda_timer_total = 0;
  }
  cuda_timer_count++;

  return 1;  // finished and ready for next step
}


#endif  // NAMD_CUDA

