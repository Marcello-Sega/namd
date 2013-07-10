
#include "common.h"
#include "charm++.h"

#include "WorkDistrib.h"
#include "ComputeMgr.h"
#include "ProxyMgr.h"
#include "ComputeNonbondedMIC.h"
#include "ComputeNonbondedMICKernel.h"
#include "LJTable.h"
#include "ObjectArena.h"
#include "SortAtoms.h"
#include <algorithm>

#include "NamdTypes.h"

// This file only used for builds with Xeon Phi (MIC) support
#ifdef NAMD_MIC

#include <offload.h>


//#define PRINT_GBIS
#undef PRINT_GBIS

#ifdef WIN32
  #define __thread __declspec(thread)
#endif


// Declare the structures that will be used for holding kernel-specific data since
//   virial data is located within that data structure (read below).
extern __thread __attribute__((target(mic))) mic_kernel_data *device__remote_kernel_data;
extern __thread __attribute__((target(mic))) mic_kernel_data *device__local_kernel_data;

#if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
  extern __thread __attribute__((target(mic))) patch_pair* device__patch_pairs;
  extern __thread __attribute__((target(mic))) int device__patch_pairs_size;
  extern __thread __attribute__((target(mic))) int device__force_lists_size;
  #if MIC_DEVICE_TRACING_DETAILED != 0
    extern __thread __attribute__((target(mic))) double* device__device_times_computes;
    extern __thread __attribute__((target(mic))) double* device__device_times_patches;
  #endif
  extern __thread __attribute__((target(mic))) double* device__device_times_start;
#endif


int mic_device_count = 0;
int mic_get_device_count() { return mic_device_count; }


void mic_errcheck(const char *msg) {
  // DMK - TODO : Add mechanism to detect errors here
}

void mic_die(const char *msg) {
  char host[128];
  #ifdef NOHOSTNAME
    sprintf(host,"physical node %d", CmiPhysicalNodeID(CkMyPe()));
  #else
    gethostname(host, 128);  host[127] = 0;
  #endif
  char devstr[128] = "";
  int devnum = mic_device_pe();
  if (devnum == CkMyPe()) {
   sprintf(devstr, " device %d", devnum);
  }
  char errmsg[1024];
  sprintf(errmsg,"MIC error on Pe %d (%s%s): %s", CkMyPe(), host, devstr, msg);
  NAMD_die(errmsg);
}


// Function used to process MIC-specific command-line parameters along with
//   global variables to told the content (or flag value) of the parameters.
char *devicelist;
static __thread int usedevicelist;
static __thread int ignoresharing;
static __thread int mergegrids;
void mic_getargs(char **argv) {
  devicelist = 0;
  usedevicelist = CmiGetArgStringDesc(argv, "+devices", &devicelist,
	"comma-delimited list of MIC device numbers such as 0,2,1,2");
  ignoresharing = CmiGetArgFlag(argv, "+ignoresharing");
  mergegrids = CmiGetArgFlag(argv, "+mergegrids");
}


// Variables for tracking which MIC-device is associated with which host PE.
static __thread int shared_gpu;
static __thread int first_pe_sharing_gpu;
static __thread int next_pe_sharing_gpu;
static __thread int devicePe;
static __thread int numPesSharingDevice;
static __thread int *pesSharingDevice;

static __thread int gpu_is_mine;
static __thread int myDevice;

int mic_device_pe() { return devicePe; }

bool mic_device_shared_with_pe(int pe) {
  for ( int i=0; i<numPesSharingDevice; ++i ) {
    if ( pesSharingDevice[i] == pe ) return true;
  }
  return false;
}

static inline bool sortop_bitreverse(int a, int b) {
  if ( a == b ) return 0; 
  for ( int bit = 1; bit; bit *= 2 ) {
    if ( (a&bit) != (b&bit) ) return ((a&bit) < (b&bit));
  }
  return 0;
}

// Projections support
#define MIC_EVENT_ID_BASE 100
#define MIC_TRACE_REMOTE(START,END) \
  do { int dev; micGetDevice(&dev); traceUserBracketEvent( \
       MIC_EVENT_ID_BASE + 2 * dev, START, END); } while (0)
#define MIC_TRACE_LOCAL(START,END) \
  do { int dev; micGetDevice(&dev); traceUserBracketEvent( \
       MIC_EVENT_ID_BASE + 2 * dev + 1, START, END); } while (0)


void mic_initialize() {

  // If this is PE 0, register the MIC-specific user events with Projections
  if ( 0 == CkMyPe() ) {
    MIC_TRACING_REGISTER_EVENTS
  }

  // Get the host name for this node
  char host[128];
  #ifdef NOHOSTNAME
    sprintf(host,"physical node %d", CmiPhysicalNodeID(CkMyPe()));
  #else
    gethostname(host, 128);  host[127] = 0;
  #endif

  int myPhysicalNodeID = CmiPhysicalNodeID(CkMyPe());
  int myRankInPhysicalNode;
  int numPesOnPhysicalNode;
  int *pesOnPhysicalNode;
  CmiGetPesOnPhysicalNode(myPhysicalNodeID,
                          &pesOnPhysicalNode,
                          &numPesOnPhysicalNode
                         );

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

  int deviceCount = _Offload_number_of_devices();
  if ( deviceCount < 1 ) {
    mic_die("No MIC devices found.");
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
    #pragma noinline
    for ( int i=0; i<deviceCount; ++i ) {
      int dev = i % deviceCount;
      #if 0
        int dev_ok = 0;
        #pragma offload target(mic: dev) in(dev)
        dev_ok = mic_check_local();
      #else
        int dev_ok = mic_check(dev);
      #endif
      if ( dev_ok ) devices[ndevices++] = dev;
      else CkPrintf("Offload to device %d on Pe %d failed.\n",dev,CkMyPe());
    }
  }

  mic_device_count = ndevices;

  if ( ! ndevices ) {
    mic_die("No usable MIC devices found.");
  }

  shared_gpu = 0;
  gpu_is_mine = 1;
  first_pe_sharing_gpu = CkMyPe();
  next_pe_sharing_gpu = CkMyPe();

  int dev;
  if ( numPesOnPhysicalNode > 1 ) {
    int myDeviceRank = myRankInPhysicalNode * ndevices / numPesOnPhysicalNode;
    dev = devices[myDeviceRank];
    devicePe = CkMyPe();
    if ( ignoresharing ) {
      pesSharingDevice = new int[1];
      pesSharingDevice[0] = CkMyPe();
      numPesSharingDevice = 1;
    } else {
      pesSharingDevice = new int[numPesOnPhysicalNode];
      devicePe = -1;
      numPesSharingDevice = 0;
      for ( int i = 0; i < numPesOnPhysicalNode; ++i ) {
        if ( i * ndevices / numPesOnPhysicalNode == myDeviceRank ) {
          int thisPe = pesOnPhysicalNode[i];
          pesSharingDevice[numPesSharingDevice++] = thisPe;
          if ( devicePe < 1 ) devicePe = thisPe;
          if ( sortop_bitreverse(thisPe,devicePe) ) devicePe = thisPe;
        }
      }
      for ( int j = 0; j < ndevices; ++j ) {
        if ( devices[j] == dev && j != myDeviceRank ) shared_gpu = 1;
      }
    }
    if ( shared_gpu && devicePe == CkMyPe() ) {
      CkPrintf("Pe %d sharing MIC device %d\n", CkMyPe(), dev);
    }
  } else {  // in case phys node code is lying
    dev = devices[CkMyPe() % ndevices];
    devicePe = CkMyPe();
    pesSharingDevice = new int[1];
    pesSharingDevice[0] = CkMyPe();
    numPesSharingDevice = 1;
  }

  if ( devicePe != CkMyPe() ) {
    CkPrintf("Pe %d physical rank %d will use MIC device of pe %d\n",
             CkMyPe(), myRankInPhysicalNode, devicePe);
    myDevice = -1;
    return;
  }

  // disable token-passing but don't submit local until remote finished
  // if shared_gpu is true, otherwise submit all work immediately
  first_pe_sharing_gpu = CkMyPe();
  next_pe_sharing_gpu = CkMyPe();

  gpu_is_mine = ( first_pe_sharing_gpu == CkMyPe() ); 

  if ( dev >= deviceCount ) {
    char buf[256];
    sprintf(buf,"Pe %d unable to bind to MIC device %d on %s because only %d devices are present",
		CkMyPe(), dev, host, deviceCount);
    NAMD_die(buf);
  }

  CkPrintf("Pe %d physical rank %d binding to MIC device %d on %s: %d procs %d threads \n",
             CkMyPe(), myRankInPhysicalNode, dev, host,
             omp_get_num_procs_target(TARGET_MIC,dev),
             omp_get_max_threads_target(TARGET_MIC,dev));

  myDevice = dev;

  // Initialize the MIC device (note, this should only be called once per device)
  mic_init_device(CkMyPe(), CkMyNode(), myDevice);

  int dev_ok = mic_check(dev);
  if ( dev_ok ) devices[ndevices++] = dev;
  if ( ! dev_ok ) {
    char errmsg[1024];
    sprintf(errmsg,"Failed binding to device %d on Pe %d", dev, CkMyPe());
    NAMD_die(errmsg);
  }
}

void mic_free() {
  if (devicePe != CkMyPe()) { return; }
  mic_free_device(myDevice);
}

static __thread ComputeNonbondedMIC* micCompute = 0;
static __thread ComputeMgr *computeMgr = 0;

void send_build_mic_force_table() {
  computeMgr->sendBuildMICForceTable();
}

void build_mic_force_table() {
  if (devicePe != CkMyPe()) { return; }
  ComputeNonbondedMIC::bind_lj_table(myDevice);
  ComputeNonbondedMIC::bind_force_table(myDevice);
  ComputeNonbondedMIC::bind_constants(myDevice);
}

void ComputeNonbondedMIC::bind_lj_table(int deviceNum) {
  if (devicePe != CkMyPe()) { return; }
  const int lj_table_dim = ljTable->get_table_dim();
  const int lj_table_size = 2 * lj_table_dim * lj_table_dim * sizeof(LJTable::TableEntry);
  const LJTable::TableEntry *lj_table_base_ptr = ljTable->get_table();
  mic_bind_lj_table(deviceNum, (char*)lj_table_base_ptr, lj_table_dim, lj_table_size);
}

void ComputeNonbondedMIC::bind_force_table(int deviceNum) {
  if (devicePe != CkMyPe()) { return; }
  mic_bind_table_four(deviceNum, mic_table_base_ptr, mic_table_n, mic_table_n_16);
}

void ComputeNonbondedMIC::bind_constants(int deviceNum) {
  if (devicePe != CkMyPe()) { return; }
  mic_bind_constants(deviceNum,
                     ComputeNonbondedUtil::cutoff2,
                     ComputeNonbondedUtil::dielectric_1,
                     ComputeNonbondedUtil::scaling,
                     ComputeNonbondedUtil::scale14,
                     ComputeNonbondedUtil::r2_delta,
                     ComputeNonbondedUtil::r2_delta_exp,
                     ComputeNonbondedUtil::commOnly
                    );
}

struct exlist_sortop {
  bool operator() (int32 *li, int32 *lj) {
    return ( li[1] < lj[1] );
  }
};

static __thread int2 *exclusionsByAtom = NULL;

void ComputeNonbondedMIC::bind_exclusions(int deviceNum) {

  // Only do this on PEs that directly control devices (simply return for PEs that do not)
  if (devicePe != CkMyPe()) { return; }

  // DMK - Info...
  //   Each atom as a globally unique index in the molecule object.  For any
  //   given atom, the exclusion list in the molecule object contains a list of
  //   the other atoms (referenced by index) that are "excluded" for the the
  //   given atom.  Since these exclusion lists tend to follow a set of patterns
  //   (i.e. water molecules "look" the same).  For this reason, we will store
  //   patterns rather than lists, and associate each atom with the appropriate
  //   pattern (e.g. 'all oxygen atoms within water molecules will exclude the
  //   hydrogen atoms next to them in the global list of atoms,' rather than
  //   both '[water_oxygen_1 excludes water_hydrogen_1 and water_hydrogen_2] and
  //   [water_oxygen_2 excludes water_hydrogen_3 and water_hydrogen_4]'). The code
  //   below creates those patterns and associates atoms with the appropriate
  //   patterns.  By using these common patterns, the memory required is reduced.
  //   To do this, the patterns use distances between the indexes of two atoms
  //   that exclude one another, rather than absolute indexes of those same atoms.
  //   A pattern is an array of bits (2-bits per entry), that ranges from negative
  //   offsets to positive offsets.  For example, for a water molecule made up of
  //   one oxygen (O1 w/ index i) and two hydrogen (H1 w/ index i+1 and H2 w/
  //   index i+2) atoms, the pattern for the oxygen atom (O1) would have 5 entries
  //   (10 bits total) and would range from offset -2 to 2 (always -X to +X).  The
  //   pattern for the first hydrogen (H1) would have 3 entries (6 bits total)
  //   and would range from offset -1 (O1) to 1 (H1).  The value of each 2-bit
  //   entry indicate the type of interaction (NORMAL, MODIFIED, or EXCLUDED).

  // Grab the molecule object, which contains information about the
  //   particle system being simulated.
  Molecule *mol = Node::Object()->molecule;

  // Get the number of atoms
  #ifdef MEM_OPT_VERSION
    int natoms = mol->exclSigPoolSize;
  #else
    int natoms = mol->numAtoms;
  #endif

  // A data structure that will store the offsets and ranges of the patterns
  //   created by the code below (see the comments below) for each atom,
  //   indicating which pattern should be used for the atom and which other
  //   atoms need to be checked against the pattern for this atom.
  exclusionsByAtom = new int2[natoms];

  // Declare some tmp variables/buffers used in the code below.
  ObjectArena<int32> listArena;
  ResizeArray<int32*> unique_lists;
  int32 **listsByAtom = new int32*[natoms];
  SortableResizeArray<int32> curExclList;
  SortableResizeArray<int32> curModList;
  SortableResizeArray<int32> curList;

  // For each atom...
  for (int i = 0; i < natoms; i++) {

    // Clear the current lists
    // NOTE: The values in these lists are offsets between the atom indexes
    curExclList.resize(0);  // EXCLUDED (full) list
    curModList.resize(0);   // MODIFIED list
    curList.resize(0);      // Combined list

    // Always exclude self (index offset to self is zero)
    curExclList.add(0);

    // Get the exclusions for this atom from the molecule object, adding each
    //   one to the current lists, converting absolut indexes to offsets as
    //   required.
    // NOTE: NAMD uses a macro called MEM_OPT_VERSION to adjust some of its data
    //   structures to be more memory efficient.  Whether or not MEM_OPT_VERSION
    //   is set changes how the exclusion information for atoms is stored in the
    //   molecule object.
    int n = 0;
    int min_index = 0;  // Furthest negative index offset
    int max_index = 0;  // Furthest positive index offset
    #if MEM_OPT_VERSION
      const ExclusionSignature *sig = mol->exclSigPool + i;
      n = sig->fullExclCnt;
      for (int j = 0; j < n; j++) {
        int index = sig->fullOffset[j];  // NOTE: Already an index offset
        if (index < min_index) { min_index = index; }
        if (index > max_index) { max_index = index; }
        curExclList.add(index);
      }
      n = sig->modExclCnt;
      for (int j = 0; j < n; j++) {
        int index = sig->modOffset[j];  // NOTE: Already an index offset
        if (index < min_index) { min_index = index; }
        if (index > max_index) { max_index = index; }
        curModList.add(index);
      }
    #else
      const int32 *full_excl = mol->get_full_exclusions_for_atom(i);
      n = full_excl[0] + 1;
      for (int j = 1; j < n; j++) {
        int index = full_excl[j] - i;  // NOTE: An absolute index, so make offset
        if (index < min_index) { min_index = index; }
        if (index > max_index) { max_index = index; }
        curExclList.add(index);
      }
      const int32 *mod_excl = mol->get_mod_exclusions_for_atom(i);
      n = mod_excl[0] + 1;
      for (int j = 1; j < n; j++) {
        int index = mod_excl[j] - i;
        if (index < min_index) { min_index = index; }
        if (index > max_index) { max_index = index; }
        curModList.add(index);
      }
    #endif
    int maxDiff = -1 * min_index;
    if (maxDiff < max_index) { maxDiff = max_index; }
    // NOTE: maxDiff = max(abs(min_index), max_index);
    curExclList.sort(); curExclList.add(-1 * maxDiff - 1);
    curModList.sort(); curModList.add(-1 * maxDiff - 1);
    // NOTE : curExclList and curModList now contain the list of offsets for
    //   atoms that are full and modified excluded (respectively) for this atom
    //   (i) in sorted order, with a '-1 * maxDiff - 1' value at the end (this
    //   last value is added for the sake of the "matching code" in the loop below).

    // Create a single list with the combined exclusion info of the two separate lists
    n = 2 * maxDiff + 1;  // NOTE: pattern ranges from -maxDiff to +maxDiff, inclusive
    curList.resize(n);
    int curExclListIndex = 0;
    int curModListIndex = 0;

    // For each entry in the combined list...
    for (int j = -1 * maxDiff; j <= maxDiff; j++) {

      // Assume normal
      int bitPattern = 0x00;

      // NOTE: curExclList and curModList are in sorted order and we are moving
      //   through the combined list in-order, so we only need to check one
      //   index in the separate lists (moving ahead in the list when a match
      //   is found).

      // Check if is actually excluded...
      if (curExclList[curExclListIndex] == j) {
        bitPattern |= 0x01;
        curExclListIndex++;
      }

      // Check if is actually modified...
      if (curModList[curModListIndex] == j) {
        bitPattern |= 0x02;
        curModListIndex++;
      }

      // Store the generated pattern entry
      curList[j + maxDiff] = bitPattern;
    }

    // NOTE: curList now contains the bit patterns (isModified and isExcluded)
    //   flags for this atom (i) based on the offsets of the other atoms.

    // Search through the list of patterns that have already been created, checking
    //   to see the current pattern already exists.  If so, just use the existing
    //   pattern.  If not, add this new, unique pattern.
    // For each existing pattern...
    int j = 0;
    for (j = 0; j < unique_lists.size(); j++) {

      // Check to see if the size matches (skip if it does not, obviously not a match)
      if (n != unique_lists[j][0]) { continue; }

      // Check each entry in the current pattern vs this existing pattern, breaking
      //   from the loop if a mismatch is found (i.e. causing the final value of k
      //   to be less than the pattern length).
      int k = 0;
      for (; k < n; k++) { // For each element in the list...
        if (unique_lists[j][k+3] != curList[k]) { break; } // ... check for mismatch
      }
      if (k == n) { break; } // If full match, stop searching (NOTE: j < unique_lists.size() will be true)
    }

    // If no match was found (i.e. j == length of existing pattern list) in the search
    //   loop above, then add the current pattern to the end of the list of patterns.
    if (j == unique_lists.size()) {

      // Create a new list with size n
      int32 *list = listArena.getNewArray(n + 3);
      list[0] = n;        // Entry [0] = length of pattern (not including list[0] or list[1])
      list[1] = maxDiff;  // Entry [1] = offset range
      // NOTE: Entry [2] will be filled in later

      // Entries [3+] contain the bit pattern entries
      for (int k = 0; k < n; k++) { list[k + 3] = curList[k]; }

      // Add this pattern to the list of patterns
      unique_lists.add(list);  // NOTE: This adds the new list at unique_lists[j]
    }

    // Note the pattern for this atom (whether it was found or added to the list of patterns)
    listsByAtom[i] = unique_lists[j];

  } // end for (i < natoms)

  // Sort the list of patterns, placing the smaller patterns first
  std::stable_sort(unique_lists.begin(), unique_lists.end(), exlist_sortop());

  // Now, we create a data structure that simply stores the patterns, using 2 bits
  //   per entry per atom.

  // Count the total bits required to store the lists (2 bits per entry in each
  //   pattern).  At the same time, note the offsets for the "self entry" for each
  //   pattern in the final list of bits for each pattern.  We want to use the
  //   "self entry" as the "start" so that we can easily create on offset in the
  //   pattern by subtracting the indexes of the two atoms and multiplying that
  //   by 2 bits during simulation.
  long int totalBits = 0;
  int nlists = unique_lists.size();
  // For each pattern in the list of patterns...
  for (int j = 0; j < nlists; j++) {
    int32 *list = unique_lists[j];  // Get the list
    int maxDiff = list[1];  // Get the range for the pattern
    list[2] = totalBits + (2 * maxDiff);  // Note the offset for the "self entry"/"start" for this pattern
    totalBits += 2 * (2 * maxDiff + 1); // Add the total bits needed for this pattern
  }

  // For each atom, record the "start" of the list (i.e. the "self" position)
  for (int i = 0; i < natoms; i++) {
    exclusionsByAtom[i].x = (listsByAtom[i])[1]; // maxDiff or range of pattern
    exclusionsByAtom[i].y = (listsByAtom[i])[2]; // "start" (self offset in bits)
  }

  // Cleanup listsByAtom (no longer required)
  delete [] listsByAtom; listsByAtom = NULL;

  // Roundup the total number of bits to a multiple of sizeof(unsigned int)
  const long int uintBitCnt = sizeof(unsigned int) * 8;
  if (totalBits & (uintBitCnt - 1)) {
    totalBits += (uintBitCnt - (totalBits & (uintBitCnt - 1)));
  }
  long int totalBytes = totalBits / 8;

  // If this is PE 0, print some info...
  if (!CkMyPe()) {
    CkPrintf("Info: Found %d unique exclusion lists needed %ld bytes\n",
             unique_lists.size(), totalBytes);
  }

  // Allocate the memory required (using 'unsigned int' as the array's data type)
  //   and initialize the memory to zero.
  long int totalUInts = totalBits / (sizeof(unsigned int) * 8);
  unsigned int *exclusion_bits = new unsigned int[totalUInts];
  memset(exclusion_bits, 0, totalBytes);  // Zero-out the data

  // Fill in the bits
  long int offset = 0;  // Offset of current list in exclusion_bits array
  // For each of the patterns...
  for (int i = 0; i < unique_lists.size(); i++) {

    // Get the pattern
    int32 *list = unique_lists[i];

    // Sanity check: Verify that the base value matches with the stored "start"
    //   values (i.e. the start of this list is where we expect it to be)
    if (offset + (2 * list[1]) != list[2]) {
      NAMD_bug("ComputeNonbondedMIC::bind_exclusions offset mismatch");
    }

    // Pack the bits from this list into the exclusion_bits array
    const int n = list[0];
    for (int j = 0; j < n; j++) {
      const int offset_j = offset + (2 * j);  // NOTE: Each entry is 2 bits
      const int offset_major = offset_j / (sizeof(unsigned int) * 8);
      const int offset_minor = offset_j % (sizeof(unsigned int) * 8); // NOTE: Reverse indexing direction relative to offset_major
      const int entry_mask = (list[j + 3] & 0x3) << offset_minor;
      exclusion_bits[offset_major] |= entry_mask;
    }

    // Move the offset forward
    offset += 2 * (2 * list[1] + 1);
  }

  // Now that the bits for all of the patterns have been placed into the
  //   exclusion_bits array, push this array down to the device for use
  //   during simulation.
  mic_bind_exclusions(deviceNum, exclusion_bits, totalUInts);

  // Cleanup
  // NOTE: That the individual lists will be free'ed as the arena is destroyed,
  //   along with the resize arrays (curModList, etc.).
  delete [] exclusion_bits;
}


// Register a "self compute" on the host with this MIC meta-compute object, creating a
//   compute_record for the self compute.
void register_mic_compute_self(ComputeID c, PatchID pid, int part, int numParts) {

  if ( ! micCompute ) NAMD_bug("register_self called early");

  //// DMK - DEBUG
  //printf("register_mic_compute_self(c:%d, pid:%d, part:%d, numParts:%d)...\n",
  //       c, pid, part, numParts
  //      );

  // Indicate that the mic compute requires the patch information
  //   associated with the given self compute
  micCompute->requirePatch(pid);

  #if MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING != 0
    numParts = MIC_SPECIFIC_COMPUTE_PARTITIONING__SELF_P1;
    for (int part = 0; part < numParts; part++) {
  #endif

    // Create a compute record within the mic compute that represents
    //   the given self compute
    ComputeNonbondedMIC::compute_record cr;
    cr.c = c;
    cr.pid[0] = pid;  cr.pid[1] = pid;
    cr.offset = 0.;
    cr.isSelf = 1;
 
    #if (MIC_ENABLE_COMPUTE_PARTITIONING != 0) || (MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING != 0)
      cr.part = part;
      cr.numParts = numParts;
    #endif 

    #if MIC_MULTIPLE_KERNELS != 0
      if (micCompute->patchRecords[pid].isLocal) {
    #endif
        micCompute->localComputeRecords.add(cr);
    #if MIC_MULTIPLE_KERNELS != 0
      } else {
        micCompute->remoteComputeRecords.add(cr);
      }
    #endif

  #if MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING != 0
    }
  #endif

  //// DMK - NOTE - We are only issuing one offload (local) at the moment, so
  ////   make sure all computes are "local"
  //if (1) {
  ////if ( micCompute->patchRecords[pid].isLocal ) {
  //  micCompute->localComputeRecords.add(cr);
  //} else {
  //  micCompute->remoteComputeRecords.add(cr);
  //}
}

void register_mic_compute_pair(ComputeID c, PatchID pid[], int t[], int part, int numParts) {

  if ( ! micCompute ) NAMD_bug("register_pair called early");

  //// DMK - DEBUG
  //printf("register_mic_compute_self(c:%d, pid:{%d,%d}, part:%d, numParts:%d, ...)...\n",
  //       c, pid[0], pid[1], part, numParts
  //      );

  // Indicate that the mic compute requires the patch information
  //   associated with the given pair compute
  micCompute->requirePatch(pid[0]);
  micCompute->requirePatch(pid[1]);

  // Calculate the offset that will need to be applied to the atom positions (which are
  //   stored as offsets from the center of their respective patches)
  int t1 = t[0];
  int t2 = t[1];
  Vector offset = micCompute->patchMap->center(pid[0])
                - micCompute->patchMap->center(pid[1]);
  offset.x += (t1%3-1) - (t2%3-1);
  offset.y += ((t1/3)%3-1) - ((t2/3)%3-1);
  offset.z += (t1/9-1) - (t2/9-1);

  #if MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING != 0

    // Calculate the Manhattan distance between the patches and use that to determine how
    //   many parts the pair compute should be broken up into
    ComputeMap *computeMap = ComputeMap::Object();
    PatchMap *patchMap = PatchMap::Object();
    int aSize = patchMap->gridsize_a();
    int bSize = patchMap->gridsize_b();
    int cSize = patchMap->gridsize_c();
    int trans0 = computeMap->trans(c, 0);
    int trans1 = computeMap->trans(c, 1);
    int index_a0 = patchMap->index_a(pid[0]) + aSize * Lattice::offset_a(trans0);
    int index_b0 = patchMap->index_b(pid[0]) + bSize * Lattice::offset_b(trans0);
    int index_c0 = patchMap->index_c(pid[0]) + cSize * Lattice::offset_c(trans0);
    int index_a1 = patchMap->index_a(pid[1]) + aSize * Lattice::offset_a(trans1);
    int index_b1 = patchMap->index_b(pid[1]) + bSize * Lattice::offset_b(trans1);
    int index_c1 = patchMap->index_c(pid[1]) + cSize * Lattice::offset_c(trans1);
    int da = index_a0 - index_a1; da *= ((da < 0) ? (-1) : (1));
    int db = index_b0 - index_b1; db *= ((db < 0) ? (-1) : (1));
    int dc = index_c0 - index_c1; dc *= ((dc < 0) ? (-1) : (1));
    int manDist = da + db + dc;

    // Create each part
    //numParts = 8 >> manDist;
    numParts = MIC_SPECIFIC_COMPUTE_PARTITIONING__PAIR_P1 - (MIC_SPECIFIC_COMPUTE_PARTITIONING__PAIR_P2 * manDist);
    if (numParts <= 0) { numParts = 1; }
    for (int part = 0; part < numParts; part++) {

  #endif

    //// DMK - DEBUG
    //printf("[DEBUG] creating pair compute %d,%d\n", c, part);

    // Create a compute record within the mic compute that represents
    //   the given pair compute
    // DMK : NOTE - We are going to use equal and opposite on MIC, so
    //   modify this code to only create a single compute record
    ComputeNonbondedMIC::compute_record cr;
    cr.c = c;
    cr.pid[0] = pid[0];  cr.pid[1] = pid[1];

    cr.offset = offset;

    #if (MIC_ENABLE_COMPUTE_PARTITIONING != 0) || (MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING != 0)
      cr.part = part;
      cr.numParts = numParts;
    #endif

    cr.isSelf = 0;

    #if MIC_MULTIPLE_KERNELS != 0
      if (micCompute->patchRecords[pid[0]].isLocal &&
          micCompute->patchRecords[pid[1]].isLocal
         ) {
    #endif
        micCompute->localComputeRecords.add(cr);
    #if MIC_MULTIPLE_KERNELS != 0
      } else {
        micCompute->remoteComputeRecords.add(cr);
      }
    #endif

  #if MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING != 0
    }
  #endif

  //// DMK - NOTE : Changing the logic here a bit.  Since we are using
  ////   equal and opposite, the concept of a compute record being only
  ////   local or remote is more restricted (i.e. a much higher fraction
  ////   of pairs will be remote, since it would be remote if either
  ////   patch is remote and since we are contributing to both patches
  ////   instead of just one of them).  So, for MIC, call local only
  ////   if both patches are local.  Note that this will probably tend
  ////   to make the remote group larger.
  //// NOTE | TODO | FIXME : We are only issuing one kernel for now
  ////   so mark all as "local"
  //if (1) {
  //  //  micCompute->patchRecords[pid[0]].isLocal &&
  //  //  micCompute->patchRecords[pid[1]].isLocal
  //  //) {
  //  micCompute->localComputeRecords.add(cr);
  //} else {
  //  micCompute->remoteComputeRecords.add(cr);
  //}
}

void unregister_mic_compute(ComputeID c) {  // static
  NAMD_bug("unregister_compute unimplemented");
}

static __thread int atomsChanged = 0;
static __thread int computesChanged = 0;

static __thread int pairlistsValid = 0;
static __thread float pairlistTolerance = 0.;
static __thread int usePairlists = 0;
static __thread int savePairlists = 0;
static __thread float plcutoff2 = 0;

static __thread ResizeArray<patch_pair> *patch_pairs_ptr;
static __thread ResizeArray<force_list> *force_lists_ptr;

#define PATCH_PAIRS_REF ResizeArray<patch_pair> &patch_pairs(*patch_pairs_ptr);
#define FORCE_LISTS_REF ResizeArray<force_list> &force_lists(*force_lists_ptr);

ComputeNonbondedMIC::ComputeNonbondedMIC(ComputeID c,
                                         ComputeMgr *mgr,
                                         ComputeNonbondedMIC *m, int idx
                                        ) : Compute(c), slaveIndex(idx) {

  #ifdef PRINT_GBIS
    CkPrintf("C.N.MIC[%d]::constructor cid=%d\n", CkMyPe(), c);
  #endif

  //// DMK - DEBUG
  //printf("[%d] ComputeNonbondedMIC::ComputeNonbondedMIC(cid:%d) - Called...\n", CkMyPe(), c);
  //fflush(NULL);

  master = m ? m : this;
  micCompute = this;
  computeMgr = mgr;
  patchMap = PatchMap::Object();
  atomMap = AtomMap::Object();
  reduction = 0;

  SimParameters *params = Node::Object()->simParameters;
  if (params->pressureProfileOn) {
    NAMD_die("pressure profile not supported in MIC");
  }

  atomsChanged = 1;
  computesChanged = 1;
  workStarted = 0;
  basePriority = PROXY_DATA_PRIORITY;
  localWorkMsg2 = new (PRIORITY_SIZE) LocalWorkMsg;

  if ( master != this ) { // I am slave
    masterPe = master->masterPe;
    master->slaves[slaveIndex] = this;
    if ( master->slavePes[slaveIndex] != CkMyPe() ) {
      NAMD_bug("ComputeNonbondedMIC slavePes[slaveIndex] != CkMyPe");
    }
    registerPatches();
    return;
  } else {  // I am a master, identify self to ComputeMgr for load balancing data
    computeMgr->sendMICPEData(CkMyPe(), 1);
  }
  masterPe = CkMyPe();

  if (myDevice >= 0) {
    bind_exclusions(myDevice);
  }

  patch_pairs_ptr = new ResizeArray<patch_pair>;
  force_lists_ptr = new ResizeArray<force_list>;

  #if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0
    micDevice = myDevice;
    exclusionsByAtom_ptr = exclusionsByAtom;
  #endif

  // DMK - DEBUG
  timestep = 0;
}


ComputeNonbondedMIC::~ComputeNonbondedMIC() { ; }

void ComputeNonbondedMIC::requirePatch(int pid) {

  computesChanged = 1;
  patch_record &pr = patchRecords.item(pid);
  if ( pr.refCount == 0 ) {
    //if ( mergegrids ) {
    //  pr.isLocal = 0;
    //} else if ( CkNumNodes() < 2 ) {
    //  pr.isLocal = 1; //1 & ( 1 ^ patchMap->index_a(pid) ^ patchMap->index_b(pid) ^ patchMap->index_c(pid) );
    //} else {
    //  pr.isLocal = 1; //( CkNodeOf(patchMap->node(pid)) == CkMyNode() );
    //}
    pr.isLocal = ( CkNodeOf(patchMap->node(pid)) == CkMyNode() ); // && ( pid > 30 );

    #if MIC_MULTIPLE_KERNELS != 0
      if ( pr.isLocal ) {
    #endif
        localActivePatches.add(pid);
    #if MIC_MULTIPLE_KERNELS != 0
      } else {
        remoteActivePatches.add(pid);
      }
    #endif

    activePatches.add(pid);
    pr.patchID = pid;
    pr.hostPe = -1;
    pr.x = NULL;
    pr.xExt = NULL;
    pr.r = NULL;
    pr.f = NULL;
    pr.intRad      = NULL;
    pr.psiSum      = NULL;
    pr.bornRad     = NULL;
    pr.dEdaSum     = NULL;
    pr.dHdrPrefix  = NULL;
  }
  pr.refCount += 1;
}

void ComputeNonbondedMIC::registerPatches() {

  SimParameters *simParams = Node::Object()->simParameters;
  int npatches = master->activePatches.size();
  int *pids = master->activePatches.begin();
  patch_record *recs = master->patchRecords.begin();
  const int myPE = CkMyPe();

  for ( int i=0; i<npatches; ++i ) {
    int pid = pids[i];
    patch_record &pr = recs[pid];
    if ( pr.hostPe == myPE ) {
      //if ( ! reduction ) {
      //  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
      //}
      hostedPatches.add(pid);

      #if MIC_MULTIPLE_KERNELS != 0
        if ( pr.isLocal ) {
      #endif
          localHostedPatches.add(pid);
      #if MIC_MULTIPLE_KERNELS != 0
        } else {
          remoteHostedPatches.add(pid);
        }
      #endif

      ProxyMgr::Object()->createProxy(pid);
      pr.p = patchMap->patch(pid);
      pr.positionBox = pr.p->registerPositionPickup(this);
      pr.forceBox = pr.p->registerForceDeposit(this);
      if (simParams->GBISOn) {
        pr.intRadBox      = pr.p->registerIntRadPickup(this);
        pr.psiSumBox      = pr.p->registerPsiSumDeposit(this);
        pr.bornRadBox     = pr.p->registerBornRadPickup(this);
        pr.dEdaSumBox     = pr.p->registerDEdaSumDeposit(this);
        pr.dHdrPrefixBox  = pr.p->registerDHdrPrefixPickup(this);
      }
    }
  }
  if ( master == this ) setNumPatches(activePatches.size());
  else setNumPatches(hostedPatches.size());
  CkPrintf("Pe %d hosts %d local and %d remote patches for pe %d\n", CkMyPe(), localHostedPatches.size(), remoteHostedPatches.size(), masterPe);
}

void ComputeNonbondedMIC::assignPatches() {

  int *pesOnNodeSharingDevice = new int[CkMyNodeSize()];
  int numPesOnNodeSharingDevice = 0;
  int masterIndex = -1;
  for ( int i=0; i<numPesSharingDevice; ++i ) {
    int pe = pesSharingDevice[i];
    if ( pe == CkMyPe() ) masterIndex = numPesOnNodeSharingDevice;
    if ( CkNodeOf(pe) == CkMyNode() ) {
      pesOnNodeSharingDevice[numPesOnNodeSharingDevice++] = pe;
    }
  }

  int npatches = activePatches.size();

  if ( npatches ) {
    reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  }

  int *count = new int[npatches];
  memset(count, 0, sizeof(int)*npatches);
  int *pcount = new int[numPesOnNodeSharingDevice];
  memset(pcount, 0, sizeof(int)*numPesOnNodeSharingDevice);
  int *rankpcount = new int[CkMyNodeSize()];
  memset(rankpcount, 0, sizeof(int)*CkMyNodeSize());
  char *table = new char[npatches*numPesOnNodeSharingDevice];
  memset(table, 0, npatches*numPesOnNodeSharingDevice);

  int unassignedpatches = npatches;

  if ( 0 ) { // assign all to device pe
    for ( int i=0; i<npatches; ++i ) {
      int pid = activePatches[i];
      patch_record &pr = patchRecords[pid];
      pr.hostPe = CkMyPe();
    }
    unassignedpatches = 0;
    pcount[masterIndex] = npatches;
  } else 

  // assign if home pe and build table of natural proxies
  for ( int i=0; i<npatches; ++i ) {
    int pid = activePatches[i];
    patch_record &pr = patchRecords[pid];
    int homePe = patchMap->node(pid);

    for ( int j=0; j<numPesOnNodeSharingDevice; ++j ) {
      int pe = pesOnNodeSharingDevice[j];
      if ( pe == homePe ) {
        pr.hostPe = pe;  --unassignedpatches;
        pcount[j] += 1;
      }
      if ( PatchMap::ObjectOnPe(pe)->patch(pid) ) {
        table[i*numPesOnNodeSharingDevice+j] = 1;
      }
    }
    if ( pr.hostPe == -1 && CkNodeOf(homePe) == CkMyNode() ) {
      pr.hostPe = homePe;  --unassignedpatches;
      rankpcount[CkRankOf(homePe)] += 1;
    }
  }

  // assign if only one pe has a required proxy
  int assignj = 0;
  for ( int i=0; i<npatches; ++i ) {
    int pid = activePatches[i];
    patch_record &pr = patchRecords[pid];
    if ( pr.hostPe != -1 ) continue;
    int c = 0;
    int lastj;
    for ( int j=0; j<numPesOnNodeSharingDevice; ++j ) {
      if ( table[i*numPesOnNodeSharingDevice+j] ) { ++c; lastj=j; }
    }
    count[i] = c;
    if ( c == 1 ) {
      pr.hostPe = pesOnNodeSharingDevice[lastj];
      --unassignedpatches;
      pcount[lastj] += 1;
    }
  }
  while ( unassignedpatches ) {
    int i;
    for ( i=0; i<npatches; ++i ) {
      if ( ! table[i*numPesOnNodeSharingDevice+assignj] ) continue;
      int pid = activePatches[i];
      patch_record &pr = patchRecords[pid];
      if ( pr.hostPe != -1 ) continue;
      pr.hostPe = pesOnNodeSharingDevice[assignj];
      --unassignedpatches;
      pcount[assignj] += 1;
      if ( ++assignj == numPesOnNodeSharingDevice ) assignj = 0;
      break;
    }
    if ( i<npatches ) continue;  // start search again
    for ( i=0; i<npatches; ++i ) {
      int pid = activePatches[i];
      patch_record &pr = patchRecords[pid];
      if ( pr.hostPe != -1 ) continue;
      if ( count[i] ) continue;
      pr.hostPe = pesOnNodeSharingDevice[assignj];
      --unassignedpatches;
      pcount[assignj] += 1;
      if ( ++assignj == numPesOnNodeSharingDevice ) assignj = 0;
      break;
    }
    if ( i<npatches ) continue;  // start search again
    if ( ++assignj == numPesOnNodeSharingDevice ) assignj = 0;
  }

  for ( int i=0; i<npatches; ++i ) {
    int pid = activePatches[i];
    patch_record &pr = patchRecords[pid];
    // CkPrintf("Pe %d patch %d hostPe %d\n", CkMyPe(), pid, pr.hostPe);
  }

  slavePes = new int[CkMyNodeSize()];
  slaves = new ComputeNonbondedMIC*[CkMyNodeSize()];
  numSlaves = 0;
  for ( int j=0; j<numPesOnNodeSharingDevice; ++j ) {
    int pe = pesOnNodeSharingDevice[j];
    int rank = pe - CkNodeFirst(CkMyNode());
    // CkPrintf("host %d sharing %d pe %d rank %d pcount %d rankpcount %d\n",
    //          CkMyPe(),j,pe,rank,pcount[j],rankpcount[rank]);
    if ( pe == CkMyPe() ) continue;
    if ( ! pcount[j] && ! rankpcount[rank] ) continue;
    rankpcount[rank] = 0;  // skip in rank loop below
    slavePes[numSlaves] = pe;
    computeMgr->sendCreateNonbondedMICSlave(pe,numSlaves);
    ++numSlaves;
  }
  for ( int j=0; j<CkMyNodeSize(); ++j ) {
    int pe = CkNodeFirst(CkMyNode()) + j;
    // CkPrintf("host %d rank %d pe %d rankpcount %d\n",
    //          CkMyPe(),j,pe,rankpcount[j]);
    if ( ! rankpcount[j] ) continue;
    if ( pe == CkMyPe() ) continue;
    slavePes[numSlaves] = pe;
    computeMgr->sendCreateNonbondedMICSlave(pe,numSlaves);
    ++numSlaves;
  }
  registerPatches();

  delete [] pesOnNodeSharingDevice;
  delete [] count;
  delete [] pcount;
  delete [] rankpcount;
  delete [] table;
}

static __thread int num_atom_records_allocated;
// static __thread int* atom_order;
static __thread atom_param* atom_params;
static __thread atom* atoms;
// static __thread float4* forces;
// static __thread float4* slow_forces;
//static __thread int num_virials;
//static __thread int num_virials_allocated;
//static __thread float *virials;
//static __thread float *slow_virials;
static __thread float *energy_gbis;

//GBIS host pointers
static __thread float *intRad0H;
static __thread float *intRadSH;
//static __thread GBReal *psiSumH; //moved into class
static __thread float *bornRadH;
//static __thread GBReal *dEdaSumH; //moved into class
static __thread float *dHdrPrefixH;

static __thread int mic_timer_count;
static __thread double mic_timer_total;
static __thread double kernel_time;
static __thread double remote_submit_time;
static __thread double local_submit_time;

// TRACING
#if MIC_TRACING != 0
  __thread double mic_tracing_offload_start_remote = 0.0;
  __thread double mic_tracing_offload_start_local = 0.0;
  static __thread double mic_tracing_polling_start = 0.0;
  static __thread int mic_tracing_polling_count = 0;
  #define MIC_TRACING_POLLING_SET_FREQ  ( 100 )
#endif

#define MIC_POLL(FN,ARG) CcdCallFnAfter(FN,ARG,0.1)

#ifdef PRINT_GBIS
#define GBISP(...) CkPrintf(__VA_ARGS__);
#else
#define GBISP(...)
#endif


#define count_limit 5000000  // DMK - NOTE : I have this set fairly high so that I can test single-thread executions on the MIC device, which can take quite some time if using printfs to display intermediate values.
static __thread int check_remote_count;
static __thread int check_local_count;

void mic_check_remote_progress(void *arg, double) {

  // If the offloaded work is finished, trigger finishWork()
  if (mic_check_remote_kernel_complete(myDevice)) {

    ((ComputeNonbondedMIC*)arg)->messageFinishWork();
    #if MIC_TRACING != 0
      double _now = CmiWallTimer();
      MIC_TRACING_RECORD(MIC_EVENT_OFFLOAD_REMOTE, mic_tracing_offload_start_remote, _now);
      MIC_TRACING_RECORD(MIC_EVENT_OFFLOAD_POLLSET, mic_tracing_polling_start, _now);
      mic_tracing_polling_count = 0;
    #endif
    check_remote_count = 0;

    //// DMK - DEBUG
    //printf("[DEBUG:%d] :: << detected remote kernel completion >>\n", CkMyPe()); fflush(NULL);

  // Otherwise, check to see if it has been a long time (if so, timeout with error)
  } else if (++check_remote_count >= count_limit) {
    char errmsg[256];
    sprintf(errmsg,
            "mic_check_remote_progress polled %d times over %f s on step %d",
            check_remote_count, CkWallTimer() - remote_submit_time,
            ((ComputeNonbondedMIC*)arg)->step
           );
    //mic_errcheck(errmsg);
    NAMD_die(errmsg);

  // Otherwise, queue another poll attempt in the future to try again later
  } else {

    #if MIC_TRACING != 0
      if (++mic_tracing_polling_count > MIC_TRACING_POLLING_SET_FREQ) {
        double now = CmiWallTimer();
        MIC_TRACING_RECORD(MIC_EVENT_OFFLOAD_POLLSET, mic_tracing_polling_start, now);
        mic_tracing_polling_start = now;
        mic_tracing_polling_count = 0;
      }
    #endif
    MIC_POLL(mic_check_remote_progress, arg);
  }
}

void mic_check_local_progress(void *arg, double) {

  // If the offloaded work is finished, trigger finishWork()
  if (mic_check_local_kernel_complete(myDevice)) {

    ((ComputeNonbondedMIC*)arg)->messageFinishWork();
    #if MIC_TRACING != 0
      double _now = CmiWallTimer();
      MIC_TRACING_RECORD(MIC_EVENT_OFFLOAD_LOCAL, mic_tracing_offload_start_local, _now);
      MIC_TRACING_RECORD(MIC_EVENT_OFFLOAD_POLLSET, mic_tracing_polling_start, _now);
      mic_tracing_polling_count = 0;
    #endif
    #if MIC_PRAGMA_TIMING_STATS != 0
      #if MIC_TRACING != 0
      #else
        double _now = CmiWallTimer();
      #endif
	pragma_timing_submitFinish(_now, CkMyPe(), CkMyRank());
    #endif
    check_local_count = 0;

    //// DMK - DEBUG
    //printf("[DEBUG:%d] :: << detected local kernel completion >>\n", CkMyPe()); fflush(NULL);

  // Otherwise, check to see if it has been a long time (if so, timeout with error)
  } else if (++check_local_count >= count_limit) {
    char errmsg[256];
    sprintf(errmsg,
            "mic_check_local_progress polled %d times over %f s on step %d",
            check_local_count, CkWallTimer() - local_submit_time,
            ((ComputeNonbondedMIC*)arg)->step
           );
    //mic_errcheck(errmsg);
    NAMD_die(errmsg);

  // Shouldn't be polling for local complete until remote complete was already detected (check for this case)
  } else if ( check_remote_count ) {
    NAMD_bug("nonzero check_remote_count in mic_check_local_progres");

  // Otherwise, queue another poll attmpt in the future to try again later
  } else {

    //#if MIC_TRACING != 0
    //  if (++mic_tracing_polling_count > MIC_TRACING_POLLING_SET_FREQ) {
    //    double now = CmiWallTimer();
    //    MIC_TRACING_RECORD(MIC_EVENT_OFFLOAD_POLLSET, mic_tracing_polling_start, now);
    //    mic_tracing_polling_start = now;
    //    mic_tracing_polling_count = 0;
    //  }
    //#endif
    MIC_POLL(mic_check_local_progress, arg);
  }
}

void ComputeNonbondedMIC::atomUpdate() { atomsChanged = 1; }

static __thread int kernel_launch_state = 0;

struct cr_sortop {
  const Lattice &l;
  cr_sortop(const Lattice &lattice) : l(lattice) { }
  bool operator() (ComputeNonbondedMIC::compute_record i,
                   ComputeNonbondedMIC::compute_record j) {
    Vector a = l.a();
    Vector b = l.b();
    Vector c = l.c();
    BigReal ri = (i.offset.x * a + i.offset.y * b + i.offset.z * c).length2();
    BigReal rj = (j.offset.x * a + j.offset.y * b + j.offset.z * c).length2();
    return ( ri < rj );
  }
};

//void ComputeNonbondedMIC::patchReady(PatchID patchID, int doneMigration, int seq) {
//  printf("[DMK-DEBUG] :: PE %d :: ComputeNonbondedMIC::patchReady(pID:%d, dM:%d, s:%d) - Called...\n",
//         CkMyPe(), patchID, doneMigration, seq
//        );
//  Compute::patchReady(patchID, doneMigration, seq);
//}

#if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0
void ComputeNonbondedMIC::patchReady(PatchID patchID, int doneMigration, int seq) {

  //// DMK - DEBUG
  //printf("[DEBUG] :: PE %d :: ComputeNonbondedMIC::patchReady(pID:%d, dM:%d, s:%d) - Called...\n", CkMyPe(), patchID, doneMigration, seq);
  //printf("[DEBUG] :: PE %d ::   &patchID:%p\n", CkMyPe(), &patchID);
  //fflush(NULL);

  // This function overrides the Compute::patchReady function so that some additional
  //   work can be done.  Compute::patchReady is called at the end of the function.
  //   When submitting atoms as the atom data arrives on the node, this function along
  //   with mic_submit_patch_data take care of pushing that data to the devices that
  //   require them.  However, this function is called by each of the host cores, so
  //   some care must be taken when doing this.  There are two items that need to be
  //   accomplished each timestep.  First, the atom data needs to be packed up for
  //   transfer to the devices (once per patch).  Second, for each device that requires
  //   the data, push that atom data to the device (once per patch per device).  A lock
  //   is used to protect flags that track what portions of this work have been done
  //   as each of the host threads are notified that the patches are ready.  The flags
  //   are mic_atomData_seq for once per patch work and mic_atomData_deviceSeq[] for
  //   once per patch per device work.

  // Get the patch information
  patch_record &pr = master->patchRecords[patchID];

  //// DMK - DEBUG
  //printf("[DEBUG] :: PE %d :: pr:%p, master:%p, master->patchRecords:%p\n", CkMyPe(), &pr, master, &(master->patchRecords));
  //fflush(NULL);

  Patch *p = pr.p;
  CudaAtom *ca = p->getCudaAtomList();
  CompAtom *a = pr.positionBox->open(); pr.x = a;
  CompAtomExt *aExt = p->getCompAtomExtInfo();
  int2 *exclusionsByAtom = master->exclusionsByAtom_ptr;
  int numAtoms = p->getNumAtoms();
  int numAtoms_16 = ((numAtoms + 15) & (~15));

  // Create references to the atomData pointers on the host for this patch (host copy of data on device)
  void* &atomData = p->mic_atomData;
  int &allocSize = p->mic_atomData_allocSize_host;
  // NOTE: Since we are using a mutex (per patch) to protect the following region,
  //   only the first thread through will reallocate the memory and copy the patch
  //   data in to the host's copy of the atomData buffer (seq number as check).
  //   The remaining threads will check the size and see that everything is fine,
  //   skipping this work.

  // Note the MIC device that this call wants the device on (the device that the master object drives)
  int micDevice = master->micDevice;
  __ASSERT(micDevice < MIC_MAX_DEVICES_PER_NODE);

  int transferFlag = 0;

  // If either the 'once per patch' or 'once per patch per device' work needs to be done, then grab
  //   the lock and dow the work
  if (p->mic_atomData_seq != seq || p->mic_atomData_deviceSeq[micDevice] != seq) {
    pthread_mutex_lock(&(p->mic_atomData_mutex));

    // Clear the flags early so other threads are more likely to skip by the checks for the work that
    //   will be done by this thread (and the lock).  If another thread does pass the check, the lock and
    //   flag values will still avoid the work being done multiple times.  This is just an optimization.
    int tmp_mic_atomData_seq = p->mic_atomData_seq;
    int tmp_mic_atomData_deviceSeq = p->mic_atomData_deviceSeq[micDevice];
    p->mic_atomData_seq = seq;
    p->mic_atomData_deviceSeq[micDevice] = seq;

    // Once per patch per timestep work
    if (tmp_mic_atomData_seq != seq) {  // Check again in case the first check passed while another thread was in the region

      // Allocate the memory as needed
      int allocFlag = 0;
      if (numAtoms_16 > allocSize || atomData == NULL) {
        if (atomData != NULL) { _mm_free(atomData); }
        int toAllocSize = (int)(numAtoms_16 * 1.1f);
        toAllocSize = ((toAllocSize + 15) & (~15));
        atomData = (char*)(_mm_malloc((sizeof(atom) + sizeof(atom_param)) * toAllocSize, MIC_ALIGN));
        __ASSERT(atomData != NULL);
        allocSize = toAllocSize;
        allocFlag = 1;
      }

      // Copy the data to the buffer that will be passed to the device(s)
      // NOTE: the number of atoms will only change when doneMigration is set
      // WARNING | NOTE : The following memcopy assumes CudaAtom and atom data structures match !!!
      atom* dev_a = (atom*)atomData;
      atom_param* dev_aExt = (atom_param*)(dev_a + numAtoms_16);
      memcpy(dev_a, a, sizeof(atom) * numAtoms);  // atoms always
      if (doneMigration || allocFlag) { // atom_params sometimes
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          for (int k = 0; k < numAtoms; k++) {

            int j = aExt[k].sortOrder;
            dev_aExt[k].vdw_type = a[j].vdwType;
            dev_aExt[k].index = aExt[j].id;

            //printf("[DEBUG] :: j:%d, k:%d, dev_aExt:%p, a:%p, aExt:%p, exclusionsByAtom:%p\n", j, k, dev_aExt, a, aExt, exclusionsByAtom);
            //fflush(NULL);
            //printf("[DEBUG] :: aExt[j:%d].id:%d\n", j, aExt[j].id);
            //fflush(NULL);

            #ifdef MEM_OPT_VERSION
              dev_aExt[k].excl_index = exclusionsByAtom[aExt[j].exclId].y;
              dev_aExt[k].excl_maxdiff = exclusionsByAtom[aExt[j].exclId].x;
            #else
              dev_aExt[k].excl_index = exclusionsByAtom[aExt[j].id].y;
              dev_aExt[k].excl_maxdiff = exclusionsByAtom[aExt[j].id].x;
            #endif
          }
        #else
          int *dev_aExt_vdwType     = ((int*)dev_aExt) + (0 * numAtoms_16);
          int *dev_aExt_index       = ((int*)dev_aExt) + (1 * numAtoms_16);
          int *dev_aExt_exclIndex   = ((int*)dev_aExt) + (2 * numAtoms_16);
          int *dev_aExt_exclMaxDiff = ((int*)dev_aExt) + (3 * numAtoms_16);
          for (int k = 0; k < numAtoms; k++) {
            int j = aExt[k].sortOrder;
            dev_aExt_vdwType[k] = a[j].vdwType;
            dev_aExt_index[k] = aExt[j].id;
            #ifdef MEM_OPT_VERSION
              dev_aExt_exclIndex[k] = exclusionsByAtom[aExt[j].exclId].y;
              dev_aExt_exclMaxDiff[k] = exclusionsByAtom[aExt[j].exclId].x;
            #else
              dev_aExt_exclIndex[k] = exclusionsByAtom[aExt[j].id].y;
              dev_aExt_exclMaxDiff[k] = exclusionsByAtom[aExt[j].id].x;
            #endif
          }
        #endif
      }
    } // end if (mic_atomData_seq != seq)

    // Once per patch per timestep per device work
    // NOTE : Within the protected region, simply flag that the transfer needs to take place
    //   and move on.  This will allow transfers to multiple MIC cards to occur in parallel
    //   without the per-patch-lock serializing them.
    if (tmp_mic_atomData_deviceSeq != seq) { transferFlag = 1; }

    pthread_mutex_unlock(&(p->mic_atomData_mutex));
  } // end if (mic_atomData_seq != seq || mic_atomData_deviceSeq[micDevice] != seq)

  // Transfer the data to the given device
  if (transferFlag != 0) {
    int allocBytes = allocSize * (sizeof(atom) + sizeof(atom_param));
    int transferBytes = numAtoms_16 * sizeof(atom);
    if (doneMigration) { transferBytes += numAtoms_16 * sizeof(atom_params); }
    mic_submit_patch_data(micDevice,
                          p->mic_atomData,
                          p->mic_atomData_prev[micDevice],
                          transferBytes,
                          allocBytes,
                          p->mic_atomData_allocSize_device[micDevice],
                          p->mic_atomData_devicePtr[micDevice]
                         );
  }

  //// DMK - DEBUG
  //printf("[DEBUG] :: PE %d :: ComputeNonbondedMIC::patchReady(...) - Finished.\n", CkMyPe());
  //fflush(NULL);

  // Call parent class patchReady function to perform the usual work
  Compute::patchReady(patchID, doneMigration, seq);
}
#endif // MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0


int ComputeNonbondedMIC::noWork() {

  //// DMK - DEBUG
  //printf("[DEBUG:%d] :: ComputeNonbondedMIC::noWork() - Called...\n", CkMyPe()); fflush(NULL);

  SimParameters *simParams = Node::Object()->simParameters;
  Flags &flags = master->patchRecords[hostedPatches[0]].p->flags;
  lattice = flags.lattice;
  doSlow = flags.doFullElectrostatics;
  doEnergy = flags.doEnergy;
  step = flags.step;

  if ( ! flags.doNonbonded ) {
    GBISP("GBIS[%d] noWork() don't do nonbonded\n",CkMyPe());
    for ( int i=0; i<hostedPatches.size(); ++i ) {
      patch_record &pr = master->patchRecords[hostedPatches[i]];
      pr.positionBox->skip();
      pr.forceBox->skip();
      //if (simParams->GBISOn) {
      //  pr.intRadBox->skip();
      //  pr.psiSumBox->skip();
      //  pr.bornRadBox->skip();
      //  pr.dEdaSumBox->skip();
      //  pr.dHdrPrefixBox->skip();
      //}
    }
    if ( master != this ) {
      computeMgr->sendNonbondedMICSlaveReady(masterPe,
			hostedPatches.size(),atomsChanged,sequence());
    }
    if ( reduction ) reduction->submit();
    return 1;
  }

  for ( int i=0; i<hostedPatches.size(); ++i ) {
    patch_record &pr = master->patchRecords[hostedPatches[i]];
    if (!simParams->GBISOn || gbisPhase == 1) {
      GBISP("GBIS[%d] noWork() P0[%d] open()\n",CkMyPe(), pr.patchID);
      // DMK - NOTE : When patches are pushed to the device individually, open called in patchReady function
      #if MIC_SUBMIT_ATOMS_ON_ARRIVAL == 0
        pr.x = pr.positionBox->open();
      #endif
      pr.xExt = pr.p->getCompAtomExtInfo();
    }

    //if (simParams->GBISOn) {
    //  if (gbisPhase == 1) {
    //    GBISP("GBIS[%d] noWork() P1[%d] open()\n",CkMyPe(),pr.patchID);
    //    pr.intRad     = pr.intRadBox->open();
    //    pr.psiSum     = pr.psiSumBox->open();
    //  } else if (gbisPhase == 2) {
    //    GBISP("GBIS[%d] noWork() P2[%d] open()\n",CkMyPe(),pr.patchID);
    //    pr.bornRad    = pr.bornRadBox->open();
    //    pr.dEdaSum    = pr.dEdaSumBox->open();
    //  } else if (gbisPhase == 3) {
    //    GBISP("GBIS[%d] noWork() P3[%d] open()\n",CkMyPe(),pr.patchID);
    //    pr.dHdrPrefix = pr.dHdrPrefixBox->open();
    //  }
    //  GBISP("opened GBIS boxes");
    //}
  }

  if ( master == this ) return 0; //work to do, enqueue as usual

  // message masterPe
  computeMgr->sendNonbondedMICSlaveReady(masterPe,
                                         hostedPatches.size(),
                                         atomsChanged,
                                         sequence()
                                        );

  //// DMK - NOTE : Skip the "remote" work for now
  //workStarted = 2; //1;
  #if MIC_MULTIPLE_KERNELS != 0
    workStarted = 1;
  #else
    workStarted = 2;
  #endif
  basePriority = COMPUTE_PROXY_PRIORITY;

  return 1;
}

void ComputeNonbondedMIC::doWork() {

  // DMK - Continue here... place the device pointer set in patchReady into the patch_pair records so that
  //   the computes can find the atom data on the device

  //// DMK - DEBUG
  //printf("[DEBUG:%d] :: ComputeNonbondedMIC::doWork() - Called...  this:%p\n", CkMyPe(), this); fflush(NULL);

  GBISP("C.N.MIC[%d]::doWork: seq %d, phase %d, workStarted %d\n", \
        CkMyPe(), sequence(), gbisPhase, workStarted);

  PATCH_PAIRS_REF;
  FORCE_LISTS_REF;

  if ( workStarted ) { //if work already started, check if finished
    if ( finishWork() ) {  // finished

      //// DMK - DEBUG
      //printf("[DEBUG] :: doWork() :: kernels all finished...\n"); fflush(NULL);

      workStarted = 0;
      basePriority = PROXY_DATA_PRIORITY;  // higher to aid overlap

      // DMK - DEBUG
      timestep++;

    } else {  // need to call again

      //// DMK - DEBUG
      //printf("[DEBUG] :: doWork() :: remote kernel finished, starting check for local...\n"); fflush(NULL);

      workStarted = 2;
      //basePriority = PROXY_RESULTS_PRIORITY;  // lower for local
      basePriority = PROXY_DATA_PRIORITY;  // lower for local
      if ( master == this && kernel_launch_state > 2 ) {
        #if MIC_TRACING != 0
          mic_tracing_polling_start = CmiWallTimer();
          mic_tracing_polling_count = 0;
        #endif
        mic_check_local_progress(this,0.);  // launches polling
      }
    }
    return;
  }

  #if MIC_TRACING != 0
    double doWork_start = CmiWallTimer();
  #endif

  //// DMK - NOTE : Skip the "remote" work for now
  //workStarted = 2; //1;
  #if MIC_MULTIPLE_KERNELS != 0
    workStarted = 1;
  #else
    workStarted = 2;
  #endif
  basePriority = COMPUTE_MIC_PRIORITY; //COMPUTE_PROXY_PRIORITY;

  Molecule *mol = Node::Object()->molecule;
  Parameters *params = Node::Object()->parameters;
  SimParameters *simParams = Node::Object()->simParameters;

  //execute only during GBIS phase 1, or if not using GBIS
  if (!simParams->GBISOn || gbisPhase == 1) {

    // bind new patches to GPU
    if ( atomsChanged || computesChanged ) {

      int npatches = activePatches.size();

      pairlistsValid = 0;
      pairlistTolerance = 0.;

      if ( computesChanged ) {

        computesChanged = 0;

	// Merge the local and remote active patch lists into a single list
        num_local_patch_records = localActivePatches.size();
        num_remote_patch_records = remoteActivePatches.size();
        npatches = num_local_patch_records + num_remote_patch_records;  // DMK
        activePatches.resize(npatches);                                 // DMK
        int *ap = activePatches.begin();
        for ( int i=0; i<num_local_patch_records; ++i ) {               // DMK
          *(ap++) = localActivePatches[i];
        }
        for ( int i=0; i<num_remote_patch_records; ++i ) {              // DMK
          *(ap++) = remoteActivePatches[i];
        }

        // sort computes by distance between patches
        #if MIC_SORT_COMPUTES != 0
          cr_sortop so(lattice);
          std::stable_sort(localComputeRecords.begin(),localComputeRecords.end(),so);
          std::stable_sort(remoteComputeRecords.begin(),remoteComputeRecords.end(),so);
        #endif

        // Merge the sorted lists of local and remote computes into a single list of computes 
        num_local_compute_records = localComputeRecords.size();
        num_remote_compute_records = remoteComputeRecords.size();
        computeRecords.resize(num_local_compute_records+num_remote_compute_records);
        compute_record *cr = computeRecords.begin();
        for ( int i=0; i<num_local_compute_records; ++i ) {
          *(cr++) = localComputeRecords[i];
        }
        for ( int i=0; i<num_remote_compute_records; ++i ) {
          *(cr++) = remoteComputeRecords[i];
        }

        // Set each patch record's localIndex and initialize force_list_size for each force_list
        force_lists.resize(npatches);
        for ( int i=0; i<npatches; ++i ) {
          patchRecords[activePatches[i]].localIndex = i;
          force_lists[i].force_list_size = 0;
        }

        // For the patch_pairs, one for each compute object...
        int ncomputes = computeRecords.size();
        patch_pairs.resize(ncomputes);
        for ( int i=0; i<ncomputes; ++i ) {

          ComputeNonbondedMIC::compute_record &cr = computeRecords[i];
          int lp1 = patchRecords[cr.pid[0]].localIndex;
          int lp2 = patchRecords[cr.pid[1]].localIndex;

          // Count the number of writers/force-contributers
          force_lists[lp1].force_list_size++;
          if (!cr.isSelf) {
            force_lists[lp2].force_list_size++;
	  }

          // Initialize the offset
          patch_pair &pp = patch_pairs[i];
          pp.offset.x = cr.offset.x;
          pp.offset.y = cr.offset.y;
          pp.offset.z = cr.offset.z;
        }

        // Set the force list information for each patch pair
        for ( int i=0; i<ncomputes; ++i ) {

          ComputeNonbondedMIC::compute_record &cr = computeRecords[i];
          patch_pair &pp = patch_pairs[i];

          int lp1 = patchRecords[cr.pid[0]].localIndex;
          pp.patch1_force_list_index = lp1;
          pp.patch1_force_list_size = force_lists[lp1].force_list_size;

          if (cr.isSelf) {
            pp.patch2_force_list_index = pp.patch1_force_list_index;
            pp.patch2_force_list_size = pp.patch1_force_list_size;
	  } else {
            int lp2 = patchRecords[cr.pid[1]].localIndex;
            pp.patch2_force_list_index = lp2;
            pp.patch2_force_list_size = force_lists[lp2].force_list_size;
	  }
        }

        CkPrintf("Pe %d has %d local and %d remote patches and %d local and %d remote computes.\n",
                 CkMyPe(), localActivePatches.size(), remoteActivePatches.size(),
                 localComputeRecords.size(), remoteComputeRecords.size()
                );

      }  // computesChanged

      //// DMK - DEBUG
      //printf("[DEBUG] :: 1.0 :: this:%p\n", this); fflush(NULL);

      // Count the number of atoms (and non-fixed atoms), recording the accumulated
      //   values as we go into the patch record and force list structures.
      int istart = 0;
      int flstart = 0;
      int vlstart = 0;
      int max_atoms_per_patch = 0;
      int i;
      for (i = 0; i < npatches; i++) {

        // If we have walked off the end of the local list, then take
        //   note of the numer of atoms so far (the local atom count).
        if ( i == localActivePatches.size() ) {
          num_local_atom_records = istart;
        }

        // Record the current offsets as the offsets as the
        // beginning of this force list.
        force_lists[i].force_list_start = flstart;
        force_lists[i].force_output_start = istart;
        force_lists[i].atom_start = istart;
        //force_lists[i].virial_list_start = vlstart;
        //force_lists[i].virial_output_start = 16*i;

        // Record the current offset as the patch record's starting offset.
        patch_record &pr = patchRecords[activePatches[i]];
        pr.localStart = istart;

        // Count number of atoms (and non-fixed atoms), recording the actual
        //   counts to the patch record and then rounding up for alignment
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
        force_lists[i].patch_size = natoms; // DMK - nfreeatoms;
        natoms = (natoms + 15) & (~15);
        nfreeatoms = (nfreeatoms + 15) & (~15);
        //if ( natoms & 15 ) { natoms += 16 - (natoms & 15); }  // Round up 16 (DMK - NOTE)
        //if ( nfreeatoms & 15 ) { nfreeatoms += 16 - (nfreeatoms & 15); }  // Round up 16 (DMK - NOTE)
        force_lists[i].patch_stride = natoms; // DMK - nfreeatoms;

        // Update the offsets by the atom counts for this patch record.
        flstart += natoms * force_lists[i].force_list_size; // DMK - nfreeatoms * force_lists[i].force_list_size;
        vlstart += 16 * force_lists[i].force_list_size;
        istart += natoms;  // already rounded up

        force_lists[i].force_list_size = 0;  // rebuild below
      }
      // Handle the case where all are local for recording local atom count.
      if ( i == localActivePatches.size() ) {
        num_local_atom_records = istart;
      }

      //// DMK - DEBUG
      //printf("[DEBUG] :: 2.0 :: this:%p\n", this); fflush(NULL);

      // Record final offsets (i.e. lengths/counts).
      num_force_records = flstart;
      num_atom_records = istart;
      num_remote_atom_records = num_atom_records - num_local_atom_records;

      // Allocate the memory for the atom and force information
      if ( num_atom_records > num_atom_records_allocated ) {
        if ( num_atom_records_allocated ) {
          _MM_FREE_WRAPPER(atom_params); //delete [] atom_params;
          _MM_FREE_WRAPPER(atoms); //delete [] atoms;
          _MM_FREE_WRAPPER(forces); //delete [] forces;
          _MM_FREE_WRAPPER(slow_forces); //delete [] slow_forces;
          //if (simParams->GBISOn) {
          //  delete [] intRad0H;//6 GBIS arrays
          //  delete [] intRadSH;
          //  delete [] psiSumH;
          //  delete [] bornRadH;
          //  delete [] dEdaSumH;
          //  delete [] dHdrPrefixH;
          //}
        }
        num_atom_records_allocated = 1.1 * num_atom_records + 1;
        atom_params = (atom_param*)_MM_MALLOC_WRAPPER(num_atom_records_allocated * sizeof(atom_param), 64, "atom_params"); //new atom_param[num_atom_records_allocated];
        atoms = (atom*)_MM_MALLOC_WRAPPER(num_atom_records_allocated * sizeof(atom), 64, "atoms"); //new atom[num_atom_records_allocated];
        if (atom_params == NULL || atoms == NULL) { NAMD_die("Unable to allocate atoms in ComputeNonbondedMIC::doWork"); }
        forces = (double4*)_MM_MALLOC_WRAPPER(num_atom_records_allocated * sizeof(double4), 64, "forces"); //new double4[num_atom_records_allocated];
        slow_forces = (double4*)_MM_MALLOC_WRAPPER(num_atom_records_allocated * sizeof(double4), 64, "slow_forces"); //new double4[num_atom_records_allocated];
        //allocate GBIS memory
        //if (simParams->GBISOn) {
        //  intRad0H = new float[num_atom_records_allocated];
        //  intRadSH = new float[num_atom_records_allocated];
        //  psiSumH = new GBReal[num_atom_records_allocated];
        //  bornRadH = new float[num_atom_records_allocated];
        //  dEdaSumH = new GBReal[num_atom_records_allocated];
        //  dHdrPrefixH = new float[num_atom_records_allocated];
        //}
        if (forces == NULL || slow_forces == NULL || (simParams->GBISOn &&
            (intRad0H == NULL || intRadSH == NULL || psiSumH == NULL ||
             bornRadH == NULL || dEdaSumH == NULL || dHdrPrefixH == NULL))) {
          NAMD_die("Unable to allocate forces in ComputeNonbondedMIC::doWork");
        }
      }

      //// DMK - DEBUG
      //printf("[DEBUG] :: 3.0 :: this:%p\n", this); fflush(NULL);

      // DMK - NOTE - Continue filling in the patch pair records
      int bfstart = 0;
      int ncomputes = computeRecords.size();
      for ( int i=0; i<ncomputes; ++i ) {

        ComputeNonbondedMIC::compute_record &cr = computeRecords[i];

        int p1 = cr.pid[0];
        int p2 = cr.pid[1];
        int lp1 = patchRecords[p1].localIndex;
        int lp2 = patchRecords[p2].localIndex;
        patch_pair &pp = patch_pairs[i];
        pp.patch1_atom_start = patchRecords[p1].localStart;
        pp.patch1_force_start = force_lists[lp1].force_list_start
                              + (force_lists[lp1].patch_stride * force_lists[lp1].force_list_size);
        pp.patch1_size = patchRecords[p1].numAtoms;
        pp.patch1_force_size = patchRecords[p1].numAtoms; 

        if (cr.isSelf) {
          pp.patch2_atom_start = pp.patch1_atom_start;
          pp.patch2_force_start = pp.patch1_force_start;
          pp.patch2_size = pp.patch1_size;
          pp.patch2_force_size = pp.patch1_force_size;
	} else {
          pp.patch2_atom_start = patchRecords[p2].localStart;
          pp.patch2_force_start = force_lists[lp2].force_list_start
                                + (force_lists[lp2].patch_stride * force_lists[lp2].force_list_size);
          pp.patch2_size = patchRecords[p2].numAtoms;
          pp.patch2_force_size = patchRecords[p2].numAtoms;
	}

        // Get a custom pairlist cutoff for each patch pair
        pp.plcutoff = ComputeNonbondedUtil::cutoff +
          patchRecords[p1].p->flags.pairlistTolerance +
          patchRecords[p2].p->flags.pairlistTolerance;

        // Record the parts
        #if MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING != 0
          pp.numParts = cr.numParts;
          pp.part = cr.part;
        #endif

        // Record the patch centers
        Vector p1_center = micCompute->patchMap->center(p1);
        pp.patch1_center_x = p1_center.x;
        pp.patch1_center_y = p1_center.y;
        pp.patch1_center_z = p1_center.z;
        Vector p2_center = micCompute->patchMap->center(p2);
        pp.patch2_center_x = p2_center.x;
        pp.patch2_center_y = p2_center.y;
        pp.patch2_center_z = p2_center.z;

        // DMK - DEBUG
        pp.p1 = p1;
        pp.p2 = p2;
        pp.cid = cr.c;

        pp.block_flags_start = bfstart;
        bfstart += ((pp.patch1_force_size + 127) >> 7) << 5;

        //pp.virial_start = force_lists[lp1].virial_list_start +
        //  16 * force_lists[lp1].force_list_size;

        force_lists[lp1].force_list_size++;
        if (!cr.isSelf) { force_lists[lp2].force_list_size++; }

      } // for ncomputes

      //// DMK - DEBUG
      //printf("[DEBUG] :: 4.0 :: this:%p\n", this); fflush(NULL);

      #if 0
        CkPrintf("Pe %d mic_bind_patch_pairs %d %d %d %d %d\n", CkMyPe(),
	  patch_pairs.size(), force_lists.size(),
          num_atom_records, num_force_records,
          max_atoms_per_patch);
      #endif

      int totalmem = patch_pairs.size() * sizeof(patch_pair) +
                  force_lists.size() * sizeof(force_list) +
                  num_force_records * sizeof(double4) +
                  num_atom_records * sizeof(atom) +
                  num_atom_records * sizeof(atom_param) +
                  num_atom_records * sizeof(double4);
      int totalcopy = num_atom_records * ( sizeof(atom) + sizeof(double4) );
      /*
      CkPrintf("Pe %d allocating %d MB of MIC memory, will copy %d kB per step\n",
  			CkMyPe(), totalmem >> 20, totalcopy >> 10);
      */

      // Push the data structures just created, such as patch pairs, to the device.
      mic_bind_patch_pairs_only(myDevice, patch_pairs.begin(), patch_pairs.size(), patch_pairs.bufSize());
      mic_bind_force_lists_only(myDevice, force_lists.begin(), force_lists.size(), force_lists.bufSize());
      mic_bind_atoms_only(myDevice, atoms, atom_params, forces, slow_forces, num_atom_records, num_atom_records_allocated);
      mic_bind_force_buffers_only(myDevice, num_force_records);

    } // atomsChanged || computesChanged

    double charge_scaling = sqrt(COULOMB * scaling * dielectric_1);

    Flags &flags = patchRecords[hostedPatches[0]].p->flags;
    float maxAtomMovement = 0.;
    float maxPatchTolerance = 0.;

    //// DMK - DEBUG
    //printf("[DEBUG] :: 5.0 :: this:%p\n", this); fflush(NULL);

    for ( int i=0; i<activePatches.size(); ++i ) {

      patch_record &pr = patchRecords[activePatches[i]];

      float maxMove = pr.p->flags.maxAtomMovement;
      if ( maxMove > maxAtomMovement ) maxAtomMovement = maxMove;

      float maxTol = pr.p->flags.pairlistTolerance;
      if ( maxTol > maxPatchTolerance ) maxPatchTolerance = maxTol;

      int start = pr.localStart;
      int n = pr.numAtoms;
      int n_16 = (n + 15) & (~15);
      const CompAtom *a = pr.x;
      const CompAtomExt *aExt = pr.xExt;

      //// DMK - DEBUG
      //printf("[DEBUG] :: 5.1 :: this:%p, atom_params:%p\n", this, atom_params); fflush(NULL);

      if ( atomsChanged ) {

        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          atom_param *ap = atom_params + start;
          for (int k = 0; k < n; k++) {
            int j = aExt[k].sortOrder;
            ap[k].vdw_type = a[j].vdwType;
            ap[k].index = aExt[j].id;
            #ifdef MEM_OPT_VERSION
              ap[k].excl_index = exclusionsByAtom[aExt[j].exclId].y;
              ap[k].excl_maxdiff = exclusionsByAtom[aExt[j].exclId].x;
            #else
              ap[k].excl_index = exclusionsByAtom[aExt[j].id].y;
              ap[k].excl_maxdiff = exclusionsByAtom[aExt[j].id].x;
            #endif
	  }
        #else
          atom_param *ap = atom_params + start;
          int *ap_vdwType = ((int*)ap) + (0 * n_16);
          int *ap_index = ((int*)ap) + (1 * n_16);
          int *ap_exclIndex = ((int*)ap) + (2 * n_16);
          int *ap_exclMaxDiff = ((int*)ap) + (3 * n_16);
          for ( int k=0; k<n; ++k ) {
            int j = aExt[k].sortOrder;
            ap_vdwType[k] = a[j].vdwType;
            ap_index[k] = aExt[j].id;
            #ifdef MEM_OPT_VERSION
              ap_exclIndex[k] = exclusionsByAtom[aExt[j].exclId].y;
              ap_exclMaxDiff[k] = exclusionsByAtom[aExt[j].exclId].x;
            #else
              ap_exclIndex[k] = exclusionsByAtom[aExt[j].id].y;
              ap_exclMaxDiff[k] = exclusionsByAtom[aExt[j].id].x;
            #endif
          }
        #endif

      } // end if (atomsChanged)

      { // start block
        const CudaAtom *ac = pr.p->getCudaAtomList();
        atom *ap = atoms + start;
        #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
          memcpy(ap, ac, sizeof(atom)*n);
        #else
          memcpy(ap, ac, sizeof(atom)*n_16);
        #endif

	//// DMK - DEBUG
        //for (int i = 0; i < n; i++) {
        //  printf("[DEBUG:PID-%04d-ATOMS-%04d] :: ap[%06d] = { x:%.4le, y:%.4le, z:%.4le, q:%.4le }\n",
        //         pr.patchID, timestep, i, ap[i].x, ap[i].y, ap[i].z, ap[i].charge
        //        );
        //}

      } // end block

    } // end for (i < activePatches.size())

    //GBISP("finished active patches\n")

    //CkPrintf("maxMovement = %f  maxTolerance = %f  save = %d  use = %d\n",
    //  maxAtomMovement, maxPatchTolerance,
    //  flags.savePairlists, flags.usePairlists);

    savePairlists = 0;
    usePairlists = 0;
    if ( flags.savePairlists ) {
      savePairlists = 1;
      usePairlists = 1;
    } else if ( flags.usePairlists ) {
      if ( ! pairlistsValid ||
           ( 2. * maxAtomMovement > pairlistTolerance ) ) {
        reduction->item(REDUCTION_PAIRLIST_WARNINGS) += 1;
      } else {
        usePairlists = 1;
      }
    }
    if ( ! usePairlists ) {
      pairlistsValid = 0;
    }
    float plcutoff = cutoff;
    if ( savePairlists ) {
      pairlistsValid = 1;
      pairlistTolerance = 2. * maxPatchTolerance;
      plcutoff += pairlistTolerance;
    }
    plcutoff2 = plcutoff * plcutoff;

    //CkPrintf("plcutoff = %f  listTolerance = %f  save = %d  use = %d\n",
    //  plcutoff, pairlistTolerance, savePairlists, usePairlists);

  } // !GBISOn || gbisPhase == 1

  //// Do GBIS
  //if (simParams->GBISOn) {
  //  //open GBIS boxes depending on phase
  //  for ( int i=0; i<activePatches.size(); ++i ) {
  //    patch_record &pr = master->patchRecords[activePatches[i]];
  //    GBISP("doWork[%d] accessing arrays for P%d\n",CkMyPe(),gbisPhase);
  //    if (gbisPhase == 1) {
  //      //Copy GBIS intRadius to Host
  //      if (atomsChanged) {
  //        float *intRad0 = intRad0H + pr.localStart;
  //        float *intRadS = intRadSH + pr.localStart;
  //        for ( int k=0; k<pr.numAtoms; ++k ) {
  //          int j = pr.xExt[k].sortOrder;
  //          intRad0[k] = pr.intRad[2*j+0];
  //          intRadS[k] = pr.intRad[2*j+1];
  //        }
  //      }
  //    } else if (gbisPhase == 2) {
  //      float *bornRad = bornRadH + pr.localStart;
  //      for ( int k=0; k<pr.numAtoms; ++k ) {
  //        int j = pr.xExt[k].sortOrder;
  //        bornRad[k] = pr.bornRad[j];
  //      }
  //    } else if (gbisPhase == 3) {
  //      float *dHdrPrefix = dHdrPrefixH + pr.localStart;
  //      for ( int k=0; k<pr.numAtoms; ++k ) {
  //        int j = pr.xExt[k].sortOrder;
  //        dHdrPrefix[k] = pr.dHdrPrefix[j];
  //      }
  //    } // end phases
  //  } // end for patches
  //} // if GBISOn

  #if MIC_TRACING != 0
    MIC_TRACING_RECORD(MIC_EVENT_FUNC_DOWORK, doWork_start, CmiWallTimer());
  #endif

  kernel_time = CkWallTimer();
  #if MIC_MULTIPLE_KERNELS != 0
    kernel_launch_state = 1;
  #else
    kernel_launch_state = 2;
  #endif
  if ( gpu_is_mine ) recvYieldDevice(-1);
}

void mic_check_remote_calc(void *arg, double) {
  if (mic_check_remote_kernel_complete(myDevice)) {
    computeMgr->sendYieldDevice(next_pe_sharing_gpu);
  } else {
    MIC_POLL(mic_check_remote_calc, arg);
  }
}

void mic_check_local_calc(void *arg, double) {
  if (mic_check_local_kernel_complete(myDevice)) {
    computeMgr->sendYieldDevice(next_pe_sharing_gpu);
  } else {
    MIC_POLL(mic_check_local_calc, arg);
  }
}

void ComputeNonbondedMIC::recvYieldDevice(int pe) {

  //// DMK - DEBUG
  //printf("[DEBUG:%d] :: ComputeNonbondedMIC::recvYieldDevice() - Called...\n", CkMyPe()); fflush(NULL);

  GBISP("C.N.MIC[%d]::recvYieldDevice: seq %d, workStarted %d, \
        gbisPhase %d, kls %d, from pe %d\n", CkMyPe(), sequence(), \
        workStarted, gbisPhase, kernel_launch_state, pe)

  mic_position3_t lata, latb, latc;
  lata.x = lattice.a().x;
  lata.y = lattice.a().y;
  lata.z = lattice.a().z;
  latb.x = lattice.b().x;
  latb.y = lattice.b().y;
  latb.z = lattice.b().z;
  latc.x = lattice.c().x;
  latc.y = lattice.c().y;
  latc.z = lattice.c().z;
  SimParameters *simParams = Node::Object()->simParameters;

  switch ( kernel_launch_state ) {

    ////////////////////////////////////////////////////////////
    // Remote
    case 1:

      #if MIC_MULTIPLE_KERNELS != 0

        GBISP("C.N.MIC[%d]::recvYieldDeviceR: case 1\n", CkMyPe())
        ++kernel_launch_state;
        gpu_is_mine = 0;
        remote_submit_time = CkWallTimer();

        if (!simParams->GBISOn || gbisPhase == 1) {

          // DMK - NOTE : GBIS is not supported in this port yet, so cause a runtime error
          //   if it has been enabled and execution has made it this far.
          if (simParams->GBISOn) {
            NAMD_die("Unsupported feature (DMK33949330)");
          }

          //// DMK - DEBUG
          //printf("[DEBUG:%d] :: << issuing remote kernel >>\n", CkMyPe()); fflush(NULL);

          #if MIC_TRACING != 0
            mic_tracing_offload_start_remote = CmiWallTimer();
          #endif

          // Issue the remote kernel
          mic_nonbonded_forces(myDevice, 1,
                               num_local_atom_records,
                               num_local_compute_records,
                               num_local_patch_records,
                               lata, latb, latc,
                               doSlow, doEnergy,
                               usePairlists, savePairlists,
                               atomsChanged
                              );

          // Start the polling check for the completion of the remote kernel
          #if MIC_TRACING != 0
            mic_tracing_polling_start = CmiWallTimer();
            mic_tracing_polling_count = 0;
          #endif
          MIC_POLL(mic_check_remote_progress, this);

          //// DMK - DEBUG - Only use one kernel for now
          ////#if MIC_TRACING != 0
          ////  mic_tracing_offload_start_remote = CmiWallTimer();
          ////#endif
          ////mic_nonbonded_forces(myDevice, 1,
          ////                     lata, latb, latc,
          ////                     doSlow, doEnergy,
          ////                     usePairlists, savePairlists,
          ////                     atomsChanged
          ////                    );
          ////MIC_POLL(mic_check_remote_progress, this);
          ////if ( shared_gpu && ! mergegrids ) {
          ////  MIC_POLL(mic_check_remote_calc, this);
          ////  break;
          ////}
        }

      #else // if MIC_MULTIPLE_KERNELS != 0

        // DMK - DEBUG
        { static int __reached_here_flag__ = 0;
          if (!__reached_here_flag__) {
            CkPrintf("[MIC] :: PE %d reached remote case with MIC_MULTIPLE_KERNELS disabled\n", CkMyPe());
            __reached_here_flag__ = 1;
	  }
        }

      #endif

      // NOTE : Fall through to next case (i.e. break not missing)

    ////////////////////////////////////////////////////////////
    // Local
    case 2:

      GBISP("C.N.MIC[%d]::recvYieldDeviceL: case 2\n", CkMyPe())
      ++kernel_launch_state;
      gpu_is_mine = 0;

      if (!simParams->GBISOn || gbisPhase == 1) {

        // DMK - NOTE : GBIS is not supported in this port yet, so cause a runtime error
        //   if it has been enabled and execution has made it this far.
        if (simParams->GBISOn) {
          NAMD_die("Unsupported feature (DMK83620583)");
        }

        // DMK - TIMING - NOTE : Only local being submitted for now, so
        //   take the local time now
        local_submit_time = CkWallTimer();

        #if MIC_TRACING != 0
          mic_tracing_offload_start_local = CmiWallTimer();
        #endif

	//// DMK - DEBUG
        //printf("[DEBUG:%d] :: << issuing local kernel >>\n", CkMyPe()); fflush(NULL);

        // Issue the local kernel
        mic_nonbonded_forces(myDevice, 0,
                             num_local_atom_records,
                             num_local_compute_records,
                             num_local_patch_records,
                             lata, latb, latc,
                             doSlow, doEnergy,
                             usePairlists, savePairlists,
                             atomsChanged
                            );

        if ( workStarted == 2 ) {
          // Start the polling check for the completion of the local kernel
          #if MIC_TRACING != 0
            mic_tracing_polling_start = CmiWallTimer();
            mic_tracing_polling_count = 0;
          #endif
          MIC_POLL(mic_check_local_progress, this);
        }

        //if ( shared_gpu && ! mergegrids ) {
        //  MIC_POLL(mic_check_local_calc, this);
        //}

      } // end if (!simParams->GBISOn || gbisPhase == 1)

    default:

      GBISP("C.N.MIC[%d]::recvYieldDevice: case default\n", CkMyPe())
      gpu_is_mine = 1;
      break;

  } // switch

  GBISP("C.N.MIC[%d]::recvYieldDevice: DONE\n", CkMyPe())
}


void ComputeNonbondedMIC::messageFinishWork() {

  mic_errcheck("at mic stream completed");

  for ( int i = 0; i < numSlaves; ++i ) {
    computeMgr->sendNonbondedMICSlaveEnqueue(slaves[i],slavePes[i],sequence(),priority(),workStarted);
  }

  WorkDistrib::messageEnqueueWork(this);
}

//dtanner
int ComputeNonbondedMIC::finishWork() {

  //// DMK - DEBUG
  //printf("[DEBUG] :: ComputeNonbondedMIC::finishWork() - Called...\n"); fflush(NULL);

  #if MIC_TRACING != 0
    double finishWork_start = CmiWallTimer();
  #endif

  GBISP("C.N.MIC[%d]::fnWork: workStarted %d, phase %d\n", \
  CkMyPe(), workStarted, gbisPhase)

  Molecule *mol = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;

  #if MIC_MULTIPLE_KERNELS != 0
    ResizeArray<int> &patches( workStarted == 1 ? remoteHostedPatches : localHostedPatches );
    mic_kernel_data* &kernel_data = (workStarted == 1) ? (device__remote_kernel_data) : (device__local_kernel_data);
  #else
    ResizeArray<int> &patches( localHostedPatches );
    mic_kernel_data* &kernel_data = device__local_kernel_data;
  #endif

  //// DMK - DEBUG
  //printf("[DEBUG:%d] :: ComputeNonbondedMIC::finishWork() - p:%d (r:%d, l:%d)\n",
  //       CkMyPe(), patches.size(), remoteHostedPatches.size(), localHostedPatches.size()
  //      ); fflush(NULL);

  // DMK - NOTE : Open the force boxes for the patches so forces can be contributed
  if ( !simParams->GBISOn || gbisPhase == 1 ) {
    for ( int i=0; i<patches.size(); ++i ) {
      patch_record &pr = master->patchRecords[patches[i]];
      GBISP("GBIS[%d] fnWork() P0[%d] force.open()\n",CkMyPe(), pr.patchID);
      pr.r = pr.forceBox->open();
    }
  } // !GBISOn || gbisPhase==1

  //// DMK - DEBUG
  //#if MIC_MULTIPLE_KERNELS != 0
  //  //if (workStarted == 2) { // NOTE: Only print after the local kernel has completed
  //#endif
  //    if (this == master) {
  //      for (int i = 0; i < num_atom_records; i++) {
  //        printf("[DEBUG:%d:FORCE-%d] :: force[%d] = { x:%.3lf, y:%.3lf, z:%.3lf }\n",
  //               CkMyPe(), workStarted, i, forces[i].x, forces[i].y, forces[i].z
  //              );
  //      }
  //    }
  //#if MIC_MULTIPLE_KERNELS != 0
  //  //}
  //#endif

  // DMK - NOTE : For each patch...
  for ( int i=0; i<patches.size(); ++i ) {

    // CkPrintf("Pe %d patch %d of %d pid %d\n",CkMyPe(),i,patches.size(),patches[i]);
    patch_record &pr = master->patchRecords[patches[i]];
    int start = pr.localStart;
    const CompAtomExt *aExt = pr.xExt;

    //// DMK - DEBUG
    //printf("[DEBUG:PID-%04d(%04d)-%04d-%01d-PATCHCONTRIB] - Contributing forces\n",
    //       patches[i], pr.patchID, timestep, workStarted
    //      );

    // DMK - NOTE : If this is the end of the timestep for this compute
    if ( !simParams->GBISOn || gbisPhase == 3 ) {

      // DMK - NOTE : Contribute the calculated forces and slow_forces from
      //   this compute to the patch.
      //int nfree = pr.numFreeAtoms;
      int nAtoms = pr.numAtoms;
      int nAtoms_16 = (nAtoms + 15) & (~15);
      pr.f = pr.r->f[Results::nbond];
      Force *f = pr.f;
      Force *f_slow = pr.r->f[Results::slow];
      const CompAtom *a = pr.x;
      const CompAtomExt *aExt = pr.xExt;
      // int *ao = atom_order + start;
      //float4 *af = master->forces + start;
      //float4 *af_slow = master->slow_forces + start;

      #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0

        double4 *af = master->forces + start;
        for (int k = 0; k < nAtoms; k++) {
          int j = aExt[k].sortOrder;
          f[j].x += af[k].x;
          f[j].y += af[k].y;
          f[j].z += af[k].z;

          //// DMK - DEBUG
          //printf("[DEBUG:PID-%04d-FORCES-%04d] :: af[%06d] = { x:%.4le, y:%.4le, z:%.4le }\n",
          //       pr.patchID, timestep, k, af[k].x, af[k].y, af[k].z
          //      );

	}
        if (doSlow) {
          double4 *af_slow = master->slow_forces + start;
          for (int k = 0; k < nAtoms; k++) {
            int j = aExt[k].sortOrder;
            f_slow[j].x += af_slow[k].x;
            f_slow[j].y += af_slow[k].y;
            f_slow[j].z += af_slow[k].z;
          }
	}

      #else

        double4 *af = master->forces + start;
        double4 *af_slow = master->slow_forces + start;
        double *af_x = ((double*)af) + (0 * nAtoms_16);
        double *af_y = ((double*)af) + (1 * nAtoms_16);
        double *af_z = ((double*)af) + (2 * nAtoms_16);
        double *af_w = ((double*)af) + (3 * nAtoms_16);
        double *af_slow_x = ((double*)af_slow) + (0 * nAtoms_16);
        double *af_slow_y = ((double*)af_slow) + (1 * nAtoms_16);
        double *af_slow_z = ((double*)af_slow) + (2 * nAtoms_16);
        double *af_slow_w = ((double*)af_slow) + (3 * nAtoms_16);

        for (int k = 0; k < nAtoms; k++) {
          int j = aExt[k].sortOrder;
          f[j].x += af_x[k];
          f[j].y += af_y[k];
          f[j].z += af_z[k];
          if (doSlow) {
            f_slow[j].x += af_slow_x[k];
            f_slow[j].y += af_slow_y[k];
            f_slow[j].z += af_slow_z[k];
	  }
	}

      #endif

      // DMK - DEBUG
      #if 0
        if (i <= 10) {
          for (int k = 0; k < nAtoms; k++) {
            printf("hfo :: i:%d, k:%d :: f = { %+.3le %+.3le %+.3le %+.3le }\n",
                   #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
                     i, k, af[k].x, af[k].y, af[k].z, af[k].w
                   #else
                     i, k, af_x[k], af_y[k], af_z[k], af_w[k]
                   #endif
		  );
	  }
        }
      #endif

      // Verify exclusion counts on 1 PE (reported as force's 'w' value)
      #if MIC_EXCL_CHECKSUM != 0
        if ( CkNumPes() == 1 ) {
          for (int k = 0; k < nAtoms; ++k) {
	  //for ( int j = 0; j < nAtoms; j++) { //int k=0; k<nfree; ++k ) {
            int j = aExt[k].sortOrder;
            #ifdef MEM_OPT_VERSION
	      int excl_expected = mol->exclSigPool[aExt[j].exclId].fullExclCnt; // DMK : Don't count self? + 1;
              excl_expected += mol->exclSigPool[aExt[j].exclId].modExclCnt; // DMK : and count modified
            #else
	      int excl_expected = mol->get_full_exclusions_for_atom(aExt[j].id)[0]; // DMK : Don't count self? + 1;
	      excl_expected += mol->get_mod_exclusions_for_atom(aExt[j].id)[0]; // DMK : and count modified
            #endif
            #if MIC_HANDCODE_FORCE_SOA_VS_AOS != 0
              int excl_actual = (int)(af[k].w);
            #else
              int excl_actual = (int)(af_w[k]);
            #endif
	    if (excl_actual != excl_expected) {
              static int __limit = 1000;
              static int __count = 0;
              if (++__count <= __limit) {
                CkPrintf("%d:%d(%d) atom %d found %d exclusions but expected %d\n",
                         i, j, k, aExt[j].id, excl_actual, excl_expected );
                if (__count == __limit) { CkPrintf("  Limit reached (not printing remaining exclusion check failures).\n"); }
              }
            } // end if excl check fail
          } // end for (j < nAtoms)
        } // end if ( CkNumPes() == 1 )
      #endif

    } // !GBISOn || gbisPhase == 3

    #if 0
      if ( i % 31 == 0 ) for ( int j=0; j<3; ++j ) {
        CkPrintf("Pe %d patch %d atom %d (%f %f %f) force %f\n", CkMyPe(), i,
	         j, pr.x[j].position.x, pr.x[j].position.y, pr.x[j].position.z,
	         af[j].w);
      }
    #endif

    //// Close Boxes depending on Phase
    //if (simParams->GBISOn) {
    //
    //  if (gbisPhase == 1) {
    //
    //    //Copy dEdaSum from Host to Patch Box
    //    GBReal *psiSumMaster = master->psiSumH + start;
    //    for ( int k=0; k<pr.numAtoms; ++k ) {
    //      int j = aExt[k].sortOrder;
    //      pr.psiSum[j] += psiSumMaster[k];
    //    }
    //    GBISP("C.N.MIC[%d]::fnWork: P1 psiSum.close()\n", CkMyPe());
    //    pr.psiSumBox->close(&(pr.psiSum));
    //
    //  } else if (gbisPhase == 2) {
    //
    //    //Copy dEdaSum from Host to Patch Box
    //    GBReal *dEdaSumMaster = master->dEdaSumH + start;
    //    for ( int k=0; k<pr.numAtoms; ++k ) {
    //      int j = aExt[k].sortOrder;
    //      pr.dEdaSum[j] += dEdaSumMaster[k];
    //    }
    //    GBISP("C.N.MIC[%d]::fnWork: P2 dEdaSum.close()\n", CkMyPe());
    //    pr.dEdaSumBox->close(&(pr.dEdaSum));
    //
    //  } else if (gbisPhase == 3) {
    //
    //    GBISP("C.N.MIC[%d]::fnWork: P3 all.close()\n", CkMyPe());
    //    pr.intRadBox->close(&(pr.intRad)); //box 6
    //    pr.bornRadBox->close(&(pr.bornRad)); //box 7
    //    pr.dHdrPrefixBox->close(&(pr.dHdrPrefix)); //box 9
    //    pr.positionBox->close(&(pr.x)); //box 0
    //    pr.forceBox->close(&(pr.r));
    //
    //  } //end phases
    //
    //} else { //not GBIS
    //
    //  GBISP("C.N.MIC[%d]::fnWork: pos/force.close()\n", CkMyPe());
      pr.positionBox->close(&(pr.x));
      pr.forceBox->close(&(pr.r));
    //}

  }  // end for (i<patches.size())

  #if 0
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
  #endif

  // DMK - NOTE : Contribute virial values
  if ( master == this && (!simParams->GBISOn || gbisPhase == 3) && workStarted == 2 ) {

    double virial_xx = device__local_kernel_data->virial_xx;
    double virial_xy = device__local_kernel_data->virial_xy;
    double virial_xz = device__local_kernel_data->virial_xz;
    double virial_yy = device__local_kernel_data->virial_yy;
    double virial_yz = device__local_kernel_data->virial_yz;
    double virial_zz = device__local_kernel_data->virial_zz;
    double fullElectVirial_xx = device__local_kernel_data->fullElectVirial_xx;
    double fullElectVirial_xy = device__local_kernel_data->fullElectVirial_xy;
    double fullElectVirial_xz = device__local_kernel_data->fullElectVirial_xz;
    double fullElectVirial_yy = device__local_kernel_data->fullElectVirial_yy;
    double fullElectVirial_yz = device__local_kernel_data->fullElectVirial_yz;
    double fullElectVirial_zz = device__local_kernel_data->fullElectVirial_zz;
    double vdwEnergy = device__local_kernel_data->vdwEnergy;
    double electEnergy = device__local_kernel_data->electEnergy;
    double fullElectEnergy = device__local_kernel_data->fullElectEnergy;
    #if MIC_MULTIPLE_KERNELS != 0
      virial_xx += device__remote_kernel_data->virial_xx;
      virial_xy += device__remote_kernel_data->virial_xy;
      virial_xz += device__remote_kernel_data->virial_xz;
      virial_yy += device__remote_kernel_data->virial_yy;
      virial_yz += device__remote_kernel_data->virial_yz;
      virial_zz += device__remote_kernel_data->virial_zz;
      fullElectVirial_xx += device__remote_kernel_data->fullElectVirial_xx;
      fullElectVirial_xy += device__remote_kernel_data->fullElectVirial_xy;
      fullElectVirial_xz += device__remote_kernel_data->fullElectVirial_xz;
      fullElectVirial_yy += device__remote_kernel_data->fullElectVirial_yy;
      fullElectVirial_yz += device__remote_kernel_data->fullElectVirial_yz;
      fullElectVirial_zz += device__remote_kernel_data->fullElectVirial_zz;
      vdwEnergy += device__remote_kernel_data->vdwEnergy;
      electEnergy += device__remote_kernel_data->electEnergy;
      fullElectEnergy += device__remote_kernel_data->fullElectEnergy;
    #endif

    //// DMK - DEBUG
    //#if MIC_MULTIPLE_KERNELS != 0
    //  printf("[DEBUG:FWl-%04d-%d] :: virial = { %.3lf %.3lf %.3lf  %.3lf %.3lf  %.3lf }\n",
    //         timestep, workStarted, device__local_kernel_data->virial_xx, device__local_kernel_data->virial_xy, device__local_kernel_data->virial_xz, device__local_kernel_data->virial_yy, device__local_kernel_data->virial_yz, device__local_kernel_data->virial_zz
    //        );
    //  printf("[DEBUG:FWr-%04d-%d] :: virial = { %.3lf %.3lf %.3lf  %.3lf %.3lf  %.3lf }\n",
    //         timestep, workStarted, device__remote_kernel_data->virial_xx, device__remote_kernel_data->virial_xy, device__remote_kernel_data->virial_xz, device__remote_kernel_data->virial_yy, device__remote_kernel_data->virial_yz, device__remote_kernel_data->virial_zz
    //        );
    //#endif
    //printf("[DEBUG:FW-%04d-%d] :: virial = { %.3lf %.3lf %.3lf  %.3lf %.3lf  %.3lf }\n",
    //       timestep, workStarted, virial_xx, virial_xy, virial_xz, virial_yy, virial_yz, virial_zz
    //      );
    //printf("[DEBUG:FW-%04d-%d] :: energy = { vdw:%.3lf elect:%.3lf fullElect:%.3lf }\n",
    //       timestep, workStarted, vdwEnergy, electEnergy, fullElectEnergy
    //      );

    // DMK - NOTE : Contribute virial
    Tensor virial_tensor;
    virial_tensor.xx = virial_xx;
    virial_tensor.xy = virial_xy;
    virial_tensor.xz = virial_xz;
    virial_tensor.yx = virial_xy;
    virial_tensor.yy = virial_yy;
    virial_tensor.yz = virial_yz;
    virial_tensor.zx = virial_xz;
    virial_tensor.zy = virial_yz;
    virial_tensor.zz = virial_zz;
    //virial_tensor.xx = kernel_data->virial_xx;
    //virial_tensor.xy = kernel_data->virial_xy;
    //virial_tensor.xz = kernel_data->virial_xz;
    //virial_tensor.yx = kernel_data->virial_xy;
    //virial_tensor.yy = kernel_data->virial_yy;
    //virial_tensor.yz = kernel_data->virial_yz;
    //virial_tensor.zx = kernel_data->virial_xz;
    //virial_tensor.zy = kernel_data->virial_yz;
    //virial_tensor.zz = kernel_data->virial_zz;
    // DMK - TODO | FIXME : GBIS support needed eventually
    ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NBOND, virial_tensor);
    if (doEnergy) {
      reduction->item(REDUCTION_LJ_ENERGY) += vdwEnergy;
      reduction->item(REDUCTION_ELECT_ENERGY) += electEnergy;
      //reduction->item(REDUCTION_LJ_ENERGY) += kernel_data->vdwEnergy;
      //reduction->item(REDUCTION_ELECT_ENERGY) += kernel_data->electEnergy;
    }
    if (doSlow) {
      Tensor virial_slow_tensor;
      virial_slow_tensor.xx = fullElectVirial_xx;
      virial_slow_tensor.xy = fullElectVirial_xy;
      virial_slow_tensor.xz = fullElectVirial_xz;
      virial_slow_tensor.yx = fullElectVirial_xy;
      virial_slow_tensor.yy = fullElectVirial_yy;
      virial_slow_tensor.yz = fullElectVirial_yz;
      virial_slow_tensor.zx = fullElectVirial_xz;
      virial_slow_tensor.zy = fullElectVirial_yz;
      virial_slow_tensor.zz = fullElectVirial_zz;
      //virial_slow_tensor.xx = kernel_data->fullElectVirial_xx;
      //virial_slow_tensor.xy = kernel_data->fullElectVirial_xy;
      //virial_slow_tensor.xz = kernel_data->fullElectVirial_xz;
      //virial_slow_tensor.yx = kernel_data->fullElectVirial_xy;
      //virial_slow_tensor.yy = kernel_data->fullElectVirial_yy;
      //virial_slow_tensor.yz = kernel_data->fullElectVirial_yz;
      //virial_slow_tensor.zx = kernel_data->fullElectVirial_xz;
      //virial_slow_tensor.zy = kernel_data->fullElectVirial_yz;
      //virial_slow_tensor.zz = kernel_data->fullElectVirial_zz;
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_SLOW, virial_slow_tensor);
      if (doEnergy) { reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += fullElectEnergy; }
      //if (doEnergy) { reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += kernel_data->fullElectEnergy; }
    }

    // DMK - TRACING
    #if (MIC_TRACING != 0) && (MIC_DEVICE_TRACING != 0)
    {
      #if MIC_DEVICE_TRACING_DETAILED != 0

        // Create compute user events
        ComputeMap *computeMap = ComputeMap::Object();
        PatchMap *patchMap = PatchMap::Object();
        int aSize = patchMap->gridsize_a();
        int bSize = patchMap->gridsize_b();
        int cSize = patchMap->gridsize_c();
        for (int i = 0; i < device__patch_pairs_size; i++) {

          // Determine the distance between the patches
          int pid0 = device__patch_pairs[i].p1;
          int pid1 = device__patch_pairs[i].p2;
          int dist = 0;
          if (pid0 != pid1) { // Is a pair, not a self
            int trans0 = computeMap->trans(device__patch_pairs[i].cid, 0);
            int trans1 = computeMap->trans(device__patch_pairs[i].cid, 1);
            int index_a0 = patchMap->index_a(pid0) + aSize * Lattice::offset_a(trans0);
            int index_b0 = patchMap->index_b(pid0) + bSize * Lattice::offset_b(trans0);
            int index_c0 = patchMap->index_c(pid0) + cSize * Lattice::offset_c(trans0);
            int index_a1 = patchMap->index_a(pid1) + aSize * Lattice::offset_a(trans1);
            int index_b1 = patchMap->index_b(pid1) + bSize * Lattice::offset_b(trans1);
            int index_c1 = patchMap->index_c(pid1) + cSize * Lattice::offset_c(trans1);
            int da = index_a0 - index_a1; da *= ((da < 0) ? (-1) : (1));
            int db = index_b0 - index_b1; db *= ((db < 0) ? (-1) : (1));
            int dc = index_c0 - index_c1; dc *= ((dc < 0) ? (-1) : (1));
            dist = da + db + dc;
	  }

          // Retrieve the start and end times for the "compute's task"
          double ueStart = device__device_times_computes[i * 2];
          double ueEnd = device__device_times_computes[i * 2 + 1];
          #if MIC_MULTIPLE_KERNELS != 0
            ueStart += mic_tracing_offload_start_remote - device__device_times_start[0];
            ueEnd += mic_tracing_offload_start_remote - device__device_times_start[0];
          #else
            ueStart += mic_tracing_offload_start_local - device__device_times_start[1];
            ueEnd += mic_tracing_offload_start_local - device__device_times_start[1];
          #endif

          //if (i == 0) { printf("[DEBUG:MIC:TRACING:COMPUTE] user event = { %.6le, %.6le } - %.6le\n", ueStart, ueEnd, ueEnd - ueStart); }
          if (dist > 7) { dist = 7; }  // NOTE: Make sure that the distance is put in the "7+" category if it is >7
          traceUserBracketEvent(MIC_EVENT_DEVICE_COMPUTE + dist, ueStart, ueEnd);
        }

        // Create patch (force reduction) user events
        for (int i = 0; i < device__force_lists_size; i++) {
          double ueStart = device__device_times_patches[i * 2];
          double ueEnd = device__device_times_patches[i * 2 + 1];
          #if MIC_MULTIPLE_KERNELS != 0
            ueStart += mic_tracing_offload_start_remote - device__device_times_start[0];
            ueEnd += mic_tracing_offload_start_remote - device__device_times_start[0];
          #else
            ueStart += mic_tracing_offload_start_local - device__device_times_start[1];
            ueEnd += mic_tracing_offload_start_local - device__device_times_start[1];
          #endif
	  //if (i == 0) { printf("[DEBUG:MIC:TRACING:PATCH] user event = { %.6le, %.6le } - %.6le\n", ueStart, ueEnd, ueEnd - ueStart); }
          traceUserBracketEvent(MIC_EVENT_DEVICE_PATCH, ueStart, ueEnd);
        }

      #endif // MIC_DEVICE_TRACING_DETAILED

      // Create phases
      double lPTime0 = device__device_times_start[3];
      double lPTime1 = device__device_times_start[5];
      double lPTime2 = device__device_times_start[7];
      double lPTime3 = device__device_times_start[9];
      #if MIC_MULTIPLE_KERNELS != 0
        double rPTime0 = device__device_times_start[2];
        double rPTime1 = device__device_times_start[4];
        double rPTime2 = device__device_times_start[6];
        double rPTime3 = device__device_times_start[8];
        rPTime0 += mic_tracing_offload_start_remote - device__device_times_start[0];
        lPTime0 += mic_tracing_offload_start_remote - device__device_times_start[0];
        rPTime1 += mic_tracing_offload_start_remote - device__device_times_start[0];
        lPTime1 += mic_tracing_offload_start_remote - device__device_times_start[0];
        rPTime2 += mic_tracing_offload_start_remote - device__device_times_start[0];
        lPTime2 += mic_tracing_offload_start_remote - device__device_times_start[0];
        rPTime3 += mic_tracing_offload_start_remote - device__device_times_start[0];
        lPTime3 += mic_tracing_offload_start_remote - device__device_times_start[0];
        traceUserBracketEvent(MIC_EVENT_DEVICE_COMPUTES, rPTime0, rPTime1);
        traceUserBracketEvent(MIC_EVENT_DEVICE_VIRIALS, rPTime1, rPTime2);
        traceUserBracketEvent(MIC_EVENT_DEVICE_PATCHES, rPTime2, rPTime3);
      #else
        lPTime0 += mic_tracing_offload_start_local - device__device_times_start[1];
        lPTime1 += mic_tracing_offload_start_local - device__device_times_start[1];
        lPTime2 += mic_tracing_offload_start_local - device__device_times_start[1];
        lPTime3 += mic_tracing_offload_start_local - device__device_times_start[1];
      #endif
      traceUserBracketEvent(MIC_EVENT_DEVICE_COMPUTES, lPTime0, lPTime1);
      traceUserBracketEvent(MIC_EVENT_DEVICE_VIRIALS, lPTime1, lPTime2);
      traceUserBracketEvent(MIC_EVENT_DEVICE_PATCHES, lPTime2, lPTime3);
    }
    #endif  // MIC_DEVICE_TRACING

    // DMK - MEMORY CHECK
    #if MIC_MALLOC_CHECK != 0
      _MM_MALLOC_VERIFY;
      #pragma offload target(mic:myDevice)
      { _MM_MALLOC_VERIFY; }
    #endif
  }

  if (workStarted == 1) { return 0; }
  //if ( workStarted == 1 && ! mergegrids &&
  //     ( localHostedPatches.size() || master == this ) ) {
  //  GBISP("not finished, call again\n");
  //  return 0;  // not finished, call again
  //}

  if ( master != this ) {  // finished
    GBISP("finished\n");
    if (simParams->GBISOn) gbisPhase = 1 + (gbisPhase % 3);//1->2->3->1...
    atomsChanged = 0;

    #if MIC_TRACING != 0
      MIC_TRACING_RECORD(MIC_EVENT_FUNC_FINISHWORK, finishWork_start, CmiWallTimer());
    #endif

    return 1;
  }

  mic_timer_total += kernel_time;

  // DMK - NOTE : If this is the end of the step for this compute
  if ( !simParams->GBISOn || gbisPhase == 3 ) {

    //// DMK - NOTE : Contribute virial
    //Tensor virial_tensor;
    //virial_tensor.xx = kernel_data->virial_xx;
    //virial_tensor.xy = kernel_data->virial_xy;
    //virial_tensor.xz = kernel_data->virial_xz;
    //virial_tensor.yx = kernel_data->virial_xy;
    //virial_tensor.yy = kernel_data->virial_yy;
    //virial_tensor.yz = kernel_data->virial_yz;
    //virial_tensor.zx = kernel_data->virial_xz;
    //virial_tensor.zy = kernel_data->virial_yz;
    //virial_tensor.zz = kernel_data->virial_zz;
    //// DMK - TODO | FIXME : GBIS support needed eventually
    //ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NBOND, virial_tensor);
    //if (doEnergy) {
    //  reduction->item(REDUCTION_LJ_ENERGY) += kernel_data->vdwEnergy;
    //  reduction->item(REDUCTION_ELECT_ENERGY) += kernel_data->electEnergy;
    //  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += kernel_data->fullElectEnergy;
    //}
    //if (doSlow) {
    //  Tensor virial_slow_tensor;
    //  virial_slow_tensor.xx = kernel_data->fullElectVirial_xx;
    //  virial_slow_tensor.xy = kernel_data->fullElectVirial_xy;
    //  virial_slow_tensor.xz = kernel_data->fullElectVirial_xz;
    //  virial_slow_tensor.yx = kernel_data->fullElectVirial_xy;
    //  virial_slow_tensor.yy = kernel_data->fullElectVirial_yy;
    //  virial_slow_tensor.yz = kernel_data->fullElectVirial_yz;
    //  virial_slow_tensor.zx = kernel_data->fullElectVirial_xz;
    //  virial_slow_tensor.zy = kernel_data->fullElectVirial_yz;
    //  virial_slow_tensor.zz = kernel_data->fullElectVirial_zz;
    //  ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_SLOW, virial_slow_tensor);
    //}

    atomsChanged = 0;
    reduction->submit();

    #if (VTUNE_INFO != 0) && (VTUNE_TASKS != 0)
      endOffloadTask();
    #endif

    mic_timer_count++;
    if ( simParams->outputCudaTiming &&
	  step % simParams->outputCudaTiming == 0 ) {

      // int natoms = mol->numAtoms; 
      // double wpa = wcount;  wpa /= natoms;

      // CkPrintf("Pe %d MIC kernel %f ms, total %f ms, wpa %f\n", CkMyPe(),
      //          kernel_time * 1.e3, time * 1.e3, wpa);

      #if 0
        float upload_ms, remote_calc_ms;
        float local_calc_ms, total_ms;
        mic_errcheck("before event timers");
        micEventElapsedTime(&upload_ms, start_upload, start_calc);
        mic_errcheck("in event timer 1");
        micEventElapsedTime(&remote_calc_ms, start_calc, end_remote_download);
        mic_errcheck("in event timer 2");
        micEventElapsedTime(&local_calc_ms, end_remote_download, end_local_download);
        mic_errcheck("in event timer 3");
        micEventElapsedTime(&total_ms, start_upload, end_local_download);
        mic_errcheck("in event timer 4");
        mic_errcheck("in event timers");

        CkPrintf("MIC EVENT TIMING: %d %f %f %f %f\n",
                 CkMyPe(), upload_ms, remote_calc_ms,
                 local_calc_ms, total_ms);
      #endif

      if ( mic_timer_count >= simParams->outputCudaTiming ) {
        mic_timer_total /= mic_timer_count;
        CkPrintf("MIC TIMING: %d  %f ms/step on node %d\n",
                 step, mic_timer_total * 1.e3, CkMyPe());
      }
      mic_timer_count = 0;
      mic_timer_total = 0;
    }

  } // !GBISOn || gbisPhase==3  

  // Next GBIS Phase
  GBISP("C.N.MIC[%d]::fnWork: incrementing phase\n", CkMyPe())
  if (simParams->GBISOn) gbisPhase = 1 + (gbisPhase % 3);//1->2->3->1...

  GBISP("C.N.MIC[%d] finished ready for next step\n",CkMyPe());

  #if MIC_TRACING != 0
    MIC_TRACING_RECORD(MIC_EVENT_FUNC_FINISHWORK, finishWork_start, CmiWallTimer());
  #endif

  return 1;  // finished and ready for next step
}

void mic_print_stats() { _mic_print_stats(myDevice); }


#endif  // NAMD_MIC
