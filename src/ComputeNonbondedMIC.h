#include "ComputeNonbondedUtil.h"
#include "ComputeHomeTuples.h"
#include "ComputeNonbondedMICKernel.h"

class ComputeMgr;

class ComputeNonbondedMICKernel;

class float4;
class double4;

int mic_device_pe();

bool mic_device_shared_with_pe(int pe);

class ComputeNonbondedMIC : public Compute, private ComputeNonbondedUtil {
  public:

  struct compute_record {
    ComputeID c;
    PatchID pid[2];
    Vector offset;
    int isSelf;
    #if (MIC_ENABLE_COMPUTE_PARTITIONING != 0) || (MIC_ENABLE_MIC_SPECIFIC_COMPUTE_PARTITIONING != 0)
      int part;   // DMK - TODO : Look into making this a 'char' to reduce memory usage
      int numParts;
    #endif
  };

  struct patch_record {
    int localIndex;
    int localStart;
    int numAtoms;
    int numFreeAtoms;
    int refCount;
    int isLocal;
    int hostPe;
    PatchID patchID;
    Patch *p;
    Box<Patch,CompAtom> *positionBox;
    Box<Patch,Results> *forceBox;
    Box<Patch,Real>   *intRadBox; //5 GBIS Boxes
    Box<Patch,GBReal> *psiSumBox;
    Box<Patch,Real>   *bornRadBox;
    Box<Patch,GBReal> *dEdaSumBox;
    Box<Patch,Real>   *dHdrPrefixBox;
    CompAtom *x;
    CompAtomExt *xExt;
    Results *r;
    Force *f;
    Real   *intRad; //5 GBIS arrays
    GBReal *psiSum;
    Real   *bornRad;
    GBReal *dEdaSum;
    Real   *dHdrPrefix;

    patch_record() { refCount = 0; }
  };


    ComputeNonbondedMIC(ComputeID c, ComputeMgr *mgr,
		ComputeNonbondedMIC *m = 0, int idx = -1);
    ~ComputeNonbondedMIC();

    void atomUpdate();
    void doWork();
    int noWork();

    void recvYieldDevice(int pe);
    LocalWorkMsg *localWorkMsg2;

    int workStarted;
    Lattice lattice;
    int doSlow, doEnergy;
    int step;
    int finishWork();  // returns true when finished, false to continue
    void messageFinishWork();

    static void bind_lj_table(int deviceNum);
    static void bind_force_table(int deviceNum);
    static void bind_constants(int deviceNum);
    static void bind_exclusions(int deviceNum);

    void build_exclusions();

    void requirePatch(int pid);
    void assignPatches();
    void registerPatches();
    ResizeArray<int> activePatches, localActivePatches, remoteActivePatches;
    ResizeArray<int> hostedPatches, localHostedPatches, remoteHostedPatches;
    ResizeArray<patch_record> patchRecords;
    ResizeArray<compute_record> computeRecords;
    ResizeArray<compute_record> localComputeRecords, remoteComputeRecords;

    int num_atom_records;
    int num_local_atom_records;
    int num_remote_atom_records;
    int num_force_records;

    int num_local_compute_records;
    int num_remote_compute_records;
    int num_local_patch_records;
    int num_remote_patch_records;

    double4 *forces;
    double4 *slow_forces;
    GBReal *psiSumH;
    GBReal *dEdaSumH;

    PatchMap *patchMap;
    AtomMap *atomMap;
    SubmitReduction *reduction;

    ComputeNonbondedMICKernel *kernel;

    ComputeNonbondedMIC *master;
    int masterPe;
    int slaveIndex;
    ComputeNonbondedMIC **slaves;
    int *slavePes;
    int numSlaves;

    #if MIC_SUBMIT_ATOMS_ON_ARRIVAL != 0
      int micDevice;
      int2 *exclusionsByAtom_ptr;
      virtual void patchReady(PatchID patchID, int doneMigration, int seq);
    #endif

    // DMK - DEBUG
    int timestep;
};


