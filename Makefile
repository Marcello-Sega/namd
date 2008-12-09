# pass version/platform information to compile
NAMD_VERSION = 2.6

default: all

# compiler flags (Win32 overrides)
COPTI = -I
COPTC = -c
COPTD = -D
COPTO = -o $(SPACE)

include Make.config

# pass version/platform information to compile
RELEASE=$(COPTD)NAMD_VERSION=\"$(NAMD_VERSION)\" $(COPTD)NAMD_PLATFORM=\"$(NAMD_PLATFORM)\" $(SCYLDFLAGS)

# directories
SRCDIR = src
DSTDIR = obj
INCDIR = inc
DPMTADIR=dpmta-2.6
DPMEDIR=dpme2
PLUGINSRCDIR= plugins/molfile_plugin/src
PLUGININCDIR= plugins/include
SBSRCDIR = sb/src

# comment/uncomment these lines for (D)PMTA routines
#DPMTAINCL=$(COPTI)$(DPMTADIR)/mpole $(COPTI)$(DPMTADIR)/src
#DPMTALIB=-L$(DPMTADIR)/mpole -L$(DPMTADIR)/src -ldpmta2 -lmpole -lpvmc
#DPMTAFLAGS=$(COPTD)DPMTA
#DPMTA=$(DPMTAINCL) $(DPMTAFLAGS)
#DPMTALIBS=$(DPMTADIR)/mpole/libmpole.a $(DPMTADIR)/src/libdpmta2.a

# comment/uncomment these lines for DPME routines
#DPMEINCL=$(COPTI)$(DPMEDIR)
#DPMELIB=-L$(DPMEDIR) -ldpme
#DPMEFLAGS=$(COPTD)DPME
#DPME=$(DPMEINCL) $(DPMEFLAGS)
#DPMELIBS= $(DPMEDIR)/libdpme.a

# to compile a memory optimized version, add -DMEM_OPT_VERSION
# to compile version that uses node aware spanning tree, add -DNODEAWARE_PROXY_SPANNINGTREE
# to the variable EXTRADEFINES
#EXTRADEFINES=-DMEM_OPT_VERSION
#EXTRADEFINES=-DREMOVE_PROXYDATAMSG_EXTRACOPY -DREMOVE_PROXYRESULTMSG_EXTRACOPY
EXTRADEFINES=-DREMOVE_PROXYRESULTMSG_EXTRACOPY 

# defaults for special cases
CXXTHREADOPTS = $(CXXOPTS) 
CXXSIMPARAMOPTS = $(CXXOPTS) 
CXXNOALIASOPTS = $(CXXOPTS) 
CUDACC = $(CXX)
CUDAOBJS =

include Make.config

# Add new source files here.

OBJS = \
	$(DSTDIR)/common.o \
	$(DSTDIR)/dcdlib.o \
	$(DSTDIR)/erf.o \
	$(DSTDIR)/fitrms.o \
	$(DSTDIR)/main.o \
	$(DSTDIR)/mainfunc.o \
	$(DSTDIR)/memusage.o \
	$(DSTDIR)/strlib.o \
	$(DSTDIR)/AlgSeven.o \
	$(DSTDIR)/AlgRecBisection.o \
	$(DSTDIR)/AlgNbor.o \
	$(DSTDIR)/AtomMap.o \
	$(DSTDIR)/BackEnd.o \
	$(DSTDIR)/BroadcastMgr.o \
	$(DSTDIR)/BroadcastClient.o \
	$(DSTDIR)/CollectionMaster.o \
	$(DSTDIR)/CollectionMgr.o \
	$(DSTDIR)/Communicate.o \
	$(DSTDIR)/Compute.o \
	$(DSTDIR)/ComputeAngles.o \
	$(DSTDIR)/ComputeBonds.o \
	$(DSTDIR)/ComputeConsForce.o \
	$(DSTDIR)/ComputeConsForceMsgs.o \
	$(DSTDIR)/ComputeCrossterms.o \
	$(DSTDIR)/ComputeCylindricalBC.o \
	$(DSTDIR)/ComputeDihedrals.o \
	$(DSTDIR)/ComputeDPME.o \
	$(DSTDIR)/ComputeDPMEMsgs.o \
	$(DSTDIR)/ComputeDPMTA.o \
	$(DSTDIR)/ComputeEField.o \
	$(DSTDIR)/ComputeEwald.o \
	$(DSTDIR)/ComputeExt.o \
	$(DSTDIR)/ComputeFullDirect.o \
	$(DSTDIR)/ComputeHomePatch.o \
	$(DSTDIR)/ComputeHomePatches.o \
	$(DSTDIR)/ComputeImpropers.o \
	$(DSTDIR)/ComputeGlobal.o \
	$(DSTDIR)/ComputeGlobalMsgs.o \
	$(DSTDIR)/ComputeGridForce.o \
	$(DSTDIR)/ComputeMap.o \
	$(DSTDIR)/ComputeMgr.o \
	$(DSTDIR)/ComputeNonbondedSelf.o \
	$(DSTDIR)/ComputeNonbondedPair.o \
	$(DSTDIR)/ComputeNonbondedUtil.o \
	$(DSTDIR)/ComputeNonbondedStd.o \
	$(DSTDIR)/ComputeNonbondedFEP.o \
	$(DSTDIR)/ComputeNonbondedLES.o \
	$(DSTDIR)/ComputeNonbondedPProf.o \
	$(DSTDIR)/ComputeNonbondedCUDA.o \
	$(DSTDIR)/ComputePatch.o \
	$(DSTDIR)/ComputePatchPair.o \
	$(DSTDIR)/ComputePme.o \
	$(DSTDIR)/ComputeRestraints.o \
	$(DSTDIR)/ComputeSphericalBC.o \
	$(DSTDIR)/ComputeStir.o \
	$(DSTDIR)/ComputeTclBC.o \
	$(DSTDIR)/ConfigList.o \
	$(DSTDIR)/Controller.o \
	$(DSTDIR)/ccsinterface.o \
	$(DSTDIR)/DataStream.o \
	$(DSTDIR)/DumpBench.o \
        $(DSTDIR)/FreeEnergyAssert.o \
        $(DSTDIR)/FreeEnergyGroup.o \
        $(DSTDIR)/FreeEnergyLambda.o \
        $(DSTDIR)/FreeEnergyLambdMgr.o \
        $(DSTDIR)/FreeEnergyParse.o \
        $(DSTDIR)/FreeEnergyRestrain.o \
        $(DSTDIR)/FreeEnergyRMgr.o \
        $(DSTDIR)/FreeEnergyVector.o \
	$(DSTDIR)/GlobalMaster.o \
	$(DSTDIR)/GlobalMasterServer.o \
	$(DSTDIR)/GlobalMasterTest.o \
	$(DSTDIR)/GlobalMasterIMD.o \
	$(DSTDIR)/GlobalMasterTcl.o \
	$(DSTDIR)/GlobalMasterSMD.o \
	$(DSTDIR)/GlobalMasterTMD.o \
	$(DSTDIR)/GlobalMasterFreeEnergy.o \
	$(DSTDIR)/GlobalMasterEasy.o \
	$(DSTDIR)/GlobalMasterMisc.o \
	$(DSTDIR)/colvarmodule.o \
	$(DSTDIR)/colvarparse.o \
	$(DSTDIR)/colvar.o \
	$(DSTDIR)/colvarvalue.o \
	$(DSTDIR)/colvarbias.o \
	$(DSTDIR)/colvarbias_abf.o \
	$(DSTDIR)/colvarbias_meta.o \
	$(DSTDIR)/colvaratoms.o \
	$(DSTDIR)/colvarcomp.o \
	$(DSTDIR)/colvarcomp_angles.o \
	$(DSTDIR)/colvarcomp_coordnums.o \
	$(DSTDIR)/colvarcomp_distances.o \
	$(DSTDIR)/colvarcomp_protein.o \
	$(DSTDIR)/colvarcomp_rotations.o \
	$(DSTDIR)/colvarproxy_namd.o \
	$(DSTDIR)/GridForceGrid.o \
        $(DSTDIR)/GromacsTopFile.o \
	$(DSTDIR)/heap.o \
	$(DSTDIR)/HomePatch.o \
	$(DSTDIR)/IMDOutput.o \
	$(DSTDIR)/InfoStream.o \
	$(DSTDIR)/LdbCoordinator.o \
	$(DSTDIR)/LJTable.o \
	$(DSTDIR)/Measure.o \
	$(DSTDIR)/MGridforceParams.o \
	$(DSTDIR)/MStream.o \
	$(DSTDIR)/MigrateAtomsMsg.o \
	$(DSTDIR)/Molecule.o \
	$(DSTDIR)/NamdCentLB.o \
	$(DSTDIR)/NamdNborLB.o \
	$(DSTDIR)/NamdState.o \
	$(DSTDIR)/NamdOneTools.o \
	$(DSTDIR)/Node.o \
	$(DSTDIR)/Output.o \
	$(DSTDIR)/Parameters.o \
	$(DSTDIR)/ParseOptions.o \
	$(DSTDIR)/Patch.o \
	$(DSTDIR)/PatchMgr.o \
	$(DSTDIR)/PatchMap.o \
	$(DSTDIR)/PDB.o \
	$(DSTDIR)/PDBData.o \
	$(DSTDIR)/PmeBase.o \
	$(DSTDIR)/PmeKSpace.o \
	$(DSTDIR)/PmeRealSpace.o \
	$(DSTDIR)/ProcessorPrivate.o \
	$(DSTDIR)/ProxyMgr.o \
	$(DSTDIR)/ProxyPatch.o \
	$(DSTDIR)/Rebalancer.o \
	$(DSTDIR)/RecBisection.o \
	$(DSTDIR)/ReductionMgr.o \
	$(DSTDIR)/RefineOnly.o \
	$(DSTDIR)/RefineTorusLB.o \
	$(DSTDIR)/ScriptTcl.o \
	$(DSTDIR)/Sequencer.o \
	$(DSTDIR)/Set.o \
	$(DSTDIR)/Settle.o \
	$(DSTDIR)/SimParameters.o \
	$(DSTDIR)/Sync.o \
	$(DSTDIR)/TclCommands.o \
	$(DSTDIR)/TorusLB.o \
	$(DSTDIR)/WorkDistrib.o \
	$(DSTDIR)/pub3dfft.o \
	$(DSTDIR)/vmdsock.o \
	$(DSTDIR)/parm.o \
	$(DSTDIR)/imd.o \
	$(DSTDIR)/CompressPsf.o \
	$(DSTDIR)/AtomsDisInfo.o \
	$(DSTDIR)/FileIO.o 


# Add new modules here and also define explicit rule below.

CIFILES = 	\
		$(INCDIR)/BroadcastMgr.decl.h \
		$(INCDIR)/BroadcastMgr.def.h \
		$(INCDIR)/CollectionMaster.decl.h \
		$(INCDIR)/CollectionMaster.def.h \
		$(INCDIR)/CollectionMgr.decl.h \
		$(INCDIR)/CollectionMgr.def.h \
		$(INCDIR)/ComputeMgr.decl.h \
		$(INCDIR)/ComputeMgr.def.h \
		$(INCDIR)/ComputePmeMgr.decl.h \
		$(INCDIR)/ComputePmeMgr.def.h \
		$(INCDIR)/ComputeExtMgr.decl.h \
		$(INCDIR)/ComputeExtMgr.def.h \
		$(INCDIR)/LdbCoordinator.decl.h \
		$(INCDIR)/LdbCoordinator.def.h \
		$(INCDIR)/NamdCentLB.decl.h \
		$(INCDIR)/NamdCentLB.def.h \
		$(INCDIR)/NamdNborLB.decl.h \
		$(INCDIR)/NamdNborLB.def.h \
		$(INCDIR)/Node.decl.h \
		$(INCDIR)/Node.def.h \
		$(INCDIR)/PatchMgr.decl.h \
		$(INCDIR)/PatchMgr.def.h \
		$(INCDIR)/ProxyMgr.decl.h \
		$(INCDIR)/ProxyMgr.def.h \
		$(INCDIR)/ReductionMgr.decl.h \
		$(INCDIR)/ReductionMgr.def.h \
		$(INCDIR)/Sync.decl.h \
		$(INCDIR)/Sync.def.h \
		$(INCDIR)/WorkDistrib.decl.h \
		$(INCDIR)/WorkDistrib.def.h \
		$(INCDIR)/main.decl.h \
		$(INCDIR)/main.def.h \
		$(INCDIR)/AtomsDisInfo.decl.h \
		$(INCDIR)/AtomsDisInfo.def.h \
		$(INCDIR)/FileIO.decl.h \
		$(INCDIR)/FileIO.def.h

# Add new source files here.

PLUGINOBJS = \
	$(DSTDIR)/dcdplugin.o \
	$(DSTDIR)/jsplugin.o \
	$(DSTDIR)/pdbplugin.o \
	$(DSTDIR)/psfplugin.o

PLUGINLIB = $(PLUGINOBJS)

CUDAOBJSRAW = \
	$(DSTDIR)/ComputeNonbondedCUDAKernel.o

$(DSTDIR)/ComputeNonbondedCUDAKernel.o: \
	$(SRCDIR)/ComputeNonbondedCUDAKernel.cu \
	$(SRCDIR)/ComputeNonbondedCUDAKernel.h
	$(CUDACC) $(COPTO)$(DSTDIR)/ComputeNonbondedCUDAKernel.o $(COPTC) $(SRCDIR)/ComputeNonbondedCUDAKernel.cu
	$(CUDACC) -ptx $(SRCDIR)/ComputeNonbondedCUDAKernel.cu
	grep global ComputeNonbondedCUDAKernel.ptx
	$(CUDACC) -cubin $(SRCDIR)/ComputeNonbondedCUDAKernel.cu
	grep reg ComputeNonbondedCUDAKernel.cubin

SBOBJS = \
	$(DSTDIR)/tcl_main.o \
	$(DSTDIR)/tcl_psfgen.o \
	$(DSTDIR)/charmm_file.o \
	$(DSTDIR)/charmm_parse_topo_defs.o \
	$(DSTDIR)/extract_alias.o \
	$(DSTDIR)/hash.o \
	$(DSTDIR)/hasharray.o \
	$(DSTDIR)/memarena.o \
	$(DSTDIR)/pdb_file.o \
	$(DSTDIR)/pdb_file_extract.o \
	$(DSTDIR)/psf_file.o \
	$(DSTDIR)/psf_file_extract.o \
	$(DSTDIR)/topo_defs.o \
	$(DSTDIR)/topo_mol.o \
	$(DSTDIR)/topo_mol_output.o \
	$(DSTDIR)/stringhash.o

# definitions for Charm routines
CHARMC = $(CHARM)/bin/charmc
CHARMXI = $(CHARM)/bin/charmc
CHARMINC = $(CHARM)/include $(COPTD)CMK_OPTIMIZE=1
CHARMLIB = $(CHARM)/lib

# Libraries we may have changed
LIBS = $(CUDAOBJS) $(PLUGINLIB) $(DPMTALIBS) $(DPMELIBS) $(TCLDLL)

# CXX is platform dependent
CXXBASEFLAGS = $(COPTI)$(CHARMINC) $(COPTI)$(SRCDIR) $(COPTI)$(INCDIR) $(DPMTA) $(DPME) $(COPTI)$(PLUGININCDIR) $(TCL) $(FFT) $(CUDA) $(CCS) $(RELEASE) $(EXTRADEFINES) $(TRACEOBJDEF)
CXXFLAGS = $(CXXBASEFLAGS) $(CXXOPTS)
CXXTHREADFLAGS = $(CXXBASEFLAGS) $(CXXTHREADOPTS)
CXXSIMPARAMFLAGS = $(CXXBASEFLAGS) $(CXXSIMPARAMOPTS)
CXXNOALIASFLAGS = $(CXXBASEFLAGS) $(CXXNOALIASOPTS)
GXXFLAGS = $(CXXBASEFLAGS) -DNO_STRSTREAM_H
CFLAGS = $(COPTI)$(SRCDIR) $(TCL) $(COPTS) $(RELEASE) $(EXTRADEFINES) $(TRACEOBJDEF)
PLUGINGCCFLAGS = $(COPTI)$(PLUGINSRCDIR) $(COPTI)$(PLUGININCDIR)
PLUGINCFLAGS = $(PLUGINGCCFLAGS) $(COPTS)
SBCFLAGS = $(COPTI)$(SBSRCDIR) $(TCL) $(COPTS) $(RELEASE) $(EXTRADEFINES) $(TRACEOBJDEF)
SBGCCFLAGS = $(COPTI)$(SBSRCDIR) $(TCL) $(RELEASE) $(EXTRADEFINES) $(TRACEOBJDEF)

# Add new executables here.

BINARIES = namd2 psfgen charmrun flipdcd flipbinpdb

# This should be rebuilt at every compile, but not on Win32.
BUILDINFO = $(DSTDIR)/buildinfo
MAKEBUILDINFO = \
	$(RM) $(BUILDINFO).C; \
	echo 'const char *namd_build_date = ' \"`date`\"\; > $(BUILDINFO).C; \
	echo 'const char *namd_build_user = ' \"$(USER)\"\; >> $(BUILDINFO).C; \
	echo 'const char *namd_build_machine = ' \"`hostname`\"\; >> $(BUILDINFO).C; \
	cat $(BUILDINFO).C; \
	$(CXX) $(CXXFLAGS) $(COPTO)$(BUILDINFO).o $(COPTC) $(BUILDINFO).C

all:	$(BINARIES)

namd2:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose -ld++-option \
	"$(COPTI)$(CHARMINC) $(COPTI)$(INCDIR) $(COPTI)$(SRCDIR) $(CXXOPTS)" \
	-module NeighborLB -module commlib -language charm++ \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	-lm -o namd2

charmrun: $(CHARM)/bin/charmrun # XXX
	$(COPY) $(CHARM)/bin/charmrun $@

windowsbinaries: namd2.exe psfgen.exe charmd.exe charmd_faceless.exe charmrun.exe

namd2.exe:  $(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose \
	-module NeighborLB -module commlib -language charm++ \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	-o namd2

charmd.exe:
	$(COPY) $(CHARM)/bin/charmd.exe charmd.exe

charmd_faceless.exe:
	$(COPY) $(CHARM)/bin/charmd_faceless.exe charmd_faceless.exe

charmrun.exe:
	$(COPY) $(CHARM)/bin/charmrun.exe charmrun.exe

psfgen:	$(DSTDIR) $(SBOBJS)
	$(CC) $(SBCFLAGS) -o psfgen $(SBOBJS) $(TCLLIB) $(TCLAPPLIB) -lm

psfgen.exe:	$(DSTDIR) $(SBOBJS)
	$(LINK) $(LINKOPTS) /out:psfgen.exe $(SBOBJS) $(TCLWINLIB) $(TCLAPPLIB)

flipdcd:	$(SRCDIR)/flipdcd.c
	$(CC) $(CFLAGS) -o $@ $(SRCDIR)/flipdcd.c || \
	echo "#!/bin/sh\necho unavailable on this platform" > $@; \
	chmod +x $@

flipbinpdb:	$(SRCDIR)/flipbinpdb.c
	$(CC) $(CFLAGS) -o $@ $(SRCDIR)/flipbinpdb.c || \
	echo "#!/bin/sh\necho unavailable on this platform" > $@; \
	chmod +x $@

fixdcd:	$(SRCDIR)/fixdcd.c
	$(CC) $(CFLAGS) -o fixdcd $(SRCDIR)/fixdcd.c

dumpdcd:	$(SRCDIR)/dumpdcd.c
	$(CC) $(CFLAGS) -o dumpdcd $(SRCDIR)/dumpdcd.c

loaddcd:	$(SRCDIR)/loaddcd.c
	$(CC) $(CFLAGS) -o loaddcd $(SRCDIR)/loaddcd.c

updatefiles:
	touch ../src/ComputeSelfTuples.h
	rm -f obj/ComputeNonbondedPair.o
	rm -f obj/ComputeNonbondedSelf.o
	rm -f obj/ComputePme.o

#To compile tracecomputes, type the command "make tracecomputes TRACEOBJDEF=-DTRACE_COMPUTE_OBJECTS"
tracecomputes: updatefiles $(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose -ld++-option \
	"$(COPTI)$(CHARMINC) $(COPTI)$(INCDIR) $(COPTI)$(SRCDIR) $(CXXOPTS)" \
	-module NeighborLB -module commlib -language charm++ \
	-tracemode projections -tracemode summary \
	$(BUILDINFO).o \
	$(OBJS) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	-lm -o namd2

projections: $(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose -ld++-option \
	"$(COPTI)$(CHARMINC) $(COPTI)$(INCDIR) $(COPTI)$(SRCDIR) $(CXXOPTS)" \
	-module NeighborLB -module commlib -language charm++ \
	-tracemode projections -tracemode summary \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	-lm -o namd2

summary: $(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose -ld++-option \
	"$(COPTI)$(CHARMINC) $(COPTI)$(INCDIR) $(COPTI)$(SRCDIR) $(CXXOPTS)" \
	-module NeighborLB -module commlib -language charm++ \
	-tracemode summary \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	-lm -o namd2

$(DPMTADIR)/mpole/libmpole.a: $(DPMTADIR)/src/libdpmta2.a

$(DPMTADIR)/src/libdpmta2.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMEDIR)/libdpme.a:
	cd $(DPMEDIR) ; $(MAKE) ; cd ..

# Unix commands

ECHO = echo
MOVE = mv
COPY = cp
RM = rm -f
LDD = ldd

include Make.config

# Explicit rules for modules.

$(INCDIR)/BroadcastMgr.def.h: $(INCDIR)/BroadcastMgr.decl.h

$(INCDIR)/BroadcastMgr.decl.h: $(SRCDIR)/BroadcastMgr.ci
	$(CHARMXI) $(SRCDIR)/BroadcastMgr.ci
	$(MOVE) BroadcastMgr.def.h $(INCDIR)
	$(MOVE) BroadcastMgr.decl.h $(INCDIR)

$(INCDIR)/CollectionMaster.def.h: $(INCDIR)/CollectionMaster.decl.h

$(INCDIR)/CollectionMaster.decl.h: $(SRCDIR)/CollectionMaster.ci
	$(CHARMXI) $(SRCDIR)/CollectionMaster.ci
	$(MOVE) CollectionMaster.def.h $(INCDIR)
	$(MOVE) CollectionMaster.decl.h $(INCDIR)

$(INCDIR)/CollectionMgr.def.h: $(INCDIR)/CollectionMgr.decl.h

$(INCDIR)/CollectionMgr.decl.h: $(SRCDIR)/CollectionMgr.ci
	$(CHARMXI) $(SRCDIR)/CollectionMgr.ci
	$(MOVE) CollectionMgr.def.h $(INCDIR)
	$(MOVE) CollectionMgr.decl.h $(INCDIR)

$(INCDIR)/ComputeMgr.def.h: $(INCDIR)/ComputeMgr.decl.h

$(INCDIR)/ComputeMgr.decl.h: $(SRCDIR)/ComputeMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputeMgr.ci
	$(MOVE) ComputeMgr.def.h $(INCDIR)
	$(MOVE) ComputeMgr.decl.h $(INCDIR)

$(INCDIR)/ComputePmeMgr.def.h: $(INCDIR)/ComputePmeMgr.decl.h

$(INCDIR)/ComputePmeMgr.decl.h: $(SRCDIR)/ComputePmeMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputePmeMgr.ci
	$(MOVE) ComputePmeMgr.def.h $(INCDIR)
	$(MOVE) ComputePmeMgr.decl.h $(INCDIR)

$(INCDIR)/ComputeExtMgr.def.h: $(INCDIR)/ComputeExtMgr.decl.h

$(INCDIR)/ComputeExtMgr.decl.h: $(SRCDIR)/ComputeExtMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputeExtMgr.ci
	$(MOVE) ComputeExtMgr.def.h $(INCDIR)
	$(MOVE) ComputeExtMgr.decl.h $(INCDIR)

$(INCDIR)/LdbCoordinator.def.h: $(INCDIR)/LdbCoordinator.decl.h

$(INCDIR)/LdbCoordinator.decl.h: $(SRCDIR)/LdbCoordinator.ci
	$(CHARMXI) $(SRCDIR)/LdbCoordinator.ci
	$(MOVE) LdbCoordinator.def.h $(INCDIR)
	$(MOVE) LdbCoordinator.decl.h $(INCDIR)

$(INCDIR)/NamdCentLB.def.h: $(INCDIR)/NamdCentLB.decl.h

$(INCDIR)/NamdCentLB.decl.h: $(SRCDIR)/NamdCentLB.ci
	$(CHARMXI) $(SRCDIR)/NamdCentLB.ci
	$(MOVE) NamdCentLB.def.h $(INCDIR)
	$(MOVE) NamdCentLB.decl.h $(INCDIR)

$(INCDIR)/NamdNborLB.def.h: $(INCDIR)/NamdNborLB.decl.h

$(INCDIR)/NamdNborLB.decl.h: $(SRCDIR)/NamdNborLB.ci
	$(CHARMXI) $(SRCDIR)/NamdNborLB.ci
	$(MOVE) NamdNborLB.def.h $(INCDIR)
	$(MOVE) NamdNborLB.decl.h $(INCDIR)

$(INCDIR)/Node.def.h: $(INCDIR)/Node.decl.h

$(INCDIR)/Node.decl.h: $(SRCDIR)/Node.ci
	$(CHARMXI) $(SRCDIR)/Node.ci
	$(MOVE) Node.def.h $(INCDIR)
	$(MOVE) Node.decl.h $(INCDIR)

$(INCDIR)/PatchMgr.def.h: $(INCDIR)/PatchMgr.decl.h

$(INCDIR)/PatchMgr.decl.h: $(SRCDIR)/PatchMgr.ci
	$(CHARMXI) $(SRCDIR)/PatchMgr.ci
	$(MOVE) PatchMgr.def.h $(INCDIR)
	$(MOVE) PatchMgr.decl.h $(INCDIR)

$(INCDIR)/ProxyMgr.def.h: $(INCDIR)/ProxyMgr.decl.h

$(INCDIR)/ProxyMgr.decl.h: $(SRCDIR)/ProxyMgr.ci
	$(CHARMXI) $(SRCDIR)/ProxyMgr.ci
	$(MOVE) ProxyMgr.def.h $(INCDIR)
	$(MOVE) ProxyMgr.decl.h $(INCDIR)

$(INCDIR)/ReductionMgr.def.h: $(INCDIR)/ReductionMgr.decl.h

$(INCDIR)/ReductionMgr.decl.h: $(SRCDIR)/ReductionMgr.ci
	$(CHARMXI) $(SRCDIR)/ReductionMgr.ci
	$(MOVE) ReductionMgr.def.h $(INCDIR)
	$(MOVE) ReductionMgr.decl.h $(INCDIR)

$(INCDIR)/Sync.def.h: $(INCDIR)/Sync.decl.h

$(INCDIR)/Sync.decl.h: $(SRCDIR)/Sync.ci
	$(CHARMXI) $(SRCDIR)/Sync.ci
	$(MOVE) Sync.def.h $(INCDIR)
	$(MOVE) Sync.decl.h $(INCDIR)

$(INCDIR)/WorkDistrib.def.h: $(INCDIR)/WorkDistrib.decl.h

$(INCDIR)/WorkDistrib.decl.h: $(SRCDIR)/WorkDistrib.ci
	$(CHARMXI) $(SRCDIR)/WorkDistrib.ci
	$(MOVE) WorkDistrib.def.h $(INCDIR)
	$(MOVE) WorkDistrib.decl.h $(INCDIR)

$(INCDIR)/main.def.h: $(INCDIR)/main.decl.h

$(INCDIR)/main.decl.h: $(SRCDIR)/main.ci
	$(CHARMXI) $(SRCDIR)/main.ci
	$(MOVE) main.def.h $(INCDIR)
	$(MOVE) main.decl.h $(INCDIR)

$(INCDIR)/AtomsDisInfo.def.h: $(INCDIR)/AtomsDisInfo.decl.h

$(INCDIR)/AtomsDisInfo.decl.h: $(SRCDIR)/AtomsDisInfo.ci
	$(CHARMXI) $(SRCDIR)/AtomsDisInfo.ci
	$(MOVE) AtomsDisInfo.def.h $(INCDIR)
	$(MOVE) AtomsDisInfo.decl.h $(INCDIR)

$(INCDIR)/FileIO.def.h: $(INCDIR)/FileIO.decl.h

$(INCDIR)/FileIO.decl.h: $(SRCDIR)/FileIO.ci
	$(CHARMXI) $(SRCDIR)/FileIO.ci
	$(MOVE) FileIO.def.h $(INCDIR)
	$(MOVE) FileIO.decl.h $(INCDIR)


DEPENDFILE = .rootdir/Make.depends

# This is a CPU killer...  Don't make depends if you don't need to.
depends: $(INCDIR) $(CIFILES) $(DSTDIR) $(DEPENDFILE)
	$(ECHO) "Creating " $(DEPENDFILE) " ..."; \
	if [ -f $(DEPENDFILE) ]; then \
	   $(MOVE) -f $(DEPENDFILE) $(DEPENDFILE).old; \
	fi; \
	touch $(DEPENDFILE); \
	for i in $(OBJS) ; do \
	      SRCFILE=$(SRCDIR)/`basename $$i .o`.C ; \
	      $(ECHO) "checking dependencies for $$SRCFILE" ; \
	      g++ -MM $(GXXFLAGS) $$SRCFILE | \
	      perl $(SRCDIR)/dc.pl $(CHARMINC) $(TCLDIR) $(FFTDIR) /usr/include /usr/local >> $(DEPENDFILE); \
	      $(ECHO) '	$$(CXX) $$(CXXFLAGS) $$(COPTO)'$$i '$$(COPTC)' \
		$$SRCFILE >> $(DEPENDFILE) ; \
	done; \
	for i in $(PLUGINOBJS) ; do \
	      BASENAME=`basename $$i .o` ; \
	      SRCFILE=$(PLUGINSRCDIR)/$$BASENAME.c ; \
	      $(ECHO) "checking dependencies for $$SRCFILE" ; \
	      gcc -MM $(PLUGINGCCFLAGS) $$SRCFILE | \
	      perl $(SRCDIR)/dc.pl /usr/include /usr/local >> $(DEPENDFILE); \
	      $(ECHO) '	$$(CC) $$(PLUGINCFLAGS) $$(COPTO)'$$i '$$(COPTC)' \
		'$$(COPTD)'VMDPLUGIN=molfile_$$BASENAME \
		$$SRCFILE >> $(DEPENDFILE) ; \
	done; \
	for i in $(SBOBJS) ; do \
	      SRCFILE=$(SBSRCDIR)/`basename $$i .o`.c ; \
	      $(ECHO) "checking dependencies for $$SRCFILE" ; \
	      gcc -MM $(SBGCCFLAGS) $$SRCFILE | \
	      perl $(SRCDIR)/dc.pl $(CHARMINC) $(TCLDIR) $(FFTDIR) /usr/include /usr/local >> $(DEPENDFILE); \
	      $(ECHO) '	$$(CC) $$(SBCFLAGS) $$(COPTO)'$$i '$$(COPTC)' \
		$$SRCFILE >> $(DEPENDFILE) ; \
	done; \
	$(RM) $(DEPENDFILE).sed; \
	sed -e "/obj\/Controller.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/Sequencer.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/ComputeFullDirect.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/ReductionMgr.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/SimParameters.o/ s/CXXFLAGS/CXXSIMPARAMFLAGS/" \
	    -e "/obj\/ComputeNonbondedStd.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    -e "/obj\/ComputeNonbondedFEP.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    -e "/obj\/ComputeNonbondedLES.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    $(DEPENDFILE) > $(DEPENDFILE).sed; \
	$(MOVE) -f $(DEPENDFILE).sed $(DEPENDFILE);

$(DEPENDFILE):
	touch $(DEPENDFILE)

include Make.depends

$(DSTDIR):
	mkdir $(DSTDIR)

$(INCDIR):
	mkdir $(INCDIR)

clean:
	rm -rf ptrepository Templates.DB SunWS_cache $(DSTDIR) $(INCDIR)

veryclean:	clean
	rm -f $(BINARIES)

RELEASE_DIR_NAME = NAMD_$(NAMD_VERSION)_$(NAMD_PLATFORM)

DOC_FILES = README.txt announce.txt license.txt notes.txt

RELEASE_FILES = flipdcd flipbinpdb psfgen charmrun namd2

WIN32_RELEASE_FILES = namd2.exe psfgen.exe charmrun.exe charmd.exe charmd_faceless.exe $(TCLDLL)

release: all
	$(ECHO) Creating release $(RELEASE_DIR_NAME)
	mkdir $(RELEASE_DIR_NAME)
	cp $(RELEASE_FILES) $(RELEASE_DIR_NAME)
	for f in $(DOC_FILES); do cp .rootdir/$$f $(RELEASE_DIR_NAME); done
	cp -r .rootdir/lib $(RELEASE_DIR_NAME)
	/bin/rm -rf $(RELEASE_DIR_NAME)/lib/CVS $(RELEASE_DIR_NAME)/lib/*/CVS
	if [ -r $(CHARM)/bin/charmd ]; then \
	  $(COPY) $(CHARM)/bin/charmd $(RELEASE_DIR_NAME); \
	fi
	if [ -r $(CHARM)/bin/charmd_faceless ]; then \
	  $(COPY) $(CHARM)/bin/charmd_faceless $(RELEASE_DIR_NAME); \
	fi
	chmod -R a+rX $(RELEASE_DIR_NAME)
	tar cf $(RELEASE_DIR_NAME).tar $(RELEASE_DIR_NAME)
	gzip $(RELEASE_DIR_NAME).tar
	echo $(CHARM)
	ls -l $(CHARM)/lib
	-for f in $(RELEASE_FILES); do echo $$f; $(LDD) $$f; done

winrelease: winall
	$(ECHO) Creating release $(RELEASE_DIR_NAME)
	mkdir $(RELEASE_DIR_NAME)
	cp $(DOC_FILES) $(WIN32_RELEASE_FILES) $(RELEASE_DIR_NAME)
	cp -r lib $(RELEASE_DIR_NAME)
	chmod -R a+rX $(RELEASE_DIR_NAME)
	echo $(CHARM)
	ls -l $(CHARM)/lib
	echo $(CHARM)
	zip -r $(RELEASE_DIR_NAME).zip $(RELEASE_DIR_NAME)

