include Makearch

# pass version/platform information to compile
NAMD_VERSION = 2.1
RELEASE=-DNAMD_VERSION=\"$(NAMD_VERSION)\" -DNAMD_PLATFORM=\"$(NAMD_PLATFORM)\"

# directories
SRCDIR = src
DSTDIR = obj
INCDIR = inc
DPMTADIR=dpmta-2.6
DPMEDIR=dpme2

# comment/uncomment these lines for (D)PMTA routines
DPMTAINCL=-I$(DPMTADIR)/mpole -I$(DPMTADIR)/src
DPMTALIB=-L$(DPMTADIR)/mpole -L$(DPMTADIR)/src -ldpmta2 -lmpole -lpvmc
DPMTAFLAGS=-DDPMTA
DPMTA=$(DPMTAINCL) $(DPMTAFLAGS)
DPMTALIBS=$(DPMTADIR)/mpole/libmpole.a $(DPMTADIR)/src/libdpmta2.a

# comment/uncomment these lines for DPME routines
DPMEINCL=-I$(DPMEDIR)
DPMELIB=-L$(DPMEDIR) -ldpme
DPMEFLAGS=-DDPME
DPME=$(DPMEINCL) $(DPMEFLAGS)
DPMELIBS= $(DPMEDIR)/libdpme.a

# Add new source files here.

OBJS = \
	$(DSTDIR)/common.o \
	$(DSTDIR)/dcdlib.o \
	$(DSTDIR)/main.o \
	$(DSTDIR)/mainfunc.o \
	$(DSTDIR)/strlib.o \
	$(DSTDIR)/AlgSeven.o \
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
	$(DSTDIR)/ComputeCylindricalBC.o \
	$(DSTDIR)/ComputeDihedrals.o \
	$(DSTDIR)/ComputeDPME.o \
	$(DSTDIR)/ComputeDPMEMsgs.o \
	$(DSTDIR)/ComputeDPMTA.o \
	$(DSTDIR)/ComputeFreeEnergy.o \
	$(DSTDIR)/ComputeFullDirect.o \
	$(DSTDIR)/ComputeHomePatch.o \
	$(DSTDIR)/ComputeHomePatches.o \
	$(DSTDIR)/ComputeImpropers.o \
	$(DSTDIR)/ComputeGlobal.o \
	$(DSTDIR)/ComputeGlobalEasy.o \
	$(DSTDIR)/ComputeGlobalMaster.o \
	$(DSTDIR)/ComputeGlobalMsgs.o \
	$(DSTDIR)/ComputeIMD.o \
	$(DSTDIR)/ComputeMap.o \
	$(DSTDIR)/ComputeMisc.o \
	$(DSTDIR)/ComputeMgr.o \
	$(DSTDIR)/ComputeNonbondedExcl.o \
	$(DSTDIR)/ComputeNonbondedSelf.o \
	$(DSTDIR)/ComputeNonbondedPair.o \
	$(DSTDIR)/ComputeNonbondedUtil.o \
	$(DSTDIR)/ComputePatch.o \
	$(DSTDIR)/ComputePatchPair.o \
	$(DSTDIR)/ComputePme.o \
	$(DSTDIR)/ComputePmeMsgs.o \
	$(DSTDIR)/ComputeRestraints.o \
	$(DSTDIR)/ComputeSMD.o \
	$(DSTDIR)/ComputeSphericalBC.o \
	$(DSTDIR)/ComputeTcl.o \
	$(DSTDIR)/ConfigList.o \
	$(DSTDIR)/Controller.o \
	$(DSTDIR)/ccsinterface.o \
        $(DSTDIR)/FreeEnergyAssert.o \
        $(DSTDIR)/FreeEnergyGroup.o \
        $(DSTDIR)/FreeEnergyLambda.o \
        $(DSTDIR)/FreeEnergyLambdMgr.o \
        $(DSTDIR)/FreeEnergyParse.o \
        $(DSTDIR)/FreeEnergyRestrain.o \
        $(DSTDIR)/FreeEnergyRMgr.o \
        $(DSTDIR)/FreeEnergyVector.o \
	$(DSTDIR)/heap.o \
	$(DSTDIR)/HomePatch.o \
	$(DSTDIR)/IMDOutput.o \
	$(DSTDIR)/InfoStream.o \
	$(DSTDIR)/LdbCoordinator.o \
	$(DSTDIR)/LJTable.o \
	$(DSTDIR)/Measure.o \
	$(DSTDIR)/MStream.o \
	$(DSTDIR)/MigrateAtomsMsg.o \
	$(DSTDIR)/Molecule.o \
	$(DSTDIR)/NamdCentLB.o \
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
	$(DSTDIR)/PmeCoulomb.o \
	$(DSTDIR)/PmeFFT.o \
	$(DSTDIR)/PmeKSpace.o \
	$(DSTDIR)/PmeRealSpace.o \
	$(DSTDIR)/ProcessorPrivate.o \
	$(DSTDIR)/ProxyMgr.o \
	$(DSTDIR)/ProxyPatch.o \
	$(DSTDIR)/Rebalancer.o \
	$(DSTDIR)/RecBisection.o \
	$(DSTDIR)/ReductionMgr.o \
	$(DSTDIR)/RefineOnly.o \
	$(DSTDIR)/ScriptTcl.o \
	$(DSTDIR)/Sequencer.o \
	$(DSTDIR)/Set.o \
	$(DSTDIR)/SimParameters.o \
	$(DSTDIR)/Sync.o \
	$(DSTDIR)/TclCommands.o \
	$(DSTDIR)/WorkDistrib.o \
	$(DSTDIR)/pub3dfft.o \
	$(DSTDIR)/vmdsock.o \
	$(DSTDIR)/imd.o

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
		$(INCDIR)/LdbCoordinator.decl.h \
		$(INCDIR)/LdbCoordinator.def.h \
		$(INCDIR)/NamdCentLB.decl.h \
		$(INCDIR)/NamdCentLB.def.h \
		$(INCDIR)/Node.decl.h \
		$(INCDIR)/Node.def.h \
		$(INCDIR)/PatchMgr.decl.h \
		$(INCDIR)/PatchMgr.def.h \
		$(INCDIR)/ProxyMgr.decl.h \
		$(INCDIR)/ProxyMgr.def.h \
		$(INCDIR)/ReductionMgr.decl.h \
		$(INCDIR)/ReductionMgr.def.h \
		$(INCDIR)/WorkDistrib.decl.h \
		$(INCDIR)/WorkDistrib.def.h \
		$(INCDIR)/main.decl.h \
		$(INCDIR)/main.def.h

# definitions for Charm routines
CHARMC = $(CHARM)/bin/charmc
CHARMXI = $(CHARM)/bin/charmc
INCLUDE = $(CHARM)/include

# Libraries we may have changed
LIBS = $(DPMTALIBS) $(DPMELIBS) fftw/libfftw_ampi.a

# CXX is platform dependent
CXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) -Ifftw $(DPMTA) $(DPME) $(TCL) $(FFT) $(CCS) $(CXXOPTS) $(RELEASE)
CXXTHREADFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) -Ifftw $(DPMTA) $(DPME) $(TCL) $(FFT) $(CCS) $(CXXTHREADOPTS) $(RELEASE)
CXXSIMPARAMFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) -Ifftw $(DPMTA) $(DPME) $(TCL) $(FFT) $(CCS) $(CXXSIMPARAMOPTS) $(RELEASE)
GXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) -Ifftw $(DPMTA) $(DPME) $(TCL) $(FFT) $(CCS) $(RELEASE)

# Add new executables here.

BINARIES = namd2 flipdcd flipbinpdb

all:	$(BINARIES)

namd2:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(INCDIR) -I$(SRCDIR) $(CXXOPTS)" \
	-language ampi \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(DPMELIB) \
	-Lfftw -lfftw_ampi \
	$(TCLLIB) \
	$(FFTLIB)

flipdcd:	$(SRCDIR)/flipdcd.c
	$(CC) -o $@ $(SRCDIR)/flipdcd.c || \
	echo "#!/bin/sh\necho unavailable on this platform" > $@; \
	chmod +x $@

flipbinpdb:	$(SRCDIR)/flipbinpdb.c
	$(CC) -o $@ $(SRCDIR)/flipbinpdb.c || \
	echo "#!/bin/sh\necho unavailable on this platform" > $@; \
	chmod +x $@

fixdcd:	$(SRCDIR)/fixdcd.c
	$(CC) -o fixdcd $(SRCDIR)/fixdcd.c

dumpdcd:	$(SRCDIR)/dumpdcd.c
	$(CC) -o dumpdcd $(SRCDIR)/dumpdcd.c

loaddcd:	$(SRCDIR)/loaddcd.c
	$(CC) -o loaddcd $(SRCDIR)/loaddcd.c

projections:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(INCDIR) -I$(SRCDIR) $(CXXOPTS)" \
	-language ampi -tracemode projections \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(DPMELIB) \
	-Lfftw -lfftw_ampi \
	$(TCLLIB) \
	$(FFTLIB)

summary:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(INCDIR) -I$(SRCDIR) $(CXXOPTS)" \
	-language ampi -tracemode summary \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(DPMELIB) \
	-Lfftw -lfftw_ampi \
	$(TCLLIB) \
	$(FFTLIB)

$(DPMTADIR)/mpole/libmpole.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMTADIR)/src/libdpmta2.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMEDIR)/libdpme.a:
	cd $(DPMEDIR) ; $(MAKE) ; cd ..

fftw/libfftw_ampi.a:
	cd fftw ; $(MAKE) ; cd ..

# Unix commands

ECHO = echo
MOVE = mv
COPY = cp
RM = rm -f

# Explicit rules for modules.

$(INCDIR)/BroadcastMgr.def.h: $(INCDIR)/BroadcastMgr.decl.h

$(INCDIR)/BroadcastMgr.decl.h: $(SRCDIR)/BroadcastMgr.ci
	$(CHARMXI) $(SRCDIR)/BroadcastMgr.ci
	$(MOVE) BroadcastMgr.decl.h BroadcastMgr.def.h $(INCDIR)

$(INCDIR)/CollectionMaster.def.h: $(INCDIR)/CollectionMaster.decl.h

$(INCDIR)/CollectionMaster.decl.h: $(SRCDIR)/CollectionMaster.ci
	$(CHARMXI) $(SRCDIR)/CollectionMaster.ci
	$(MOVE) CollectionMaster.decl.h CollectionMaster.def.h $(INCDIR)

$(INCDIR)/CollectionMgr.def.h: $(INCDIR)/CollectionMgr.decl.h

$(INCDIR)/CollectionMgr.decl.h: $(SRCDIR)/CollectionMgr.ci
	$(CHARMXI) $(SRCDIR)/CollectionMgr.ci
	$(MOVE) CollectionMgr.decl.h CollectionMgr.def.h $(INCDIR)

$(INCDIR)/ComputeMgr.def.h: $(INCDIR)/ComputeMgr.decl.h

$(INCDIR)/ComputeMgr.decl.h: $(SRCDIR)/ComputeMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputeMgr.ci
	$(MOVE) ComputeMgr.decl.h ComputeMgr.def.h $(INCDIR)

$(INCDIR)/ComputePmeMgr.def.h: $(INCDIR)/ComputePmeMgr.decl.h

$(INCDIR)/ComputePmeMgr.decl.h: $(SRCDIR)/ComputePmeMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputePmeMgr.ci
	$(MOVE) ComputePmeMgr.decl.h ComputePmeMgr.def.h $(INCDIR)

$(INCDIR)/LdbCoordinator.def.h: $(INCDIR)/LdbCoordinator.decl.h

$(INCDIR)/LdbCoordinator.decl.h: $(SRCDIR)/LdbCoordinator.ci
	$(CHARMXI) $(SRCDIR)/LdbCoordinator.ci
	$(MOVE) LdbCoordinator.decl.h LdbCoordinator.def.h $(INCDIR)

$(INCDIR)/NamdCentLB.def.h: $(INCDIR)/NamdCentLB.decl.h

$(INCDIR)/NamdCentLB.decl.h: $(SRCDIR)/NamdCentLB.ci
	$(CHARMXI) $(SRCDIR)/NamdCentLB.ci
	$(MOVE) NamdCentLB.decl.h NamdCentLB.def.h $(INCDIR)

$(INCDIR)/Node.def.h: $(INCDIR)/Node.decl.h

$(INCDIR)/Node.decl.h: $(SRCDIR)/Node.ci
	$(CHARMXI) $(SRCDIR)/Node.ci
	$(MOVE) Node.decl.h Node.def.h $(INCDIR)

$(INCDIR)/PatchMgr.def.h: $(INCDIR)/PatchMgr.decl.h

$(INCDIR)/PatchMgr.decl.h: $(SRCDIR)/PatchMgr.ci
	$(CHARMXI) $(SRCDIR)/PatchMgr.ci
	$(MOVE) PatchMgr.decl.h PatchMgr.def.h $(INCDIR)

$(INCDIR)/ProxyMgr.def.h: $(INCDIR)/ProxyMgr.decl.h

$(INCDIR)/ProxyMgr.decl.h: $(SRCDIR)/ProxyMgr.ci
	$(CHARMXI) $(SRCDIR)/ProxyMgr.ci
	$(MOVE) ProxyMgr.decl.h ProxyMgr.def.h $(INCDIR)

$(INCDIR)/ReductionMgr.def.h: $(INCDIR)/ReductionMgr.decl.h

$(INCDIR)/ReductionMgr.decl.h: $(SRCDIR)/ReductionMgr.ci
	$(CHARMXI) $(SRCDIR)/ReductionMgr.ci
	$(MOVE) ReductionMgr.decl.h ReductionMgr.def.h $(INCDIR)

$(INCDIR)/Sync.def.h: $(INCDIR)/Sync.decl.h

$(INCDIR)/Sync.decl.h: $(SRCDIR)/Sync.ci
	$(CHARMXI) $(SRCDIR)/Sync.ci
	$(MOVE) Sync.decl.h Sync.def.h $(INCDIR)

$(INCDIR)/WorkDistrib.def.h: $(INCDIR)/WorkDistrib.decl.h

$(INCDIR)/WorkDistrib.decl.h: $(SRCDIR)/WorkDistrib.ci
	$(CHARMXI) $(SRCDIR)/WorkDistrib.ci
	$(MOVE) WorkDistrib.decl.h WorkDistrib.def.h $(INCDIR)

$(INCDIR)/main.def.h: $(INCDIR)/main.decl.h

$(INCDIR)/main.decl.h: $(SRCDIR)/main.ci
	$(CHARMXI) $(SRCDIR)/main.ci
	$(MOVE) main.decl.h main.def.h $(INCDIR)

DEPENDFILE = .rootdir/Make.depends

# make depends is ugly!  The problem: we have obj/file.o and want src/file.C.
# Solution: heavy use of basename and awk.
# This is a CPU killer...  Don't make depends if you don't need to.
depends: $(INCDIR) $(CIFILES) $(DSTDIR) $(DEPENDFILE)
	$(ECHO) "Creating " $(DEPENDFILE) " ..."; \
	if [ -f $(DEPENDFILE) ]; then \
	   $(MOVE) -f $(DEPENDFILE) $(DEPENDFILE).old; \
	fi; \
	touch $(DEPENDFILE); \
	for i in $(OBJS) ; do \
	      $(ECHO) "checking dependencies for" \
	        `basename $$i | awk -F. '{print $$1".C"}'` ; \
	      g++ -MM $(GXXFLAGS) \
	        $(SRCDIR)/`basename $$i | awk -F. '{print $$1".C"}'` | \
	      perl $(SRCDIR)/dc.pl $(INCLUDE) /usr/include /usr/local >> $(DEPENDFILE); \
	      $(ECHO) '	$$(CXX) $$(CXXFLAGS)' -o $$i -c \
	        $(SRCDIR)/`basename $$i | awk -F. '{print $$1".C"}'` \
		>> $(DEPENDFILE) ; \
	done; \
	$(RM) $(DEPENDFILE).sed; \
	sed -e "/obj\/Controller.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/Sequencer.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/ComputeFullDirect.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/ReductionMgr.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/SimParameters.o/ s/CXXFLAGS/CXXSIMPARAMFLAGS/" \
	    $(DEPENDFILE) > $(DEPENDFILE).sed; \
	$(MOVE) -f $(DEPENDFILE).sed $(DEPENDFILE);

$(DEPENDFILE):
	touch $(DEPENDFILE)

include	$(DEPENDFILE)

$(DSTDIR):
	mkdir $(DSTDIR)

$(INCDIR):
	mkdir $(INCDIR)

clean:
	rm -rf ptrepository Templates.DB SunWS_cache $(DSTDIR) $(INCDIR)
	cd $(DPMTADIR) ; $(MAKE) clean ; cd ..
	cd $(DPMEDIR) ; $(MAKE) clean ; cd ..
	cd fftw ; $(MAKE) clean ; cd ..

veryclean:	clean
	rm -f $(BINARIES)

RELEASE_DIR_NAME = NAMD_$(NAMD_VERSION)_$(NAMD_PLATFORM)
RELEASE_FILES = .rootdir/README.txt \
		.rootdir/announce.txt \
		.rootdir/license.txt \
		.rootdir/notes.txt \
		namd2 flipdcd flipbinpdb

release: all
	$(ECHO) Creating release $(RELEASE_DIR_NAME)
	mkdir $(RELEASE_DIR_NAME)
	cp $(RELEASE_FILES) $(RELEASE_DIR_NAME)
	if [ -r conv-host ]; then \
	   $(COPY) conv-host $(RELEASE_DIR_NAME); \
	   $(ECHO) "group main" > $(RELEASE_DIR_NAME)/nodelist; \
	   $(ECHO) " host localhost" >> $(RELEASE_DIR_NAME)/nodelist; \
	fi
	chmod -R a+rX $(RELEASE_DIR_NAME)
	tar cf $(RELEASE_DIR_NAME).tar $(RELEASE_DIR_NAME)
	gzip $(RELEASE_DIR_NAME).tar

