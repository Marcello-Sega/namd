include Make.charm
include Makearch

# directories
SRCDIR = src
DSTDIR = obj
INCDIR = inc
DPMTADIR=dpmta-2.6
PVMDIR=pvm3
DPMEDIR=dpme2

# comment/uncomment these lines for (D)PMTA routines
DPMTAINCL=-I$(DPMTADIR)/mpole -I$(DPMTADIR)/src -I$(PVMDIR)
DPMTALIB=-L$(DPMTADIR)/mpole -L$(DPMTADIR)/src -ldpmta2 -lmpole -L$(PVMDIR) -lpvmc
DPMTAFLAGS=-DDPMTA
DPMTA=$(DPMTAINCL) $(DPMTAFLAGS)
DPMTALIBS=$(DPMTADIR)/mpole/libmpole.a $(DPMTADIR)/src/libdpmta2.a $(PVMDIR)/libpvmc.a

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
	$(DSTDIR)/strlib.o \
	$(DSTDIR)/AlgSeven.o \
	$(DSTDIR)/AtomMap.o \
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
	$(DSTDIR)/ComputeGeneral.o \
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
        $(DSTDIR)/defmain.o \
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
	$(DSTDIR)/HBondParam.o \
	$(DSTDIR)/IMDOutput.o \
	$(DSTDIR)/Inform.o \
	$(DSTDIR)/InfoStream.o \
	$(DSTDIR)/LdbCoordinator.o \
	$(DSTDIR)/LJTable.o \
	$(DSTDIR)/MStream.o \
	$(DSTDIR)/MigrateAtomsMsg.o \
	$(DSTDIR)/Molecule.o \
	$(DSTDIR)/Namd.o \
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
	$(DSTDIR)/SMD.o \
	$(DSTDIR)/SMDMsgs.o \
	$(DSTDIR)/TclCommands.o \
	$(DSTDIR)/TestController.o \
	$(DSTDIR)/TestSequencer.o \
	$(DSTDIR)/VoidTree.o \
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
		$(INCDIR)/LdbCoordinator.decl.h \
		$(INCDIR)/LdbCoordinator.def.h \
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
LIBS = $(DPMTALIBS) $(DPMELIBS)

# CXX is platform dependent
CXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(DPME) $(TCL) $(FFT) $(CCS) $(CXXOPTS)
CXXTHREADFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(DPME) $(TCL) $(FFT) $(CCS) $(CXXTHREADOPTS)
CXXSIMPARAMFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(DPME) $(TCL) $(FFT) $(CCS) $(CXXSIMPARAMOPTS)
GXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(DPME) $(TCL) $(FFT) $(CCS)

# Add new executables here.

BINARIES = namd2 flipdcd flipbinpdb

all:	$(BINARIES)

namd2:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(INCDIR) -I$(SRCDIR) $(CXXOPTS)" \
	-language charm++ \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(TCLLIB) \
	$(FFTLIB)

flipdcd:	$(SRCDIR)/flipdcd.c
	$(CC) -o flipdcd $(SRCDIR)/flipdcd.c

flipbinpdb:	$(SRCDIR)/flipbinpdb.c
	$(CC) -o flipbinpdb $(SRCDIR)/flipbinpdb.c

fixdcd:	$(SRCDIR)/fixdcd.c
	$(CC) -o fixdcd $(SRCDIR)/fixdcd.c

dumpdcd:	$(SRCDIR)/dumpdcd.c
	$(CC) -o dumpdcd $(SRCDIR)/dumpdcd.c

loaddcd:	$(SRCDIR)/loaddcd.c
	$(CC) -o loaddcd $(SRCDIR)/loaddcd.c

projections:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(INCDIR) -I$(SRCDIR) $(CXXOPTS)" \
	-language charm++ -tracemode projections \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(TCLLIB) \
	$(FFTLIB)

$(DPMTADIR)/mpole/libmpole.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMTADIR)/src/libdpmta2.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMEDIR)/libdpme.a:
	cd $(DPMEDIR) ; $(MAKE) ; cd ..

$(PVMDIR)/libpvmc.a:
	cd $(PVMDIR) ; $(MAKE) ; cd ..

# Unix commands

ECHO = echo
MOVE = mv
COPY = cp
RM = rm -f

# Explicit rules for modules.

MOVECIFILES = $(MOVE) *.decl.h *.def.h $(INCDIR)

$(INCDIR)/BroadcastMgr.decl.h:	$(SRCDIR)/BroadcastMgr.ci
	$(CHARMXI) $(SRCDIR)/BroadcastMgr.ci
	$(MOVECIFILES)

$(INCDIR)/BroadcastMgr.def.h:	$(SRCDIR)/BroadcastMgr.ci
	$(CHARMXI) $(SRCDIR)/BroadcastMgr.ci
	$(MOVECIFILES)

$(INCDIR)/CollectionMaster.decl.h:	$(SRCDIR)/CollectionMaster.ci
	$(CHARMXI) $(SRCDIR)/CollectionMaster.ci
	$(MOVECIFILES)

$(INCDIR)/CollectionMaster.def.h:	$(SRCDIR)/CollectionMaster.ci
	$(CHARMXI) $(SRCDIR)/CollectionMaster.ci
	$(MOVECIFILES)

$(INCDIR)/CollectionMgr.decl.h:	$(SRCDIR)/CollectionMgr.ci
	$(CHARMXI) $(SRCDIR)/CollectionMgr.ci
	$(MOVECIFILES)

$(INCDIR)/CollectionMgr.def.h:	$(SRCDIR)/CollectionMgr.ci
	$(CHARMXI) $(SRCDIR)/CollectionMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ComputeMgr.decl.h:	$(SRCDIR)/ComputeMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputeMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ComputeMgr.def.h:	$(SRCDIR)/ComputeMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputeMgr.ci
	$(MOVECIFILES)

$(INCDIR)/LdbCoordinator.decl.h:	$(SRCDIR)/LdbCoordinator.ci
	$(CHARMXI) $(SRCDIR)/LdbCoordinator.ci
	$(MOVECIFILES)

$(INCDIR)/LdbCoordinator.def.h:	$(SRCDIR)/LdbCoordinator.ci
	$(CHARMXI) $(SRCDIR)/LdbCoordinator.ci
	$(MOVECIFILES)

$(INCDIR)/Node.decl.h:	$(SRCDIR)/Node.ci
	$(CHARMXI) $(SRCDIR)/Node.ci
	$(MOVECIFILES)

$(INCDIR)/Node.def.h:	$(SRCDIR)/Node.ci
	$(CHARMXI) $(SRCDIR)/Node.ci
	$(MOVECIFILES)

$(INCDIR)/PatchMgr.decl.h:	$(SRCDIR)/PatchMgr.ci
	$(CHARMXI) $(SRCDIR)/PatchMgr.ci
	$(MOVECIFILES)

$(INCDIR)/PatchMgr.def.h:	$(SRCDIR)/PatchMgr.ci
	$(CHARMXI) $(SRCDIR)/PatchMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ProxyMgr.decl.h:	$(SRCDIR)/ProxyMgr.ci
	$(CHARMXI) $(SRCDIR)/ProxyMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ProxyMgr.def.h:	$(SRCDIR)/ProxyMgr.ci
	$(CHARMXI) $(SRCDIR)/ProxyMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ReductionMgr.decl.h:	$(SRCDIR)/ReductionMgr.ci
	$(CHARMXI) $(SRCDIR)/ReductionMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ReductionMgr.def.h:	$(SRCDIR)/ReductionMgr.ci
	$(CHARMXI) $(SRCDIR)/ReductionMgr.ci
	$(MOVECIFILES)

$(INCDIR)/WorkDistrib.decl.h:	$(SRCDIR)/WorkDistrib.ci
	$(CHARMXI) $(SRCDIR)/WorkDistrib.ci
	$(MOVECIFILES)

$(INCDIR)/WorkDistrib.def.h:	$(SRCDIR)/WorkDistrib.ci
	$(CHARMXI) $(SRCDIR)/WorkDistrib.ci
	$(MOVECIFILES)

$(INCDIR)/main.decl.h:	$(SRCDIR)/main.ci
	$(CHARMXI) $(SRCDIR)/main.ci
	$(MOVECIFILES)

$(INCDIR)/main.def.h:	$(SRCDIR)/main.ci
	$(CHARMXI) $(SRCDIR)/main.ci
	$(MOVECIFILES)

DEPENDFILE = Make.depends

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
	    -e "/obj\/TestController.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/TestSequencer.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
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
	cd $(PVMDIR) ; $(MAKE) clean ; cd ..
	cd $(DPMEDIR) ; $(MAKE) clean ; cd ..

veryclean:	clean
	rm -f $(BINARIES)

NAMD_VERSION = 2.1b3
RELEASE_DIR_NAME = NAMD_$(NAMD_VERSION)_$(NAMD_ARCH)
RELEASE_FILES = .rootdir/README.txt \
		.rootdir/announce.txt \
		.rootdir/license.txt \
		.rootdir/notes.txt \
		namd2 flipdcd flipbinpdb

release: all
	$(ECHO) Creating release $(RELEASE_DIR_NAME)
	mkdir $(RELEASE_DIR_NAME)
	cp $(RELEASE_FILES) $(RELEASE_DIR_NAME)
	if [ -e conv-host ]; then \
	   $(COPY) conv-host $(RELEASE_DIR_NAME); \
	   $(ECHO) "group main" > $(RELEASE_DIR_NAME)/nodelist; \
	   $(ECHO) " host localhost" >> $(RELEASE_DIR_NAME)/nodelist; \
	fi
	tar cf $(RELEASE_DIR_NAME).tar $(RELEASE_DIR_NAME)
	gzip $(RELEASE_DIR_NAME).tar

