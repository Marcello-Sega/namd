include Makearch

#####
# Directories
#####

# source directory
SRCDIR = src
# destination directory (binaries) -- currently, MUST be .
DSTDIR = obj
# temp include directory for .top.h and .bot.h files
INCDIR = inc

#####
# definitions for Charm routines
#####

# CHARM is platform dependent
CHARMC = $(CHARM)/bin/charmc $(PURIFY)
CHARMXI = $(CHARM)/bin/charmc $(PURIFY)


#####
# definitions for (D)PMTA routines
#####

DPMTADIR=dpmta-2.6
DPMTAINCL=-I$(DPMTADIR)/mpole -I$(DPMTADIR)/src
DPMTALIB=-L$(DPMTADIR)/mpole -L$(DPMTADIR)/src -ldpmta2 -lmpole
DPMTAFLAGS=-DDPMTA
DPMTA=$(DPMTAINCL) $(DPMTAFLAGS)
DPMTALIBS=$(DPMTADIR)/mpole/libmpole.a $(DPMTADIR)/src/libdpmta2.a


#####
# definitions for DPME routines
#####

DPMEDIR=dpme2
DPMEINCL=-I$(DPMEDIR)
DPMELIB=-L$(DPMEDIR) -ldpme
DPMEFLAGS=-DDPME
DPME=$(DPMEINCL) $(DPMEFLAGS)
DPMELIBS= $(DPMEDIR)/libdpme.a


######
## definitions for PVM routines
######

PVMDIR=pvm3
PVMLIB=-L$(PVMDIR) -lpvmc
PVM=-I$(PVMDIR)
PVMLIBS=pvm3/libpvmc.a


######
## Libraries we may have changed
######

LIBS = $(DPMTALIBS) $(PVMLIBS) $(DPMELIBS)


# CXX is platform dependent
INCLUDE = $(CHARM)/include
CXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(DPME) $(PVM) $(MDCOMM) $(TCL) $(CXXOPTS) $(NOWARN) $(NAMDFLAGS)
CXXTHREADFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(DPME) $(PVM) $(MDCOMM) $(TCL) $(CXXTHREADOPTS) $(NOWARN) $(NAMDFLAGS)
GXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(DPME) $(PVM) $(MDCOMM) $(TCL) $(NOWARN) $(NAMDFLAGS)

DEPENDFILE = Make.depends

# Add new source files here.

OBJS = \
	$(DSTDIR)/common.o \
	$(DSTDIR)/dcdlib.o \
	$(DSTDIR)/main.o \
	$(DSTDIR)/mdcomm.o \
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
	$(DSTDIR)/ComputeHomePatches.o \
	$(DSTDIR)/ComputeImpropers.o \
	$(DSTDIR)/ComputeGeneral.o \
	$(DSTDIR)/ComputeGlobal.o \
	$(DSTDIR)/ComputeGlobalMaster.o \
	$(DSTDIR)/ComputeGlobalMsgs.o \
	$(DSTDIR)/ComputeMap.o \
	$(DSTDIR)/ComputeMgr.o \
	$(DSTDIR)/ComputeMDComm.o \
	$(DSTDIR)/ComputeNonbondedExcl.o \
	$(DSTDIR)/ComputeNonbondedSelf.o \
	$(DSTDIR)/ComputeNonbondedPair.o \
	$(DSTDIR)/ComputeNonbondedUtil.o \
	$(DSTDIR)/ComputePatch.o \
	$(DSTDIR)/ComputePatchPair.o \
	$(DSTDIR)/ComputeRestraints.o \
	$(DSTDIR)/ComputeSMD.o \
	$(DSTDIR)/ComputeSphericalBC.o \
	$(DSTDIR)/ComputeTcl.o \
	$(DSTDIR)/ConfigList.o \
	$(DSTDIR)/Controller.o \
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
	$(DSTDIR)/Inform.o \
	$(DSTDIR)/InfoStream.o \
	$(DSTDIR)/IntTree.o \
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
	$(DSTDIR)/ProcessorPrivate.o \
	$(DSTDIR)/ProxyMgr.o \
	$(DSTDIR)/ProxyPatch.o \
	$(DSTDIR)/Rebalancer.o \
	$(DSTDIR)/RecBisection.o \
	$(DSTDIR)/ReductionMgr.o \
	$(DSTDIR)/RefineOnly.o \
	$(DSTDIR)/Sequencer.o \
	$(DSTDIR)/Set.o \
	$(DSTDIR)/SimParameters.o \
	$(DSTDIR)/SMD.o \
	$(DSTDIR)/SMDMsgs.o \
	$(DSTDIR)/TclCommands.o \
	$(DSTDIR)/TestController.o \
	$(DSTDIR)/TestSequencer.o \
	$(DSTDIR)/VoidTree.o \
	$(DSTDIR)/WorkDistrib.o

# Add new modules here and also define explicit rule below.

CIFILES = 	\
		$(INCDIR)/BroadcastMgr.top.h \
		$(INCDIR)/BroadcastMgr.bot.h \
		$(INCDIR)/CollectionMaster.top.h \
		$(INCDIR)/CollectionMaster.bot.h \
		$(INCDIR)/CollectionMgr.top.h \
		$(INCDIR)/CollectionMgr.bot.h \
		$(INCDIR)/ComputeMgr.top.h \
		$(INCDIR)/ComputeMgr.bot.h \
		$(INCDIR)/LdbCoordinator.top.h \
		$(INCDIR)/LdbCoordinator.bot.h \
		$(INCDIR)/Node.top.h \
		$(INCDIR)/Node.bot.h \
		$(INCDIR)/PatchMgr.top.h \
		$(INCDIR)/PatchMgr.bot.h \
		$(INCDIR)/ProxyMgr.top.h \
		$(INCDIR)/ProxyMgr.bot.h \
		$(INCDIR)/ReductionMgr.top.h \
		$(INCDIR)/ReductionMgr.bot.h \
		$(INCDIR)/WorkDistrib.top.h \
		$(INCDIR)/WorkDistrib.bot.h \
		$(INCDIR)/main.top.h \
		$(INCDIR)/main.bot.h

# Add new executables here.

BINARIES = namd2 flipdcd flipbinpdb

all:	$(BINARIES)

namd2:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(SRCDIR) $(CXXOPTS) " \
	-language charm++ \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(PVMLIB) \
	$(MDCOMMLIB) \
	$(TCLLIB)

flipdcd:	$(SRCDIR)/flipdcd.c
	$(CC) -o flipdcd $(SRCDIR)/flipdcd.c

flipbinpdb:	$(SRCDIR)/flipbinpdb.c
	$(CC) -o flipbinpdb $(SRCDIR)/flipbinpdb.c

# Now sit back, have a coke, and relax.

projections:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	cd $(PVMDIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..
	cd $(DPMTADIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..
	cd $(DPMEDIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(SRCDIR) $(CXXOPTS) " \
	-language charm++ -tracemode projections \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(PVMLIB) \
	$(MDCOMMLIB) \
	$(TCLLIB)

# Now sit back, have a coke, and relax.

$(DPMTADIR)/mpole/libmpole.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMTADIR)/src/libdpmta2.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMEDIR)/libdpme.a:
	cd $(DPMEDIR) ; $(MAKE) ; cd ..

pvm3/libpvmc.a:
	cd $(PVMDIR) ; $(MAKE) ; cd ..

# Explicit rules for modules.

MOVECIFILES = $(MOVE) $(SRCDIR)/*.top.h $(SRCDIR)/*.bot.h $(INCDIR)

$(INCDIR)/BroadcastMgr.top.h:	$(SRCDIR)/BroadcastMgr.ci
	$(CHARMXI) $(SRCDIR)/BroadcastMgr.ci
	$(MOVECIFILES)

$(INCDIR)/BroadcastMgr.bot.h:	$(SRCDIR)/BroadcastMgr.ci
	$(CHARMXI) $(SRCDIR)/BroadcastMgr.ci
	$(MOVECIFILES)

$(INCDIR)/CollectionMaster.top.h:	$(SRCDIR)/CollectionMaster.ci
	$(CHARMXI) $(SRCDIR)/CollectionMaster.ci
	$(MOVECIFILES)

$(INCDIR)/CollectionMaster.bot.h:	$(SRCDIR)/CollectionMaster.ci
	$(CHARMXI) $(SRCDIR)/CollectionMaster.ci
	$(MOVECIFILES)

$(INCDIR)/CollectionMgr.top.h:	$(SRCDIR)/CollectionMgr.ci
	$(CHARMXI) $(SRCDIR)/CollectionMgr.ci
	$(MOVECIFILES)

$(INCDIR)/CollectionMgr.bot.h:	$(SRCDIR)/CollectionMgr.ci
	$(CHARMXI) $(SRCDIR)/CollectionMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ComputeMgr.top.h:	$(SRCDIR)/ComputeMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputeMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ComputeMgr.bot.h:	$(SRCDIR)/ComputeMgr.ci
	$(CHARMXI) $(SRCDIR)/ComputeMgr.ci
	$(MOVECIFILES)

$(INCDIR)/LdbCoordinator.top.h:	$(SRCDIR)/LdbCoordinator.ci
	$(CHARMXI) $(SRCDIR)/LdbCoordinator.ci
	$(MOVECIFILES)

$(INCDIR)/LdbCoordinator.bot.h:	$(SRCDIR)/LdbCoordinator.ci
	$(CHARMXI) $(SRCDIR)/LdbCoordinator.ci
	$(MOVECIFILES)

$(INCDIR)/Node.top.h:	$(SRCDIR)/Node.ci
	$(CHARMXI) $(SRCDIR)/Node.ci
	$(MOVECIFILES)

$(INCDIR)/Node.bot.h:	$(SRCDIR)/Node.ci
	$(CHARMXI) $(SRCDIR)/Node.ci
	$(MOVECIFILES)

$(INCDIR)/PatchMgr.top.h:	$(SRCDIR)/PatchMgr.ci
	$(CHARMXI) $(SRCDIR)/PatchMgr.ci
	$(MOVECIFILES)

$(INCDIR)/PatchMgr.bot.h:	$(SRCDIR)/PatchMgr.ci
	$(CHARMXI) $(SRCDIR)/PatchMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ProxyMgr.top.h:	$(SRCDIR)/ProxyMgr.ci
	$(CHARMXI) $(SRCDIR)/ProxyMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ProxyMgr.bot.h:	$(SRCDIR)/ProxyMgr.ci
	$(CHARMXI) $(SRCDIR)/ProxyMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ReductionMgr.top.h:	$(SRCDIR)/ReductionMgr.ci
	$(CHARMXI) $(SRCDIR)/ReductionMgr.ci
	$(MOVECIFILES)

$(INCDIR)/ReductionMgr.bot.h:	$(SRCDIR)/ReductionMgr.ci
	$(CHARMXI) $(SRCDIR)/ReductionMgr.ci
	$(MOVECIFILES)

$(INCDIR)/WorkDistrib.top.h:	$(SRCDIR)/WorkDistrib.ci
	$(CHARMXI) $(SRCDIR)/WorkDistrib.ci
	$(MOVECIFILES)

$(INCDIR)/WorkDistrib.bot.h:	$(SRCDIR)/WorkDistrib.ci
	$(CHARMXI) $(SRCDIR)/WorkDistrib.ci
	$(MOVECIFILES)

$(INCDIR)/main.top.h:	$(SRCDIR)/main.ci
	$(CHARMXI) $(SRCDIR)/main.ci
	$(MOVECIFILES)

$(INCDIR)/main.bot.h:	$(SRCDIR)/main.ci
	$(CHARMXI) $(SRCDIR)/main.ci
	$(MOVECIFILES)

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
	rm -rf ptrepository Templates.DB $(DSTDIR) $(INCDIR)
	cd $(DPMTADIR) ; $(MAKE) clean ; cd ..
	cd $(PVMDIR) ; $(MAKE) clean ; cd ..
	cd $(DPMEDIR) ; $(MAKE) clean ; cd ..

veryclean:	clean
	rm -f $(BINARIES)

accesslist:
	cvs admin -aari,brunner,jim,milind,nealk .
