include Makearch

#####
# Directories
#####

# source directory
SRCDIR = src
# destination directory (binaries) -- currently, MUST be .
DSTDIR = obj
# temp include directory for cifiles
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
DPMEINCL=-I$(DPMEDIR)/include
DPMELIB=-L$(DPMEDIR) -ldpme2 -lmpole
DPMEFLAGS=-DDPME
DPME=$(DPMEINCL) $(DPMEFLAGS)
#DPMELIBS= dpme2/libdpme2.a


######
## definitions for PVM routines
######

PVMDIR=pvm3
PVMLIB=-L$(PVMDIR) -lpvmc
PVM=-I$(PVMDIR)
PVMLIBS=pvm3/libpvmc.a

######
## definitions for TCL interface
######
##
## MOVED TO Makearch.* !!!
##
######
##
#TCLDIR=/usr/local
#TCLINCL=-I$(TCLDIR)/include
#TCLLIB=-L$(TCLDIR)/lib -ltclx -ltcl
#TCLFLAGS=-DNAMD_TCL
#TCL=$(TCLINCL) $(TCLFLAGS)

######
## Libraries we may have changed
######

LIBS = $(DPMTALIBS) $(PVMLIBS) $(DPMELIBS)


# CXX is platform dependent
INCLUDE = $(CHARM)/include
CXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(PVM) $(TCL) $(CXXOPTS) $(NOWARN) $(NAMDFLAGS)
GXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(PVM) $(TCL) $(NOWARN) $(NAMDFLAGS)

.SUFFIXES: 	.ci

DEPENDFILE = Make.depends

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
	$(DSTDIR)/ComputeDPMTA.o \
	$(DSTDIR)/ComputeFullDirect.o \
	$(DSTDIR)/ComputeHomePatches.o \
	$(DSTDIR)/ComputeImpropers.o \
	$(DSTDIR)/ComputeGeneral.o \
	$(DSTDIR)/ComputeGlobal.o \
	$(DSTDIR)/ComputeGlobalMaster.o \
	$(DSTDIR)/ComputeGlobalMsgs.o \
	$(DSTDIR)/ComputeMap.o \
	$(DSTDIR)/ComputeMgr.o \
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
	$(DSTDIR)/VoidTree.o \
	$(DSTDIR)/WorkDistrib.o

#Compute.ci - nolonger necessary(?)

INTERFACES = main.ci Node.ci WorkDistrib.ci PatchMgr.ci \
		ComputeMgr.ci ProxyMgr.ci ReductionMgr.ci \
		CollectionMgr.ci CollectionMaster.ci BroadcastMgr.ci \
		LdbCoordinator.ci

namd2:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(SRCDIR) $(CXXOPTS) " \
	-language charm++ \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(PVMLIB) \
	$(TCLLIB)

# Now sit back, have a coke, and relax.

projections:	$(INCDIR) $(DSTDIR) $(OBJS) $(LIBS)
	cd $(PVMDIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..
	cd $(DPMTADIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..
	$(CHARMC) -verbose -ld++-option \
	"-I$(INCLUDE) -I$(SRCDIR) $(CXXOPTS) " \
	-language charm++ -tracemode projections \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(PVMLIB) \
	$(TCLLIB)

# Now sit back, have a coke, and relax.

$(DPMTADIR)/mpole/libmpole.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMTADIR)/src/libdpmta2.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

dpme2/libdpme2.a:
	cd $(DPMEDIR) ; $(MAKE) ; cd ..

pvm3/libpvmc.a:
	cd $(PVMDIR) ; $(MAKE) ; cd ..

cifiles:	$(INCDIR) $(DSTDIR)
	for i in $(INTERFACES); do \
	   $(CHARMXI) $(SRCDIR)/$$i; \
	done;
	$(MOVE) $(SRCDIR)/*.top.h $(INCDIR)
	$(MOVE) $(SRCDIR)/*.bot.h $(INCDIR)

# make depends is ugly!  The problem: we have obj/file.o and want src/file.C.
# Solution: heavy use of basename and awk.
# This is a CPU killer...  Don't make depends if you don't need to.
depends: cifiles $(DSTDIR) $(DEPENDSFILE)
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
	done;

Make.depends:
	touch $(DEPENDSFILE)

include	$(DEPENDFILE)

#$(INTERFACES:.ci=.top.h):	$(INCDIR) $$(@:.top.h=.ci)
#	$(CHARMXI) $?
#	$(MOVE) $(SRCDIR)/*.top.h $(INCDIR)

#$(INTERFACES:.ci=.bot.h):	$(INCDIR) $$(@:.bot.h=.ci)
#	$(CHARMXI) $?
#	$(MOVE) $(SRCDIR)/*.bot.h $(INCDIR)

$(DSTDIR):
	mkdir $(DSTDIR)

$(INCDIR):
	mkdir $(INCDIR)

clean:
	rm -rf ptrepository
	rm -rf $(DSTDIR)
	rm -f namd2
	cd $(DPMTADIR) ; $(MAKE) clean ; cd ..
	cd $(PVMDIR) ; $(MAKE) clean ; cd ..
	cd $(DPMEDIR) ; $(MAKE) clean ; cd ..

veryclean:	clean
	rm -rf $(INCDIR)
	rm -f *.depends
	# allow for the makefile to continue to work
	touch $(DEPENDFILE)

accesslist:
	cvs admin -aari,brunner,jim,milind,nealk .
