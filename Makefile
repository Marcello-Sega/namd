
#####
# Directories
#####
# source directory
SRCDIR = src
# destination directory (binaries) -- currently, MUST be .
DSTDIR = obj
# temp include directory for cifiles
INCDIR = inc
# Libraries we may have changed
LIBS = dpmta2/mpole/libmpole.a dpmta2/src/libdpmta2.a pvm3/libpvmc.a

#####
# Flags
#####
# uncomment for no warnings during compile
#NOWARN=-w
# uncomment for purify during compile
#PURIFY=-purify

#####
# definitions for Charm routines
#####
CHARM = /Projects/l1/namd.2.0/charm
CHARMC = $(CHARM)/bin/charmc $(PURIFY)
CHARMXI = $(CHARM)/bin/charmc $(PURIFY)

#####
# definitions for (D)PMTA routines
#####
DPMTADIR=dpmta2
DPMTAINCL=-I$(DPMTADIR)/include
DPMTALIB=-L$(DPMTADIR) -ldpmta2 -lmpole
DPMTAFLAGS=-DDPMTA
DPMTA=$(DPMTAINCL) $(DPMTAFLAGS)
######
## definitions for PVM routines
######
PVMDIR=pvm3
PVMLIB=-L$(PVMDIR) -lpvmc
PVM=-I$(PVMDIR)

#####
## Choose your favorite Compiler options
####
#CXXOPTS = -O
#CXXOPTS = +O3 -G -z 
#CXXOPTS = -O +DAK460 +DSK460
CXXOPTS = -g

CXX = CC -Aa -D_HPUX_SOURCE
INCLUDE = /Projects/l1/namd.2.0/charm/include
CXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(PVM) $(CXXOPTS) $(NOWARN)
GXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(DPMTA) $(PVM) $(NOWARN)

.SUFFIXES: 	.ci

DEPENDFILE = Make.depends
ECHO = echo
MOVE = mv
COPY = cp

OBJS = \
	$(DSTDIR)/common.o \
	$(DSTDIR)/dcdlib.o \
	$(DSTDIR)/main.o \
	$(DSTDIR)/strlib.o \
	$(DSTDIR)/AtomMap.o \
	$(DSTDIR)/CollectionMaster.o \
	$(DSTDIR)/CollectionMgr.o \
	$(DSTDIR)/Communicate.o \
	$(DSTDIR)/CommunicateConverse.o \
	$(DSTDIR)/Compute.o \
	$(DSTDIR)/ComputeAngles.o \
	$(DSTDIR)/ComputeBonds.o \
	$(DSTDIR)/ComputeDihedrals.o \
	$(DSTDIR)/ComputeDPMTA.o \
	$(DSTDIR)/ComputeFullDirect.o \
	$(DSTDIR)/ComputeHomePatches.o \
	$(DSTDIR)/ComputeImpropers.o \
	$(DSTDIR)/ComputeGeneral.o \
	$(DSTDIR)/ComputeMap.o \
	$(DSTDIR)/ComputeMgr.o \
	$(DSTDIR)/ComputeNonbondedExcl.o \
	$(DSTDIR)/ComputeNonbondedSelf.o \
	$(DSTDIR)/ComputeNonbondedPair.o \
	$(DSTDIR)/ComputeNonbondedUtil.o \
	$(DSTDIR)/ComputePatch.o \
	$(DSTDIR)/ComputePatchPair.o \
	$(DSTDIR)/ConfigList.o \
	$(DSTDIR)/Controller.o \
	$(DSTDIR)/HomePatch.o \
	$(DSTDIR)/HBondParam.o \
	$(DSTDIR)/Inform.o \
	$(DSTDIR)/InfoStream.o \
	$(DSTDIR)/IntTree.o \
	$(DSTDIR)/LJTable.o \
	$(DSTDIR)/Message.o \
	$(DSTDIR)/MessageManager.o \
	$(DSTDIR)/MessageQueue.o \
	$(DSTDIR)/MigrateAtomsMsg.o \
	$(DSTDIR)/Molecule.o \
	$(DSTDIR)/Namd.o \
	$(DSTDIR)/NamdState.o \
	$(DSTDIR)/Node.o \
	$(DSTDIR)/Output.o \
	$(DSTDIR)/Parameters.o \
	$(DSTDIR)/ParseOptions.o \
	$(DSTDIR)/Patch.o \
	$(DSTDIR)/PatchMgr.o \
	$(DSTDIR)/PatchMap.o \
	$(DSTDIR)/PDB.o \
	$(DSTDIR)/PDBData.o \
	$(DSTDIR)/ProxyMgr.o \
	$(DSTDIR)/ProxyPatch.o \
	$(DSTDIR)/ReductionMgr.o \
	$(DSTDIR)/Sequencer.o \
	$(DSTDIR)/SimParameters.o \
	$(DSTDIR)/VoidTree.o \
	$(DSTDIR)/WorkDistrib.o

INTERFACES = main.ci Node.ci WorkDistrib.ci PatchMgr.ci Compute.ci \
		ComputeMgr.ci ProxyMgr.ci ReductionMgr.ci \
		CollectionMgr.ci CollectionMaster.ci

TEMPLATES = \
	$(SRCDIR)/ComputeHomeTuples.C \
	$(SRCDIR)/PositionBox.C \
	$(SRCDIR)/PositionOwnerBox.C \
	$(SRCDIR)/Templates/Box.C \
	$(SRCDIR)/Templates/OwnerBox.C \
	$(SRCDIR)/Templates/ResizeArray.C \
	$(SRCDIR)/Templates/SortableResizeArray.C \
	$(SRCDIR)/Templates/SortedArray.C \
	$(SRCDIR)/Templates/UniqueSortedArray.C

namd2:	$(INCDIR) $(DSTDIR) $(OBJS) $(TEMPLATES) $(LIBS)
	cd $(PVMDIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..
	cd $(DPMTADIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..
	$(CHARMC) -verbose -ld++-option \
	"-I $(INCLUDE) -I $(SRCDIR) $(CXXOPTS) " \
	-language charm++ \
	-o namd2 $(OBJS) \
	$(DPMTALIB) \
	$(PVMLIB) \
	# Now sit back, have a coke, and relax.

dpmta2/mpole/libmpole.a:
	cd $(DPMTADIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..

dpmta2/src/libdpmta2.a:
	cd $(DPMTADIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..

pvm3/libpvmc.a:
	cd $(PVMDIR) ; $(MAKE) CHARM=$(CHARM) ; cd ..

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
	      $(SRCDIR)/dc.pl $(INCLUDE) /usr/include /usr/local >> $(DEPENDFILE); \
	      $(ECHO) "\t$$(CXX) $$(CXXFLAGS) -o "$$i" -c" \
	        "$$(SRCDIR)/"`basename $$i | awk -F. '{print $$1".C"}'` \
		>> $(DEPENDFILE) ; \
	done;

Make.depends:
	touch $(DEPENDSFILE)

include	$(DEPENDFILE)

$(INTERFACES:.ci=.top.h):	$(INCDIR) $$(@:.top.h=.ci)
	$(CHARMXI) $?
	$(MOVE) $(SRCDIR)/*.top.h $(INCDIR)

$(INTERFACES:.ci=.bot.h):	$(INCDIR) $$(@:.bot.h=.ci)
	$(CHARMXI) $?
	$(MOVE) $(SRCDIR)/*.bot.h $(INCDIR)

$(DSTDIR):
	mkdir $(DSTDIR)

$(INCDIR):
	mkdir $(INCDIR)

clean:
	rm -rf ptrepository
	rm -rf $(DSTDIR)
	rm -f namd2
#	cd $(DPMTADIR) ; $(MAKE) clean ; cd ..

veryclean:	clean
	rm -rf $(INCDIR)
	rm -f *.depends
	# allow for the makefile to continue to work
	touch $(DEPENDFILE)

accesslist:
	cvs admin -aari,brunner,jim,milind,nealk .
