CHARMC = /Projects/l1/namd.2.0/charm/bin/charmc
CHARMXI = /Projects/l1/namd.2.0/charm/bin/charmc

# source directory
SRCDIR = src
# destination directory (binaries) -- currently, MUST be .
DSTDIR = .
# temp include directory for cifiles
INCDIR = inc

CXXOPTS = -g
CXX = CC -Aa -D_HPUX_SOURCE
INCLUDE = /Projects/l1/namd.2.0/charm/include
CXXFLAGS = -I$(INCLUDE) -I$(SRCDIR) -I$(INCDIR) $(CXXOPTS)  -w

.SUFFIXES: 	.ci

DEPENDFILE = Make.depends
ECHO = echo
MOVE = mv

OBJS = \
	main.o \
	AtomMap.o \
	common.o strlib.o \
	CollectionMaster.o CollectionMgr.o \
	Communicate.o CommunicateConverse.o \
	Compute.o \
	ComputeAngles.o ComputeBonds.o ComputeDihedrals.o ComputeImpropers.o \
	ComputeGeneral.o ComputeMap.o ComputeMgr.o \
	ComputeNonbondedExcl.o ComputeNonbondedSelf.o \
	ComputeNonbondedPair.o ComputeNonbondedUtil.o \
	ComputePatch.o ComputePatchPair.o \
	ConfigList.o \
	Controller.o \
	HomePatch.o \
	Inform.o InfoStream.o IntTree.o \
	LJTable.o \
	Message.o MessageManager.o MessageQueue.o \
	Molecule.o \
	Namd.o NamdState.o \
	Node.o \
	ParseOptions.o Parameters.o \
	SimParameters.o PDB.o PDBData.o \
	Patch.o PatchMgr.o PatchMap.o \
	ProxyMgr.o ProxyPatch.o \
	ReductionMgr.o \
	Sequencer.o \
	WorkDistrib.o

CXXFILES = $(OBJS:.o=.C)

INTERFACES = main.ci Node.ci WorkDistrib.ci PatchMgr.ci Compute.ci \
		ComputeMgr.ci ProxyMgr.ci ReductionMgr.ci \
		CollectionMgr.ci CollectionMaster.ci

namd2:	$(INCDIR) $(DSTDIR) $(OBJS)
	$(CHARMC) -ld++-option "-I $(INCLUDE) -I $(SRCDIR)" -g -language charm++ -o namd2 $(OBJS)

cifiles:	$(INCDIR) $(DSTDIR)
	for i in $(INTERFACES); do \
	   $(CHARMXI) $(SRCDIR)/$$i; \
	done;
	$(MOVE) $(SRCDIR)/*.top.h $(INCDIR)
	$(MOVE) $(SRCDIR)/*.bot.h $(INCDIR)

depends: $(DEPENDSFILE)
	$(ECHO) "Creating " $(DEPENDFILE) " ..."; \
	if [ -f $(DEPENDFILE) ]; then \
	   $(MOVE) -f $(DEPENDFILE) $(DEPENDFILE).old; \
	fi; \
	touch $(DEPENDFILE); \
	for i in ZZZ $(CXXFILES) ; do \
	   if [ "$$i" != "ZZZ" -a -f $(SRCDIR)/$$i ]; then \
	      $(ECHO) "checking dependencies for $$i ..."; \
	      g++ -MM $(CXXFLAGS) $(SRCDIR)/$$i |  \
	      $(SRCDIR)/dc.pl $(INCLUDE) /usr/include >> $(DEPENDFILE);\
	   fi; \
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
	rm -f *.o
	rm -f $(SRCDIR)/*.o
	rm -rf ptrepository

veryclean:
	rm -f *.o
	rm -rf ptrepository
	rm -f $(INCDIR)/*.top.h $(INCDIR)/*.bot.h *.depends

