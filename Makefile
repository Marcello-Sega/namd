CHARMC = /Projects/l1/namd.2.0/charm/bin/charmc
CHARMXI = /Projects/l1/namd.2.0/charm/bin/charmc

CXXOPTS = -g
CXX = CC -Aa -D_HPUX_SOURCE
INCLUDE = /Projects/l1/namd.2.0/charm/include
CXXFLAGS = -I$(INCLUDE) $(CXXOPTS)  -w

.SUFFIXES: 	.ci

DEPENDFILE = Make.depends
ECHO = echo
MOVE = mv

OBJS = \
	main.o Message.o Molecule.o PDB.o PDBData.o \
	ConfigList.o Inform.o InfoStream.o Parameters.o common.o \
	strlib.o SimParameters.o ParseOptions.o Namd.o \
	NamdState.o WorkDistrib.o Node.o PatchMap.o ComputeMap.o \
	PatchMgr.o Patch.o HomePatch.o Sequencer.o Compute.o \
	ComputeGeneral.o AtomMap.o ComputePatchPair.o \
	ComputePatch.o ComputeNonbondedUtil.o LJTable.o \
	ComputeNonbondedSelf.o ComputeNonbondedPair.o \
	ComputeAngles.o ComputeDihedrals.o ComputeImpropers.o \
	ComputeBonds.o ComputeNonbondedExcl.o ComputeMgr.o \
	ProxyMgr.o ProxyPatch.o CommunicateConverse.o \
	Communicate.o IntTree.o MessageQueue.o MessageManager.o \
	ReductionMgr.o

CXXFILES = $(OBJS:.o=.C)

INTERFACES = main.ci Node.ci WorkDistrib.ci PatchMgr.ci Compute.ci \
		ComputeMgr.ci ProxyMgr.ci

namd2:	$(OBJS)
	$(CHARMC) -ld++-option "-I $(INCLUDE)" -g -language charm++ -o namd2 $(OBJS)

cifiles:
	for i in $(INTERFACES); do \
	   $(CHARMXI) $$i; \
	done;

depends: $(DEPENDSFILE)
	$(ECHO) "Creating " $(DEPENDFILE) " ..."; \
	if [ -f $(DEPENDFILE) ]; then \
	   $(MOVE) -f $(DEPENDFILE) $(DEPENDFILE).old; \
	fi; \
	touch $(DEPENDFILE); \
	for i in ZZZ $(CXXFILES) ; do \
	   if [ "$$i" != "ZZZ" -a -f $$i ]; then \
	      $(ECHO) "checking dependencies for $$i ..."; \
	      g++ -MM $(CXXFLAGS) $$i |  \
	      dc.pl $(INCLUDE) /usr/include >> $(DEPENDFILE);\
	   fi; \
	done;

Make.depends:
	touch $(DEPENDSFILE)

include	$(DEPENDFILE)

$(INTERFACES:.ci=.top.h):	$$(@:.top.h=.ci)
	$(CHARMXI) $?

$(INTERFACES:.ci=.bot.h):	$$(@:.bot.h=.ci)
	$(CHARMXI) $?


clean:
	rm -f *.o
	rm -rf ptrepository

veryclean:
	rm -f *.o
	rm -rf ptrepository
	rm -f *.top.h *.bot.h *.depends
