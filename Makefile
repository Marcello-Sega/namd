CHARMC = /Projects/l1/namd.2.0/charm/bin/charmc
CHARMXI = /Projects/l1/namd.2.0/charm/bin/charmc

CXX = CC -Aa -D_HPUX_SOURCE -g
CXXFLAGS = -I/Projects/l1/namd.2.0/charm/include

.SUFFIXES: 	.ci

DEPENDFILE = Make.depends
ECHO = echo
MOVE = mv

OBJS = \
	main.o Message.o Molecule.o PDB.o PDBData.o ConfigList.o Inform.o Parameters.o common.o strlib.o SimParameters.o ParseOptions.o

CXXFILES = $(OBJS:.o=.C)

INTERFACES = main.ci

namd2:	$(OBJS)
	$(CHARMC) -g -language charm++ -o namd2 $(OBJS)

cifiles:
	for i in $(INTERFACES); do \
	   $(CHARMXI) $$i; \
	done;

depends: 
	$(ECHO) "Creating " $(DEPENDFILE) " ..."; \
	if [ -f $(DEPENDFILE) ]; then \
	   $(MOVE) -f $(DEPENDFILE) $(DEPENDFILE).old; \
	fi; \
	touch $(DEPENDFILE); \
	for i in ZZZ $(CXXFILES) ; do \
	   if [ "$$i" != "ZZZ" -a -f $$i ]; then \
	      $(ECHO) "checking dependencies for $$i ..."; \
	      g++ -MM $(CXXFLAGS) $$i |grep -v "/usr/include" >> $(DEPENDFILE);\
	   fi; \
	done;

include	$(DEPENDFILE)

$(INTERFACES:.ci=.top.h):	$$(@:.top.h=.ci)
	$(CHARMXI) $?

$(INTERFACES:.ci=.bot.h):	$$(@:.bot.h=.ci)
	$(CHARMXI) $?


clean:
	rm -rf *.o
