
TCLDIR=/Projects/namd2/tcl/solaris
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -ltcl8.3 -ldl
TCLAPPLIB=-lsocket -lnsl
TCLFLAGS=-DNAMD_TCL -DUSE_NON_CONST
TCL=$(TCLINCL) $(TCLFLAGS)

