
TCLDIR=/Projects/namd2/tcl/solaris
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -ltcl8.1
TCLAPPLIB=-lsocket -lnsl
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

