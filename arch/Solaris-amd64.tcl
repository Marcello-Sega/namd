
TCLDIR=/Projects/namd2/tcl/solx86_64
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -ltcl8.4 -lnsl -lsocket
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

