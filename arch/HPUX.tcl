
TCLDIR=/Projects/namd2/tcl/hpux
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -ltcl
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

