
TCLDIR=/Projects/namd2/tcl/alpha
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -ltcl
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

