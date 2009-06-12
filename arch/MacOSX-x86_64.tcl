
TCLDIR=/Projects/namd2/tcl/macosx-x86
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -lnamdtcl8.4
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

