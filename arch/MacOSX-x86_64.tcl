
TCLDIR=/Projects/namd2/tcl/macosx-x86_64
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -lnamdtcl8.4 -framework CoreFoundation
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

