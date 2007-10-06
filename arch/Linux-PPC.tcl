
TCLDIR=/Projects/namd2/tcl/linux-ppc
TCLINCL=-I$(TCLDIR)/include -I$(HOME)/tcl/include
TCLLIB=-L$(TCLDIR)/lib -L$(HOME)/tcl/lib -lnamdtcl8.3
TCLFLAGS=-DNAMD_TCL -DUSE_NON_CONST
TCL=$(TCLINCL) $(TCLFLAGS)

