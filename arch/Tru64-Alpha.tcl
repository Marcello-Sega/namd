
TCLDIR=/Projects/namd2/tcl/alpha
TCLINCL=-I$(TCLDIR)/include -I$(HOME)/tcl/include
TCLLIB=-L$(TCLDIR)/lib -L$(HOME)/tcl/lib -ltcl8.3
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

