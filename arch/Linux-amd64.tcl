
TCLDIR=/Projects/namd2/tcl/linux-amd64
TCLINCL=-I$(TCLDIR)/include -I$(HOME)/tcl/include
TCLLIB=-L$(TCLDIR)/lib -L$(HOME)/tcl/lib -ltcl8.3 -ldl
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

