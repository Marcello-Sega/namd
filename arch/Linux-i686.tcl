
TCLDIR=/Projects/namd2/tcl/linux
TCLINCL=-I$(TCLDIR)/include -I$(HOME)/tcl-linux/include
TCLLIB=-L$(TCLDIR)/lib -L$(HOME)/tcl-linux/lib -ltcl8.3
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

