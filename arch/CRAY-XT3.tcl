TCLDIR=$(HOME)/tcl/cray-xt3
TCLINCL=-I$(TCLDIR)/include -I$(HOME)/tcl/include
TCLLIB=-L$(TCLDIR)/lib -L$(HOME)/tcl/lib -ltcl8.3 -ldl
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)
