# config for ASCI Red

# compiler options
CXX = ciCC -DPLACEMENT_NEW
CXXOPTS = -O
CXXTHREADOPTS = $(CXXOPTS)
CXXSIMPARAMOPTS = $(CXXOPTS)
CC = cicc
COPTS = -O

# definitions for Charm routines
CHARMARCH = paragon-red
CHARM = $(CHARMBASE)/$(CHARMARCH)

# definitions for TCL interface
#TCLDIR=/usr/local
#TCLINCL=-I$(TCLDIR)/include
#TCLLIB=-L$(TCLDIR)/lib -ltcl
#TCLFLAGS=-DNAMD_TCL
#TCL=$(TCLINCL) $(TCLFLAGS)

# definitions for FFT (used by PME)
#FFTDIR=/usr/local
#FFTINCL=-I$(FFTDIR)/include
#FFTLIB=-L$(FFTDIR)/lib -lrfftw -lfftw
#FFTFLAGS=-DNAMD_FFTW
#FFT=$(FFTINCL) $(FFTFLAGS)

# definitions for CCS (used by BioCoRE)
#CCS=-DNAMD_CCS

