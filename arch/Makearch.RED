##########################################
##########################################
# Makefile architecture HPUX10
##########################################
##########################################

#####
## Choose your favorite Compiler options
####
CC = cicc
CFLAGS = -O
CXX = ciCC
CXXOPTS = -O
CXXTHREADOPTS = $(CXXOPTS)

#####
# Flags
#####
# uncomment for no warnings during compile
#NOWARN=-w
# uncomment for purify during compile
# PURIFY=-purify

DFLAGS = -Ae -Wp,-H500000
PVM_ARCH = HPPA

#####
# definitions for Charm routines
#####
CHARMARCH = paragon-red
CHARM = $(CHARMBASE)/$(CHARMARCH)

######
## definitions for TCL interface
######

TCLDIR=/usr/local
TCLINCL=-I$(TCLDIR)/include
TCLLIB=-L$(TCLDIR)/lib -ltclx -ltcl
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

######
## definitions for MDCOMM interface
######

#MDCOMMDIR=/home/rkufrin/rapp
#MDCOMMINCL=-I${MDCOMMDIR}/include
#MDCOMMLIB=-L$(MDCOMMDIR)/lib/HPPA -lrapp
#MDCOMMFLAGS=-DMDCOMM
#MDCOMM=$(MDCOMMINCL) $(MDCOMMFLAGS)

#####
# Unix commands
#####
ECHO = echo
MOVE = mv
COPY = cp
RM = rm -f
AR = xar
