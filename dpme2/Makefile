include ../Make.charm
include ../Makearch

#**************************************************************
# This is the makefile to build dpme_test and other support libs
# Nour Toukmaji/Dan Paul - Duke University 1995 Copyright.
#**************************************************************

MAKE   = make
SHELL  = /bin/sh
# CC     = gcc # -Wall  sun :,sgi: cc -fullwarn -g3 
CC     = $(CXX)
CFLAGS = $(CXXOPTS)
AR     = ar
RM     = rm -f

# specify PVM environmental variables.  PVM_HOME is the location
# where the executables will be placed.
#
# SUN
#PVM_ARCH = SUN4SOL2       #DSIL
# SGI 
#PVM_ARCH = SGI64
#don't forget to change the LIBS

#PVM_HOME = $(HOME)/pvm3
#PVM_ARCH = HPUX
PVM_HOME = ../pvm3

# Specify compile time flags
# -DTIMEME: activates timing of different sections of program to measure performance.
# -DVERBOSE: prints messages that indicates the status of a run.
# -DVIRIAL: activates virial calculation.
#  -DDPME_DEBUG(2,3): print debug particle data; memory alloc; pvm prints to screen

#  suggested: for the HPPA machines, add "-Ae" 
#  for the SGI machines, add "-Woff all" 
#

# Must set the WORKINGDIR and DATADIR in dpme_pvm.h

#DFLAGS = -D$(PVM_ARCH) -DVIRIAL -DTIMEME -DVERBOSE -DDPME_DEBUG(2,3)
#DFLAGS = -D$(PVM_ARCH)  -DTIMEME  -DVERBOSE -DVIRIAL
DFLAGS = -D$(PVM_ARCH)  -DVIRIAL
# tell the compiler where to locate the PVM include
# files
#
IFLAGS = -I$(PVM_HOME) -I../src -I../dpmta2 -I$(CHARM)/include

# library flags - set the flags to point to the directories
#   containing the PVM 
#
LDIR = -L. -L$(PVMDIR)

#  tell the linker where to which libraries to link in.
#  this will vary by operating system
#
# solaris :
LIBS = -lpvmc -lm 

# otherwise :
#LIBS = -lpvm3 -lgpvm3  -lm
###################################################################
#        END  OF SYSTEM CONFIGURATION SECTION   
###################################################################

####################################################################
DEBUG= 
MYCFLAGS=   $(DFLAGS) $(IFLAGS) $(DEBUG) #-O3
LFLAGS = $(LDIR)
#FFLAGS= -O +U77 -lf -lm -lF77
MYCC = $(CC) $(MYCFLAGS)
####################################################################
objects = dpme2_recip_sumW.o dpme2_paralibsW.o  erfcfun.o  dpme2_utilityW.o
Dpme_Recip = charge_grid2.o fftcalls.o bspline2.o pmesh_kspace2.o
Dpubfft = cfftb.o cfftf1.o passb.o passb4.o passf2.o \
	  passf5.o cfftb1.o cffti.o passb2.o passb5.o passf3.o \
	  cfftf.o cffti1.o  passb3.o  passf.o passf4.o pub3dfft.o

all:	$(objects) $(Dpme_Recip) $(Dpubfft)
	ar -vr libdpme.a $(objects) $(Dpme_Recip) $(Dpubfft)

.c.o:
	$(CC) $(CFLAGS) $(DFLAGS) $(IFLAGS) -c -o $*.o  $*.c

dpme_test:    $(objects) dpme2_testW.o
	$(CC) $(CFLAGS)    -o dpme_test  $(objects) dpme2_testW.o $(libobjects) $(LDLIBS) $(LFLAGS)

clean:
	rm -f *.o *.a dpme_test

