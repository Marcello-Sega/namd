#**************************************************************
# This is the makefile to build dpme_test and other support libs
# Nour Toukmaji/Dan Paul - Duke University 1995 Copyright.
#**************************************************************

MAKE   = make
SHELL  = /bin/sh
CC     =gcc  -Wall # SGI: cc -fullwarn -g3 -lmalloc
AR     = ar
RM     = rm -f

# specify PVM environmental variables.  PVM_HOME is the location
# where the executables will be placed.
#
# SUN
#PVM_ARCH = SUN4SOL2       #DSIL
PVM_ROOT = ../pvm3
# SGI 
#PVM_ARCH = SGI64
#PVM_ROOT = /usr/opt/pvm/pvm3
#don't forget to change the LIBS

#PVM_HOME = $(HOME)/pvm3
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
DFLAGS = -D$(PVM_ARCH)  -DTIMEME   -DVIRIAL -DVERBOSE  
# tell the compiler where to locate the PVM include
# files
#
IFLAGS = -I$(PVM_HOME)
CFLAGS = $(IFLAGS)

# library flags - set the flags to point to the directories
#   containing the PVM 
#
LDIR = -L. -L../pvm3

#  tell the linker where to which libraries to link in.
#  this will vary by operating system
#
# solaris :
# LIBS = -lpvmc -ldpme2 -lsocket -lnsl  -lm
# otherwise :
LIBS = -lpvmc -ldpme2 -lm 
###################################################################
#        END  OF SYSTEM CONFIGURATION SECTION   
###################################################################

####################################################################
DEBUG= 
MYCFLAGS=   $(DFLAGS) $(IFLAGS) $(DEBUG) #-O3
LFLAGS = $(LDIR) $(LIBS)
#FFLAGS= -O +U77 -lf -lm -lF77
MYCC = $(CC) $(MYCFLAGS)
####################################################################
testobjects = dpme2_test.o
objects = dpme2_recip_sum.o dpme2_paralibs.o erfcfun.o dpme2_utility.o \
	  bspline2.o charge_grid2.o fftcalls.o pmesh_kspace2.o \
	  cfftb.o cfftb1.o cfftf.o cfftf1.o cffti.o cffti1.o \
	  passb.o passb2.o passb3.o passb4.o passb5.o \
	  passf.o passf2.o passf3.o passf4.o passf5.o \
	  pub3dfft.o
headers = dpme2.h dpme2_pvm.h dpme2def.h f2c.h protype.h
libobjects =  libdpme2.a

all:	libdpme2.a

.c.o:	$(headers)
	$(CC) $(CFLAGS) -c $*.c
	ar r libdpme2.a $*.o

libdpme2.a:	$(objects)

dpme_test:    $(objects)  $(testobjects) 
	$(CC) $(CFLAGS) -o dpme2_test $(objects) $(testobjects) $(LDLIBS) $(LFLAGS)

clean:
	-rm *.o *.a

