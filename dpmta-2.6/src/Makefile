# 
# Makefile to build dpmta_test and dpmta_slave
#
# w. t. rankin
#
# Copyright (c) 1995 Duke University
# All rights reserved
#
# this makefile is set up to compile for a network of sun workstations
# running solaris 2.3.  the compile flags below wil need to be modifed
# for other platforms.
#
#
# RCS Info: $Id: Makefile,v 1.3 1997/09/29 23:58:33 jim Exp $
#
# Revision History:
# $Log: Makefile,v $
# Revision 1.3  1997/09/29 23:58:33  jim
# Incorporated changes from version 2.6.1 of DPMTA.
#   - fixes for bad handling of empty/invalid multipoles when
#     using large processor sets.
#   - moved functions that provide data mapping to processors.  master
#     and slave routines now call the same function in dpmta_distmisc.c
# Also, switched pvmc.h back to pvm3.h.
#
# Revision 1.2  1997/09/28 10:24:07  milind
# Fixed Makefiles to include arch params from a single file.
#
# Revision 1.1  1997/09/05 19:41:47  jim
# Original distribution.
#
# Revision 2.45  1997/05/13  17:53:16  wrankin
# use basic types defined in mpe.h instead of redeclaring them
#
# Revision 2.44  1997/05/12  17:59:03  wrankin
# updated compiler declarations
#
# Revision 2.43  1997/02/26  20:43:23  wrankin
# updated timing measurements and placed routines in single module
#
# Revision 2.42  1997/02/26  16:54:23  wrankin
# added support for serial non-pvm compilation of DPMTA
#
# Revision 2.41  1997/01/13  22:10:01  wrankin
# added supprt for dpmta_test5.c which does MIMD with PDB file input
#
# Revision 2.40  1996/11/01  02:26:41  wrankin
# modifications to support cray t3d compilation
# version update for 2.5 release
#
# Revision 2.39  1996/10/29  19:34:11  wrankin
# changes to makefile to support other platforms
# fix for multi-master communication interface
# new global dumpdata routine to support direct pbc/macro verification
#
# Revision 2.38  1996/10/28  23:01:49  wrankin
# additions to test routines to provide processing paramenters as
#   part of output file.
# dpmta_direct will now perform macroscopic computations, reading
#   in processing parameters from particle position file.
#
# Revision 2.37  1996/10/18  17:04:26  wrankin
# improved PVM message passing structures.
# made current and cleaned up T3D library code
# added additional test codes
# implement new heirarchical make structure
#
# Revision 2.36  1996/09/24  18:41:10  wrankin
# many changes for support of version 2.5
#  - non cubic cells
#  - non power-of-2 processors
#  - fixed the resize code
#  - new virial interface
#  - changes to macroscopic code (still not working)
#  - code cleanup
#  - new test program that reads PDB file
#  - update code for T3D support
#
# Revision 2.35  1996/08/20  17:12:31  wrankin
# Enhancement to allow any number of processes to be utilized.
# Formally the number of precesses was limited to a power of 2.
#
# Revision 2.34  1996/02/29  21:13:13  wrankin
# New relaease: 2.4 (.1)
#    - simplified calling structure for initialization
#    - macroscopic periodic code
#    - fixed PBC calculation - all particles are now stored in the
#      cell as positions relative to the cell center.
#    - virial preasure tensor computed
#    - fix to allow particles on the outer cube boundary to be
#      included. (UNC fix)
#    - fix to order reception of particle data during the distributed
#      calling sequence (UIUC fix)
#    - removed M2L code that didn't use transfer matrices
#    - early hooks in to perform interaction list sorting.
#    - fixed LJ scaling factor for 1/r^12 potential.
#    - cleaned up the LJ interface.
#    - and of course, my continued efforts to ANSI-fy this beast.
#
# Revision 2.33  1996/02/12  15:20:48  wrankin
# Added Macroscopic Assemblies code.
#
# Revision 2.32  1996/02/06  00:55:02  wrankin
# update for version 2.3 release.
#
# Revision 2.31  1995/12/08  23:00:13  wrankin
# preliminary release of DPMTA 2.3
#   - added working Periodic Boundary Conditions (PDC) flag
#   - added code for Virial Computation (not working yet)
#   - general cleanup of some modules
#
# Revision 2.30  1995/11/29  22:28:57  wrankin
# addition of periodic boundary code
# beginning work on virial code.
#
# Revision 2.29  1995/10/02  20:58:34  wrankin
# changes to support creation and use of MPE x-fer matrix for
#   LJ multipole calculation
#
# Revision 2.28  1995/10/01  21:45:39  wrankin
# implemented LJ multipole calculations - compile time option only
#
# Revision 2.27  1995/07/26  01:52:55  wrankin
# updated Makefiles
# slight performance improvement for direct calculations
# clean up test code example
# fixed t3dlib for new multipole code
#
# Revision 2.26  1995/07/17 01:13:31  wrankin
# updates to support SGI Power Challenge (IRIX 6.0.1)
# cleanup of Makefile
# initial work on T3D port (not yet supported in this release)
#
# Revision 2.25  1995/07/10  02:47:08  wrankin
# multipole processing code modified to use precomputed mpe transfer
#   matrices.
#
# in addition, all multippole calculation routines have been removed
#   from this source distribution and will be placed in a separate
#   multipole library to be supplied externally.
#
# Revision 2.24  1995/07/01  03:26:51  wrankin
# initial cut at precomputation of mpe transfer functions.
# this works for vanilla M2L.  it does not work for FFT enhancements.
#
# Revision 2.23  1995/06/27  14:20:12  wrankin
# Implementation of relative interaction list algorithm.
# Interaction list generation code removed from master library.
# Cell Table restructured as a 2-D table indexed by level/cell
# All direct computations are now performed in a single pass.
# General cleanup and removal of old files.
#
# Revision 2.22  1995/05/15  19:33:57  wrankin
# removded dependancy declaration for dpmta_mastilist.c
# added dependancy declaration for dpmta_mastmkil.c
#
# Revision 2.21  1995/04/26  15:48:35  wrankin
# updates for version 2.1.1
#
# Revision 2.20  1995/03/06  08:01:47  wrankin
# added dpmta_slvmkiil.c to source files for dpmta_slave
#
# Revision 2.19  1995/02/24  21:53:11  wrankin
# builds with new interaction list generator (dpmta_mastmkil.c)
# "make install" will now rebuild sources, if needed.
# "make clean" will remove libdpmta.a
#
# Revision 2.18  1995/01/13  17:15:19  wrankin
# added notes for HPPA compilation
# removed automatic installation step
#
# Revision 2.17  1995/01/02  19:50:06  wrankin
# set library definitions for SUN4SOL2 compilation
#
# Revision 2.16  1995/01/01  14:56:59  wrankin
# changed ARCH to PVM_ARCH in install section
#
# Revision 2.15  1995/01/01  14:46:37  wrankin
# removed automatic compilation by platform
# user must now edit Makefile for site-specific information
# removed cray-t3d compilation entirely
#
# Revision 2.14  1994/11/30  18:46:41  wrankin
# fixed declarations for HP ATM compilation
#
# Revision 2.13  1994/11/30  17:52:53  wrankin
# added dpmta.h dependancies
#
# Revision 2.12  1994/11/30  16:07:25  wrankin
# added support for additional test files
#
# Revision 2.11  1994/11/28  04:32:58  wrankin
# compilation now produces master interface library (libdpmta.a) and
#   uses thew library when building test example.
# added include file dependancies
# removed useless libraries from compilation lines
#
# Revision 2.10  1994/10/26  03:12:18  wrankin
# updated source file lists
# dpmta_master has been replaced by dpmta_test
# removed the 'mil' tool, no longer supported
#
# Revision 2.9  1994/10/19  00:27:41  wrankin
# added O2 optimization to compile line
#
# Revision 2.8  1994/10/14  04:25:26  wrankin
# modified file lists for new multipole library
#
# Revision 2.7  1994/10/09  18:47:52  wrankin
# added module dpmta_fft.c to compile list
# added DATADUMP for default solaris build
#
# Revision 2.6  1994/08/25  19:17:27  wrankin
# removed T3DBUG flag from all but "cray" make flags
#
# Revision 2.5  1994/08/25  06:09:19  wrankin
# removed T3DBUG flag from Solaris compilation
# added TIME flag to HP compilation
# default Make automatically installs executables
#
# Revision 2.4  1994/08/11  18:22:34  wrankin
# cleaned up compilation architecture options
# we don't need to edit the makefile anymore
#
# Revision 2.3  1994/06/18  06:15:53  wrankin
# support for optional ilist generation if no file specified
#
# Revision 2.2  1994/06/16  07:37:31  wrankin
# Added -lgpvm3 library support needed for PVM 3.3
# Added compilation of mil (Make Interaction List) program
#
# Revision 2.1  1994/06/02  12:52:27  wrankin
# No change.  Update version number for Release 2
#
# Revision 1.9  1994/06/02  03:40:43  wrankin
# added optimizer switched to compiler command line
#
# Revision 1.8  1994/05/27  15:27:01  wrankin
# added settings for sun4sol2
#
# Revision 1.7  1994/05/27  07:59:27  wrankin
# added recursive make to allow for different CFLAGS defs
#  for master and slave
# cray and linux flags set, still need to add sun4sol2
#
# Revision 1.6  1994/03/15  16:24:07  wrankin
# update to add new modules
# added header file dependancies
#
# Revision 1.5  1994/02/22  03:58:04  wrankin
# added dpmta_slvpart.c compilation
#
# Revision 1.4  1994/02/16  05:01:51  wrankin
# added interaction list modules
#
# Revision 1.3  1994/02/09  21:32:38  wrankin
# include particle generation code into master compile
#
# Revision 1.2  1994/02/09  17:38:03  wrankin
# changes for new filenames
#
# Revision 1.1  1994/02/01  19:32:31  wrankin
# Initial revision
#
################################################################

include ../../Makearch

# tell the compiler where to locate the PVM include
# files and the mpole include file

IFLAGS = -I$(PVM_ROOT)/include/ -I../mpole/


# library flags - set the flags to point to the directories
#   containing the PVM and DPMTA libraries.

TLDIR = -L$(PVM_ROOT)/lib/$(PVM_ARCH) -L.
SLDIR = -L$(PVM_ROOT)/lib/$(PVM_ARCH) -L. -L../mpole

#
# libraries to add
#
TLIBS = -lpvm3 -lgpvm3 -ldpmta -lm $(OSLIBS)
SLIBS = -lpvm3 -lgpvm3 -lmpole -lm $(OSLIBS)



################################################################
#
#  nothing below this line should need to be modified
#
################################################################

CFLAGS = -O3 -D$(PVM_ARCH) $(DFLAGS) $(IFLAGS)

MLFLAGS = 
SLFLAGS = $(SLDIR)
TLFLAGS = $(TLDIR)


MSRC  = dpmta_mastiface.c dpmta_mastinit.c dpmta_mastiter.c dpmta_timer.c \
	dpmta_distmisc.c
MOBJS = $(MSRC:.c=.o)

SSRC  = dpmta_slave.c dpmta_slvmkcell.c dpmta_slvmcalc.c dpmta_slvpcalc.c \
        dpmta_slvcomm.c  dpmta_slvclean.c dpmta_slvstart.c \
	dpmta_slvmkiil.c  dpmta_slvmkil.c dpmta_slvmkhl.c \
	dpmta_slvcompute.c dpmta_slvmacro.c dpmta_slvscale.c \
	dpmta_timer.c \
	dpmta_distmisc.c
SOBJS = $(SSRC:.c=.o)


TSRC  = dpmta_test.c dpmta_test2.c dpmta_test2b.c dpmta_test3.c \
        dpmta_test4.c dpmta_test5.c dpmta_direct.c dpmta_error.c 
TSRC2 = dpmta_testdump.c
TOBJS = $(TSRC:.c=.o) $(TSRC2:.c=.o)
TEXEC = $(TSRC:.c=)

MINCS = dpmta_pvm.h dpmta_cell.h dpmta.h dpmta_version.h dpmta_timer.h
SINCS = dpmta_pvm.h dpmta_cell.h dpmta_slave.h dpmta_version.h dpmta_timer.h

################################################################

all:
	$(MAKE) dpmta_slave CFLAGS="$(CFLAGS)" LDLIBS="$(SLIBS)" \
                LDFLAGS="$(SLFLAGS)" CC=$(CC)
	$(MAKE) libdpmta.a CFLAGS="$(CFLAGS)" CC=$(CC)
	$(MAKE) dpmta_test CFLAGS="$(CFLAGS)" LDLIBS="$(TLIBS)" \
                LDFLAGS="$(TLFLAGS)" CC=$(CC)
	$(MAKE) dpmta_test2 CFLAGS="$(CFLAGS)" LDLIBS="$(TLIBS)" \
                LDFLAGS="$(TLFLAGS)" CC=$(CC)
	$(MAKE) dpmta_test2b CFLAGS="$(CFLAGS)" LDLIBS="$(TLIBS)" \
                LDFLAGS="$(TLFLAGS)" CC=$(CC)
	$(MAKE) dpmta_test3 CFLAGS="$(CFLAGS)" LDLIBS="$(TLIBS)" \
                LDFLAGS="$(TLFLAGS)" CC=$(CC)
	$(MAKE) dpmta_test4 CFLAGS="$(CFLAGS)" LDLIBS="$(TLIBS)" \
                LDFLAGS="$(TLFLAGS)" CC=$(CC)
	$(MAKE) dpmta_test5 CFLAGS="$(CFLAGS)" LDLIBS="$(TLIBS)" \
                LDFLAGS="$(TLFLAGS)" CC=$(CC)
	$(MAKE) dpmta_direct CFLAGS="$(CFLAGS)" LDLIBS="$(TLIBS)" \
                LDFLAGS="$(TLFLAGS)" CC=$(CC)
	$(MAKE) dpmta_error CFLAGS="$(CFLAGS)" LDLIBS="$(TLIBS)" \
                LDFLAGS="$(TLFLAGS)" CC=$(CC)



dpmta_slave : $(SOBJS) ../mpole/libmpole.a

#	$(CC) $(LDFLAGS) -o $@ $? $(LDLIBS)

libdpmta.a : $(MOBJS)
	$(RM) libdpmta.a
	$(AR) r libdpmta.a $(MOBJS)

dpmta_test : dpmta_test.o dpmta_testdump.o libdpmta.a

dpmta_test2 : dpmta_test2.o dpmta_testdump.o libdpmta.a

dpmta_test2b : dpmta_test2b.o libdpmta.a

dpmta_test3 : dpmta_test3.o dpmta_testdump.o libdpmta.a

dpmta_test4 : dpmta_test4.o dpmta_testdump.o libdpmta.a

dpmta_test5 : dpmta_test5.o dpmta_testdump.o libdpmta.a

dpmta_direct : dpmta_direct.o

dpmta_error : dpmta_error.o

install: 
	cp dpmta_slave $(PVM_HOME)/bin/$(PVM_ARCH)
	cp $(TEXEC) $(PVM_HOME)/bin/$(PVM_ARCH)
	cp libdpmta.a $(PVM_HOME)/lib/$(PVM_ARCH)

clean:
	rm -f libdpmta.a dpmta_slave $(MOBJS) $(SOBJS) $(TOBJS) $(TEXEC)
	rm -f libdpmta2.a dpmta_t3dlib.o
	rm -f libspmta.a dpmta_serial.o


################################################################
#
#  include file dependencies
#

dpmta_test.c : dpmta.h

dpmta_test2.c : dpmta.h

dpmta_test2b.c : dpmta.h

dpmta_test3.c : dpmta.h 

dpmta_test4.c : dpmta.h 

dpmta_direct.c : dpmta.h 

dpmta_mastiface.c : dpmta.h dpmta_pvm.h

dpmta_mastinit.c : dpmta.h dpmta_pvm.h 

dpmta_mastiter.c : dpmta.h dpmta_pvm.h 

dpmta_slave.c : dpmta_cell.h dpmta_pvm.h

dpmta_slvmkcell.c : dpmta_cell.h dpmta_pvm.h dpmta_slave.h

dpmta_slvmcalc.c :  dpmta_cell.h dpmta_pvm.h dpmta_slave.h dpmta_slvmkil.h

dpmta_slvpcalc.c : dpmta_cell.h dpmta_pvm.h dpmta_slave.h

dpmta_slvmisc.c : dpmta_cell.h dpmta_pvm.h dpmta_slave.h

dpmta_slvcomm.c : dpmta_cell.h dpmta_pvm.h dpmta_slave.h

dpmta_slvclean.c : dpmta_cell.h dpmta_slave.h

dpmta_slvstart.c : dpmta_cell.h dpmta_slave.h

dpmta_slvmkhl.c : dpmta_cell.h dpmta_slave.h dpmta_slvmkil.h

dpmta_slvmkil.c : dpmta_cell.h dpmta_slave.h 

dpmta_slvmkiil.c : dpmta_cell.h dpmta_pvm.h dpmta_slave.h

dpmta_slvcompute.c : dpmta_cell.h dpmta_slave.h

dpmta_slvscale : dpmta_cell.h dpmta_slave.h 
