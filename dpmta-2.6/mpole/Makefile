# 
# Makefile to build multipole library
#
# w. t. rankin
#
# Copyright (c) 1995 Duke University
# All rights reserved
#
# RCS Info: $Id: Makefile,v 1.4 1998/08/23 20:59:22 jim Exp $
#
# Revision History:
# $Log: Makefile,v $
# Revision 1.4  1998/08/23 20:59:22  jim
# Fixed compile option passing.
#
# Revision 1.3  1997/09/28 10:24:03  milind
# Fixed Makefiles to include arch params from a single file.
#
# Revision 1.2  1997/09/12 22:56:21  jim
# Modifications to work with converse pvm.
#
# Revision 1.1  1997/09/05 19:41:26  jim
# Original distribution.
#
# Revision 1.9  1997/05/13  17:52:11  wrankin
# fixed irritation bug seen when math.h is not included
# added a missing prototype to mpe.h
#
# Revision 1.8  1997/05/09  20:14:55  wrankin
# added routines to de-allocate global arrays created in [C,LJ]init()
# added routines to free up multipole expansion matrices
# added LJ prototypes and ansi-fied more procedures
#
# Revision 1.7  1996/10/18  17:00:51  wrankin
# implement new heirarchical make structure
#
# Revision 1.6  1996/09/24  19:05:48  wrankin
# added more test codes.
#
# Revision 1.5  1996/09/19  18:38:02  wrankin
# Makefile now leaves libmpole in current directory
#
# Revision 1.4  1995/09/15  15:08:51  wrankin
# Added Lennard-Jones multipole routines to library.
#
# Revision 1.3  1995/07/19  01:55:14  wrankin
# updates to reflect move of multipole library code to a new directory
#
# Revision 1.2  1995/07/17 01:31:43  wrankin
# Makefile updates to support SGI Power Challenge platform
#
# Revision 1.1.1.1  1995/07/10  13:11:46  wrankin
# Initial release of the Multipole Library
# Based upon W. Elliott's MDMA codes
# Implements Colomb Potentials only (no LJ-potentials yet)
#
#
################################################################

include ../../Makearch

#  for the HPPA machines, add "-Ae" and "-Wp,-H500000" to CFLAGS
#  for the SGI machines, add "-Woff all" to CFLAGS

CFLAGS = $(COPTS) $(DFLAGS) -I$(CHARM)/include

SRC  = mpe_fft.c mpe_misc.c mpe_mpoleC.c mpe_allocC.c \
       mpe_mpoleLJ.c mpe_allocLJ.c
OBJS = $(SRC:.c=.o)

TSRC = mpe_testC.c mpe_testLJ.c
TOBJS = $(TSRC:.c=.o)

################################################################

default : all

all : libmpole.a

libmpole.a : $(OBJS)
	$(RM) libmpole.a
	$(AR) rs libmpole.a $(OBJS)

install : all

clean:
	$(RM) libmpole.a $(OBJS) $(TOBJS)

testLJ : mpe_testLJ.o libmpole.a
	$(CC) -o mpe_testLJ mpe_testLJ.o -L. -lmpole -lm

testC : mpe_testC.o libmpole.a
	$(CC) -o mpe_testC mpe_testC.o -L. -lmpole -lm

################################################################
#
#  include file dependencies
#

mpe_fft.c : mpe.h mpe_fftC.h
mpe_misc.c : mpe.h
mpe_mpoleC.c : mpe.h mpe_legendre.h mpe_fftC.h
mpe_mpoleLJ.c : mpe.h
