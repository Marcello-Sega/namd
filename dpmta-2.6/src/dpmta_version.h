/*
*  dpmta_version.h - version number for distributed pmta implementation
*    declarations
*
*  w. t. rankin
*
*  Copyright (c) 1994 Duke University
*  All rights reserved
*
*/

/* $Id: dpmta_version.h,v 1.1 1997/09/05 19:42:18 jim Exp $
 *
 * revision history:
 *
 * $Log: dpmta_version.h,v $
 * Revision 1.1  1997/09/05 19:42:18  jim
 * Original distribution.
 *
 * Revision 2.14  1997/05/09  20:18:44  wrankin
 * new release version 2.6 to include parallel-piped code
 *
 * Revision 2.13  1997/04/28  15:58:35  wrankin
 * update for latest release
 *
 * Revision 2.12  1997/04/08  19:11:21  wrankin
 * Fixed syntax problems in serial and LJ code.
 *
 * Revision 2.11  1997/03/26  20:48:39  wrankin
 * new release -
 *  - serial version of dpmta
 *  - memory leak fixed in interaction list generation
 *  - improved timing code
 *
 * Revision 2.10  1997/02/07  18:15:23  wrankin
 * updates to docs and version file for new release (2.5.1)
 *
 * Revision 2.9  1996/11/01  02:26:50  wrankin
 * modifications to support cray t3d compilation
 * version update for 2.5 release
 *
 * Revision 2.8  1996/02/06  00:55:06  wrankin
 * update for version 2.3 release.
 *
 * Revision 2.7  1995/06/27  14:20:35  wrankin
 * Implementation of relative interaction list algorithm.
 * Interaction list generation code removed from master library.
 * Cell Table restructured as a 2-D table indexed by level/cell
 * All direct computations are now performed in a single pass.
 * General cleanup and removal of old files.
 *
 * Revision 2.6  1995/06/05  19:39:17  wrankin
 * to support release of version 2.1.2
 *
 * Revision 2.5  1995/04/09  22:22:49  wrankin
 * additional functionality and maintenance release
 * - slaves now generate the inverse interaction list
 * - T3D library functions implemented
 *
 * Revision 2.4  1995/02/24  21:24:13  wrankin
 * update for new version 2.1.0
 *
 * Revision 2.3  1995/01/13  17:08:07  wrankin
 * updated release number for new release (2.0.4)
 *   - dynamic particle allocation implemented
 *   - general cleanup of code
 *
 * Revision 2.2  1995/01/01  15:00:37  wrankin
 * updates for fixes at UIUC, new Makefile
 *
 * Revision 2.1  1994/12/10  16:41:49  wrankin
 * updated version for new release
 *
 * Revision 1.2  1994/12/06  19:56:03  wrankin
 * first set of release patches.
 *
 * Revision 1.1  1994/12/06  19:55:09  wrankin
 * Initial revision
 *
 *
 *
*/

#define DPMTA_VERSION 2.6
