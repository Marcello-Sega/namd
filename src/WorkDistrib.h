/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************/
/* DESCRIPTION:								   */
/*									   */
/***************************************************************************/

#ifndef _WORKDISTRIB_H
#define _WORKDISTRIB_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "NamdTypes.h"

class Node;
class InitMsg;

const Position sysDimX = 20.;
const Position sysDimY = 20.;
const Position sysDimZ = 20.;
const Position patchSize = 4.;

enum { maxPatchDepends = 126 };

class WorkDistrib : public groupmember
{
private:
  Node *node;
  void mapPatches(void);
  void mapComputes(void);
  void mapAngleComputes(void);
  void mapElectComputes(void);

public:
  WorkDistrib(InitMsg *msg);
  ~WorkDistrib(void);

  void parentNode(Node *inode);

  void buildMapsFromScratch(void); 
};

#endif /* WORKDISTRIB_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: WorkDistrib.h,v $
 *	$Author: brunner $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/08/16 21:41:11 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: WorkDistrib.h,v $
 * Revision 1.3  1996/08/16 21:41:11  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 21:16:04  ari
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 20:34:23  brunner
 * Initial revision
 *
 * Revision 1.6  1996/08/03 20:08:09  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/07/01 16:25:30  brunner
 * *** empty log message ***
 *
 * Revision 1.4  1996/06/14 18:40:54  brunner
 * *** empty log message ***
 *
 * Revision 1.3  1996/06/11 22:38:04  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/06/04 16:12:18  brunner
 * Create dummy patches
 *
 * Revision 1.1  1996/05/30 20:16:09  brunner
 * Initial revision
 *
 * Revision 1.1  1996/05/30 20:11:09  brunner
 * Initial revision
 *
 ***************************************************************************/
