/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: MessageManager.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/12/06 19:54:12 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *	MessageManager is used to maintain queues of received messages.
 * We only need this .C file because the HP CC compiler won't inline
 * for loops, thus the contructor must be in a .C.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: MessageManager.C,v $
 * Revision 1.1  1996/12/06 19:54:12  ari
 * Initial revision
 *
 * Revision 1.3  1995/10/09 04:50:56  hazen
 * Updated memory allocation to use C++ new/delete
 *
 * Revision 1.2  1995/03/20  13:50:54  nelson
 * Fixed bug in number of tags calculations
 *
 * Revision 1.1  95/03/18  13:06:13  13:06:13  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/MessageManager.C,v 1.1 1996/12/06 19:54:12 ari Exp $";

#include "MessageManager.h"
#include "common.h"
typedef MessageQueue *MessageQueuePtr;

MessageManager::MessageManager()
{
      int i;		//  Loop counter

      numMessages = 0;
      numTags = (MAXTAGVALUE-MINTAGVALUE) + 1;  //  Figure out number of tags
      queues  = new MessageQueuePtr[numTags];
      if (queues == NULL)
      {
	NAMD_die("memory allocation failed in MessageManager::MessageManager");
      }

      //  Loop through and create the necessary MessageQueues
      for (i=0; i<numTags; i++)
      {
	queues[i] = new MessageQueue();

	if (queues[i] == NULL)
	{
		NAMD_die("memory allocation failed in MessageManager::MessageManager");
	}
      }
}
