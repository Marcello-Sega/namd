/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************/
/* DESCRIPTION:							           */	
/*									   */
/***************************************************************************/

#ifndef _NODE_H
#define _NODE_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.h"

class Molecule;
class Parameters;
class SimParameters;
class ConfigList;
class WorkDistrib;

class Node : public groupmember
{
private:
  Molecule *molecule;
  Parameters *parameters;
  SimParameters *simParameters;
  ConfigList *configList;

  WorkDistrib *workDistrb;

public:
  // Charm Entry point - distributed contructor
  Node(NodeInitMsg *msg);
  ~Node(void);

  int myid(void);		   
  
  // Charm Entry point - Read in system data, get all ready to simulate
  void startup(InitMsg *initmsg);  

  // Run for the number of steps specified in the sim_parameters
  void run(void);                  
};

#endif /* _NODE_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 1996/08/16 21:19:34 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.h,v $
 * Revision 1.1  1996/08/16 21:19:34  ari
 * Initial revision
 *
 *
 ***************************************************************************/
