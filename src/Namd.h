//-*-c++-*-
/***************************************************************************/
/*              (C) Copyright 1996,1997 The Board of Trustees of the       */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION: Used with single instance owned by main()
 *
 ***************************************************************************/

#ifndef _NAMD_H
#define _NAMD_H

#include "charm++.h"
#include "converse.h"

class Node;
class NamdState;

class Namd {
public:
  // Constructor for startup node
  Namd(void);             

  // Shut down slaves and clean up
  ~Namd();                

  // read in various input files by invoking
  // proper classes (parameters, molecule etc)
  void startup(char *);   

  // last call of system
  static void namdDone(void);

  // Emergency bailout 
  static void die() { CmiAbort("NAMD ABORTING DUE TO Namd::die().\n"); }
  static void startTimer() {
	cmiWallStart = CmiWallTimer(); 
	cmiCpuStart = CmiCpuTimer();
	if ( ! cmiFirstStart ) {
		cmiFirstStart = 1;
		cmiWallFirstStart = cmiWallStart;
		cmiCpuFirstStart = cmiCpuStart;
	}
  }

  static float cmiWallStart;
  static float cmiCpuStart;

private:
  static int cmiFirstStart;
  static float cmiWallFirstStart;
  static float cmiCpuFirstStart;

  Node *node;
  int nodeGroup;
  int workDistribGroup;
  int patchMgrGroup;

  NamdState *namdState;

};

#endif /* _NAMD_H */

