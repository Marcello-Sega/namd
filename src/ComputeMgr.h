/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef COMPUTEMGR_H
#define COMPUTEMGR_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"
#include "main.h"

#include "NamdTypes.h"
#include "BOCgroup.h"

#include "Templates/ResizeArray.h"
#include "Compute.h"

class ComputeMgr : public groupmember
{
public:

  ComputeMgr(InitMsg *);
  ~ComputeMgr();
  void createComputes(ComputeMap *map);
  void enqueueWork(Compute *);

private:

  class ComputeElem {
  public:
    ComputeID   cid;
    Compute *c;

    void * operator new(size_t size) { return ::operator new(size); }
    void operator delete(void* ptr) { ::operator delete(ptr); }

    operator<(ComputeElem e) { return (cid < e.cid); }
    operator==(ComputeElem e) { return (cid == e.cid); }

    ComputeElem(ComputeID id=-1, Compute *compute=NULL) : 
      cid(id), c(compute) {};
    ~ComputeElem() { };
    ComputeElem& operator=(const ComputeElem &e) { 
      cid = e.cid; c = e.c;  // Do not delete c!  This op only used to shuffle
                             // we delete the c here only when the Compute is 
		             // moved off!
      return(*this);
    };
  };

  int numComputes;

  typedef ResizeArray<ComputeElem> ComputeList;
  typedef ResizeArray<int> ComputeIndex;

  // global patch number to local patch table conversion table
  ComputeIndex computeIndex;

  // an array of compute pointers residing on this node
  ComputeList computeList;

  int workDistribGroup;
};

#endif /* COMPUTEMGR_H */
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeMgr.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.3 $	$Date: 1996/12/16 23:16:56 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeMgr.h,v $
 * Revision 1.3  1996/12/16 23:16:56  jim
 * eliminated warning about multiple new's
 *
 * Revision 1.2  1996/11/30 00:34:05  jim
 * added some includes, now uses InitMsg
 *
 * Revision 1.1  1996/11/27 20:19:59  jim
 * Initial revision
 *
 *
 ***************************************************************************/
