/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*                                                                         */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

/***************************************************************************
 The order of execution is expected to be:
   (mode 0) 1. register() and subscribe()
   ------------------ (processing barrier)
   (mode 1) 2. submit() and request()
   ------------------
   (mode 2) 3. unregister() and unsubscribe()
 Doing this out-of-order will cause errors.
 ***************************************************************************/

static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/ReductionMgr.C,v 1.1 1996/12/18 19:53:39 nealk Exp $";

#include <stdlib.h>
#include <stdio.h>

#include "chare.h"
#include "ckdefs.h"
#include "c++interface.h"

#include "InfoStream.h"
#include "ReductionMgr.h"
#define DEBUGM
#include "Debug.h"

/*******************************************
 ReductionMgr::ReductionMgr(): init object
 *******************************************/
ReductionMgr::ReductionMgr()
  {
    #if PANIC > 0
    panicMode = 0;
    #endif
    data = NULL;
    for(i=0; i<REDUCTION_MAX_RESERVED; i++)
    {
    	numSubscribed[i] = 0;
	maxData[i] = 0;
    }
    return *this;
  }

/*******************************************
 ReductionMgr::~ReductionMgr(): free object
 *******************************************/
ReductionMgr::~ReductionMgr()
  {
    ReductionMgrData *nextdata;
    while(data != NULL)
    {
      nextdata = data->next;
      delete data;
      data = nextdata;
    }
  }

/*******************************************
 ReductionMgr::createdata(): create a blank
 data element.
 *******************************************/
ReductionMgrData *ReductionMgr::createdata(int seq)
{
  ReductionMgrData *data;
  data = new ReductionMgrData;
  data->sequence = seq;
  data->dataToSend = REDUCTION_MAX_RESERVED;
  data->next = NULL;
  for(int i=0; i<REDUCTION_MAX_RESERVED; i++)
  {
      data->numData[i] = maxData[i];
      data->tagData[i] = 0;
  }
  return(data);
} /* ReductionMgr::createdata() */

/*******************************************
  ReductionMgr::register(): increase counter
  to a reduction tag.
  ASSUMPTION: this function will only be called
  before data has been depositied.
  (un)register to submit data for reduction
  may cause an error if reductions are active
 *******************************************/
void	ReductionMgr::register(ReductionTag tag)
{
  #if PANIC > 0
  if (panicMode != 0)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die(""Panic due to wrong mode");
  }
  panicMode = 0;
  #endif

  maxData[tag]++;
} /* ReductionMgr::register() */

/*******************************************
 ReductionMgr::unregister(): 
 of counters.
 ASSUMPTION: this function will only be called
 after data has been depositied.
 (un)register to submit data for reduction
 may cause an error if reductions are active
 *******************************************/
void	ReductionMgr::unregister(ReductionTag tag)
{
  #if PANIC > 0
  if (panicMode < 1)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die(""Panic due to wrong mode");
  }
  panicMode = 2;
  #endif

  maxData[tag]--;
} /* ReductionMgr::unregister() */

/*******************************************
 ReductionMgr::remove(): remove a sequence
 of counters.
 ASSUMPTION: this function will only be called
 after data has been depositied.
 (un)register to submit data for reduction
 may cause an error if reductions are active
 *******************************************/
void	ReductionMgr::remove(int seq)
{
  if (!data) return;
  ReductionMgrData *currentdata = data;
  ReductionMgrData *previousdata = NULL;

  for(int i=data->sequence; i<seq; i++)
  {
	#if PANIC > 0
	if (!currentdata)
	  {
	  iout << iERRORF << "Oh no, Mr. Bill!  Cannot find " << seq
	       << " to remove!\n" << endi;
	  NAMD_die("Cannot remove sequence.");
	  }
	#endif
	previousdata = currentdata;
	currentdata = currentdata->next;
  }
  #if PANIC > 0
  if (currentdata->sequence != seq)
  {
    iout << iERRORF << "Yikes! Missed sequence " << seq << " and found "
	 << currentdata->sequence << " instead!\n" << endi;
    NAMD_die("Remove without register");
  }
  #endif

  // delete data
  if (!previousdata)
    {
      // head of queue
      data = currentdata->next;
    }
  else
    {
      // body of queue
      previousdata = currentdata->next;
    }
  delete currentdata;
} /* ReductionMgr::remove() */

/*******************************************
 ReductionMgr::find(): find data for reduction.
 If sequence does not exist, then create it.
 *******************************************/
ReductionMgrData *	ReductionMgr::find(int seq)
{
  ReductionMgrData *current=data;
  ReductionMgrData *previous=NULL;

  #if PANIC > 0
  if (current && (seq < current->sequence))
  {
    iout << iERRORF << "Searching for data that has been removed. "
	 << "want " << seq << " but queue starts at " << current->sequence
	 << "\n" << endi;
    NAMD_die("Searching for data that has been removed");
  }
  #endif

  // find the sequence
  while(current && (current->sequence < seq))
  {
    previous = current;
    current = current->next;
  }

  if (!current)
    {
    // no sequences available (create)
    if (!previous)
	{
	data = createdata(seq);
	current = data;
	}
    // it's a middle sequence (insert)
    else
	{
	previous->next = createdata(seq);
	previous->next->next = current;
	current = previous->next;
	}
    }
  // current now contains the sequence
  return(current);
} /* ReductionMgr::find() */

/*******************************************
 ReductionMgr::submit(): submit data for reduction.
 more == 1 signals immediate submission of other data
 *******************************************/
ReductionMgr::submit(int seq, ReductionTag tag, BigReal data, int more=0)
{
  #if PANIC > 0
  if (panicMode > 1)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die(""Panic due to wrong mode");
  }
  panicMode = 1;
  #endif

  ReductionMgrData *current=find(seq);

  // add to tag
  current->tagData[tag] += data;
  current->numData[tag]--;	/* expect 1 less */

  #if PANIC > 0
  if (current->numData[tag] < 0)
  {
    iout << iERRORF << "Too many submits!\n" << endi;
    NAMD_die("Too many submits");
  }
  #endif

  if (current->numData == 0)
  {
    // all done here!
    // send the data[tag] to main collector

    // one less data to send (delete if all done)
    current->dataToSend--;
    if ((CMyPe() != 0) && (current->dataToSend == 0))
	remove(seq);
  }
} /* ReductionMgr::submit() */

/*******************************************
 ReductionMgr::submit(): pass on submitting data
 *******************************************/
ReductionMgr::submit(int seq, ReductionTag tag, int more=0)
{
  #if PANIC > 0
  if (panicMode > 1)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die(""Panic due to wrong mode");
  }
  panicMode = 1;
  #endif
} /* ReductionMgr::submit() */

/*******************************************
 ReductionMgr::require(): get data.
 Note: this suspends until this data is ready
 and should be called only from Sequencer thread
 *******************************************/
void	ReductionMgr::require(int seq, ReductionTag tag, BigReal &data)
{
  #if PANIC > 0
  if (panicMode > 1)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die(""Panic due to wrong mode");
  }
  panicMode = 1;
  #endif

  // 1. find the correct sequence
  ReductionMgrData *current = find(seq);
  if (current->numData[tag] > 0)
  {
    // suspend thread until numData is 0...
    // ...then return value
  }
  data = current->tagData[tag];
  current->numRequire[tag]--;
  if (current->numRequire[tag] == 0)
	remove(seq);	// free it.
} /* ReductionMgr::require() */

/*******************************************
 *******************************************/

/*******************************************
 ReductionMgr::subscribe(): Allow a process to
 require the data.
 *******************************************/
void	ReductionMgr::subscribe(ReductionTag tag)
{
  #if PANIC > 0
  if (panicMode != 0)
  {
    iout << iERRORF << "Panic due to wrong mode: " << panicMode << "\n" << endi;
    NAMD_die(""Panic due to wrong mode");
  }
  panicMode = 0;
  #endif

  numSubscribed[tag]++;
} /* ReductionMgr::subscribe() */

/*******************************************
 ReductionMgr::unsubscribe(): Allow a process
 to no-longer require the data.
 *******************************************/
void	ReductionMgr::unsubscribe(ReductionTag tag)
{
  #if PANIC > 0
  panicMode = 2;
  #endif

  numSubscribed[tag]--;
  #if PANIC > 0
  if (numSubscribed < 0)
  {
    iout << iERRORF << "Too many unsubscribed. " << numSubscribed[tag]
	 << " for tag==" << tag << "\n" << endi;
    NAMD_die("Too many unsubscribed.");
  }
  #endif
} /* ReductionMgr::unsubscribe() */

