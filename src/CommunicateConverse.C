/***************************************************************************/
/*         (C) Copyright 1995,1996,1997 The Board of Trustees of the       */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *
 * CommunicateConverse - Converse version of the Communicate object.  Allows
 * user to establish id's for available nodes, establish connections, and
 * send/receive data.
 * 
 ***************************************************************************/
static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/CommunicateConverse.C,v 1.1004 1997/03/06 22:05:56 ari Exp $";

#include <iostream.h>
#include <string.h>
#include <unistd.h>
#include "CommunicateConverse.h"
#include "common.h"
#include "MessageQueue.h"
#include "MessageManager.h"
#include "pvmc.h"
// #define DEBUGM
#include "Debug.h"

#define MemCopy memcpy

#define Align1(x) (x)
#define Align4(x) (((x)%4)? ((((x)+4)/4)*4): (x))
#define Align8(x) (((x)%8)? ((((x)+8)/8)*8): (x))
#define Align(x)  Align1(x)

///////////////////////////  static data  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

CpvStaticDeclare(CmmTable, CsmMessages);

static void CsmHandler(void *msg)
{
  CmiGrabBuffer((char **)&msg);
  int *m = (int *) msg;
  // sending node acts as a single tag
  CmmPut(CpvAccess(CsmMessages), 1, &(m[2]), m);
}

static void mycpy(void *d, const void *s, int bytes)
{
        double *dst = (double *) d, *src = (double *) s;
        unsigned char *cdst, *csrc;

        while(bytes>8)
        {
                *dst++ = *src++;
                bytes -= 8;
        }
        cdst = (unsigned char *) dst;
        csrc = (unsigned char *) src;
        while(bytes)
        {
                *cdst++ = *csrc++;
                bytes--;
        }
}

///////////////////////////  private routines  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

// determines the size of the message
int CommunicateConverse::size_message(Message *msg) 
{
  Message::Types item_type;    //  Type of current item
  int item_size;
  int nitems = msg->items();    //  Number of items in Message
  int retsize=Align(sizeof(int));      // space for number of items

  msg->reset(); //current item = beginning of message
  for(int i=0; i < nitems; i++)
  {
    item_type = msg->type();
    retsize += Align(sizeof(short));
    item_size = msg->size();
    retsize += Align(sizeof(int));
    //inefficient right now.. maybe table would be a better idea
    switch(item_type) {
      case Message::CHAR : 
        retsize += Align(item_size*sizeof(char));
        break;
      case Message::SHORT : 
        retsize += Align(item_size*sizeof(short));
        break;
      case Message::USHORT :
        retsize += Align(item_size*sizeof(unsigned short));
        break;
      case Message::INT : 
        retsize += Align(item_size*sizeof(int));
        break;
      case Message::UINT : 
        retsize += Align(item_size*sizeof(unsigned int));
        break;
      case Message::LONG : 
        retsize += Align(item_size*sizeof(long));
        break;
      case Message::ULONG :
        retsize += Align(item_size*sizeof(unsigned long));
        break;
      case Message::FLOAT :
        retsize += Align(item_size*sizeof(float));
        break;
      case Message::DOUBLE :
        retsize += Align(item_size*sizeof(double));
        break;
      case Message::UNKNOWN :
        break;
    }
    msg->skip();  // move on to the next item in the message
  }
  return retsize;
}

// take the data from a Message object and pack it into the current send buffer.
// each message is packed in this order:
//  number of items  (int)
//    type of item 1  (short)
//    size of item 1, in number of elements  (int)
//    item 1 data  (various)
//    ...
//    type of item N  (short)
//    size of item N, in number of elements  (int)
//    item N data  (various)
char *CommunicateConverse::pack_message(Message *msg, char *mcmsg) 
{
  Message::Types item_type;    //  Type of current item
  int item_size;
  int nitems = msg->items();    //  Number of items in Message
  short stype;

  MemCopy(mcmsg, &nitems, sizeof(int));
  mcmsg += Align(sizeof(int));
  msg->reset();
  for(int i=0; i < nitems; i++)
  {
    item_type = msg->type();
    stype = (short) item_type;
    MemCopy(mcmsg, &stype, sizeof(short));
    mcmsg += Align(sizeof(short));
    item_size = msg->size();
    MemCopy(mcmsg, &item_size, sizeof(int));
    mcmsg += Align(sizeof(int));
    switch(item_type) {
      case Message::CHAR : 
        MemCopy(mcmsg, (char *)(msg->item()), item_size*sizeof(char));
        mcmsg += Align(item_size*sizeof(char));
        break;
      case Message::SHORT : 
        MemCopy(mcmsg, (short *)(msg->item()), item_size*sizeof(short));
        mcmsg += Align(item_size*sizeof(short));
        break;
      case Message::USHORT :
        MemCopy(mcmsg, (unsigned short *)(msg->item()), 
               item_size*sizeof(unsigned short));
        mcmsg += Align(item_size*sizeof(unsigned short));
        break;
      case Message::INT : 
        MemCopy(mcmsg, (int *)(msg->item()), item_size*sizeof(int));
        mcmsg += Align(item_size*sizeof(int));
        break;
      case Message::UINT : 
        MemCopy(mcmsg, (unsigned int *)(msg->item()), 
               item_size*sizeof(unsigned int));
        mcmsg += Align(item_size*sizeof(unsigned int));
        break;
      case Message::LONG : 
        MemCopy(mcmsg, (long *)(msg->item()), item_size*sizeof(long));
        mcmsg += Align(item_size*sizeof(long));
        break;
      case Message::ULONG :
        MemCopy(mcmsg, (unsigned long *)(msg->item()), 
               item_size*sizeof(unsigned long));
        mcmsg += Align(item_size*sizeof(unsigned long));
        break;
      case Message::FLOAT :
        MemCopy(mcmsg, (float *)(msg->item()), item_size*sizeof(float));
        mcmsg += Align(item_size*sizeof(float));
        break;
      case Message::DOUBLE :
        MemCopy(mcmsg, (double *)(msg->item()), item_size*sizeof(double));
        mcmsg += Align(item_size*sizeof(double));
        break;
      case Message::UNKNOWN :
        break;
    }
    msg->skip();  // move on to the next item in the message
  }
  return mcmsg;
}

// unpack the current receive message into individual Messages and MsgList
// structures, and look for a message with the given tag.  Returns a
// MsgList struct for the first matching tag; other MsgList structs for the
// other messages in the big message are added to the received msg queue
// messages are in this format:
//  number of individual messages    (int)
//  sending node (0 ... N-1)    (int)
//    data for message 1 (see comments for pack_message() above)
//    ...
//    data for message N
Message* CommunicateConverse::unpack_message(char *mcmsg, int &tag, int &node) 
{
  int messages, sender;
  Message *newmsg, *retmsg = NULL;
  int nitems, item_size, msg_tag;
  
  // get the number of messages, and sending node.
  MemCopy(&messages, mcmsg, sizeof(int));
  mcmsg += Align(sizeof(int));
  MemCopy(&sender, mcmsg, sizeof(int));
  mcmsg += Align(sizeof(int));

  for(int i=0; i < messages; i++) 
  {
    // get the tag of the message
    MemCopy(&msg_tag, mcmsg, sizeof(int));
    mcmsg += Align(sizeof(int));
    // get number of items in this message
    MemCopy(&nitems, mcmsg, sizeof(int));
    mcmsg += Align(sizeof(int));
    DebugM(1,"Received message with tag " << msg_tag
         << " from " << sender 
         << " containing " << nitems 
         << " Items\n");

    // create data structure for this message
    if(nitems > 0) 
    {
      newmsg = new Message;
      if ( newmsg == NULL ) {
        NAMD_die("new failed in CommunicateConverse::unpack_message");
      }

      // get all the items and add to the message
      for(int j=0; j < nitems; j++) 
      {
        short stype;
        MemCopy(&stype, mcmsg, sizeof(short));
        mcmsg += Align(sizeof(short));
        MemCopy(&item_size, mcmsg, sizeof(int));
        mcmsg += Align(sizeof(int));
  
  // for each item, allocate storage, and give to Message object as
  // storage that should NOT be copied, and that SHOULD be freed up
  // when the Message is deleted.
        switch(stype) {
          case Message::CHAR : {
            char *data = new char[item_size];
            if ( data == NULL ) {
              NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(char));
            mcmsg += Align(item_size*sizeof(char));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
          break;
          case Message::SHORT : {
            short *data = new short[item_size];
            if ( data == NULL ) {
              NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(short));
            mcmsg += Align(item_size*sizeof(short));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
            break;
          case Message::USHORT : {
            unsigned short *data = new unsigned short[item_size];
            if ( data == NULL ) {
               NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(unsigned short));
            mcmsg += Align(item_size*sizeof(unsigned short));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
            break;
          case Message::INT : {
            int *data = new int[item_size];
            if ( data == NULL ) {
               NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(int));
            mcmsg += Align(item_size*sizeof(int));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
            break;
          case Message::UINT : {
            unsigned int *data = new unsigned int[item_size];
            if ( data == NULL ) {
              NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(unsigned int));
            mcmsg += Align(item_size*sizeof(unsigned int));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
            break;
          case Message::LONG : {
            long *data = new long[item_size];
            if ( data == NULL ) {
              NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(long));
            mcmsg += Align(item_size*sizeof(long));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
            break;
          case Message::ULONG : {
            unsigned long *data = new unsigned long[item_size];
            if ( data == NULL ) {
              NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(unsigned long));
            mcmsg += Align(item_size*sizeof(unsigned long));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
            break;
          case Message::FLOAT : {
            float *data = new float[item_size];
            if ( data == NULL ) {
              NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(float));
            mcmsg += Align(item_size*sizeof(float));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
            break;
          case Message::DOUBLE : {
            double *data = new double[item_size];
            if ( data == NULL ) {
              NAMD_die("new failed in CommunicateConverse::unpack_message");
            }
            MemCopy(data, mcmsg, item_size*sizeof(double));
            mcmsg += Align(item_size*sizeof(double));
            (newmsg)->put(item_size, data, FALSE, TRUE);
          }
            break;
          case Message::UNKNOWN: break;
        }
      }

      newmsg->reset();
  
      // check and see if this is the message we're looking for
      if(!retmsg && ( (tag == msg_tag) || (tag < 0) ) ) 
      {
        retmsg = newmsg;
        tag = msg_tag;
        node = sender;
      } 
      else 
      {    
  // add the message to the received queue
        Communicate::remoteMsgs->add_message(newmsg, msg_tag, sender);
      }
    }
  }
  return retmsg;
}

        
///////////////////////////  protected routines  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

// send the entire list of MsgList messages as a single message.  Return
// number of messages sent.
// The destination node is the node of the first MsgList object;
// all MsgList objects following that one until the end of list, or the
// node number changes, are sent.
int CommunicateConverse::do_send_queue(int nodenum) 
{
  int sent = 0, nmsg = 0;
  Message *curmsg, *first;
  int tag, node, total_size = Align(CmiMsgHeaderSizeBytes + 2*sizeof(int));
  char *mcmsg, *curptr;

  if (Communicate::sendQueues[nodenum].num() == 0)
    return 0;

  if(!node_available(nodenum)) 
  {
    ErrorStatus = NOSEND;
    return 0;
  }
  first = curmsg = Communicate::sendQueues[nodenum].get_head(tag, node); 
  if(first == NULL)
    return 0;
  do{
    total_size += size_message(curmsg) + Align(sizeof(int));
    nmsg++;
    Communicate::sendQueues[nodenum].add_msg(curmsg, tag, node);
    curmsg = Communicate::sendQueues[nodenum].get_head(tag, node);
  } while(curmsg != first) ;
  mcmsg = (char *)CmiAlloc(total_size);
  curptr = mcmsg + Align(CmiMsgHeaderSizeBytes);
  MemCopy(curptr, &nmsg, sizeof(int));
  curptr += Align(sizeof(int));
  MemCopy(curptr, &thisNode, sizeof(int));
  curptr += Align(sizeof(int));
  do{
    MemCopy(curptr, &tag, sizeof(int));
    curptr += Align(sizeof(int));
    curptr = pack_message(curmsg,curptr);
    delete curmsg;
    curmsg = Communicate::sendQueues[nodenum].get_head(tag, node);
  } while(curmsg != NULL) ;
  CmiSetHandler(mcmsg, CsmHandlerIndex);
  CmiSyncSendAndFree(node, total_size, mcmsg);
  return nmsg;
}

int CommunicateConverse::do_send_msg(Message *msg, int node, int tag, int delmsg) 
{
  int one=1;

  if (msg == NULL)
    return 0;
  if(!node_available(node)) 
  {
    ErrorStatus = NOSEND;
    return 0;
  }
  int cursize = size_message(msg);
  cursize += Align(3*sizeof(int));
  char *mcmsg = (char *)CmiAlloc(CmiMsgHeaderSizeBytes+cursize);
  char *curptr = mcmsg+CmiMsgHeaderSizeBytes;
  MemCopy(curptr, &one, sizeof(int));
  curptr += Align(sizeof(int));
  MemCopy(curptr, &thisNode, sizeof(int));
  curptr += Align(sizeof(int));
  MemCopy(curptr, &tag, sizeof(int));
  curptr += Align(sizeof(int));
  pack_message(msg,curptr);
  CmiSetHandler(mcmsg, CsmHandlerIndex);
  CmiSyncSendAndFree(node, CmiMsgHeaderSizeBytes+cursize, mcmsg);
  if (delmsg) delete msg;
  return 1;
}
  
Message *CommunicateConverse::do_receive(int &node, int &tag) 
{
  Message *newmsg = NULL;
  int nomoremessages = 0;
  int rtag, itag;
  char *msg;
  
  while ( (newmsg == NULL) && !nomoremessages)
  {
    itag = ((node<0 || node>=TotalNodes) ? (CmmWildCard) : node);
    CmiDeliverMsgs(0);
    msg = (char *)CmmGet(CpvAccess(CsmMessages), 1, &itag, &rtag);
    if(msg) {
      newmsg = unpack_message(msg+CmiMsgHeaderSizeBytes, tag, node);
      CmiFree(msg);
    } else {
      nomoremessages = 1;
    }
  }
  return newmsg;
}


///////////////////////////  public routines  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

// class constructor.
CommunicateConverse::CommunicateConverse(int argc, char *argv[], int dbg) 
{
  // general initialization
  debug(dbg);

  // initialize converse
  // ConverseInit(argv);
  CsmHandlerIndex = CmiRegisterHandler((CmiHandler) CsmHandler);
  CpvInitialize(CmmTable, CsmMessages);
  CpvAccess(CsmMessages) = CmmNew();

  // Initialize PVM-converse library
  pvmc_init();

  // initialize superclass
  TotalNodes = CmiNumPes();
  thisNode = CmiMyPe();

  // Allocate the sendQueues for messages held in the WAIT state
  Communicate::sendQueues = new MessageQueue[TotalNodes];
  if ( Communicate::sendQueues == NULL )
  {
    NAMD_die("new failed in CommunicateConverse::CommunicateConverse");
  }
}


// class destructor
CommunicateConverse::~CommunicateConverse(void) 
{
  ConverseExit();
}


// add a node to the list of communications.  Return node number if OK,
//  (-1) if problem.
int CommunicateConverse::add_node(void *id) 
{
  if(id);
  if(debug())
    cout << "CommunicateConverse: error, cannot add new node.  Currently "
         << TotalNodes << " nodes." << endl;
  return (-1);
}

int *CommunicateConverse::get_tids()
{
  return ((int *)NULL);
}

/***************************************************************************
 * RCS INFORMATION:
 *
 *  $RCSfile: CommunicateConverse.C,v $
 *  $Author: ari $  $Locker:  $    $State: Exp $
 *  $Revision: 1.1004 $  $Date: 1997/03/06 22:05:56 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: CommunicateConverse.C,v $
 * Revision 1.1004  1997/03/06 22:05:56  ari
 * Removed Compute.ci
 * Comments added - more code cleaning
 *
 *
 ***************************************************************************/
