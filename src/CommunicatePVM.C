/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: CommunicatePVM.C,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.777 $	$Date: 1997/01/17 19:35:30 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * CommunicatePVM - PVM version of the Communicate object.  Allows
 * user to establish id's for available nodes, establish connections, and
 * send/receive data.
 *
 * When created, this object checks to see if other nodes have been
 * spawned yet; if no, this spawns them and confirms everthing is running
 * well.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: CommunicatePVM.C,v $
 * Revision 1.777  1997/01/17 19:35:30  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/12/06 19:52:20  ari
 * Initial revision
 *
 * Revision 1.26  1996/01/28 21:50:08  jean
 * Attempting to make stable RCS without Mark Nelson's
 * fma/pairlist decoupling
 *
 * Revision 1.27  1995/12/04 20:46:57  brunner
 * More communication stats - message sizes
 *
 * Revision 1.26  1995/11/22 12:10:55  brunner
 * Addded simple message accounting information.
 * for Communicate:print_comm_stats
 *
 * Revision 1.25  95/10/24  14:53:24  14:53:24  nelson (Mark T. Nelson)
 * Added exclusion of pvm_setopt call for SGI power challenge
 * 
 * Revision 1.24  95/10/07  01:32:35  01:32:35  hazen (Brett Hazen)
 * Memory Allocation error-checking added
 * 
 * Revision 1.23  1995/09/19  15:55:09  nelson
 * Added ifdef's for PVMe for IBM SP1 and SP2
 *
 * Revision 1.22  95/09/19  09:09:20  09:09:20  nelson (Mark T. Nelson)
 * Fixed problem where is SPECIFY_NUM_NODES was set and the number
 * of nodes wasn't specified on the command line namd core dumped
 * 
 * Revision 1.21  95/06/20  11:34:51  11:34:51  nelson (Mark T. Nelson)
 * Added explicit int declaration to get rid of warning messages
 * 
 * Revision 1.20  95/05/23  12:34:12  12:34:12  nelson (Mark T. Nelson)
 * Added code to support T3D with much help from Tom MacFarland from IPP Rechenzentrum
 * 
 * Revision 1.19  95/04/10  11:27:36  11:27:36  nelson (Mark T. Nelson)
 * Changed handling of previously static variables
 * 
 * Revision 1.18  95/03/20  11:53:52  11:53:52  nelson (Mark T. Nelson)
 * Removed numActive flag
 * 
 * Revision 1.17  95/03/18  02:42:11  02:42:11  nelson (Mark T. Nelson)
 * Reworked extensively to improve performance
 * 
 * Revision 1.16  95/03/09  10:43:06  10:43:06  nelson (Mark T. Nelson)
 * Added code for Exemplar port
 * 
 * Revision 1.15  95/03/08  14:42:42  14:42:42  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.14  95/02/27  14:52:09  14:52:09  nelson (Mark T. Nelson)
 * Changed so that PVM errors terminate namd if the reach a max number
 * 
 * Revision 1.13  95/01/30  14:43:12  14:43:12  nelson (Mark T. Nelson)
 * Changed do_receive in vain attempt to find memory leak on master node
 * 
 * Revision 1.12  94/12/14  16:13:42  16:13:42  nelson (Mark T. Nelson)
 * Added get_tids and changed do_receive to work with FMA
 * 
 * Revision 1.11  94/11/22  13:40:30  13:40:30  nelson (Mark T. Nelson)
 * Change way that do_receive works so that it retrieves messages until
 * it either finds one that matches, or runs out of messages
 * 
 * Revision 1.10  94/10/24  09:25:58  09:25:58  nelson (Mark T. Nelson)
 * Made minor changes to MsgList definitions
 * 
 * Revision 1.9  94/09/03  18:10:53  18:10:53  billh (Bill Humphrey)
 * Changed to put data in Message object without a second copy of the data
 * (thus, we allocate storage, and put the allocated storage into the Message
 * instead of copying the data and deleting the allocated space).
 * 
 * Revision 1.8  94/08/03  21:55:34  21:55:34  nelson (Mark T. Nelson)
 * Added unsigned short, int and long
 * 
 * Revision 1.7  94/07/26  16:52:26  16:52:26  billh (Bill Humphrey)
 * (Hopefully) fixed problem with broadcast ... now just loops through nodes
 * and sends copy of message to each node, instead of having specific
 * broadcast routine in child class.
 * 
 * Revision 1.6  94/07/22  14:03:57  14:03:57  billh (Bill Humphrey)
 * Changed pvm_advise to pvm_setopt due to new PVM 3.3 version.
 * 
 * Revision 1.5  94/07/03  01:22:54  01:22:54  billh (Bill Humphrey)
 * Made Communicate enum's and MsgItem struct public members of
 * Communicate, instead of global.
 * 
 * Revision 1.4  94/07/03  00:15:00  00:15:00  billh (Bill Humphrey)
 * New format ... user creates a Message object, and given that to a
 * Communicate object; these can be kept in a list to be all sent at once
 * or sent as soon as requested.
 * 
 ***************************************************************************/
static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/CommunicatePVM.C,v 1.777 1997/01/17 19:35:30 ari Exp $";

#include <iostream.h>
#include <string.h>
#include <unistd.h>
#include "CommunicatePVM.h"
#include "common.h"
#include "MessageQueue.h"
#include "MessageManager.h"

#ifdef SPECIFY_NUM_NODES
#include <ctype.h>
#endif

#if NAMD_DEBUG

#define COMM_DEBUG(m) { if(debug()) { cout << m; } }

#else

#define COMM_DEBUG(m)

#endif

// defines
#define COMM_HOSTS_TAG  0
#define COMM_DIE_TAG 	1
#define NAMDPVMTAG	1003    //  Tag used by PVM to transmit all messages

//  Maximum number of errors we will accept before giving up and dying
#define MAX_ERRS	100

typedef char *charptr;

///////////////////////////  static data  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

///////////////////////////  private routines  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

// take the data from a Message object and pack it into the current send buffer.
// each message is packed in this order:
//	tag		(int)
//	number of items	(int)
//		type of item 1	(short)
//		size of item 1, in number of elements	(int)
//		item 1 data	(various)
//		...
//		type of item N	(short)
//		size of item N, in number of elements	(int)
//		item N data	(various)
void CommunicatePVM::pack_message(Message *msg, int tag) 
{
  Message::Types item_type;		//  Type of current item
  int nitems = msg->items();		//  Number of items in Message
  int item_size;			//  Size of the current item
  short stype;				//  Short int version of type

  COMM_DEBUG("CommunicatePVM: Preparing message with tag " << tag \
	     << " and " << nitems << " items." << endl);

  pvm_pkint(&tag, 1, 1);	// tag
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
  pvm_pkint(&nitems, 1, 1);	// number of items
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
  
  // scan through the message items, and pack their data
  msg->reset();		// reset current item to beginnning
  for(int i=0; i < nitems; i++) 
  {
    item_type = msg->type();
    item_size = msg->size();
    stype = (short)item_type;
    pvm_pkshort(&stype, 1, 1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(short);
#endif /* CHEAP_PERFORMANCE */
    pvm_pkint(&item_size, 1, 1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
    
    switch(item_type) {
      case Message::CHAR : pvm_pkbyte((char *)(msg->item()), item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(char)*item_size;
#endif /* CHEAP_PERFORMANCE */
      			   break;
      case Message::SHORT : pvm_pkshort((short *)(msg->item()), item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(short)*item_size;
#endif /* CHEAP_PERFORMANCE */
      			    break;
      case Message::USHORT : pvm_pkushort((unsigned short *)(msg->item()), 
					   item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(unsigned short)*item_size;
#endif /* CHEAP_PERFORMANCE */
			    break;
      case Message::INT : pvm_pkint((int *)(msg->item()), item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(int)*item_size;
#endif /* CHEAP_PERFORMANCE */
      			  break;
      case Message::UINT : pvm_pkuint((unsigned int *)(msg->item()), 
					   item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(unsigned int)*item_size;
#endif /* CHEAP_PERFORMANCE */
			    break;
      case Message::LONG : pvm_pklong((long *)(msg->item()), item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(long)*item_size;
#endif /* CHEAP_PERFORMANCE */
      			   break;
      case Message::ULONG : pvm_pkulong((unsigned long *)(msg->item()), 
					   item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(unsigned long)*item_size;
#endif /* CHEAP_PERFORMANCE */
			    break;
      case Message::FLOAT : pvm_pkfloat((float *)(msg->item()), item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(float)*item_size;
#endif /* CHEAP_PERFORMANCE */
      			    break;
      case Message::DOUBLE : pvm_pkdouble((double *)(msg->item()),item_size,1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_sent+=sizeof(double)*item_size;
#endif /* CHEAP_PERFORMANCE */
      			     break;
      case Message::UNKNOWN: break;
    }
    msg->skip();	// move on to the next item in the message
  }
}


// unpack the current receive message into individual Messages and MsgList
// structures, and look for a message with the given tag.  Returns a
// MsgList struct for the first matching tag; other MsgList structs for the
// other messages in the big message are added to the received msg queue
// messages are in this format:
//	number of individual messages		(int)
//	sending node (0 ... N-1)		(int)
//		data for message 1 (see comments for pack_message() above)
//		...
//		data for message N
Message* CommunicatePVM::unpack_message(int &tag, int &node) 
{
  int messages, sender;
  Message *newmsg, *retmsg = NULL;
  int nitems, item_size, msg_tag;
  
  // get the number of messages, and sending node.
  pvm_upkint(&messages, 1, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
  pvm_upkint(&sender, 1, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */


  COMM_DEBUG("CommunicatePVM: Unpacking big message containing " \
         << messages << " messages" \
         << " from node " << sender << "; want tag " << tag << endl);

  for(int i=0; i < messages; i++) 
  {
    pvm_upkint(&msg_tag, 1, 1);		// get the tag of the message
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
    pvm_upkint(&nitems, 1, 1);		// get number of items in this message
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */

    // create data structure for this message
    if(nitems > 0) 
    {
      newmsg = new Message;
      if ( newmsg == NULL )
      {
	NAMD_die("new failed in CommunicatePVM::unpack_message");
      }

      // get all the items and add to the message
      short stype;
      for(int j=0; j < nitems; j++) 
      {
        pvm_upkshort(&stype, 1, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(short);
#endif /* CHEAP_PERFORMANCE */
        pvm_upkint(&item_size, 1, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
	
	// for each item, allocate storage, and give to Message object as
	// storage that should NOT be copied, and that SHOULD be freed up
	// when the Message is deleted.
        switch(stype) {
          case Message::CHAR : {
	       char *data = new char[item_size];
	       if ( data == NULL ) {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upkbyte(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(char)*item_size;
#endif /* CHEAP_PERFORMANCE */
	       (newmsg)->put(item_size, data, FALSE, TRUE);
	    }
      	    break;
          case Message::SHORT : {
	       short *data = new short[item_size];
	       if ( data == NULL )
	       {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upkshort(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(short)*item_size;
#endif /* CHEAP_PERFORMANCE */
	       (newmsg)->put(item_size, data, FALSE, TRUE);
	    }
      	    break;
          case Message::USHORT : {
	       unsigned short *data = new unsigned short[item_size];
	       if ( data == NULL )
	       {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upkushort(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(unsigned short)*item_size;
#endif /* CHEAP_PERFORMANCE */
	       (newmsg)->put(item_size, data, FALSE, TRUE);
	    }
      	    break;
          case Message::INT : {
	       int *data = new int[item_size];
	       if ( data == NULL )
	       {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upkint(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int)*item_size;
#endif /* CHEAP_PERFORMANCE */
	       (newmsg)->put(item_size, data, FALSE, TRUE);
	    }
      	    break;
          case Message::UINT : {
	       unsigned int *data = new unsigned int[item_size];
	       if ( data == NULL )
	       {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upkuint(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(unsigned int)*item_size;
#endif /* CHEAP_PERFORMANCE */
	       (newmsg)->put(item_size, data, FALSE, TRUE);
	    }
      	    break;
          case Message::LONG : {
	       long *data = new long[item_size];
	       if ( data == NULL )
	       {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upklong(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(long)*item_size;
#endif /* CHEAP_PERFORMANCE */
	       (newmsg)->put(item_size, data, FALSE, TRUE);
	    }
      	    break;
          case Message::ULONG : {
	       unsigned long *data = new unsigned long[item_size];
	       if ( data == NULL )
	       {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upkulong(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(unsigned long)*item_size;
#endif /* CHEAP_PERFORMANCE */
	       (newmsg)->put(item_size, data, FALSE, TRUE);
	    }
      	    break;
          case Message::FLOAT : {
	       float *data = new float[item_size];
	       if ( data == NULL )
	       {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upkfloat(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(float)*item_size;
#endif /* CHEAP_PERFORMANCE */
	       (newmsg)->put(item_size, data, FALSE, TRUE);
	    }
      	    break;
          case Message::DOUBLE : {
	       double *data = new double[item_size];
	       if ( data == NULL )
	       {
		 NAMD_die("new failed in CommunicatePVM::unpack_message");
	       }
	       pvm_upkdouble(data, item_size, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(double)*item_size;
#endif /* CHEAP_PERFORMANCE */
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
int CommunicatePVM::do_send_queue(int nodenum) 
{
  int errstat, newbufid, sent = 0;
  Message *curmsg;
  int tag, node;

  if (Communicate::sendQueues[nodenum].num() == 0)
	return 0;

  sent = Communicate::sendQueues[nodenum].num();

  // create new message, and sent it off
  if(node_available(nodenum)) 
  {
#ifdef T3D
    if((newbufid = pvm_initsend(PvmDataRaw)) >= 0) 
#else
    if((newbufid = pvm_initsend(PvmDataDefault)) >= 0) 
#endif
    {
      // take data out of the Message objects, and pack into the current buf
      // the first item must be the number of messages, the second the
      // sending node
      pvm_pkint(&sent, 1, 1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
      pvm_pkint(&thisNode, 1, 1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */

      // pack all the individual messages
      while(curmsg = Communicate::sendQueues[nodenum].get_head(tag, node)) 
      {
          pack_message(curmsg, tag);
	  delete curmsg;
      }
      
      // finally, send the complete message
      errstat = pvm_send(tids[nodenum], NAMDPVMTAG);
#ifdef CHEAP_PERFORMANCE
      num_remote_sends++;
#endif
    }

    ErrorStatus = (errstat == 0 ? NOERROR : NOSEND);

    if(ErrorStatus == NOERROR) 
    {      
      COMM_DEBUG("CommunicatePVM: Sent message in buffer " << newbufid \
	     << " containing " << sent << " messages" \
             << " from " << myHost << " to " << node \
	     << " with tag " << tag << endl);
    } 
    else 
    {
      COMM_DEBUG("Communicate PVM: Error " << errstat
           << " sending message in buffer " << newbufid
 	   << " containing " << sent << " messages"
           << " from " << myHost << " to " << node
	   << " with tag " << tag << endl);
      sent = 0;
    }
  } 
  else 
  {
    ErrorStatus = NOSEND;
    sent = 0;
  }

  return sent;
}

  
Message *CommunicatePVM::do_receive(int &node, int &tag) 
{
  int bufid, size, checknode, checktag;
  Message *newmsg = NULL;
  static int numErrors = 0;
  Bool nomoremessages=FALSE;

  while ( (newmsg == NULL) && !(nomoremessages) )
  {
  	checknode = (node < 0 || node >= TotalNodes ? (-1) : tids[node]);
  	checktag = NAMDPVMTAG;

  	// clear out any current message
  	if((bufid = pvm_getrbuf()) > 0)
    	   pvm_freebuf(bufid);

  	// get the message, if available
  	bufid = pvm_probe(checknode, checktag);

  	// check if there really was a message
  	if(bufid > 0) 
  	{		
    		// yes, there was ...
    		// get info about message
    		if(pvm_bufinfo(bufid, &size, &checktag, &checknode) < 0) 
    		{
    			cout << "CommunicatePVM: error, receive shows msg, but "
 		             << "cannot get info about this buffer." << endl;
      			cout << "   id = " << bufid << endl;
    		} 
    		else if(size <= 0) 
    		{
      			cout << "CommunicatePVM: error, received message has size "
           		     << size << endl;
    		} 
    		else 
    		{
			// we have a valid message.  break it up into Messages,
			// and look for a message with a matching tag.
      			bufid = pvm_recv(checknode, checktag);
#ifdef CHEAP_PERFORMANCE
			num_remote_recvs++;
#endif
      			newmsg = unpack_message(tag, node);
      			pvm_freebuf(bufid);
    		}
 	 } 
  	 else if(bufid < 0) 
    	 {
    		cout << "CommunicatePVM: error, receive cannot receive msg "
	             << "from node " << node << ", tag " << tag << endl;

    		numErrors++;

    		if (numErrors > MAX_ERRS)
    		{
			NAMD_die("Maximum number of PVM receive errors exceeded.  PVM is hosed!!");
    		}
  	}
	else
	{
		nomoremessages=TRUE;
	}
   }

  // return message object
  return newmsg;
}


///////////////////////////  public routines  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// class constructor.  argc and argv are the arguments for use in spawning
// new tasks
CommunicatePVM::CommunicatePVM(int argc, char *argv[], int dbg) 
{
  int i, mytid, parent_tid, testnumhosts, added;
  int *child_ready, reported, rep_host, rep_tid; 
  char *exnmptr;
  charptr *spawnargv = NULL;
#ifdef CRAY
  int my_pe;		//  Processor element # for T3D
#endif

//  Declare the PVM host and task pointers.  Take into account
//  machines with older versions of PVM
#if defined(EXEMPLAR) || defined(PVMe)
  //  This is the old hostinfo pointer for the Exemplar and PVMe
  struct hostinfo *hostp;
#ifdef PVMe
  //  The old taskinfo pointer for PVMe
  struct taskinfo *taskp;
#endif  // PVMe
#else
  //  Newer pointers for everyone else
  struct pvmhostinfo *hostp;
  struct pvmtaskinfo *taskp;
#endif  // EXEMPLAR || PVMe

#ifdef SPECIFY_NUM_NODES
  int len;
#else
  int pvmd_tid;
#endif

  //  Initialize useful info
  execName = NULL;
  maxHosts = (-1);
  myHost = 0;
  tids = NULL;

  // general initialization
  debug(dbg);

  // this is the first PVM communication object

  // make copy of executable name
  exnmptr = strrchr(argv[0],'/');
  if(!exnmptr)
    exnmptr = argv[0];
  else
    exnmptr++;
  execName = new char[strlen(exnmptr) + 1];
  if ( execName == NULL )
  {
    NAMD_die("new failed in CommunicatePVM::unpack_message");
  }
  strcpy(execName,exnmptr);

  // initialize; assume as first just a single node available
  maxHosts = TotalNodes = 1;
  myHost = thisNode = 0;
  parent_tid = (-1);

  // determine my tid number, and if this process is the parent process
  mytid = pvm_mytid();
  if(mytid >= 0) 
  {
	// PVM is indeed running, get info about the virtual machine.
#ifdef T3D
	//  This section is specific to the T3D
	parent_tid = pvm_gettid(0,0);

	// get information about the virtual machine.
	// used to set number of hosts, and to spawn processes if this is the
	// master node.
	TotalNodes=pvm_gsize(0);
	my_pe=pvm_get_PE(pvm_mytid());
	if(TotalNodes==1) parent_tid=(-1);
	if(my_pe==0) parent_tid=(-1);
#else
	//  This is the generic PVM code for all other platforms
    	parent_tid = pvm_parent();
#ifndef SGI64
    	pvm_setopt(PvmRoute, PvmRouteDirect);
#endif

    	// get information about the virtual machine.
    	// used to set number of hosts, and to spawn processes if this is the
    	// master node.
    	if(pvm_config(&TotalNodes, NULL, &hostp) != 0) 
    	{
        	// error getting config info; go back to being single-node machine
		COMM_DEBUG("Error: cannot get virtual machine configuration.\n" \
      	       		<< "Going to single-node status." << endl);
		TotalNodes = 1;
		parent_tid = (-1);
    	}
#endif
  }

#ifdef SPECIFY_NUM_NODES
  //  This is a machine like the Exemplar, where there is one host, but you
  //  want to start multiple processes.  For this, we need the user to specify
  //  how many nodes and get this value from the command line
  if (argc != 3)
  {
	//  NOTE:  We *DON'T* want to use NAMD_die to exit here because that
	//  	   routine will try and delete the Communicate object, and
	//   	   since we are in the constructor, core dumps insue . . .
	cerr    << "ABNORMAL TERMINATION - " 
		<< "Usage: namd <config_file> <num_nodes>"
		<< endl;
	exit(-1);
  }

  len = strlen(argv[2]);

  for (i=0; i<len; i++)
  {
	if (!isdigit(argv[2][i]))
	{
		//  NOTE:  We *DON'T* want to use NAMD_die to exit here because that
		//  	   routine will try and delete the Communicate object, and
		//   	   since we are in the constructor, core dumps insue . . .
		cerr    << "ABNORMAL TERMINATION - " 
			<< "Number of nodes must be a number!!!"
			<< endl;
		exit(-1);
	}
  }

  TotalNodes = atoi(argv[2]);
#endif	// SPECIFY_NUM_NODES

  tids = new int[TotalNodes];
  if ( tids == NULL )
  {
    NAMD_die("new failed in CommunicatePVM::CommunicatePVM");
  }
#ifdef T3D
  if(my_pe == 0)
#else
  if(parent_tid < 0) 
#endif
  {	
	// this is the parent process; spawn others
      	tids[0] = mytid;
      	added = 1;
      	myHost = 0;
      	if(TotalNodes > 1) 
	{
#ifndef T3D
		if(argc > 1) 
		{
  			spawnargv = new charptr[argc];
			if ( spawnargv == NULL )
			{
			  NAMD_die("new failed in CommunicatePVM::CommunicatePVM");
  }
  			for(i=0; i < (argc - 1); i++)
			{
    				spawnargv[i] = argv[i+1];
			}

  			spawnargv[argc - 1] = NULL;
		}
#endif	//  Not T3D

#ifdef SPECIFY_NUM_NODES
		//  This is for machines like the Exemplar where all the
		//  processes are explicitly spawned on one machine

		//  Start all the jobs on one host
		for (i=1; i<TotalNodes; i++)
		{
    			pvm_spawn(execName, spawnargv, PvmTaskHost, 
			   (hostp[0]).hi_name, 1, tids + added);
    			if(tids[added] < 0)
			{
      				cout << "\nCommunicatePVM: error, cannot start job # "
				     << i;
			}
			else
			{
				added++;
			}
		}
#elif T3D
		//  For the T3D, the jobs are automatically spawned,
		//  so all we need to do is get the tids
	        for(i=1; i < TotalNodes; i++)
	        {
			*(tids+added)=pvm_gettid(0,i);
			added++;
		}
#else
		//  This is a "normal" pvm setup, so start a job on
		//  each host in the virtual machine

		//  Start jobs on each host in the machine

		//  Find master host so we don't start another job there
		pvm_tasks(0, NULL, &taskp);
		pvmd_tid = (taskp[0]).ti_host;

		// go to each host and spawn one process
#ifndef PVMe
		for(i=0; i < TotalNodes; i++) 
#else
		for(i=1; i < TotalNodes; i++)
#endif
		{
	  		if(pvmd_tid != (hostp[i]).hi_tid) 
			{
	      			COMM_DEBUG("Parent starting job on host '"  \
				    	     << (hostp[i]).hi_name \
	           		    	     << "' ...");

	    			pvm_spawn(execName, spawnargv, PvmTaskHost, 
				   (hostp[i]).hi_name, 1, tids + added);

	    			if(tids[added] < 0)
				{
	      				cout << "\nCommunicatePVM: error, cannot start job on "
	           			     << (hostp[i]).hi_name << endl;
				}
	    			else 
				{
	        			COMM_DEBUG(" new tid=" << tids[added]  \
						   << " (" << added+1 \
	             			     	   << " now added)" << endl);

	      				added++;
	    			}
	  		}
		}
#endif

		// adjust maxHosts to the actual number of hosts which responded
		maxHosts = TotalNodes = added;
      
#ifdef T3D
		pvm_initsend(PvmDataRaw);
#else
		pvm_initsend(PvmDataDefault);
#endif
		pvm_pkint(&TotalNodes,1,1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
		pvm_pkint(tids, TotalNodes, 1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int)*TotalNodes;
#endif /* CHEAP_PERFORMANCE */
		pvm_mcast(tids+1, TotalNodes-1, COMM_HOSTS_TAG);
#ifdef T3D
		pvm_initsend(PvmDataRaw);
#else
		pvm_initsend(PvmDataDefault);
#endif

		// wait for the spawned processes to report back that they're ready
		child_ready = new int[TotalNodes];
		if ( child_ready == NULL )
		{
		  NAMD_die("new failed in CommunicatePVM::CommunicatePVM");
		}
		for(i=0; i < TotalNodes; child_ready[i++] = FALSE);
		COMM_DEBUG("CommunicatePVM: Parent process waiting for children " \
		       << "to report back ..." << endl);
		reported = 1;		// since the parent is already ready
		while(reported < TotalNodes) 
		{
			if(pvm_nrecv((-1),COMM_HOSTS_TAG) > 0) 
			{
#ifdef CHEAP_PERFORMANCE
			  num_remote_recvs++;
#endif
		    		pvm_upkint(&rep_host,1,1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
		    		pvm_upkint(&rep_tid,1,1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
		    		if(rep_host >= 0 && rep_host < TotalNodes
		    				 && !(child_ready[rep_host])) 
				{
		      			child_ready[rep_host] = TRUE;
		      			reported++;
		     			COMM_DEBUG("Child " << rep_host << " (tid = " << rep_tid \
		  			        << " ) reported ready." << endl);
		    		} 
				else 
				{
		 		     cout << "Error with child " << rep_host
		     		          << " reporting to parent.\n"
		       		          << "   rep_host = " << rep_host << "\n"
		        		  << "   rep_tid  = " << rep_tid << "\n"
		      		          << "   child_ready[] = " 
					  << child_ready[rep_host] << endl;
		    		}
	  		}
		}

		delete [] child_ready;
      	}
  } 
  else 
  {			// this is a child process; get data from pops
      pvm_recv(parent_tid, COMM_HOSTS_TAG);
      pvm_upkint(&testnumhosts,1,1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
      if(testnumhosts != TotalNodes) 
      {
	COMM_DEBUG("For tid " << mytid << ": Adjusting # of nodes from " \
	       << TotalNodes << " to " << testnumhosts << endl);
	maxHosts = TotalNodes = testnumhosts;
      }
      pvm_upkint(tids, TotalNodes, 1);
#ifdef CHEAP_PERFORMANCE
	remote_bytes_rcvd+=sizeof(int)*TotalNodes;
#endif /* CHEAP_PERFORMANCE */
      for(i=1; i < TotalNodes; i++)
	if(mytid == tids[i]) { myHost = i; break; }
#ifdef T3D
      pvm_initsend(PvmDataRaw);
#else
      pvm_initsend(PvmDataDefault);
#endif
      pvm_pkint(&myHost,1,1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
      pvm_pkint(&mytid,1,1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
      pvm_send(parent_tid,COMM_HOSTS_TAG);
#ifdef CHEAP_PERFORMANCE
      num_remote_sends++;
#endif

#ifdef T3D
      pvm_initsend(PvmDataRaw);
#else
      pvm_initsend(PvmDataDefault);
#endif
  }
    
  // set the base class variable which indicates which host this is
  thisNode = myHost;

  // Allocate the sendQueues for messages held in the WAIT state
  Communicate::sendQueues = new MessageQueue[TotalNodes];
  if ( Communicate::sendQueues == NULL )
  {
    NAMD_die("new failed in CommunicatePVM::CommunicatePVM");
  }
}


// class destructor
CommunicatePVM::~CommunicatePVM(void) 
{
  // assume the user has made sure all the other nodes are ready to die ...
  // if this is the last active CommunicatePVM object, clean up
  if(TotalNodes > 1) 
  {
    if(myHost == 0) 
    {		
#ifndef T3D
      // this is the parent process; kill others
      COMM_DEBUG("CommunicatePVM: Parent process killing all spawned processes " \
             << "before exiting ..." << endl);

      for(int i=1; i < TotalNodes; i++)
        pvm_kill(tids[i]);

      pvm_exit();		// and then exit the parent from pvm
#endif

      delete [] tids;
      delete [] execName;

    } 
    else 
    {
      delete [] tids;
      delete [] execName;
#ifdef T3D
      exit(1);
#else
      while(TRUE)
        sleep(1);	 // just sit and wait for the end
#endif
    }
  }
}


// add a node to the list of communications.  Return node number if OK,
//  (-1) if problem.
int CommunicatePVM::add_node(void *id) 
{
  if(id);
  if(debug())
    cout << "CommunicatePVM: error, cannot add new node.  Currently "
         << TotalNodes << " nodes." << endl;
  return (-1);
}

int *CommunicatePVM::get_tids()
{
	return (tids);
}

int CommunicatePVM::do_send_msg(Message *msg, int node, int tag, int delmsg) 
{
  int errstat, newbufid, sent = 1;

  if (msg == NULL)
	return 0;

  // create new message, and sent it off
  if(node_available(node)) 
  {
#ifdef T3D
    if((newbufid = pvm_initsend(PvmDataRaw)) >= 0) 
#else
    if((newbufid = pvm_initsend(PvmDataDefault)) >= 0) 
#endif
    {
      // take data out of the Message objects, and pack into the current buf
      // the first item must be the number of messages, the second the
      // sending node
      pvm_pkint(&sent, 1, 1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */
      pvm_pkint(&thisNode, 1, 1);
#ifdef CHEAP_PERFORMANCE
      remote_bytes_sent+=sizeof(int);
#endif /* CHEAP_PERFORMANCE */

      // pack all the individual messages
      pack_message(msg, tag);
      
      // finally, send the complete message
      errstat = pvm_send(tids[node], NAMDPVMTAG);
#ifdef CHEAP_PERFORMANCE
      num_remote_sends++;
#endif

    }

    ErrorStatus = (errstat == 0 ? NOERROR : NOSEND);

    if(ErrorStatus == NOERROR) 
    {      
      if (delmsg)
	delete msg;

      COMM_DEBUG("CommunicatePVM: Sent message in buffer " << newbufid \
	     << " containing " << sent << " messages" \
             << " from " << myHost << " to " << node \
	     << " with tag " << tag << endl);
    } 
    else 
    {
      cout << "Communicate PVM: Error " << errstat
           << " sending message in buffer " << newbufid
 	   << " containing " << sent << " messages"
           << " from " << myHost << " to " << node
	   << " with tag " << tag << endl;
      sent = 0;
    }
  } 
  else 
  {
    ErrorStatus = NOSEND;
    sent = 0;
  }

  return sent;
}
