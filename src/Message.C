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
 *      $RCSfile: Message.C,v $
 *      $Author: ari $        $Locker:  $             $State: Exp $
 *      $Revision: 1.1000 $        $Date: 1997/02/06 15:58:38 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * Message - contains a list of items that comprise a set of data to be
 * sent, or received from, another node in a parallel architecture.  A person
 * creates a Message object, loads it with data to be sent, and gives it
 * to a Communicate object.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Message.C,v $
 * Revision 1.1000  1997/02/06 15:58:38  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:46  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:36:21  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.11  1995/10/09 04:44:05  hazen
 * Memory allocation error-checking added
 *
 * Revision 1.10  1995/05/13  10:44:06  nelson
 * Made change to allow BigReal to be mapped to float rather than double
 *
 * Revision 1.9  95/04/10  13:34:46  13:34:46  nelson (Mark T. Nelson)
 * Removed output operator because of gcc complaints
 * 
 * Revision 1.8  95/03/08  14:33:20  14:33:20  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.7  94/11/22  12:00:10  12:00:10  nelson (Mark T. Nelson)
 * Added function get_reallist_by_ref
 * 
 * Revision 1.6  94/09/03  18:10:19  18:10:19  billh (Bill Humphrey)
 * Inlined many functions, added 'don't copy, but delete storage' option
 * to 'put' routines to help performance.
 * 
 * Revision 1.5  94/09/03  09:58:15  09:58:15  nelson (Mark T. Nelson)
 * Added functions get_vectorlist_by_ref and
 * get_inliat_by_ref
 * 
 * Revision 1.4  94/09/02  16:40:18  16:40:18  nelson (Mark T. Nelson)
 * changed error messages in get
 * 
 * Revision 1.3  94/08/03  21:55:13  21:55:13  nelson (Mark T. Nelson)
 * Added unsigned short, int, and long
 * 
 * Revision 1.2  94/07/22  14:17:44  14:17:44  billh (Bill Humphrey)
 * Added optional third argument to put routines which put an array
 * of data; if third argument is FALSE, the array pointer is stored,
 * not copied, and so no storage is created/deleted for that item.
 * 
 * Revision 1.1  94/07/03  00:18:51  00:18:51  billh (Bill Humphrey)
 * Initial revision
 * 
 ***************************************************************************/
static char ident[] = "@(#)$Header: /home/cvs/namd/cvsroot/namd2/src/Attic/Message.C,v 1.1000 1997/02/06 15:58:38 ari Exp $";

#include <string.h>
#include "Message.h"


// delete the given MsgItem
void Message::del_item(MsgItem *m) {

  if(!m)
    return;

  // change the current message, if necessary
  if(m == curr) {
    if(m->next)
      curr = m->next;
    else
      curr = m->prev;
  }

  // see if this is the first or last message
  if(m == head)
    head = m->next;
  if(m == tail)
    tail = m->prev;

  // unlink this message    
  if(m->prev)
    (m->prev)->next = m->next;

  if(m->next)
    (m->next)->prev = m->prev;

  // delete the item stored in this MsgItem, and the MsgItem itself
  if(m->needDelete && m->item)
    delete [] ((char *)(m->item));  
  delete m;
  Items--;
}


// clear the message; delete's all it's items
Message& Message::clear(void) {
  while(head)
    del_item(head);

  return *this;
}


// general put routine; called by other cases
// arguments are the data item, it's type, its size, 
// size of each element (in bytes), and whether to copy the data,
// or just store a pointer to the data
Message& Message::putmsg(void *data, Types t, int s, int ts, int copy,
	int delstor) {
  MsgItem *m = new MsgItem;
  if ( m == NULL )
  {
    NAMD_die("Memory allocation failed in Message::putmsg");
  }

  // make sure we have some data of the proper size
  if(!data || s < 0 || ts < 1)
    return *this;

  // initialize new MsgItem
  m->type = t;
  m->size = s;

  // save the total storage of this item, in bytes
  if(s > 0)
    m->bytesize = s * ts;
  else
    m->bytesize = ts;

  // allocate storage for the data, and copy it
  m->needDelete = copy;
  if(copy) {
    m->item = (void *)(new char[m->bytesize]);
    if ( m->item == NULL )
    {
      NAMD_die("Memory allocation failed in Message::putmsg");
    }
    memcpy(m->item, data, m->bytesize);
    m->needDelete = copy;
  } else {
    m->item = (void *)data;
    m->needDelete = delstor;
  }

  // put MsgItem at end of linked list
  if(Items == 0) {
    curr = head = tail = m;
    m->prev = m->next = NULL;
  } else {
    m->prev = tail;
    m->next = NULL;
    tail->next = m;
    tail = m;
  }
  
  Items++;
  return *this;
}

// get the current message data by reference
Vector *Message::get_vectorlist_by_ref()

{
  Vector *tmp;

  // check to make sure there is an item and it is the correct type
#ifdef SHORTREALS
  if (!curr || curr->type != FLOAT)
#else
  if (!curr || curr->type != DOUBLE)
#endif
  {
	if (!curr)
	{
		cerr << 
		  "Message: get_vectorlist_by_ref past end of message" << endl;
	}
	else
	{
		cerr << 
		  "Message: get_vectorlist_by_ref data type mismatch" << endl;
	}

	return(NULL);
  } 

  tmp = (Vector *) curr->item;

  skip();

  return(tmp);
}

// get the current message data by reference
int *Message::get_intlist_by_ref()

{
  int *tmp;

  // check to make sure there is an item and it is the correct type
  if (!curr || curr->type != INT)
  {
	if (!curr)
	{
		cerr << 
		  "Message: get_intlist_by_ref past end of message" << endl;
	}
	else
	{
		cerr << 
		  "Message: get_intlist_by_ref data type mismatch" << endl;
	}

	return(NULL);
  } 

  tmp = (int *) curr->item;

  skip();

  return(tmp);
}

// get the current message data by reference
Real *Message::get_reallist_by_ref()

{
  float *tmp;

  // check to make sure there is an item and it is the correct type
  if (!curr || curr->type != FLOAT)
  {
	if (!curr)
	{
		cerr << 
		  "Message: get_reallist_by_ref past end of message" << endl;
	}
	else
	{
		cerr << 
		  "Message: get_reallist_by_ref data type mismatch" << endl;
	}

	return(NULL);
  } 

  tmp = (float *) curr->item;

  skip();

  return(tmp);
}

// print out the contents of the message
/*
ostream& operator<<(ostream &o, Message &msg) {
  o << "Message with " << msg.Items << " items:";
  int i = 0;
  Message::MsgItem *mi = msg.head;
  while(mi) {
    o << "\ntype=" << mi->type << ", size=" << mi->size
      << ", bytesize=" << mi->bytesize << "\ndata[0]=";
    switch(mi->type) {
      case Message::CHAR : o << ((char *)(mi->item))[0]; break;
      case Message::SHORT : o << ((short *)(mi->item))[0]; break;
      case Message::USHORT : o << ((unsigned short *)(mi->item))[0]; break;
      case Message::INT : o << ((int *)(mi->item))[0]; break;
      case Message::UINT : o << ((unsigned int *)(mi->item))[0]; break;
      case Message::LONG : o << ((long *)(mi->item))[0]; break;
      case Message::ULONG : o << ((unsigned long *)(mi->item))[0]; break;
      case Message::FLOAT : o << ((float *)(mi->item))[0]; break;
      case Message::DOUBLE : o << ((double *)(mi->item))[0]; break;
      case Message::UNKNOWN : o << "unknown"; break;
    }
    mi = mi->next;
  }
  
  return o;
}
*/
