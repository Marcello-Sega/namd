//-*-c++-*-
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
 *      $RCSfile: Message.h,v $
 *      $Author: ari $        $Locker:  $             $State: Exp $
 *      $Revision: 1.778 $        $Date: 1997/01/28 00:30:48 $
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
 * $Log: Message.h,v $
 * Revision 1.778  1997/01/28 00:30:48  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:22  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:22  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.2  1996/12/12 20:14:50  milind
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.10  1996/04/27 20:42:24  billh
 * Added insertion  operator for Vector, and changed how Vector.h included.
 *
 * Revision 1.9  1995/03/08 14:34:24  nelson
 * Added copyright
 *
 * Revision 1.8  95/01/19  15:16:42  15:16:42  brunner (Robert Brunner)
 * Changed put(float d) to store d in a temporary variable, to work
 * around faulty type coercion produced by objectcenter compiler.  This
 * is not necessary for the standard compiler.
 * 
 * Revision 1.7  94/11/22  11:59:57  11:59:57  nelson (Mark T. Nelson)
 * Added function get_reallist_by_ref
 * 
 * Revision 1.6  94/09/03  18:10:52  18:10:52  billh (Bill Humphrey)
 * Inlined many functions, added 'don't copy, but delete storage' option
 * to 'put' routines to help performance.
 * 
 * Revision 1.5  94/09/03  09:57:47  09:57:47  nelson (Mark T. Nelson)
 * Added functins get_vectorlist_by_ref and 
 * get_intlist_by_ref
 * 
 * Revision 1.4  94/08/04  13:15:15  13:15:15  nelson (Mark T. Nelson)
 * Added Vector class
 * 
 * Revision 1.3  94/08/03  21:54:59  21:54:59  nelson (Mark T. Nelson)
 * Added unsigned shorts, ints, and longs
 * 
 * Revision 1.2  94/07/22  14:18:45  14:18:45  billh (Bill Humphrey)
 * Added optional third argument to put routines which put an array
 * of data; if third argument is FALSE, the array pointer is stored,
 * not copied, and so no storage is created/deleted for that item.
 * 
 * Revision 1.1  94/07/03  00:18:51  00:18:51  billh (Bill Humphrey)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef MESSAGE_H
#define MESSAGE_H

#include <string.h>
#include <iostream.h>
#include "common.h"

// forward references
class Message;
class Vector;


class Message {

public:
  // enum with types of data objects
  enum Types { CHAR, SHORT, USHORT, INT, UINT, LONG, ULONG, FLOAT, DOUBLE, UNKNOWN };

  // struct for MsgItem linked list
  typedef struct MsgItem_linklist {
    Types type;			// what type of data is this
    int size;			// if 0, scalar item; if > 0, array of len size
    int bytesize;		// size of item storage, in bytes
    void *item;			// pointer to the item; must be new/delete
    int needDelete;		// do we need to delete this item storage?
    struct MsgItem_linklist *next;	// next item in linked list
    struct MsgItem_linklist *prev;	// prev item in linked list
  } MsgItem;

private:
  // first and last elements of the linked list, and current item
  // (items are processed from head to tail)
  MsgItem *head, *tail, *curr;
  
  // items in list
  int Items;

  // return the Nth MsgItem
  MsgItem *msg_item(int N) {
    MsgItem *retval = NULL;
    if(N == (-1))
      retval = curr;
    else if(N == (Items - 1))
      retval = tail;
    else if(N < 0 || N >= Items)
      ;
    else {
      MsgItem *m = head;
      for(int i = 0; i < N; i++)
        m = m->next;
      retval = m;
    }
    return retval;
  }

  // delete the given MsgItem
  void del_item(MsgItem *);

  // delete the Nth MsgItem
  void del_item(int N) {
    del_item(msg_item(N));
  }
  
  // general put routine; called by other cases
  // arguments are the item, it's type, its size, size of 1 element,
  // whether the item should be copied (T) or just the pointer saved (F),
  // and whether the storage should be freed when the Message is deleted
  // (T or F, ignored if the message is being copied).
  Message& putmsg(void *, Types, int, int, int = TRUE, int = TRUE);
  
  // general get routine; called by other cases
  // arguments are the 
  Message& getmsg(void *data, Types t) {
    // check to see if there is an item, and it is of the correct type
    if(!curr || curr->type != t || !data) { 
	if (!curr)
	{
		iout << iERROR << "Message: no more items in Message\n" << endi; 
	}
	else if (curr->type != t)
	{
		iout << iERROR << "Message: data type mismatch in get\n" << endi;
	}
	else
	{
		iout << iERROR << "Message: !data\n" << endi;
	}
		
	return *this;
    }

    // copy the data into the given location
    memcpy(data, curr->item, curr->bytesize);

    // change the current item to the next item
    skip();
  
    return *this;
  }

public:
  Message(void) {
    head = tail = curr = NULL;
    Items = 0;
  }
  ~Message(void) {
    // delete all items in this message
    clear();
  }

  // print out the contents of the message
  friend ostream& operator<<(ostream &o, Message &msg);

  // clear the message; delete's all it's items
  Message& clear(void);

  // return number of items in this message
  int items(void) { return Items; }
  
  // return size of item N, in number of elements.  If N is outside range
  // of 0 ... Items() - 1, or there is no current item, return 0.
  // If N == (-1), return info about current item.
  int size(int N = (-1)) {
    MsgItem *m = msg_item(N);
    int retval = 0;
  
    // return 1 if a scalar (m->size == 0), or the array length if an array
    if(m)
      retval = (m->size > 0 ? m->size : 1);
    else
      iout << iERROR << "Message: cannot get size of item " << N << "\n" << endi;
    return retval;
  }

  // return type of item N.  If N is outside range
  // of 0 ... Items() - 1, or there is no current item, return 0.
  // If N == (-1), return info about current item.
  Types type(int N = (-1)) {
    MsgItem *m = msg_item(N);
    Types retval = UNKNOWN;
    if(m)
      retval = m->type;
    else
      iout << iERROR << "Message: cannot get type of item " << N << "\n" << endi;
    return retval;
  }

  // return, via a pointer to void (ugly, but necessary), the message item
  // If N == (-1), return info about current item.
  void *item(int N = (-1)) {
    MsgItem *m = msg_item(N);
    void *retptr = NULL;
    if(m)
      retptr = m->item;
    else
      iout << iERROR << "Message: cannot get pointer to item " << N << "\n" << endi;
    return retptr;
  }

  // delete the Nth item ... this changes the number of items in
  // the message
  Message& del(int N = (-1)) {
    del_item(msg_item(N));
    return *this;
  }

  //
  // routine to affect the status of the current item
  //

  // reset the current item to the beginning
  Message& reset(void) {
    curr = head;
    return *this;
  }
  
  // move the current item on to the next one, if available
  Message& skip(void) {
    if(curr)
      curr = curr->next;
    return *this;
  }
  
  // move the current item back to the previous one, if available
  Message& back(void) {
    if(curr)
      curr = curr->prev;
    return *this;
  }
  
  // set the current item to be the Nth item
  Message& current(int N) {
    curr = msg_item(N);
    return *this;
  }

  //
  // routines to get and put data out of/into the message.
  // NOTES for 'put' routines:
  // For scalar values, there is one argument: the scalar item.
  // For vector values, there are four arguments:
  //	n = size of vector, in number of elements
  //	d = pointer to data
  //	copy = flag whether the data should be copied into new location (T)
  //		or else just the given pointer stored (F)
  //	delstor = If copy=FALSE, this flag tells whether to delete the
  //		storage used by the data (T), or do not touch it (F).
  // Note: the copy and delstor are used to increase performance.  They
  // should be used in the following circumstances:
  //	a. Data to be sent/rec is in a location that will not change before the
  // msg is used, but should not be affected after the Message is deleted.
  // For this case, use copy=F, delstor=F.
  //	b. Data for a message has already had space allocated for it, and so
  // does not need to be copied.  Also, since the space has been allocated
  // already, it must be freed when the msg is sent.  For this case, use
  // copy=F, delstor=T.
  //	c. For data that is in a volatile location, it must be copied to new
  // storage, so use copy=T, and don't specify any value for delstor (it will
  // be ignored if copy=T).

  // specific put routines, for different data types
  Message& put(char *d, int copy=TRUE, int delstor=FALSE) {
    // null-terminated string
    return putmsg((void *)d, CHAR, strlen(d) + 1, sizeof(char), copy, delstor);
  }

  Message& put(char d) {			// character
    return putmsg((void *)(&d), CHAR, 0, sizeof(char));
  }

  Message& put(int n, char *d, int copy=TRUE, int delstor=FALSE) {
    // char array (no terminator)
    return putmsg((void *)d, CHAR, n, sizeof(char), copy, delstor);
  }

  Message& put(short d) {			// short
    return putmsg((void *)(&d), SHORT, 0, sizeof(short));
  }

  Message& put(int n, short *d, int copy=TRUE, int delstor=FALSE) {
    // short array
    return putmsg((void *)d, SHORT, n, sizeof(short), copy, delstor);
  }

  Message& put(unsigned short d) {		//  unsigned short
    return putmsg((void *)(&d), USHORT, 0, sizeof(unsigned short));
  }

  Message& put(int n, unsigned short *d, int copy=TRUE, int delstor=FALSE) {
    // unsigined short array
    return putmsg((void *)d, USHORT, n, sizeof(unsigned short), copy, delstor);
  }

  Message& put(int d) {				// int
    return putmsg((void *)(&d), INT, 0, sizeof(int));
  }

  Message& put(int n, int *d, int copy=TRUE, int delstor=FALSE) {
    // int array
    return putmsg((void *)d, INT, n, sizeof(int), copy, delstor);
  }

  Message& put(unsigned int d) {		//  unsigned int
    return putmsg((void *)(&d), UINT, 0, sizeof(unsigned int));
  }

  Message& put(int n, unsigned int *d, int copy=TRUE, int delstor=FALSE) {
    // unsigned int array
    return putmsg((void *)d, UINT, n, sizeof(unsigned int), copy, delstor);
  }

  Message& put(long d) {			// long
    return putmsg((void *)(&d), LONG, 0, sizeof(long));
  }

  Message& put(int n, long *d, int copy=TRUE, int delstor=FALSE) {
    // long array
    return putmsg((void *)d, LONG, n, sizeof(long), copy, delstor);
  }

  Message& put(unsigned long d) {		//  unsigned long
    return putmsg((void *)(&d), ULONG, 0, sizeof(unsigned long));
  }

  Message& put(int n, unsigned long *d, int copy=TRUE, int delstor=FALSE) {
    // unsigined long array
    return putmsg((void *)d, ULONG, n, sizeof(unsigned long), copy, delstor);
  }

  Message& put(float d) {			// float
    float d2=d;
    return putmsg((void *)(&d2), FLOAT, 0, sizeof(float));
  }

  Message& put(int n, float *d, int copy=TRUE, int delstor=FALSE) {
    // float array
    return putmsg((void *)d, FLOAT, n, sizeof(float), copy, delstor);
  }

  Message& put(double d) {			// double
    return putmsg((void *)(&d), DOUBLE, 0, sizeof(double));
  }

  Message& put(int n, double *d, int copy=TRUE, int delstor=FALSE) {
    // double array
    return putmsg((void *)d, DOUBLE, n, sizeof(double), copy, delstor);
  }

  //  NOTE:  For the Vector class, we just kind of fake it out by
  //         telling it that each vector is really just an array
  //         of three BigReals.  This will work great unless someone
  //	     adds virtual functions to the Vector class.  So don't
  //	     do that!!
  Message& put(Vector *v, int copy=TRUE, int delstor=FALSE) {
    // single vector (which is still an array, but only of size 3)
    return put(3, (BigReal *) v, copy, delstor);
  }

  Message& put(int n, Vector *v, int copy=TRUE, int delstor=FALSE) {
    // vector array
    return put(3*n, (BigReal *) v, copy, delstor);
  }

  // specific get routines, for different data types
  Message& get(char &d) { return getmsg((void *)(&d), CHAR); }
  Message& get(char *d) { return getmsg((void *)d, CHAR); }
  Message& get(short &d) { return getmsg((void *)(&d), SHORT); }
  Message& get(short *d) { return getmsg((void *)d, SHORT); }
  Message& get(unsigned short &d) { return getmsg((void *)(&d), USHORT); }
  Message& get(unsigned short *d) { return getmsg((void *)d, USHORT); }
  Message& get(int &d) { return getmsg((void *)(&d), INT); }
  Message& get(int *d) { return getmsg((void *)d, INT); }
  Message& get(unsigned int &d) { return getmsg((void *)(&d), UINT); }
  Message& get(unsigned int *d) { return getmsg((void *)d, UINT); }
  Message& get(long &d) { return getmsg((void *)(&d), LONG); }
  Message& get(long *d) { return getmsg((void *)d, LONG); }
  Message& get(unsigned long &d) { return getmsg((void *)(&d), ULONG); }
  Message& get(unsigned long *d) { return getmsg((void *)d, ULONG); }
  Message& get(float &d) { return getmsg((void *)(&d), FLOAT); }
  Message& get(float *d) { return getmsg((void *)d, FLOAT); }
  Message& get(double &d) { return getmsg((void *)(&d), DOUBLE); }
  Message& get(double *d) { return getmsg((void *)d, DOUBLE); }

  //  Again, for the vector class, fake it out by telling it that it is
  //  an array of BigReals
  Message& get(Vector *d) { return get((BigReal *) d); }

  //  These routines are used to get data by reference rather than
  //  copying them.  Only the most likely data types are here . . .
  Vector *get_vectorlist_by_ref();
  int *get_intlist_by_ref();
  Real *get_reallist_by_ref();

};

#endif
