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
 *	$RCSfile: ConfigList.h,v $
 *	$Author: ari $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:12 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Read in a configuration file of the form:
 *       keyword = information\n
 *-or-   keyword information\n
 *-or-   keyword = {\n line 0\n line 1\n ... line n\n}
 * and produces a database which can return a linked list of strings (char *)
 * to all the information fields associated with that keyword.
 * 
 *    A "word" is a seqeunce of characters that are not white space (see
 * isspace(3C) ).  The equals sign ('=') is optional (though if there is more
 * more than one equals sign, then the 2nd is not ignored).  The "information"
 * field may contain more than one word.  White space is ignored except that
 * white space between multiple words in the information field is maintained.
 * Everything on the line at and beyond a pound sign ('#') is ignored.  Hence
 * a data file can be:
 *   fullname = George   Washington # the first president of the U.S.
 *   fullname = Martha Washington   # his second wife
 * Once that is read in, all data associated with "name" can be retreived as
 *  StringList *strList = configFile.find("fullname");
 *  for (StringList *tmp=strList; tmp!=NULL; tmp = tmp -> next)
 *      cout << tmp->data << '\n';
 * Note:
 *   The returned StringList * is NOT new'ed.  Do NOT free it.
 *   Keywords are case INsensitive
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ConfigList.h,v $
 * Revision 1.1001  1997/03/19 11:54:12  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:58:20  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:30:31  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777.2.1  1997/01/27 22:45:11  ari
 * Basic Atom Migration Code added.
 * Added correct magic first line to .h files for xemacs to go to C++ mode.
 * Compiles and runs without migration turned on.
 *
 * Revision 1.777  1997/01/17 19:36:05  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.8  1995/11/11 17:25:27  dalke
 * added a
 *   keyword {
 * this is a line
 * this is another line
 * }
 *
 * option for reading multiple lines without specify the keyword
 * each time.
 *
 * Revision 1.7  95/10/07  01:49:57  01:49:57  hazen (Brett Hazen)
 * Memory Allocation error-checking added
 * 
 * Revision 1.6  1995/06/20  10:36:32  nelson
 * Fixed bug where delete [] was causing core dump on SGI's
 *
 * Revision 1.5  95/04/04  18:44:56  18:44:56  dalke (Andrew Dalke)
 * fixed const references
 * 
 * Revision 1.4  1995/03/08  14:46:33  nelson
 * Added copyright
 *
 * Revision 1.3  94/11/01  21:51:44  21:51:44  dalke (Andrew Dalke)
 * Fixed a delete problem
 * 
 * Revision 1.2  94/10/08  04:09:12  04:09:12  dalke (Andrew Dalke)
 * added missing .h file
 * 
 * Revision 1.1  94/06/27  20:46:09  20:46:09  dalke (Andrew Dalke)
 * Initial revision
 * 
 ***************************************************************************/

// This header introduces two names to the global name space
// They are:
//    StringList -- a linked list of char *.  It has an initilizer and
//      destructor to handle the char * by alloc'ing new space
//    ConfigList -- its constructor opens a file containing lines of the
//      format "keyword = data" and creates a listing of all the data
//      for each keyword.  This data can be retreived with "find(char *)"

#ifndef CONFIGLIST_H
#define CONFIGLIST_H
#include "common.h"
#include <string.h>

typedef struct StringList {
  char *data;
  StringList *next;
  StringList(char *newdata) {  // take a string, and copy it
     data = new char[strlen(newdata)+1];
     if ( data == NULL )
     {
       NAMD_die("new failed in struct StringList");
     }
     strcpy( data, newdata);
     next = NULL;
  }
  ~StringList( void) {  // just clear out my info
    delete [] data;
    data = NULL;
    next = NULL;
  }
} StringList;

class ConfigList {
  public:
    class ConfigListNode 
    {
    private:
    public:
      char *name;
      StringList *data;
      ConfigListNode *next;
      ConfigListNode( ConfigListNode *newnext, char *newname, 
                                          StringList *newdata) {
        name = new char[strlen(newname)+1];  // create space for the name
	if ( name == NULL )
	{
	  NAMD_die("new failed in ConfigListNode::ConfigListNode");
	}
        strcpy((char *) name, newname);      // and copy the new name
        data = newdata;
        next = newnext;
      }
      ~ConfigListNode( void) 
      {
        delete [] name;                  // delete the new'ed name
        name = NULL;
        next = NULL;
        StringList *curr=data, *next=NULL;        // go through the string list

        while ( curr!=NULL ) 
	{
          next = curr -> next;           // and delete each element
          delete curr;
	  curr = next;
        }
      }
    };
 private:
    ConfigListNode *theList;
       // copy the information into a String, as appropriate
       // this is really a "push"
    void add_element( char *s1, int len1, char *s2, int len2);
    ConfigListNode *find_key_word( const char *keyword) const;
    Bool isokay;
  public:
    ConfigList( const char *filename);
    Bool okay( void) { return isokay; }
    ~ConfigList( void);
    StringList *find( const char *name) const;   //search for values by name
       // NOTE: do not change or delete this information.  It is not new'ed
       //   and any changed you make will be permanent.
       
    ConfigListNode *head( void) const { return theList;  } // return everything
       // NOTE:  you _REALLY_ not not want to change the information
       //   you get from this (unless you really want to)
};


#endif // CONFIGLIST_H


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1001 $	$Date: 1997/03/19 11:54:12 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ConfigList.h,v $
 * Revision 1.1001  1997/03/19 11:54:12  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
