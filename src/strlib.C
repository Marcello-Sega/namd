/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:
 *	strlib contains a bunch of functions that are useful for basic 
 * file input and string manipulation.  See strlib.h for a list and 
 * description of the functions that are available.
 *
 ***************************************************************************/

#include "strlib.h"

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_read_line				*/
/*									*/
/*   INPUTS:								*/
/*	fd - FILE to read line from					*/
/*	buf - buffer to put line into					*/
/*									*/
/*   OUTPUTS:								*/
/*	this function returns a 0 if the line is read successfully or   */
/*   a -1 if an EOF is encountered					*/
/*									*/
/*	NAMD_read_line reads in a line from an open file.  It places    */
/*   the line that is read in into the buffer passed in via buf.  If    */
/*   an EOF is encountered, a -1 is returned.  If an EOF is encountered */
/*   in the middle of a line, the program will terminate abnormally.    */
/*   Also, any comments of the form {blah blah blah} will automatically */
/*   be skipped by NAMD_read_line, even if usch comments span lines.    */
/*   Also, the string will be left justified.  That is, any spaces at   */
/*   the begining of the line will be removed.				*/
/*									*/
/************************************************************************/

int NAMD_read_line(FILE *fd, char *buf)

{
	int i=0;	//  Current position in buf
	int c;		//  Character read in from file

	/*  Loop and read characters until we get either an EOF or a    */
	/*  newline							*/
	while ( ((c=fgetc(fd)) != EOF) && (c != '\n') )
	{
		/*  If we encounter a bracketed comment, skip it.  This */
		/*  basically means read EVERYTHING until the next } and*/
		/*  throw it into the big bit bucket			*/
		if (c == '{')
		{
			while ( ((c=fgetc(fd)) != EOF) && (c!='}') )
			{
			}

			if (c==EOF)
			{
				char err_msg[512];

				sprintf(err_msg, "ABNORMAL EOF FOUND - buffer=*%s*\n", 
				   buf);
				NAMD_die(err_msg);
			}

			continue;
		}

		/*  Also, use this little piece of logic to remove      */
		/*  any leading spaces from the line			*/
		if ((i>0) || !isspace(c))
		{
			buf[i] = c;
	
			i++;
		}
	}

	/*  NULL terminate the string					*/
	buf[i]=STRINGNULL;

	/*  Check for an EOF in the middle of a line			*/
	if ((c==EOF) && (i!=0))
	{
		char err_msg[512];

		sprintf(err_msg, "ABNORMAL EOF FOUND - buffer=*%s*\n", 
		   buf);
		NAMD_die(err_msg);
	}

	/*  Return the appropriate value				*/
	if (c==EOF)
		return(-1);
	else
		return(0);
}
/*			END OF FUNCTION NAMD_read_line			*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_remove_comment			*/
/*									*/
/*   INPUTS:								*/
/*	str - String to remove comment from				*/
/*									*/
/*	This function removes comments from the end of a line that	*/
/*   are of the form:							*/
/*									*/
/*	sample line		! This is a comment			*/
/*									*/
/************************************************************************/

void NAMD_remove_comment(char *str)

{
	int i=0;

	while ( (str[i] != STRINGNULL) && (str[i] != '!') )
	{
		i++;
	}

	str[i] = STRINGNULL;
}
/*			END OF FUNCTION NAMD_remove_comment		*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_truncate				*/
/*									*/
/*   INPUTS:								*/
/*	str - string to truncate					*/
/*									*/
/*      NAMD_truncate will remove any trailing spaces from a string.    */
/*   i.e.  "AB DF FG     "  would be truncated to "AB DF FG".		*/
/*									*/
/************************************************************************/

void NAMD_truncate(char *str)

{
	int i;		//  Current position in str

	i=strlen(str);

	/*  Loop from the back of the string until we find a non-space  */
	for (i--; i>=0 && isspace(str[i]); i--)
	{
	}
	
	str[i+1]=STRINGNULL;
}
/*			END OF FUNCTION NAMD_truncate			*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_find_word				*/
/*									*/
/*   INPUTS:								*/
/*	source - the string to be searched in				*/
/*	search - the string to be searched for				*/
/*									*/
/*   OUTPUTS:								*/
/*	a 1 is returned if the string search is found within the string */
/*   source, otherwise a 0 is returned.					*/
/*									*/
/*	NAMD_find_word searches for one string within another.  It is   */
/*   usually used to determine if a word appears within a given line.   */
/*   If the word is found, a 1 is returned.  Otherwise, 0 is returned.  */
/*   Case is ignored while doing this search.				*/
/*									*/
/************************************************************************/

int NAMD_find_word(char *source, char *search)

{
	int i=0;		//  Position inside source
	int search_length;	//  Length of string search
	int source_length;	//  Length of string source
	int found=0;		//  Flag 1-> found the value

	search_length=strlen(search);
	source_length=strlen(source);

	/*  While we haven't found it and we haven't readched the      */
	/*  point where our current position plus the length of the    */
	/*  string we are looking for is longer than the string itself,*/
	/*  check the next search_length characters for a match	       */
	while (!found && ((i+search_length)<=(source_length)))
	{
		found = (strncasecmp(source+i, search, search_length)==0);

		i++;
	}

	return(found);
}
/*			END OF FUNCTION NAMD_find_word			*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_blank_str				*/
/*									*/
/*   INPUTS:								*/
/*	str - string to test						*/
/*									*/
/*   OUTPUTS:								*/
/*	a 1 is returned if the string is blank, otherwise a 0 is        */
/*   returned								*/
/*									*/
/*	NAMD_blank_str tests to see if a string is blank.  That is,     */
/*   contains only characters where isspace() is true			*/
/*									*/
/************************************************************************/

int NAMD_blank_string(char *str)
{
	int i;		// Current position in str
	int blank=1;	// Flag 1-> string is blank
	int len;	// length of the string str

	len=strlen(str);

	for (i=0; i<len && blank; i++)
	{
		blank = isspace(str[i]);
	}

	return(blank);
}
/*			END OF FUNCTION NAMD_blank_string		*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_find_first_word			*/
/*									*/
/*   INPUTS:								*/
/*	source - string to obtain word from				*/
/*	word - buffer to place word into				*/
/*									*/
/*   OUTPUTS:								*/
/*	word is returned containing the first word of source		*/
/*									*/
/*	This function finds the first word in a string.  The first word */
/*   is defined to be the first set of continuous non-space charaters   */
/*   in a string.  So in the string "   AB14^  FDGFD GFDG"  the first   */
/*   word would be "AB14^".  The word is returned in the string pointed */
/*   to by word.							*/
/*									*/
/************************************************************************/

void NAMD_find_first_word(char *source, char *word)

{
	int i=0;	//  Position within source
	int j=0;	//  Position within word

	/*  Skip leading white space					*/
	while ( (source[i] != STRINGNULL) && isspace(source[i]))
		i++;

	/*  Copy the word						*/
	while ( (source[i] != STRINGNULL) && !isspace(source[i]))
	{
		word[j]=source[i];
		i++;
		j++;
	}

	word[j]=STRINGNULL;

	return;
}
/*			END OF FUNCTION NAMD_find_first_word		*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_read_int 				*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor to read from				*/
/*	msg - string indicating what we are trying to read		*/
/*									*/
/*   OUTPUTS:								*/
/*	the value of the next integer in the file is returned		*/
/*									*/
/*	NAMD_read_int is used to read in integer lists from the .psf    */
/*   file.  It will read the next integer it finds in the file passed   */
/*   to it.  If an alpha character is encountered, the program          */
/*   terminates.  If an EOF is encountered, the program terminates.     */
/*   The string msg is used to indicate what we were trying to read     */
/*   in any error messages.						*/
/*									*/
/************************************************************************/

int NAMD_read_int(FILE *fd, char *msg)

{
	int i;			//  Loop counter
	int c;			//  Character read in from file
	char tmp_string[11];	//  Temporary string for integer

	/*  Skip white space				*/
	while ( ((c=fgetc(fd)) == '\n') || isspace(c) )
	{
	}

	/*  Check to make sure we didn't hit EOF	*/
	if (c==EOF)
	{
		char err_msg[128];

		sprintf(err_msg, "EOF ENCOUNTERED READING %s FROM PSF FILE", msg);
		NAMD_die(err_msg);
	}

	/*  Now read in the integer itself		*/
	i=0;
	
	while (!isspace(c))
	{
		/*  Check to make sure we only get #'s  */
		if (!isdigit(c))
		{
			char err_msg[128];

			sprintf(err_msg, "ALPHA CHARCTER ENCOUNTERED WHILE READING %s FROM PSF FILE", msg);
			NAMD_die(err_msg);
		}

		tmp_string[i] = c;
		i++;

		c=fgetc(fd);

		/*  Check to make sure we didn't hit EOF*/
		if (c==EOF)
		{
			char err_msg[128];

			sprintf(err_msg, "EOF ENCOUNTERED WHILE READING %s FROM PSF FILE", msg);
			NAMD_die(err_msg);
		}
	}

	tmp_string[i]=STRINGNULL;

	/*  Convert the string to an integer and return its value	*/
	return(atoi(tmp_string));
}
/*			END OF FUNCTION NAMD_read_int			*/

/************************************************************************/
/*									*/
/*			FUNCTION NAMD_pad				*/
/*									*/
/*	This function pads a string with leading spaces to a specified  */
/*   length.								*/
/*									*/
/************************************************************************/

void NAMD_pad(char *str, size_t length)

{
	char tmp_string[128];
	size_t i;

	if (strlen(str) >= length)
		return;

	for (i=0; i<(length-strlen(str)); i++)
	{
		tmp_string[i] = ' ';
	}

	tmp_string[i] = STRINGNULL;

	strcat(str, tmp_string);
}
/*			END OF FUNCTION NAMD_pad			*/

//  For AIX, implement the library routines strcasecmp and strncasecmp, since
//  for some bizarre reason, IBM hasn't seen fit to!

//  The implementation is a modified version of the strcmp shown
//  in the K&R C manual

#ifdef SP2
int strcasecmp(const char s[], const char t[])

{
	int i=0;

	while (tolower(s[i]) == tolower(t[i]))
		if (s[i++] == STRINGNULL)
			return(0);

	return(tolower(s[i]) - tolower(t[i]));
}

int strncasecmp(const char s[], const char t[], int n)

{
	int i=0;

	while (tolower(s[i]) == tolower(t[i]))
		if ( (s[i++] == STRINGNULL) || (i==n) )
			return(0);

	return(tolower(s[i]) - tolower(t[i]));
}
#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: strlib.C,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1003 $	$Date: 1998/10/24 19:58:10 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: strlib.C,v $
 * Revision 1.1003  1998/10/24 19:58:10  jim
 * Eliminated warnings generated by g++ -Wall.
 *
 * Revision 1.1002  1997/07/09 21:26:46  milind
 * Ported NAMD2 to SP3. The SP specific code is within #ifdef SP2
 * and #endif's.
 *
 * Revision 1.1001  1997/03/19 11:55:04  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 * Revision 1.1000  1997/02/06 15:59:38  ari
 * Resetting CVS to merge branches back into the main trunk.
 * We will stick to main trunk development as suggested by CVS manual.
 * We will set up tags to track fixed points of development/release
 * as suggested by CVS manual - all praise the CVS manual.
 *
 * Revision 1.778  1997/01/28 00:31:39  ari
 * internal release uplevel to 1.778
 *
 * Revision 1.777  1997/01/17 19:37:15  ari
 * Internal CVS leveling release.  Start development code work
 * at 1.777.1.1.
 *
 * Revision 1.1  1996/08/06 20:38:38  ari
 * Initial revision
 *
 * Revision 1.7  1995/09/21 17:45:51  billh
 * Use `\0' instead of NULL when terminating a string (or comparing a
 * character to the string terminator character ... yes they're both 0,
 * but NULL is a void *, and '\0' is a char)
 *
 * Revision 1.6  95/04/10  11:30:50  11:30:50  nelson (Mark T. Nelson)
 * Added strcasecmp and strncasecmp for AIX
 * 
 * Revision 1.5  95/03/08  14:46:18  14:46:18  nelson (Mark T. Nelson)
 * Added copyright
 * 
 * Revision 1.4  95/01/26  14:14:12  14:14:12  nelson (Mark T. Nelson)
 * Added function NAMD_remove_comment for charmm22 parameters
 * 
 * Revision 1.3  94/10/28  12:49:44  12:49:44  nelson (Mark T. Nelson)
 * Fixed bug in NAMD_pad
 * 
 * Revision 1.2  94/09/30  09:09:59  09:09:59  nelson (Mark T. Nelson)
 * added NAMD_pad function
 * 
 * Revision 1.1  94/06/22  15:05:44  15:05:44  nelson (Mark T. Nelson)
 * Initial revision
 * 
 ***************************************************************************/
