#ifndef DEBUG_H
#define DEBUG_H

#ifndef MIN_DEBUG_LEVEL
  #define MIN_DEBUG_LEVEL 0
#endif
#ifndef MAX_DEBUG_LEVEL
  #define MAX_DEBUG_LEVEL 10
#endif
#ifndef STDERR_LEVEL
  /* anything >= this error level goes to stderr */
  #define STDERR_LEVEL 5
#endif

#include <stdio.h>
#include <stdarg.h>

/*****************************************************************
 *  DebugM(): function to display a debug message.
 *  Messages have different levels.  The low numbers are low severity
 *  while the high numbers are really important.  Very high numbers
 *  are sent to stderr rather than stdout.
 *  The default severity scale is from 0 to 10.
 *     0 = plain message
 *     4 = important message
 *     5 = warning (stderr)
 *     10 = CRASH BANG BOOM error (stderr)
 *  The remaining args are like printf: a format string and some args.
 *  This function can be turned off by compiling without the DEBUGM flag
 *  No parameters to this function should have a side effect!
 *  No functions should be passed as parameters!  (including inline)
 *****************************************************************/
 #ifdef DEBUGM

  inline void DebugM(const int level, const char *format, ... )
  {
	if ((level >= MIN_DEBUG_LEVEL) && (level <= MAX_DEBUG_LEVEL))
	{
	  va_list args;
	  va_start(args, format);
	  if (level >= STDERR_LEVEL)
	  {
	    fprintf(stderr,"%s %d: ",__FILE__,__LINE__);
	    vfprintf(stderr,format,args);
	  }
	  else
	  {
	    fprintf(stdout,"%s %d: ",__FILE__,__LINE__);
	    vfprintf(stdout,format,args);
	  }
	  va_end(args);
	}
  }

 #else
  /* make a void function.  Even dumb compilers will remove this as long
     as the parameters do not have side effects and are not functions. */
  inline void DebugM(...) {;}

 #endif /* DEBUGM */

#endif /* DEBUG_H */

