//-*-c++-*-
#ifndef INFORMCHARM_H
#define INFORMCHARM_H

#define MAX_MSG_SIZE    10000

class Inform {

private :

   int On;
   int n;
   char name[16];
   char buffer[MAX_MSG_SIZE]; 
   char tmp_str[MAX_MSG_SIZE];

public  :

   Inform(const char *);

   ~Inform(void);

   void on(int o) {On = o;}

   int  on(void)  {return On;}
 
   Inform& operator<<(const char *);
   Inform& operator<<(char);
   Inform& operator<<(short);
   Inform& operator<<(int);
   Inform& operator<<(unsigned int);
   Inform& operator<<(long);
   Inform& operator<<(float);
   Inform& operator<<(double);
   Inform& operator<<(Inform& (*f)(Inform&)) {return f(*this);}
   void display_message();

};

extern Inform& sendmsg(Inform&);

#endif


/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1002 $	$Date: 1999/07/06 20:32:44 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Inform.h,v $
 * Revision 1.1002  1999/07/06 20:32:44  jim
 * Eliminated warnings from new generation of picky compilers.
 *
 * Revision 1.1001  1997/03/19 11:54:20  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
