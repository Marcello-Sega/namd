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

   Inform(char *);

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
