#include <stdio.h>
#include "common.h"
#include "Inform.h"
#include "string.h"

Inform::Inform(char *str)
{
   strncpy(name,str,15);
   name[15] = '\0';
}

Inform::~Inform(){}


Inform& Inform::operator<<(const char *str)
{

    sprintf(tmp_str,"%s",str);
    n = MAX_MSG_SIZE - strlen(buffer)-1; 
    strncat(buffer,tmp_str,n);
    return *this;
}

Inform& Inform::operator<<(char c)
{
    sprintf(tmp_str,"%c",c);
    n = MAX_MSG_SIZE - strlen(buffer)-1;
    strncat(buffer,tmp_str,n); 
    return *this;
}

Inform& Inform::operator<<(short i)
{
    sprintf(tmp_str,"%d",i);
    n = MAX_MSG_SIZE - strlen(buffer)-1;
    strncat(buffer,tmp_str,n);
    return *this;
}

Inform& Inform::operator<<(int i)
{
    sprintf(tmp_str,"%d",i);
    n = MAX_MSG_SIZE - strlen(buffer)-1;
    strncat(buffer,tmp_str,n); 
    return *this;
}


Inform& Inform::operator<<(unsigned int i)
{
    sprintf(tmp_str,"%d",i);
    n = MAX_MSG_SIZE - strlen(buffer)-1;
    strncat(buffer,tmp_str,n);
    return *this;
}

Inform& Inform::operator<<(long i)
{
    sprintf(tmp_str,"%ld",i);
    n = MAX_MSG_SIZE - strlen(buffer)-1;
    strncat(buffer,tmp_str,n);
    return *this;
}

Inform& Inform::operator<<(float f)
{
    sprintf(tmp_str,"%f",f);
    n = MAX_MSG_SIZE - strlen(buffer)-1;
    strncat(buffer,tmp_str,n);
    return *this;
}

Inform& Inform::operator<<(double f)
{
    sprintf(tmp_str,"%lf",f);
    n = MAX_MSG_SIZE - strlen(buffer)-1;
    strncat(buffer,tmp_str,n);
    return *this;
}


void Inform::display_message()
{
    printf("Node%d:%s>%s\n",0,name,buffer);
    buffer[0] = '\0';
}



Inform& sendmsg(Inform& inform)
{
    inform.display_message();
    return inform;
}


