/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

extern "C" {
#include "converse.h"
}

#include "ProcessorPrivate.h"

extern void _initCharm(int, char**);

void charm_init(int argc, char **argv)
{
  ProcessorPrivateInit();
  _initCharm(argc, argv);
}

int main(int argc, char **argv)
{
  ConverseInit(argc, argv, charm_init,0,0);
}

