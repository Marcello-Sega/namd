/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "memusage.h"
#include "converse.h"
#include "common.h"
#include "BackEnd.h"
#include "InfoStream.h"
#include "Broadcasts.h"

#include "NamdState.h"
#include "Node.h"
#ifdef WIN32
#include <direct.h>
#define CHDIR _chdir
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#else
#include <unistd.h>
#define CHDIR chdir
#define PATHSEP '/'
#define PATHSEPSTR "/"
#endif
#include <sys/stat.h>
#include "ConfigList.h"
#include "ScriptTcl.h"

int main(int argc, char **argv) {
  BackEnd::init(argc,argv);
  for(argc = 0; argv[argc]; ++argc);
  if ( argc < 2 ) {
    NAMD_die("No simulation config file specified on command line.");
  }

  char *confFile = argv[argc-1];
  char *currentdir=confFile;
  char *tmp;
  for(tmp=confFile;*tmp;++tmp); // find final null
  for( ; tmp != confFile && *tmp != PATHSEP; --tmp); // find last '/'
  if ( tmp != confFile )
  {
    *tmp = 0; confFile = tmp + 1;
    if ( CHDIR(currentdir) ) NAMD_die("chdir() failed!");
    iout << iINFO << "Changed directory to " << currentdir << "\n" << endi;
  }
  else if ( *tmp == PATHSEP ) // config file in / is odd, but it might happen
    if ( CHDIR(PATHSEPSTR) ) NAMD_die("chdir() failed!");
  currentdir = NULL;

  iout << iINFO << "Configuration file is " << confFile << "\n" << endi;

  struct stat statBuf;
  if (stat(confFile, &statBuf)) {
    NAMD_die("Simulation config file is not accessible.");
  }

  ScriptTcl *script = new ScriptTcl;
  Node::Object()->setScript(script);
  script->run(confFile,0);

  BackEnd::exit();
  return 0;
}

