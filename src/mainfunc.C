/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "converse.h"
#include "common.h"
#include "BackEnd.h"

#include "NamdState.h"
#include "Node.h"
#include <sys/stat.h>
#include "ConfigList.h"
#include "ScriptTcl.h"

int main(int argc, char **argv) {
  BackEnd::init(argc,argv);
  if ( argc < 2 ) {
    NAMD_die("No simulation config file specified on command line.");
  }

  char *confFile = argv[argc-1];
  char *currentdir=confFile;
  char *tmp;
  for(tmp=confFile;*tmp;++tmp); // find final null
  for( ; tmp != confFile && *tmp != '/'; --tmp); // find last '/'
  if ( tmp != confFile )
  {
    *tmp = 0; confFile = tmp + 1;
    if ( chdir(currentdir) ) NAMD_die("chdir() failed!");
    iout << iINFO << "Changed directory to " << currentdir << "\n" << endi;
  }
  else if ( *tmp == '/' ) // config file in / is odd, but it might happen
    if ( chdir("/") ) NAMD_die("chdir() failed!");
  currentdir = NULL;

  iout << iINFO << "Configuration file is " << confFile << "\n" << endi;

  struct stat statBuf;
  if (stat(confFile, &statBuf)) {
    NAMD_die("Simulation config file is not accessible.");
  }

  ConfigList *configList;
#ifdef NAMD_TCL
  configList = new ConfigList;  // empty, will be filled by Tcl

  ScriptTcl *script = new ScriptTcl;
  Node::Object()->setScript(script);
  script->run(confFile,configList);
#else
  if ( NULL == confFile || NULL == (configList = new ConfigList(confFile)) ) {
    NAMD_die("Simulation config file is empty.");
  }
  NamdState *state = new NamdState;
  state->configListInit(configList);
  Node::Object()->saveMolDataPointers(state);
  Node::messageStartUp();
  BackEnd::suspend();
#endif

  BackEnd::exit();
}

