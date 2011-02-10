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
#if defined(WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define CHDIR _chdir
#define GETCWD _getcwd
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#else
#include <unistd.h>
#define CHDIR chdir
#define GETCWD getcwd
#define PATHSEP '/'
#define PATHSEPSTR "/"
#endif
#include <sys/stat.h>
#include "ConfigList.h"
#include "ScriptTcl.h"


void after_backend_init(int argc, char **argv);

#ifdef MEM_OPT_VERSION
//record the working directory when reading the configuration file
//for Parallel IO Input --Chao Mei
char *gWorkDir = NULL;
#endif

int main(int argc, char **argv) {
  BackEnd::init(argc,argv);
  after_backend_init(argc, argv);
  return 0;
}

void after_backend_init(int argc, char **argv){
  ScriptTcl *script = new ScriptTcl;
  Node::Object()->setScript(script);

  for(argc = 0; argv[argc]; ++argc);
  if ( argc < 2 ) {
    NAMD_die("No simulation config file specified on command line.");
  }
  char *origcwd = GETCWD(0,0);
#ifdef NAMD_TCL
  for(int i = 1; i < argc; ++i) {
  if ( strstr(argv[i],"--") == argv[i] ) {
    char buf[1024];
    if ( i + 1 == argc ) {
      sprintf(buf, "missing argument for command line option %s", argv[i]);
      NAMD_die(buf);
    }
    sprintf(buf, "%s %s", argv[i]+2, argv[i+1]);
    iout << iINFO << "Command-line argument is --" << buf << "\n" << endi;
    script->eval(buf);
    ++i;
    continue;
  }
  char *confFile = argv[i];
#else
  char *confFile = argv[argc-1];
#endif
  iout << iINFO << "Configuration file is " << confFile << "\n" << endi;

  char *currentdir=confFile;
  char *tmp;
  for(tmp=confFile;*tmp;++tmp); // find final null
  for( ; tmp != confFile && *tmp != PATHSEP; --tmp); // find last '/'
#if defined(WIN32) && !defined(__CYGWIN__)
  if (tmp == confFile) {
    // in case this is under cygwin, search for '/' as well
    for(tmp=confFile;*tmp;++tmp); // find final null
    for( ; tmp != confFile && *tmp != '/'; --tmp); // find last '/'
  }
#endif
  if ( CHDIR(origcwd) ) NAMD_err(origcwd);
  if ( tmp != confFile )
  {
    *tmp = 0; confFile = tmp + 1;
    if ( CHDIR(currentdir) ) NAMD_err(currentdir);
    iout << iINFO << "Changed directory to " << currentdir << "\n" << endi;
    currentdir = GETCWD(0,0);
  }
  else{
      if ( *tmp == PATHSEP ){ // config file in / is odd, but it might happen
          if ( CHDIR(PATHSEPSTR) ) NAMD_err(PATHSEPSTR);
      }else{ // just a config file name, so the path is the current working path
          char tmpcurdir[3];
          tmpcurdir[0] = '.';
          tmpcurdir[1] = PATHSEP;
          tmpcurdir[2] = 0;
          currentdir = tmpcurdir;
          iout << iINFO << "Working in the current directory " << origcwd << "\n" << endi;
      }
  }

#ifdef MEM_OPT_VERSION
    int dirlen = strlen(currentdir);
    gWorkDir = new char[dirlen+1];
    gWorkDir[dirlen]=0;
    memcpy(gWorkDir, currentdir, dirlen);
#endif

  currentdir = NULL;

  struct stat statBuf;
  if (stat(confFile, &statBuf)) {
    NAMD_die("Simulation config file is not accessible.");
  }

#ifdef NAMD_TCL
  script->load(confFile);
#else
  script->run(confFile);
#endif

#ifdef NAMD_TCL
}
  script->run();
#endif

  BackEnd::exit();
}

