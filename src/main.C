/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"
#include "InfoStream.h"
#include "memusage.h"

#include "main.decl.h"
#include "main.h"

#ifndef WIN32

#define TBSOFT_TRACK_HOST   "130.126.120.106" /* www.ks.uiuc.edu */
#define TBSOFT_TRACK_PORT   3141              /* UDP port 3141   */
#define TBSOFT_TRACK_MAXLEN 1024              /* maximum message length */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <pwd.h>


int send_dgram(const char *host_addr, int port, const char *buf, int buflen) {
  struct sockaddr_in addr;
  int sockfd;

#ifndef NOHOSTNAME
  if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
    return -1;
  } 

  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_port = htons(port);
  addr.sin_addr.s_addr = inet_addr(host_addr);

  sendto(sockfd, buf, buflen, 0, (struct sockaddr *)&addr, sizeof(addr));

  close(sockfd);
#endif

  return 0;
}                     


int tbsoft_sendusage(const char *program, 
                     const char *versionnum,
                     const char *platform,
                     const char *numcpus,
                     const char *miscinfo) {

#ifndef NOHOSTNAME
  iout << iINFO << "Sending usage information to NAMD developers via UDP.\n" << endi;

  char sendbuf[TBSOFT_TRACK_MAXLEN];
  char host[128];
  struct passwd *pw;
  char user[128];

  memset(sendbuf, 0, sizeof(sendbuf));

  gethostname(host, 128);  host[127] = 0;
  pw = getpwuid(getuid());
  if ( pw && pw->pw_name ) {
    strncpy(user, pw->pw_name, 127);  user[127] = 0;
  } else {
    sprintf(user,"%d",getuid());
  }

  sprintf(sendbuf, "1 %s  %s  %s  %s  %s  %s  %s", 
    program, versionnum, platform, numcpus, miscinfo, host, user);
  iout << iINFO << "Sent data is: " << sendbuf << "\n" << endi;
  send_dgram(TBSOFT_TRACK_HOST, TBSOFT_TRACK_PORT, sendbuf, strlen(sendbuf));

#endif
  return 0;
}

#endif

class main : public Chare
{
public:
  main(CkArgMsg *)
  {
    // print banner
    iout << iINFO << "NAMD " << NAMD_VERSION << " for " << NAMD_PLATFORM
         << "\n"
#if 0
         << iWARN << "\n"
         << iWARN << "          ***  UNRELEASED EXPERIMENTAL VERSION  ***\n"
         << iWARN << "\n"
#endif
#ifdef SCYLD_NOTICE

         << iINFO << "\n"
         << iINFO << "NAMD is a parallel, object-oriented molecular dynamics\n"
         << iINFO << "code designed for high-performance simulation of large\n"
         << iINFO << "biomolecular systems.  NAMD is distributed free of\n"
         << iINFO << "charge and includes source code.  For more information\n" 
         << iINFO << "please visit http://www.ks.uiuc.edu/Research/namd/\n"
         << iINFO << "\n"
         << iINFO << "******************************************************\n"
         << iINFO << "This version of NAMD may be distributed only as a part\n"
         << iINFO << "of the Scyld Beowulf CDROM and all other distribution\n"
         << iINFO << "is prohibited.  Any use of this software is bound by\n"
         << iINFO << "the terms of the NAMD License, which is available at\n"
         << iINFO << "http://www.ks.uiuc.edu/Research/namd/license.html\n"
         << iINFO << "The NAMD development team will not provide support for\n"
         << iINFO << "any version of NAMD unless you have first registered\n"
         << iINFO << "and downloaded the latest version of NAMD available at\n"
         << iINFO << "http://www.ks.uiuc.edu/Research/namd/\n"
         << iINFO << "******************************************************\n"
         << iINFO << "\n"
#else
         << iINFO << "Please visit http://www.ks.uiuc.edu/Research/namd/\n"
         << iINFO << "and send feedback or bug reports to namd@ks.uiuc.edu\n"
#endif
         << endi;

#ifndef WIN32
    char numcpus[512];
    sprintf(numcpus,"%d",CmiNumPes());
    tbsoft_sendusage("NAMD",NAMD_VERSION,NAMD_PLATFORM,numcpus,"");
#endif

    iout << iINFO << "Running on " << CmiNumPes() << " processors.\n" << endi;
    iout << iINFO << (memusage()/1024) << " kB of memory in use.\n" << endi;

  }

};

#include "main.def.h"

