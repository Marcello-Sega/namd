/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"
#include "InfoStream.h"

#include "main.decl.h"
#include "main.h"


#define TBSOFT_TRACK_HOST   "130.126.120.106" /* www.ks.uiuc.edu */
#define TBSOFT_TRACK_PORT   3141              /* UDP port 3141   */
#define TBSOFT_TRACK_MAXLEN 1024              /* maximum message length */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <arpa/inet.h>
#include <fcntl.h>              
#include <sys/types.h>
#include <unistd.h>
#include <sys/socket.h>


int send_dgram(const char *host_addr, int port, const char *buf, int buflen) {
  struct sockaddr_in addr;
  int sockfd;

  if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
    return -1;
  } 

  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_port = htons(port);
  addr.sin_addr.s_addr = inet_addr(host_addr);

  sendto(sockfd, buf, buflen, 0, (struct sockaddr *)&addr, sizeof(addr));

  close(sockfd);
  return 0;
}                     


int tbsoft_sendusage(const char *program, 
                     const char *versionnum,
                     const char *platform,
                     const char *numcpus,
                     const char *miscinfo) {

  char sendbuf[TBSOFT_TRACK_MAXLEN];
  char host[1024];
  char user[1024];

  memset(sendbuf, 0, sizeof(sendbuf));

  gethostname(host, 1023);
#if defined(__linux)
  strcpy(user, getlogin());
#else
  cuserid(user);
#endif

  sprintf(sendbuf, "1 %s  %s  %s  %s  %s  %s  %s", 
    program, versionnum, platform, numcpus, miscinfo, host, user);
  send_dgram(TBSOFT_TRACK_HOST, TBSOFT_TRACK_PORT, sendbuf, strlen(sendbuf));

  return 0;
}


class main : public Chare
{
public:
  main(CkArgMsg *)
  {
    // print banner
    iout << iINFO << "NAMD " << NAMD_VERSION << " for " << NAMD_PLATFORM
#ifdef NAMD_FFTW
         << " (includes FFTW, do not distribute)"
#endif
         << "\n"
#if 1
         << iWARN << "\n"
         << iWARN << "          ***  UNRELEASED EXPERIMENTAL VERSION  ***\n"
         << iWARN << "\n"
#endif
         << iINFO << "Please complete the registration form at\n"
         << iINFO << "http://www.ks.uiuc.edu/Research/namd/download.html\n"
         << iINFO << "and send feedback or bug reports to namd@ks.uiuc.edu\n"
         << endi;

    iout << iINFO
      << "Reporting usage information via UDP to NAMD developers.\n" << endi;
    char numcpus[512];
    sprintf(numcpus,"%d",CmiNumPes());
    tbsoft_sendusage("NAMD",NAMD_VERSION,NAMD_PLATFORM,numcpus,"");
  }

};

#include "main.def.h"

