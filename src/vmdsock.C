/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#if defined(DUMMY_VMDSOCK)

#include "vmdsock.h"

void * vmdsock_create(void) { return 0; }
int  vmdsock_connect(void *v, const char *host, int port) { return 0; }
int vmdsock_bind(void * v, int port) { return 0; }
int vmdsock_listen(void * v) { return 0; }
int  vmdsock_accept(void * v) { return 0; }
int  vmdsock_write(void * v, const void *buf, int len) { return 0; }
int  vmdsock_read(void * v, void *buf, int len) { return 0; }
void vmdsock_destroy(void * v) { return; }
int vmdsock_selread(void *v, int sec) { return 0; }
int vmdsock_selwrite(void *v, int sec) { return 0; }

#else

#define VMDSOCKINTERNAL

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _AIX
#include <strings.h>
#endif
#include <arpa/inet.h>
#include <fcntl.h>

#include <unistd.h>   /* for Linux */
#include <sys/socket.h>
#include <netdb.h>
#include <errno.h>

#include "vmdsock.h"

void * vmdsock_create(void) {
  vmdsocket * s;

  s = (vmdsocket *) malloc(sizeof(vmdsocket));
  if (s != NULL)
    memset(s, 0, sizeof(vmdsocket)); 

  if ((s->sd = socket(PF_INET, SOCK_STREAM, 0)) == -1) {
    printf("Failed to open socket.");
    free(s);
    return NULL;
  }

  return (void *) s;
}

int  vmdsock_connect(void *v, const char *host, int port) {
  vmdsocket *s = (vmdsocket *) v;
  char address[1030];
  struct hostent *h;

  h=gethostbyname(host);
  if (h == NULL) 
    return -1;
  sprintf(address, "%d.%d.%d.%d",
    (unsigned char) h->h_addr_list[0][0],
    (unsigned char) h->h_addr_list[0][1],
    (unsigned char) h->h_addr_list[0][2],
    (unsigned char) h->h_addr_list[0][3]);

  memset(&(s->addr), 0, sizeof(s->addr)); 
  s->addr.sin_family = PF_INET;
  s->addr.sin_addr.s_addr = inet_addr(address);
  s->addr.sin_port = htons(port);  

  return connect(s->sd, (struct sockaddr *) &s->addr, sizeof(s->addr)); 
}

int vmdsock_bind(void * v, int port) {
  vmdsocket *s = (vmdsocket *) v;
  memset(&(s->addr), 0, sizeof(s->addr)); 
  s->addr.sin_family = PF_INET;
  s->addr.sin_port = htons(port);

  return bind(s->sd, (struct sockaddr *) &s->addr, sizeof(s->addr));
}

int vmdsock_listen(void * v) {
  vmdsocket *s = (vmdsocket *) v;
  return listen(s->sd, 5);
}

int  vmdsock_accept(void * v) {
  int rc;
  vmdsocket *s = (vmdsocket *) v;
#if defined __linux__ || defined _AIX
  socklen_t len;
#else
  int len;
#endif
  len = sizeof(s->addr);
  rc = accept(s->sd, (struct sockaddr *) &s->addr, &len);
  if (rc >= 0) s->sd = rc;
  return rc;
}

int  vmdsock_write(void * v, const void *buf, int len) {
  vmdsocket *s = (vmdsocket *) v;
  return write(s->sd, buf, len);
}

int  vmdsock_read(void * v, void *buf, int len) {
  vmdsocket *s = (vmdsocket *) v;
  return read(s->sd, buf, len);
}

void vmdsock_destroy(void * v) {
  vmdsocket * s = (vmdsocket *) v;
  if (s == NULL)
    return;

  close(s->sd);
  free(s);  
}

int vmdsock_selread(void *v, int sec) {
  vmdsocket *s = (vmdsocket *)v;
  // struct fd_set rfd;
  fd_set rfd;
  struct timeval tv;
  int rc;
 
  FD_ZERO(&rfd);
  FD_SET(s->sd, &rfd);
  memset((void *)&tv, 0, sizeof(struct timeval));
  tv.tv_sec = sec;
  do {
    rc = select(s->sd+1, &rfd, NULL, NULL, &tv);
  } while (rc < 0 && errno == EINTR);
  return rc;

}
  
int vmdsock_selwrite(void *v, int sec) {
  vmdsocket *s = (vmdsocket *)v;
  // struct fd_set wfd;
  fd_set wfd;
  struct timeval tv;
  int rc;
 
  FD_ZERO(&wfd);
  FD_SET(s->sd, &wfd);
  memset((void *)&tv, 0, sizeof(struct timeval));
  tv.tv_sec = sec;
  do {
    rc = select(s->sd + 1, NULL, &wfd, NULL, &tv);
  } while (rc < 0 && errno == EINTR);
  return rc;
}

#endif
