
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "imd.h"
#include "vmdsock.h"

static const char *headers[] = {
  "NRGS",  
  "EXIT",
  "FCOO",
  "KILL",
  "MDCO",
  "PAUS",
  "TEXT",
  "TRTE",
  "XXXX"
};

void imd_setheader(IMDHeader *header, IMDHeaderType type, int length, int sz) {
  sprintf(header->length, "%d", length);
  sprintf(header->size, "%d", sz);
  strcpy(header->type, headers[type]); 
}

void imd_getheader(const IMDHeader *header, IMDHeaderType *type, int *length,
                   int *sz) {
  *length=atoi(header->length);
  *sz = atoi(header->size); 
  if (!strncmp(header->type, headers[IMD_ENERGIES], 4)) 
    *type = IMD_ENERGIES; 
  else if (!strncmp(header->type, headers[IMD_EXIT], 4)) 
    *type = IMD_EXIT; 
  else if (!strncmp(header->type, headers[IMD_FCOORDS], 4)) 
    *type = IMD_FCOORDS; 
  else if (!strncmp(header->type, headers[IMD_KILL], 4)) 
    *type = IMD_KILL; 
  else if (!strncmp(header->type, headers[IMD_MDCOMM], 4)) 
    *type = IMD_MDCOMM; 
  else if (!strncmp(header->type, headers[IMD_PAUSE], 4)) 
    *type = IMD_PAUSE; 
  else if (!strncmp(header->type, headers[IMD_TEXT], 4)) 
    *type = IMD_TEXT; 
  else if (!strncmp(header->type, headers[IMD_TRATE], 4)) 
    *type = IMD_TRATE; 
  else
    *type = IMD_XXXX; 
}  

int imd_sendheader(void *s, IMDHeaderType type, int length, int size) {
  IMDHeader header;
  imd_setheader(&header, type, length, size);
  return imd_writen(s, (char *)&header, sizeof(IMDHeader));
}

int imd_readheader(void * s, IMDHeaderType *type, int *length, int *size) {
  IMDHeader header;
  int rc = imd_readn(s, (char *)&header, sizeof(IMDHeader));
  if (rc < sizeof(IMDHeader)) return -1;
  imd_getheader(&header, type, length, size);
  return rc;
}

ssize_t imd_readn(void *s, char *ptr, size_t n) {
  size_t nleft;
  ssize_t nread;
  
  nleft = n;
  while (nleft > 0) {
    if ((nread = vmdsock_read(s, ptr, nleft)) < 0) {
      if (errno == EINTR)
        nread = 0;         /* and call read() again */
      else
        return -1;
    } else if (nread == 0)
      break;               /* EOF */
    nleft -= nread;
    ptr += nread;
  }
  return n-nleft;
}

ssize_t imd_writen(void *s, const char *ptr, size_t n) {
  size_t nleft;
  ssize_t nwritten;
  
  nleft = n;
  while (nleft > 0) {
    if ((nwritten = vmdsock_write(s, ptr, nleft)) <= 0) {
      if (errno == EINTR)
        nwritten = 0;
      else
        return -1;
    }
    nleft -= nwritten;
    ptr += nwritten;
  }
  return n;
}

