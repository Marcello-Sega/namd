
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "imd.h"
#include "vmdsock.h"

void imd_setheader(IMDHeader *header, IMDHeaderType type, int length, int sz) {
  sprintf(header->length, "%d", length);
  sprintf(header->size, "%d", sz);
  switch (type) {
    case TEXT:
      strcpy(header->type, "TEXT");
      break;
    case MDCOMM:
      strcpy(header->type, "MDCO");
      break;
    case FCOORDS:
      strcpy(header->type, "FCOO");
      break;
    case EXIT:
      strcpy(header->type, "EXIT");
      break;
    default:
      strcpy(header->type, "XXXX");
  }
}

void imd_getheader(const IMDHeader *header, IMDHeaderType *type, int *length,
                   int *sz) {
  *length=atoi(header->length);
  *sz = atoi(header->size); 
  if (!strncmp(header->type, "TEXT", 4))
    *type = TEXT;
  else if (!strncmp(header->type, "MDCO",4))
    *type = MDCOMM;
  else if (!strncmp(header->type, "FCOO", 4))
    *type = FCOORDS;
  else if (!strncmp(header->type, "EXIT", 4))
    *type = EXIT;
  else            
    *type = ERROR;   
} 

int imd_sendheader(void *s, IMDHeaderType type, int length, int size) {
  IMDHeader header;
  imd_setheader(&header, type, length, size);
  return imd_blockwrite(s, (char *)&header, sizeof(IMDHeader));
}

int imd_readheader(void * s, IMDHeaderType *type, int *length, int *size) {
  IMDHeader header;
  int rc = imd_blockread(s, (char *)&header, sizeof(IMDHeader));
  if (rc < 0) return rc;
  imd_getheader(&header, type, length, size);
  return rc;
}

int imd_blockwrite(void *s, const char *buf, int ntotal) {
  int rc;
  int ndone = 0;
  if (s == NULL) return -1;
  do {
    rc = vmdsock_write(s, buf+ndone, ntotal-ndone);
    if (rc < 0) break;
    ndone += rc;
  } while (ndone < ntotal);
  return ndone;
}

int imd_blockread(void *s, char *buf, int ntotal) {
  int rc;
  int ndone = 0;
  if (s == NULL) return -1;
  do {
    rc = vmdsock_read(s, buf+ndone, ntotal-ndone);
    if (rc < 0) break;
    ndone += rc;
  } while (ndone < ntotal);
  return ndone;
}

 
