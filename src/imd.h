/* 
 *
 * IMDHeader: The idea behind the messages is that first you send a header 
 * message consisting of a fixed number of chars; that message tells you what
 * kinf of message you're getting next, including the all the type information
 * needed to byteswap, sort into structs, allocate memory, etc.
 */

#include <sys/types.h>  /* For types size_t, ssize_t, etc.  */

typedef struct {
  char type[8];      /* Type of message  */
  char length[8];    /* Some length parameter, not necessarily the number  */
                     /* of bytes in the data message.  */
  char size[8];      /* The number of bytes in the following data message  */
} IMDHeader;

enum IMDHeaderType {
  IMD_ENERGIES,  /* NRGS */
  IMD_ERROR,     /* XXXX */
  IMD_EXIT,      /* EXIT */
  IMD_FCOORDS,   /* FCOO */
  IMD_HANDSHAKE, /* HAND */
  IMD_KILL,      /* KILL */
  IMD_MDCOMM,    /* MDCO */
  IMD_PAUSE,     /* PAUS */
  IMD_TEXT,      /* TEXT */
  IMD_TRATE      /* TRTE */
};
typedef enum IMDHeaderType IMDHeaderType;
 
extern int  imd_sendheader(void *, IMDHeaderType,int, int);
extern int  imd_readheader(void *, IMDHeaderType *, int *, int *);

/* readn and writen take from Stevens 2nd ed., p. 78. */
extern ssize_t imd_readn(void *, char *, size_t);
extern ssize_t imd_writen(void *, const char *, size_t);
 
/*
 * Data structures sent by NAMD to VMD
 *
 */

typedef struct {
  float T;
  float Etot;
  float Epot;
  float Evdw;
  float Eelec;
  float Ebond;
  float Eangle;
  float Edihe;
  float Eimpr;
} IMDEnergies;

