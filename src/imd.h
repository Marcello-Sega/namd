
/* 
 *
 * IMDHeader: The idea behind the messages is that first you send a header 
 * message consisting of a fixed number of chars; that message tells you what
 * kinf of message you're getting next, including the all the type information
 * needed to byteswap, sort into structs, allocate memory, etc.
 */

typedef struct {
  char type[8];      /* Type of message  */
  char length[8];    /* Some length parameter, not necessarily the number  */
                     /* of bytes in the data message.  */
  char size[8];      /* The number of bytes in the following data message  */
} IMDHeader;

typedef enum {
  TEXT,      /* TEXT */
  FCOORDS,   /* FCOO */
  MDCOMM,    /* MDCO */
  ENERGIES,  /* NRGS */
  EXIT,      /* EXIT */
  ERROR      /* XXXX */
} IMDHeaderType;
  
extern void imd_setheader(IMDHeader *, IMDHeaderType, int length, int size);
extern void imd_getheader(const IMDHeader *,IMDHeaderType *, int *, int *);
extern int  imd_sendheader(void *, IMDHeaderType,int, int);
extern int  imd_readheader(void *, IMDHeaderType *, int *, int *);

extern int  imd_blockwrite(void *, const char *, int);
extern int  imd_blockread( void *, char *, int);

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

