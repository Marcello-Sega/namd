
#ifndef IMD_H__
#define IMD_H__

#include <limits.h>

#if ( INT_MAX == 2147483647 )
typedef int     int32;
#else
typedef short   int32;
#endif

enum IMDType {
  IMD_DISCONNECT,
  IMD_ENERGIES, 
  IMD_FCOORDS,   
  IMD_HANDSHAKE, 
  IMD_KILL,      
  IMD_MDCOMM,    
  IMD_PAUSE,
  IMD_TRATE,
  IMD_IOERROR
};

typedef struct {
  int32 tstep;
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

// Send simple messages - these consist of a header with no subsequent data
extern int   imd_disconnect(void *);
extern int   imd_pause(void *);
extern int   imd_kill(void *);
extern int   imd_handshake(void *);
extern int   imd_trate(void *, int32);

// Send data
extern int   imd_send_mdcomm(void *, int32, const int32 *, const float *);
extern int   imd_send_energies(void *, const IMDEnergies *);
extern int   imd_send_fcoords(void *, int32, const float *);

/// Receive header and data 
extern IMDType imd_recv_header(void *, int32 *);
extern int imd_recv_mdcomm(void *, int32, int32 *, float *);
extern int imd_recv_energies(void *, IMDEnergies *);
extern int imd_recv_fcoords(void *, int32, float *);
 
#endif

