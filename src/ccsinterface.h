#include <stdio.h>
#include <string.h>
#include "converse.h"
#include <errno.h>

#if CMK_WHEN_PROCESSOR_IDLE_USLEEP
#include <sys/types.h>
#include <sys/time.h>
#endif

#if CMK_TIMER_USE_TIMES
#include <sys/times.h>
#include <limits.h>
#include <unistd.h>
#endif

#if CMK_TIMER_USE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#endif

#ifdef NAMDCCS

#ifdef __cplusplus
extern "C" {
#endif
void CApplicationInit(void);
void CApplicationDepositData(char *data);
void CApplicationDepositNode0Data(char *data);
#ifdef __cplusplus
}
#endif

#endif
