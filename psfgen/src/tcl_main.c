
#if defined(NAMD_TCL) || ! defined(NAMD_VERSION)

#include <tcl.h>

extern int Psfgen_Init(Tcl_Interp *);

int main(int argc, char *argv[]) {
#ifdef NAMD_VERSION
  printf("PSFGEN from NAMD %s for %s\n",NAMD_VERSION,NAMD_PLATFORM);
#endif
  Tcl_Main(argc, argv, Psfgen_Init);
  return 0;
}

#else

#include <stdio.h>

int main(int argc, char **argv) {
  fprintf(stderr,"%s unavailable on this platform (no Tcl)\n",argv[0]);
  exit(-1);
}

#endif

