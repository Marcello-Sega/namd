
#include <stdio.h>
#include "charmm_file.h"

#define TOKLEN 100
#define BUFLEN 200

int main(int argc, char **argv) {

  char *tok[TOKLEN];
  char sbuf[BUFLEN];
  int ntok;
  int itok;

  while ( ntok = charmm_get_tokens(tok,TOKLEN,sbuf,BUFLEN,stdin) ) {
    for ( itok = 0; itok < ntok; ++itok ) {
      printf("{%s} ", tok[itok]);
    }
    printf("\n");
  }

}

