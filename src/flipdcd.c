/********************************************************************/
/*                                                                  */
/*   Copyright 1998, Jim Phillips and the University of Illinois.   */
/*                                                                  */
/********************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>

#ifndef MAP_FILE
#define MAP_FILE 0
#endif

int main(int argc, char *argv[]) {

int fd;
struct stat statbuf;
int i, j, n, isbig, itmp;
char b[8];
char *d;

if ( argc != 2 ) {
  fprintf(stderr,"This program flips byte-ordering of DCD files.\n");
  fprintf(stderr,"Usage: %s <filename>\n",argv[0]);
  exit(-1);
}

if ( ( fd = open(argv[1], O_RDWR) ) < 0 ) {
  fprintf(stderr,"Can't open %s for updating.\n",argv[1]);
  exit(-1);
}

if ( fstat(fd,&statbuf) < 0 ) {
  fprintf(stderr,"Can't stat %s.\n",argv[1]);
  exit(-1);
}

n = statbuf.st_size;

if ( n <= 104 ) {
  fprintf(stderr,"%s is not in DCD format.\n",argv[1]);
  exit(-1);
}

if ( n % 4 ) {
  fprintf(stderr,"%s is not in DCD format.\n",argv[1]);
  exit(-1);
}
if ( ( d = mmap(0,n,PROT_READ|PROT_WRITE,MAP_FILE|MAP_SHARED,fd,0) )
							== (caddr_t) -1 ) {
  fprintf(stderr,"Can't mmap %s.\n",argv[1]);
  exit(-1);
}

if ( d[0] == 84 ) {
  isbig = 0;
  fprintf(stderr,"%s was little-endian, will be big-endian.\n",argv[1]);
}
else if ( d[3] == 84 ) {
  isbig = 1;
  fprintf(stderr,"%s was big-endian, will be little-endian.\n",argv[1]);
}
else {
  fprintf(stderr,"%s is not in DCD format.\n",argv[1]);
  exit(-1);
}

#define FLIPFOUR {for(j=0;j<4;++j)b[j]=d[j];for(j=3;j>=0;--j,++d)*d=b[j];n-=4;}
#define FLIPEIGHT {for(j=0;j<8;++j)b[j]=d[j];for(j=7;j>=0;--j,++d)*d=b[j];n-=8;}
#define SKIPFOUR {d+=4;n-=4;}
#define SKIP(X) {d+=(X);n-=(X);}
#define READINT(X) { X=0; if (isbig) { for(j=0;j<4;++j,X<<8) X+=d[j]; } \
	else { for(j=3;j>=0;--j,X<<8) X+=d[j]; } }


FLIPFOUR;  /* 84 */
SKIPFOUR;  /* "CORD" */
for ( i = 0; i < 9; ++i ) FLIPFOUR;
FLIPEIGHT;  /* DELTA */
for ( i = 0; i < 9; ++i ) FLIPFOUR;
FLIPFOUR;  /* 84 */
FLIPFOUR;  /* TITLE SIZE */
READINT(itmp); FLIPFOUR;  /* NTITLE */
if ( n <= (80*itmp + 4) ) {
  fprintf(stderr,"%s is too short for DCD format.\n",argv[1]);
  exit(-1);
}
SKIP(80*itmp);
FLIPFOUR;  /* TITLE SIZE */

while ( n ) FLIPFOUR;  /* WHATEVER UNTIL END */

exit(0);

}

