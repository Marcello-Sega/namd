
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "charmm_file.h"
#include "topo_defs.h"
#include "charmm_parse_topo_defs.h"

#define TOKLEN 100
#define BUFLEN 200

static char * parse_atom(char *aref, int *res, int *rel) {
  if ( isdigit(*aref) ) { *res = *aref - '1'; ++aref; }
  else { *res = 0; }
  if ( *aref == '-' ) { *rel = -1; ++aref; }
  else if ( *aref == '+' ) { *rel = 1; ++aref; }
  else if ( *aref == '#' ) { *rel = 2; ++aref; }
  else { *rel = 0; }
  return aref;
}

static void null_print_msg(void *v, const char *s) {
  printf("%s\n", s);
}

static void debug_msg(const char *s) {
  ;
}

int charmm_parse_topo_defs(topo_defs *defs, FILE *file, void *v,
				void (*print_msg)(void *,const char *)) {

  char *tok[TOKLEN];
  char sbuf[BUFLEN];
  char msgbuf[BUFLEN];
  int ntok;
  int itok;
  int first;
  char *s1, *s2, *s3, *s4;
  int i1, i2, i3, i4; 
  int j1, j2, j3, j4; 

  first = 1;

  if ( ! defs ) return -1;
  if ( ! file ) return -2;
  if ( print_msg == 0 ) print_msg = null_print_msg;

  while ( (ntok = charmm_get_tokens(tok,TOKLEN,sbuf,BUFLEN,file)) ) {
    if ( ! tok[0][0] ) print_msg(v,tok[1]);
    else if ( first ) {
      if ( ntok == 2 ) {
        sprintf(msgbuf,"Created by CHARMM version %s %s",tok[0],tok[1]);
        print_msg(v,msgbuf);
      } else {
        print_msg(v,"ERROR!  Unusual CHARMM version record.");
      }
      first = 0;
    }
    else if ( ! strncmp("ACCE",tok[0],4) ) {
      debug_msg("Recognized acceptor statement.");
    }
    else if ( ! strncmp("DONO",tok[0],4) ) {
      debug_msg("Recognized donor statement.");
    }
    else if ( ! strncmp("BOND",tok[0],4) ||
		! strncmp("DOUB",tok[0],4) ||
		! strncmp("TRIP",tok[0],4) ) {
      debug_msg("Recognized bond statement.");
      if ( ntok < 3 || (ntok-1)%2 ) {
        print_msg(v,"ERROR!  Failed to parse bond statement.");
      } else {
        for ( itok = 1; itok < ntok; itok += 2 ) {
          s1 = parse_atom(tok[itok],&i1,&j1);
          s2 = parse_atom(tok[itok+1],&i2,&j2);
          if ( topo_defs_bond(defs,0,0,s1,i1,j1,s2,i2,j2) )
            print_msg(v,"ERROR!  Failed to parse bond statement.");
        }
      }
    }
    else if ( ! strncmp("ANGL",tok[0],4) ||
		! strncmp("THET",tok[0],4) ) {
      debug_msg("Recognized angle statement.");
      if ( ntok < 4 || (ntok-1)%3 ) {
        print_msg(v,"ERROR!  Failed to parse angle statement.");
      } else {
        for ( itok = 1; itok < ntok; itok += 3 ) {
          s1 = parse_atom(tok[itok],&i1,&j1);
          s2 = parse_atom(tok[itok+1],&i2,&j2);
          s3 = parse_atom(tok[itok+2],&i3,&j3);
          if ( topo_defs_angle(defs,0,0,s1,i1,j1,s2,i2,j2,s3,i3,j3) )
            print_msg(v,"ERROR!  Failed to parse angle statement.");
        }
      }
    }
    else if ( ! strncmp("DIHE",tok[0],4) ) {
      debug_msg("Recognized dihedral statement.");
      if ( ntok < 5 || (ntok-1)%4 ) {
        print_msg(v,"ERROR!  Failed to parse dihedral statement.");
      } else {
        for ( itok = 1; itok < ntok; itok += 4 ) {
          s1 = parse_atom(tok[itok],&i1,&j1);
          s2 = parse_atom(tok[itok+1],&i2,&j2);
          s3 = parse_atom(tok[itok+2],&i3,&j3);
          s4 = parse_atom(tok[itok+3],&i4,&j4);
          if (topo_defs_dihedral(defs,0,0,s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4))
            print_msg(v,"ERROR!  Failed to parse dihedral statement.");
        }
      }
    }
    else if ( ! strncmp("IMPH",tok[0],4) ||
		! strncmp("IMPR",tok[0],4) ) {
      debug_msg("Recognized improper statement.");
      if ( ntok < 5 || (ntok-1)%4 ) {
        print_msg(v,"ERROR!  Failed to parse improper statement.");
      } else {
        for ( itok = 1; itok < ntok; itok += 4 ) {
          s1 = parse_atom(tok[itok],&i1,&j1);
          s2 = parse_atom(tok[itok+1],&i2,&j2);
          s3 = parse_atom(tok[itok+2],&i3,&j3);
          s4 = parse_atom(tok[itok+3],&i4,&j4);
          if (topo_defs_improper(defs,0,0,s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4))
            print_msg(v,"ERROR!  Failed to parse improper statement.");
        }
      }
    }
    else if ( ! strncmp("DECL",tok[0],4) ) {
      debug_msg("Recognized atom declaration statement.");
    }
    else if ( ! strncmp("ATOM",tok[0],4) ) {
      debug_msg("Recognized atom statement.");
      if ( ntok < 4 ) {
        print_msg(v,"ERROR!  Failed to parse atom statement.");
      } else if ( ntok > 4 ) {
        print_msg(v,"ERROR!  Explicit exclustions not supported, atom ignored.");
      } else {
        s1 = parse_atom(tok[1],&i1,&j1);
        if ( topo_defs_atom(defs,0,0, s1,i1,j1,tok[2],atof(tok[3])) ) {
          print_msg(v,"ERROR!  Failed to parse atom statement.");
        }
      }
    }
    else if ( ! strncmp("MASS",tok[0],4) ) {
      debug_msg("Recognized mass statement.");
      if ( ntok < 4 || topo_defs_type(defs,tok[2],atof(tok[3]),atoi(tok[1])) ) {
        print_msg(v,"ERROR!  Failed to parse mass statement.");
      }
    }
    else if ( ! strncmp("AUTO",tok[0],4) ) {
      debug_msg("Recognized autogenerate statement.");
      for ( itok = 1; itok < ntok; itok += 1 ) {
        if ( ! strncmp("ANGL",tok[itok],4) ) {
          topo_defs_auto_angles(defs,1);
        } else if ( ! strncmp("DIHE",tok[itok],4) ) {
          topo_defs_auto_dihedrals(defs,1);
        } else {
          print_msg(v,"ERROR!  Failed to parse autogenerate statement.");
        }
      }
    }
    else if ( ! strncmp("DEFA",tok[0],4) ) {
      debug_msg("Recognized default patch statement.");
      if ( ntok < 3 || (ntok-1)%2 ) {
        print_msg(v,"ERROR!  Failed to parse default patching statement.");
      } else {
        i1 = i2 = 0;
        for ( itok = 1; itok < ntok; itok += 2 ) {
          if ( ! strncmp("FIRS",tok[itok],4) ) {
              i1 = topo_defs_default_patching_first(defs,tok[itok+1]);
          } else if ( ! strncmp("LAST",tok[itok],4) ) {
              i2 = topo_defs_default_patching_last(defs,tok[itok+1]);
          } else {
            print_msg(v,"ERROR!  Failed to parse default patching statement.");
          }
        }
        if ( i1 || i2 )
          print_msg(v,"ERROR!  Failed to parse default patching statement.");
      }
    }
    else if ( ! strncmp("BILD",tok[0],4) ||
		! strncmp("IC",tok[0],4) ) {
      debug_msg("Recognized internal coordinate statement.");
      if ( ntok < 10 ) {
        print_msg(v,"ERROR!  Failed to parse internal coordinate statement.");
      } else {
        s1 = parse_atom(tok[1],&i1,&j1);
        s2 = parse_atom(tok[2],&i2,&j2);
        s3 = tok[3] + ( tok[3][0] == '*' ? 1 : 0 );
        s3 = parse_atom(s3,&i3,&j3);
        s4 = parse_atom(tok[4],&i4,&j4);
        if ( topo_defs_conformation(defs,0,0,
		s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4,
		atof(tok[5]),atof(tok[6]),atof(tok[7]),
		(tok[3][0]=='*'?1:0),atof(tok[8]),atof(tok[9])) )
          print_msg(v,"ERROR!  Failed to parse internal coordinate statement.");
      }
    }
    else if ( ! strncmp("DELE",tok[0],4) ) {
      debug_msg("Recognized delete statement.");
      if ( ntok < 2 ) {
        print_msg(v,"ERROR!  Failed to parse delete statement.");
      } else {
        if ( ! strncmp("ATOM",tok[1],4) ) {
          if ( ntok < 3 ) {
            print_msg(v,"ERROR!  Failed to parse delete atom statement.");
          } else {
            s1 = parse_atom(tok[2],&i1,&j1);
            if ( topo_defs_atom(defs,0,1, s1,i1,j1,"DEL",0.0) ) {
              print_msg(v,"ERROR!  Failed to parse delete atom statement.");
            }
          }
        } else if ( ! strncmp("ACCE",tok[1],4) ) {
          ;
        } else if ( ! strncmp("DONO",tok[1],4) ) {
          ;
        } else if ( ! strncmp("BOND",tok[1],4) ||
		! strncmp("DOUB",tok[1],4) ||
		! strncmp("TRIP",tok[1],4) ) {
          if ( ntok < 4 || (ntok-2)%2 ) {
            print_msg(v,"ERROR!  Failed to parse delete bond statement.");
          } else {
            for ( itok = 2; itok < ntok; itok += 2 ) {
              s1 = parse_atom(tok[itok],&i1,&j1);
              s2 = parse_atom(tok[itok+1],&i2,&j2);
              if ( topo_defs_bond(defs,0,1,s1,i1,j1,s2,i2,j2) )
                print_msg(v,"ERROR!  Failed to parse delete bond statement.");
            }
          }
        } else if ( ! strncmp("ANGL",tok[1],4) ||
		! strncmp("THET",tok[1],4) ) {
          if ( ntok < 5 || (ntok-2)%3 ) {
            print_msg(v,"ERROR!  Failed to parse delete angle statement.");
          } else {
            for ( itok = 2; itok < ntok; itok += 3 ) {
              s1 = parse_atom(tok[itok],&i1,&j1);
              s2 = parse_atom(tok[itok+1],&i2,&j2);
              s3 = parse_atom(tok[itok+2],&i3,&j3);
              if ( topo_defs_angle(defs,0,1,s1,i1,j1,s2,i2,j2,s3,i3,j3) )
                print_msg(v,"ERROR!  Failed to parse delete angle statement.");
            }
          }
        } else if ( ! strncmp("DIHE",tok[1],4) ) {
          if ( ntok < 6 || (ntok-2)%4 ) {
            print_msg(v,"ERROR!  Failed to parse delete dihedral statement.");
          } else {
            for ( itok = 2; itok < ntok; itok += 4 ) {
              s1 = parse_atom(tok[itok],&i1,&j1);
              s2 = parse_atom(tok[itok+1],&i2,&j2);
              s3 = parse_atom(tok[itok+2],&i3,&j3);
              s4 = parse_atom(tok[itok+3],&i4,&j4);
              if (topo_defs_dihedral(defs,0,1,s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4))
                print_msg(v,"ERROR!  Failed to parse delete dihedral statement.");
            }
          }
        } else if ( ! strncmp("IMPH",tok[1],4) ||
		! strncmp("IMPR",tok[1],4) ) {
          if ( ntok < 6 || (ntok-2)%4 ) {
            print_msg(v,"ERROR!  Failed to parse delete improper statement.");
          } else {
            for ( itok = 2; itok < ntok; itok += 4 ) {
              s1 = parse_atom(tok[itok],&i1,&j1);
              s2 = parse_atom(tok[itok+1],&i2,&j2);
              s3 = parse_atom(tok[itok+2],&i3,&j3);
              s4 = parse_atom(tok[itok+3],&i4,&j4);
              if (topo_defs_improper(defs,0,1,s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4))
                print_msg(v,"ERROR!  Failed to parse delete improper statement.");
            }
          }
        } else if ( ! strncmp("BILD",tok[1],4) ||
		! strncmp("IC",tok[1],4) ) {
          if ( ntok < 6 ) {
            print_msg(v,"ERROR!  Failed to parse delete internal coordinate statement.");
          } else {
            s1 = parse_atom(tok[2],&i1,&j1);
            s2 = parse_atom(tok[3],&i2,&j2);
            s3 = tok[4] + ( tok[4][0] == '*' ? 1 : 0 );
            s3 = parse_atom(s3,&i3,&j3);
            s4 = parse_atom(tok[5],&i4,&j4);
            if ( topo_defs_conformation(defs,0,1,
		s1,i1,j1,s2,i2,j2,s3,i3,j3,s4,i4,j4,
		0,0,0,(tok[4][0]=='*'?1:0),0,0) )
              print_msg(v,"ERROR!  Failed to parse delete internal coordinate statement.");
          }
        } else {
          print_msg(v,"ERROR!  Failed to parse delete statement.");
        }
      }
    }
    else if ( ! strncmp("END",tok[0],4) ) {
      debug_msg("Recognized file end statement.");
      if ( topo_defs_end(defs) ) {
        print_msg(v,"ERROR!  Failed to parse file end statement.");
      }
    }
    else if ( ! strncmp("GROU",tok[0],4) ) {
      debug_msg("Recognized group statement.");
    }
    else if ( ! strncmp("PATC",tok[0],4) ) {
      debug_msg("Recognized patching statement.");
      if ( ntok < 3 || (ntok-1)%2 ) {
        print_msg(v,"ERROR!  Failed to parse patching statement.");
      } else {
        i1 = i2 = 0;
        for ( itok = 1; itok < ntok; itok += 2 ) {
          if ( ! strncmp("FIRS",tok[itok],4) ) {
              i1 = topo_defs_patching_first(defs,0,tok[itok+1]);
          } else if ( ! strncmp("LAST",tok[itok],4) ) {
              i2 = topo_defs_patching_last(defs,0,tok[itok+1]);
          } else {
            print_msg(v,"ERROR!  Failed to parse patching statement.");
          }
        }
        if ( i1 || i2 )
          print_msg(v,"ERROR!  Failed to parse patching statement.");
      }
    }
    else if ( ! strncmp("RESI",tok[0],4) ) {
      debug_msg("Recognized residue statement.");
      if ( ntok < 2 || topo_defs_residue(defs,tok[1],0) ) {
        print_msg(v,"ERROR!  Failed to parse residue statement.");
      }
    }
    else if ( ! strncmp("PRES",tok[0],4) ) {
      debug_msg("Recognized patch residue statement.");
      if ( ntok < 2 || topo_defs_residue(defs,tok[1],1) ) {
        print_msg(v,"ERROR!  Failed to parse patch residue statement.");
      }
    }
    else {
      sprintf(msgbuf,"ERROR!  FAILED TO RECOGNIZE %s",tok[0]);
      print_msg(v,msgbuf);
    }

  }

  topo_defs_end(defs);

  return 0;

}

