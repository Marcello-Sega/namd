
#include "vmdplugin.h"

extern int molfile_dcdplugin_init(void);
extern int molfile_dcdplugin_register(void *, vmdplugin_register_cb);
extern int molfile_dcdplugin_fini(void);
/*
extern int molfile_jsplugin_init(void);
extern int molfile_jsplugin_register(void *, vmdplugin_register_cb);
extern int molfile_jsplugin_fini(void);
*/

#include "molfile_plugin.h"

static int register_cb(void *v, vmdplugin_t *p) {
  *((vmdplugin_t **)v) = p;
}

#include <stdio.h>
#include <string.h>
#include <errno.h>

int main(int argc, char **argv) {

  molfile_timestep_t frame;
  molfile_plugin_t *plugin;
  char *output_root;
  char *filename;
  int num_replicas;
  int runs_per_frame;
  FILE **hist_in;
  FILE **hist_out;
  void **traj_in;
  void **traj_out;
  int natoms=MOLFILE_NUMATOMS_UNKNOWN;
  int i, i_run;

  molfile_dcdplugin_init();
  molfile_dcdplugin_register(&plugin, register_cb);

  if ( argc < 4 ) {
    fprintf(stderr, "args: <job_output_root> <num_replicas> <runs_per_frame>\n");
    exit(-1);
  }
  output_root = argv[1];
  num_replicas = atoi(argv[2]);
  runs_per_frame = atoi(argv[3]);

  filename = (char*) malloc(strlen(output_root)+100);
  hist_in = (FILE**) malloc(num_replicas*sizeof(FILE*));
  hist_out = (FILE**) malloc(num_replicas*sizeof(FILE*));
  traj_in = (void**) malloc(num_replicas*sizeof(FILE*));
  traj_out = (void**) malloc(num_replicas*sizeof(FILE*));

  for ( i=0; i<num_replicas; ++i ) {
    char *root_end;
    if ( strstr(output_root,"%s") ) {
      char istr[10];
      sprintf(istr,"%d",i);
      sprintf(filename,output_root,istr);
    } else {
      sprintf(filename,output_root,i);
    }
    root_end = filename + strlen(filename);

    sprintf(root_end,".%d.history",i);
    hist_in[i] = fopen(filename,"r");
    if ( ! hist_in[i] ) {
      fprintf(stderr, "error opening input file %s: %s\n",
					filename, strerror(errno));
      exit(-1);
    }
    sprintf(root_end,".%d.dcd",i);
    traj_in[i] = plugin->open_file_read(filename,"dcd",&natoms);
    if ( ! traj_in[i] ) {
      fprintf(stderr, "error opening input file %s: %s\n",
					filename, strerror(errno));
      exit(-1);
    }
    sprintf(root_end,".%d.sort.history",i);
    hist_out[i] = fopen(filename,"w");
    if ( ! hist_out[i] ) {
      fprintf(stderr, "error opening output file %s: %s\n",
					filename, strerror(errno));
      exit(-1);
    }
    sprintf(root_end,".%d.sort.dcd",i);
    traj_out[i] = plugin->open_file_write(filename,"dcd",natoms);
    if ( ! traj_out[i] ) {
      fprintf(stderr, "error opening output file %s: %s\n",
					filename, strerror(errno));
      exit(-1);
    }
  }

  frame.coords = (float*) malloc(3*natoms*sizeof(float));
  frame.velocities = (float*) NULL;

#define LINE_MAX 10000

  i_run = 0;
  for ( ; 1; ++i_run ) { /* loop until read fails */
    char line[LINE_MAX];
    for ( i=0; i<num_replicas; ++i ) {
      char *r;
      char sav;
      int f1,f2;
      int rc;
      int rep_id = -1;
      r = fgets(line, LINE_MAX, hist_in[i]);
      if ( ! r ) { break; }
      sscanf(line, "%*s %n%d%n", &f1, &rep_id, &f2);
      if ( rep_id < 0 || rep_id >= num_replicas ) {
        fprintf(stderr,"Invalid replica ID for replica %d at line %d: %s",
							i, i_run, line);
        exit(-1);
      }
      sav = line[f1];
      line[f1] = 0;
      fprintf(hist_out[rep_id],"%s%d%s",line,i,line+f2);
      line[f1] = sav;
      if ( i_run % runs_per_frame ) continue;
      rc = plugin->read_next_timestep(traj_in[i],natoms,&frame);
      if ( rc == MOLFILE_SUCCESS ) {
        plugin->write_timestep(traj_out[rep_id],&frame);
      } else {
        fprintf(stderr,"Unable to read frame for replica %d at line %d: %s",
							i, i_run, line);
        break;
      }
    }
    if ( i < num_replicas ) {
      printf("Processed %d runs.\n",i_run);
      if ( i ) fprintf(stderr,"Uneven input lengths for replica %d at line %d: %s",
							i, i_run, line);
      break;
    }
  }

  free(frame.coords);

  for ( i=0; i<num_replicas; ++i ) {
    fclose(hist_in[i]);
    plugin->close_file_read(traj_in[i]);
    fclose(hist_out[i]);
    plugin->close_file_write(traj_out[i]);
  }

  molfile_dcdplugin_fini();
  exit(0);
}

