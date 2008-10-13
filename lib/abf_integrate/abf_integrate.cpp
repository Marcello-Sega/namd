
#include "abf_data.h"
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <cmath>


// Note: this reads in *FE gradients* as output by colvars_abf (as opposed to average force)

char * parse_cl (int argc, char* argv[], unsigned int *nsteps, double *temp, int *meta, double *hill, double *hill_fact);

int main (int argc, char* argv[]) {
  char		*data_file;
  char		*out_file;
  unsigned int  step, nsteps, total;
  int		*pos, *dpos, *newpos;
  unsigned int  *histogram;
  const double  *grad, *newgrad;
  unsigned int  offset, newoffset;
  int		not_accepted;
  double	dA;
  double	temp;
  double	mbeta;
  int		meta;
  double	hill, hill_fact;
  int		phase;

  // Setting default values
  nsteps = 10000000;
  temp = 500;
  meta = 1;
  hill = 0.01;
  hill_fact = 0.5;

  if ( ! (data_file = parse_cl (argc, argv, &nsteps, &temp, &meta, &hill, &hill_fact)) ) {
    std::cerr << "Syntax: " << argv[0] << " <filename> [-n <nsteps>] [-t <temp>] [-m [0|1] (metadynamics)]"
	  " [-h <hill_height>] [-f <variable_hill_factor>]\n";
    exit(1);
  }

  if ( meta ) {
	std::cout << "\nUsing metadynamics-style sampling with hill height: " << hill << "\n";
	if ( hill_fact ) {
	      std::cout << "Varying hill height by factor " << hill_fact << "\n";
	}
  } else {
	std::cout << "\nUsing unbiased MC sampling\n";
  }
	std::cout << "Sampling " << nsteps << " steps at temperature " << temp << "\n\n";
 
  // Inverse temperature in (kcal/mol)-1 
  mbeta = -1/(0.001987*temp);

  ABFdata data (data_file);

  srand (time(NULL));

  pos = new int [data.Nvars];
  dpos = new int [data.Nvars];
  newpos = new int [data.Nvars];


  for (int i = 0; i<data.Nvars; i++) {
    pos[i] = rand() % data.sizes[i];
  }

  total = 0;
  for (step = 1; step <= nsteps; step++) {

    offset = data.offset(pos);
    data.histogram[offset]++;

    if (meta) {
      if ( hill_fact && ((10*step)%nsteps == 0) && ((phase=(10*step)/nsteps) >= 5) ) {
	hill *= hill_fact;
	std::cerr << "Changing hill height to " << hill << "\n";
      }
      data.bias[offset] += hill;
    }

    grad = data.gradients + offset * data.Nvars;

    not_accepted = 1;
    while (not_accepted) {
      dA = 0.0;
      total++;
      for (int i = 0; i<data.Nvars; i++) {
	dpos[i] = rand() % 3 - 1;
	// std::cout << dpos[i] << " ";

	newpos[i] = pos[i] + dpos[i];
	data.wrap (newpos[i], i);
	if (newpos[i] == pos[i]) dpos[i] = 0;

	if ( dpos[i] ) {
	  dA += grad[i] * dpos[i] * data.widths[i];
	  // usefulness of the interpolation below depends on
	  // where the grid points are for the histogram wrt to the gradients
	  // If done, it has to be done in all directions
	  // the line below is useless
	  //dA += 0.5 * (newgrad[i] + grad[i]) * dpos[i] * data.widths[i];
	}
      }

      if ( meta ) {
	newoffset = data.offset(newpos);
	dA += data.bias[newoffset] - data.bias[offset];
      }

      if (((float)rand())/RAND_MAX < exp(mbeta * dA)) {
	// Accept move
	for (int i = 0; i<data.Nvars; i++) {
	  pos[i] = newpos[i];
	  not_accepted = 0;
	}
      }
    }
  }
  std::cout << "Run " << total << " total iterations; acceptance ratio is "
	    << ((float)nsteps)/total << "\n";

  out_file = new char [strlen (data_file) + 8];
  if (meta) {
    sprintf (out_file, "%s.pmf", data_file);
    std::cout << "Writing bias data to file " << out_file << "\n";
    data.write_bias(out_file);
  }
  sprintf (out_file, "%s.histo", data_file);
  std::cout << "Writing histogram data to file " << out_file << "\n";
  data.write_histogram (out_file);

  // Computing deviation between gradients differentiated from pmf
  // and input data
  double *dev = data.deviation;
  grad = data.gradients;
  double estimate;

  for (int i=0; i<data.Nvars; i++)
    pos[i] = 0;
  for (offset = 0; offset < data.scalar_dim; offset++) {
    for (int i=data.Nvars-1; i>0; i--) {
      if (pos[i] == data.sizes[i]) {
	pos[i] = 0;
	pos[i-1]++;
      }
    }
    for (int i=0; i<data.Nvars; i++)
      newpos[i] = pos[i];
    
    for (int i=0; i<data.Nvars; i++) {
      // TODO: handle the borders of the grid properly
      newpos[i] = pos[i] - 1;
      data.wrap (newpos[i], i);
      newoffset = data.offset (newpos);
      if (meta) {
	estimate = 0.5 * (data.bias[newoffset] - data.bias[offset]) / data.widths[i];
      } else {
	// Note: computing a 300K estimate, arbitrarily
	if (data.histogram[offset] && data.histogram[newoffset])
	  estimate = 0.5 * 0.592 * log (double (data.histogram[offset]) / double (data.histogram[newoffset])) / data.widths[i];
	else
	  estimate = 0.0;
      }
      newpos[i] = pos[i] + 1;
      data.wrap (newpos[i], i);
      newoffset = data.offset (newpos);
      if (meta) {
	estimate += 0.5 * (data.bias[offset] - data.bias[newoffset]) / data.widths[i];
      } else {
	// Note: computing a 300K estimate, arbitrarily
	if (data.histogram[offset] && data.histogram[newoffset])
	  estimate += 0.5 * 0.592 * log (double (data.histogram[newoffset]) / double (data.histogram[offset])) / data.widths[i];
      }
      //dev[i] = estimate;
      dev[i] = grad[i] - estimate;
    }

    pos[data.Nvars-1]++;   // move on to next position
    dev += data.Nvars;
    grad += data.Nvars;
  }

  sprintf (out_file, "%s.dev", data_file);
  std::cout << "Writing FE gradient deviation to file " << out_file << "\n\n";
  data.write_deviation (out_file);

  delete[] pos;
  delete[] dpos;
  delete[] newpos;
  delete[] out_file;
  exit(0);
}


char * parse_cl (int argc, char* argv[], unsigned int *nsteps, double *temp, int *meta, double *hill, double *hill_fact)
{
    char    *filename = NULL;
    float   f_temp, f_hill;
    // "Syntax: " << argv[0] << " <filename> [-n <nsteps>] [-t <temp>] [-m [0|1] (metadynamics)] [-h <hill_height>]\n";
    if (argc < 2) {
      return NULL;
    }

    for (int i=2; i+1 < argc; i+=2) {
      if ( argv[i][0] != '-' ) {
	return NULL;
      }
      switch (argv[i][1]) {
	case 'n':
	  if ( sscanf (argv[i+1], "%u", nsteps) != 1) return NULL;
	  break;
	case 't':
	  if ( sscanf (argv[i+1], "%lf", temp) != 1) return NULL;
	  break;
	case 'm':
	  if ( sscanf (argv[i+1], "%u", meta) != 1) return NULL;
	  break;
	case 'h':
	  if ( sscanf (argv[i+1], "%lf", hill) != 1) return NULL;
	  break;
	case 'f':
	  if ( sscanf (argv[i+1], "%lf", hill_fact) != 1) return NULL;
	  break;
	default:
	  return NULL;
      }
    }
    return argv[1];
}
