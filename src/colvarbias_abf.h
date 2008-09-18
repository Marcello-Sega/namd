#ifndef COLVARBIAS_ABF_H
#define COLVARBIAS_ABF_H

#include <vector>
#include <list>
#include <sstream>
#include <iomanip>
//#include <cmath>

#include "colvarbias.h"

typedef cvm::real* gradient_t;


/// ABF bias
class colvarbias_abf : public colvarbias {

public:

  colvarbias_abf (std::string const &conf, char const *key);
  ~colvarbias_abf ();

  void update ();

private:

  /// Filenames for human-readable gradient/sample count output (updated at restartFreq and end-of-run)
  std::string	gradients_out_name;
  std::string	samples_out_name;
  std::string	pmf_out_name;		// only used for 1D calculations

  /// Base filename(s) for reading previous gradient data (replaces data from restart file)
  std::vector<std::string> input_prefix;

  bool		apply_bias;
  bool		hide_Jacobian;
  size_t	full_samples;
  size_t	min_samples;
  /// frequency for updating output files (default: same as restartFreq?)
  int		output_freq;

  // below: stuff for force histograms
  //float	fMax; 
  //float	df; 

  // Internal data and methods

  std::vector<int>  bin, prev_bin;
  gradient_t	    prev_force;

  /// n-dim grid of free energy gradients
  colvar_grid_gradient  *gradients;
  /// n-dim grid of number of samples
  colvar_grid_count     *samples;

  /// Write human-readable FE gradients and sample count
  void		  write_gradients_samples ();
  std::ofstream	  gradients_os;  /// Stream for writing FE gradients to disk
  std::ofstream	  samples_os;    /// Stream for writing sampling to disk
  std::ofstream	  pmf_os;		 /// Stream for writing 1D pmf to disk

  /// Read human-readable FE gradients and sample count (if not using restart)
  void		  read_gradients_samples ();

  std::istream& read_restart  (std::istream&);
  std::ostream& write_restart (std::ostream&);
};


/// Histogram "bias" (does as the name says)
class colvarbias_histogram : public colvarbias {

public:

  colvarbias_histogram (std::string const &conf, char const *key);
  ~colvarbias_histogram ();

  void update ();

private:

  /// n-dim histogram
  colvar_grid_count    *grid;
  std::vector<int>  bin;
  std::string	  out_name;

  int		  output_freq;  
  void		  write_grid ();
  std::ofstream	  grid_os;  /// Stream for writing grid to disk

  std::istream& read_restart  (std::istream&);
  std::ostream& write_restart (std::ostream&);
};

#endif
