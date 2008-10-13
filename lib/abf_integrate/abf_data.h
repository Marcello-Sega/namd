/// \file integrate.h General headers for ABF_integrate

#include <iostream>
#include <vector>

/// Free energy gradients class
class ABFdata {

    protected:
	/// Sizes of (i-1) dimension blocks
	/// computed as Prod_(j<i) sizes[j]
	int *blocksizes;
	/// Minimum values of each variable
	double *mins;
	unsigned int vec_dim;

    public:
	int Nvars;
	/// Free energy gradients (vector field)
	double *gradients;
	/// Bin widths
	double *widths;

	unsigned int scalar_dim;
	unsigned int *histogram;

	/// History-dependent bias
	double *bias;

	/// Deviation between starting free energy gradient and that computed
	/// from MtD bias or histogram in standard MC
	double *deviation;

	void write_histogram ( const char *fileName );
	void write_bias ( const char *fileName );
	void write_deviation ( const char *fileName );

	/// Grid sizes
	int *sizes;

	/// Flag stating if each variable is periodic
	int *PBC;

	/// Constructor: reads from a file
	ABFdata ( const char *gradFileName );
       ~ABFdata ();

	/// \brief Returns an offset for scalar fields based on a n-index.
	/// multiply by Nvars to get an offset in a Nvars-vector field
	unsigned int offset (const int*);

	inline void wrap (int &pos, int i);
};

inline void ABFdata::wrap (int &pos, int i)
{
  if ( PBC[i] ) {
	if (pos == -1) {
	  pos = sizes[i] - 1;
	}
	if (pos == sizes[i]) {
	  pos = 0;
	}
  } else {
	// No PBC
	if (pos == -1) {
	  pos = 0;
	}
	if (pos == sizes[i]) {
	  pos = sizes[i] - 1;
	}
  }
  return;
}
