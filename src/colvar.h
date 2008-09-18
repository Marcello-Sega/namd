#ifndef COLVAR_H
#define COLVAR_H

#include <iostream>
#include <iomanip>
#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"


/// \brief A collective variable (main class); needs at least one
/// colvar::cvc object to be defined and calculate a \link
/// colvarvalue \endlink object
///
/// This class parses the configuration, defines the behaviour and
/// stores the value (\link colvar::x \endlink) and all related data
/// of a collective variable.  How the value is calculated is defined
/// in \link colvar::cvc \endlink and its derived classes.  The
/// \link colvar \endlink object contains pointers to multiple \link
/// colvar::cvc \endlink derived objects, which can be combined
/// together into one collective variable.  This makes possible to
/// implement new collective variables at runtime based on the
/// existing ones.  Currently, this possibility is limited to a
/// polynomial, using the coefficients cvc::sup_coeff and the
/// exponents cvc::sup_np.  In case of non-scalar variables,
/// only exponents equal to 1 are accepted.
///
/// Please note that most of its members are \link colvarvalue
/// \endlink objects, i.e. they can handle different data types
/// together, and must all be set to the same type of colvar::x by
/// using the colvarvalue::type() member function before using them
/// together in assignments or other operations; this is usually done
/// automatically in the constructor.  If you add a new member of
/// \link colvarvalue \endlink type, you should also add its
/// initialization line in the \link colvar \endlink constructor.

class colvar : public colvarparse {

public:

  /// Name
  std::string name;

  /// Type of value
  colvarvalue::Type type() const;

  /// \brief Current value (previously obtained from calc() or read_traj())
  colvarvalue const & value() const;

  /// \brief Current velocity (previously obtained from calc() or read_traj())
  colvarvalue const & velocity() const;

  /// \brief Current system force (previously obtained from calc() or
  /// read_traj()).  Note: this is calculated using the atomic forces
  /// from the last simulation step.
  ///
  /// Total atom forces are read from the MD program, the total force
  /// acting on the collective variable is calculated summing those
  /// from all definition elements, the bias and walls forces are
  /// subtracted.
  colvarvalue const & system_force() const;

  /// \brief

  /// \brief Typical fluctuation amplitude for this collective
  /// variable (e.g. local width of a free energy basin)
  ///
  /// In metadynamics calculations, \link colvarbias_meta \endlink,
  /// this value is used to calculate the width of a gaussian.  In ABF
  /// calculations, \link colvarbias_abf \endlink, it is used to
  /// calculate the grid spacing in the direction of this collective
  /// variable.
  cvm::real width;

  /// \brief True if this \link colvar \endlink is a linear
  /// superposition of \link cvc \endlink elements
  bool b_linear;

  /// \brief True if all \link cvc \endlink objects are capable
  /// of calculating inverse gradients
  bool b_inverse_gradients;

  /// \brief True if all \link cvc \endlink objects are capable
  /// of calculating Jacobian forces
  bool b_Jacobian_force;

  /// \brief Options controlling the behaviour of the colvar during
  /// the simulation, which are set from outside the cvcs
  enum task {
    /// \brief Gradients are calculated and temporarily stored, so
    /// that external forces can be applied
    task_gradients,
    /// \brief Calculate the velocity with finite differences
    task_fdiff_velocity,
    /// \brief The system force is calculated, projecting the atomic
    /// forces on the inverse gradients
    task_system_force,
    /// \brief The variable has a harmonic restraint around a moving
    /// center with fictitious mass; bias forces will be applied to
    /// the center
    task_extended_lagrangian,
    /// \brief Compute analytically the "force" arising from the
    /// largest entropy component (for example, that from the angular
    /// states orthogonal to a distance vector)
    task_Jacobian_force,
    /// \brief Report the Jacobian force as part of the system force
    /// if disabled, apply a correction internally to cancel it
    task_report_Jacobian_force,
    /// \brief Output the value to the trajectory file (on by default)
    task_output_value,
    /// \brief Output the velocity to the trajectory file
    task_output_velocity,
    /// \brief Output the applied force to the trajectory file
    task_output_applied_force,
    /// \brief Output the system force to the trajectory file
    task_output_system_force,
    /// \brief Provide a discretization of the values of the colvar to
    /// be used by the biases or in analysis
    task_grid,
    /// \brief A restraining potential (|x-xb|^2) is applied to
    /// discourage the colvar from going below the lower boundary
    task_lower_wall,
    /// \brief A restraining potential (|x-xb|^2) is applied to
    /// discourage the colvar from going above the upper boundary
    task_upper_wall,
    /// \brief Number of possible tasks
    task_ntot
  };

  /// Tasks performed by this colvar
  bool tasks[task_ntot];

protected:


  /*
    extended:
    H = H_{0} + \sum_{i} 1/2*\lambda*(S_i(x(t))-s_i(t))^2 \\
    + \sum_{i} 1/2*m_i*(ds_i(t)/dt)^2 \\
    + \sum_{t'<t} W * exp (-1/2*\sum_{i} (s_i(t')-s_i(t))^2/(\delta{}s_i)^2) \\
    + \sum_{w} (\sum_{i}c_{w,i}s_i(t) - D_w)^M/(\Sigma_w)^M

    normal:
    H = H_{0} + \sum_{t'<t} W * exp (-1/2*\sum_{i} (S_i(x(t'))-S_i(x(t)))^2/(\delta{}S_i)^2) \\
    + \sum_{w} (\sum_{i}c_{w,i}S_i(t) - D_w)^M/(\Sigma_w)^M

    output:
    H = H_{0}   (only output S(x), no forces)

    Here:
    S(x(t)) = x
    s(t)    = xr
    DS = Ds = delta
  */


  /// Value of the colvar
  colvarvalue x;

  /// Cached reported value (x may be manipulated)
  colvarvalue x_reported;

  /// Finite-difference velocity
  colvarvalue v_fdiff;

  inline colvarvalue fdiff_velocity (colvarvalue const &xold,
                                     colvarvalue const &xnew)
  {
    // using the gradient of the square distance to calculate the
    // velocity (non-scalar variables automatically taken into
    // account)
    return ( ( (cvm::dt > 0.0) ? (1.0/cvm::dt) : 1.0 ) *
             0.5 * dist2_lgrad (xnew, xold) );
  }

  /// Cached reported velocity
  colvarvalue v_reported;

  // Options for task_extended_lagrangian
  /// Restraint center
  colvarvalue xr;
  /// Velocity of the restraint center
  colvarvalue vr;
  /// Mass of the restraint center
  cvm::real ext_mass;
  /// Restraint force constant
  cvm::real ext_force_k;
  /// \brief Harmonic restraint force
  colvarvalue fr;

  /// \brief Jacobian force, when task_Jacobian_force is enabled
  colvarvalue fj;

  /// Cached reported system force
  colvarvalue ft_reported;

public:


  /// \brief Bias force; reset_bias_force() should be called before
  /// the biases are updated
  colvarvalue fb;

  /// \brief Total \em applied force; fr (if task_extended_lagrangian
  /// is defined), fb (if biases are applied) and the walls' forces
  /// (if defined) contribute to it
  colvarvalue f;

  /// \brief Total force, as derived from the atomic trajectory;
  /// should equal the total system force plus \link f \endlink
  colvarvalue ft;


  /// Period, if defined
  cvm::real period;

  /// \brief Lower boundary value
  colvarvalue lower_boundary;
  /// \brief True if the lower boundary has been explicitly set
  bool        b_lower_boundary;
  /// \brief Force constant for the lower boundary potential (|x-xb|^2)
  cvm::real   lower_wall_k;

  /// \brief Upper boundary value
  colvarvalue upper_boundary;
  /// \brief True if the upper boundary has been explicitly set
  bool        b_upper_boundary;
  /// \brief Force constant for the upper boundary potential (|x-xb|^2)
  cvm::real   upper_wall_k;


  /// \brief Use the two boundaries and the width to report which bin
  /// the current value is in
  int current_bin_scalar() const;

  /// \brief Use the lower boundary and the width to report which bin
  /// the provided value is in
  int value_to_bin_scalar (colvarvalue const &val) const;

  /// \brief Same as the standard version, but uses another value
  /// instead of the lower boundary
  int value_to_bin_scalar (colvarvalue const &val,
                           colvarvalue const &offset) const;

  /// \brief Use the two boundaries and the width to report the
  /// central value corresponding to a bin index
  colvarvalue bin_to_value_scalar (int const &i_bin) const;

  /// \brief Same as the standard version, but uses another value
  /// instead of the lower boundary
  colvarvalue bin_to_value_scalar (int const &i_bin,
                                   colvarvalue const &offset) const;

  /// \brief Is the interval defined by the two boundaries periodic?
  bool periodic_boundaries() const;


  /// Constructor
  colvar (std::string const &conf);

  /// Enable the specified task
  void enable (colvar::task const &t);

  /// Disable the specified task
  void disable (colvar::task const &t);

  /// Destructor
  ~colvar();


  /// \brief Calculate the colvar value and all the other requested
  /// quantities
  void calc();

  /// \brief Calculate just the value, and store it in the argument
  void calc_value (colvarvalue &xn);

  /// Set the total biasing force to zero
  void reset_bias_force();

  /// Add to the total force from biases
  void add_bias_force (colvarvalue const &force);

  /// \brief Collect all forces on this colvar, integrate internal
  /// equations of motion of internal degrees of freedom; see also
  /// colvar::communicate_forces()
  void update();

  /// \brief Communicate forces (previously calculated in
  /// colvar::update()) to the external degrees of freedom
  void communicate_forces();


  /// \brief Use the internal metrics (as from \link cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  cvm::real dist2 (colvarvalue const &x1,
                   colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  colvarvalue dist2_lgrad (colvarvalue const &x1,
                           colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  colvarvalue dist2_rgrad (colvarvalue const &x1,
                           colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link cvc
  /// \endlink objects) to compare colvar values
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  cvm::real compare (colvarvalue const &x1,
                     colvarvalue const &x2) const;




  /// Read the analysis tasks
  void parse_analysis (std::string const &conf);
  /// Perform analysis tasks
  void analyse();


  /// Read the value from a collective variable trajectory file
  std::istream & read_traj (std::istream &is);
  /// Output formatted values to the trajectory file
  std::ostream & write_traj (std::ostream &os);
  /// Write a label to the trajectory file (comment line)
  std::ostream & write_traj_label (std::ostream &os);


  /// Read the collective variable from a restart file
  std::istream & read_restart (std::istream &is);
  /// Write the collective variable to a restart file
  std::ostream & write_restart (std::ostream &os);


protected:

  /// Previous value (to calculate velocities during analysis)
  colvarvalue            x_old;

  /// Time series of values and velocities used in auto-correlation
  /// functions (ACF) and running averages
  std::list< std::list<colvarvalue> > x_store, v_store;
  /// Time series of values and velocities used in auto-correlation
  /// functions (ACF) and running averages (standard deviations)
  std::list< std::list<colvarvalue> >::iterator x_store_p, v_store_p;

  /// Length of autocorrelation function (ACF)
  size_t                 acf_length;
  /// After how many steps the ACF starts
  size_t                 acf_offset;
  /// How many timesteps separate two ACF values
  size_t                 acf_stride;
  /// Number of frames for each ACF point
  size_t                 acf_nframes;
  /// Normalize the ACF to a maximum value of 1?
  bool                   acf_normalize;
  /// ACF values
  std::vector<cvm::real> acf;
  /// Name of the file to write the ACF
  std::string            acf_outfile;

  /// Type of autocorrelation function (ACF)
  enum acf_type_e {
    /// Unset type
    acf_notset,
    /// Velocity ACF, scalar product between v(0) and v(t)
    acf_vel,
    /// Coordinate ACF, scalar product between x(0) and x(t)
    acf_coor,
    /// \brief Coordinate ACF, second order Legendre polynomial
    /// between x(0) and x(t) (does not work with scalar numbers)
    acf_p2coor
  };

  /// Type of autocorrelation function (ACF)
  acf_type_e             acf_type;

  /// \brief Velocity ACF, scalar product between v(0) and v(t)
  void calc_vel_acf (std::list<colvarvalue> &v_store,
                     colvarvalue const      &v);

  /// \brief Coordinate ACF, scalar product between x(0) and x(t)
  /// (does not work with scalar numbers)
  void calc_coor_acf (std::list<colvarvalue> &x_store,
                      colvarvalue const      &x);

  /// \brief Coordinate ACF, second order Legendre polynomial between
  /// x(0) and x(t) (does not work with scalar numbers)
  void calc_p2coor_acf (std::list<colvarvalue> &x_store,
                        colvarvalue const      &x);

  /// Calculate the auto-correlation function (ACF)
  void calc_acf();
  /// Save the ACF to a file
  void write_acf (std::ostream &os);

  /// Length of running average series
  size_t         runave_length;
  /// Timesteps to skip between two values in the running average series
  size_t         runave_stride;
  /// Name of the file to write the running average
  std::ofstream  runave_os;
  /// Current value of the running average
  colvarvalue    runave;
  /// Current value of the square deviation from the running average
  cvm::real      runave_variance;

  /// Calculate the running average and its standard deviation
  void calc_runave();


public:


  // collective variable definition base class
  class cvc;

  // currently available collective variable definitions

  // scalar colvar components
  class distance;
  class distance_z;
  class distance_xy;
  class angle;
  class dihedral;
  class coordnum;
  class h_bond;
  class rmsd;
  class gyration;
  class alpha_dihedrals;
  class alpha_angles;

  // non-scalar components
  class distance_vec;
  class distance_dir;
  class orientation;
  class rotation;

protected:

  /// \brief Array of \link cvc \endlink objects
  std::vector<cvc *> cvcs;

public:

  inline size_t n_components () const {
    return cvcs.size();
  }
};


inline colvar * cvm::colvar_p (std::string const &name)
{
  for (std::vector<colvar *>::iterator cvi = cvm::colvars.begin();
       cvi != cvm::colvars.end();
       cvi++) {
    if ((*cvi)->name == name) {
      return (*cvi);
    }
  }
  return NULL;
}


inline colvarvalue::Type colvar::type() const
{
  return x.type();
}


inline colvarvalue const & colvar::value() const
{
  return x_reported;
}


inline colvarvalue const & colvar::velocity() const
{
  return v_reported;
}


inline colvarvalue const & colvar::system_force() const
{
  return ft_reported;
}


inline int colvar::current_bin_scalar() const
{
  return this->value_to_bin_scalar (this->value());
}

inline int colvar::value_to_bin_scalar (colvarvalue const &val) const
{
  return (int) ::floor ( (val.real_value - lower_boundary.real_value) / width );
}

inline int colvar::value_to_bin_scalar (colvarvalue const &val,
                                              colvarvalue const &offset) const
{
  return (int) ::floor ( (val.real_value - offset.real_value) / width );
}

inline colvarvalue colvar::bin_to_value_scalar (int const &i_bin) const
{
  return lower_boundary.real_value + width * (0.5 + i_bin);
}

inline colvarvalue colvar::bin_to_value_scalar (int const &i_bin,
                                                colvarvalue const &offset) const
{
  return offset.real_value + width * (0.5 + i_bin);
}


inline void colvar::add_bias_force (colvarvalue const &force)
{
  fb += force;
}

inline void colvar::reset_bias_force() {
  fb.reset();
}




/// \brief Grid of values of a function of several collective
/// variables \param T The data type
template <class T> class colvar_grid {

protected:

  /// Number of dimensions
  size_t nd;

  /// Number of points along each dimension
  std::vector<int> nx;

  /// Cumulative number of points along each dimension
  std::vector<int> nxc;

  /// \brief Multiplicity of each datum (allow the binning of
  /// non-scalar types)
  size_t mult;

  /// Total number of grid points
  size_t nt;

  /// Low-level array of values
  std::vector<T> data;

  /// Newly read data (used for count grids, when adding several grids read from disk)
  std::vector<size_t> new_data;

  /// Colvars collected in this grid
  std::vector<colvar *> cv;

  /// Get the low-level index corresponding to an index
  inline size_t address (std::vector<int> const &ix)
  {
    size_t addr = 0;
    for (size_t i = 0; i < nd; i++) {
      addr += ix[i]*nxc[i];
      if (cvm::debug()) {
        if (ix[i] >= nx[i])
          cvm::fatal_error ("Error: exceeding bounds in colvar_grid.\n");
      }
    }
    return addr;
  }

public:

  /// True if this is a count grid related to another grid of data
  bool has_parent_data;

  /// Return the number of colvars
  inline size_t number_of_colvars()
  {
    return nd;
  }

  /// Return the number of points in the i-th direction, if provided, or
  /// the total number
  inline size_t number_of_points (int const icv = -1)
  {
    if (icv < 0) {
      return nt;
    } else {
      return nx[icv];
    }
  }

  /// Return the multiplicity of the type used
  inline size_t multiplicity()
  {
    return mult;
  }

  /// \brief Allocate data (allow initialization also after construction)
  void create (std::vector<int> const &nx_i,
               T const &t = T(),
               size_t const &mult_i = 1)
  {
    mult = mult_i;
    nd = nx_i.size();
    nxc.resize (nd);
    nx = nx_i;

    nt = mult;
    for (int i = nd-1; i >= 0; i--) {
      if (nx_i[i] <= 0)
        cvm::fatal_error ("Error: providing an invalid number of points, "+
                          cvm::to_str (nx_i[i])+".\n");
      nxc[i] = nt;
      nt *= nx[i];
    }

    data.reserve (nt);
    data.assign (nt, t);
  }

  /// Default constructor
  colvar_grid()
  {
    nd = nt = 0;
  }

  /// Destructor
  virtual ~colvar_grid()
  {}

  /// \brief Constructor from explicit grid sizes \param nx_i Number
  /// of grid points along each dimension \param t Initial value for
  /// the function at each point (optional) \param mult_i Multiplicity
  /// of each value
  colvar_grid (std::vector<int> const &nx_i,
               T const &t = T(),
               size_t const &mult_i = 1)
  {
    this->create (nx_i, t, mult_i);
  }

  /// \brief Constructor from a vector of colvars
  colvar_grid (std::vector<colvar *> const &colvars,
               T const &t = T(),
               size_t const &mult_i = 1,
               size_t const &bins_scale = 1)
    : cv (colvars)
  {
    std::vector<int> nx_i;

    if (cvm::debug())
      cvm::log ("Allocating a grid for "+cvm::to_str (colvars.size())+
                " collective variables.\n");

    for (size_t i =  0; i < cv.size(); i++) {

      if (cv[i]->type() != colvarvalue::type_scalar) {
        cvm::fatal_error ("Colvar grids can only be automatically "
                          "constructed for scalar variables.\n");
      }

      if (cv[i]->width <= 0.0) {
        cvm::fatal_error ("Tried to initialize a grid on a"
                          "variable with negative or zero width.\n");
      }

      if (!cv[i]->b_lower_boundary || !cv[i]->b_upper_boundary) {
        cvm::fatal_error ("Tried to initialize a count grid on "
                          "variable with undefined boundaries.\n");
      }

      double integer;
      cvm::real const fract = ::modf (( cv[i]->upper_boundary -
                                        cv[i]->lower_boundary ) / cv[i]->width,
                                      &integer);

      int n = ((int) integer) * bins_scale;
      if (fract > 1.0E-10) {
        cvm::log ("Warning: colvar interval ("+
                  cvm::to_str (cv[i]->lower_boundary, cvm::cv_width, cvm::cv_prec)+" - "+
                  cvm::to_str (cv[i]->upper_boundary, cvm::cv_width, cvm::cv_prec)+
                  ") is not commensurate to its width ("+
                  cvm::to_str (cv[i]->width, cvm::cv_width, cvm::cv_prec)+").\n");
        n += bins_scale;
      }

      if (cvm::debug())
        cvm::log ("Number of points is "+cvm::to_str (n)+
                  " for the colvar no. "+cvm::to_str (i+1)+".\n");

      nx_i.push_back (n);
    }

    create (nx_i, t, mult_i);
  }

  /// Set the value at the point with index ix
  inline void set_value (std::vector<int> const &ix,
                         T const &t,
                         size_t const &imult = 0)
  {
    data[this->address (ix)+imult] = t;
  }

  /// \brief Get the binned value indexed by ix, or the first of them
  /// if the multiplicity is larger than 1
  inline T const & value (std::vector<int> const &ix,
                          size_t const &imult = 0)
  {
    return data[this->address (ix) + imult];
  }



  /// \brief Add a constant to all elements (fast loop)
  inline void add_constant (T const &t)
  {
    for (size_t i = 0; i < nt; i++) 
      data[i] += t;
  }

  /// \brief Multiply all elements by a scalar constant (fast loop)
  inline void multiply_constant (cvm::real const &a)
  {
    for (size_t i = 0; i < nt; i++) 
      data[i] *= a;
  }

  /// \brief Add data from another grid of the same type \param grid
  /// Reference to the other grid \param grid_boundaries Lower
  /// boundaries of the other grid (the width is assumed to be the
  /// same)
  void add_grid (colvar_grid<T> &grid,
                 std::vector<colvarvalue> const &grid_boundaries)
  {
    if (grid.multiplicity() != this->multiplicity())
      cvm::fatal_error ("Error: trying to merge two grids with values of "
                        "different multiplicities.\n");

    std::vector<colvarvalue> const &gb = grid_boundaries;
    std::vector<int> ixn = this->new_index();
    for (std::vector<int> ix = grid.new_index();
         grid.index_ok (ix); grid.incr (ix)) {
      for (size_t i = 0; i < nd; i++) {
        ixn[i] =
          cv[i]->value_to_bin_scalar (cv[i]->bin_to_value_scalar (ix[i],
                                                                  gb[i]));
      }
      if (! this->index_ok (ixn))
        continue;
      for (size_t im = 0; im < mult; im++) {
        this->set_value (ixn, (this->value (ixn, im) +
                               grid.value (ix, im)), im);
      }
    }
  }

  /// \brief Return the value suitable for output purposes (so that it
  /// may be rescaled or manipulated without changing it permanently)
  virtual inline T value_output (std::vector<int> const &ix,
                                 size_t const &imult = 0)
  {
    return value (ix, imult);
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (it may have been rescaled or
  /// manipulated)
  virtual inline void value_input (std::vector<int> const &ix,
                                   T const &t,
                                   size_t const &imult = 0,
                                   bool add = false)
  {
    if ( add )
      data[address (ix) + imult] += t;
    else
      data[address (ix) + imult] = t;
  }

  //   /// Get the pointer to the binned value indexed by ix
  //   inline T const *value_p (std::vector<int> const &ix)
  //   {
  //     return &(data[address (ix)]);
  //   }

  /// \brief Get the index corresponding to the "first" bin, to be
  /// used as the initial value for an index in looping
  inline std::vector<int> const new_index() const
  {
    return std::vector<int> (nd, 0);
  }

  /// \brief Check that the index is within range in each of the
  /// dimensions
  inline bool index_ok (std::vector<int> const &ix) const
  {
    for (size_t i = 0; i < nd; i++) {
      if ( (ix[i] < 0) || (ix[i] >= int (nx[i])) )
        return false;
    }
    return true;
  }

  /// \brief Increment the index, in a way that will make it loop over
  /// the whole nd-dimensional array
  inline void incr (std::vector<int> &ix) const
  {
    for (int i = ix.size()-1; i >= 0; i--) {

      ix[i]++;

      if (ix[i] >= nx[i]) {

        if (i > 0) {
          ix[i] = 0;
          continue;
        } else {
          // this is the last iteration, a non-valid index is being
          // set for the outer index, which will be caught by
          // index_ok()
          ix[0] = nx[0];
          return;
        }
      } else {
        return;
      }
    }
  }

  /// \brief Write the grid data without labels, as they are
  /// represented in memory
  /// \param buf_size Number of values per line
  void write_raw (std::ostream &os, size_t const &buf_size = 3)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    std::vector<int> ix = new_index();
    size_t count = 0;
    for ( ; index_ok (ix); incr (ix)) {
      for (size_t imult = 0; imult < mult; imult++) {
        os << " "
           << std::setw (w) << std::setprecision (p)
           << value_output (ix, imult);
        if (((++count) % buf_size) == 0)
          os << "\n";
      }
    }
    // write a final newline only if buffer is not empty
    if ((count % buf_size) != 0)
      os << "\n";
  }

  /// \brief Read data written by colvar_grid::write_raw()
  void read_raw (std::istream &is)
  {
    std::vector<int> ix = new_index();
    size_t count = 0;
    for ( ; index_ok (ix); incr (ix)) {
      for (size_t imult = 0; imult < mult; imult++) {
        T new_value;
        if (is >> new_value) {
          value_input (ix, new_value, imult);
          count++;
        }
      }
    }

    if (count != number_of_points())
      cvm::fatal_error ("Error: number of points read from file ("+
                        cvm::to_str (count)+") is wrong (should be "+
                        cvm::to_str (number_of_points())+").\n");
  }

  /// \brief Write the grid in a format which is both human readable
  /// and suitable for visualization e.g. with gnuplot
  void write_multicol (std::ostream &os)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    // Data in the header: nColvars, then for each
    // xiMin, dXi, nPoints, periodic

    os << std::setw (2) << "# " << nd << "\n";
    for (size_t i = 0; i < nd; i++) {
      os << "# "
         << std::setw (10) << cv[i]->lower_boundary
         << std::setw (10) << cv[i]->width
         << std::setw (10) << nx[i] << "  "
         << cv[i]->periodic_boundaries() << "\n";
    }

    for (std::vector<int> ix = new_index(); index_ok (ix); incr (ix) ) {

      if (ix.back() == 0) {
        // add a new line to delimit the record
        os << "\n";
      }

      for (size_t i = 0; i < nd; i++) {
        os << " "
           << std::setw (w) << std::setprecision (p)
           << cv[i]->bin_to_value_scalar (ix[i]);
      }
      os << " ";
      for (size_t imult = 0; imult < mult; imult++) {
        os << " "
           << std::setw (w) << std::setprecision (p)
           << value_output (ix, imult);
      }
      os << "\n";
    }
  }

  /// \brief Read a grid written by colvar_grid::write_multicol()
  /// Adding data if add is true, replacing if false
  /// If adding data to a count grid, keep newly read data separately
  /// (taken care of by value_input)
  void read_multicol (std::istream &is, bool add = false)
  {
    // Data in the header: nColvars, then for each
    // xiMin, dXi, nPoints, periodic

    std::string   hash;
    cvm::real             lower, width, x;
    size_t        n, periodic;
    bool          remap;
    std::vector<T>        new_value;
    std::vector<int>      nx_read;
    std::vector<int>      bin;

    if ( cv.size() != nd ) {
      cvm::fatal_error ("Cannot grid file: missing reference to colvars.");
    }

    if ( !(is >> hash) || (hash != "#") ) {
      cvm::fatal_error ("Error reading grid at position "+
                        cvm::to_str (is.tellg())+" in stream (read \"" + hash + "\")\n");
    }

    is >> n;
    if ( n != nd ) {
      cvm::fatal_error ("Error reading grid: wrong number of collective variables.\n");
    }

    nx_read.resize (n);
    bin.resize (n);
    new_value.resize (mult);

    if ( this->has_parent_data && add ) {
      new_data.resize (data.size());
    }

    remap = false;
    for (size_t i = 0; i < nd; i++ ) {
      if ( !(is >> hash) || (hash != "#") ) {
        cvm::fatal_error ("Error reading grid at position "+
                          cvm::to_str (is.tellg())+" in stream (read \"" + hash + "\")\n");
      }

      is >> lower >> width >> nx_read[i] >> periodic;

      if ( (::fabs (lower - cv[i]->lower_boundary.real_value) > 1.0e-10) ||
           (::fabs (width - cv[i]->width ) > 1.0e-10) ||
           (nx_read[i] != nx[i]) ) {
        cvm::log ("Warning: reading from different grid definition (colvar "
                  + cvm::to_str (i+1) + "); remapping data on new grid.\n");
        remap = true;
      }
    }

    if ( remap ) {
      // re-grid data
      while (is.good()) {
	bool end_of_file = false;

        for (size_t i = 0; i < nd; i++ ) {
          if ( !(is >> x) ) end_of_file = true;
          bin[i] = cv[i]->value_to_bin_scalar (x);
        }
	if (end_of_file) break;	

        for (size_t imult = 0; imult < mult; imult++) {
          is >> new_value[imult];
        }

        if ( index_ok(bin) ) {
          for (size_t imult = 0; imult < mult; imult++) {
            value_input (bin, new_value[imult], imult, add);
          }
        }
      }
    } else {
      // do not re-grid the data but assume the same grid is used
      for (std::vector<int> ix = new_index(); index_ok (ix); incr (ix) ) {
        for (size_t i = 0; i < nd; i++ ) {
          is >> x;
        }
        for (size_t imult = 0; imult < mult; imult++) {
          is >> new_value[imult];
          value_input (ix, new_value[imult], imult, add);
        }
      }
    }
    return;
  }

};



/// \brief Colvar_grid derived class to hold counters in discrete
/// n-dim colvar space
class colvar_grid_count : public colvar_grid<size_t>
{
public:

  /// Default constructor
  colvar_grid_count();

  /// Destructor
  virtual inline ~colvar_grid_count()
  {}

  /// Constructor
  colvar_grid_count (std::vector<int> const &nx_i,
                     size_t const           &def_count = 0);

  /// Constructor from a vector of colvars
  colvar_grid_count (std::vector<colvar *>  &colvars,
                     size_t const           &def_count = 0,
                     size_t const           &bins_scale = 1);

  /// Increment the counter at given position
  inline void incr_count (std::vector<int> const &ix)
  {
    ++(data[this->address (ix)]);
  }

  /// \brief Get the binned count indexed by ix from the newly read data 
  inline size_t const & new_count (std::vector<int> const &ix,
                                   size_t const &imult = 0)
  {
    return new_data[address (ix) + imult];
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (it may have been rescaled or
  /// manipulated)
  virtual inline void value_input (std::vector<int> const &ix,
                                   size_t const &t,
                                   size_t const &imult = 0,
                                   bool add = false)
  {
    if ( add ) {
      data[address (ix)] += t;
      // Save newly read data for inputting parent grid
      if ( this->has_parent_data )
	new_data[address (ix)] = t;
    } else {
      data[address (ix)] = t;
    }
  }
};


/// Class for accumulating a scalar function on a grid
class colvar_grid_scalar : public colvar_grid<cvm::real>
{
public:

  /// \brief Provide the associated sample count by which each binned value
  /// should be divided
  colvar_grid_count *samples;

  /// Default constructor
  colvar_grid_scalar();

  /// Destructor
  virtual inline ~colvar_grid_scalar()
  {}

  /// Constructor from specific sizes arrays
  colvar_grid_scalar (std::vector<int> const &nx_i);

  /// Constructor from a vector of colvars
  colvar_grid_scalar (std::vector<colvar *> &colvars,
                      size_t const          &bins_scale = 1);

  /// Accumulate the value
  inline void acc_value (std::vector<int> const &ix,
                         cvm::real const &new_value,
                         size_t const &imult = 0)
  {
    // only legal value of imult here is 0
    data[address (ix)] += new_value;
    if (samples)
      samples->incr_count (ix);
  }

  /// \brief Return the value of the function at ix divided by its
  /// number of samples (if the count grid is defined)
  virtual cvm::real value_output (std::vector<int> const &ix,
                                  size_t const &imult = 0)
  {
    if (imult > 0)
      cvm::fatal_error ("Error: trying to access a component "
                        "larger than 1 in a scalar data grid.\n");
    if (samples)
      return (samples->value (ix) > 0) ?
        (data[address (ix)] / cvm::real (samples->value (ix))) :
        0.0;
    else
      return data[address (ix)];
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (it may have been rescaled or
  /// manipulated)
  virtual void value_input (std::vector<int> const &ix,
                            cvm::real const &new_value,
                            size_t const &imult = 0,
                            bool add = false)
  {
    if (imult > 0)
      cvm::fatal_error ("Error: trying to access a component "
                        "larger than 1 in a scalar data grid.\n");
    if (add) {
      if (samples)
        data[address (ix)] += new_value * samples->new_count (ix);
      else
        data[address (ix)] += new_value;
    } else {
      if (samples)
        data[address (ix)] = new_value * samples->value (ix);
      else
        data[address (ix)] = new_value;
    }
  }

  /// \brief Return the highest value
  inline cvm::real maximum_value()
  {
    cvm::real max = -1.0E+99;
    for (size_t i = 0; i < nt; i++) {
      if (data[i] > max) max = data[i];
    }
    return max;
  }

  /// \brief Return the lowest value
  inline cvm::real minimum_value()
  {
    cvm::real min = 1.0E+99;
    for (size_t i = 0; i < nt; i++) {
      if (data[i] > min) min = data[i];
    }
    return min;
  }

};



/// Class for accumulating the gradient of a scalar function on a grid
class colvar_grid_gradient : public colvar_grid<cvm::real>
{
public:

  /// \brief Provide the sample count by which each binned value
  /// should be divided
  colvar_grid_count *samples;

  /// Default constructor
  colvar_grid_gradient();

  /// Destructor
  virtual inline ~colvar_grid_gradient()
  {}

  /// Constructor from specific sizes arrays
  colvar_grid_gradient (std::vector<int> const &nx_i);

  /// Constructor from a vector of colvars
  colvar_grid_gradient (std::vector<colvar *>  &colvars,
                        size_t const           &bins_scale = 1);

  /// \brief Accumulate the gradient
  inline void acc_grad (std::vector<int> const &ix, cvm::real const *grads) {
    for (size_t imult = 0; imult < mult; imult++) {
      data[address (ix) + imult] += grads[imult];
    }
    if (samples)
      samples->incr_count (ix);
  }

  /// \brief Accumulate the gradient based on the force (i.e. sums the
  /// opposite of the force)
  inline void acc_force (std::vector<int> const &ix, cvm::real const *forces) {
    for (size_t imult = 0; imult < mult; imult++) {
      data[address (ix) + imult] -= forces[imult];
    }
    if (samples)
      samples->incr_count (ix);
  }

  /// \brief Return the value of the function at ix divided by its
  /// number of samples (if the count grid is defined)
  virtual inline cvm::real value_output (std::vector<int> const &ix,
                                         size_t const &imult = 0)
  {
    if (samples)
      return (samples->value (ix) > 0) ?
        (data[address (ix) + imult] / cvm::real (samples->value (ix))) :
        0.0;
    else
      return data[address (ix) + imult];
  }

  /// \brief Get the value from a formatted output and transform it
  /// into the internal representation (it may have been rescaled or
  /// manipulated)
  virtual inline void value_input (std::vector<int> const &ix,
                                   cvm::real const &new_value,
                                   size_t const &imult = 0,
                                   bool add = false)
  {
    if (add) {
      if (samples)
        data[address (ix) + imult] += new_value * samples->new_count (ix);
      else
        data[address (ix) + imult] += new_value;
    } else {
      if (samples)
        data[address (ix) + imult] = new_value * samples->value (ix);
      else
        data[address (ix) + imult] = new_value;
    }
  }

  /// \brief If the grid is 1-dimensional, integrate it and write the
  /// integral to a file
  void write_1D_integral (std::ostream &os);
};



#endif


// Emacs
// Local Variables:
// mode: C++
// End:
