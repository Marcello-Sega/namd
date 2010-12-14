#ifndef COLVARBIAS_META_H
#define COLVARBIAS_META_H

#include <vector>
#include <list>
#include <sstream>
#include <fstream>

#include "colvarbias.h"
#include "colvargrid.h"

/// Metadynamics bias (implementation of \link colvarbias \endlink)
class colvarbias_meta : public colvarbias {

public:

  /// Communication between different replicas
  enum Communication {
    /// One replica (default)
    single_replica,
    /// Hills added concurrently by several replicas
    multiple_replicas
  };

  /// Communication between different replicas
  Communication comm;

  /// Constructor
  colvarbias_meta (std::string const &conf, char const *key);

  /// Destructor
  virtual ~colvarbias_meta();
  
  virtual cvm::real update();

  /// Perform analysis
  virtual void analyse();

  virtual std::istream & read_restart (std::istream &is);

  virtual std::ostream & write_restart (std::ostream &os);

  class hill;
  typedef std::list<hill>::iterator hill_iter;

protected:

  /// Parse analysis tasks and options
//   virtual void parse_analysis (std::string const &conf);

  /// \brief List of hills used on this bias (total); if a grid is
  /// employed, these don't need to be updated at every time step
  std::list<hill> hills;

  /// \brief Iterator to the first of the "newest" hills (when using
  /// grids, those who haven't been mapped yet)
  hill_iter new_hills_begin;

  /// \brief List of hills used on this bias that are on the boundary
  /// edges; these are updated regardless of whether hills are used
  std::list<hill> hills_off_grid;

  /// \brief Same as new_hills_begin, but for the off-grid ones
  hill_iter new_hills_off_grid_begin;

  /// Regenerate the hills_off_grid list
  void recount_hills_off_grid (hill_iter h_first, hill_iter h_last,
                               colvar_grid_scalar *ge);

  /// Read a hill from a file
  std::istream & read_hill (std::istream &is);

  /// \brief Add a new hill; if a .hills trajectory is written,
  /// write it there; if there is more than one replica, communicate
  /// it to the others
  virtual std::list<hill>::const_iterator create_hill (hill const &h);

  /// \brief Remove a previously saved hill (returns an iterator for
  /// the next hill in the list)
  virtual std::list<hill>::const_iterator delete_hill (std::list<hill>::iterator &h);

  /// \brief Calculate the values of the hills, incrementing
  /// bias_energy
  virtual void calc_hills (hill_iter  h_first,
                           hill_iter  h_last,
                           cvm::real &energy,
                           std::vector<colvarvalue> const &values = std::vector<colvarvalue> (0));

  /// \brief Calculate the forces acting on the i-th colvar,
  /// incrementing colvar_forces[i]; must be called after calc_hills
  /// each time the values of the colvars are changed
  virtual void calc_hills_force (size_t const &i,
                                 hill_iter h_first,
                                 hill_iter h_last,
                                 std::vector<colvarvalue> &forces,
                                 std::vector<colvarvalue> const &values = std::vector<colvarvalue> (0));

  /// Height of new hills
  cvm::real  hill_weight;

  /// \brief Bin the hills on grids of energy and forces, and use them
  /// to force the colvars (as opposed to deriving the hills analytically)
  bool       use_grids;

  /// \brief Rebin the hills upon restarting
  bool       rebin_grids;

  /// \brief Should the grids be expanded if necessary?
  bool       expand_grids;

  /// \brief How often the hills should be projected onto the grids
  size_t     grids_freq;

  /// \brief Whether to keep the hills in the restart file (e.g. to do
  /// meaningful accurate rebinning afterwards)
  bool       keep_hills;

  /// \brief Dump the free energy surface (.pmf file) every restartFrequency
  bool       dump_fes;

  /// \brief Keep the free energy surface files at different
  /// iterations, appending a step number to each
  bool       dump_fes_save;

  /// Hill energy, cached on a grid
  colvar_grid_scalar    *hills_energy;

  /// Hill forces, cached on a grid
  colvar_grid_gradient  *hills_energy_gradients;

  /// Project the selected hills onto grids
  void project_hills (hill_iter h_first, hill_iter h_last,
                      colvar_grid_scalar *ge, colvar_grid_gradient *gf,
                      cvm::real const scale_factor = 1.0);



  /// \brief width of a hill
  ///
  /// The local width of each collective variable, multiplied by this
  /// number, provides the hill width along that direction
  cvm::real  hill_width;

  /// \brief Number of simulation steps between two hills
  size_t     new_hill_freq;

  /// Write the hill logfile
  bool           b_hills_traj;
  /// Logfile of hill management (creation and deletion)
  std::ofstream  hills_traj_os;


  /// Identifier for this replica
  std::string    replica;

  /// \brief Add this replica to the registry (called only when \link
  /// comm \endlink != \link single_replica \endlink)
  virtual void register_replica_file (std::string const &new_file);

  /// \brief Read the names all replica output files which have been
  /// added, skipping the one which is currently being written to
  /// (called only when \link colvarbias_meta::comm \endlink != \link
  /// colvarbias_meta::single_replica \endlink)
  virtual void update_replica_files_registry();

  /// \brief Read new data from replicas' files (called only when \link
  /// colvarbias_meta::comm \endlink != \link
  /// colvarbias_meta::single_replica \endlink)
  virtual void read_replica_files();

  /// \brief Frequency at which output files from other replicas are
  /// checked
  size_t                 replica_update_freq;

  /// List of hill files from all the replicas
  std::string            replica_files_registry;

  /// Replicas file names
  std::list<std::string> replica_files;

  /// Positions in replica files (files are reopened at these positions)
  std::list<size_t>      replica_files_pos;

  /// \brief File to contain the hills created by this replica in this run
  std::ofstream          replica_out_file;

  /// \brief File to contain the hills created by this replica in this run (name)
  std::string            replica_out_file_name;


  // Analysis

  /// \brief Select only hills after this step
  size_t                 free_energy_begin;

  /// \brief Select only hills before this step
  size_t                 free_energy_end;

  /// \brief Make the free energy surface be larger than or equal to zero
  bool                   shift_fes;

  /// \brief Subtract this value from the free energy (after 
  cvm::real              free_energy_offset;

  /// \brief Free energy surface output file
  std::string            free_energy_file;

  /// \brief Free energy gradients output file
  std::string            free_energy_gradients_file;

  /// \brief Boltzmann weights output file
  std::string            boltzmann_weights_file;

  /// \brief Boltzmann counts output file.
  /// 
  /// These are discretized Boltzmann weights, useful e.g. to
  /// initialize the samples of an ABF calculation.
  /// boltzmann_weights_scale must be much larger than 1 for a proper
  /// discretization
  std::string            boltzmann_counts_file;

  // weight = scale * exp (-(fe-offset)/(kB*temp))

  /// \brief Multiply the exponential weights by this constant
  cvm::real              boltzmann_weights_scale;

  /// \brief Temperature (in K) for Boltzmann weights
  cvm::real              boltzmann_weights_temp;

  /// Free energy values
  colvar_grid_scalar    *free_energy;

  /// Free energy gradients
  colvar_grid_gradient  *free_energy_gradients;

  /// Boltzmann populations (real version)
  colvar_grid_scalar    *boltzmann_weights;

  /// Boltzmann populations (integer version)
  colvar_grid_count     *boltzmann_counts;

};




/// \brief A hill for the metadynamics bias
class colvarbias_meta::hill {

protected:

  /// Value of the hill function (ranges between 0 and 1)
  cvm::real hill_value;

  /// Scale factor, which could be modified at runtime (default: 1)
  cvm::real sW;

  /// Maximum height in energy of the hill
  cvm::real W;

  /// Center of the hill in the collective variable space
  std::vector<colvarvalue>  centers;

  /// Widths of the hill in the collective variable space
  std::vector<cvm::real>    widths;

public:

  friend class colvarbias_meta;

  /// Time step at which this hill was added
  size_t      it;

  /// Identity of the replica who added this hill (only in multiple replica simulations)
  std::string replica;

  /// \brief Runtime constructor: data are read directly from
  /// collective variables \param weight Weight of the hill \param
  /// cv Pointer to the array of collective variables involved \param
  /// replica (optional) Identity of the replica which creates the
  /// hill
  inline hill (cvm::real             const &W_in,
               std::vector<colvar *>       &cv,
               cvm::real             const &hill_width,
               std::string           const &replica_in = "")
    : sW (1.0),
      W (W_in),
      centers (cv.size()),
      widths (cv.size()),
      it (cvm::it),
      replica (replica_in)
  {
    for (size_t i = 0; i < cv.size(); i++) {
      centers[i].type (cv[i]->type());
      centers[i] = cv[i]->value();
      widths[i] = cv[i]->width * hill_width;
    }
    if (cvm::debug()) 
      cvm::log ("New hill, applied to "+cvm::to_str (cv.size())+
                " collective variables, with reference values "+
                cvm::to_str (centers)+", widths "+
                cvm::to_str (widths)+" and weight "+
                cvm::to_str (W)+".\n");
  }

  /// \brief General constructor: all data are explicitly passed as
  /// arguments (used for instance when reading hills saved on a
  /// file) \param it Time step of creation of the hill \param
  /// weight Weight of the hill \param centers Center of the hill
  /// \param widths Width of the hill around centers \param replica
  /// (optional) Identity of the replica which creates the hill
  inline hill (size_t                    const &it_in,
               cvm::real                 const &W_in,
               std::vector<colvarvalue>  const &centers_in,
               std::vector<cvm::real>    const &widths_in,
               std::string               const &replica_in = "")
    : sW (1.0),
      W (W_in),
      centers (centers_in),
      widths (widths_in),
      it (it_in),
      replica (replica_in)
  {}

  /// Copy constructor
  inline hill (colvarbias_meta::hill const &h)
    : sW (1.0),
      W (h.W),
      centers (h.centers),
      widths (h.widths),
      it (h.it),
      replica (h.replica)
  {}

  /// Destructor
  inline ~hill()
  {}
  
  /// Get the energy
  inline cvm::real energy()
  {
    return W * sW * hill_value;
  }

  /// Get the energy using another hill weight
  inline cvm::real energy (cvm::real const &new_weight)
  {
    return new_weight * sW * hill_value;
  }

  /// Get the current hill value
  inline cvm::real const &value()
  {
    return hill_value;
  }

  /// Set the hill value as specified
  inline void value (cvm::real const &new_value)
  {
    hill_value = new_value;
  }

  /// Get the weight
  inline cvm::real weight()
  {
    return W * sW;
  }

  /// Scale the weight with this factor (by default 1.0 is used)
  inline void scale (cvm::real const &new_scale_fac)
  {
    sW = new_scale_fac;
  }

  /// Get the center of the hill
  inline std::vector<colvarvalue> & center()
  {
    return centers;
  }

  /// Get the i-th component of the center
  inline colvarvalue & center (size_t const &i)
  {
    return centers[i];
  }

  /// Comparison operator
  inline friend bool operator < (hill const &h1, hill const &h2)
  {
    if (h1.it < h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator <= (hill const &h1, hill const &h2)
  {
    if (h1.it <= h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator > (hill const &h1, hill const &h2)
  {
    if (h1.it > h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator >= (hill const &h1, hill const &h2)
  {
    if (h1.it >= h2.it) return true;
    else return false;
  }

  /// Comparison operator
  inline friend bool operator == (hill const &h1, hill const &h2)
  {
    if ( (h1.it >= h2.it) && (h1.replica == h2.replica) ) return true;
    else return false;
  }

  /// Represent the hill ina string suitable for a trajectory file
  std::string output_traj();

  /// Write the hill to an output stream
  inline friend std::ostream & operator << (std::ostream &os,
                                            hill const &h);

};


#endif


// Emacs
// Local Variables:
// mode: C++
// End:
