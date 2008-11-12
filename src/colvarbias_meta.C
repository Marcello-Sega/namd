#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

// used to set the absolute path of a replica file
#if defined (WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define CHDIR ::_chdir
#define GETCWD ::_getcwd
#define PATHSEP "\\"
#else
#include <unistd.h>
#define CHDIR ::chdir
#define GETCWD ::getcwd
#define PATHSEP "/"
#endif


#include "colvar.h"
#include "colvarbias_meta.h"


cvm::real const colvarbias_meta::hills_overlap = 0.61;
cvm::real const colvarbias_meta::hills_downscaling = 0.90;


colvarbias_meta::colvarbias_meta (std::string const &conf, char const *key)
  : colvarbias (conf, key),
    comm (single_replica),
    fill_mode (fill_normal),
    free_energy (NULL), free_energy_gradients (NULL),
    boltzmann_weights (NULL), boltzmann_counts (NULL)
{
  if (cvm::n_abf_biases > 0)
    cvm::log ("Warning: ABF and metadynamics running at the "
              "same time can give inconsistent results.\n");

  get_keyval (conf, "hillWeight", hill_weight, 0.1);
  if (hill_weight == 0.0)
    cvm::log ("Warning: zero weight specified, "
              "this bias will have no effect.\n");

  get_keyval (conf, "newHillFrequency", new_hill_freq, 1000);

  get_keyval (conf, "hillWidth", hill_width, 3.0);

  // use the grids?
  get_keyval (conf, "useGrids", use_grids, true);
  if (use_grids) {

    get_keyval (conf, "expandGrids", expand_grids, true);
    get_keyval (conf, "gridsUpdateFrequency", grids_freq, new_hill_freq);
    get_keyval (conf, "dumpFreeEnergyFile", dump_fes, true);
    get_keyval (conf, "saveFreeEnergyFile", dump_fes_save, false);

    for (size_t i = 0; i < colvars.size(); i++) {
      colvars[i]->enable (colvar::task_grid);
    }

    hills_energy = new colvar_grid_scalar (colvars);
    hills_energy_gradients = new colvar_grid_gradient (colvars);
  }

  bool b_replicas = false;
  get_keyval (conf, "multipleReplicas", b_replicas, false);

  if (b_replicas) {

    if (use_grids)
      cvm::fatal_error ("Error: calculations with a grid and multiple replicas "
                        "is currently not supported, set useGrids to \"no\".\n");

    // by default, the replicas are loosely coupled
    comm = multiple_replicas;

    get_keyval (conf, "replicaID", replica, std::string (""));
    if (!replica.size())
      cvm::fatal_error ("Error: you must define an id for this replica "
                        "when using more than one.\n");

    get_keyval (conf, "replicaFilesRegistry",
                replica_files_registry, 
                (this->name+".replica_files.txt"));

    get_keyval (conf, "replicaUpdateFrequency",
                replica_update_freq, new_hill_freq);

    //     if (key_lookup (conf, "replicaExchange")) {
    //       cvm::fatal_error ("Error: sorry, replica exchange not yet implemented.\n");
    //       comm = exchange_replicas;
    //     }      

  }

  if (comm != single_replica) {

    char *pwd = new char[321];
    if (GETCWD (pwd, 320) == NULL)
      cvm::fatal_error ("Error: cannot read the current working directory.\n");
    replica_out_file_name =
      (std::string (pwd)+std::string (PATHSEP)+
       cvm::output_prefix+".colvars."+this->name+
       "."+replica+".hills");
    delete pwd;

    replica_out_file.open (replica_out_file_name.c_str());
    if (!replica_out_file.good())
      cvm::fatal_error ("Error: in opening hills output file \""+
                        replica_out_file_name+"\".\n");
    register_replica_file (replica_out_file_name);

  }

  get_keyval (conf, "writeHillsTrajectory", b_hills_traj, false);
  if (b_hills_traj) {

    std::string const traj_file_name (cvm::output_prefix+
                                      ".colvars."+this->name+
                                      ( (comm != single_replica) ?
                                        ("."+replica) :
                                        ("") )+
                                      ".hills.traj");
    hills_traj_os.open (traj_file_name.c_str());
    if (!hills_traj_os.good())
      cvm::fatal_error ("Error: in opening hills output file \"" +
                        traj_file_name + "\".\n");
  }
              
  // define the fill mode

  if (get_keyval (conf, "hillsTimeCutoff", nt_cutoff, 1000, parse_silent)) {
    fill_mode = strip_older;
  }

  if ( (fill_mode != fill_normal) && (comm != single_replica) )
    cvm::fatal_error ("Error: cannot post-process hills "
                      "when multiple replicas use them concurrently.\n");

  //   if (cvm::b_analysis)
  //     parse_analysis (conf);

  if (cvm::debug())
    cvm::log ("Done initializing the metadynamics bias \""+this->name+"\".\n");
}


void colvarbias_meta::parse_analysis (std::string const &conf)
{
  // at some point, all should be done by the grids

  get_keyval (conf, "freeEnergyFile", free_energy_file);

  get_keyval (conf, "freeEnergyGradientsFile", free_energy_gradients_file);

  get_keyval (conf, "freeEnergyFirstStep",  free_energy_begin, 0);
  get_keyval (conf, "freeEnergyLastStep",   free_energy_end,   0);

  get_keyval (conf, "shiftFreeEnergy",  shift_fes, true);

  get_keyval (conf, "freeEnergyOffset", free_energy_offset, 0.0);

  get_keyval (conf, "boltzmannWeightsFile", boltzmann_weights_file);

  get_keyval (conf, "boltzmannCountsFile", boltzmann_counts_file);

  get_keyval (conf, "boltzmannWeightsTemp", boltzmann_weights_temp, 300.0);

  if (get_keyval (conf, "boltzmannWeightsScale",
                  boltzmann_weights_scale, 1.0)) {
    if (boltzmann_counts_file.size() && (boltzmann_weights_scale <= 1.0)) {
      cvm::log ("Warning: boltzmann_weights_scale is too small "
                "for a proper discretization: Boltzmann counts "
                "will be inaccurate.\n");
    }
  }

  if (free_energy_file.size() ||
      boltzmann_weights_file.size() ||
      boltzmann_counts_file.size()) {
    // all these require the free energy to be read and shifted
    if (cvm::debug())
      cvm::log ("Allocating free energy grid.\n");
    free_energy = new colvar_grid_scalar (colvars);
  } else
    free_energy = NULL;
  
  if (free_energy_gradients_file.size()) {
    free_energy_gradients = new colvar_grid_gradient (colvars);
  } else 
    free_energy_gradients = NULL;

  if (boltzmann_weights_file.size()) {
    boltzmann_weights = new colvar_grid_scalar (colvars);
  } else
    boltzmann_weights = NULL;

  if (boltzmann_counts_file.size()) {
    boltzmann_counts = new colvar_grid_count (colvars);
  } else
    boltzmann_counts = NULL;
}


colvarbias_meta::~colvarbias_meta()
{
  if (hills_energy) {
    delete hills_energy;
    hills_energy = NULL;
  }

  if (hills_energy_gradients) {
    delete hills_energy_gradients;
    hills_energy_gradients = NULL;
  }

  if (free_energy) {
    delete free_energy;
    free_energy = NULL;
  }

  if (free_energy_gradients) {
    delete free_energy_gradients;
    free_energy_gradients = NULL;
  }

  if (boltzmann_weights) {
    delete boltzmann_weights;
    boltzmann_weights = NULL;
  }

  if (boltzmann_counts) {
    delete boltzmann_counts;
    boltzmann_counts = NULL;
  }

  if (replica_out_file.good())
    replica_out_file.close();

  if (hills_traj_os.good())
    hills_traj_os.close();
}



// **********************************************************************
// Hill management member functions
// **********************************************************************

std::list<colvarbias_meta::hill>::const_iterator 
colvarbias_meta::create_hill (colvarbias_meta::hill const &h)
{
  hills.push_back (h);

  // output to trajectory
  if (hills_traj_os.good()) {
    hills_traj_os << (hills.back()).output_traj();
    if (cvm::debug()) {
      hills_traj_os.flush();
    }
  }

  // output to the replica file
  if (comm != single_replica) {
    if (replica_out_file.good()) {
      replica_out_file << hills.back();
      // the buffer should be always emptied to keep the other
      // replicas as much in sync as possible
      replica_out_file.flush();
    } else {
      cvm::fatal_error ("Error in writing to the hills database for replica \""+
                        replica+"\".\n");
    }
  }

  return hills.end();
}


std::list<colvarbias_meta::hill>::const_iterator
colvarbias_meta::delete_hill (hill_iter &h)
{
  if (cvm::debug()) {
    cvm::log ("Deleting hill from the metadynamics bias \""+this->name+
              "\", with step number "+
              cvm::to_str (h->it)+(h->replica.size() ?
                                   ", replica id \""+h->replica :
                                   "")+".\n");
  }

  if (hills_traj_os.good()) {
    // output to the trajectory 
    hills_traj_os << "# DELETED this hill: "
                  << (hills.back()).output_traj()
                  << "\n";
    if (cvm::debug())
      hills_traj_os.flush();
  }

  return hills.erase (h);
}


void colvarbias_meta::update()
{
  if (cvm::debug())
    cvm::log ("Updating the metadynamics bias \""+this->name+"\".\n");

  if (use_grids) {

    // first of all, expand the grids, if defined

    std::vector<int> curr_bin (colvars.size());
    set_bin (curr_bin);
    if (cvm::debug())
      cvm::log ("Coordinates on the current grid: "+
                cvm::to_str (curr_bin)+".\n");

    if (expand_grids) {

      if (cvm::debug())
        cvm::log ("Checking if the grid is big enough.\n");

      bool new_grids = false;
      int const min_buffer =
        (3 * (size_t) ::floor (hill_width)) + 1;
      std::vector<colvarvalue> old_lower_boundaries (colvars.size());
      for (size_t i = 0; i < colvars.size(); i++) {

        old_lower_boundaries[i].type (colvars[i]->type());
        old_lower_boundaries[i] = colvars[i]->lower_boundary;

        bool new_lb = false, new_ub = false;

        int np = hills_energy->number_of_points (i);

        // if the boundaries cover a period for a periodic colvar, it
        // does not make sense to expand the grid
        if (colvars[i]->periodic_boundaries())
          continue;

        // only update the boundary if there is no wall there
        if (colvars[i]->lower_wall_k == 0.0)
          if (curr_bin[i] < min_buffer) {
            int const extra_points = (2*min_buffer - curr_bin[i]);
            colvars[i]->lower_boundary -=
              extra_points * colvars[i]->width;
            curr_bin[i] += extra_points;
            np += extra_points;
            new_lb = true;
            cvm::log ("Setting a new lower boundary for colvar \""+
                      colvars[i]->name+"\", at "+
                      cvm::to_str (colvars[i]->lower_boundary)+".\n");
          }

        // only update the boundary if there is no wall there
        if (colvars[i]->upper_wall_k == 0.0)
          if (curr_bin[i] > np-min_buffer-1) {
            int const extra_points = (curr_bin[i]-np+1+2*min_buffer);
            colvars[i]->upper_boundary +=
              extra_points * colvars[i]->width;
            new_ub = true;
            cvm::log ("Setting a new upper boundary for colvar \""+
                      colvars[i]->name+"\", at "+
                      cvm::to_str (colvars[i]->upper_boundary)+".\n");
          }

        if (new_lb || new_ub)
          new_grids = true;
      }

      if (new_grids) {

        // remap everything into new grids
        cvm::log ("Expanding the grids for the metadynamics bias \""+
                  name+"\".\n");

        colvar_grid_scalar *new_hills_energy =
          new colvar_grid_scalar (colvars);
        new_hills_energy->add_grid (*hills_energy,
                                    old_lower_boundaries);

        colvar_grid_gradient *new_hills_energy_gradients =
          new colvar_grid_gradient (colvars);
        new_hills_energy_gradients->add_grid (*hills_energy_gradients,
                                              old_lower_boundaries);

        delete hills_energy;
        delete hills_energy_gradients;

        hills_energy = new_hills_energy;
        hills_energy_gradients = new_hills_energy_gradients;

        set_bin (curr_bin);
        if (cvm::debug())
          cvm::log ("Coordinates on the new grid: "+
                    cvm::to_str (curr_bin)+".\n");
      }
    }
  }

  // add a new hill if the required time interval has passed
  if ((cvm::it % new_hill_freq) == 0) {

    if (cvm::debug())
      cvm::log ("Adding a new hill under the bias \""+
                this->name+"\", at step "+cvm::to_str (cvm::it)+".\n");

    switch (comm) {

    case single_replica:
      create_hill (hill (hill_weight, colvars, hill_width));
      break;

    case multiple_replicas:
      create_hill (hill (hill_weight, colvars, hill_width, replica));
      break;

    case exchange_replicas:
      cvm::fatal_error ("Error: sorry, replica exchange not yet implemented.\n");
      break;
    }
  }

  // syncronise with the other replicas if specified
  if (comm != single_replica) {
    if ((cvm::it % replica_update_freq) == 0) {
      replica_out_file.flush();
      // read all replica files, except those for this replica; hills
      // from previous run were already loaded in memory by restarts
      update_replica_files_registry();
      read_replica_files();
    }
  }

  // remove the older hills, if specified
  if (fill_mode != fill_normal)
    remove_hills();

  // calculate the biasing energy and forces
  colvar_energy = 0.0;
  for (size_t i = 0; i < colvars.size(); i++) {
    colvar_forces[i].reset();
  }


  if (use_grids) {

    std::vector<int> curr_bin (colvars.size());
    set_bin (curr_bin);
    if (cvm::debug())
      cvm::log ("Current coordinates on the grid: "+
                cvm::to_str (curr_bin)+".\n");

    if ((cvm::step_relative() % grids_freq) == 0) {
      // map all the current gaussians to the grids (and delete them
      // afterwards)
      project_hills (hills.begin(), hills.end(),
                     hills_energy,  hills_energy_gradients);
    }

    // add the energy and the forces
    colvar_energy += hills_energy->value (curr_bin);
    for (size_t i = 0; i < colvars.size(); i++) {
      cvm::real const *f = &(hills_energy_gradients->value (curr_bin));
      colvar_forces[i].real_value += -1.0 * f[i]; // the gradients are stored
    }
  }

  // calculate the current value of each hill and add it to colvar_energy
  calc_hills (hills.begin(), hills.end(), colvar_energy);
  
  // calculate the current derivatives of each hill and add them to
  // colvar_forces
  for (size_t i = 0; i < colvars.size(); i++) {
    calc_hills_force (i, hills.begin(), hills.end(), colvar_forces);
  }

  if (cvm::debug()) 
    cvm::log ("Hills energy = "+cvm::to_str (colvar_energy)+
              ", hills forces = "+cvm::to_str (colvar_forces)+".\n");
}


void colvarbias_meta::remove_hills()
{
  switch (fill_mode) {

  case strip_older:
    // remove the oldest hills
    if (hills.size() && nt_cutoff)
      while ( (cvm::it - (hills.begin())->it) > nt_cutoff ) {
        hill_iter h = hills.begin();
        delete_hill (h); 
      }
    break;

  case fill_normal:
  default:
    break;
  }
}


void colvarbias_meta::scale_hills()
{
  // time of the first hill (not used with more than one replica)
  //   size_t const it_start =
  //     (hills.size() ? (hills.front()).it : cvm::it_restart);

  cvm::real energy = 0.0;
  for (hill_iter h = hills.begin(); h != hills.end(); h++) {
    energy += h->energy();
  }
  int scale_count = 1;
  while (energy > energy_bound_max) {
    if (cvm::debug())
      cvm::log ("Rescaling hills in \""+this->name+"\" bias, cycle no. "+
                cvm::to_str (scale_count++)+".\n");
    energy = 0.0;
    for (hill_iter h = hills.begin(); h != hills.end(); h++) {
      if (h->value() > colvarbias_meta::hills_overlap) {
        h->scale (colvarbias_meta::hills_downscaling);
      }
      energy += h->energy();
    }
  }
}


void colvarbias_meta::calc_hills (colvarbias_meta::hill_iter      h_first,
                                  colvarbias_meta::hill_iter      h_last,
                                  cvm::real                      &energy,
                                  std::vector<colvarvalue> const &colvar_values)
{
  std::vector<colvarvalue> curr_values (colvars.size());
  for (size_t i = 0; i < colvars.size(); i++) {
    curr_values[i].type (colvars[i]->type());
  }

  if (colvar_values.size()) {
    for (size_t i = 0; i < colvars.size(); i++) {
      curr_values[i] = colvar_values[i];
    }
  } else {
    for (size_t i = 0; i < colvars.size(); i++) {
      curr_values[i] = colvars[i]->value();
    }
  }

  for (hill_iter h = h_first; h != h_last; h++) {
      
    // calculate the gaussian exponent
    cvm::real cv_sqdev = 0.0;
    for (size_t i = 0; i < colvars.size(); i++) {
      colvarvalue const &s  = curr_values[i];
      colvarvalue const &s0 = h->s0[i];
      cvm::real const   &ds = h->ds[i];
      cv_sqdev += (colvars[i]->dist2 (s, s0)) / (ds*ds);
    }

    // calculate the gaussian
    h->value (::exp (-0.5*cv_sqdev));
    energy += h->energy();
  }
}


void colvarbias_meta::calc_hills_force (size_t const &i,
                                        colvarbias_meta::hill_iter      h_first,
                                        colvarbias_meta::hill_iter      h_last,
                                        std::vector<colvarvalue>       &forces,
                                        std::vector<colvarvalue> const &values)
{
  // Retrieve the value of the colvar
  colvarvalue s (values.size() ? values[i].type() : colvars[i]->type());
  s = (values.size() ? values[i] : colvars[i]->value());

  // do the type check only once (all colvarvalues in the hills series
  // were already saved with their types matching those in the
  // colvars)

  switch (colvars[i]->type()) {

  case colvarvalue::type_scalar:
    for (hill_iter h = h_first; h != h_last; h++) {
      colvarvalue const &s0 = h->s0[i];
      cvm::real const   &ds = h->ds[i];
      forces[i].real_value += 
        ( h->weight() * h->value() * (0.5 / (ds*ds)) * 
          (colvars[i]->dist2_lgrad (s, s0)).real_value );
    }
    break;

  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
    for (hill_iter h = h_first; h != h_last; h++) {
      colvarvalue const &s0 = h->s0[i];
      cvm::real const   &ds = h->ds[i];
      forces[i].rvector_value +=
        ( h->weight() * h->value() * (0.5 / (ds*ds)) *
          (colvars[i]->dist2_lgrad (s, s0)).rvector_value );
    }
    break;

  case colvarvalue::type_quaternion:
    for (hill_iter h = h_first; h != h_last; h++) {
      colvarvalue const &s0 = h->s0[i];
      cvm::real const   &ds = h->ds[i];
      forces[i].quaternion_value +=
        ( h->weight() * h->value() * (0.5 / (ds*ds)) *
          (colvars[i]->dist2_lgrad (s, s0)).quaternion_value );
    }
    break;

  case colvarvalue::type_notset:
    break;
  }
}



// **********************************************************************
// multiple replicas member functions
// **********************************************************************

void colvarbias_meta::register_replica_file (std::string const &new_file)
{
  std::ifstream reg_is (replica_files_registry.c_str());
  if (reg_is.good()) {
    std::string s;
    while (reg_is >> s)
      if (new_file == s) {
        // already found this replica in registry
        return;
      }
  }
  reg_is.close();

  std::ofstream reg_os (replica_files_registry.c_str(), std::ios::app);
  reg_os << new_file << "\n";
  reg_os.close();
}


void colvarbias_meta::update_replica_files_registry()
{
  if (cvm::debug())
    cvm::log ("Updating the list of replica files, currently containing "+
              cvm::to_str (replica_files.size())+" elements.\n");

  // update the list of replica files by reading those created after restarting
  std::ifstream reg_is (replica_files_registry.c_str());

  if (reg_is.good()) {

    std::string s ("");
    while ((reg_is >> s) && s.size()) {

      if (s == replica_out_file_name) {
        // this is a file for this replica, and skip_this_replica is true 
        if (cvm::debug())
          cvm::log ("Skipping this job's replica file, \""+s+"\".\n");
        s = "";
        continue;
      }

      bool already_loaded = false;
      std::list<std::string>::iterator rsi = replica_files.begin();
      for ( ; rsi != replica_files.end(); rsi++) {
        if (s == *rsi) {
          // this file was already added to the list
          if (cvm::debug())
            cvm::log ("Skipping the file already loaded, \""+(*rsi)+"\".\n");
          already_loaded = true;;
        }
      }

      if ( (replica_files.size() == 0) || (!already_loaded) ) {
        // add this file to the registry
        if (cvm::debug())
          cvm::log ("New replica file found: \""+s+"\".\n");
        replica_files.push_back (s);
        replica_files_pos.push_back (0);
      }

    }

    s = "";
  } else {
    cvm::fatal_error ("Error: cannot read the replica files registry, \""+
                      replica_files_registry+"\".\n");
  }

  reg_is.close();

  if (cvm::debug())
    cvm::log ("The list of replica files now contains "+
              cvm::to_str (replica_files.size())+" elements.\n");
}


void colvarbias_meta::read_replica_files()
{
  // read the hills from the replica files, starting from the
  // positions recorded last time
  std::list<std::string>::iterator rsi = replica_files.begin();
  std::list<size_t>::iterator      pi  = replica_files_pos.begin();
  for ( ; rsi != replica_files.end(); pi++, rsi++) {

    if (cvm::debug()) 
      cvm::log ("Checking for new hills in the replica file \""+
                (*rsi)+"\".\n");

    std::ifstream is ((*rsi).c_str());
    if (is.good()) {

      is.seekg ((*pi), std::ios::beg);
      while (read_hill (is)) {

        if (cvm::debug())
          cvm::log ("Read a previously saved hill under the "
                    "metadynamics bias \""+
                    this->name+"\", created by replica "+
                    (hills.back()).replica+
                    " at step "+
                    cvm::to_str ((hills.back()).it)+".\n");

        // output to trajectory
        if (hills_traj_os.good()) {
          hills_traj_os << (hills.back()).output_traj();
          if (cvm::debug()) {
            hills_traj_os.flush();
          }
        }
      }
      is.clear();
      // store it for the next read
      *pi = is.tellg();
      if (cvm::debug())
        cvm::log ("Stopped reading at position "+
                  cvm::to_str (*pi)+".\n");
    } else {
      cvm::fatal_error ("Error: cannot read from file \""+(*rsi)+"\".\n");
    }

    is.close();
  }

}


// **********************************************************************
// input member functions
// **********************************************************************

std::istream & colvarbias_meta::read_restart (std::istream& is)
{
  size_t const start_pos = is.tellg();

  cvm::log ("Restarting metadynamics bias \""+
            this->name+"\".\n");
  std::string key, brace, conf;
  if ( !(is >> key)   || !(key == "metadynamics") ||
       !(is >> brace) || !(brace == "{") ||
       !(is >> colvarparse::read_block ("configuration", conf)) ) {

    cvm::log ("Error: in reading restart configuration for metadynamics bias \""+
              this->name+"\" at position "+
              cvm::to_str (start_pos)+" in stream.\n");
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  std::string name = "";
  if ( colvarparse::get_keyval (conf, "name", name, std::string (""), colvarparse::parse_silent) &&
       (name != this->name) )
    cvm::fatal_error ("Error: in the restart file, the "
                      "\"metadynamics\" block has a wrong name: different system?\n");
  if (name.size() == 0) {
    cvm::fatal_error ("Error: \"metadynamics\" block in the restart file "
                      "has no identifiers.\n");
  }

  if (this->comm == single_replica) {

    if (use_grids) {

      if (expand_grids) {
        // reallocate the grid, the boundaries of the colvars may have
        // been changed
        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy = new colvar_grid_scalar (colvars);
        hills_energy_gradients = new colvar_grid_gradient (colvars);
      }

      {
        size_t const pos = is.tellg();
        if ( !(is >> key) || !(key == "hills_energy")) {
          cvm::log ("Error: in reading restart configuration for metadynamics bias \""+
                    this->name+"\", expected \"hills_energy\" at position "+
                    cvm::to_str (pos)+" in stream.\n");
          is.clear();
          is.seekg (start_pos, std::ios::beg);
          is.setstate (std::ios::failbit);
          return is;
        }
        hills_energy->read_raw (is);
      }      

      {
        size_t const pos = is.tellg();
        if ( !(is >> key) || !(key == "hills_energy_gradients")) {
          cvm::log ("Error: in reading restart configuration for metadynamics bias \""+
                    this->name+"\", expected \"hills_energy_gradients\" at position "+
                    cvm::to_str (pos)+" in stream.\n");
          is.clear();
          is.seekg (start_pos, std::ios::beg);
          is.setstate (std::ios::failbit);
          return is;
        }
        hills_energy_gradients->read_raw (is);
      }      

    } else {

      // with one replica, the restart file contains all hills
      while (this->read_hill (is)) {
        if (cvm::debug()) 
          cvm::log ("Read a previously saved hill under the "
                    "metadynamics bias \""+
                    this->name+"\", created at step "+
                    cvm::to_str ((hills.back()).it)+".\n");
      }
      is.clear();
    }

  } else {
    // otherwise, read hill from the replica files
    cvm::log ("Multiple replicas have been defined for the "
              "metadynamics bias \""+this->name+"\", skipping "
              "the restart file and reading hills from the replica "
              "files, whose list is contained in \""+
              this->replica_files_registry+"\".\n");
    this->update_replica_files_registry();
    this->read_replica_files();
  }

  is >> brace;
  if (brace != "}") {
    cvm::fatal_error ("Error: corrupt restart information for metadynamics bias \""+
                      this->name+"\": no matching brace at position "+
                      cvm::to_str (is.tellg())+" in the restart file.\n");
    is.setstate (std::ios::failbit);
  }
  return is;
}  


std::istream & colvarbias_meta::read_hill (std::istream &is)
{
  if (!is) return is; // do nothing if failbit is set

  size_t const start_pos = is.tellg();

  std::string data;
  if ( !(is >> read_block ("hill", data)) ) {
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  std::string h_replica = "";

  int h_it;
  get_keyval (data, "step", h_it, 0, parse_silent);
  cvm::real h_weight;
  get_keyval (data, "weight", h_weight, hill_weight, parse_silent);

  std::vector<cvm::real> h_ds (colvars.size());
  get_keyval (data, "ds", h_ds, std::vector<cvm::real> (colvars.size(), 1.0), parse_silent);
 
  std::vector<colvarvalue> h_s0 (colvars.size());
  for (size_t i = 0; i < colvars.size(); i++) {
    h_s0[i].type ((colvars[i]->value()).type()); 
  }
  {
    // it is safer to read colvarvalue objects one at a time; XXX
    // todo: check to remove it later
    std::string s0_input;
    key_lookup (data, "s0", s0_input);
    std::istringstream s0_is (s0_input);
    for (size_t i = 0; i < colvars.size(); i++) {
      s0_is >> h_s0[i];
    }
  }
  
  if (comm != single_replica)
    get_keyval (data, "replica", h_replica, replica, parse_silent);

  hills.push_back ( hill (h_it, h_weight, h_s0, h_ds, h_replica) );
  return is;
}


// **********************************************************************
// output functions
// **********************************************************************

std::ostream & colvarbias_meta::write_restart (std::ostream& os)
{
  os << "metadynamics {\n"
     << "  configuration {\n"
    //      << "    id " << this->id << "\n"
     << "    name " << this->name << "\n";
  if (this->comm != single_replica)
    os << "    replica " << this->replica << "\n";
  os << "  }\n";

  if (this->comm == single_replica) {

    if (use_grids) {

      os << "  hills_energy\n";
      hills_energy->write_raw (os);
      os << "  hills_energy_gradients\n";
      hills_energy_gradients->write_raw (os);

      if (dump_fes) {
        cvm::real const max = hills_energy->maximum_value();
        hills_energy->add_constant (-1.0 * max);
        hills_energy->multiply_constant (-1.0);
        // if this is the only free energy integrator, the pmf file
        // name is general, otherwise there is a label with the bias
        // name
        std::string const fes_file_name =
          ((cvm::n_meta_biases == 1) && (cvm::n_abf_biases == 0)) ?
          std::string (cvm::output_prefix+
                       (dump_fes_save ? "."+cvm::to_str (cvm::step_absolute()) : "")+
                       ".pmf") :
          std::string (cvm::output_prefix+"."+this->name+
                       (dump_fes_save ? "."+cvm::to_str (cvm::step_absolute()) : "")+
                       ".pmf");
        std::ofstream fes_os (fes_file_name.c_str());
        hills_energy->write_multicol (fes_os);
        hills_energy->multiply_constant (-1.0);
        hills_energy->add_constant (max);
      }

    } else {

      for (std::list<hill>::const_iterator h = this->hills.begin();
           h != this->hills.end();
           h++) {
        os << *h;
      }
    }
  } else {
    cvm::log ("Not writing hills added by the metadynamics bias \""+
              this->name+"\", they are already saved in the replica file.\n");
  }

  os << "}\n\n";
  return os;
}  


std::string colvarbias_meta::hill::output_traj()
{
  std::ostringstream os;
  os.setf (std::ios::fixed, std::ios::floatfield);
  os << std::setw (cvm::it_width) << it << " ";

  os.setf (std::ios::scientific, std::ios::floatfield);

  os << "  ";
  for (size_t i = 0; i < s0.size(); i++) {
    os << " ";
    os << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)  << s0[i];
  }

  os << "  ";
  for (size_t i = 0; i < ds.size(); i++) {
    os << " ";
    os << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width) << ds[i];
  }

  os << "  ";
  os << std::setprecision (cvm::en_prec)
     << std::setw (cvm::en_width) << W << "\n";

  return os.str();
}  


std::ostream & operator << (std::ostream &os, colvarbias_meta::hill const &h)
{
  os.setf (std::ios::scientific, std::ios::floatfield);

  os << "hill {\n";
  os << "  step " << std::setw (cvm::it_width) << h.it << "\n";
  os << "  weight  "
     << std::setprecision (cvm::en_prec)
     << std::setw (cvm::en_width)
     << h.W << "\n";

  if (h.replica.size())
    os << "  replica " << h.replica << "\n";

  os << "  s0 ";
  for (size_t i = 0; i < (h.s0).size(); i++) {
    os << " "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << h.s0[i];
  }
  os << "\n";

  os << "  ds ";
  for (size_t i = 0; i < (h.ds).size(); i++) {
    os << " "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << h.ds[i];
  }
  os << "\n";

  os << "}\n";

  return os;
}


void colvarbias_meta::project_hills (colvarbias_meta::hill_iter  h_first,
                                     colvarbias_meta::hill_iter  h_last,
                                     colvar_grid_scalar         *ge,
                                     colvar_grid_gradient       *gg)
{
  if (cvm::debug())
    cvm::log ("Projecting hills.\n");

  std::vector<colvarvalue> colvar_values (colvars.size());
  std::vector<cvm::real> colvar_forces_scalar (colvars.size());

  std::vector<int> ge_ix = ge->new_index();
  std::vector<int> gg_ix = gg->new_index();
  cvm::real colvar_energy_here = 0.0;
  std::vector<colvarvalue> colvar_forces_here (colvars.size());
  for ( ;
        ge->index_ok (ge_ix) && gg->index_ok (gg_ix);
        ge->incr (ge_ix), gg->incr (gg_ix) ) {

    for (size_t i = 0; i < colvars.size(); i++) {
      colvar_values[i] = colvars[i]->bin_to_value_scalar (ge_ix[i]);
    }

    colvar_energy_here = 0.0;
    calc_hills (h_first, h_last, colvar_energy_here, colvar_values);
    // increment the hill energy
    ge->acc_value (ge_ix, colvar_energy_here);

    // increment the hill forces
    
    for (size_t i = 0; i < colvars.size(); i++) {
      colvar_forces_here[i].reset();
      calc_hills_force (i, h_first, h_last, colvar_forces_here, colvar_values);
      colvar_forces_scalar[i] = colvar_forces_here[i].real_value;
    }
    gg->acc_force (gg_ix, &(colvar_forces_scalar.front()));
  }

  hills.erase (hills.begin(), hills.end());
}



// **********************************************************************
// analysis member functions
// **********************************************************************


void colvarbias_meta::analyse()
{
  if ((cvm::step_relative() == 0) || cvm::debug())
    cvm::log ("Performing analysis for the metadynamics bias \""+
              this->name+"\".\n");

  cvm::log ("Sorting hills according to their timestep...\n");
  hills.sort();
  cvm::log ("Sorting done.\n");


  if (cvm::step_relative() == 0) {

    if (free_energy || free_energy_gradients) {

      cvm::log ("Calculating the free energy surface of the metadynamics bias \""+
                this->name+"\".\n");

      if (!free_energy_end) {
        free_energy_end = (hills.back()).it;
      }

      hill_iter h_first;
      for (h_first = hills.begin();
           (h_first != hills.end()) && (h_first->it < free_energy_begin);
           h_first++) {
      }

      hill_iter h_last;
      for (h_last = h_first;
           (h_last != hills.end()) && (h_last->it <= free_energy_end);
           h_last++) {
      }

      if (h_first != h_last) {

        size_t np_total = 0;

        std::vector<int> fe_ix;
        colvar_grid_scalar *fe = free_energy;
        if (fe) {
          fe_ix = fe->new_index();
          np_total = fe->number_of_points();
          cvm::log ("Calculating free energy values on a grid "
                    "of "+cvm::to_str (np_total)+" points.\n");
        }

        std::vector<int> fg_ix;
        colvar_grid_gradient *fg = free_energy_gradients;
        if (fg) {
          fg_ix = fg->new_index();
          np_total = fg->number_of_points();
          cvm::log ("Calculating free energy gradients on a grid "
                    "of "+cvm::to_str (np_total)+" points.\n");
        }

        std::vector<int> bw_ix;
        colvar_grid_scalar *bw = boltzmann_weights;
        if (bw) {
          bw_ix = bw->new_index();
          np_total = bw->number_of_points();
          cvm::log ("Calculating Boltzmann weights on a grid "
                    "of "+cvm::to_str (np_total)+" points.\n");
        }

        std::vector<int> bc_ix;
        colvar_grid_count *bc = boltzmann_counts;
        if (bc) {
          bc_ix = bc->new_index();
        }

        std::vector<colvarvalue> cv_values (colvars.size());
        for (size_t i = 0; i < colvars.size(); i++) {
          cv_values[i].type (colvars[i]->type());
        }


        // loop over all grids in the same sweep
        cvm::real free_energy_minimum = 1.0E+12;
        for (size_t np = 0;
             ( (fe ? fe->index_ok (fe_ix) : true) &&
               (fg ? fg->index_ok (fg_ix) : true) ); np++) {

          for (size_t i = 0; i < colvars.size(); i++) {
            cv_values[i] = colvars[i]->bin_to_value_scalar (fe_ix[i]);
          }

          colvar_energy = 0.0;
          calc_hills (h_first, h_last, colvar_energy, cv_values);

          if (fe) {
            cvm::real const fe_here = -1.0*colvar_energy;
            if (fe_here < free_energy_minimum)
              free_energy_minimum = fe_here;
            // introduce the offset already, but the minimum is
            // calculated from the hill energy
            fe->set_value (fe_ix, fe_here - free_energy_offset);
            fe->incr (fe_ix);
          }

          if (fg) {

            for (size_t i = 0; i < colvars.size(); i++) {
              colvar_forces[i].reset();
              calc_hills_force (i, h_first, h_last, colvar_forces, cv_values);
              colvar_forces[i] *= -1.0;
              fg->set_value (fg_ix, colvar_forces[i].real_value, i);
            }
            fg->incr (fg_ix);
          }

#if defined (COLVARS_STANDALONE)
          std::cerr.setf (std::ios::fixed, std::ios::floatfield); 
          std::cerr << std::setw (6) << std::setprecision (2)
                    << 100.0 * double (np) / double (np_total)
                    << "% done.\r";
#endif
        }
#if defined (COLVARS_STANDALONE)
        std::cerr << "100.00% done.\n";
#endif

        if (shift_fes && fe) {
          fe_ix = fe->new_index();
          for ( ; fe->index_ok (fe_ix); fe->incr (fe_ix)) {
            fe->set_value (fe_ix, fe->value (fe_ix) - free_energy_minimum);
          }
        }


        if (bw || bc) {

          fe_ix = fe->new_index();
          if (bw) bw_ix = bw->new_index();
          if (bc) bc_ix = bc->new_index();

          for ( ; fe->index_ok (fe_ix); ) {

            cvm::real const fe_here = fe->value (fe_ix);

            cvm::real const boltzmann_weight = 
              boltzmann_weights_scale * 
              ( (fe_here > 0.0) ? 
                ::exp (-1.0 * fe_here / 
                       (cvm::boltzmann() * boltzmann_weights_temp)) :
                1.0 );

            if (bw) {
              bw->set_value (bw_ix, boltzmann_weight); 
              bw->incr (bw_ix);
            }

            if (bc) {
              bc->set_value (bc_ix, size_t (boltzmann_weight));
              bc->incr (bc_ix);
            }

            fe->incr (fe_ix);
          }
        }

        if (free_energy_file.size()) {
          std::ofstream os (free_energy_file.c_str());
          os.setf (std::ios::fixed, std::ios::floatfield);
          os << std::setw (cvm::cv_width) << std::setprecision (cvm::cv_prec);
          free_energy->write_multicol (os);
        }

        if (free_energy_gradients_file.size()) {
          std::ofstream os (free_energy_gradients_file.c_str());
          os.setf (std::ios::fixed, std::ios::floatfield);
          os << std::setw (cvm::cv_width) << std::setprecision (cvm::cv_prec);
          free_energy_gradients->write_multicol (os);
        }

        if (boltzmann_weights_file.size()) {
          std::ofstream os (boltzmann_weights_file.c_str());
          os.setf (std::ios::fixed, std::ios::floatfield);
          os << std::setw (cvm::cv_width) << std::setprecision (cvm::cv_prec);
          boltzmann_weights->write_multicol (os);
        }

        if (boltzmann_counts_file.size()) {
          std::ofstream os (boltzmann_counts_file.c_str());
          os.setf (std::ios::fixed, std::ios::floatfield);
          os << std::setw (cvm::cv_width) << std::setprecision (cvm::cv_prec);
          boltzmann_counts->write_multicol (os);
        }
      
      } else {
        cvm::log ("Warning: no hills found within the requested interval.\n");
      }

    } // end of free energy plotting

  } // end of first-step analysis
}
