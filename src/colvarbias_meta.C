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


colvarbias_meta::colvarbias_meta (std::string const &conf, char const *key)
  : colvarbias (conf, key),
    comm (single_replica),
    new_hills_begin (hills.end())
//     free_energy (NULL), free_energy_gradients (NULL),
//     boltzmann_weights (NULL), boltzmann_counts (NULL)
{
  if (cvm::n_abf_biases > 0)
    cvm::log ("Warning: ABF and metadynamics running at the "
              "same time can give inconsistent results.\n");

  get_keyval (conf, "hillWeight", hill_weight, 0.1);
  if (hill_weight == 0.0)
    cvm::log ("Warning: zero weight specified, "
              "this bias will have no effect.\n");

  get_keyval (conf, "newHillFrequency", new_hill_freq, 1000);

  get_keyval (conf, "hillWidth", hill_width, ::sqrt (2.0 * PI) / 2.0);

  {
    bool b_replicas = false;
    get_keyval (conf, "multipleReplicas", b_replicas, false);
    if (b_replicas)
      comm = multiple_replicas;
    else 
      comm = single_replica;
  }

  get_keyval (conf, "useGrids", use_grids, true);
  if (use_grids && (comm != single_replica)) {
    cvm::log ("Warning: calculations with a grid and "
              "multiple replicas is currently not supported: "
              "setting useGrids to \"no\".\n");
    use_grids = false;
  }

  if (use_grids) {
    get_keyval (conf, "gridsUpdateFrequency", grids_freq, new_hill_freq);

    get_keyval (conf, "rebinGrids", rebin_grids, false);

    expand_grids = false;
    for (size_t i = 0; i < colvars.size(); i++) {
      if (colvars[i]->expand_boundaries) {
        expand_grids = true;
        cvm::log ("Will expand the metadynamics grid when the colvar \""+
                  colvars[i]->name+"\" approaches the boundaries.\n");
      }
    }

    //    get_keyval (conf, "expandGrids", expand_grids, false);

    get_keyval (conf, "keepHills", keep_hills, false);
    get_keyval (conf, "dumpFreeEnergyFile", dump_fes, true);
    get_keyval (conf, "saveFreeEnergyFile", dump_fes_save, false);

    for (size_t i = 0; i < colvars.size(); i++) {
      colvars[i]->enable (colvar::task_grid);
    }

    hills_energy           = new colvar_grid_scalar   (colvars);
    hills_energy_gradients = new colvar_grid_gradient (colvars);
  }

  if (comm != single_replica) {

    get_keyval (conf, "replicaID", replica, std::string (""));
    if (!replica.size())
      cvm::fatal_error ("Error: you must define an id for this replica "
                        "when using more than one.\n");

    get_keyval (conf, "replicaFilesRegistry",
                replica_files_registry, 
                (this->name+".replica_files.txt"));

    get_keyval (conf, "replicaUpdateFrequency",
                replica_update_freq, new_hill_freq);

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
              
  if (cvm::debug())
    cvm::log ("Done initializing the metadynamics bias \""+this->name+"\".\n");

  save_delimiters = false;
}


// void colvarbias_meta::parse_analysis (std::string const &conf)
// {
//   // at some point, all of this should be done in the standard run

//   get_keyval (conf, "freeEnergyFile", free_energy_file);

//   get_keyval (conf, "freeEnergyGradientsFile", free_energy_gradients_file);

//   get_keyval (conf, "freeEnergyFirstStep",  free_energy_begin, 0);
//   get_keyval (conf, "freeEnergyLastStep",   free_energy_end,   0);

//   get_keyval (conf, "shiftFreeEnergy",  shift_fes, true);

//   get_keyval (conf, "freeEnergyOffset", free_energy_offset, 0.0);

//   get_keyval (conf, "boltzmannWeightsFile", boltzmann_weights_file);

//   get_keyval (conf, "boltzmannCountsFile", boltzmann_counts_file);

//   get_keyval (conf, "boltzmannWeightsTemp", boltzmann_weights_temp, 300.0);

//   if (get_keyval (conf, "boltzmannWeightsScale",
//                   boltzmann_weights_scale, 1.0)) {
//     if (boltzmann_counts_file.size() && (boltzmann_weights_scale <= 1.0)) {
//       cvm::log ("Warning: boltzmann_weights_scale is too small "
//                 "for a proper discretization: Boltzmann counts "
//                 "will be inaccurate.\n");
//     }
//   }

//   if (free_energy_file.size() ||
//       boltzmann_weights_file.size() ||
//       boltzmann_counts_file.size()) {
//     // all these require the free energy to be read and shifted
//     if (cvm::debug())
//       cvm::log ("Allocating free energy grid.\n");
//     free_energy = new colvar_grid_scalar (colvars);
//   } else
//     free_energy = NULL;
  
//   if (free_energy_gradients_file.size()) {
//     free_energy_gradients = new colvar_grid_gradient (colvars);
//   } else 
//     free_energy_gradients = NULL;

//   if (boltzmann_weights_file.size()) {
//     boltzmann_weights = new colvar_grid_scalar (colvars);
//   } else
//     boltzmann_weights = NULL;

//   if (boltzmann_counts_file.size()) {
//     boltzmann_counts = new colvar_grid_count (colvars);
//   } else
//     boltzmann_counts = NULL;
// }


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

//   if (free_energy) {
//     delete free_energy;
//     free_energy = NULL;
//   }

//   if (free_energy_gradients) {
//     delete free_energy_gradients;
//     free_energy_gradients = NULL;
//   }

//   if (boltzmann_weights) {
//     delete boltzmann_weights;
//     boltzmann_weights = NULL;
//   }

//   if (boltzmann_counts) {
//     delete boltzmann_counts;
//     boltzmann_counts = NULL;
//   }

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
  hill_iter const hills_end = hills.end();
  hills.push_back (h);
  if (new_hills_begin == hills_end) {
    // set the beginning of the "new" hills
    new_hills_begin = hills.end();
    new_hills_begin--;
  }

  if (use_grids) {

    // also add it to the list of hills that are off-grid, which
    // receive special treatment (i.e. they are computed analytically)

    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries (h.centers);
    if (min_dist < (3.0 * ::floor (hill_width)) + 1.0) {
      hills_off_grid.push_back (hills.back());
    }
  }

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

  if (use_grids && hills_off_grid.size()) {
    for (hill_iter hoff = hills_off_grid.begin();
         hoff != hills_off_grid.end(); hoff++) {
      if (*h == *hoff) {
        hills_off_grid.erase (hoff);
        break;
      }
    }
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

    std::vector<int> curr_bin = hills_energy->get_colvars_index();

    if (expand_grids) {

      // first of all, expand the grids, if specified
      if (cvm::debug())
        cvm::log ("Current coordinates on the grid: "+
                  cvm::to_str (curr_bin)+".\n");

      if (cvm::debug())
        cvm::log ("Checking if the grid is big enough.\n");

      bool changed_grids = false;
      int const min_buffer =
        (3 * (size_t) ::floor (hill_width)) + 1;

      std::vector<int>         new_sizes (hills_energy->sizes());
      std::vector<colvarvalue> new_lower_boundaries (hills_energy->lower_boundaries);
      std::vector<colvarvalue> new_upper_boundaries (hills_energy->upper_boundaries);

      for (size_t i = 0; i < colvars.size(); i++) {

        if (! colvars[i]->expand_boundaries)
          continue;

        cvm::real &new_lb   = new_lower_boundaries[i].real_value;
        cvm::real &new_ub   = new_upper_boundaries[i].real_value;
        int       &new_size = new_sizes[i];
        bool changed_lb = false, changed_ub = false;

        if (curr_bin[i] < min_buffer) {
          int const extra_points = (min_buffer - curr_bin[i]);
          new_lb -= extra_points * colvars[i]->width;
          new_size += extra_points;
          // changed offset in this direction => the pointer needs to
          // be changed, too
          curr_bin[i] += extra_points;

          changed_lb = true;
          cvm::log ("Metadynamics"+((cvm::n_meta_biases > 1) ? " \""+name+"\"" : "")+
                    ": new lower boundary for colvar \""+
                    colvars[i]->name+"\", at "+
                    cvm::to_str (new_lower_boundaries[i])+".\n");
        }

        if (curr_bin[i] > new_size - min_buffer - 1) {
          int const extra_points = (curr_bin[i] - (new_size - 1) + min_buffer);
          new_ub += extra_points * colvars[i]->width;
          new_size += extra_points;

          changed_ub = true;
          cvm::log ("Metadynamics"+((cvm::n_meta_biases > 1) ? " \""+name+"\"" : "")+
                    ": new upper boundary for colvar \""+
                    colvars[i]->name+"\", at "+
                    cvm::to_str (new_upper_boundaries[i])+".\n");
        }

        if (changed_lb || changed_ub)
          changed_grids = true;
      }

      if (changed_grids) {

        // map everything into new grids
        if (cvm::debug())
          cvm::log ("Expanding the grids for the metadynamics bias \""+
                    name+"\".\n");

        colvar_grid_scalar *new_hills_energy =
          new colvar_grid_scalar (*hills_energy);
        colvar_grid_gradient *new_hills_energy_gradients =
          new colvar_grid_gradient (*hills_energy_gradients);

        // supply new boundaries to the new grids

        new_hills_energy->lower_boundaries = new_lower_boundaries;
        new_hills_energy->upper_boundaries = new_upper_boundaries;
        new_hills_energy->create (new_sizes, 0.0, 1);

        new_hills_energy_gradients->lower_boundaries = new_lower_boundaries;
        new_hills_energy_gradients->upper_boundaries = new_upper_boundaries;
        new_hills_energy_gradients->create (new_sizes, 0.0, colvars.size());

        new_hills_energy->map_grid (*hills_energy);
        new_hills_energy_gradients->map_grid (*hills_energy_gradients);

        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy = new_hills_energy;
        hills_energy_gradients = new_hills_energy_gradients;

        curr_bin = hills_energy->get_colvars_index();
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
    }
  }

  // syncronise with the other replicas if specified
  if (comm != single_replica) {
    if ((cvm::it % replica_update_freq) == 0) {
      // the buffer should be always emptied to keep the other
      // replicas as much in sync as possible
      replica_out_file.flush();
      // read all replica files, except those for this replica; hills
      // from previous run were already loaded in memory by restarts
      update_replica_files_registry();
      read_replica_files();
    }
  }

  // calculate the biasing energy and forces
  colvar_energy = 0.0;
  for (size_t i = 0; i < colvars.size(); i++) {
    colvar_forces[i].reset();
  }


  if (use_grids) {

    // get the forces from the grid

    if ((cvm::step_relative() % grids_freq) == 0) {
      // map the most recent gaussians to the grids
      project_hills (new_hills_begin, hills.end(),
                     hills_energy,    hills_energy_gradients);
      new_hills_begin = hills.end();
    }

    std::vector<int> curr_bin = hills_energy->get_colvars_index();
    if (cvm::debug())
      cvm::log ("Current coordinates on the grid: "+
                cvm::to_str (curr_bin)+".\n");

    if (hills_energy->index_ok (curr_bin)) {

      // within the grid: add the energy and the forces from there

      colvar_energy += hills_energy->value (curr_bin);
      for (size_t i = 0; i < colvars.size(); i++) {
        cvm::real const *f = &(hills_energy_gradients->value (curr_bin));
        colvar_forces[i].real_value += -1.0 * f[i]; // the gradients
                                                    // are stored, not
                                                    // the forces
      }

    } else {

      // we're off the grid, computing analytically only the hills
      // within range

      calc_hills (hills_off_grid.begin(), hills_off_grid.end(), colvar_energy);
      for (size_t i = 0; i < colvars.size(); i++) {
        calc_hills_force (i, hills_off_grid.begin(), hills_off_grid.end(), colvar_forces);
      }
    }

  } else {

    // calculate the current value of each hill and add it to colvar_energy
    calc_hills (hills.begin(), hills.end(), colvar_energy);
  
    // calculate the current derivatives of each hill and add them to
    // colvar_forces
    for (size_t i = 0; i < colvars.size(); i++) {
      calc_hills_force (i, hills.begin(), hills.end(), colvar_forces);
    }
  }

  if (cvm::debug()) 
    cvm::log ("Hills energy = "+cvm::to_str (colvar_energy)+
              ", hills forces = "+cvm::to_str (colvar_forces)+".\n");
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
      
    // compute the gaussian exponent
    cvm::real cv_sqdev = 0.0;
    for (size_t i = 0; i < colvars.size(); i++) {
      colvarvalue const &x  = curr_values[i];
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      cv_sqdev += (colvars[i]->dist2 (x, center)) / (half_width*half_width);
    }

    // compute the gaussian
    if (cv_sqdev > 23.0) {
      // set it to zero if the exponent is more negative than log(1.0E-05)
      h->value (0.0);
    } else {
      h->value (::exp (-0.5*cv_sqdev));
    }
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
  colvarvalue x (values.size() ? values[i].type() : colvars[i]->type());
  x = (values.size() ? values[i] : colvars[i]->value());

  // do the type check only once (all colvarvalues in the hills series
  // were already saved with their types matching those in the
  // colvars)

  switch (colvars[i]->type()) {

  case colvarvalue::type_scalar:
    for (hill_iter h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].real_value += 
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) * 
          (colvars[i]->dist2_lgrad (x, center)).real_value );
    }
    break;

  case colvarvalue::type_vector:
  case colvarvalue::type_unitvector:
    for (hill_iter h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].rvector_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (colvars[i]->dist2_lgrad (x, center)).rvector_value );
    }
    break;

  case colvarvalue::type_quaternion:
    for (hill_iter h = h_first; h != h_last; h++) {
      if (h->value() == 0.0) continue;
      colvarvalue const &center = h->centers[i];
      cvm::real const    half_width = 0.5 * h->widths[i];
      forces[i].quaternion_value +=
        ( h->weight() * h->value() * (0.5 / (half_width*half_width)) *
          (colvars[i]->dist2_lgrad (x, center)).quaternion_value );
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
        // this is the file from this same job, skip it
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
            cvm::log ("Skipping a replica file already loaded, \""+(*rsi)+"\".\n");
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

  // read hills from the other replicas' files; for each file, resume
  // the position recorded last time
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
        // the boundaries of the colvars may have been changed (note:
        // this second reallocation may be deleted when the new
        // restart format for the grids has kicked in)
        delete hills_energy;
        delete hills_energy_gradients;
        hills_energy = new colvar_grid_scalar (colvars);
        hills_energy_gradients = new colvar_grid_gradient (colvars);
      }

      if ( !(is >> key) ||
           !(key == std::string ("hills_energy")) ||
           !(hills_energy->read_restart (is)) ) {
        cvm::log ("Error: in reading restart information for metadynamics bias \""+
                  this->name+"\".\n");
        is.clear();
        is.seekg (start_pos, std::ios::beg);
        is.setstate (std::ios::failbit);
        return is;
      }

      if ( !(is >> key) ||
           !(key == std::string ("hills_energy_gradients")) || 
           !(hills_energy_gradients->read_restart (is)) ) {
        cvm::log ("Error: in reading restart information for metadynamics bias \""+
                  this->name+"\".\n");
        is.clear();
        is.seekg (start_pos, std::ios::beg);
        is.setstate (std::ios::failbit);
        return is;
      }

    }

    // read the hills explicitly (if there are any)
    while (this->read_hill (is)) {
      if (cvm::debug()) 
        cvm::log ("Read a previously saved hill under the "
                  "metadynamics bias \""+
                  this->name+"\", created at step "+
                  cvm::to_str ((hills.back()).it)+".\n");
    }
    is.clear();
    new_hills_begin = hills.end();

    if (rebin_grids) {

      // allocate new grids, based on the CURRENT boundaries and
      // widths, as per the restarts, and project the grids from the
      // restart file onto them
      colvar_grid_scalar   *new_hills_energy =
        new colvar_grid_scalar (colvars);
      colvar_grid_gradient *new_hills_energy_gradients =
        new colvar_grid_gradient (colvars);
      
      if (keep_hills && (hills.size() > 0)) {

        // if there are hills, recompute the new grids from them
        cvm::log ("rebinGrids and keepHills defined, recomputing "
                  "analytically the energy and force grids; this may take a while...\n");
        project_hills (hills.begin(), hills.end(),
                       new_hills_energy, new_hills_energy_gradients);
        cvm::log ("rebinning done.\n");

      } else {

        // otherwise, use the grids in the restart file

        cvm::log ("rebinGrids defined, mapping energy and forces grid to new grids.\n");
        new_hills_energy->map_grid (*hills_energy);
        new_hills_energy_gradients->map_grid (*hills_energy_gradients);
      }

      delete hills_energy;
      delete hills_energy_gradients;
      hills_energy = new_hills_energy;
      hills_energy_gradients = new_hills_energy_gradients;

      if (hills.size())
        recount_hills_off_grid (hills.begin(), hills.end(), hills_energy);
    }

  } else {

    // read all hills from the replica files
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
    cvm::fatal_error ("Error: corrupted restart information for metadynamics bias \""+
                      this->name+"\": no matching brace at position "+
                      cvm::to_str (is.tellg())+" in the restart file.\n");
    is.setstate (std::ios::failbit);
  }

  if (cvm::debug())
    cvm::log ("colvarbias_meta::read_restart() done\n");

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

  std::vector<colvarvalue> h_centers (colvars.size());
  for (size_t i = 0; i < colvars.size(); i++) {
    h_centers[i].type ((colvars[i]->value()).type()); 
  }
  {
    // it is safer to read colvarvalue objects one at a time; XXX
    // todo: check to remove it later
    std::string centers_input;
    key_lookup (data, "centers", centers_input);
    std::istringstream centers_is (centers_input);
    for (size_t i = 0; i < colvars.size(); i++) {
      centers_is >> h_centers[i];
    }
  }

  std::vector<cvm::real> h_widths (colvars.size());
  get_keyval (data, "widths", h_widths, std::vector<cvm::real> (colvars.size(), 1.0), parse_silent);
  
  if (comm != single_replica)
    get_keyval (data, "replica", h_replica, replica, parse_silent);

  hill_iter const hills_end = hills.end();
  hills.push_back (hill (h_it, h_weight, h_centers, h_widths, h_replica));
  if (new_hills_begin == hills_end) {
    // set the beginning of the "new" hills
    new_hills_begin = hills.end();
    new_hills_begin--;
  }

  if (use_grids) {

    // add this also to the list of hills that are off-grid, which will
    // be computed analytically
    cvm::real const min_dist =
      hills_energy->bin_distance_from_boundaries ((hills.back()).centers);
    if (min_dist < (3.0 * ::floor (hill_width)) + 1.0) {
      hills_off_grid.push_back (hills.back());
    }
  }

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
      hills_energy->write_restart (os);
      os << "  hills_energy_gradients\n";
      hills_energy_gradients->write_restart (os);

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

        // restore the grid to original values
        hills_energy->multiply_constant (-1.0);
        hills_energy->add_constant (max);
      }

    }

    if ( (!use_grids) || keep_hills ) {

      // write all the hills in memory
      for (std::list<hill>::const_iterator h = this->hills.begin();
           h != this->hills.end();
           h++) {
        os << *h;
      }

    } else {

      // write just those that need to be computed analytically
      for (std::list<hill>::const_iterator h = this->hills_off_grid.begin();
           h != this->hills_off_grid.end();
           h++) {
        os << *h;
      }
    }

  } else {
    cvm::log ("Hills from the metadynamics bias \""+
              this->name+"\" have been previously saved in the file \""+
              replica_out_file_name+"\".\n");
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
  for (size_t i = 0; i < centers.size(); i++) {
    os << " ";
    os << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)  << centers[i];
  }

  os << "  ";
  for (size_t i = 0; i < widths.size(); i++) {
    os << " ";
    os << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width) << widths[i];
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
  os << "  weight   "
     << std::setprecision (cvm::en_prec)
     << std::setw (cvm::en_width)
     << h.W << "\n";

  if (h.replica.size())
    os << "  replica   " << h.replica << "\n";

  os << "  centers ";
  for (size_t i = 0; i < (h.centers).size(); i++) {
    os << " "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << h.centers[i];
  }
  os << "\n";

  os << "  widths  ";
  for (size_t i = 0; i < (h.widths).size(); i++) {
    os << " "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << h.widths[i];
  }
  os << "\n";

  os << "}\n";

  return os;
}


void colvarbias_meta::project_hills (colvarbias_meta::hill_iter  h_first,
                                     colvarbias_meta::hill_iter  h_last,
                                     colvar_grid_scalar         *he,
                                     colvar_grid_gradient       *hg,
                                     cvm::real const scale_factor)
{
  if (cvm::debug())
    cvm::log ("Projecting hills.\n");

  std::vector<colvarvalue> colvar_values (colvars.size());
  std::vector<cvm::real> colvar_forces_scalar (colvars.size());

  std::vector<int> he_ix = he->new_index();
  std::vector<int> hg_ix = (hg != NULL) ? hg->new_index() : std::vector<int> (0);
  cvm::real hills_energy_here = 0.0;
  std::vector<colvarvalue> hills_forces_here (colvars.size(), 0.0);
  // loop over the points of the grid
  for ( ;
        he->index_ok (he_ix) && ((hg != NULL) ? hg->index_ok (hg_ix) : true);
        ) {

    for (size_t i = 0; i < colvars.size(); i++) {
      colvar_values[i] = hills_energy->bin_to_value_scalar (he_ix[i], i);
    }

    // loop over the hills and increment the energy grid locally
    hills_energy_here = 0.0;
    calc_hills (h_first, h_last, hills_energy_here, colvar_values);
    he->acc_value (he_ix, scale_factor * hills_energy_here);

    if (hg != NULL) {
      for (size_t i = 0; i < colvars.size(); i++) {
        hills_forces_here[i].reset();
        calc_hills_force (i, h_first, h_last, hills_forces_here, colvar_values);
        colvar_forces_scalar[i] = scale_factor * hills_forces_here[i].real_value;
      }
      hg->acc_force (hg_ix, &(colvar_forces_scalar.front()));
    }

    he->incr (he_ix);
    if (hg != NULL)
      hg->incr (hg_ix);
  }

  if (! keep_hills) {
    hills.erase (hills.begin(), hills.end());
  }
}


void colvarbias_meta::recount_hills_off_grid (colvarbias_meta::hill_iter  h_first,
                                              colvarbias_meta::hill_iter  h_last,
                                              colvar_grid_scalar         *he)
{
  hills_off_grid.clear();

  for (hill_iter h = h_first; h != h_last; h++) {
    cvm::real const min_dist = hills_energy->bin_distance_from_boundaries (h->centers);
    if (min_dist < (3.0 * ::floor (hill_width)) + 1.0) {
      hills_off_grid.push_back (*h);
    }
  }
}


void colvarbias_meta::analyse()
{}


// void colvarbias_meta::analyse()
// {
//   if ((cvm::step_relative() == 0) || cvm::debug())
//     cvm::log ("Performing analysis for the metadynamics bias \""+
//               this->name+"\".\n");

//   cvm::log ("Sorting hills according to their timestep...\n");
//   hills.sort();
//   cvm::log ("Sorting done.\n");


//   if (cvm::step_relative() == 0) {

//     if (free_energy || free_energy_gradients) {

//       cvm::log ("Calculating the free energy surface of the metadynamics bias \""+
//                 this->name+"\".\n");

//       if (!free_energy_end) {
//         free_energy_end = (hills.back()).it;
//       }

//       hill_iter h_first;
//       for (h_first = hills.begin();
//            (h_first != hills.end()) && (h_first->it < free_energy_begin);
//            h_first++) {
//       }

//       hill_iter h_last;
//       for (h_last = h_first;
//            (h_last != hills.end()) && (h_last->it <= free_energy_end);
//            h_last++) {
//       }

//       if (h_first != h_last) {

//         size_t np_total = 0;

//         std::vector<int> fe_ix;
//         colvar_grid_scalar *fe = free_energy;
//         if (fe) {
//           fe_ix = fe->new_index();
//           np_total = fe->number_of_points();
//           cvm::log ("Calculating free energy values on a grid "
//                     "of "+cvm::to_str (np_total)+" points.\n");
//         }

//         std::vector<int> fg_ix;
//         colvar_grid_gradient *fg = free_energy_gradients;
//         if (fg) {
//           fg_ix = fg->new_index();
//           np_total = fg->number_of_points();
//           cvm::log ("Calculating free energy gradients on a grid "
//                     "of "+cvm::to_str (np_total)+" points.\n");
//         }

//         std::vector<int> bw_ix;
//         colvar_grid_scalar *bw = boltzmann_weights;
//         if (bw) {
//           bw_ix = bw->new_index();
//           np_total = bw->number_of_points();
//           cvm::log ("Calculating Boltzmann weights on a grid "
//                     "of "+cvm::to_str (np_total)+" points.\n");
//         }

//         std::vector<int> bc_ix;
//         colvar_grid_count *bc = boltzmann_counts;
//         if (bc) {
//           bc_ix = bc->new_index();
//         }

//         std::vector<colvarvalue> cv_values (colvars.size());
//         for (size_t i = 0; i < colvars.size(); i++) {
//           cv_values[i].type (colvars[i]->type());
//         }


//         // loop over all grids in the same sweep
//         cvm::real free_energy_minimum = 1.0E+12;
//         for (size_t np = 0;
//              ( (fe ? fe->index_ok (fe_ix) : true) &&
//                (fg ? fg->index_ok (fg_ix) : true) ); np++) {

//           for (size_t i = 0; i < colvars.size(); i++) {
//             cv_values[i] = colvars[i]->bin_to_value_scalar (fe_ix[i]);
//           }

//           colvar_energy = 0.0;
//           calc_hills (h_first, h_last, colvar_energy, cv_values);

//           if (fe) {
//             cvm::real const fe_here = -1.0*colvar_energy;
//             if (fe_here < free_energy_minimum)
//               free_energy_minimum = fe_here;
//             // introduce the offset already, but the minimum is
//             // calculated from the hill energy
//             fe->set_value (fe_ix, fe_here - free_energy_offset);
//             fe->incr (fe_ix);
//           }

//           if (fg) {

//             for (size_t i = 0; i < colvars.size(); i++) {
//               colvar_forces[i].reset();
//               calc_hills_force (i, h_first, h_last, colvar_forces, cv_values);
//               colvar_forces[i] *= -1.0;
//               fg->set_value (fg_ix, colvar_forces[i].real_value, i);
//             }
//             fg->incr (fg_ix);
//           }

// #if defined (COLVARS_STANDALONE)
//           std::cerr.setf (std::ios::fixed, std::ios::floatfield); 
//           std::cerr << std::setw (6) << std::setprecision (2)
//                     << 100.0 * double (np) / double (np_total)
//                     << "% done.\r";
// #endif
//         }
// #if defined (COLVARS_STANDALONE)
//         std::cerr << "100.00% done.\n";
// #endif

//         if (shift_fes && fe) {
//           fe_ix = fe->new_index();
//           for ( ; fe->index_ok (fe_ix); fe->incr (fe_ix)) {
//             fe->set_value (fe_ix, fe->value (fe_ix) - free_energy_minimum);
//           }
//         }


//         if (bw || bc) {

//           fe_ix = fe->new_index();
//           if (bw) bw_ix = bw->new_index();
//           if (bc) bc_ix = bc->new_index();

//           for ( ; fe->index_ok (fe_ix); ) {

//             cvm::real const fe_here = fe->value (fe_ix);

//             cvm::real const boltzmann_weight = 
//               boltzmann_weights_scale * 
//               ( (fe_here > 0.0) ? 
//                 ::exp (-1.0 * fe_here / 
//                        (cvm::boltzmann() * boltzmann_weights_temp)) :
//                 1.0 );

//             if (bw) {
//               bw->set_value (bw_ix, boltzmann_weight); 
//               bw->incr (bw_ix);
//             }

//             if (bc) {
//               bc->set_value (bc_ix, size_t (boltzmann_weight));
//               bc->incr (bc_ix);
//             }

//             fe->incr (fe_ix);
//           }
//         }

//         if (free_energy_file.size()) {
//           std::ofstream os (free_energy_file.c_str());
//           os.setf (std::ios::fixed, std::ios::floatfield);
//           os << std::setw (cvm::cv_width) << std::setprecision (cvm::cv_prec);
//           free_energy->write_multicol (os);
//         }

//         if (free_energy_gradients_file.size()) {
//           std::ofstream os (free_energy_gradients_file.c_str());
//           os.setf (std::ios::fixed, std::ios::floatfield);
//           os << std::setw (cvm::cv_width) << std::setprecision (cvm::cv_prec);
//           free_energy_gradients->write_multicol (os);
//         }

//         if (boltzmann_weights_file.size()) {
//           std::ofstream os (boltzmann_weights_file.c_str());
//           os.setf (std::ios::fixed, std::ios::floatfield);
//           os << std::setw (cvm::cv_width) << std::setprecision (cvm::cv_prec);
//           boltzmann_weights->write_multicol (os);
//         }

//         if (boltzmann_counts_file.size()) {
//           std::ofstream os (boltzmann_counts_file.c_str());
//           os.setf (std::ios::fixed, std::ios::floatfield);
//           os << std::setw (cvm::cv_width) << std::setprecision (cvm::cv_prec);
//           boltzmann_counts->write_multicol (os);
//         }
      
//       } else {
//         cvm::log ("Warning: no hills found within the requested interval.\n");
//       }

//     } // end of free energy plotting

//   } // end of first-step analysis
// }
