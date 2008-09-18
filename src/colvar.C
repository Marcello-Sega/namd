#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"


// XX TODO make the acf code handle forces as well as values and velocities


colvar::colvar (std::string const &conf)
{
  cvm::log ("Initializing a new collective variable.\n");
  
  get_keyval (conf, "name", this->name,
              (std::string ("colvar")+cvm::to_str (cvm::colvars.size()+1)));

  for (std::vector<colvar *>::iterator cvi = cvm::colvars.begin();
       cvi < cvm::colvars.end();
       cvi++) {
    if ((*cvi)->name == this->name)
      cvm::fatal_error ("Error: this colvar has the same name, \""+this->name+
                        "\", of another colvar.\n");
  }

  for (size_t i = 0; i < task_ntot; i++) {
    tasks[i] = false;
  }

  // read the configuration and set up corresponding instances, for
  // each type of forcepar implemented
#define initialize_forcepars(def_desc,def_config_key,def_class_name)    \
  {                                                                     \
    size_t def_count = 0;                                               \
    std::string def_conf = "";                                          \
    size_t pos = 0;                                                     \
    while ( this->key_lookup (conf,                                     \
                              def_config_key,                           \
                              def_conf,                                 \
                              pos) ) {                                  \
      if (!def_conf.size()) continue;                                   \
      cvm::log ("Initializing "                                         \
                "a new \""+std::string (def_desc)+"\" component"+       \
                (cvm::debug() ? ", with configuration:\n"+def_conf      \
                 : ".\n"));                                             \
      cvm::increase_depth();                                            \
      cvc *cvdp = new colvar::def_class_name (def_conf);          \
      cvdp->check_keywords (def_conf, def_config_key);                  \
      cvm::decrease_depth();                                            \
      if (cvdp != NULL)                                                 \
        cvcs.push_back (cvdp);                                          \
      else cvm::fatal_error ("Error: in allocating definition \""       \
                             def_config_key"\".\n");                    \
      if ( ! cvcs.back()->name.size())                                  \
        cvcs.back()->name = std::string (def_config_key)+               \
          (cvm::to_str (++def_count));                                  \
      if (cvm::debug())                                                 \
        cvm::log ("Done initializing a definition of type \""+          \
                  std::string (def_desc)+"\""+                          \
                  (cvm::debug() ?                                       \
                   ", named \""+cvcs.back()->name+"\""                  \
                   : "")+".\n");                                        \
      def_conf = "";                                                    \
      if (cvm::debug())                                                 \
        cvm::log ("Parsed "+cvm::to_str (cvcs.size())+                  \
                  " definitions at this time.\n");                      \
    }                                                                   \
  }


  initialize_forcepars ("distance",         "distance",       distance);
  initialize_forcepars ("distance vector",  "distanceVec",    distance_vec);
  initialize_forcepars ("distance vector "
                        "direction",        "distanceDir",    distance_dir);
  initialize_forcepars ("distance projection "
                        "on an axis",       "distanceZ",      distance_z);
  initialize_forcepars ("distance projection "
                        "on a plane",       "distanceXY",     distance_xy);

  initialize_forcepars ("coordination "
                        "number",           "coordnum",       coordnum);

  initialize_forcepars ("angle",            "angle",          angle);
  initialize_forcepars ("dihedral",         "dihedral",       dihedral);

  initialize_forcepars ("hydrogen bond",    "hBond",          h_bond);

  initialize_forcepars ("alpha helix",      "alphaDihedrals", alpha_dihedrals);
  initialize_forcepars ("alpha helix",      "alpha",          alpha_angles);

  initialize_forcepars ("absolute "
                        "orientation",      "orientation",    orientation);

  initialize_forcepars ("RMSD",             "rmsd",           rmsd);

  initialize_forcepars ("radius of "
                        "gyration",         "gyration",       gyration);

  if (!cvcs.size())
    cvm::fatal_error ("Error: no definitions were provided "
                      "for this collective variable.\n");

  cvm::log ("All components initialized.\n");


  // this is set false if any of the definitions has an exponent
  // different from 1 in the polynomial
  b_linear = true;

  // these are set false if any of the cvcs has them false
  b_inverse_gradients = true;
  b_Jacobian_force = true;

  this->period = 0.0;

  // check the available features of each cvc
  for (size_t i = 0; i < cvcs.size(); i++) {

    if ((cvcs[i])->sup_np != 1) {
      if (cvm::debug() && b_linear)
        cvm::log ("Warning: You are using a non-linear polynomial "
                  "superposition to define this collective variable, "
                  "some biasing methods may be unavailable.\n");
      b_linear = false;

      if ((cvcs[i])->sup_np < 0) {
        cvm::log ("Warning: you chose a negative exponent in the superposition; "
                  "if you apply forces, the simulation may become unstable "
                  "in case the colvar definition \""+
                  (cvcs[i])->function_type+"\" approaches zero.\n");
      }
    }

    if ((cvcs[i])->period > 0.0) {
      if (!b_linear)
        cvm::fatal_error ("Error: cannot use a non-linear superposition with "
                          "periodic variables.\n");

      if (this->period > 0.0) {
        if (::fabs (this->period - (cvcs[i])->period) > 1.0E-07)
          cvm::fatal_error ("Error: combining periodic collective variables "
                            "definitions with different periods.\n");
      } else 
        this->period = (cvcs[i])->period;
    }

    if (! (cvcs[i])->b_inverse_gradients)
      b_inverse_gradients = false;

    if (! (cvcs[i])->b_Jacobian_derivative)
      b_Jacobian_force = false;

    for (size_t j = i; j < cvcs.size(); j++) {
      if ( (cvcs[i])->type() != (cvcs[j])->type() ) {
        cvm::fatal_error ("ERROR: you are definining this collective variable "
                          "by using components of different types, \""+
                          colvarvalue::type_desc[(cvcs[i])->type()]+
                          "\" and \""+
                          colvarvalue::type_desc[(cvcs[j])->type()]+
                          "\". "
                          "You must use the same type in order to "
                          " sum them together.\n");
      }
    }
  }

  {
    colvarvalue::Type const value_type = (cvcs[0])->type();
    if (cvm::debug())
      cvm::log ("This collective variable is a "+
                colvarvalue::type_desc[value_type]+", corresponding to "+
                cvm::to_str (colvarvalue::dof_num[value_type])+
                " internal degrees of freedom.\n");
    x.type (value_type);
    x_reported.type (value_type);
  }

  get_keyval (conf, "width", width, 1.0);
  if (width <= 0.0)
    cvm::fatal_error ("Error: \"width\" must be positive.\n");

  lower_boundary.type (this->type());
  if (get_keyval (conf, "lowerBoundary", lower_boundary, lower_boundary)) {
    if (this->type() != colvarvalue::type_scalar) {
      cvm::fatal_error ("Error: you requested to define a lower boundary, but "
                        "this colvar does not calculate a scalar number.\n");
    } 
    cvm::log ("Lower boundary defined.\n");
    b_lower_boundary = true;
  } else {
    cvm::log ("Lower boundary was not defined.\n");
    b_lower_boundary = false;
  }


  upper_boundary.type (this->type());
  if (get_keyval (conf, "upperBoundary", upper_boundary, upper_boundary)) {
    if (this->type() != colvarvalue::type_scalar) {
      cvm::fatal_error ("Error: you requested to define an upper boundary, but "
                        "this colvar does not calculate a scalar number.\n");
    }
    cvm::log ("Upper boundary defined.\n");
    b_upper_boundary = true;
  } else {
    cvm::log ("Upper boundary was not defined.\n");
    b_upper_boundary = false;
  }


  // TODO: this test breaks down for variables with spatial periodicity
  // (ex: distanceZ) where the (period > 0.0) criterion is not valid
  
/*
  if (b_lower_boundary && b_upper_boundary && (!(period > 0.0))) {
    if (this->compare (lower_boundary, upper_boundary) > 0.0) {
      cvm::fatal_error ("Error: the lower boundary, "+
                        cvm::to_str (lower_boundary)+", is larger "
                        "than the upper boundary, "+
                        cvm::to_str (upper_boundary)+".\n");
    }
  }
*/

  if (b_lower_boundary)
    if (get_keyval (conf, "lowerWallConstant", lower_wall_k, 0.0)) {
      cvm::log ("Applying a harmonic lower wall at "+
                cvm::to_str (lower_boundary)+", with coefficient "+
                cvm::to_str (lower_wall_k)+".\n");
      enable (task_lower_wall);
    }


  if (b_upper_boundary) 
    if (get_keyval (conf, "upperWallConstant", upper_wall_k, 0.0)) {
      cvm::log ("Applying a harmonic upper wall at "+
                cvm::to_str (upper_boundary)+", with coefficient "+
                cvm::to_str (upper_wall_k)+".\n");
      enable (task_upper_wall);
    }


  {
    bool b_extended_lagrangian;
    get_keyval (conf, "extendedLagrangian", b_extended_lagrangian, false);

    if (b_extended_lagrangian) {

      cvm::log ("Enabling the extended lagrangian term for colvar \""+
                this->name+"\".\n");
    
      enable (task_extended_lagrangian);

      xr.type (this->type());
      vr.type (this->type());
      fr.type (this->type());

      get_keyval (conf, "extendedForceConstant", ext_force_k, 1.0);
      if (ext_force_k <= 0.0)
        cvm::fatal_error ("Error: \"extended_force_constant\" must be positive.\n");

      get_keyval (conf, "extendedFictitiousMass",   ext_mass, 1.0);
      if (ext_mass <= 0.0)
        cvm::fatal_error ("Error: \"extended_fictitious_mass\" must be positive.\n");
    }
  }

  //{
  //  bool b_jacobian_force;
  //  get_keyval (conf, "JacobianForce", b_jacobian_force, false);
  //  if (b_jacobian_force) {
  //    enable (task_Jacobian_force);
  //  }
  //}

  {
    bool b_output_value;
    get_keyval (conf, "outputValue", b_output_value, true);
    if (b_output_value) {
      enable (task_output_value);
    }
  }

  {
    bool b_output_velocity;
    get_keyval (conf, "outputVelocity", b_output_velocity, false);
    if (b_output_velocity) {
      enable (task_output_velocity);
    }
  }

  {
    bool b_output_system_force;
    get_keyval (conf, "outputSystemForce", b_output_system_force, false);
    if (b_output_system_force) {
      enable (task_output_system_force);
    }
  }

  {
    bool b_output_applied_force;
    get_keyval (conf, "outputAppliedForce", b_output_applied_force, false);
    if (b_output_applied_force) {
      enable (task_output_applied_force);
    }
  }

  if (cvm::b_analysis)
    parse_analysis (conf);

  if (cvm::debug())
    cvm::log ("Done initializing collective variable \""+this->name+"\".\n");
}


void colvar::parse_analysis (std::string const &conf) {

  //   if (cvm::debug())
  //     cvm::log ("Parsing analysis flags for collective variable \""+
  //               this->name+"\".\n");

  bool b_runave = false;
  if (get_keyval (conf, "runningAverage", b_runave)) {

    //     cvm::log ("Calculating running average.\n");

    get_keyval (conf, "runAveLength", runave_length, 1000);
    get_keyval (conf, "runAveStride", runave_stride, 1);

    std::string runave_outfile;
    get_keyval (conf, "runAveOutputFile", runave_outfile,
                std::string (this->name+".runave.dat"));

    size_t const this_cv_width = x.output_width (cvm::cv_width);
    runave_os.open (runave_outfile.c_str());
    runave_os << "# " << cvm::wrap_string ("step", cvm::it_width-2)
              << "  "
              << cvm::wrap_string ("running average", this_cv_width)
              << " "
              << cvm::wrap_string ("running stddev", this_cv_width)
              << "\n";
  } else 
    runave_length = 0;
  
  bool b_acf = false;
  if (get_keyval (conf, "acf", b_acf)) {

    //     cvm::log ("Calculating auto-correlation function.\n");

    std::string acf_type_str;
    get_keyval (conf, "acfType", acf_type_str, to_lower_cppstr (std::string ("velocity")));
    if (acf_type_str == to_lower_cppstr (std::string ("coordinate"))) {
      //       cvm::log ("Calculating coordinate ACF for \""+this->name+"\".\n");
      acf_type = acf_coor;
    } else if (acf_type_str == to_lower_cppstr (std::string ("velocity"))) {
      //       cvm::log ("Calculating velocity ACF for \""+this->name+"\".\n");
      acf_type = acf_vel;
      enable (task_fdiff_velocity);
    } else if (acf_type_str == to_lower_cppstr (std::string ("coordinate_p2"))) {
      //       cvm::log ("Calculating second order Legendre polynomial "
      //                 "coordinate ACF for \""+this->name+"\".\n");
      acf_type = acf_p2coor;
    } else {
      cvm::fatal_error ("Unknown type of ACF for \""+this->name+
                        "\", (\""+acf_type_str+"\").\n");
    }
    
    get_keyval (conf, "acfNormalize", acf_normalize, true);
    get_keyval (conf, "acfLength", acf_length, 1000);
    get_keyval (conf, "acfOffset", acf_offset, 0);
    get_keyval (conf, "acfStride", acf_stride, 1);
    get_keyval (conf, "acfOutputFile", acf_outfile,
                std::string (this->name+"."+acf_type_str+".acf.dat"));

  } else {
    acf_length = 0;
  }

  if (runave_length && acf_length)
    cvm::fatal_error ("Error: sorry, cannot calculate running_average "
                     "and acf at the same time.\n");
}


void colvar::enable (colvar::task const &t)
{
  switch (t) {

  case task_output_system_force:
    enable (task_system_force);
    break;

  case task_report_Jacobian_force:
    enable (task_Jacobian_force);
    enable (task_system_force);
    if (cvm::debug())
      cvm::log ("Adding the Jacobian force to the system force, "
                "rather than applying its opposite silently.\n");
    break;

  case task_Jacobian_force:
    enable (task_gradients);

    if (!b_Jacobian_force) 
      cvm::fatal_error ("Error: colvar \""+this->name+
                        "\" does not have Jacobian forces implemented.\n");
    if (!b_linear) 
      cvm::fatal_error ("Error: colvar \""+this->name+
                        "\" must be defined as a linear superposition "
                        "to calculate the Jacobian force.\n");
    if (cvm::debug())
      cvm::log ("Enabling calculation of the Jacobian force "
                "on this colvar.\n");
    fj.type (this->type());
    break;

  case task_system_force:
    if (!b_inverse_gradients)
      cvm::fatal_error ("Error: one or more of the definitions of "
                        "colvar \""+this->name+
                        "\" is unable to calculate system forces.\n");
    ft.type (this->type());
    ft_reported.type (this->type());
    break;

  case task_output_applied_force:
  case task_lower_wall:
  case task_upper_wall:
    // all of the above require gradients
    enable (task_gradients);
    break;

  case task_fdiff_velocity:
    x_old.type (this->type());
    v_fdiff.type (this->type());
    v_reported.type (this->type());
    break;

  case task_output_velocity:
    enable (task_fdiff_velocity);
    break;

  case task_grid:
    if (this->type() != colvarvalue::type_scalar)
      cvm::fatal_error ("Cannot calculate a grid for collective variable, \""+
                        this->name+"\" because it has a non scalar value.\n");
    break;

  case task_extended_lagrangian:
    v_reported.type (this->type());
    break;

  case task_output_value:
  case task_ntot:
    break;

  case task_gradients:
    f.type  (this->type());
    fb.type (this->type());
    break;
  }

  tasks[t] = true;
}


void colvar::disable (colvar::task const &t)
{
  // check dependencies
  switch (t) {
  case task_gradients:
    disable (task_upper_wall);
    disable (task_lower_wall);
    disable (task_output_applied_force);
    disable (task_system_force);
    disable (task_Jacobian_force);
    break;

  case task_system_force:
    disable (task_output_system_force);
    break;

  case task_Jacobian_force:
    disable (task_report_Jacobian_force);
    break;

  case task_fdiff_velocity:
    disable (task_output_velocity);
    break;

  case task_extended_lagrangian:
  case task_report_Jacobian_force:
  case task_output_value:
  case task_output_velocity:
  case task_output_applied_force:
  case task_output_system_force:
  case task_grid:
  case task_lower_wall:
  case task_upper_wall:
  case task_ntot:
    break;
  }

  tasks[t] = false;
}


colvar::~colvar()
{
  if (cvm::b_analysis) {

    if (acf.size()) {
      cvm::log ("Writing acf to file \""+acf_outfile+"\".\n");

      std::ofstream acf_os (acf_outfile.c_str());
      if (! acf_os.good())
        cvm::fatal_error ("Cannot open file \""+acf_outfile+"\".\n");
      write_acf (acf_os);
      acf_os.close();
    }

    if (runave_os.good()) {
      runave_os.close();
    }
  }

  for (size_t i = 0; i < cvcs.size(); i++) {
    delete cvcs[i];
  }
}  



// ******************** CALC FUNCTIONS ********************


void colvar::calc()
{
  if (cvm::debug())
    cvm::log ("Calculating colvar \""+this->name+"\".\n");

  // calculate the value of the colvar

  x.reset();
  if (x.type() == colvarvalue::type_scalar) {

    // scalar variable, polynomial superposition allowed
    for (size_t i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      (cvcs[i])->calc_value();
      cvm::decrease_depth();
      if (cvm::debug())
        cvm::log ("Colvar definition no. "+cvm::to_str (i+1)+
                  " within colvar \""+this->name+"\" has value "+
                  cvm::to_str ((cvcs[i])->value(),
                               cvm::cv_width, cvm::cv_prec)+".\n");
      x += (cvcs[i])->sup_coeff *
        ( ((cvcs[i])->sup_np != 1) ?
          ::pow ((cvcs[i])->value().real_value, (cvcs[i])->sup_np) :
          (cvcs[i])->value().real_value );
    } 
  } else {

    for (size_t i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      (cvcs[i])->calc_value();
      cvm::decrease_depth();
      if (cvm::debug())
        cvm::log ("Colvar definition no. "+cvm::to_str (i+1)+
                  " within colvar \""+this->name+"\" has value "+
                  cvm::to_str ((cvcs[i])->value(),
                               cvm::cv_width, cvm::cv_prec)+".\n");
      x += (cvcs[i])->sup_coeff * (cvcs[i])->value();
    } 
  }

  if (cvm::debug())
    cvm::log ("Colvar \""+this->name+"\" has value "+
              cvm::to_str (x, cvm::cv_width, cvm::cv_prec)+".\n");

  if (tasks[task_gradients]) {
    // calculate the gradients
    for (size_t i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      (cvcs[i])->calc_gradients();
      cvm::decrease_depth();
    }
  }

  if (tasks[task_system_force] && (cvm::step_relative() > 0)) {
    // get from the cvcs the system forces from the PREVIOUS step
    ft.reset();
    for (size_t i = 0; i < cvcs.size(); i++) {
      (cvcs[i])->calc_force_invgrads();
      // linear superposition is assumed
      cvm::increase_depth();
      ft += (cvcs[i])->system_force() / ((cvcs[i])->sup_coeff * cvm::real (cvcs.size()));
      cvm::decrease_depth();
    }
  }

  if (tasks[task_report_Jacobian_force]) {
    // add the Jacobian force to the system force, and don't apply the
    // correction internally: biases such as colvarbias_abf will handle it
    ft += fj;
  }

  if (tasks[task_fdiff_velocity]) {
    // calculate the velocity by finite differences
    if (cvm::step_relative() == 0)
      x_old = x;
    else {
      v_fdiff = fdiff_velocity (x_old, x);
      v_reported = v_fdiff;
    }
  }

  if (tasks[task_extended_lagrangian]) {

    // initialize the restraint center in the first step to the value
    // just calculated from the cvcs; XX TODO: put it in the
    // restart information
    if (cvm::step_relative() == 0)
      xr = x;

    // report the restraint center as "value"
    x_reported = xr;
    v_reported = vr;
    // the "system force" with the extended lagrangian is just the
    // harmonic term
    ft_reported = fr;

  } else {

    x_reported = x;
    ft_reported = ft;
  }

  if (cvm::debug())
    cvm::log ("Done calculating colvar \""+this->name+"\".\n");
}


void colvar::update()
{
  if (cvm::debug())
    cvm::log ("Updating colvar \""+this->name+"\".\n");


  // set to zero the applied force
  f.reset();

  // add the biases' force, which at this point should already have
  // been summed over each bias using this colvar
  f += fb;


  if (tasks[task_lower_wall] || tasks[task_upper_wall]) {

    // Wall force
    colvarvalue fw (this->type());

    // with walls, the colvar is assumed to be always a scalar; but
    // the cvc::compare() function should be used anyway
    // (not sure about this comment - Jerome)

    // if the two are applied concurrently, decide which is the closer
    // (with periodicity, the colvar may be outside of both boundaries
    // at the same time)

    if ( (!tasks[task_upper_wall]) ||
         ( this->dist2 (x, lower_boundary) <
           this->dist2 (x, (upper_boundary +
                            colvarvalue (1.0E-12))) ) ) {
      // if the two boundaries are equal, they are artificially made
      // different just above machine precision

      cvm::real const grad = this->dist2_lgrad (x, lower_boundary);
      if (grad < 0.0) {
        fw = -0.5 * lower_wall_k * grad;
        if (cvm::debug())
          cvm::log ("Applying a lower boundary force ("+
                    cvm::to_str (fw)+") to \""+this->name+"\".\n");
        f += fw;
      }

    } else {

      cvm::real const grad = this->dist2_lgrad (x, upper_boundary);
      if (grad > 0.0) {
        fw = -0.5 * upper_wall_k * grad;
        if (cvm::debug())
          cvm::log ("Applying an upper boundary force ("+
                    cvm::to_str (fw)+") to \""+this->name+"\".\n");
        f += fw;
      }
    }
  }

  if (tasks[task_Jacobian_force]) {
    
    for (size_t i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      (cvcs[i])->calc_Jacobian_derivative();
      cvm::decrease_depth();
    }

    // JH - here we could compute the dot product of the cvc inverse gradient
    // with the colvar gradient, and renormalize.

    fj.reset();
    for (size_t i = 0; i < cvcs.size(); i++) {
      // linear superposition is assumed
      fj += 1.0 / ( cvm::real (cvcs.size()) *  cvm::real ((cvcs[i])->sup_coeff) ) *
        (cvcs[i])->Jacobian_derivative();
    }
    fj *= cvm::boltzmann() * cvm::temperature();

    // the Jacobian force has not been added to the system force, so
    // it's now subtracted from the applied force
    if (! tasks[task_report_Jacobian_force]) 
      f -= fj;
  }


  if (tasks[task_extended_lagrangian]) {

    // the total force is applied on the fictitious mass, while the
    // atoms only feel the harmonic force
    fr   = f;
    fr  += (-0.5 * ext_force_k) * this->dist2_lgrad (xr, x);
    f    = (-0.5 * ext_force_k) * this->dist2_rgrad (xr, x);

    // leap frog
    vr  += (0.5 * cvm::dt) * fr / ext_mass;
    xr  += cvm::dt * vr;
    // if the colvarvalue is set to a type with constraints, apply
    // them
    xr.apply_constraints();
    vr  += (0.5 * cvm::dt) * fr / ext_mass; 
  }


  if (tasks[task_fdiff_velocity]) {
    // set it for the next step
    x_old = x;
  }

  if (cvm::debug())
    cvm::log ("Done updating colvar \""+this->name+"\".\n");
}


void colvar::communicate_forces()
{
  if (cvm::debug())
    cvm::log ("Communicating forces from colvar \""+this->name+"\".\n"); 

  if (x.type() == colvarvalue::type_scalar) {

    for (size_t i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      (cvcs[i])->apply_force (f * (cvcs[i])->sup_coeff * 
                              cvm::real ((cvcs[i])->sup_np) *
                              (::pow ((cvcs[i])->value().real_value,
                                       (cvcs[i])->sup_np-1)) );
      cvm::decrease_depth();
    }

  } else {

    for (size_t i = 0; i < cvcs.size(); i++) {
      cvm::increase_depth();
      (cvcs[i])->apply_force (f * (cvcs[i])->sup_coeff);
      cvm::decrease_depth();
    }
  }

  if (cvm::debug())
    cvm::log ("Done communicating forces from colvar \""+this->name+"\".\n"); 
}



// ******************** METRIC FUNCTIONS ********************
// Use the metrics defined by \link cvc \endlink objects

bool colvar::periodic_boundaries() const
{
  if (this->type() != colvarvalue::type_scalar) {
    cvm::fatal_error ("Error: this colvar is not a scalar value "
                      "and cannot produce a grid.\n");
  }

  if ( (!b_lower_boundary) || (!b_upper_boundary) ) {
    cvm::fatal_error ("Error: requesting to histogram the "
                      "collective variable \""+this->name+"\", but a "
                      "pair of lower and upper boundaries must be "
                      "defined.\n");
  }

  if (period > 0.0) {
    if ( ((::sqrt (this->dist2 (lower_boundary,
                                upper_boundary))) / this->width)
         < 1.0E-06 ) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}


cvm::real colvar::dist2 (colvarvalue const &x1,
                         colvarvalue const &x2) const
{
  return (cvcs[0])->dist2 (x1, x2);
}

colvarvalue colvar::dist2_lgrad (colvarvalue const &x1,
                                 colvarvalue const &x2) const
{
  return (cvcs[0])->dist2_lgrad (x1, x2);
}

colvarvalue colvar::dist2_rgrad (colvarvalue const &x1,
                                 colvarvalue const &x2) const
{
  return (cvcs[0])->dist2_rgrad (x1, x2);
}

cvm::real colvar::compare (colvarvalue const &x1,
                           colvarvalue const &x2) const
{
  return (cvcs[0])->compare (x1, x2);
}



// ******************** INPUT FUNCTIONS ********************

std::istream & colvar::read_restart (std::istream &is)
{
  size_t const start_pos = is.tellg();

  std::string conf;
  if ( !(is >> colvarparse::read_block ("colvar", conf)) ) {
    // this is not a colvar block
    is.clear();
    is.seekg (start_pos, std::ios::beg);
    is.setstate (std::ios::failbit);
    return is;
  }

  {
    std::string check_name = "";
    if ( (get_keyval (conf, "name", check_name,
                      std::string (""), colvarparse::parse_silent)) &&
         (check_name != name) )  {
      cvm::fatal_error ("Error: the state file does not match the "
                        "configuration file, at colvar \""+name+"\".\n");
    }
    if (check_name.size() == 0) {
      cvm::fatal_error ("Error: Collective variable in the "
                        "restart file without any identifier.\n");
    }
  }

  if ( !(get_keyval (conf, "x", x,
                     colvarvalue (x.type()), colvarparse::parse_silent)) ) {
    cvm::log ("Error: restart file does not contain "
              "the value of the colvar \""+
              name+"\" .\n");
  } else {
    cvm::log ("Restarting collective variable \""+name+"\" from value: "+
              cvm::to_str (x)+"\n");
  }

  if (tasks[colvar::task_extended_lagrangian]) {

    if ( !(get_keyval (conf, "extended_x", xr,
                       colvarvalue (x.type()), colvarparse::parse_silent)) &&
         !(get_keyval (conf, "extended_v", vr,
                       colvarvalue (x.type()), colvarparse::parse_silent)) ) {
      cvm::log ("Error: restart file does not contain "
                "\"extended_x\" or \"extended_v\" for the colvar \""+
                name+"\", but you requested \"extendedLagrangian\".\n");
    }
  }

  if (tasks[task_extended_lagrangian]) {
    x_reported = xr;
  } else {
    x_reported = x;
  }

  if (tasks[task_output_velocity]) {

    if ( !(get_keyval (conf, "v", v_fdiff,
                       colvarvalue (x.type()), colvarparse::parse_silent)) ) {
      cvm::log ("Error: restart file does not contain "
                "the velocity for the colvar \""+
                name+"\", but you requested \"outputVelocity\".\n");
    }

    if (tasks[task_extended_lagrangian]) {
      v_reported = vr;
    } else {
      v_reported = v_fdiff;
    }
  }

  if (b_lower_boundary) {
    if (get_keyval (conf, "lowerBoundary", lower_boundary,
                    lower_boundary, parse_silent)) {
      cvm::log ("Reading a new lowerBoundary in the "
                "state file for the colvar \""+name+"\".\n");
    }
  }

  if (b_upper_boundary) {
    if (get_keyval (conf, "upperBoundary", upper_boundary,
                    upper_boundary, parse_silent)) {
      cvm::log ("Reading an new upperBoundary in the "
                "state file for the colvar \""+name+"\".\n");
    }
  }

  return is;
}


std::istream & colvar::read_traj (std::istream &is)
{
  size_t const start_pos = is.tellg();

  if (tasks[task_output_value]) {

    if (!(is >> x)) {
      cvm::log ("Error: in reading the value of colvar \""+
                this->name+"\" from trajectory.\n");
      is.clear();
      is.seekg (start_pos, std::ios::beg);
      is.setstate (std::ios::failbit);
      return is;
    }

    if (tasks[task_extended_lagrangian]) {
      is >> xr;
      x_reported = xr;
    } else {
      x_reported = x;
    }
  }

  if (tasks[task_output_velocity]) {

    is >> v_fdiff;

    if (tasks[task_extended_lagrangian]) {
      is >> vr;
      v_reported = vr;
    } else {
      v_reported = v_fdiff;
    }
  }

  if (tasks[task_output_system_force]) {

    is >> ft;

    if (tasks[task_extended_lagrangian]) { 
      is >> fr;
      ft_reported = fr;
    } else {
      ft_reported = ft;
    }
  }

  if (tasks[task_output_applied_force]) {
    is >> f;
  }

  return is;
}


// ******************** OUTPUT FUNCTIONS ********************

std::ostream & colvar::write_restart (std::ostream &os) {

  os << "colvar {\n"
     << "  name " << name << "\n"
     << "  x "
     << std::setprecision (cvm::cv_prec)
     << std::setw (cvm::cv_width)
     << x << "\n";

  if (tasks[task_output_velocity]) {
    os << "  v "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << x << "\n";
  }

  if (tasks[task_extended_lagrangian]) {
    os << "  extended_x "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << xr << "\n"
       << "  extended_v "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << vr << "\n";
  }

  if (b_lower_boundary)
    os << "  lowerBoundary "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << lower_boundary << "\n";

  if (b_upper_boundary)
    os << "  upperBoundary "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width)
       << upper_boundary << "\n";

  os << "}\n\n";

  return os;
}  


std::ostream & colvar::write_traj_label (std::ostream & os)
{
  size_t const this_cv_width = x.output_width (cvm::cv_width);

  os << " ";

  if (tasks[task_output_value]) {

    os << " "
       << cvm::wrap_string (this->name, this_cv_width);

    if (tasks[task_extended_lagrangian]) {
      // restraint center
      os << " r_"
         << cvm::wrap_string (this->name, this_cv_width-2);
    }
  }

  if (tasks[task_output_velocity]) {

    os << " v_"
       << cvm::wrap_string (this->name, this_cv_width-2);

    if (tasks[task_extended_lagrangian]) {
      // restraint center
      os << " vr_"
         << cvm::wrap_string (this->name, this_cv_width-3);
    }
  }

  if (tasks[task_output_system_force]) {

    os << " fs_"
       << cvm::wrap_string (this->name, this_cv_width-2);

    if (tasks[task_extended_lagrangian]) {
      // restraint center
      os << " fr_"
         << cvm::wrap_string (this->name, this_cv_width-3);
    }
  }

  if (tasks[task_output_applied_force]) {
    os << " fa_"
       << cvm::wrap_string (this->name, this_cv_width-2);
  }

  return os;
}


std::ostream & colvar::write_traj (std::ostream &os)
{
  os << " ";

  if (tasks[task_output_value]) {

    if (tasks[task_extended_lagrangian]) {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << x;
    }

    os << " "
       << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
       << x_reported;
  }

  if (tasks[task_output_velocity]) {

    if (tasks[task_extended_lagrangian]) {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << v_fdiff;
    }

    os << " "
       << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
       << v_reported;
  }

  if (tasks[task_output_system_force]) {

    if (tasks[task_extended_lagrangian]) {
      os << " "
         << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
         << ft;
    }

    os << " "
       << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
       << ft_reported;
  }

  if (tasks[task_output_applied_force]) {
    os << " "
       << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
       << f;
  }

  return os;
}



// ******************** ANALYSIS FUNCTIONS ********************

void colvar::analyse()
{
  // either acf or runave one should manage the store of past values,
  // not both
  if (acf_length > 0) {
    calc_acf();
  } else if (runave_length > 0) {
    calc_runave();
  }
}


inline void store_add_value (size_t const &store_length,
                             std::list<colvarvalue> &store,
                             colvarvalue const &new_value)
{
  store.push_front (new_value);
  if (store.size() > store_length)
    store.pop_back();
}

inline void store_incr (std::list< std::list<colvarvalue> >           &store,
                        std::list< std::list<colvarvalue> >::iterator &store_p)
{
  if ((++store_p) == store.end()) 
    store_p = store.begin();
}


void colvar::calc_acf()
{
  // using here an acf_stride-long list of vectors for either
  // coordinates (x_store) or velocities (v_store); each vector can
  // contain up to acf_length values, which are contiguous in memory
  // representation but separated by acf_stride in the time series;
  // the pointer to each vector is changed at every step

  if (! (x_store.size() || v_store.size()) ) {

    // first-step operations

    acf_nframes = 0;

    if (cvm::debug())
      cvm::log ("Colvar \""+this->name+"\": initializing ACF calculation.\n");

    if (acf.size() < acf_length+1)
      acf.resize (acf_length+1, 0.0);

    switch (acf_type) {

    case acf_vel:
      // allocate space for velocities to be stored
      for (size_t i = 0; i < acf_stride; i++) {
        v_store.push_back (std::list<colvarvalue>());
        //        (v_store.back()).reserve (acf_length);
      }
      v_store_p = v_store.begin();
      break;

    case acf_coor:
    case acf_p2coor:
      // allocate space for coordinates to be stored
      for (size_t i = 0; i < acf_stride; i++) {
        x_store.push_back (std::list<colvarvalue>());
        //        (x_store.back()).reserve (acf_length);
      }
      x_store_p = x_store.begin();
      break;

    default:
      break;
    }

  } else {

    switch (acf_type) {

    case acf_vel:

      if (tasks[task_fdiff_velocity]) {
        // only way to calculate velocities for the moment
        v_reported = v_fdiff = fdiff_velocity (x_old, x);
      }

      calc_vel_acf ((*v_store_p), this->velocity());
      // store this value in the list of previous ones
      store_add_value (acf_length+acf_offset, *v_store_p, this->velocity());
      // if stride is larger than one, cycle the store lists
      store_incr (v_store, v_store_p);
      break;

    case acf_coor:

      calc_coor_acf ((*x_store_p), x);
      store_add_value (acf_length+acf_offset, *x_store_p, x);
      store_incr (x_store, x_store_p);
      break;

    case acf_p2coor:

      calc_p2coor_acf ((*x_store_p), x);
      store_add_value (acf_length+acf_offset, *x_store_p, x);
      store_incr (x_store, x_store_p);
      break;

    default:
      break;
    }
  }

  if (tasks[task_fdiff_velocity]) {
    // set it for the next step
    x_old = x;
  }
}


void colvar::calc_vel_acf (std::list<colvarvalue> &v_store,
                           colvarvalue const      &v)
{
  // loop over stored velocities and add to the ACF, but only the
  // length is sufficient to hold an entire row of ACF values
  if (v_store.size() >= acf_length+acf_offset) {
    std::list<colvarvalue>::iterator  vs_i = v_store.begin();
    std::vector<cvm::real>::iterator acf_i = acf.begin();

    for (size_t i = 0; i < acf_offset; i++)
      vs_i++;

    // current vel with itself
    *(acf_i++) += v.norm2();

    // inner products of previous velocities with current (acf_i and
    // vs_i are updated)
    colvarvalue::inner_opt (v, vs_i, v_store.end(), acf_i); 

    acf_nframes++;
  }
}


void colvar::calc_coor_acf (std::list<colvarvalue> &x_store,
                            colvarvalue const      &x)
{
  // same as above but for coordinates
  if (x_store.size() >= acf_length+acf_offset) {
    std::list<colvarvalue>::iterator  xs_i = x_store.begin();
    std::vector<cvm::real>::iterator acf_i = acf.begin();

    for (size_t i = 0; i < acf_offset; i++)
      xs_i++;

    *(acf_i++) += x.norm2();

    colvarvalue::inner_opt (x, xs_i, x_store.end(), acf_i); 

    acf_nframes++;
  }
}


void colvar::calc_p2coor_acf (std::list<colvarvalue> &x_store,
                              colvarvalue const      &x)
{
  // same as above but with second order Legendre polynomial instead
  // of just the scalar product
  if (x_store.size() >= acf_length+acf_offset) {
    std::list<colvarvalue>::iterator  xs_i = x_store.begin();
    std::vector<cvm::real>::iterator acf_i = acf.begin();

    for (size_t i = 0; i < acf_offset; i++)
      xs_i++;

    // value of P2(0) = 1
    *(acf_i++) += 1.0;

    colvarvalue::p2leg_opt (x, xs_i, x_store.end(), acf_i); 

    acf_nframes++;
  }
}


void colvar::write_acf (std::ostream &os)
{
  if (!acf_nframes)
    cvm::log ("Warning: ACF was not calculated (insufficient frames).\n");
  os.setf (std::ios::scientific, std::ios::floatfield);
  os << "# Autocorrelation function for collective variable \""
     << this->name << "\"\n";
  // one frame is used for normalization, the statistical sample is
  // hence decreased
  os << "# nframes = " << (acf_normalize ?
                           acf_nframes - 1 :
                           acf_nframes) << "\n";

  cvm::real const acf_norm = acf.front() / cvm::real (acf_nframes);
  std::vector<cvm::real>::iterator acf_i;
  size_t it = acf_offset;
  for (acf_i = acf.begin(); acf_i != acf.end(); acf_i++) {

    (*acf_i) /= cvm::real (acf_nframes);

    if (acf_normalize)
      (*acf_i) /= acf_norm;

    os << std::setw (cvm::it_width) << acf_stride * (it++) << " "
       << std::setprecision (cvm::cv_prec)
       << std::setw (cvm::cv_width) << *acf_i << "\n";
  }
}


void colvar::calc_runave()
{
  size_t const store_length = runave_length-1;

  if (!x_store.size()) {

    runave.type (x.type());

    // first-step operations

    if (cvm::debug())
      cvm::log ("Colvar \""+this->name+
                "\": initializing running average calculation.\n");

    acf_nframes = 0;

    x_store.push_back (std::list<colvarvalue>());
    //    (x_store.back()).reserve (store_length);
    x_store_p = x_store.begin();

  } else {

    if ( (cvm::step_relative() % runave_stride) == 0) {

      if ((*x_store_p).size() >= store_length) {

        runave = x;
        for (std::list<colvarvalue>::iterator xs_i = (*x_store_p).begin();
             xs_i != (*x_store_p).end(); xs_i++) {
          runave += (*xs_i);
        }
        runave *= 1.0 / cvm::real (runave_length);
        runave.apply_constraints();

        runave_variance = 0.0;
        runave_variance += this->dist2 (x, runave);
        for (std::list<colvarvalue>::iterator xs_i = (*x_store_p).begin();
             xs_i != (*x_store_p).end(); xs_i++) {
          runave_variance += this->dist2 (x, (*xs_i));
        }
        runave_variance *= 1.0 / cvm::real (store_length);

        runave_os << std::setw (cvm::it_width) << cvm::step_relative()
                  << "  "
                  << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
                  << runave << " "
                  << std::setprecision (cvm::cv_prec) << std::setw (cvm::cv_width)
                  << ::sqrt (runave_variance) << "\n";
      }

      store_add_value (store_length, *x_store_p, x);
    }
  }

}



colvar_grid_count::colvar_grid_count()
  : colvar_grid<size_t>()
{}

colvar_grid_count::colvar_grid_count (std::vector<int> const &nx_i,
                                      size_t const           &def_count)
  : colvar_grid<size_t> (nx_i, def_count)
{}

colvar_grid_count::colvar_grid_count (std::vector<colvar *>  &colvars,
                                      size_t const           &def_count,
                                      size_t const           &bins_scale)
  : colvar_grid<size_t> (colvars, def_count, bins_scale)
{}


colvar_grid_scalar::colvar_grid_scalar()
  : colvar_grid<cvm::real>(), samples (NULL)
{}

colvar_grid_scalar::colvar_grid_scalar (std::vector<int> const &nx_i)
  : colvar_grid<cvm::real> (nx_i, 0.0), samples (NULL)
{}

colvar_grid_scalar::colvar_grid_scalar (std::vector<colvar *> &colvars,
                                        size_t const          &bins_scale)
  : colvar_grid<cvm::real> (colvars, 0.0, bins_scale), samples (NULL)
{}



colvar_grid_gradient::colvar_grid_gradient()
  : colvar_grid<cvm::real>(), samples (NULL)
{}

colvar_grid_gradient::colvar_grid_gradient (std::vector<int> const &nx_i)
  : colvar_grid<cvm::real> (nx_i, 0.0, nx_i.size()), samples (NULL)
{}

colvar_grid_gradient::colvar_grid_gradient (std::vector<colvar *> &colvars,
                                            size_t const          &bins_scale)
  : colvar_grid<cvm::real> (colvars, 0.0, colvars.size(), bins_scale), samples (NULL)
{}


void colvar_grid_gradient::write_1D_integral (std::ostream &os)
{
  cvm::real bin, min, integral;
  std::vector<cvm::real> int_vals;
  
  os << "#       xi            A(xi)\n";

  if ( cv.size() != 1 ) {
    cvm::fatal_error ("Cannot write integral for multi-dimensional gradient grids.");
  }

  integral = 0.0;
  int_vals.push_back ( 0.0 );
  bin = 0.0;
  min = 0.0;

  for (std::vector<int> ix = new_index(); index_ok (ix); incr (ix), bin += 1.0 ) {
	  
    if (samples) {
      size_t const samples_here = samples->value (ix);
      if (samples_here)
        integral += value (ix) / cvm::real (samples_here) * cv[0]->width;
    } else {
      integral += value (ix) * cv[0]->width;
    }

    if ( integral < min ) min = integral;
    int_vals.push_back ( integral );
  }

  bin = 0.0;
  for ( int i = 0; i < nx[0]; i++, bin += 1.0 ) {
    os << std::setw(10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
       << std::setw(16) << std::setprecision (6) << int_vals[i] - min << "\n";
  }

  os << std::setw(10) << cv[0]->lower_boundary.real_value + cv[0]->width * bin << " "
     << std::setw(16) << std::setprecision (6) << int_vals[nx[0]] - min << "\n";

  return;
}




// quaternion_grid::quaternion_grid (std::vector<colvar *>      const &cv_i,
//                                   std::vector<std::string>   const &grid_str)
// {
//   cv = cv_i;

//   std::istringstream is (grid_str[0]);
//   is >> grid_size;

//   min.assign (3, -1.0);
//   max.assign (3,  1.0);
//   np.assign  (3, grid_size);
//   dx.assign  (3, 2.0/(cvm::real (grid_size)));

//   // assumes a uniform grid in the three directions; change
//   // get_value() if you want to use different sizes
//   cvm::log ("Allocating quaternion grid ("+cvm::to_str (np.size())+" dimensional)...");
//   data.create (np, 0.0);
//   cvm::log ("done.\n");
//   if (cvm::debug()) cvm::log ("Grid size = "+data.size());
// }



// Emacs
// Local Variables:
// mode: C++
// End:
