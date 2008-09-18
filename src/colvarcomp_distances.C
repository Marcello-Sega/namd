#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"



/// \file cvc_distance.cpp \brief Collective variables
/// determining various type of distances between two groups


colvar::distance::distance (std::string const &conf)
  : cvc (conf)
{
  function_type = "distance";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  if (get_keyval (conf, "oneSiteSystemForce", b_1site_force, false)) {
    cvm::log ("Computing system force on group 1 only");
  }
  parse_group (conf, "group1", group1);
  parse_group (conf, "group2", group2);
  x.type (colvarvalue::type_scalar);
}

colvar::distance::distance()
  : cvc ()
{
  function_type = "distance";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  b_1site_force = false;
  x.type (colvarvalue::type_scalar);
}

void colvar::distance::calc_value()
{
  group1.reset_atoms_data();
  group2.reset_atoms_data();

  group1.read_positions();
  group2.read_positions();

  dist_v = cvm::position_distance (group1.center_of_mass(),
                                   group2.center_of_mass());
  x.real_value = dist_v.norm();
}

void colvar::distance::calc_gradients()
{
  cvm::rvector const u = dist_v.unit();
  group1.set_weighted_gradient (-1.0 * u);
  group2.set_weighted_gradient (       u);

  if (b_debug_gradients) {

    if (x_old.type() == colvarvalue::type_notset)
      x_old.type (x.type());

    colvarvalue const dx_real = 0.5 * this->dist2_lgrad (x, x_old);

    colvarvalue dx_grad (x.type());
    dx_grad += fdiff_change (group1);
    dx_grad += fdiff_change (group2);
    
    cvm::log ("cvdiff (grad) = "+cvm::to_str (dx_grad, cvm::cv_width, cvm::cv_prec)+
              ", cvdiff (real)  = "+cvm::to_str (dx_real, cvm::cv_width, cvm::cv_prec)+"\n");

    x_old = x;
  }
}

void colvar::distance::calc_force_invgrads()
{
  group1.read_system_forces();
  if ( b_1site_force ) {
    ft.real_value = -1.0 * (group1.system_force() * dist_v.unit());
  } else {
    group2.read_system_forces();
    ft.real_value = 0.5 * ((group2.system_force() - group1.system_force()) * dist_v.unit());
  }
}

void colvar::distance::calc_Jacobian_derivative()
{
  jd.real_value = x.real_value ? (2.0 / x.real_value) : 0.0;
}

void colvar::distance::apply_force (colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force (force);

  if (!group2.noforce)
    group2.apply_colvar_force (force);
}




colvar::distance_vec::distance_vec (std::string const &conf)
  : distance (conf)
{
  function_type = "distance_vec";
  x.type (colvarvalue::type_vector);
}

colvar::distance_vec::distance_vec()
  : distance()
{
  function_type = "distance_vec";
  x.type (colvarvalue::type_vector);
}

void colvar::distance_vec::calc_value()
{
  group1.reset_atoms_data();
  group2.reset_atoms_data();

  group1.read_positions();
  group2.read_positions();

  cvm::rvector const dist_v = cvm::position_distance (group1.center_of_mass(),
                                                      group2.center_of_mass());
  x.rvector_value = dist_v;
}

void colvar::distance_vec::calc_gradients()
{ 
  // gradients are not stored: a 3x3 matrix for each atom would be
  // needed to store just the identity matrix
}

void colvar::distance_vec::apply_force (colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_force (-1.0 * force.rvector_value);

  if (!group2.noforce)
    group2.apply_force (       force.rvector_value);
}



colvar::distance_z::distance_z (std::string const &conf)
{
  function_type = "distance_z";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  x.type (colvarvalue::type_scalar);

  parse_group (conf, "main", main);
  parse_group (conf, "ref", ref1);
  // this group is optional
  parse_group (conf, "ref2", ref2, true);
 
  if (ref2.size()) {
    cvm::log ("Using axis joining the centers of mass of groups \"ref\" and \"ref2\"");
    fixed_axis = false;
    if (key_lookup (conf, "axis"))
      cvm::log ("Warning: explicit axis definition will be ignored!");
  } else {
    if (get_keyval (conf, "axis", axis, cvm::rvector (0.0, 0.0, 1.0))) {
      if (axis.norm2() == 0.0)
        cvm::fatal_error ("Axis vector is zero!");
      axis = axis.unit();
    }
    fixed_axis = true;
  }

  if (get_keyval (conf, "oneSiteSystemForce", b_1site_force, false)) {
    cvm::log ("Computing system force on group \"main\" only");
  }
}

colvar::distance_z::distance_z()
{
  function_type = "distance_z";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  x.type (colvarvalue::type_scalar);
}

void colvar::distance_z::calc_value()
{
  main.reset_atoms_data();
  ref1.reset_atoms_data();

  main.read_positions();
  ref1.read_positions();

  if (fixed_axis) {
    dist_v = cvm::position_distance (ref1.center_of_mass(),
                                     main.center_of_mass());
  } else {
    ref2.reset_atoms_data();
    ref2.read_positions();

    dist_v = cvm::position_distance (0.5 * (ref1.center_of_mass() + ref2.center_of_mass()),
                                      main.center_of_mass());
    axis = cvm::position_distance (ref1.center_of_mass(), ref2.center_of_mass());
    axis_norm = axis.norm();
    axis = axis.unit();
  }
  x.real_value = axis * dist_v;
}

void colvar::distance_z::calc_gradients()
{
  if (fixed_axis) {
    ref1.set_weighted_gradient (-1.0 * axis);
    main.set_weighted_gradient (       axis);
  } else {
    ref1.set_weighted_gradient ( 1.0 / axis_norm * (
      cvm::position_distance (ref2.center_of_mass(), main.center_of_mass()) - x.real_value * axis ));
    ref2.set_weighted_gradient ( 1.0 / axis_norm * (
      cvm::position_distance (main.center_of_mass(), ref1.center_of_mass()) + x.real_value * axis ));
    main.set_weighted_gradient ( axis );
  }
}

void colvar::distance_z::calc_force_invgrads()
{
  main.read_system_forces();

  if (fixed_axis && !b_1site_force) {
    ref1.read_system_forces();
    ft.real_value = 0.5 * ((main.system_force() - ref1.system_force()) * axis);
  } else {
    ft.real_value = main.system_force() * axis;
  }
}

void colvar::distance_z::calc_Jacobian_derivative()
{
  jd.real_value = 0.0;
}

void colvar::distance_z::apply_force (colvarvalue const &force)
{
  if (!ref1.noforce)
    ref1.apply_colvar_force (force.real_value);

  if (ref2.size() && !ref2.noforce)
    ref2.apply_colvar_force (force.real_value);

  if (!main.noforce)
    main.apply_colvar_force (force.real_value);
}



colvar::distance_xy::distance_xy (std::string const &conf)
  : distance_z (conf)
{
  function_type = "distance_xy";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  x.type (colvarvalue::type_scalar);
}

colvar::distance_xy::distance_xy()
  : distance_z()
{
  function_type = "distance_xy";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  x.type (colvarvalue::type_scalar);
}

void colvar::distance_xy::calc_value()
{
  ref1.reset_atoms_data();
  main.reset_atoms_data();

  ref1.read_positions();
  main.read_positions();

  dist_v = cvm::position_distance (ref1.center_of_mass(),
                                   main.center_of_mass());
  if (!fixed_axis) {
    ref2.reset_atoms_data();
    ref2.read_positions();

    v12 = cvm::position_distance (ref1.center_of_mass(), ref2.center_of_mass());
    axis_norm = v12.norm();
    axis = v12.unit();
  }

  dist_v_ortho = dist_v - (dist_v * axis) * axis;
  x.real_value = dist_v_ortho.norm();
}

void colvar::distance_xy::calc_gradients()
{
  // Intermediate quantity (r_P3 / r_12 where P is the projection
  // of 3 (main) on the plane orthogonal to 12, containing 1 (ref1))
  cvm::real A;
  cvm::real x_inv;

  if (x.real_value == 0.0) return;
  x_inv = 1.0 / x.real_value;

  if (fixed_axis) {
    ref1.set_weighted_gradient (-1.0 * x_inv * dist_v_ortho);
    main.set_weighted_gradient (       x_inv * dist_v_ortho);
  } else {
    v13 = cvm::position_distance (ref1.center_of_mass(), main.center_of_mass());
    A = (dist_v * axis) / axis_norm;

    ref1.set_weighted_gradient ( (A - 1.0) * x_inv * dist_v_ortho);
    ref2.set_weighted_gradient ( -A        * x_inv * dist_v_ortho);
    main.set_weighted_gradient (      1.0  * x_inv * dist_v_ortho);
  }
}

void colvar::distance_xy::calc_force_invgrads()
{
  main.read_system_forces();

  if (fixed_axis && !b_1site_force) {
    ref1.read_system_forces();
    ft.real_value = 0.5 / x.real_value * ((main.system_force() - ref1.system_force()) * dist_v_ortho);
  } else {
    ft.real_value = 1.0 / x.real_value * main.system_force() * dist_v_ortho;
  }
}

void colvar::distance_xy::calc_Jacobian_derivative()
{
  jd.real_value = x.real_value ? (1.0 / x.real_value) : 0.0;
}

void colvar::distance_xy::apply_force (colvarvalue const &force)
{
  if (!ref1.noforce)
    ref1.apply_colvar_force (force.real_value);

  if (ref2.size() && !ref2.noforce)
    ref2.apply_colvar_force (force.real_value);

  if (!main.noforce)
    main.apply_colvar_force (force.real_value);
}



colvar::distance_dir::distance_dir (std::string const &conf)
  : distance (conf)
{
  function_type = "distance_dir";
  x.type (colvarvalue::type_unitvector);
}


colvar::distance_dir::distance_dir()
  : distance()
{
  function_type = "distance_dir";
  x.type (colvarvalue::type_unitvector);
}


void colvar::distance_dir::calc_value()
{
  group1.reset_atoms_data();
  group2.reset_atoms_data();

  group1.read_positions();
  group2.read_positions();

  dist_v = cvm::position_distance (group1.center_of_mass(),
                                   group2.center_of_mass());
  x.rvector_value = dist_v.unit();
}


void colvar::distance_dir::calc_gradients()
{
  // same as for distance_vec; gradients are computed on the fly in apply_force()
}


void colvar::distance_dir::apply_force (colvarvalue const &force)
{
  // remove the radial force component
  cvm::real const iprod = force.rvector_value * x.rvector_value;
  cvm::rvector const force_tang = (1.0 - iprod) * force.rvector_value;

  if (!group1.noforce)
    group1.apply_force (-1.0 * force_tang);

  if (!group2.noforce)
    group2.apply_force (       force_tang);
}




colvar::gyration::gyration (std::string const &conf)
  : cvc (conf)
{
  function_type = "gyration";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  parse_group (conf, "atoms", atoms);
  x.type (colvarvalue::type_scalar);
}


colvar::gyration::gyration()
{
  function_type = "gyration";
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  x.type (colvarvalue::type_scalar);
}


void colvar::gyration::calc_value()
{
  atoms.reset_atoms_data();
  atoms.read_positions();
  atoms.apply_translation ((-1.0) * atoms.center_of_geometry());

  x.real_value = 0.0;
  for (cvm::atom_iter ai = atoms.begin(); ai != atoms.end(); ai++) {
    x.real_value += (ai->mass/atoms.total_mass) * (ai->pos).norm2();
  }
  x.real_value = ::sqrt (x.real_value);
}


void colvar::gyration::calc_gradients()
{
  cvm::real const drdx = 1.0/(cvm::real (atoms.size()) * x.real_value);
  for (cvm::atom_iter ai = atoms.begin(); ai != atoms.end(); ai++) {
    ai->grad = drdx * ai->pos;
  }
}


void colvar::gyration::calc_force_invgrads()
{
  atoms.read_system_forces();

  cvm::real const dxdr = 1.0/x.real_value;
  ft.real_value = 0.0;

  for (cvm::atom_iter ai = atoms.begin(); ai != atoms.end(); ai++) {
    ft.real_value += dxdr * ai->pos * ai->system_force;
  }
}


void colvar::gyration::calc_Jacobian_derivative()
{
  jd = x.real_value ? (3.0 * cvm::real (atoms.size()) - 4.0) / x.real_value : 0.0;
}


void colvar::gyration::apply_force (colvarvalue const &force)
{
  if (!atoms.noforce)
    atoms.apply_colvar_force (force.real_value);
}

