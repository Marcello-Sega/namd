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




// colvar::distance_vec::distance_vec (std::string const &conf)
//   : distance (conf)
// {
//   function_type = "distance_vec";
//   x.type (colvarvalue::type_vector);
// }

// colvar::distance_vec::distance_vec()
//   : distance()
// {
//   function_type = "distance_vec";
//   x.type (colvarvalue::type_vector);
// }

// void colvar::distance_vec::calc_value()
// {
//   group1.reset_atoms_data();
//   group2.reset_atoms_data();

//   group1.read_positions();
//   group2.read_positions();

//   cvm::rvector const dist_v = cvm::position_distance (group1.center_of_mass(),
//                                                       group2.center_of_mass());
//   x.rvector_value = dist_v;
// }

// void colvar::distance_vec::calc_gradients()
// { 
//   // gradients are not stored: a 3x3 matrix for each atom would be
//   // needed to store just the identity matrix
// }

// void colvar::distance_vec::apply_force (colvarvalue const &force)
// {
//   if (!group1.noforce)
//     group1.apply_force (-1.0 * force.rvector_value);

//   if (!group2.noforce)
//     group2.apply_force (       force.rvector_value);
// }



colvar::distance_z::distance_z (std::string const &conf)
  : cvc (conf)
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



colvar::min_distance::min_distance (std::string const &conf)
  : distance (conf)
{
  function_type = "min_distance";
  x.type (colvarvalue::type_scalar);

  get_keyval (conf, "smoothing", smoothing, (1.0 * cvm::unit_angstrom()));
}

colvar::min_distance::min_distance()
  : distance()
{
  function_type = "min_distance";
  x.type (colvarvalue::type_scalar);
}

void colvar::min_distance::calc_value()
{
  group1.reset_atoms_data();
  group2.reset_atoms_data();

  group1.read_positions();
  group2.read_positions();

  x.real_value = 0.0;

  bool zero_dist = false;

  for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++) {
    for (cvm::atom_iter ai2 = group2.begin(); ai2 != group2.end(); ai2++) {
      cvm::rvector const dv = cvm::position_distance (ai1->pos, ai2->pos);
      cvm::real const d = dv.norm();
      if (d > 0.0)
        x.real_value += ::exp (smoothing / d);
      else
        zero_dist = true;
    }
  }

  x.real_value = zero_dist ? 0.0 : smoothing/(::log (x.real_value));
}

void colvar::min_distance::calc_gradients()
{
  if (x.real_value > 0.0) {
    cvm::real const sum = ::exp (smoothing/x.real_value);
    cvm::real const dxdsum = -1.0 *
      (x.real_value/smoothing) * (x.real_value/smoothing) *
      (1.0 / sum);

    for (cvm::atom_iter ai1 = group1.begin(); ai1 != group1.end(); ai1++) {
      for (cvm::atom_iter ai2 = group2.begin(); ai2 != group2.end(); ai2++) {
        cvm::rvector const dv = cvm::position_distance (ai1->pos, ai2->pos);
        cvm::real const d = dv.norm();
        if (d > 0.0) {
          cvm::rvector const dvu = dv / dv.norm();
          ai1->grad += dxdsum * ::exp (smoothing / d) *
            smoothing * (-1.0/(d*d)) * (-1.0) * dvu;
          ai2->grad += dxdsum * ::exp (smoothing / d) *
            smoothing * (-1.0/(d*d)) * dvu;
        }
      }
    }
  }
}

void colvar::min_distance::apply_force (colvarvalue const &force)
{
  if (!group1.noforce)
    group1.apply_colvar_force (force.real_value);

  if (!group2.noforce)
    group2.apply_colvar_force (force.real_value);
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
  // gradients are computed on the fly within apply_force()
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



colvar::rmsd::rmsd (std::string const &conf)
  : orientation (conf)
{
  // TODO: refuse to start if hideJacobian is on

  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  function_type = "rmsd";
  x.type (colvarvalue::type_scalar);

  ref_pos_sum2 = 0.0;
  for (size_t i = 0; i < ref_pos.size(); i++) {
    ref_pos_sum2 += ref_pos[i].norm2();
  }
}

  
void colvar::rmsd::calc_value()
{
  atoms.reset_atoms_data();
  atoms.read_positions();

  atoms_cog = atoms.center_of_geometry();
  rot.calc_optimal_rotation (ref_pos, atoms.positions_shifted (-1.0 * atoms_cog));

  cvm::real group_pos_sum2 = 0.0;
  for (size_t i = 0; i < atoms.size(); i++) {
    group_pos_sum2 += (atoms[i].pos - atoms_cog).norm2();
  }

  // value of the RMSD (Coutsias et al)
  cvm::real const MSD = 1.0/(cvm::real (atoms.size())) *
    ( group_pos_sum2 + ref_pos_sum2 - 2.0 * rot.lambda );

  x.real_value = (MSD > 0.0) ? ::sqrt (MSD) : 0.0;
}


void colvar::rmsd::calc_gradients()
{
  cvm::real const drmsddx2 = (x.real_value > 0.0) ?
    0.5 / (x.real_value * cvm::real (atoms.size())) :
    0.0;

  for (size_t ia = 0; ia < atoms.size(); ia++) {
    atoms[ia].grad  = (drmsddx2 * 2.0 * (atoms[ia].pos - atoms_cog -
                                         rot.q.rotate (ref_pos[ia])));
  }
}


void colvar::rmsd::apply_force (colvarvalue const &force)
{
  if (!atoms.noforce)
    atoms.apply_colvar_force (force.real_value);
}


void colvar::rmsd::calc_force_invgrads()
{
  atoms.read_system_forces();
  ft.real_value = 0.0;
    
  // Note: gradient square norm is 1/N_atoms
          
  for (size_t ia = 0; ia < atoms.size(); ia++) {
    ft.real_value += atoms[ia].grad * atoms[ia].system_force;
  }
  ft.real_value *= atoms.size();
}


void colvar::rmsd::calc_Jacobian_derivative()
{
  // divergence of the back-rotated target coordinates
  cvm::real divergence = 0.0;
 
  // gradient of the rotation matrix
  cvm::matrix2d <cvm::rvector, 3, 3> grad_rot_mat;

  // gradients of products of 2 quaternion components 
  cvm::rvector g11, g22, g33, g01, g02, g03, g12, g13, g23;

  for (size_t ia = 0; ia < atoms.size(); ia++) {

    // Gradient of optimal quaternion wrt current Cartesian position
    cvm::vector1d< cvm::rvector, 4 >      &dq = rot.dQ0_2[ia];

    g11 = 2.0 * (rot.q)[1]*dq[1];
    g22 = 2.0 * (rot.q)[2]*dq[2];
    g33 = 2.0 * (rot.q)[3]*dq[3];
    g01 = (rot.q)[0]*dq[1] + (rot.q)[1]*dq[0];
    g02 = (rot.q)[0]*dq[2] + (rot.q)[2]*dq[0];
    g03 = (rot.q)[0]*dq[3] + (rot.q)[3]*dq[0];
    g12 = (rot.q)[1]*dq[2] + (rot.q)[2]*dq[1];
    g13 = (rot.q)[1]*dq[3] + (rot.q)[3]*dq[1];
    g23 = (rot.q)[2]*dq[3] + (rot.q)[3]*dq[2];

    // Gradient of the rotation matrix wrt current Cartesian position
    grad_rot_mat[0][0] = -2.0 * (g22 + g33); 
    grad_rot_mat[1][0] =  2.0 * (g12 + g03); 
    grad_rot_mat[2][0] =  2.0 * (g13 - g02); 
    grad_rot_mat[0][1] =  2.0 * (g12 - g03); 
    grad_rot_mat[1][1] = -2.0 * (g11 + g33); 
    grad_rot_mat[2][1] =  2.0 * (g01 + g23); 
    grad_rot_mat[0][2] =  2.0 * (g02 + g13); 
    grad_rot_mat[1][2] =  2.0 * (g23 - g01); 
    grad_rot_mat[2][2] = -2.0 * (g11 + g22); 

    cvm::atom_pos &y = ref_pos[ia]; 

    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        divergence += grad_rot_mat[i][j][i] * y[j];
      }
    }
  }

  jd.real_value = x.real_value > 0.0 ? (3.0 * atoms.size() - 4.0 - divergence) / x.real_value : 0.0;
}



colvar::logmsd::logmsd (std::string const &conf)
  : orientation (conf)
{
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  function_type = "logmsd";
  x.type (colvarvalue::type_scalar);

  ref_pos_sum2 = 0.0;
  for (size_t i = 0; i < ref_pos.size(); i++) {
    ref_pos_sum2 += ref_pos[i].norm2();
  }
}

  
void colvar::logmsd::calc_value()
{
  atoms.reset_atoms_data();
  atoms.read_positions();

  if (cvm::debug())
    cvm::log ("colvar::logmsd: current com: "+
              cvm::to_str (atoms.center_of_mass())+"\n");

  atoms_cog = atoms.center_of_geometry();
  rot.calc_optimal_rotation (ref_pos, atoms.positions_shifted (-1.0 * atoms_cog));

  cvm::real group_pos_sum2 = 0.0;
  for (size_t i = 0; i < atoms.size(); i++) {
    group_pos_sum2 += (atoms[i].pos-atoms_cog).norm2();
  }

  // value of the MSD (Coutsias et al)
  MSD = 1.0/(cvm::real (atoms.size())) *
    ( group_pos_sum2 + ref_pos_sum2 - 2.0 * rot.lambda );

  x.real_value = (MSD > 0.0) ? ::log(MSD) : 0.0;
}


void colvar::logmsd::calc_gradients()
{
  cvm::real fact = (MSD > 0.0) ? 2.0/(cvm::real (atoms.size()) * MSD) : 0.0;

  for (size_t ia = 0; ia < atoms.size(); ia++) {
    atoms[ia].grad = fact * (atoms[ia].pos - atoms_cog - rot.dL0_2[ia]);
  }
}


void colvar::logmsd::apply_force (colvarvalue const &force)
{
  if (!atoms.noforce)
    atoms.apply_colvar_force (force.real_value);
}


void colvar::logmsd::calc_force_invgrads()
{
  atoms.read_system_forces();
  ft.real_value = 0.0;
    
  // Note: gradient square norm is 4.0 / (N_atoms * E)
          
  for (size_t ia = 0; ia < atoms.size(); ia++) {
    ft.real_value += atoms[ia].grad * atoms[ia].system_force;
  }
  ft.real_value *= atoms.size() * MSD / 4.0;
}


void colvar::logmsd::calc_Jacobian_derivative()
{
  // divergence of the back-rotated target coordinates
  cvm::real divergence = 0.0;
 
  // gradient of the rotation matrix
  cvm::matrix2d <cvm::rvector, 3, 3> grad_rot_mat;

  // gradients of products of 2 quaternion components 
  cvm::rvector g11, g22, g33, g01, g02, g03, g12, g13, g23;
 
  for (size_t ia = 0; ia < atoms.size(); ia++) {

    // Gradient of optimal quaternion wrt current Cartesian position
    cvm::vector1d< cvm::rvector, 4 >      &dq = rot.dQ0_2[ia];

    g11 = 2.0 * (rot.q)[1]*dq[1];
    g22 = 2.0 * (rot.q)[2]*dq[2];
    g33 = 2.0 * (rot.q)[3]*dq[3];
    g01 = (rot.q)[0]*dq[1] + (rot.q)[1]*dq[0];
    g02 = (rot.q)[0]*dq[2] + (rot.q)[2]*dq[0];
    g03 = (rot.q)[0]*dq[3] + (rot.q)[3]*dq[0];
    g12 = (rot.q)[1]*dq[2] + (rot.q)[2]*dq[1];
    g13 = (rot.q)[1]*dq[3] + (rot.q)[3]*dq[1];
    g23 = (rot.q)[2]*dq[3] + (rot.q)[3]*dq[2];

    // Gradient of the rotation matrix wrt current Cartesian position
    // Note: we are only going to use "diagonal" terms: grad_rot_mat[i][j][i]
    grad_rot_mat[0][0] = -2.0 * (g22 + g33); 
    grad_rot_mat[1][0] =  2.0 * (g12 + g03); 
    grad_rot_mat[2][0] =  2.0 * (g13 - g02); 
    grad_rot_mat[0][1] =  2.0 * (g12 - g03); 
    grad_rot_mat[1][1] = -2.0 * (g11 + g33); 
    grad_rot_mat[2][1] =  2.0 * (g01 + g23); 
    grad_rot_mat[0][2] =  2.0 * (g02 + g13); 
    grad_rot_mat[1][2] =  2.0 * (g23 - g01); 
    grad_rot_mat[2][2] = -2.0 * (g11 + g22); 

    cvm::atom_pos &y = ref_pos[ia]; 

    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        divergence += grad_rot_mat[i][j][i] * y[j];
      }
    }
  }

  jd.real_value = (3.0 * atoms.size() - 3.0 - divergence) / 2.0;
}



colvar::eigenvector::eigenvector (std::string const &conf)
  : cvc (conf)
{
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  function_type = "eigenvector";
  x.type (colvarvalue::type_scalar);

  parse_group (conf, "atoms", atoms);

  if (get_keyval (conf, "refPositions", ref_pos, ref_pos)) {
    cvm::log ("Using reference positions from input file.\n");
    if (ref_pos.size() != atoms.size()) {
      cvm::fatal_error ("Error: reference positions do not "
                        "match the number of atom indexes.\n");
    }
  }

  {
    std::string file_name;
    if (get_keyval (conf, "refPositionsFile", file_name)) {

      std::string file_col;
      get_keyval (conf, "refPositionsCol", file_col, std::string ("O"));

      double file_col_value;
      bool found = get_keyval (conf, "refPositionsColValue", file_col_value, 0.0);
      if (found && !file_col_value)
        cvm::fatal_error ("Error: refPositionsColValue, "
                          "if provided, must be non-zero.\n");

      ref_pos.resize (atoms.size());
      cvm::load_coords (file_name.c_str(), ref_pos, file_col, file_col_value);
    }
  }

  // now load the eigenvector
  if (get_keyval (conf, "vector", eigenvec, eigenvec)) {
    cvm::log ("Using reference positions from input file.\n");
    if (eigenvec.size() != atoms.size()) {
      cvm::fatal_error ("Error: reference positions do not "
                        "match the number of atom indexes.\n");
    }
  }

  {
    std::string file_name;
    if (get_keyval (conf, "vectorFile", file_name)) {

      std::string file_col;
      get_keyval (conf, "vectorCol", file_col, std::string ("O"));

      double file_col_value;
      bool found = get_keyval (conf, "vectorColValue", file_col_value, 0.0);
      if (found && !file_col_value)
        cvm::fatal_error ("Error: eigenvectorColValue, "
                          "if provided, must be non-zero.\n");

      eigenvec.resize (atoms.size());
      cvm::load_coords (file_name.c_str(), eigenvec, file_col, file_col_value);
    }
  }

  if (!ref_pos.size() || !eigenvec.size()) {
    cvm::fatal_error ("Error: must define both reference "
                      "coordinates and eigenvector.\n");
  }

  cvm::rvector center (0.0, 0.0, 0.0);
  eigenvec_invnorm2 = 0.0;

  for (size_t i = 0; i < atoms.size(); i++) {
    center += eigenvec[i];
  }

  cvm::log ("Subtracting sum of eigenvector components: " + cvm::to_str (center) + "\n");

  for (size_t i = 0; i < atoms.size(); i++) {
    eigenvec[i] = eigenvec[i] - center;
    eigenvec_invnorm2 += eigenvec[i].norm2();
  }
  eigenvec_invnorm2 = 1.0 / eigenvec_invnorm2;

  // request derivatives of optimal rotation wrt 2nd group
  // for Jacobian
  atoms.rot.request_group1_gradients(atoms.size());
  atoms.rot.request_group2_gradients(atoms.size());
}

  
void colvar::eigenvector::calc_value()
{
  atoms.reset_atoms_data();
  atoms.read_positions();

  x.real_value = 0.0;
  for (size_t i = 0; i < atoms.size(); i++) {
    x.real_value += (atoms[i].pos - ref_pos[i]) * eigenvec[i];
  }
}


void colvar::eigenvector::calc_gradients()
{
  for (size_t i = 0; i < atoms.size(); i++) {
    atoms[i].grad = eigenvec[i];
  }

  // WARNING: either the gradient must be rotated, or the automatic
  // rotation of forces within the atom group should be disabled.
}


void colvar::eigenvector::apply_force (colvarvalue const &force)
{
  if (!atoms.noforce)
    atoms.apply_colvar_force (force.real_value);
}


void colvar::eigenvector::calc_force_invgrads()
{
  atoms.read_system_forces();
  ft.real_value = 0.0;
    
  for (size_t ia = 0; ia < atoms.size(); ia++) {
    ft.real_value += eigenvec_invnorm2 * atoms[ia].grad *
      atoms[ia].system_force;
  }
}


void colvar::eigenvector::calc_Jacobian_derivative()
{
  // gradient of the rotation matrix
  cvm::matrix2d <cvm::rvector, 3, 3> grad_rot_mat;
  cvm::quaternion &quat0 = atoms.rot.q;

  // gradients of products of 2 quaternion components 
  cvm::rvector g11, g22, g33, g01, g02, g03, g12, g13, g23;

  cvm::atom_pos x_relative; 
  cvm::real sum = 0.0;

  for (size_t ia = 0; ia < atoms.size(); ia++) {

    // Gradient of optimal quaternion wrt current Cartesian position
    // trick: d(R^-1)/dx = d(R^t)/dx = (dR/dx)^t
    // we can just transpose the derivatives of the direct matrix
    cvm::vector1d< cvm::rvector, 4 >      &dq_1 = atoms.rot.dQ0_1[ia];

    g11 = 2.0 * quat0[1]*dq_1[1];
    g22 = 2.0 * quat0[2]*dq_1[2];
    g33 = 2.0 * quat0[3]*dq_1[3];
    g01 = quat0[0]*dq_1[1] + quat0[1]*dq_1[0];
    g02 = quat0[0]*dq_1[2] + quat0[2]*dq_1[0];
    g03 = quat0[0]*dq_1[3] + quat0[3]*dq_1[0];
    g12 = quat0[1]*dq_1[2] + quat0[2]*dq_1[1];
    g13 = quat0[1]*dq_1[3] + quat0[3]*dq_1[1];
    g23 = quat0[2]*dq_1[3] + quat0[3]*dq_1[2];

    // Gradient of the inverse rotation matrix wrt current Cartesian position
    // (transpose of the gradient of the direct rotation)
    grad_rot_mat[0][0] = -2.0 * (g22 + g33); 
    grad_rot_mat[0][1] =  2.0 * (g12 + g03); 
    grad_rot_mat[0][2] =  2.0 * (g13 - g02); 
    grad_rot_mat[1][0] =  2.0 * (g12 - g03); 
    grad_rot_mat[1][1] = -2.0 * (g11 + g33); 
    grad_rot_mat[1][2] =  2.0 * (g01 + g23); 
    grad_rot_mat[2][0] =  2.0 * (g02 + g13); 
    grad_rot_mat[2][1] =  2.0 * (g23 - g01); 
    grad_rot_mat[2][2] = -2.0 * (g11 + g22); 

    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        sum += grad_rot_mat[i][j][i] * eigenvec[ia][j];
      }
    }
  }

  jd.real_value = sum * sqrt (eigenvec_invnorm2); 
}
