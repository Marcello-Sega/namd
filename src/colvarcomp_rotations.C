#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"


// XX todo implement rotation_angle


colvar::orientation::orientation (std::string const &conf)
  : cvc (conf)
{
  function_type = "orientation";
  parse_group (conf, "atoms", atoms);
  x.type (colvarvalue::type_quaternion);

  ref_pos.reserve (atoms.size());

  {
    std::string ref_pos_str = "";
    if (key_lookup (conf, "refPositions", ref_pos_str)) {

      std::istringstream rcis (ref_pos_str);
      for (size_t i = 0; i < atoms.size(); i++) {
        cvm::atom_pos x (0.0, 0.0, 0.0);
        if (! (rcis >> x)) {
          cvm::fatal_error ("Error: could not read all reference "
                            "coordinates for RMSD calculation.\n");
        }
        ref_pos.push_back (x);
      }
    } 
  }

  {
    std::string file_name;
    if (get_keyval (conf, "refPositionsFile", file_name)) {
      std::string file_col;
      get_keyval (conf, "refPositionsCol", file_col, std::string ("O"));
      ref_pos.resize (atoms.size());
      cvm::load_coords (file_name.c_str(), ref_pos, file_col);
    }
  }

  if (!ref_pos.size()) {
    cvm::fatal_error ("Error: must define a set of "
                      "reference coordinates.\n");
  } else {

    cvm::log ("Centering the reference coordinates: it is "
              "assumed that each atom is the closest "
              "periodic image to the center of mass.\n");
    cvm::rvector com (0.0, 0.0, 0.0);
    for (size_t i = 0; i < ref_pos.size(); i++) {
      com += atoms[i].mass * ref_pos[i];
    }
    com /= atoms.total_mass;
    if (cvm::debug())
      cvm::log ("colvar::orientation: ref_pos com: "+cvm::to_str (com)+"\n");
    for (size_t i = 0; i < ref_pos.size(); i++) {
      ref_pos[i] -= com;
    }
  }


  if (get_keyval (conf, "overlapQuaternion", ref_quat,
                  cvm::quaternion ())) {
    cvm::log ("Outputting quaternions with a positive scalar product to "+
              cvm::to_str (ref_quat)+".\n");
  }

  // initialize rot member data
  if (!atoms.noforce) {
    rot.dS_2.resize  (atoms.size(), cvm::matrix2d<cvm::rvector, 4, 4>());
    rot.dL0_2.resize (atoms.size(), cvm::rvector());
    rot.dQ0_2.resize (atoms.size(), cvm::vector1d<cvm::rvector, 4>());
  }
}

  
colvar::orientation::orientation()
  : cvc ()
{
  function_type = "orientation";
  x.type (colvarvalue::type_quaternion);
}



void colvar::orientation::calc_value()
{
  atoms.reset_atoms_data();
  atoms.read_positions();

  if (cvm::debug())
    cvm::log ("colvar::orientation: current com: "+
              cvm::to_str (atoms.center_of_mass())+"\n");

  if (!atoms.b_center)
    atoms.apply_translation (-1.0 * atoms.center_of_mass());
  //   if (cvm::debug()) {
  //     for (size_t i = 0; i < atoms.size(); i++) {
  //       cvm::log ("position "+cvm::to_str (i+1)+", Delta (current, ref): "+
  //                 cvm::to_str (atoms[i].pos - ref_pos[i])+"\n");
  //     }
  //   }

  rot.calc_optimal_rotation (ref_pos, atoms.positions(),
                             lambda0, quat0,
                             rot.dS_1,  rot.dS_2,
                             rot.dL0_1, rot.dL0_2,
                             rot.dQ0_1, rot.dQ0_2);

  if (quat0.inner (ref_quat) >= 0.0) {
    x.quaternion_value = quat0;
  } else {
    x.quaternion_value = -1.0 * quat0;
  }
}


void colvar::orientation::calc_gradients()
{
  // gradients have already been calculated (they are relatively cheap
  // with respect to the matrix diagonalization)


  if (b_debug_gradients) {

    colvarvalue dx_grad (x.type());

    if (atoms.old_pos.size()) {
      for (size_t i = 0; i < atoms.size(); i++) {
        cvm::rvector const &pold = atoms.old_pos[i];
        cvm::rvector const &p = atoms[i].pos;
        for (size_t iq = 0; iq < 4; iq++) {
          dx_grad.quaternion_value[iq] +=
            (p - pold) * rot.dQ0_2[i][iq];
        }
      }
    }

    // save for next step
    atoms.old_pos = atoms.positions();
    
    if (cvm::step_relative() == 0) {
      x_old.type (x.type());
    } else {
      colvarvalue const dx_real = 0.5 * this->dist2_lgrad (x, x_old);
      cvm::log ("Gradients check: difference in Dx = "+
                cvm::to_str (dx_grad-dx_real, cvm::cv_width, cvm::cv_prec)+"\n");
    }

    x_old = x;
  }
}


void colvar::orientation::apply_force (colvarvalue const &force)
{
  cvm::quaternion const &FQ = force.quaternion_value;

  if (!atoms.noforce) {
    for (size_t ia = 0; ia < atoms.size(); ia++) {
      for (size_t i = 0; i < 4; i++) {
        atoms[ia].apply_force (FQ[i] * rot.dQ0_2[ia][i]);
      }
    }
  }
}


colvar::rmsd::rmsd (std::string const &conf)
  : orientation (conf)
{
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  function_type = "rmsd";
  x.type (colvarvalue::type_scalar);

  ref_pos_sum2 = 0.0;
  for (size_t i = 0; i < ref_pos.size(); i++) {
    ref_pos_sum2 += ref_pos[i].norm2();
  }
}

  
colvar::rmsd::rmsd()
  : orientation()
{
  b_inverse_gradients = true;
  b_Jacobian_derivative = true;
  function_type = "rmsd";
  x.type (colvarvalue::type_scalar);
}


void colvar::rmsd::calc_value()
{
  atoms.reset_atoms_data();
  atoms.read_positions();

  if (cvm::debug())
    cvm::log ("colvar::rmsd: current com: "+
              cvm::to_str (atoms.center_of_mass())+"\n");

  if (!atoms.b_center)
    atoms.apply_translation (-1.0 * atoms.center_of_mass());
  //   if (cvm::debug()) {
  //     for (size_t i = 0; i < atoms.size(); i++) {
  //       cvm::log ("position "+cvm::to_str (i+1)+", Delta (current, ref): "+
  //                 cvm::to_str (atoms[i].pos - ref_pos[i])+"\n");
  //     }
  //   }

  rot.calc_optimal_rotation (ref_pos, atoms.positions(), 
                             lambda0, quat0, 
                             rot.dS_1,  rot.dS_2,
                             rot.dL0_1, rot.dL0_2,
                             rot.dQ0_1, rot.dQ0_2);

  cvm::real group_pos_sum2 = 0.0;
  for (size_t i = 0; i < atoms.size(); i++) {
    group_pos_sum2 += atoms[i].pos.norm2();
  }

  if (cvm::debug())
    cvm::log ("colvar::rmsd: ref_pos_sum2 = "+cvm::to_str (ref_pos_sum2, cvm::cv_width, cvm::cv_prec)+
              ", group_pos_sum2 = "+cvm::to_str (group_pos_sum2, cvm::cv_width, cvm::cv_prec)+
              ", lambda = "+cvm::to_str (lambda0, cvm::cv_width, cvm::cv_prec)+"\n");

  // value of the RMSD (Coutsias et al)
  cvm::real const MSD = 1.0/(cvm::real (atoms.size())) *
    ( group_pos_sum2 + ref_pos_sum2 - 2.0 * lambda0 );

  x.real_value = (MSD > 0.0) ? ::sqrt (MSD) : 0.0;

  if (cvm::debug()) {
    cvm::real MSDe = 0.0;
    for (size_t i = 0; i < atoms.size(); i++) {
      MSDe += (atoms[i].pos - quat0.rotate (ref_pos[i])).norm2();
    }
    MSDe /= cvm::real (atoms.size());
    cvm::log ("colvar::rmsd: explicit value = "+
              cvm::to_str (::sqrt (MSDe), cvm::cv_width, cvm::cv_prec)+"\n");
  }
}


void colvar::rmsd::calc_gradients()
{
  cvm::real const drdx2 = (x.real_value > 0.0) ?
    0.5 * (1.0/x.real_value) * 1.0/(cvm::real (atoms.size())) :
    0.0;
  cvm::real const drdl  = (x.real_value > 0.0) ?
    0.5 * (1.0/x.real_value) * 1.0/(cvm::real (atoms.size())) * (-2.0) :
    0.0;

  for (size_t ia = 0; ia < atoms.size(); ia++) {
    atoms[ia].grad = (drdx2 * (2.0*atoms[ia].pos) +
                      drdl  * rot.dL0_2[ia]);
  }

  if (b_debug_gradients) {

    colvarvalue dx_grad (x.type());
    dx_grad += fdiff_change (atoms);
    
    if (cvm::step_relative() == 0) {
      x_old.type (x.type());
    } else {
      colvarvalue const dx_real = 0.5 * this->dist2_lgrad (x, x_old);
      cvm::log ("Gradients check: difference in Dx = "+
                cvm::to_str (dx_grad-dx_real, cvm::cv_width, cvm::cv_prec)+"\n");
    }

    x_old = x;
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
    
  // Note: gradient norm is 1/sqrt(N_atoms)
          
  for (size_t ia = 0; ia < atoms.size(); ia++) {
    ft.real_value += atoms.size() * atoms[ia].grad *
      atoms[ia].system_force;
  }
}


void colvar::rmsd::calc_Jacobian_derivative()
{
  jd.real_value = x.real_value ? (3.0 * atoms.size() - 6.0) / x.real_value : 0.0;
}
