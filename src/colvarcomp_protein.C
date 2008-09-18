#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"


//////////////////////////////////////////////////////////////////////
// alpha component
//////////////////////////////////////////////////////////////////////

colvar::alpha_angles::alpha_angles (std::string const &conf)
  : cvc (conf)
{
  if (cvm::debug())
    cvm::log ("Initializing alpha_angles object.\n");

  function_type = "alpha_angles";
  x.type (colvarvalue::type_scalar);

  std::vector<int> residues;
  {
    std::string residues_conf = "";
    key_lookup (conf, "residueRange", residues_conf);
    if (residues_conf.size()) {
      std::istringstream is (residues_conf);
      int initial, final;
      char dash;
      if ( (is >> initial) && (initial > 0) &&
           (is >> dash) && (dash == '-') &&
           (is >> final) && (final > 0) ) {
        for (int rnum = initial; rnum <= final; rnum++) {
          residues.push_back (rnum);
        }
      }
    } else {
      cvm::fatal_error ("Error: no residues defined in \"residueRange\".\n");
    }
  }

  if (residues.size() < 5) {
    cvm::fatal_error ("Error: not enough residues defined in \"residueRange\".\n");
  }

  std::string segment_id;
  get_keyval (conf, "psfSegID", segment_id, std::string ("MAIN"));

  std::string const &sid    = segment_id;
  std::vector<int> const &r = residues;


  get_keyval (conf, "hBondCoeff", hb_coeff, 0.5);
  if ( (hb_coeff < 0.0) || (hb_coeff > 1.0) ) {
    cvm::fatal_error ("Error: hBondCoeff must be defined between 0 and 1.\n");
  }


  get_keyval (conf, "angleRef", theta_ref, 88.0, parse_silent);
  get_keyval (conf, "angleTol", theta_tol, 15.0, parse_silent);

  if (hb_coeff < 1.0) {

    for (size_t i = 0; i < residues.size()-2; i++) {
      theta.push_back (new colvar::angle (cvm::atom (r[i  ], "CA", sid),
                                          cvm::atom (r[i+1], "CA", sid),
                                          cvm::atom (r[i+2], "CA", sid)));
    }

  } else {
    cvm::log ("The hBondCoeff specified will disable the Calpha-Calpha-Calpha angle terms.\n");
  }

  {
    cvm::real r0;
    size_t en, ed;
    get_keyval (conf, "hBondCutoff",   r0, (3.3 * cvm::unit_angstrom()));
    get_keyval (conf, "hBondExpNumer", en, 6);
    get_keyval (conf, "hBondExpDenom", ed, 8);

    if (hb_coeff > 0.0) {

      for (size_t i = 0; i < residues.size()-4; i++) {
        hb.push_back (new colvar::h_bond (cvm::atom (r[i  ], "O",  sid),
                                          cvm::atom (r[i+4], "N",  sid),
                                          r0, en, ed));
      }

    } else {
      cvm::log ("The hBondCoeff specified will disable the hydrogen bond terms.\n");
    }
  }

  if (cvm::debug())
    cvm::log ("Done initializing alpha_angles object.\n");
}


colvar::alpha_angles::alpha_angles()
  : cvc ()
{
  function_type = "alpha_angles";
  x.type (colvarvalue::type_scalar);
}


void colvar::alpha_angles::calc_value()
{
  x.real_value = 0.0;

  if (theta.size()) {

    cvm::real const theta_norm = 
      (1.0-hb_coeff) / cvm::real (theta.size());

    for (size_t i = 0; i < theta.size(); i++) {

      (theta[i])->calc_value();

      cvm::real const t = ((theta[i])->value().real_value-theta_ref)/theta_tol;
      cvm::real const f = ( (1.0 - ::pow (t, (int) 2)) /
                            (1.0 - ::pow (t, (int) 4)) );

      x.real_value += theta_norm * f;

      if (cvm::debug())
        cvm::log ("Calpha-Calpha angle no. "+cvm::to_str (i+1)+" in \""+
                  this->name+"\" has a value of "+
                  (cvm::to_str ((theta[i])->value().real_value))+
                  " degrees, f = "+cvm::to_str (f)+".\n");
    }
  }

  if (hb.size()) {

    cvm::real const hb_norm =
      hb_coeff / cvm::real (hb.size());

    for (size_t i = 0; i < hb.size(); i++) {
      (hb[i])->calc_value();
      x.real_value += hb_norm * (hb[i])->value().real_value;
      if (cvm::debug())
        cvm::log ("Hydrogen bond no. "+cvm::to_str (i+1)+" in \""+
                  this->name+"\" has a value of "+
                  (cvm::to_str ((hb[i])->value().real_value))+".\n");
    }
  }
}


void colvar::alpha_angles::calc_gradients()
{
  for (size_t i = 0; i < theta.size(); i++) 
    (theta[i])->calc_gradients();

  for (size_t i = 0; i < hb.size(); i++)
    (hb[i])->calc_gradients();
}


void colvar::alpha_angles::apply_force (colvarvalue const &force)
{

  if (theta.size()) {

    cvm::real const theta_norm = 
      (1.0-hb_coeff) / cvm::real (theta.size());
    
    for (size_t i = 0; i < theta.size(); i++) {

      cvm::real const t = ((theta[i])->value().real_value-theta_ref)/theta_tol;
      cvm::real const f = ( (1.0 - ::pow (t, (int) 2)) /
                            (1.0 - ::pow (t, (int) 4)) );

      cvm::real const dfdt =
        1.0/(1.0 - ::pow (t, (int) 4)) * 
        ( (-2.0 * t) + (-1.0*f)*(-4.0 * ::pow (t, (int) 3)) );

      (theta[i])->apply_force (theta_norm * 
                               dfdt * (1.0/theta_tol) *
                               force.real_value );
    }
  }

  if (hb.size()) {

    cvm::real const hb_norm =
      hb_coeff / cvm::real (hb.size());

    for (size_t i = 0; i < hb.size(); i++) {
      (hb[i])->apply_force (0.5 * hb_norm * force.real_value);
    }
  }
}



//////////////////////////////////////////////////////////////////////
// alphaDihedrals component
//////////////////////////////////////////////////////////////////////


// Restraint function for the dihedral
inline cvm::real dih_func (cvm::real const &w,
                           cvm::real const &w0)
{
  return 0.5*(1.0 + ::cos (w-w0));
}

// Derivative
inline cvm::real dih_deriv (cvm::real const &w,
                            cvm::real const &w0)
{
  return (-0.5 * ::sin (w-w0));
}


colvar::alpha_dihedrals::alpha_dihedrals (std::string const &conf)
  : cvc (conf)
{
  if (cvm::debug())
    cvm::log ("Initializing alpha_dihedrals object.\n");

  function_type = "alpha_dihedrals";
  x.type (colvarvalue::type_scalar);

  std::vector<int> residues;

  get_keyval (conf, "residues", residues, std::vector<int>());

  phi.reserve (residues.size());
  psi.reserve (residues.size());
  hb.reserve  (residues.size());

  get_keyval (conf, "phi_ref", phi_ref, -57.8, parse_silent);
  get_keyval (conf, "psi_ref", psi_ref, -47.0, parse_silent);

  if (residues.size() < 5) {
    cvm::fatal_error ("Error: not enough residues defined.\n");
  }

  std::string segment_id;
  get_keyval (conf, "segment_id", segment_id, std::string ("MAIN"));

  std::string const &sid    = segment_id;
  std::vector<int> const &r = residues;
  for (size_t i = 0; i < residues.size()-1; i++) {

    phi.push_back (new colvar::dihedral (cvm::atom (r[i  ], "C",  sid),
                                         cvm::atom (r[i+1], "N",  sid),
                                         cvm::atom (r[i+1], "CA", sid),
                                         cvm::atom (r[i+1], "C",  sid)));

    psi.push_back (new colvar::dihedral (cvm::atom (r[i  ], "N",  sid),
                                         cvm::atom (r[i  ], "CA", sid),
                                         cvm::atom (r[i  ], "C",  sid),
                                         cvm::atom (r[i+1], "N",  sid)));
  }

  for (size_t i = 0; i < residues.size()-4; i++) {
    hb.push_back (new colvar::h_bond (cvm::atom (r[i  ], "O",  sid),
                                      cvm::atom (r[i+4], "N",  sid),
                                      3.3 * cvm::unit_angstrom(),
                                      6, 8));
  }

  if (cvm::debug())
    cvm::log ("Done initializing alpha_dihedrals object.\n");
}


colvar::alpha_dihedrals::alpha_dihedrals()
  : cvc ()
{
  function_type = "alpha_dihedrals";
  x.type (colvarvalue::type_scalar);
}


void colvar::alpha_dihedrals::calc_value()
{
  x.real_value = 0.0;

  for (size_t i = 0; i < phi.size(); i++) {

    (phi[i])->calc_value();
    (psi[i])->calc_value();

    x.real_value += 0.5 *
      dih_func (((phi[i])->value()).real_value, phi_ref) *
      dih_func (((psi[i])->value()).real_value, psi_ref);

    if (cvm::debug())
      cvm::log ("Phi dihedral no. "+cvm::to_str (i+1)+" in \""+
                this->name+"\" has a value of "+
                (cvm::to_str ((phi[i])->value().real_value))+
                " degrees.\n");

    if (cvm::debug())
      cvm::log ("Psi dihedral no. "+cvm::to_str (i+1)+" in \""+
                this->name+"\" has a value of "+
                (cvm::to_str ((psi[i])->value().real_value))+
                " degrees.\n");
  }

  for (size_t i = 0; i < hb.size(); i++) {
    (hb[i])->calc_value();
    x.real_value += 0.5 * (hb[i])->value().real_value;
    if (cvm::debug())
      cvm::log ("Hydrogen bond no. "+cvm::to_str (i+1)+" in \""+
                this->name+"\" has a value of "+
                (cvm::to_str ((hb[i])->value().real_value))+".\n");
  }
}


void colvar::alpha_dihedrals::calc_gradients()
{
  for (size_t i = 0; i < phi.size(); i++) {
    (phi[i])->calc_gradients();
    (psi[i])->calc_gradients();
  }
  for (size_t i = 0; i < hb.size(); i++) {
    (hb[i])->calc_gradients();
  }
}


void colvar::alpha_dihedrals::apply_force (colvarvalue const &force)
{
  for (size_t i = 0; i < phi.size(); i++) {

    (phi[i])->apply_force ( 0.5 *
                            dih_func (((psi[i])->value()).real_value, psi_ref) *
                            dih_deriv (((phi[i])->value()).real_value, phi_ref) *
                            force.real_value );

    (psi[i])->apply_force ( 0.5 *
                            dih_func (((phi[i])->value()).real_value, phi_ref) *
                            dih_deriv (((psi[i])->value()).real_value, psi_ref) *
                            force.real_value );
  }

  for (size_t i = 0; i < hb.size(); i++) {
    (hb[i])->apply_force (0.5 * force.real_value);
  }
}



