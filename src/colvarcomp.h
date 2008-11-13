#ifndef COLVARDEF_H
#define COLVARDEF_H

#include <fstream>
#include <cmath>


#include "colvarmodule.h"
#include "colvar.h"
#include "colvaratoms.h"


/// \brief Colvar component (base class); most implementations of
/// \link cvc \endlink utilize one or more \link
/// colvarmodule::atom \endlink or \link colvarmodule::atom_group
/// \endlink objects to access atoms.
///
/// A \link cvc \endlink object (or an object of a
/// cvc-derived class) specifies how to calculate a collective
/// variable, its gradients and other related physical quantities
/// which do not depend only on the numeric value (the \link colvar
/// \endlink class already serves this purpose).
///
/// No restriction is set to what kind of calculation a \link
/// cvc \endlink object performs (usually calculate an
/// analytical function of atomic coordinates).  The only constraint
/// is that the value calculated is implemented as a \link colvarvalue
/// \endlink object.  This serves to provide a unique way to calculate
/// scalar and non-scalar collective variables, and specify if and how
/// they can be combined together by the parent \link colvar \endlink
/// object.
///
/// If you wish to implement a new collective variable component, you
/// should write your own class by inheriting directly from \link
/// cvc \endlink, or one of its derived classes (for instance,
/// \link distance \endlink is frequently used as it can provide
/// useful data and function members for any variable based on two
/// atom groups).
/// 
/// The cvm::atom and cvm::atom_group classes are available to
/// transparently communicate with the simulation program.  However,
/// there is no need to use them as long as all the degrees of freedom
/// associated to the cvc are properly evolved by a simple call
/// to e.g. apply_force().

class colvar::cvc
  : public colvarparse
{
public:

  /// \brief The name of the object (helps to identify this
  /// cvc instance when debugging)
  std::string name;

  /// \brief Description of the type of collective variable
  /// 
  /// Normally this string is set by the parent \link colvar \endlink
  /// object within its constructor, when all \link cvc \endlink
  /// objects are initialized; therefore the main "config string"
  /// constructor does not need to define it.  If a \link cvc
  /// \endlink is initialized and/or a different constructor is used,
  /// this variable definition should be set within the constructor.
  std::string function_type;

  /// \brief Type of \link colvarvalue \endlink that this cvc
  /// provides
  colvarvalue::Type type() const;

  /// \brief Coefficient in the polynomial superposition (default: 1.0)
  cvm::real sup_coeff;
  /// \brief Exponent in the polynomial superposition (default: 1)
  int       sup_np;

  /// \brief Period of this cvc value, (default: 0.0, non periodic)
  cvm::real period;

  /// \brief Constructor
  ///
  /// At least one constructor which reads a string should be defined
  /// for every class inheriting from cvc \param conf Contents
  /// of the configuration file pertaining to this \link cvc
  /// \endlink
  cvc (std::string const &conf);

  /// \brief Within the constructor, make a group parse its own
  /// options from the provided configuration string
  void parse_group (std::string const &conf, 
                    char const *group_key,
                    cvm::atom_group &group,
                    bool optional = false);

  /// \brief Default constructor (used when \link cvc \endlink
  /// objects are declared within other ones)
  cvc();

  /// Destructor
  virtual ~cvc();

  /// \brief If true, calc_gradients() will calculate
  /// finite-difference gradients alongside the analytical ones and
  /// report their differences
  bool b_debug_gradients;

  /// \brief When b_debug_gradients is true, this function can be used
  /// to calculate the estimated change in the value using the change
  /// in the atomic coordinates times the atomic gradients
  colvarvalue fdiff_change (cvm::atom_group &group);

  /// \brief If this flag is false (default), inverse gradients
  /// (derivatives of atom coordinates with respect to x) are
  /// unavailable; it should be set to true by the constructor of each
  /// derived object capable of calculating them
  bool b_inverse_gradients;

  /// \brief If this flag is false (default), the Jacobian derivative
  /// (divergence of the inverse gradients) is unavailable; it should
  /// be set to true by the constructor of each derived object capable
  /// of calculating it
  bool b_Jacobian_derivative;

  /// \brief Calculate the variable
  virtual void calc_value() = 0;

  /// \brief Calculate the atomic gradients, to be reused later in
  /// order to apply forces
  virtual void calc_gradients() = 0;

  /// \brief Calculate the total force from the system using the
  /// inverse atomic gradients
  virtual void calc_force_invgrads();

  /// \brief Calculate the divergence of the inverse atomic gradients
  virtual void calc_Jacobian_derivative();


  /// \brief Return the previously calculated value
  virtual colvarvalue value() const;

  /// \brief Return the previously calculated system force
  virtual colvarvalue system_force() const;

  /// \brief Return the previously calculated divergence of the
  /// inverse atomic gradients
  virtual colvarvalue Jacobian_derivative() const;

  /// \brief Apply the collective variable force, by communicating the
  /// atomic forces to the simulation program (\b Note: the \link ft
  /// \endlink member is not altered by this function)
  ///
  /// Note: multiple calls to this function within the same simulation
  /// step will add the forces altogether \param cvforce The
  /// collective variable force, usually coming from the biases and
  /// eventually manipulated by the parent \link colvar \endlink
  /// object
  virtual void apply_force (colvarvalue const &cvforce) = 0;

  /// \brief Square distance between x1 and x2 (can be redefined to
  /// transparently implement constraints, symmetries and
  /// periodicities)
  /// 
  /// colvar::cvc::dist2() and the related functions are
  /// declared as "const" functions, but not "static", because
  /// additional parameters defining the metrics (e.g. the
  /// periodicity) may be specific to each colvar::cvc object.
  ///
  /// If symmetries or periodicities are present, the
  /// colvar::cvc::dist2() should be redefined to return the
  /// "closest distance" value and colvar::cvc::dist2_lgrad(),
  /// colvar::cvc::dist2_rgrad() to return its gradients.
  ///
  /// If constraints are present (and not already implemented by any
  /// of the \link colvarvalue \endlink types), the
  /// colvar::cvc::dist2_lgrad() and
  /// colvar::cvc::dist2_rgrad() functions should be redefined
  /// to provide a gradient which is compatible with the constraint,
  /// i.e. already deprived of its component normal to the constraint
  /// hypersurface.
  /// 
  /// Finally, another useful application, if you are performing very
  /// many operations with these functions, could be to override the
  /// \link colvarvalue \endlink member functions and access directly
  /// its member data.  For instance: to define dist2(x1,x2) as
  /// (x2.real_value-x1.real_value)*(x2.real_value-x1.real_value) in
  /// case of a scalar \link colvarvalue \endlink type.
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;

  /// \brief Gradient (with respect to x1) of the square distance (can
  /// be redefined to transparently implement constraints, symmetries
  /// and periodicities)
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;

  /// \brief Gradient (with respect to x2) of the square distance (can
  /// be redefined to transparently implement constraints, symmetries
  /// and periodicities)
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;

  /// \brief Return a positive number if x2>x1, zero if x2==x1,
  /// negative otherwise (can be redefined to transparently implement
  /// constraints, symmetries and periodicities) \b Note: \b it \b
  /// only \b works \b with \b scalar \b variables, otherwise raises
  /// an error
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;

protected:

  /// \brief Cached value
  colvarvalue x;

  /// \brief Value at the previous step
  colvarvalue x_old;

  /// \brief Calculated system force (\b Note: this is calculated from
  /// the total atomic forces read from the program, subtracting fromt
  /// the "internal" forces of the system the "external" forces from
  /// the colvar biases)
  colvarvalue ft;

  /// \brief Calculated Jacobian derivative (divergence of the inverse
  /// gradients): serves to calculate the phase space correction
  colvarvalue jd;
};




inline colvarvalue::Type colvar::cvc::type() const
{
  return x.type();
}

inline colvarvalue colvar::cvc::value() const
{
  return x;
}

inline colvarvalue colvar::cvc::system_force() const
{
  return ft;
}

inline colvarvalue colvar::cvc::Jacobian_derivative() const
{
  return jd;
}


inline cvm::real colvar::cvc::dist2 (colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  return x1.dist2 (x2);
}

inline colvarvalue colvar::cvc::dist2_lgrad (colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return x1.dist2_grad (x2);
}

inline colvarvalue colvar::cvc::dist2_rgrad (colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return x2.dist2_grad (x1);
}

inline cvm::real colvar::cvc::compare (colvarvalue const &x1,
                                       colvarvalue const &x2) const
{
  if (this->type() == colvarvalue::type_scalar) {
    return cvm::real (x1 - x2);
  } else {
    cvm::fatal_error ("Error: you requested an operation which requires "
                      "comparison between two non-scalar values.\n");
    return 0.0;
  }
}



/// \brief Colvar component: distance between the centers of mass of
/// two groups (colvarvalue::type_scalar type, range [0:*))
///
/// This class also serves as the template for many collective
/// variables with two atom groups: in this case, the
/// distance::distance() constructor should be called on the same
/// configuration string, to make the same member data and functions
/// available to the derived object
class colvar::distance
  : public colvar::cvc
{
protected:
  /// First atom group
  cvm::atom_group  group1;
  /// Second atom group
  cvm::atom_group  group2;
  /// Vector distance, cached to be recycled
  cvm::rvector     dist_v;
  /// Compute system force on first site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  bool b_1site_force;
public:
  distance (std::string const &conf);
  distance();
  virtual inline ~distance() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: distance vector between centers of mass
/// of two groups (\link colvarvalue::type_vector \endlink type,
/// range (-*:*)x(-*:*)x(-*:*))
class colvar::distance_vec
  : public colvar::distance
{
public:
  distance_vec (std::string const &conf);
  distance_vec();
  virtual inline ~distance_vec() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force (colvarvalue const &force);
  /// Redefined to handle the box periodicity
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  /// Redefined to handle the box periodicity
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  /// Redefined to handle the box periodicity
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  /// Redefined to handle the box periodicity
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: distance unit vector (direction) between
/// centers of mass of two groups (colvarvalue::type_unitvector type,
/// range [-1:1]x[-1:1]x[-1:1])
class colvar::distance_dir
  : public colvar::distance
{
public:
  distance_dir (std::string const &conf);
  distance_dir();
  virtual inline ~distance_dir() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force (colvarvalue const &force);
//   virtual cvm::real dist2 (colvarvalue const &x1,
//                            colvarvalue const &x2) const;
//   virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
//                                    colvarvalue const &x2) const;
//   virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
//                                    colvarvalue const &x2) const;
//   virtual cvm::real compare (colvarvalue const &x1,
//                              colvarvalue const &x2) const;
};


/// Colvar component: projection of the distance vector along
/// an axis (colvarvalue::type_scalar type, range (-*:*))
class colvar::distance_z
  : public colvar::cvc
{
protected:
  /// Main atom group
  cvm::atom_group  main;
  /// Reference atom group
  cvm::atom_group  ref1;
  /// Optional, second ref atom group
  cvm::atom_group  ref2;
  /// Compute system force on one site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  bool b_1site_force;
  /// Vector on which the distance vector is projected
  cvm::rvector axis;
  /// Norm of the axis
  cvm::real axis_norm;
  /// Vector distance, cached to be recycled
  cvm::rvector     dist_v;
  /// Flag: using a fixed axis vector?
  bool fixed_axis;
public:
  distance_z (std::string const &conf);
  distance_z();
  virtual inline ~distance_z() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: projection of the distance vector along
/// a fixed axis (colvarvalue::type_scalar type, range (-*:*))
class colvar::distance_xy
  : public colvar::distance_z
{
protected:
  /// Components of the distance vector orthogonal to the axis
  cvm::rvector dist_v_ortho;
  /// Vector distances
  cvm::rvector v12, v13;
public:
  distance_xy (std::string const &conf);
  distance_xy();
  virtual inline ~distance_xy() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: Radius of gyration of an atom group
/// (colvarvalue::type_scalar type, range [0:*))
class colvar::gyration
  : public colvar::cvc
{
protected:
  /// Atoms involved
  cvm::atom_group atoms;
public:
  /// Constructor
  gyration (std::string const &conf);
  gyration();
  virtual inline ~gyration() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: angle between the centers of mass of
/// three groups (colvarvalue::type_scalar type, range [0:PI])
class colvar::angle
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group group1;
  /// Atom group
  cvm::atom_group group2;
  /// Atom group
  cvm::atom_group group3;

  /// Inter site vectors
  cvm::rvector r21, r23;
  /// Inter site vector norms
  cvm::real r21l, r23l;

public:

  /// Initialize by parsing the configuration
  angle (std::string const &conf);
  /// \brief Initialize the three groups after three atoms
  angle (cvm::atom const &a1, cvm::atom const &a2, cvm::atom const &a3);
  angle();
  virtual inline ~angle() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: dihedral between the centers of mass of
/// four groups (colvarvalue::type_scalar type, range [-PI:PI])
class colvar::dihedral
  : public colvar::cvc
{
protected:

  /// Atom group
  cvm::atom_group group1;
  /// Atom group
  cvm::atom_group group2;
  /// Atom group
  cvm::atom_group group3;
  /// Atom group
  cvm::atom_group group4;
  /// Inter site vectors
  cvm::rvector r12, r23, r34;

  /// Compute system force on first site only to avoid unwanted
  /// coupling to other colvars (see e.g. Ciccotti et al., 2005)
  bool b_1site_force;
  
public:

  /// Initialize by parsing the configuration
  dihedral (std::string  const &conf);
  /// \brief Initialize the four groups after four atoms
  dihedral (cvm::atom const &a1, cvm::atom const &a2, cvm::atom const &a3, cvm::atom const &a4);
  dihedral();
  virtual inline ~dihedral() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force (colvarvalue const &force);

  /// Redefined to handle the 2*PI periodicity
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  /// Redefined to handle the 2*PI periodicity
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: coordination number between two groups
/// (colvarvalue::type_scalar type, range [0:N1*N2])
class colvar::coordnum
  : public colvar::distance
{
protected:
  /// \brief "Cutoff" for isotropic calculation (default)
  cvm::real     r0;
  /// \brief "Cutoff vector" for anisotropic calculation
  cvm::rvector  r0_vec;
  /// \brief Wheter dist/r0 or \vec{dist}*\vec{1/r0_vec} should ne be
  /// used
  bool b_anisotropic;
  /// Integer exponent of the function numerator
  int en;
  /// Integer exponent of the function denominator
  int ed;
  /// \brief If true, group2 will be treated as a single atom
  /// (default: loop over all pairs of atoms in group1 and group2)
  bool b_group2_center_only;
public:
  /// Constructor
  coordnum (std::string const &conf);
  coordnum();
  virtual inline ~coordnum() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force (colvarvalue const &force);
  template<bool b_gradients>
  /// \brief Calculate a coordination number through the function
  /// (1-x**n)/(1-x**m), x = |A1-A2|/r0 \param r0 "cutoff" for the
  /// coordination number \param exp_num \i n exponent \param exp_den
  /// \i m exponent \param A1 atom \param A2 atom
  static cvm::real switching_function (cvm::real const &r0,
                                       int const &exp_num, int const &exp_den,
                                       cvm::atom &A1, cvm::atom &A2);
  
  template<bool b_gradients>
  /// \brief Calculate a coordination number through the function
  /// (1-x**n)/(1-x**m), x = |(A1-A2)*(r0_vec)^-|1 \param r0_vec
  /// vector of different cutoffs in the three directions \param
  /// exp_num \i n exponent \param exp_den \i m exponent \param A1
  /// atom \param A2 atom
  static cvm::real switching_function (cvm::rvector const &r0_vec,
                                       int const &exp_num, int const &exp_den,
                                       cvm::atom &A1, cvm::atom &A2);

  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: hydrogen bond, defined as the product of
/// a colvar::coordnum and 1/2*(1-cos((180-ang)/ang_tol))
/// (colvarvalue::type_scalar type, range [0:1])
class colvar::h_bond 
  : public colvar::cvc
{
protected:
  /// Atoms involved in the component
  cvm::atom     acceptor, donor;
  /// \brief "Cutoff" distance between acceptor and donor
  cvm::real     r0;
  /// Integer exponent of the function numerator
  int en;
  /// Integer exponent of the function denominator
  int ed;
public:
  h_bond (std::string const &conf);
  /// Constructor for atoms already allocated
  h_bond (cvm::atom const &acceptor,
          cvm::atom const &donor,
          cvm::real r0, int en, int ed);
  h_bond();
  virtual inline ~h_bond() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force (colvarvalue const &force);

  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: alpha helix content of a contiguous
/// segment of 5 or more residues, implemented as a sum of phi/psi
/// dihedral angles and hydrogen bonds (colvarvalue::type_scalar type,
/// range [0:1])
class colvar::alpha_dihedrals
  : public colvar::cvc
{
protected:

  /// Alpha-helical reference phi value
  cvm::real phi_ref;

  /// Alpha-helical reference psi value
  cvm::real psi_ref;

  /// List of phi dihedral angles
  std::vector<dihedral *> phi;

  /// List of psi dihedral angles
  std::vector<dihedral *> psi;

  /// List of hydrogen bonds
  std::vector<h_bond *>   hb;

public:

  alpha_dihedrals (std::string const &conf);
  alpha_dihedrals();
  virtual inline ~alpha_dihedrals() {}
  virtual void calc_value();
  virtual void calc_gradients(); 
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: alpha helix content of a contiguous
/// segment of 5 or more residues, implemented as a sum of Ca-Ca-Ca
/// angles and hydrogen bonds (colvarvalue::type_scalar type, range
/// [0:1])
class colvar::alpha_angles
  : public colvar::cvc
{
protected:

  /// Reference Calpha-Calpha angle (default: 88 degrees)
  cvm::real theta_ref;

  /// Tolerance on the Calpha-Calpha angle
  cvm::real theta_tol;

  /// List of Calpha-Calpha angles
  std::vector<angle *> theta;

  /// List of hydrogen bonds
  std::vector<h_bond *>   hb;

  /// Contribution of the hb terms 
  cvm::real hb_coeff;

public:

  alpha_angles (std::string const &conf);
  alpha_angles();
  virtual inline ~alpha_angles() {}
  void calc_value();
  void calc_gradients();
  void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: orientation in space of an atom group,
/// with respect to a set of reference coordinates
/// (colvarvalue::type_quaternion type, range
/// [-1:1]x[-1:1]x[-1:1]x[-1:1])
class colvar::orientation
  : public colvar::cvc
{
protected:

  cvm::atom_group            atoms;
  std::vector<cvm::atom_pos> ref_pos;

  cvm::rotation              rot;
  cvm::real                  lambda0;
  cvm::quaternion            quat0;
  cvm::quaternion            ref_quat;

public:

  orientation (std::string const &conf);
  orientation();
  virtual inline ~orientation() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};


/// \brief Colvar component: root mean square deviation (RMSD) of a
/// group with respect to a set of reference coordinates; uses \link
/// colvar::orientation \endlink to calculate the rotation matrix
/// (colvarvalue::type_scalar type, range [0:*))
class colvar::rmsd
  : public colvar::orientation
{
protected:

  /// Sum of the squares of ref_coords
  cvm::real                  ref_pos_sum2;

public:

  /// Constructor
  rmsd (std::string const &conf);
  rmsd();
  virtual inline ~rmsd() {}
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void calc_force_invgrads();
  virtual void calc_Jacobian_derivative();
  virtual void apply_force (colvarvalue const &force);
  virtual cvm::real dist2 (colvarvalue const &x1,
                           colvarvalue const &x2) const;
  virtual colvarvalue dist2_lgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual colvarvalue dist2_rgrad (colvarvalue const &x1,
                                   colvarvalue const &x2) const;
  virtual cvm::real compare (colvarvalue const &x1,
                             colvarvalue const &x2) const;
};



// metrics functions for cvc implementations with a periodicity

inline cvm::real colvar::dihedral::dist2 (colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  // This version is asymptotically accurate for small distances,
  // other choices are possible
  cvm::real const diff = x1.real_value - x2.real_value;
  return (180.0*180.0/M_PI/M_PI) * 2.0 * (1.0 - ::cos (diff * (M_PI/180.00)));
}

inline colvarvalue colvar::dihedral::dist2_lgrad (colvarvalue const &x1,
                                                  colvarvalue const &x2) const
{
  cvm::real const diff = x1.real_value - x2.real_value;
  return colvarvalue ( (180.0/M_PI) * 2.0 * ::sin ( diff * (M_PI/180.0)));
}

inline colvarvalue colvar::dihedral::dist2_rgrad (colvarvalue const &x1,
                                                  colvarvalue const &x2) const
{
  cvm::real const diff = x1.real_value - x2.real_value;
  return colvarvalue ( (180.0/M_PI) * (-2.0) * ::sin ( diff * (M_PI/180.0)));
}

inline cvm::real colvar::dihedral::compare (colvarvalue const &x1,
                                            colvarvalue const &x2) const
{
  return dist2_lgrad (x1, x2);
}


// distance between three dimensional vectors

inline cvm::real colvar::distance_vec::dist2 (colvarvalue const &x1,
                                              colvarvalue const &x2) const
{
  return cvm::position_dist2 (x1, x2);
}

inline colvarvalue colvar::distance_vec::dist2_lgrad (colvarvalue const &x1,
                                                      colvarvalue const &x2) const
{
  return cvm::position_dist2_lgrad (x1, x2);
}

inline colvarvalue colvar::distance_vec::dist2_rgrad (colvarvalue const &x1,
                                                      colvarvalue const &x2) const
{
  return cvm::position_dist2_lgrad (x2, x1);
}

inline cvm::real colvar::distance_vec::compare (colvarvalue const &x1,
                                                colvarvalue const &x2) const
{
  cvm::fatal_error ("Error: cannot compare() two distance vectors.\n");
  return 0.0;
}


inline cvm::real colvar::distance_z::dist2 (colvarvalue const &x1,
                                            colvarvalue const &x2) const
{
  cvm::atom_pos const x1_vec (cvm::rvector (x1.real_value * axis));
  cvm::atom_pos const x2_vec (cvm::rvector (x2.real_value * axis));
  return cvm::position_dist2 (x1_vec, x2_vec);
}

inline colvarvalue colvar::distance_z::dist2_lgrad (colvarvalue const &x1,
                                                    colvarvalue const &x2) const
{
  cvm::atom_pos const x1_vec (cvm::rvector (x1.real_value * axis));
  cvm::atom_pos const x2_vec (cvm::rvector (x2.real_value * axis));
  return axis * cvm::position_dist2_lgrad (x1_vec, x2_vec);
}

inline colvarvalue colvar::distance_z::dist2_rgrad (colvarvalue const &x1,
                                                    colvarvalue const &x2) const
{
  cvm::atom_pos const x1_vec (cvm::rvector (x1.real_value * axis));
  cvm::atom_pos const x2_vec (cvm::rvector (x2.real_value * axis));
  return axis * cvm::position_dist2_lgrad (x2_vec, x1_vec);
}

inline cvm::real colvar::distance_z::compare (colvarvalue const &x1,
                                              colvarvalue const &x2) const
{
  return dist2_lgrad (x1, x2);
}



// simple definitions, they are here only for optimization purposes


// definitions assuming the scalar type

#define simple_scalar_dist_functions(TYPE)                              \
                                                                        \
  inline cvm::real colvar::TYPE::dist2 (colvarvalue const &x1,          \
                                        colvarvalue const &x2) const    \
  {                                                                     \
    return ::pow (x1.real_value - x2.real_value, int (2));              \
  }                                                                     \
                                                                        \
  inline colvarvalue colvar::TYPE::dist2_lgrad (colvarvalue const &x1,  \
                                                colvarvalue const &x2) const \
  {                                                                     \
    return 2.0 * (x1.real_value - x2.real_value);                       \
  }                                                                     \
                                                                        \
  inline colvarvalue colvar::TYPE::dist2_rgrad (colvarvalue const &x1,  \
                                                colvarvalue const &x2) const \
  {                                                                     \
    return this->dist2_lgrad (x2, x1);                                  \
  }                                                                     \
                                                                        \
  inline cvm::real colvar::TYPE::compare (colvarvalue const &x1,        \
                                          colvarvalue const &x2) const  \
  {                                                                     \
    return this->dist2_lgrad (x1, x2);                                  \
  }                                                                     \
                                                                        \

  simple_scalar_dist_functions (distance)
  simple_scalar_dist_functions (distance_xy)
  simple_scalar_dist_functions (angle)
  simple_scalar_dist_functions (coordnum)
  simple_scalar_dist_functions (h_bond)
  simple_scalar_dist_functions (rmsd)
  simple_scalar_dist_functions (gyration)
  simple_scalar_dist_functions (alpha_dihedrals)
  simple_scalar_dist_functions (alpha_angles)

// generic definitions


// distance between three dimensional vectors

inline cvm::real colvar::orientation::dist2 (colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return x1.quaternion_value.dist2 (x2);
}

inline colvarvalue colvar::orientation::dist2_lgrad (colvarvalue const &x1,
                                                     colvarvalue const &x2) const
{
  return x1.quaternion_value.dist2_grad (x2);
}

inline colvarvalue colvar::orientation::dist2_rgrad (colvarvalue const &x1,
                                                     colvarvalue const &x2) const
{
  return x2.quaternion_value.dist2_grad (x1);
}

inline cvm::real colvar::orientation::compare (colvarvalue const &x1,
                                                colvarvalue const &x2) const
{
  cvm::fatal_error ("Error: cannot compare() two orientation vectors.\n");
  return 0.0;
}



#endif


// Emacs
// Local Variables:
// mode: C++
// End: