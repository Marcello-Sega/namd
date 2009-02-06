#ifndef COLVARTYPES_H
#define COLVARTYPES_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923
#endif

// ----------------------------------------------------------------------
/// Linear algebra functions and data types used in the collective
/// variables implemented so far
// ----------------------------------------------------------------------


/// 1-dimensional vector of real numbers with three components
class colvarmodule::rvector {

public:

  cvm::real x, y, z;
     
  inline rvector()
    : x (0.0), y (0.0), z (0.0)
  {}

  inline rvector (cvm::real const &x_i,
                  cvm::real const &y_i,
                  cvm::real const &z_i)
    : x (x_i), y (y_i), z (z_i)
  {}

  inline rvector (cvm::real v)
    : x (v), y (v), z (v)
  {}

  inline cvm::real & operator [] (int const &i) {
    return (i == 0) ? x : (i == 1) ? y : (i == 2) ? z : x;
  }

  inline cvm::real const & operator [] (int const &i) const {
    return (i == 0) ? x : (i == 1) ? y : (i == 2) ? z : x;
  }


  inline cvm::rvector & operator = (cvm::real const &v) 
  {
    x = v;
    y = v;
    z = v;
    return *this;
  }

  inline void operator += (cvm::rvector const &v) 
  {
    x += v.x;
    y += v.y;
    z += v.z;
  }

  inline void operator -= (cvm::rvector const &v) 
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
  }

  inline void operator *= (cvm::real const &v) 
  {
    x *= v;
    y *= v;
    z *= v;
  }

  inline void operator /= (cvm::real const& v) 
  {
    x /= v;
    y /= v;
    z /= v;
  }

  inline cvm::real norm2() const
  {
    return (x*x + y*y + z*z);
  }

  inline cvm::real norm() const
  {
    return ::sqrt (this->norm2());
  }

  inline cvm::rvector unit() const
  {
    return cvm::rvector (x, y, z)/this->norm();
  }

  static inline size_t output_width (size_t const &real_width)
  {
    return 3*real_width + 10;
  }


  static inline cvm::rvector outer (cvm::rvector const &v1, cvm::rvector const &v2) 
  {
    return cvm::rvector ( v1.y*v2.z - v2.y*v1.z,
                          -v1.x*v2.z + v2.x*v1.z,
                          v1.x*v2.y - v2.x*v1.y);
  }

  friend inline cvm::rvector operator - (cvm::rvector const &v) 
  {
    return cvm::rvector (-v.x, -v.y, -v.z);
  }

  friend inline int operator == (cvm::rvector const &v1, cvm::rvector const &v2) 
  {
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
  }

  friend inline int operator != (cvm::rvector const &v1, cvm::rvector const &v2) 
  {
    return (v1.x != v2.x) || (v1.y != v2.y) || (v1.z != v2.z);
  }

  friend inline cvm::rvector operator + (cvm::rvector const &v1, cvm::rvector const &v2) 
  {
    return cvm::rvector (v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
  }
  friend inline cvm::rvector operator - (cvm::rvector const &v1, cvm::rvector const &v2) 
  {
    return cvm::rvector (v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
  }

  friend inline cvm::real operator * (cvm::rvector const &v1, cvm::rvector const &v2) 
  {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }

  friend inline cvm::rvector operator * (cvm::real const &a, cvm::rvector const &v) 
  {
    return cvm::rvector (a*v.x, a*v.y, a*v.z);
  }

  friend inline cvm::rvector operator * (cvm::rvector const &v, cvm::real const &a) 
  {
    return cvm::rvector (a*v.x, a*v.y, a*v.z);
  }

  friend inline cvm::rvector operator / (cvm::rvector const &v, cvm::real const &a) 
  {
    return cvm::rvector (v.x/a, v.y/a, v.z/a);
  }
     

};



/// \brief 1-dimensional vector of real numbers with four components and
/// a quaternion algebra
class colvarmodule::quaternion {

public:

  cvm::real q0, q1, q2, q3;

  /// Constructor from a 3-d vector
  inline quaternion (cvm::real const &x, cvm::real const &y, cvm::real const &z)
    : q0 (0.0), q1 (x), q2 (y), q3 (z)
  {}

  /// Constructor component by component
  inline quaternion (cvm::real const qv[4])
    : q0 (qv[0]), q1 (qv[1]), q2 (qv[2]), q3 (qv[3])
  {}

  /// Constructor component by component
  inline quaternion (cvm::real const &q0i,
                     cvm::real const &q1i,
                     cvm::real const &q2i,
                     cvm::real const &q3i)
    : q0 (q0i), q1 (q1i), q2 (q2i), q3 (q3i)
  {}

  /// \brief Default constructor
  inline quaternion() { reset(); }

  /// \brief Set all components to zero (null quaternion)
  inline void reset() { q0 = q1 = q2 = q3 = 0.0; }

  /// \brief Set the q0 component to 1 and the others to 0 (quaternion
  /// representing no rotation)
  inline void reset_rotation() { q0 = 1.0; q1 = q2 = q3 = 0.0; }

  /// Tell the number of characters required to print a quaternion, given that of a real number
  static inline size_t output_width (size_t const &real_width)
  {
    return 4*real_width + 13;
  }

  /// \brief Formatted output operator
  friend std::ostream & operator << (std::ostream &os, cvm::quaternion const &q);
  /// \brief Formatted input operator
  friend std::istream & operator >> (std::istream &is, cvm::quaternion &q);

  /// Access the quaternion as a 4-d array (return a reference)
  inline cvm::real & operator [] (int const &i) {
    switch (i) {
    case 0:
      return this->q0;
    case 1:
      return this->q1;
    case 2:
      return this->q2;
    case 3:
      return this->q3;
    default:
      cvm::fatal_error ("Error: incorrect quaternion component.\n");
      return q0;
    }
  }

  /// Access the quaternion as a 4-d array (return a value)
  inline cvm::real operator [] (int const &i) const {
    switch (i) {
    case 0:
      return this->q0;
    case 1:
      return this->q1;
    case 2:
      return this->q2;
    case 3:
      return this->q3;
    default:
      // deliberately return an invalid pointer
      cvm::fatal_error ("Error: trying to access an out-of-bounds "
                        "location in a quaternion.\n");
      return this->q0;
    }
  }

  /// Square norm of the quaternion
  inline cvm::real norm2() const
  {
    return q0*q0 + q1*q1 + q2*q2 + q3*q3;
  }

  /// Return the conjugate quaternion 
  inline cvm::quaternion conjugate() const
  {
    return cvm::quaternion (q0, -q1, -q2, -q3);
  }

  inline void operator *= (cvm::real const &a)
  {
    q0 *= a; q1 *= a; q2 *= a; q3 *= a;
  }

  inline void operator /= (cvm::real const &a)
  {
    q0 /= a; q1 /= a; q2 /= a; q3 /= a;
  }

  inline void set_positive()
  {
    if (q0 > 0.0) return;
    q0 = -q0;
    q1 = -q1;
    q2 = -q2;
    q3 = -q3;
  }

  inline void operator += (cvm::quaternion const &h)
  {
    q0+=h.q0; q1+=h.q1; q2+=h.q2; q3+=h.q3;
  }
  inline void operator -= (cvm::quaternion const &h)
  {
    q0-=h.q0; q1-=h.q1; q2-=h.q2; q3-=h.q3;
  }

  /// Promote a 3-vector to a quaternion
  static inline cvm::quaternion promote (cvm::rvector const &v)
  {
    return cvm::quaternion (0.0, v.x, v.y, v.z);
  }
  /// Return the vector component
  inline cvm::rvector get_rvector() const 
  {
    return cvm::rvector (q1, q2, q3);
  }


  friend inline cvm::quaternion operator + (cvm::quaternion const &h, cvm::quaternion const &q)
  {
    return cvm::quaternion (h.q0+q.q0, h.q1+q.q1, h.q2+q.q2, h.q3+q.q3);
  }

  friend inline cvm::quaternion operator - (cvm::quaternion const &h, cvm::quaternion const &q)
  {
    return cvm::quaternion (h.q0-q.q0, h.q1-q.q1, h.q2-q.q2, h.q3-q.q3);
  }

  /// \brief Provides the quaternion product.  \b NOTE: for inner
  /// product use: \code h.inner (q); \endcode
  friend inline cvm::quaternion operator * (cvm::quaternion const &h, cvm::quaternion const &q)
  {
    return cvm::quaternion (h.q0*q.q0 - h.q1*q.q1 - h.q2*q.q2 - h.q3*q.q3,
                            h.q0*q.q1 + h.q1*q.q0 + h.q2*q.q3 - h.q3*q.q2,
                            h.q0*q.q2 + h.q2*q.q0 + h.q3*q.q1 - h.q1*q.q3,
                            h.q0*q.q3 + h.q3*q.q0 + h.q1*q.q2 - h.q2*q.q1);
  }

  friend inline cvm::quaternion operator * (cvm::real const &c,
                                            cvm::quaternion const &q)
  {
    return cvm::quaternion (c*q.q0, c*q.q1, c*q.q2, c*q.q3);
  }
  friend inline cvm::quaternion operator * (cvm::quaternion const &q,
                                            cvm::real const &c)
  {
    return cvm::quaternion (q.q0*c, q.q1*c, q.q2*c, q.q3*c);
  }
  friend inline cvm::quaternion operator / (cvm::quaternion const &q,
                                            cvm::real const &c)
  {
    return cvm::quaternion (q.q0/c, q.q1/c, q.q2/c, q.q3/c);
  }


  /// \brief Rotate v through this quaternion (put it in the rotated
  /// reference frame)
  inline cvm::rvector rotate (cvm::rvector const &v) const
  {
    return ((*this) * promote (v) * ((*this).conjugate())).get_rvector();
  }

  /// \brief Rotate Q2 through this quaternion (put it in the rotated
  /// reference frame)
  inline cvm::quaternion rotate (cvm::quaternion const &Q2) const
  {
    //    return ((*this) * q * ((*this).conjugate()));
    cvm::rvector const vq_rot = this->rotate (Q2.get_rvector());
    return cvm::quaternion (Q2.q0, vq_rot.x, vq_rot.y, vq_rot.z);
  }

  /// \brief Return the square cosine of a quaternion rotation
  /// relative to another
  inline cvm::real cos2 (cvm::quaternion const &q) const
  {
    // get a vector orthogonal to both axes of rotation (*this and q)
    cvm::rvector const op = (cvm::rvector::outer (this->get_rvector(), q.get_rvector()));
    cvm::real const  opl2 = op.norm2();
    // rotate it with both quaternions and get the normalized inner product
    return ( (this->rotate (op)) * (q.rotate (op)) ) / (opl2);
  }


  /// \brief Square distance from q2 on the 4-dimensional sphere
  /// (square of the angle of the shortest geodesic)
  inline cvm::real dist2 (cvm::quaternion const &Q2) const
  {
    cvm::real const cos_theta = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;

    if (::pow (cos_theta, int (2)) > 1.0-1.0e-10) return 0.0;

    cvm::real const theta = ::acos (cos_theta);

    // get the minimum distance: x and -x are the same quaternion
    if (cos_theta > 0.0)
      return ::pow (theta, int (2));
    else
      return ::pow (M_PI-theta, int (2));
  }

  /// Gradient of the distance from q2
  inline cvm::quaternion dist2_grad (cvm::quaternion const &Q2) const
  {
    cvm::real const cos_theta = this->q0*Q2.q0 + this->q1*Q2.q1 + this->q2*Q2.q2 + this->q3*Q2.q3;
    double const theta = ::acos ( (cos_theta > 1.0) ? 1.0 : cos_theta );
    double const sin_theta = ::sin (theta);

    if (::pow (cos_theta, int (2)) > 1.0-1.0e-10) {
      // give the null element, avoid the fpe
      return cvm::quaternion();
    }

    cvm::quaternion const
      grad1 ((-1.0)*sin_theta*Q2.q0 + cos_theta*(this->q0-cos_theta*Q2.q0)/sin_theta,
             (-1.0)*sin_theta*Q2.q1 + cos_theta*(this->q1-cos_theta*Q2.q1)/sin_theta,
             (-1.0)*sin_theta*Q2.q2 + cos_theta*(this->q2-cos_theta*Q2.q2)/sin_theta,
             (-1.0)*sin_theta*Q2.q3 + cos_theta*(this->q3-cos_theta*Q2.q3)/sin_theta);

    if (cos_theta > 0.0) {
      return 2.0*theta*grad1;
    }
    else {
      return -2.0*(M_PI-theta)*grad1;
    }
  }

  /// \brief Choose the closest between Q2 and -Q2.  Not required in
  /// dist2() and dist2_grad()
  inline void match (cvm::quaternion &Q2) const
  {
    // same stuff as above

    cvm::real const cos_theta = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;

    cvm::real const theta = ::acos (cos_theta);

    if (theta > M_PI_2) Q2 *= -1.0;
  }

  /// \brief Inner product (as a 4-d vector) with Q2; requires match()
  /// if the largest overlap is looked for
  inline cvm::real inner (cvm::quaternion const &Q2) const
  {
    cvm::real const prod = this->q0*Q2.q0 + this->q1*Q2.q1 +
      this->q2*Q2.q2 + this->q3*Q2.q3;
    // assume the largest overlap between the two
    //    return (prod > -prod) ? prod : -prod;
    return prod;
  }


};


/// \brief Arbitrary size array (one dimensions) suitable for linear
/// algebra operations (i.e. for floating point numbers it can be used
/// with library functions)
template <class T, size_t const length> class colvarmodule::vector1d
{
protected:

  /// Underlying C-array
  T *array;

public:

  /// Length of the array
  inline size_t size()
  {
    return length;
  }
 
  /// Default constructor
  inline vector1d()
  {
    array = new T[length];
    reset();
  }

  /// Constructor from a 1-d C array
  inline vector1d (T const *v)
  {
    array = new T[length];
    for (size_t i = 0; i < length; i++) {
      array[i] = v[i];
    }
  }

  /// Copy constructor
  inline vector1d (vector1d const &v)
  {
    array = new T[length];
    for (size_t i = 0; i < length; i++) {
      array[i] = v.array[i];
    }
  }

  /// Destructor
  inline ~vector1d() {
    delete [] array;
  }

  /// Set all elements to zero
  inline void reset()
  {
    T null_value (0.0);
    for (size_t i = 0; i < length; i++) {
      array[i] = null_value;
    }
  }

  /// Return the 1-d C array
  inline T *c_array() { return array; }

  /// Return the 1-d C array
  inline operator T *() { return array; }

  /// Inner product
  inline friend T const operator * (vector1d<T, length> const &v1,
                                    vector1d<T, length> const &v2)
  {
    T prod (0.0);
    for (size_t i = 0; i < length; i++) {
      prod += v1.array[i] * v2.array[i];
    }
    return prod;
  }

  /// Formatted output 
  friend std::ostream & operator << (std::ostream &os,
                                     vector1d<T, length> const &v)
  {
    std::streamsize const w = os.width();
    std::streamsize const p = os.precision();

    os << "( ";
    for (size_t i = 0; i < length-1; i++) {
      os.width (w); os.precision (p);
      os << v.array[i] << " , ";
    }
    os.width (w); os.precision (p);
    os << v.array[length-1] << " )";
    return os;
  }

};



/// \brief Arbitrary size array (two dimensions) suitable for linear
/// algebra operations (i.e. for floating point numbers it can be used
/// with library functions)
template <class T,
          size_t const outer_length,
          size_t const inner_length> class colvarmodule::matrix2d
{
protected:

  /// Underlying C array
  T **array;

public:

  /// Allocation routine, used by all constructors
  inline void alloc() {
    array = new T * [outer_length];
    for (size_t i = 0; i < outer_length; i++) {
      array[i] = new T [inner_length];
    }
    reset();
  }

  /// Default constructor
  inline matrix2d()
  {
    this->alloc();
    reset();
  }
 
  /// Constructor from a 2-d C array
  inline matrix2d (T const **m)
  {
    this->alloc();
    for (size_t i = 0; i < outer_length; i++) {
      for (size_t j = 0; j < inner_length; j++)
        array[i][j] = m[i][j];
    }
  }

  /// Copy constructor
  inline matrix2d (matrix2d const &m)
  {
    this->alloc();
    for (size_t i = 0; i < outer_length; i++) {
      for (size_t j = 0; j < inner_length; j++)
        array[i][j] = m.array[i][j];
    }
  }

  /// Set all elements to zero
  inline void reset (T const &null_value = T())
  {
    for (size_t i = 0; i < outer_length; i++) {
      for (size_t j = 0; j < inner_length; j++) {
        array[i][j] = null_value;
      }
    }
  }

  /// Destructor
  inline ~matrix2d() {
    for (size_t i = 0; i < outer_length; i++) {
      delete [] array[i];
    }
    delete [] array;
  }

  /// Return the 2-d C array
  inline T **c_array() { return array; }

  /// Return the 2-d C array
  inline operator T **() { return array; }

};


/// \brief 2-dimensional array of real numbers with three components
/// along each dimension (works with colvarmodule::rvector)
class colvarmodule::rmatrix
  : public colvarmodule::matrix2d<colvarmodule::real, 3, 3> {
private:

public:

  /// Return the xx element
  inline cvm::real &xx() { return array[0][0]; }
  /// Return the xy element
  inline cvm::real &xy() { return array[0][1]; }
  /// Return the xz element
  inline cvm::real &xz() { return array[0][2]; }
  /// Return the yx element
  inline cvm::real &yx() { return array[1][0]; }
  /// Return the yy element
  inline cvm::real &yy() { return array[1][1]; }
  /// Return the yz element
  inline cvm::real &yz() { return array[1][2]; }
  /// Return the zx element
  inline cvm::real &zx() { return array[2][0]; }
  /// Return the zy element
  inline cvm::real &zy() { return array[2][1]; }
  /// Return the zz element
  inline cvm::real &zz() { return array[2][2]; }

  /// Constructor from a 2-d C array
  inline rmatrix (cvm::real const **m) 
    : cvm::matrix2d<cvm::real, 3, 3> (m) 
  {}

  /// Default constructor
  inline rmatrix() 
    : cvm::matrix2d<cvm::real, 3, 3>()
  {}

  /// Constructor component by component 
  inline rmatrix (cvm::real const &xxi,
                  cvm::real const &xyi,
                  cvm::real const &xzi,
                  cvm::real const &yxi,
                  cvm::real const &yyi,
                  cvm::real const &yzi,
                  cvm::real const &zxi,
                  cvm::real const &zyi,
                  cvm::real const &zzi) 
    : cvm::matrix2d<cvm::real, 3, 3>()
  {
    this->xx() = xxi; this->xy() = xyi; this->xz() = xzi;
    this->yx() = yxi; this->yy() = yyi; this->yz() = yzi;
    this->zx() = zxi; this->zy() = xyi; this->zz() = zzi;
  }

  /// Destructor
  inline ~rmatrix()
  {}    

  /// Return the determinant
  inline cvm::real determinant() {
    return
      xx() * (yy()*zz() - zy()*yz())
      - yx() * (xy()*zz() - zy()*xz())
      + zx() * (xy()*yz() - yy()*xz());
  }

  friend cvm::rvector operator * (cvm::rmatrix &m, cvm::rvector &r);

  inline cvm::rmatrix transpose() {
    return cvm::rmatrix (this->xx(), this->yx(), this->zx(),
                         this->xy(), this->yy(), this->zy(),
                         this->xz(), this->yz(), this->zz());
  }

  //   inline cvm::rmatrix const operator * (cvm::rmatrix const &m1, cvm::rmatrix const &m2) {
  //     return cvm::rmatrix (m1.xx()*m2.xx() + m1.xy()*m2.yx() + m1.xz()*m2.yz(),
  //                     m1.xx()*m2.xy() + m1.xy()*m2.yy() + m1.xz()*m2.zy(),
  //                     m1.xx()*m2.xz() + m1.xy()*m2.yz() + m1.xz()*m2.zz(),
  //                     m1.yx()*m2.xx() + m1.yy()*m2.yx() + m1.yz()*m2.yz(),
  //                     m1.yx()*m2.xy() + m1.yy()*m2.yy() + m1.yz()*m2.yy(),
  //                     m1.yx()*m2.xz() + m1.yy()*m2.yz() + m1.yz()*m2.yz(),
  //                     m1.zx()*m2.xx() + m1.zy()*m2.yx() + m1.zz()*m2.yz(),
  //                     m1.zx()*m2.xy() + m1.zy()*m2.yy() + m1.zz()*m2.yy(),
  //                     m1.zx()*m2.xz() + m1.zy()*m2.yz() + m1.zz()*m2.yz());
  //   }

};

                    
inline cvm::rvector operator * (cvm::rmatrix &m,
                                cvm::rvector &r)
{
  return cvm::rvector (m.xx()*r.x + m.xy()*r.y + m.xz()*r.z,
                       m.yx()*r.x + m.yy()*r.y + m.yz()*r.z,
                       m.zx()*r.x + m.zy()*r.y + m.zz()*r.z);
}


/// Numerical recipes diagonalization
void nr_jacobi (cvm::real **a, int n, cvm::real d[], cvm::real **v, int *nrot);

/// Eigenvector sort
void nr_eigsrt (cvm::real d[], cvm::real **v, int n);

/// Transpose the matrix
void nr_transpose (cvm::real **v, int n);



/// \brief A rotation between two sets of coordinates (for the moment
/// a wrapper for colvarmodule::quaternion)
class colvarmodule::rotation
{
public:

  /// \brief Positions to superimpose: the calculated rotation brings
  /// pos1 to superimpose pos2
  std::vector< cvm::atom_pos > pos1, pos2;

  /// Derivatives of S
  std::vector< cvm::matrix2d<cvm::rvector, 4, 4> > dS_1,  dS_2;
  /// Derivatives of leading eigenvalue
  std::vector< cvm::rvector >                      dL0_1, dL0_2;
  /// Derivatives of leading eigenvector
  std::vector< cvm::vector1d<cvm::rvector, 4> >    dQ0_1, dQ0_2;

  /// Default constructor
  inline rotation()
    : degeneracy (false)
  {}

  /// Constructor after a quaternion
  inline rotation (cvm::quaternion const &qi)
    : degeneracy (false)
  {
    q = qi;
  }

  /// Constructor after an axis of rotation and an angle
  inline rotation (cvm::real const &angle, cvm::rvector const &axis)
    : degeneracy (false)
  {
    cvm::rvector const an = axis.unit();
    cvm::real const sina = ::sin (angle/2.0);
    q = cvm::quaternion (::cos (angle/2.0), sina * an.x,
                         sina * an.y, sina * an.z);
  }

  /// Destructor
  inline ~rotation()
  {}

  /// Return the rotated vector
  inline cvm::rvector rotate (cvm::rvector const &v) const
  {
    return q.rotate (v);
  }

  /// Return the inverse of this rotation
  inline cvm::rotation inverse()
  {
    return cvm::rotation (this->q.conjugate());
  }

  /// \brief Calculate the optimal rotation
  ///
  /// The method used is defined in:
  /// Coutsias EA, Seok C, Dill KA.
  /// Using quaternions to calculate RMSD.
  /// J Comput Chem. 25(15):1849-57 (2004)
  /// DOI: 10.1002/jcc.20110  PubMed: 15376254
  void calc_optimal_rotation (std::vector<atom_pos> const &pos1,
                              std::vector<atom_pos> const &pos2,
                              cvm::real                   &l0,
                              cvm::quaternion             &q0);
  
  /// \brief Calculate the optimal rotation and its derivatives
  ///
  /// The method used is defined in:
  /// Coutsias EA, Seok C, Dill KA.
  /// Using quaternions to calculate RMSD.
  /// J Comput Chem. 25(15):1849-57 (2004)
  /// DOI: 10.1002/jcc.20110  PubMed: 15376254
  void calc_optimal_rotation (std::vector<atom_pos> const &pos1,
                              std::vector<atom_pos> const &pos2,
                              cvm::real                   &l0,
                              cvm::quaternion             &q0,
                              std::vector< cvm::matrix2d<cvm::rvector, 4, 4> > &dS_1,
                              std::vector< cvm::matrix2d<cvm::rvector, 4, 4> > &dS_2,
                              std::vector< cvm::rvector >                      &dL0_1,
                              std::vector< cvm::rvector >                      &dL0_2,
                              std::vector< cvm::vector1d<cvm::rvector, 4> >    &dQ0_1,
                              std::vector< cvm::vector1d<cvm::rvector, 4> >    &dQ0_2);

  /// Build the overlap matrix S
  void build_matrix (std::vector<cvm::atom_pos> const &pos1,
                     std::vector<cvm::atom_pos> const &pos2,
                     cvm::matrix2d<real, 4, 4>        &S);

  /// Diagonalize the overlap matrix S
  void diagonalize_matrix (cvm::matrix2d<cvm::real, 4, 4> &S,
                           cvm::real                       S_eigval[4],
                           cvm::matrix2d<cvm::real, 4, 4> &S_eigvec);

  /// \brief The calculated rotation
  cvm::quaternion q;

protected:

  /// \brief If a structure becomes too distorted, the eigenvalue may
  /// be degenerate, in this case warn the user
  bool degeneracy;
};


#endif

// Emacs
// Local Variables:
// mode: C++
// End:
