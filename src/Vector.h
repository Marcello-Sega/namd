/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1995 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*								   	   */
/***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Vector.h,v $
 *	$Author: nealk $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1996/11/04 19:44:42 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Implement a vector class.  This means we can do vector manipulation
 * with something like v3=v1+v2.  I don't know how fast it is compared
 * to the old fashioned way.  We may have to change things once it gets
 * working.
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Vector.h,v $
 * Revision 1.2  1996/11/04 19:44:42  nealk
 * Replaced old vector.h with new one.
 *
 * Revision 1.15  1996/10/17 15:17:11  nealk
 * Removed a negation from all cross pruoduct functions... :-)
 *
 * Revision 1.14  1996/10/14 14:02:01  nealk
 * Added negation "-" to vector class.
 *
 * Revision 1.13  1996/10/12 00:11:00  brunner
 * Made several functions const, since they don't modify the object
 *
 * Revision 1.12  1996/09/25 19:20:51  nealk
 * Added "inline" and it seems to make a small difference.
 *
 * Revision 1.11  1996/09/24 14:54:33  nealk
 * Added cross with scalar since this is a very common use for cross products.
 *
 * Revision 1.10  1996/09/17 16:40:19  nealk
 * Modified +=, *=, and -= to not return a value (void).  This makes them
 * as efficient as explicit declarations.  Note: they cannot be used to
 * return values (but the return isn't used in any of the code anyway).
 * Also added /=.
 *
 * Revision 1.9  1996/04/27 20:42:24  billh
 * Added insertion  operator for Vector, and changed how Vector.h included.
 *
 * Revision 1.8  1996/04/17 16:42:07  billh
 * Added v1 = v2, v1 = const, v1 *= const operators
 *
 * Revision 1.7  1995/04/04  18:40:48  dalke
 * added a string "set"
 *
 * Revision 1.6  1995/03/08  14:29:57  nelson
 * Added copyright
 *
 * Revision 1.5  94/09/04  20:18:50  20:18:50  nelson (Mark T. Nelson)
 * add function div and fixed length function
 * 
 * Revision 1.4  94/08/09  13:20:02  13:20:02  nelson (Mark T. Nelson)
 * Added add_const routine
 * 
 * Revision 1.3  94/07/07  16:47:45  16:47:45  dalke (Andrew Dalke)
 * changed data type from Real to BigReal
 * 
 * Revision 1.2  94/07/06  11:55:57  11:55:57  dalke (Andrew Dalke)
 * Added unit and length function, also added two-element operations
 * to minimize copies.
 * 
 * Revision 1.1  94/07/06  11:13:22  11:13:22  dalke (Andrew Dalke)
 * Initial revision
 * 
 ***************************************************************************/

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream.h>
#include <math.h>
#include <stdio.h>
#include "common.h"
#include "Inform.h"

class Vector {
   public:
     BigReal x,y,z;
     
     inline Vector( void) {         // default is to create a 0 vector
       x = y = z = 0.0;
     }
     inline Vector( const Vector &v2) { // Vector x = another_vector
       x = v2.x;
       y = v2.y;
       z = v2.z;
     }
     inline Vector( BigReal newx, BigReal newy, BigReal newz) {
       x = newx;
       y = newy;
       z = newz;
     }
     ~Vector( void) {
     }

     inline BigReal &operator[](int i) {
       return i==0 ? x
             :i==1 ? y
             :i==2 ? z
             :(NAMD_die("Vector reference out of bounds."), x);

     }

     //  v1 = v2;
     inline Vector& operator=(const Vector &v2) {
       x = v2.x;
       y = v2.y;
       z = v2.z;
       return *this;
     }

     //  v1 = const;
     inline Vector& operator=(const BigReal &v2) {
       x = v2;
       y = v2;
       z = v2;
       return *this;
     }

     //  v1 += v2;
     inline void operator+=(const Vector &v2) {
       x += v2.x;
       y += v2.y;
       z += v2.z;
     }

     // v1 -= v2;
     inline void operator-=(const Vector &v2) {
       x -= v2.x;
       y -= v2.y;
       z -= v2.z;
     }

     // v1 *= const
     inline void operator*=(const BigReal &v2) {
       x *= v2;
       y *= v2;
       z *= v2;
     }

     // v1 /= const
     inline void operator/=(const BigReal& v2) {
       x /= v2;
       y /= v2;
       z /= v2;
     }

     inline friend int operator == (const Vector& v1, const Vector& v2) {
       return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
     }
     inline friend int operator != (const Vector& v1, const Vector& v2) {
       // return !(v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
       return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
     }

     // addition of two vectors
     inline friend Vector operator+(const Vector& v1, const Vector& v2) {
       return Vector( v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
     }

     // negation
     inline friend Vector operator-(const Vector &v1) {
       return Vector( -v1.x, -v1.y, -v1.z);
     }

     // subtraction
     inline friend Vector operator-(const Vector &v1, const Vector &v2) {
       return Vector( v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
     }
     // inner ("dot") product
     inline friend BigReal operator*(const Vector &v1, const Vector &v2) {
       return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
     }
     // scalar product
     inline friend Vector operator*(const BigReal &f, const Vector &v1) {
       return Vector(f*v1.x, f*v1.y, f*v1.z);
     }
     // scalar product
     inline friend Vector operator*(const Vector &v1, const BigReal &f) {
       return Vector(f*v1.x, f*v1.y, f*v1.z);
     }
     // division by a scalar
     inline friend Vector operator/(const Vector &v1, const BigReal &f) {
//       if (!f)
//         NAMD_die("Division by 0 on a vector operation.");
       return Vector(v1.x/f, v1.y/f, v1.z/f);
     }
     
     // return the norm
     inline BigReal length(void) const {
       return sqrt(x*x+y*y+z*z);
     }
     
     inline BigReal length2(void) const {
       return (x*x + y*y + z*z);
     }

     // return the unit vector in the same direction
     inline Vector unit(void) const {
       return Vector(x, y, z)/length();
     }
     
     
     // one cross product  v3 = cross(v1, v2)
     inline friend Vector cross(const Vector &v1, const Vector &v2) {
       return Vector( v1.y*v2.z-v2.y*v1.z,
                     // -v1.x*v2.z+v2.x*v1.z,
                      v2.x*v1.z-v1.x*v2.z,
                      v1.x*v2.y-v2.x*v1.y  );
     }

     // multiplying a cross product by a scalar is very common
     // one cross product  v3 = k*cross(v1, v2)
     inline friend Vector cross(const Real &k, const Vector &v1, const Vector &v2) {
       return Vector( k*(v1.y*v2.z-v2.y*v1.z),
                      // k*(-v1.x*v2.z+v2.x*v1.z),
                      k*(v2.x*v1.z-v1.x*v2.z),
                      k*(v1.x*v2.y-v2.x*v1.y) );
     }

     inline friend Vector cross(const BigReal &k, const Vector &v1, const Vector &v2) {
       return Vector( k*(v1.y*v2.z-v2.y*v1.z),
                      // k*(-v1.x*v2.z+v2.x*v1.z),
                      k*(v2.x*v1.z-v1.x*v2.z),
                      k*(v1.x*v2.y-v2.x*v1.y) );
     }

     // print out
     friend ostream& operator<<(ostream& strm, const Vector &v1) {
       strm << "( "<< v1.x << ", " << v1.y << ", " << v1.z << ')';
       return strm;
     }

     // print out to Inform object
     friend Inform& operator<<(Inform& strm, const Vector &v1) {
       strm << "( "<< v1.x << ", " << v1.y << ", " << v1.z << ")";
       return strm;
     }

     // add a vector to this vector
     void add(const Vector &v1) {
       x+=v1.x; y+=v1.y; z+=v1.z;
     }

     // subtract the vector from this one
     void sub(const Vector &v1) {
       x-=v1.x; y-=v1.y; z-=v1.z;
     }

     // add a constant factor to each element of a vector
     void add_const(BigReal c)
     {
	x+=c;
	y+=c;
	z+=c;
     }


     // A = A x B  -- why do you want this function, anyway?
     void cross(const Vector &v2) {
       BigReal xx =  y*v2.z-v2.y*z;
       // BigReal yy = -x*v2.z+v2.x*z;
       BigReal yy = v2.x*z-x*v2.z;
       z =  x*v2.y-v2.x*y;
       y=yy;
       x=xx;
     }
     // rescale everything by a scalar -- V = a*V
     void mult(BigReal f) {
       x*=f; y*=f; z*=f;
     }

     // divide each element by a scalar
     void div(BigReal f) {
       x/=f; y/=f; z/=f;
     }

     // returns (*this) * V2
     BigReal dot(const Vector &v2) {
       return x*v2.x + y*v2.y + z*v2.z;
     }

     // set the vector based on a string.  If bad, return FALSE
     // the string can be in the form "x y z" or "x, y, z"
     Bool set(const char *s) {
	double a[3];    // read into doubles, since I don't know what
	char tmp[100];  // a "BigReal" is in real life
	// cheap way to get commas, etc.  a poor regex
       int i=sscanf(s, "%lf%99[ \t,]%lf%99[ \t,]%lf%99s",
                    a, tmp, a+1, tmp, a+2, tmp);
       if (i != 5) return FALSE;
       const char *t = s;       // now count commas (for "1,,,,2,  , 3")
       int flg = 0;                 // and check for "1 2,,3"
       i = 0;
       for (;*t;t++) {
          if (*t == ',') { 
             if (flg == 0) {   // expecting non-whitespace
                return FALSE;  //    so error
             }
             flg = 0;          // now expect non-whitespace
             i++;              // and increment comma counter
          }
          else if (*t != ' ' && *t != '\t') {  // got non-whitespace
             flg = 1;          // so the next can be whitespace or commas
          }
       }
       if (i == 0 || i == 2) {  // allow "1 2 3" or "1, 2,3" forms
          x = a[0]; y = a[1]; z = a[2];
          return TRUE;
       }
       return FALSE;
     }
};
//#define TEST_VECTOR_CLASS
#ifdef TEST_VECTOR_CLASS
main()
{
  Vector v1(1.1,2.2, 3.3);
  Vector v2(-1, 55, 32.1);
  Vector v3(v1+2*v2);
  Vector v4;
  cout << v1 << "  " << v2 << "  " << v3 << "  " << v4 << '\n';
  cout << v1*v2 << "  "  << v3-v1-2*v2 <<"  "<< v2 * v3 <<"  "<< v3*v2 <<'\n';
  v4 = v3*5 - v2/4;
  cout << v4 << "  " << v3*5.0 - v2/4.0 << '\n';
  cout << v4[0] << "  "  << v4[1] << "  " << v4[2] << '\n';
//  cout.flush();
//  cout << v4[3];
  cout << cross(v1, v2) << '\n';
  cout << v1 << '\n';  
  v1 += v2;
  cout << v1 << '\n';
  v1 -= v2;
  cout << v1 << '\n';
  {
    Vector v1(1.0, 2.0, 3.0);  // some more examples, but I was too lazy to
    Vector v2 = v1.unit();     // fix the names
    cout << v2 << '\n';
    cout << v2.dot(v1) << '\n';
    cout << v1.length() << '\n';
    v1.mult(-1);
    v1.add(v2*14);
    cout << v1 << '\n';
  }
}
#endif

#endif
