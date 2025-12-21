/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Enrico Skoruppa (Physics of Life, TU Dresden, Dresden)
------------------------------------------------------------------------- */

#ifndef LMP_LAMATH_H
#define LMP_LAMATH_H

#include <cmath>
#include <algorithm>  // for std::max
#include "lammps.h"

namespace lamath {

inline double norm(const double x[3]) noexcept {
  return std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

// -----------------------------------------------------------------
// Matrix-vector product  y = A * x
// A  : double[3][3]  (row-major, as in LAMMPS)
// x,y: double[3]
// -----------------------------------------------------------------
inline void mul(const double A[3][3],
                const double x[3],
                double       y[3]) noexcept
{
  y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
  y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
  y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
}

// -----------------------------------------------------------------
// Matrix-matrix product  C = A * B
// All arguments: double[3][3]  (row-major)
// Safe when C aliases neither A nor B.
// -----------------------------------------------------------------
inline void mul(const double A[3][3],
                const double B[3][3],
                double       C[3][3]) noexcept
{
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
  C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];

  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
  C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];

  C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
  C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
  C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
}

inline void add(const double x[3],
                const double y[3],
                double       z[3]) noexcept
{
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
}

inline void signflip( const double x[3],
                      double y[3]) noexcept
{
  y[0] = -x[0];
  y[1] = -x[1];
  y[2] = -x[2];
}


// -----------------------------------------------------------------
// Transpose
// -----------------------------------------------------------------
inline void transpose(  const double T[3][3],
                        double Ttransp[3][3]) noexcept
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Ttransp[i][j] = T[j][i];  
        }
    }
}


// -----------------------------------------------------------------
// Compute determinant of a 3×3 matrix
// mat: double[3][3] in row-major order
// -----------------------------------------------------------------
inline double det(const double mat[3][3]) noexcept {
    return 
        mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) 
      - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) 
      + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

inline bool is_positive_definite(double M[6][6]) noexcept {
  constexpr int N = 6;
  double L[N][N] = {0};

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = 0.0;
      for (int k = 0; k < j; ++k)
        sum += L[i][k] * L[j][k];

      if (i == j) {
        double val = M[i][i] - sum;
        if (val <= 0.0) return false;
        L[i][j] = std::sqrt(val);
      } else {
        L[i][j] = (1.0 / L[j][j]) * (M[i][j] - sum);
      }
    }
  }
  return true;
}

inline bool is_positive_definite(double M[3][3]) noexcept {
  constexpr int N = 3;
  double L[N][N] = {0};

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = 0.0;
      for (int k = 0; k < j; ++k)
        sum += L[i][k] * L[j][k];

      if (i == j) {
        double val = M[i][i] - sum;
        if (val <= 0.0) return false;
        L[i][j] = std::sqrt(val);
      } else {
        L[i][j] = (1.0 / L[j][j]) * (M[i][j] - sum);
      }
    }
  }
  return true;
}


} // namespace lam

#endif
