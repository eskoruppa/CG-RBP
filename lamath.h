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
   Contributing author: Enrico Skoruppa (School of Physics and Astronomy, University of Edinburgh, Edinburgh)
------------------------------------------------------------------------- */

#ifndef LMP_LAMATH_H
#define LMP_LAMATH_H

#include <cmath>
#include <algorithm>  // for std::max
#include "lammps.h"

namespace lamath {

static inline void print_mat3_lammps(FILE *fp, const char *name, double M[3][3])
{
  fprintf(fp, "%s =\n", name);
  for (int i = 0; i < 3; i++)
    fprintf(fp, "  % .4e  % .4e  % .4e\n",
            M[i][0], M[i][1], M[i][2]);
  fprintf(fp, "\n");
}

inline double norm(const double x[3]) noexcept {
  return std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

inline void mul(const double A[3][3],
                const double x[3],
                double       y[3]) noexcept
{
  y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
  y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
  y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
}

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

inline void mul_Atx(const double A[3][3], const double x[3], double y[3]) noexcept {
    y[0] = A[0][0]*x[0] + A[1][0]*x[1] + A[2][0]*x[2];
    y[1] = A[0][1]*x[0] + A[1][1]*x[1] + A[2][1]*x[2];
    y[2] = A[0][2]*x[0] + A[1][2]*x[1] + A[2][2]*x[2];
}

inline void mul_AtB(const double A[3][3],
                const double B[3][3],
                double       C[3][3]) noexcept
{
  C[0][0] = A[0][0]*B[0][0] + A[1][0]*B[1][0] + A[2][0]*B[2][0];
  C[0][1] = A[0][0]*B[0][1] + A[1][0]*B[1][1] + A[2][0]*B[2][1];
  C[0][2] = A[0][0]*B[0][2] + A[1][0]*B[1][2] + A[2][0]*B[2][2];

  C[1][0] = A[0][1]*B[0][0] + A[1][1]*B[1][0] + A[2][1]*B[2][0];
  C[1][1] = A[0][1]*B[0][1] + A[1][1]*B[1][1] + A[2][1]*B[2][1];
  C[1][2] = A[0][1]*B[0][2] + A[1][1]*B[1][2] + A[2][1]*B[2][2];

  C[2][0] = A[0][2]*B[0][0] + A[1][2]*B[1][0] + A[2][2]*B[2][0];
  C[2][1] = A[0][2]*B[0][1] + A[1][2]*B[1][1] + A[2][2]*B[2][1];
  C[2][2] = A[0][2]*B[0][2] + A[1][2]*B[1][2] + A[2][2]*B[2][2];
}


inline void add(const double x[3],
                const double y[3],
                double       z[3]) noexcept
{
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  z[2] = x[2] + y[2];
}

inline void add_to( double       x[3],
                    const double y[3]) noexcept
{
  x[0] += y[0];
  x[1] += y[1];
  x[2] += y[2];
}

inline void subtract( const double x[3],
                      const double y[3],
                      double       z[3]) noexcept
{
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  z[2] = x[2] - y[2];
}

inline void subtract_from( double       x[3],
                           const double y[3]) noexcept
{
  x[0] -= y[0];
  x[1] -= y[1];
  x[2] -= y[2];
}

inline void signflip( const double x[3],
                      double y[3]) noexcept
{
  y[0] = -x[0];
  y[1] = -x[1];
  y[2] = -x[2];
}

inline void transpose(  const double T[3][3],
                        double Ttransp[3][3]) noexcept
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Ttransp[i][j] = T[j][i];  
        }
    }
}

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

// inline bool is_positive_definite(double M[3][3]) noexcept {
//   constexpr int N = 3;
//   double L[N][N] = {0};

//   for (int i = 0; i < N; ++i) {
//     for (int j = 0; j <= i; ++j) {
//       double sum = 0.0;
//       for (int k = 0; k < j; ++k)
//         sum += L[i][k] * L[j][k];

//       if (i == j) {
//         double val = M[i][i] - sum;
//         if (val <= 0.0) return false;
//         L[i][j] = std::sqrt(val);
//       } else {
//         L[i][j] = (1.0 / L[j][j]) * (M[i][j] - sum);
//       }
//     }
//   }
//   return true;
// }

inline bool is_positive_definite( const double M[3][3],
                                  double sym_tol = 1e-12,
                                  double pd_tol  = 1e-15,
                                  bool symmetrize = false) noexcept
{
  // Finite check
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (!std::isfinite(M[i][j])) return false;
    }
  }

  // Scale-aware symmetry tolerance:
  // compare |M_ij - M_ji| against sym_tol * max(1, max(|M_ij|,|M_ji|))
  double A[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      A[i][j] = M[i][j];
    }
  }

  // Symmetry check (and optional symmetrization)
  for (int i = 0; i < 3; ++i) {
    for (int j = i + 1; j < 3; ++j) {
      const double a = A[i][j];
      const double b = A[j][i];
      const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
      const double diff  = std::abs(a - b);

      if (diff > sym_tol * scale) return false;

      if (symmetrize) {
        const double s = 0.5 * (a + b);
        A[i][j] = s;
        A[j][i] = s;
      }
    }
  }

  // Cholesky with a defensible diagonal threshold.
  // Make pd_tol scale-aware too: pd_tol * max(1, max diagonal magnitude)
  double max_diag = 0.0;
  for (int i = 0; i < 3; ++i) max_diag = std::max(max_diag, std::abs(A[i][i]));
  const double diag_eps = pd_tol * std::max(1.0, max_diag);

  double L[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = 0.0;
      for (int k = 0; k < j; ++k)
        sum += L[i][k] * L[j][k];

      if (i == j) {
        const double val = A[i][i] - sum;
        // val must be strictly positive for SPD. Allow a tiny epsilon for roundoff.
        if (!(val > diag_eps)) return false;
        L[i][j] = std::sqrt(val);
      } else {
        // L[j][j] should be > 0 if we passed the diagonal checks, but guard anyway
        if (!(L[j][j] > 0.0)) return false;
        L[i][j] = (A[i][j] - sum) / L[j][j];
      }
    }
  }

  return true;
}

} // namespace lamath
#endif
