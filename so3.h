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

#ifndef LMP_SO3_H
#define LMP_SO3_H

#include <cmath>
#include <algorithm>  // for std::max
#include "lamath.h"

namespace so3 {

inline void hat(const double v[3], double M[3][3]) noexcept {
    M[0][0] =   0.0;   M[0][1] = -v[2];  M[0][2] =  v[1];
    M[1][0] =  v[2];   M[1][1] =   0.0;  M[1][2] = -v[0];
    M[2][0] = -v[1];   M[2][1] =  v[0];  M[2][2] =   0.0;
}

inline void triads2rotmat(  const double Ti[3][3],
                            const double Tj[3][3],
                            double       R[3][3]) noexcept
{
    // R[a][b] = sum_k Ti[k][a] * Tj[k][b]
    R[0][0] = Ti[0][0]*Tj[0][0] + Ti[1][0]*Tj[1][0] + Ti[2][0]*Tj[2][0];
    R[0][1] = Ti[0][0]*Tj[0][1] + Ti[1][0]*Tj[1][1] + Ti[2][0]*Tj[2][1];
    R[0][2] = Ti[0][0]*Tj[0][2] + Ti[1][0]*Tj[1][2] + Ti[2][0]*Tj[2][2];

    R[1][0] = Ti[0][1]*Tj[0][0] + Ti[1][1]*Tj[1][0] + Ti[2][1]*Tj[2][0];
    R[1][1] = Ti[0][1]*Tj[0][1] + Ti[1][1]*Tj[1][1] + Ti[2][1]*Tj[2][1];
    R[1][2] = Ti[0][1]*Tj[0][2] + Ti[1][1]*Tj[1][2] + Ti[2][1]*Tj[2][2];

    R[2][0] = Ti[0][2]*Tj[0][0] + Ti[1][2]*Tj[1][0] + Ti[2][2]*Tj[2][0];
    R[2][1] = Ti[0][2]*Tj[0][1] + Ti[1][2]*Tj[1][1] + Ti[2][2]*Tj[2][1];
    R[2][2] = Ti[0][2]*Tj[0][2] + Ti[1][2]*Tj[1][2] + Ti[2][2]*Tj[2][2];
}

inline void euler2rotmat(const double omega[3], double R[3][3]) noexcept
{
    // Compute θ = ||omega||
    double theta2 = omega[0]*omega[0] + omega[1]*omega[1] + omega[2]*omega[2];
    if (theta2 < 1e-12) {
        // Small angle: use series R ≈ I + [ω]_× + ½ [ω]_×²
        double hatOmega[3][3], hatOmega2[3][3];
        hat(omega, hatOmega);
        lamath::mul(hatOmega, hatOmega, hatOmega2);

        // R = I + hatOmega + 0.5*hatOmega2
        R[0][0] = 1.0 + hatOmega[0][0] + 0.5*hatOmega2[0][0];
        R[0][1] =       hatOmega[0][1] + 0.5*hatOmega2[0][1];
        R[0][2] =       hatOmega[0][2] + 0.5*hatOmega2[0][2];

        R[1][0] =       hatOmega[1][0] + 0.5*hatOmega2[1][0];
        R[1][1] = 1.0 + hatOmega[1][1] + 0.5*hatOmega2[1][1];
        R[1][2] =       hatOmega[1][2] + 0.5*hatOmega2[1][2];

        R[2][0] =       hatOmega[2][0] + 0.5*hatOmega2[2][0];
        R[2][1] =       hatOmega[2][1] + 0.5*hatOmega2[2][1];
        R[2][2] = 1.0 + hatOmega[2][2] + 0.5*hatOmega2[2][2];
        return;
    }

    double theta = std::sqrt(theta2);
    double ux = omega[0]/theta;
    double uy = omega[1]/theta;
    double uz = omega[2]/theta;

    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);
    double one_minus = 1.0 - cos_t;

    // Fill R = I*cosθ + (1-cosθ) u u^T + [u]_× sinθ
    R[0][0] = cos_t + ux*ux * one_minus;
    R[0][1] = ux*uy * one_minus - uz*sin_t;
    R[0][2] = ux*uz * one_minus + uy*sin_t;

    R[1][0] = uy*ux * one_minus + uz*sin_t;
    R[1][1] = cos_t + uy*uy * one_minus;
    R[1][2] = uy*uz * one_minus - ux*sin_t;

    R[2][0] = uz*ux * one_minus - uy*sin_t;
    R[2][1] = uz*uy * one_minus + ux*sin_t;
    R[2][2] = cos_t + uz*uz * one_minus;
}

inline void rotmat2euler(const double R[3][3], double omega[3]) noexcept
{
    // Compute cosθ = (trace(R)-1)/2
    double tr = R[0][0] + R[1][1] + R[2][2];
    double cos_t = (tr - 1.0) * 0.5;
    // Clamp to [-1,1]
    if      (cos_t >  1.0) cos_t =  1.0;
    else if (cos_t < -1.0) cos_t = -1.0;
    double theta = std::acos(cos_t);

    if (theta < 1e-8) {
        // θ ≈ 0 → omega = ½ vee(R - R^T)
        omega[0] = 0.5 * (R[2][1] - R[1][2]);
        omega[1] = 0.5 * (R[0][2] - R[2][0]);
        omega[2] = 0.5 * (R[1][0] - R[0][1]);
        return;
    }

    double sin_t = std::sin(theta);
    if (std::fabs(sin_t) > 1e-6) {
        // General case:  omega = (θ / (2 sinθ)) * vee(R - R^T)
        double coeff = 0.5 * theta / sin_t;
        omega[0] = coeff * (R[2][1] - R[1][2]);
        omega[1] = coeff * (R[0][2] - R[2][0]);
        omega[2] = coeff * (R[1][0] - R[0][1]);
        return;
    }

    // θ ≈ π:  need stable extraction of axis from diagonal
    // Find max diagonal element among R[0][0],R[1][1],R[2][2]
    double diag0 = (R[0][0] + 1.0) * 0.5;  // = ux^2
    double diag1 = (R[1][1] + 1.0) * 0.5;  // = uy^2
    double diag2 = (R[2][2] + 1.0) * 0.5;  // = uz^2
    if (diag0 >= diag1 && diag0 >= diag2) {
        double ux = std::sqrt(diag0);
        double uy = (R[0][1] + R[1][0]) * (0.5/ux);
        double uz = (R[0][2] + R[2][0]) * (0.5/ux);
        omega[0] = ux * theta;
        omega[1] = uy * theta;
        omega[2] = uz * theta;
    }
    else if (diag1 >= diag0 && diag1 >= diag2) {
        double uy = std::sqrt(diag1);
        double ux = (R[0][1] + R[1][0]) * (0.5/uy);
        double uz = (R[1][2] + R[2][1]) * (0.5/uy);
        omega[0] = ux * theta;
        omega[1] = uy * theta;
        omega[2] = uz * theta;
    }
    else {
        double uz = std::sqrt(diag2);
        double ux = (R[0][2] + R[2][0]) * (0.5/uz);
        double uy = (R[1][2] + R[2][1]) * (0.5/uz);
        omega[0] = ux * theta;
        omega[1] = uy * theta;
        omega[2] = uz * theta;
    }
}

inline void triads2euler(const double Ti[3][3],
                         const double Tj[3][3],
                         const double S[3][3],
                         double       Omega[3]) noexcept
{
    double R[3][3], D[3][3];
    triads2rotmat(Ti, Tj, R);
    lamath::mul_AtB(S, R, D);
    rotmat2euler(D, Omega);
}


inline void rotmat2quat(const double R[3][3], double q[4]) noexcept {

  double trace = R[0][0] + R[1][1] + R[2][2];

  if (trace > 0.0) {
    double S = std::sqrt(trace + 1.0) * 2.0;
    q[0] = 0.25 * S;
    q[1] = (R[1][2] - R[2][1]) / S;
    q[2] = (R[2][0] - R[0][2]) / S;
    q[3] = (R[0][1] - R[1][0]) / S;
  } else {
    if ((R[0][0] > R[1][1]) && (R[0][0] > R[2][2])) {
      double S = std::sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2.0;
      q[0] = (R[1][2] - R[2][1]) / S;
      q[1] = 0.25 * S;
      q[2] = (R[1][0] + R[0][1]) / S;
      q[3] = (R[2][0] + R[0][2]) / S;
    } else if (R[1][1] > R[2][2]) {
      double S = std::sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2.0;
      q[0] = (R[2][0] - R[0][2]) / S;
      q[1] = (R[1][0] + R[0][1]) / S;
      q[2] = 0.25 * S;
      q[3] = (R[2][1] + R[1][2]) / S;
    } else {
      double S = std::sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2.0;
      q[0] = (R[0][1] - R[1][0]) / S;
      q[1] = (R[2][0] + R[0][2]) / S;
      q[2] = (R[2][1] + R[1][2]) / S;
      q[3] = 0.25 * S;
    }
  }

  double len = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  if (len > 0.0) {
    q[0] /= len;
    q[1] /= len;
    q[2] /= len;
    q[3] /= len;
  }
}


inline void leftJacobianInverseTransposed(const double Omega[3],
                                          double        JlInv[3][3]) noexcept
{
    const double x = Omega[0], y = Omega[1], z = Omega[2];
    const double theta2 = x*x + y*y + z*z;

    // series expansion for small rotation angles
    if (theta2 < 1e-8)
    {
        double hatO[3][3];
        hat(Omega, hatO);

        const double outer[3][3] = {
            { x*x, x*y, x*z },
            { y*x, y*y, y*z },
            { z*x, z*y, z*z }
        };

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                JlInv[i][j] = (i == j ? 1.0 : 0.0)
                            - 0.5 * hatO[i][j]
                            + (1.0 / 12.0) * outer[i][j];
        return;
    }

    // general case
    const double theta = std::sqrt(theta2);
    const double ux = x / theta, uy = y / theta, uz = z / theta;

    const double half     = 0.5 * theta;
    const double cot_half = std::cos(half) / std::sin(half);

    const double factor_I   = half * cot_half;
    const double factor_uu  = 1.0 - factor_I;
    const double factor_hat = half;

    const double uuT[3][3] = {
        { ux*ux, ux*uy, ux*uz },
        { uy*ux, uy*uy, uy*uz },
        { uz*ux, uz*uy, uz*uz }
    };

    double hatU[3][3];
    const double u[3] = { ux, uy, uz };
    hat(u, hatU);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            JlInv[i][j] = factor_I   * (i == j ? 1.0 : 0.0)
                        + factor_uu  * uuT[i][j]
                        + factor_hat * hatU[i][j];
}

} // namespace so3

#endif
