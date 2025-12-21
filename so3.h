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

#ifndef LMP_SO3_H
#define LMP_SO3_H

#include <cmath>
#include <algorithm>  // for std::max
#include "lamath.h"

namespace so3 {


// “Hat” map: v → [v]_×
inline void hat(const double v[3], double M[3][3]) noexcept {
    M[0][0] =   0.0;   M[0][1] = -v[2];  M[0][2] =  v[1];
    M[1][0] =  v[2];   M[1][1] =   0.0;  M[1][2] = -v[0];
    M[2][0] = -v[1];   M[2][1] =  v[0];  M[2][2] =   0.0;
}


// -----------------------------------------------------------------
// Ts2R:  Compute R = Ti^T * Tj  (two triads → relative rotation)
// Ti, Tj: double[3][3] (row-major each is an orthonormal triad)
// R:      double[3][3] (output relative rotation)
// -----------------------------------------------------------------
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

// -----------------------------------------------------------------
// euler2rotmat:  Rodrigues’ formula (Eq. A1)
//   Input:   omega[3] = axis-angle vector  (||omega|| = θ)
//   Output:  R[3][3] rotation matrix
// -----------------------------------------------------------------
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

// -----------------------------------------------------------------
// rotmat2euler:  Logarithm map (Eq. A2)
//   Input:   R[3][3] rotation matrix
//   Output:  omega[3] = axis-angle vector (||omega|| = θ ∈ [0,π])
// -----------------------------------------------------------------
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

// -----------------------------------------------------------------
// triads2euler: Compute Omega = vec(log(S^T * Ti^T * Tj))
// Inputs:  Ti, Tj, S ∈ SO(3) (each as double[3][3])
// Output: Omega[3] = axis-angle of  D = S^T * Ti^T * Tj
// -----------------------------------------------------------------
inline void triads2euler(const double Ti[3][3],
                         const double Tj[3][3],
                         const double S[3][3],
                         double       Omega[3]) noexcept
{
    double R[3][3], St[3][3], D[3][3];

    // 1) R = Ti^T * Tj
    triads2rotmat(Ti, Tj, R);

    // 2) St = S^T  (transpose of S)
    St[0][0] = S[0][0];  St[0][1] = S[1][0];  St[0][2] = S[2][0];
    St[1][0] = S[0][1];  St[1][1] = S[1][1];  St[1][2] = S[2][1];
    St[2][0] = S[0][2];  St[2][1] = S[1][2];  St[2][2] = S[2][2];

    // 3) D = St * R
    lamath::mul(St, R, D);

    // 4) Omega = rotmat2euler(D)
    rotmat2euler(D, Omega);
}


// // -----------------------------------------------------------------
// // Quaternion to Rotation matrix
// // -----------------------------------------------------------------

// inline void quat2rotmat(
//     const double *quat, 
//     double mat[3][3]) noexcept
// {
//   double w2 = quat[0]*quat[0];
//   double i2 = quat[1]*quat[1];
//   double j2 = quat[2]*quat[2];
//   double k2 = quat[3]*quat[3];
//   double twoij = 2.0*quat[1]*quat[2];
//   double twoik = 2.0*quat[1]*quat[3];
//   double twojk = 2.0*quat[2]*quat[3];
//   double twoiw = 2.0*quat[1]*quat[0];
//   double twojw = 2.0*quat[2]*quat[0];
//   double twokw = 2.0*quat[3]*quat[0];

//   mat[0][0] = w2+i2-j2-k2;
//   mat[1][0] = twoij-twokw;
//   mat[2][0] = twojw+twoik;

//   mat[0][1] = twoij+twokw;
//   mat[1][1] = w2-i2+j2-k2;
//   mat[2][1] = twojk-twoiw;

//   mat[0][2] = twoik-twojw;
//   mat[1][2] = twojk+twoiw;
//   mat[2][2] = w2-i2-j2+k2;
// }

// -----------------------------------------------------------------
// Rotation matrix to Quaternion 
// -----------------------------------------------------------------

inline void rotmat2quat(const double R[3][3], double q[4]) noexcept {

  // Compute the trace of R:
  double trace = R[0][0] + R[1][1] + R[2][2];

  if (trace > 0.0) {
    // If trace is positive, use the “trace” formula:
    double S = std::sqrt(trace + 1.0) * 2.0; // S = 4*qw
    q[0] = 0.25 * S;
    q[1] = (R[1][2] - R[2][1]) / S;
    q[2] = (R[2][0] - R[0][2]) / S;
    q[3] = (R[0][1] - R[1][0]) / S;
  } else {
    // Otherwise, find the largest diagonal element
    if ((R[0][0] > R[1][1]) && (R[0][0] > R[2][2])) {
      // R[0][0] is the largest
      double S = std::sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2.0; // S = 4*qx
      q[0] = (R[1][2] - R[2][1]) / S;
      q[1] = 0.25 * S;
      q[2] = (R[1][0] + R[0][1]) / S;
      q[3] = (R[2][0] + R[0][2]) / S;
    } else if (R[1][1] > R[2][2]) {
      // R[1][1] is the largest
      double S = std::sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2.0; // S = 4*qy
      q[0] = (R[2][0] - R[0][2]) / S;
      q[1] = (R[1][0] + R[0][1]) / S;
      q[2] = 0.25 * S;
      q[3] = (R[2][1] + R[1][2]) / S;
    } else {
      // R[2][2] is the largest
      double S = std::sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2.0; // S = 4*qz
      q[0] = (R[0][1] - R[1][0]) / S;
      q[1] = (R[2][0] + R[0][2]) / S;
      q[2] = (R[2][1] + R[1][2]) / S;
      q[3] = 0.25 * S;
    }
  }

  // Optionally, normalize to guard against numerical drift:
  double len = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  if (len > 0.0) {
    q[0] /= len;
    q[1] /= len;
    q[2] /= len;
    q[3] /= len;
  }
}



// -----------------------------------------------------------------
// leftJacobian: compute J_l(Omega) for SO(3)
//   Input:   Omega[3]  (axis-angle vector, ||Omega|| = theta)
//   Output:  Jl[3][3]  (left-Jacobian matrix)
// -----------------------------------------------------------------
inline void leftJacobian(const double Omega[3], double Jl[3][3]) noexcept {
    // Compute theta = ||Omega||
    double x = Omega[0], y = Omega[1], z = Omega[2];
    double theta2 = x*x + y*y + z*z;

    if (theta2 < 1e-8) {
        // Use series for small theta:
        // Jl ≈ I - (1/2)[Ω]_× + (1/6)(Ω Ω^T)
        double hatO[3][3];
        hat(Omega, hatO);
        double outer[3][3] = {
            { x*x,   x*y,   x*z },
            { y*x,   y*y,   y*z },
            { z*x,   z*y,   z*z }
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Jl[i][j] = (i==j ? 1.0 : 0.0)
                         - 0.5 * hatO[i][j]
                         + (1.0/6.0) * outer[i][j];
            }
        }
        return;
    }

    double theta = std::sqrt(theta2);
    double ux = x / theta;
    double uy = y / theta;
    double uz = z / theta;

    double sin_t = std::sin(theta);
    double cos_t = std::cos(theta);
    double s_over_t = sin_t / theta;
    double one_minus_s_over_t = 1.0 - s_over_t;
    double one_minus_c_over_t = (1.0 - cos_t) / theta;

    // Build u u^T
    double uuT[3][3] = {
        { ux*ux,  ux*uy,  ux*uz },
        { uy*ux,  uy*uy,  uy*uz },
        { uz*ux,  uz*uy,  uz*uz }
    };

    // Build [u]_×
    double hatU[3][3];
    double u[3] = { ux, uy, uz };
    hat(u, hatU);

    // Jl = (sinθ/θ) I + (1 - sinθ/θ) (u u^T) + ((1 - cosθ)/θ) [u]_×
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Jl[i][j] = s_over_t * (i==j ? 1.0 : 0.0)
                     + one_minus_s_over_t * uuT[i][j]
                     + one_minus_c_over_t * hatU[i][j];
        }
    }
}

// -----------------------------------------------------------------
// leftJacobianInverse: compute J_l^{-1}(Omega) for SO(3)
//   Input:   Omega[3]  (axis-angle vector, ||Omega|| = theta)
//   Output:  JlInv[3][3]  (inverse left-Jacobian)
// -----------------------------------------------------------------
inline void leftJacobianInverse(const double Omega[3], double JlInv[3][3]) noexcept {
    // Compute theta = ||Omega||
    double x = Omega[0], y = Omega[1], z = Omega[2];
    double theta2 = x*x + y*y + z*z;

    if (theta2 < 1e-8) {
        // Use series for small theta:
        // JlInv ≈ I + (1/2)[Ω]_× + (1/12)(Ω Ω^T)
        double hatO[3][3];
        hat(Omega, hatO);
        double outer[3][3] = {
            { x*x,   x*y,   x*z },
            { y*x,   y*y,   y*z },
            { z*x,   z*y,   z*z }
        };
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                JlInv[i][j] = (i==j ? 1.0 : 0.0)
                            + 0.5 * hatO[i][j]
                            + (1.0/12.0) * outer[i][j];
            }
        }
        return;
    }

    double theta = std::sqrt(theta2);
    double ux = x / theta;
    double uy = y / theta;
    double uz = z / theta;

    double half = 0.5 * theta;
    double cot_half = std::cos(half) / std::sin(half);

    // Build u u^T
    double uuT[3][3] = {
        { ux*ux,  ux*uy,  ux*uz },
        { uy*ux,  uy*uy,  uy*uz },
        { uz*ux,  uz*uy,  uz*uz }
    };

    // Build [u]_×
    double hatU[3][3];
    double u[3] = { ux, uy, uz };
    hat(u, hatU);

    // JlInv = (θ/2)cot(θ/2) * I
    //         + [1 - (θ/2)cot(θ/2)] * (u u^T)
    //         - (θ/2) [u]_×
    double factor_I = half * cot_half;
    double factor_uu = 1.0 - factor_I;
    double factor_hat = half;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            JlInv[i][j] = factor_I * (i==j ? 1.0 : 0.0)
                        + factor_uu * uuT[i][j]
                        - factor_hat * hatU[i][j];
        }
    }
}


// -----------------------------------------------------------------
// leftJacobianInverseTransposed : compute (J_l^{-1}(Ω))ᵀ for SO(3)
//   Input : Omega[3]  (axis-angle vector, ‖Ω‖ = θ)
//   Output: JlInv[3][3]  ((inverse left-Jacobian)ᵀ)
// -----------------------------------------------------------------
inline void leftJacobianInverseTransposed(const double Omega[3],
                                          double        JlInv[3][3]) noexcept
{
    const double x = Omega[0], y = Omega[1], z = Omega[2];
    const double theta2 = x*x + y*y + z*z;

    /* ----------  tiny angles: use series expansion  ---------- */
    if (theta2 < 1e-8)
    {
        // (J_l^{-1})ᵀ ≈ I − ½[Ω]_× + 1⁄12 ΩΩᵀ
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

    /* ----------  general case  ---------- */
    const double theta = std::sqrt(theta2);
    const double ux = x / theta, uy = y / theta, uz = z / theta;

    const double half     = 0.5 * theta;
    const double cot_half = std::cos(half) / std::sin(half);

    const double factor_I   = half * cot_half;
    const double factor_uu  = 1.0 - factor_I;
    const double factor_hat = half;                 // +θ/2 for the transpose!

    /* u uᵀ */
    const double uuT[3][3] = {
        { ux*ux, ux*uy, ux*uz },
        { uy*ux, uy*uy, uy*uz },
        { uz*ux, uz*uy, uz*uz }
    };

    /* [u]× */
    double hatU[3][3];
    const double u[3] = { ux, uy, uz };
    hat(u, hatU);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            JlInv[i][j] = factor_I   * (i == j ? 1.0 : 0.0)
                        + factor_uu  * uuT[i][j]
                        + factor_hat * hatU[i][j];
}


inline void rightJacobianInverse(const double Omega[3],
                                 double        JrInv[3][3]) noexcept
{
    // Calculate the squared magnitude of the rotation vector
    const double theta_sq = Omega[0] * Omega[0] + Omega[1] * Omega[1] + Omega[2] * Omega[2];

    // Define a small tolerance for detecting near-zero rotation
    constexpr double epsilon = 1e-9;

    // Case 1: The rotation angle is very small.
    // The Jacobian is the identity matrix.
    if (theta_sq < epsilon) {
        JrInv[0][0] = 1.0; JrInv[0][1] = 0.0; JrInv[0][2] = 0.0;
        JrInv[1][0] = 0.0; JrInv[1][1] = 1.0; JrInv[1][2] = 0.0;
        JrInv[2][0] = 0.0; JrInv[2][1] = 0.0; JrInv[2][2] = 1.0;
        return;
    }

    // Calculate the rotation angle (magnitude of the vector)
    const double theta = std::sqrt(theta_sq);

    // Construct the skew-symmetric cross-product matrix [Omega]x
    // Using K for the skew-symmetric matrix as is common.
    const double K[3][3] = {
        {0.0,      -Omega[2],  Omega[1]},
        {Omega[2],  0.0,      -Omega[0]},
        {-Omega[1], Omega[0],  0.0     }
    };

    // Calculate the square of the skew-symmetric matrix, K^2
    double K_sq[3][3] = {{0.0}};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                K_sq[i][j] += K[i][k] * K[k][j];
            }
        }
    }

    // --- Calculate the scalar coefficients for the formula ---
    // Jr_inv = I + (1/2)*K + C*K^2

    // Coefficient for K
    constexpr double B = 0.5;

    // Numerically stable coefficient for K^2
    // C = (1/theta^2) - ( (1 + cos(theta)) / (2*theta*sin(theta)) )
    double C;
    // Use Taylor series for C when theta is small to avoid floating point errors
    if (std::abs(theta) < 1e-5) {
        C = 1.0/12.0 + theta_sq/720.0 + theta_sq*theta_sq/30240.0;
    } else {
        C = (1.0 / theta_sq) * (1.0 - (theta * std::sin(theta)) / (2.0 * (1.0 - std::cos(theta))));
    }

    // --- Combine the terms to form the final matrix ---
    // JrInv = I + B*K + C*K^2
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            JrInv[i][j] = B * K[i][j] + C * K_sq[i][j];
        }
        JrInv[i][i] += 1.0; // Add the identity component
    }
}




} // namespace so3

#endif
