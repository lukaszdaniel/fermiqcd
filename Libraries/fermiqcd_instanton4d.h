/////////////////////////////////////////////////////////////////
/// @file fermiqcd_instanton4d.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// class for building a single instanton gauge configuration
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_INSTANTON4D_
#define FERMIQCD_INSTANTON4D_

#include <vector>
#include "mdp_global_vars.h"
#include "mdp_epsilon.h"
#include "mdp_matrix.h"
#include "mdp_site.h"
#include "fermiqcd_gamma_matrices.h"

namespace MDP
{
  /**
   * @brief SU(N) instanton gauge field in 4D Euclidean space.
   *
   * This class constructs a continuum instanton solution embedded in an SU(N)
   * gauge group. The solution is built using generalized 't Hooft symbols
   * and an SU(2) subgroup embedded into SU(N).
   *
   * @details
   * The gauge field has the form:
   * \f[
   * A_\mu(x) = \frac{2\,q}{(x - x_0)^2 + \lambda^2} \, \eta_{\mu\nu} (x - x_0)_\nu
   * \f]
   * where:
   * - \f$q\f$ is proportional to the inverse coupling (stored in charge),
   * - \f$\lambda\f$ is the instanton size,
   * - \f$\eta_{\mu\nu}\f$ are SU(N)-embedded 't Hooft symbols,
   * - the field is centered at position \f$x_0\f$.
   *
   * The construction embeds an SU(2) instanton into SU(N) by selecting
   * two color indices (sub_i, sub_j).
   */
  class Instanton4D
  {
  public:
    std::array<mdp_real, 4> p; // Instanton center position in 4D
    mdp_suint nc;              // Number of colors (SU(N))
    mdp_suint sub_i, sub_j;    // Indices defining SU(2) subgroup embedding
    mdp_real charge;           // Charge parameter (proportional to 1/g)
    mdp_real lambda;           // Instanton size (radius parameter)
    mdp_matrix eta[4][4];      // Embedded 't Hooft symbols

    /**
     * @brief Construct an SU(N) instanton configuration.
     *
     * @param nc_ Number of colors
     * @param sub_i_ First index of SU(2) subgroup embedding
     * @param sub_j_ Second index of SU(2) subgroup embedding
     * @param charge_ Charge parameter (typically 1/g)
     * @param lambda_ Instanton size parameter
     * @param p_ 4D position of the instanton center
     *
     * @details
     * This constructor builds the embedded 't Hooft symbols eta[mu][nu]
     * using rotated Pauli matrices. The SU(2) structure is inserted into
     * the full SU(N) matrix by acting only on the (sub_i, sub_j) subspace.
     */
    Instanton4D(const std::array<mdp_real, 4> &pos, mdp_real lambda_, mdp_real charge_, mdp_suint nc_, mdp_suint sub_i_, mdp_suint sub_j_,
                mdp_real alpha = 0.0, mdp_real beta = 0.0, mdp_real gamma = 0.0)
    {
      nc = nc_;
      sub_i = sub_i_;
      sub_j = sub_j_;
      lambda = lambda_;
      charge = charge_;
      p = pos;

      mdp_matrix U(2, 2);

      U(0, 0) = std::cos(beta / 2.0) * exp(I * (alpha + gamma) / 2.0);
      U(0, 1) = std::sin(beta / 2.0) * exp(I * (alpha - gamma) / 2.0);
      U(1, 0) = -std::sin(beta / 2.0) * exp(-I * (alpha - gamma) / 2.0);
      U(1, 1) = std::cos(beta / 2.0) * exp(-I * (alpha + gamma) / 2.0);

      mdp_matrix sigma_rot[4];
      for (mdp_uint a = 1; a < 4; a++)
      {
        sigma_rot[a] = U * sigma[a] * hermitian(U);
      }
#if 1 // build eta using matrix operations
      mdp_matrix T(nc, nc);

      // Construct embedded 't Hooft symbols
      for (mdp_suint mu = 0; mu < 4; mu++)
      {
        for (mdp_suint nu = 0; nu < 4; nu++)
        {
          eta[mu][nu].dimension(nc, nc);
          eta[mu][nu] = 0;

          for (mdp_uint a = 1; a < 4; a++)
          {
            T = 0;

            // Embed SU(2) generator into SU(N)
            T(sub_i, sub_i) = sigma_rot[a](0, 0);
            T(sub_i, sub_j) = sigma_rot[a](0, 1);
            T(sub_j, sub_i) = sigma_rot[a](1, 0);
            T(sub_j, sub_j) = sigma_rot[a](1, 1);

            if (a > 0 && mu > 0 && nu > 0)
            {
              eta[mu][nu] += epsilon(a, mu, nu) * T;
            }
            else if (mu == 0 && a > 0 && a == nu)
            {
              eta[mu][nu] -= T; // T[a]*(-1)*delta(a,nu)
            }
            else if (nu == 0 && a > 0 && a == mu)
            {
              eta[mu][nu] += T; // T[a]*delta(a,mu)
            }
          }
        }
      }
#else // build eta element by element
      for (mdp_suint mu = 0; mu < 4; mu++)
      {
        for (mdp_suint nu = 0; nu < 4; nu++)
        {
          eta[mu][nu].dimension(nc, nc);
          eta[mu][nu] = 0;

          mdp_matrix block(2, 2);
          block = 0;

          for (mdp_suint a = 1; a < 4; a++)
          {
            mdp_real coeff = 0.0;

            if (mu == 0 && nu == 0)
            {
              coeff = 0.0;
            }
            else if (mu == 0)
            {
              coeff = (a == nu) ? -1.0 : 0.0;
            }
            else if (nu == 0)
            {
              coeff = (a == mu) ? +1.0 : 0.0;
            }
            else
            {
              coeff = epsilon(a, mu, nu);
            }

            block += coeff * sigma_rot[a];
          }

          eta[mu][nu](sub_i, sub_i) = block(0, 0);
          eta[mu][nu](sub_i, sub_j) = block(0, 1);
          eta[mu][nu](sub_j, sub_i) = block(1, 0);
          eta[mu][nu](sub_j, sub_j) = block(1, 1);
        }
      }
#endif
    }

    /**
     * @brief Evaluate the gauge field A_mu(x).
     *
     * @param x Lattice site
     * @param mu Direction index
     *
     * @return mdp_matrix Gauge field matrix A_mu(x)
     *
     * @details
     * Computes the continuum instanton gauge field at lattice site x:
     * \f[
     * A_\mu(x) = \frac{2\,q}{(x - x_0)^2 + \lambda^2}
     * \sum_\nu \eta_{\mu\nu} (x - x_0)_\nu
     * \f]
     */
    mdp_matrix operator()(mdp_site &x, mdp_suint mu) const
    {
      mdp_int v[4];

      // Displacement from instanton center
      for (mdp_suint nu = 0; nu < 4; nu++)
        v[nu] = x(nu) - p[nu];

      // Squared distance
      mdp_real d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];

      mdp_matrix A(nc, nc);
      A = 0;

      // Linear combination with 't Hooft symbols
      for (mdp_suint nu = 0; nu < 4; nu++)
      {
        A += eta[mu][nu] * v[nu];
      }

      return (2.0 * charge / (d2 + lambda * lambda)) * A;
    }
  };
} // namespace MDP

#endif /* FERMIQCD_INSTANTON4D_ */
