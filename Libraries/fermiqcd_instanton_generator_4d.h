/////////////////////////////////////////////////////////////////
/// @file fermiqcd_instanton_generator_4d.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains an instanton class
///
/// Distributed under GPL2 license
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_INSTANTON_GENERATOR_4D_
#define FERMIQCD_INSTANTON_GENERATOR_4D_

#include "mdp_matrix.h"
#include "mdp_lattice.h"
#include "fermiqcd_gauge_field.h"
#include "fermiqcd_gamma_matrices.h"

namespace MDP
{
  /**
   * @brief Represents a single (anti-)instanton in 4D Euclidean space.
   *
   * This class stores the position, size (rho), and topological charge
   * of an instanton configuration.
   *
   * @note based on code from James Hetrick
   */
  class SingleInstanton4D
  {
  public:
    mdp_real x[4];   // Position of the instanton in 4D space-time
    mdp_real rho;    // Instanton size (radius parameter)
    mdp_sint charge; // Topological charge (+1 for instanton, -1 for anti-instanton)

    /**
     * @brief Construct a new SingleInstanton4D object.
     *
     * @param x0 Coordinate in direction 0
     * @param x1 Coordinate in direction 1
     * @param x2 Coordinate in direction 2
     * @param x3 Coordinate in direction 3
     * @param rho Instanton size parameter
     * @param charge Topological charge (+1 or -1)
     */
    SingleInstanton4D(mdp_real x0, mdp_real x1, mdp_real x2, mdp_real x3, mdp_real rho, mdp_sint charge)
    {
      x[0] = x0;
      x[1] = x1;
      x[2] = x2;
      x[3] = x3;
      this->rho = rho;
      this->charge = charge;
    }
  };

  /**
   * @brief Generates SU(2) gauge fields from a set of 4D instantons.
   *
   * This class constructs gauge links corresponding to a superposition
   * of instanton and anti-instanton configurations on a 4D lattice.
   */
  class InstantonGenerator4D
  {
  private:
    const mdp_lattice *lattice;
    mdp_matrix tau[4]; // Pauli matrices (tau basis)
    mdp_matrix bar[4]; // Conjugate matrices

    /**
     * @brief Initialize tau and bar matrices.
     *
     * Sets up the SU(2) algebra basis using Pauli matrices.
     */
    void init_tau_and_bar()
    {
      tau[0] = sigma[0]; // indentity
      tau[1] = I * sigma[1];
      tau[2] = I * sigma[2];
      tau[3] = I * sigma[3];
      bar[0] = sigma[0]; // indentity
      bar[1] = -I * sigma[1];
      bar[2] = -I * sigma[2];
      bar[3] = -I * sigma[3];
    }

    /**
     * @brief Compute the singular instanton gauge field A_mu at a point.
     *
     * Implements the continuum instanton solution in singular gauge.
     *
     * @param xl Lattice coordinates (continuous approximation)
     * @param mu Direction index
     * @param instanton Instanton parameters
     *
     * @return mdp_matrix Gauge field contribution A_mu
     */
    mdp_matrix make_singular_instanton(mdp_real xl[4], mdp_suint mu, const SingleInstanton4D &instanton)
    {
      mdp_matrix A;
      mdp_real x[4];

      for (mdp_suint i = 0; i < 4; i++)
      {
        x[i] = std::min(xl[i] - instanton.x[i],
                        (*lattice).size(i) - xl[i] + instanton.x[i]);
      }

      mdp_real x2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];

      mdp_real rho2 = std::pow(instanton.rho, 2);

      if (instanton.charge > 0)
      {
        A = (-0.5 * rho2 / (x2 + rho2) / x2) *
            (bar[mu] * (x[0] * tau[0] + x[1] * tau[1] + x[2] * tau[2] + x[3] * tau[3]) - (x[0] * bar[0] + x[1] * bar[1] + x[2] * bar[2] + x[3] * bar[3]) * tau[mu]);
      }
      else
      {
        A = (-0.5 * rho2 / (x2 + rho2) / x2) *
            (tau[mu] * (x[0] * bar[0] + x[1] * bar[1] + x[2] * bar[2] + x[3] * bar[3]) - (x[0] * tau[0] + x[1] * tau[1] + x[2] * tau[2] + x[3] * tau[3]) * bar[mu]);
      }

      return A;
    }

    /**
     * @brief Construct an SU(2) link variable from instanton fields.
     *
     * Performs a path-ordered exponential of the gauge field along
     * a lattice link using a simple discretization.
     *
     * @param xn Lattice site
     * @param mu Direction index
     * @param instantons List of instantons
     *
     * @return mdp_matrix SU(2) link matrix
     */
    mdp_matrix make_su2_link(mdp_site xn, mdp_suint mu, const std::vector<SingleInstanton4D> &instantons)
    {
      mdp_matrix A(2, 2);
      mdp_matrix P = sigma[0];
      mdp_real x[4], a[3], norm_a;
      constexpr mdp_suint steps = 5;
      mdp_real dx = 1.0 / steps;

      for (mdp_suint i = 0; i < 4; i++)
      {
        x[i] = xn(i);
      }

      for (mdp_suint m = 0; m < steps; m++)
      {
        x[mu] = (0.5 + m) / steps + xn(mu);

        for (size_t k = 0; k < instantons.size(); k++)
        {
          A += make_singular_instanton(x, mu, instantons[k]);
        }

        a[0] = imag(A(0, 1)) * dx;
        a[1] = real(A(0, 1)) * dx;
        a[2] = imag(A(0, 0)) * dx;

        norm_a = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

        if (norm_a > mdp_precision)
        {
          P *= std::sin(norm_a) / norm_a * (a[0] * tau[1] + a[1] * tau[2] + a[2] * tau[3]) + std::cos(norm_a) * tau[0];
        }
      }

      return P;
    }

  public:
    /**
     * @brief Generate a gauge field configuration from instantons.
     *
     * Fills the lattice gauge field U with SU(2) link variables derived
     * from the provided instanton ensemble.
     *
     * @param U Gauge field to be populated
     * @param instantons Vector of instanton configurations
     */
    void generate(gauge_field &U, const std::vector<SingleInstanton4D> &instantons)
    {
      init_tau_and_bar();
      lattice = &U.lattice();
      mdp_site x(U.lattice());
      mdp_matrix A;

      if (U.ndim() != 4)
        throw std::string("instantons only in 4D");

      forallsites(x)
      {
        for (mdp_suint mu = 0; mu < U.ndim(); mu++)
        {
          A = make_su2_link(x, mu, instantons);
          U(x, mu) = 1;
          for (mdp_suint i = 0; i < 2; i++)
          {
            for (mdp_suint j = 0; j < 2; j++)
            {
              U(x, mu, i, j) = A(i, j);
            }
          }
        }
      }
    }
  };
} // namespace MDP

#endif /* FERMIQCD_INSTANTON_GENERATOR_4D_ */
