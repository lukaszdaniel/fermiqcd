/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_fixing.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Gauge fixing stuff
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_GAUGE_FIXING_
#define FERMIQCD_GAUGE_FIXING_

#include "fermiqcd_gauge_field.h"

namespace MDP
{
  /**
   * @brief Structure storing statistics from the gauge fixing procedure.
   *
   * This structure collects diagnostic information about the iterative gauge
   * fixing algorithm, including convergence behaviour and final physical
   * observables estimated during the procedure.
   */
  class gaugefixing_stats
  {
  public:
    /// Maximum number of gauge fixing iterations allowed.
    mdp_uint max_steps;

    /// Target precision threshold for convergence.
    mdp_real target_precision;

    /// Number of iterations actually performed.
    mdp_uint steps;

    /// Final measured gauge fixing precision.
    mdp_real precision;

    /// Final value of the diagnostic gauge "action".
    mdp_real action;
  };

  /**
   * @brief Collection of gauge fixing algorithms for lattice gauge fields.
   *
   * This class implements iterative gauge fixing procedures for lattice gauge
   * theories, including Coulomb and Landau gauges. The algorithms are based on
   * SU(2) subgroup overrelaxation (Cabibbo–Marinari scheme) applied to SU(N)
   * gauge fields.
   *
   * The class is non-instantiable and acts purely as a namespace for static
   * algorithms.
   *
   * @example
   * @verbatim
   *    gauge_field U(lattice,nc);
   *    gaugefixing_stats stats;
   *    U.load("myfield");
   *    stats = GaugeFixing::fix(U, GaugeFixing::Coulomb, 100);
   *    U.save("myfield_gaugefixed");
   * @endverbatim
   */
  class GaugeFixing
  {
  private:
    GaugeFixing() = delete;
    ~GaugeFixing() = delete;
    GaugeFixing(const GaugeFixing &) = delete;
    GaugeFixing &operator=(const GaugeFixing &) = delete;

  public:
    /// Coulomb gauge fixing mode identifier.
    static constexpr mdp_suint Coulomb = 0;

    /// Landau gauge fixing mode identifier.
    static constexpr mdp_suint Landau = 10;

    /**
     * @brief Applies a single SU(2) subgroup gauge fixing update (overrelaxation step).
     *
     * This function performs a local gauge transformation on a selected SU(2)
     * subgroup (i, j) embedded in SU(N). The transformation is constructed from
     * local staples and applied in an overrelaxation-like manner.
     *
     * The update acts differently on even/odd lattice sites to preserve
     * consistency of link updates:
     * - forward links: left multiplication
     * - backward links: right multiplication with Hermitian conjugate
     *
     * @param U Gauge field to be updated.
     * @param mu Direction of the link being updated.
     * @param parity Lattice parity (EVEN or ODD) of updated sites.
     * @param i First color index of SU(2) subgroup.
     * @param j Second color index of SU(2) subgroup.
     * @param overrelaxation_boost Optional parameter controlling relaxation strength.
     *
     * @note This is a local update and does not itself ensure global convergence.
     */
    static void hit(gauge_field &U,
                    mdp_suint mu,
                    mdp_parity parity,
                    mdp_suint i, mdp_suint j,
                    mdp_real overrelaxation_boost = 1)
    {

      // This also works for twisted boundary even if I use mdp_field<> */

      static mdp_real a0, a1, a2, a3, b, c, d;
      static mdp_real a0_sq, ai_sq;
      static mdp_complex x0, x1;
      mdp_suint nc = U.nc();
      mdp_parity opposite_parity = EVENODD;
      mdp_complex_vector_field W(U.lattice(), 4);
      mdp_matrix U_up(nc, nc), U_dw(nc, nc);
      mdp_matrix A;
      mdp_site x(U.lattice());
      mdp_site y(U.lattice());

      switch (parity)
      {
      case EVEN:
        opposite_parity = ODD;
        break;
      case ODD:
        opposite_parity = EVEN;
        break;
      default:
        break;
      }

      forallsitesofparity(x, parity)
      {
        a0 = a1 = a2 = a3 = 0;
        for (mdp_suint nu = 0; nu < U.ndim(); nu++)
          if (nu != mu)
          {
            U_up = U(x, nu);
            U_dw = U(x - nu, nu);
            a0 +=
                +real(U_dw(i, i)) + real(U_up(i, i)) + real(U_dw(j, j)) + real(U_up(j, j));
            a1 +=
                +imag(U_dw(j, i)) - imag(U_up(j, i)) + imag(U_dw(i, j)) - imag(U_up(i, j));
            a2 +=
                -real(U_dw(j, i)) + real(U_up(j, i)) + real(U_dw(i, j)) - real(U_up(i, j));
            a3 +=
                +imag(U_dw(i, i)) - imag(U_up(i, i)) - imag(U_dw(j, j)) + imag(U_up(j, j));
          }
        ai_sq = a1 * a1 + a2 * a2 + a3 * a3;
        a0_sq = a0 * a0;
        b = (overrelaxation_boost * a0_sq + ai_sq) / (a0_sq + ai_sq);
        c = std::sqrt(a0_sq + b * b * ai_sq);
        d = b / c;
        a0 /= c;
        a1 *= d;
        a2 *= d;
        a3 *= d;
        W(x, 0) = mdp_complex(a0, a3);
        W(x, 1) = mdp_complex(a2, a1);
        W(x, 2) = mdp_complex(-a2, a1);
        W(x, 3) = mdp_complex(a0, -a3);
      }
      W.update();
      forallsitesofparity(x, parity)
      {
        for (mdp_suint nu = 0; nu < U.ndim(); nu++)
          for (mdp_suint k = 0; k < U.nc(); k++)
          {
            x0 = U(x, nu, i, k);
            x1 = U(x, nu, j, k);
            U(x, nu, i, k) = W(x, 0) * x0 + W(x, 1) * x1;
            U(x, nu, j, k) = W(x, 2) * x0 + W(x, 3) * x1;
          }
      }
      forallsitesofparity(x, opposite_parity)
      {
        for (mdp_suint nu = 0; nu < U.ndim(); nu++)
        {
          y = x + nu;
#ifndef TWISTED_BOUNDARY
          for (mdp_suint k = 0; k < U.nc(); k++)
          {
            x0 = U(x, nu, k, i);
            x1 = U(x, nu, k, j);
            U(x, nu, k, i) = conj(W(y, 0)) * x0 + conj(W(y, 1)) * x1;
            U(x, nu, k, j) = conj(W(y, 2)) * x0 + conj(W(y, 3)) * x1;
          }
#else
          if (in_block(y))
            for (mdp_suint k = 0; k < U.nc(); k++)
            {
              x0 = U(x, nu, k, i);
              x1 = U(x, nu, k, j);
              U(x, nu, k, i) = conj(W(y, 0)) * x0 + conj(W(y, 1)) * x1;
              U(x, nu, k, j) = conj(W(y, 2)) * x0 + conj(W(y, 3)) * x1;
            }
          else
          {
            A = mdp_identity(nc);
            A(i, i) = W(y, 0);
            A(i, j) = W(y, 1);
            A(j, i) = W(y, 2);
            A(j, j) = W(y, 3);
            twist_boundary(A, y);
            U(x, nu) = U(x, nu) * hermitian(A);
          }
#endif
        }
      }
      U.update();
    }

    /**
     * @brief Removes residual Z(N) center symmetry on a given lattice direction.
     *
     * This function fixes the global Z(3) (or Z(N)) ambiguity present in lattice
     * gauge configurations due to toroidal boundary conditions. It scans lattice
     * slices and applies center-phase rotations that maximize the trace alignment.
     *
     * @param U Gauge field to be center-fixed.
     * @param mu Direction along which Z(N) fixing is performed.
     *
     * @note This step is optional and typically applied after gauge fixing
     *       convergence has been reached.
     */
    static void z3_fix(gauge_field &U, mdp_suint mu)
    {
      mdp_suint i = 0;
      mdp_site x(U.lattice());
      mdp_matrix A;
      mdp_complex phase[3] = {mdp_complex(1, 0),
                              exp(2.0 * Pi * I / 3.0),
                              exp(4.0 * Pi * I / 3.0)};
      mdp_real alpha[3];
      for (mdp_uint t = 0; t < U.lattice().size(mu) - 1; t++)
      {
        A = mdp_zero(U.nc());
        forallsites(x)
        {
          if (x(mu) == t)
            A += U(x, mu);
        }
        mdp.add(A);
        alpha[0] = real(trace(A * conj(phase[0])));
        alpha[1] = real(trace(A * conj(phase[1])));
        alpha[2] = real(trace(A * conj(phase[2])));
        if (alpha[0] >= alpha[1] && alpha[0] >= alpha[2])
          i = 0;
        if (alpha[1] > alpha[0] && alpha[1] >= alpha[2])
          i = 1;
        if (alpha[2] > alpha[0] && alpha[2] >= alpha[1])
          i = 2;
        forallsites(x)
        {
          if (x(mu) == t)
            U(x, mu) *= conj(phase[i]);
          else if (x(mu) == t + 1)
            U(x, mu) *= phase[i];
        }
        U.update();
      }
    }

    /**
     * @brief Iterative gauge fixing procedure (Landau or Coulomb gauge).
     *
     * This function performs iterative gauge fixing of a lattice gauge field
     * using SU(2) subgroup overrelaxation (Cabibbo–Marinari method). The
     * algorithm alternates updates over even/odd sublattices and SU(2)
     * embeddings until convergence is reached or a maximum number of steps
     * is exceeded.
     *
     * The procedure minimizes a gauge-fixing functional (depending on the
     * chosen gauge type) by iteratively applying local gauge transformations.
     *
     * Convergence is monitored using a lattice-wide "precision" measure based
     * on deviations from the gauge condition.
     *
     * Optionally, a Z(N) center fixing step can be applied after convergence.
     *
     * @param U Gauge field to be gauge fixed (modified in place).
     * @param mu Gauge direction or gauge type selector
     *           (e.g. GaugeFixing::Landau or GaugeFixing::Coulomb).
     * @param max_steps Maximum number of gauge fixing iterations.
     * @param target_precision Convergence threshold for stopping criterion.
     * @param overrelaxation_boost Controls strength of overrelaxation updates.
     * @param z3 If true, applies Z(3) center fixing after convergence.
     *
     * @return gaugefixing_stats Structure containing convergence information:
     *         number of steps, final precision, and diagnostic action.
     *
     * @note The "action" reported here is a diagnostic quantity and is not
     *       necessarily the Wilson gauge action.
     *
     * @warning Convergence depends on initial configuration and may stall if
     *          the field is far from satisfying the gauge condition.
     */
    static gaugefixing_stats fix(gauge_field &U,
                                 mdp_suint mu = 0,
                                 mdp_uint max_steps = 1,
                                 mdp_real target_precision = 1e-5,
                                 mdp_real overrelaxation_boost = 1,
                                 bool z3 = false)
    {
      gaugefixing_stats stats;
      mdp_uint step = 0;
      mdp_site x(U.lattice());
      mdp_real action = 0;
      mdp_real precision = 0;
      mdp_matrix M(U.nc(), U.nc());

      stats.max_steps = max_steps;
      stats.target_precision = target_precision;

      mdp << "step\taction\tprecision\n";

      for (step = 0; step < max_steps; step++)
      {

        for (mdp_parity parity : {EVEN, ODD})
        {
          for (mdp_suint i = 0; i < U.nc() - 1; i++)
            for (mdp_suint j = i + 1; j < U.nc(); j++)
            {
              hit(U, mu, parity, i, j, overrelaxation_boost);
            }
        }
        // ******** What is the gauge action ? *********
        action = 0;
        precision = 0;
        forallsites(x)
        {
          M = 0;
          for (mdp_suint nu = 0; nu < U.ndim(); nu++)
            if (nu != mu)
            {
              M += U(x, nu) - U(x - nu, nu);
              action += real(trace(U(x, nu) + U(x - nu, nu)));
            }
          M = (M - trace(M) * (1.0 / U.nc()));
          M = M - hermitian(M);
          for (mdp_suint i = 0; i < U.nc(); i++)
            for (mdp_suint j = 0; j < U.nc(); j++)
              precision += std::pow(abs(M(i, j)), 2);
        }
        mdp.add(precision);
        mdp.add(action);

        precision = std::sqrt(precision / (U.nc() * U.nc() * U.lattice().global_volume()));
        action = action / (2.0 * U.nc() * U.lattice().global_volume());

        mdp << step << "\t" << action << "\t" << precision << "\n";

        if (step != 0 && precision < target_precision)
          break;
        // ***********************************************
      }
      stats.steps = step;
      stats.precision = precision;
      stats.action = action;

      if (z3)
        z3_fix(U, mu);

      mdp << "steps = " << step
          << ", action = " << action
          << ", precision = " << precision << "\n";

      return stats;
    }
  };
} // namespace MDP

#endif /* FERMIQCD_GAUGE_FIXING_ */
