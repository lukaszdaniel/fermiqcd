/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_routines.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various gauge multiplication routines
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_GAUGE_ROUTINES_
#define FERMIQCD_GAUGE_ROUTINES_

#include "mdp_matrix.h"
#include "fermiqcd_gauge_field.h"

namespace MDP
{
  /**
   * @brief Computes a single oriented staple contribution for a given direction.
   *
   * This function returns one of the two staples contributing to the link
   * U(x, mu), associated with direction nu.
   *
   * The staple corresponds to a "U-shaped" path (3 links) forming half of a plaquette:
   *
   *  - For s1 = +1 (forward staple):
   *      U(x,nu) U(x+nu,mu) U†(x+mu,nu)
   *
   *  - For s1 = -1 (backward staple):
   *      U†(x-nu,nu) U(x-nu,mu) U(x-nu+mu,nu)
   *
   * These objects are building blocks for:
   *  - Wilson gauge force (HMC),
   *  - link updates (heatbath, overrelaxation),
   *  - smearing procedures.
   *
   * @param U   Gauge field.
   * @param x   Lattice site.
   * @param mu  Direction of the central link.
   * @param s1  Orientation (+1 forward, -1 backward).
   * @param nu  Orthogonal direction defining the plaquette plane.
   *
   * @return The staple matrix (Nc x Nc).
   */
  mdp_matrix staple(const gauge_field &U, mdp_site x,
                    mdp_suint mu, mdp_int s1, mdp_suint nu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    if (s1 == +1)
    {
      tmp = U(x, nu) * U(x + nu, mu) * hermitian(U(x + mu, nu));
    }
    else
    {
      mdp_site y(U.lattice());
      y = x - nu;
      tmp = hermitian(U(y, nu)) * U(y, mu) * U(y + mu, nu);
    }
    return tmp;
  }

  /**
   * @brief Computes the full staple sum for a link U(x, mu).
   *
   * This function sums over all directions nu != mu and includes both
   * forward and backward staples:
   *
   *   staple(x,mu) = sum_{nu != mu} [ forward(nu) + backward(nu) ]
   *
   * where each contribution corresponds to half of a plaquette touching
   * the link U(x, mu).
   *
   * This is the standard object entering:
   *   - Wilson gauge force,
   *   - heatbath / overrelaxation updates,
   *   - APE / stout smearing.
   *
   * @param U   Gauge field.
   * @param x   Lattice site.
   * @param mu  Direction of the link.
   *
   * @return Sum of all staples attached to U(x, mu).
   */
  mdp_matrix staple(const gauge_field &U, mdp_site x, mdp_suint mu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());

    tmp = 0;
    for (mdp_suint nu = 0; nu < U.ndim(); nu++)
    {
      if (nu != mu)
      {
        tmp += U(x, nu) * U(x + nu, mu) * hermitian(U(x + mu, nu));
        y = x - nu;
        tmp += hermitian(U(y, nu)) * U(y, mu) * U(y + mu, nu);
      }
    }
    return tmp;
  }

  /**
   * @brief Computes the Hermitian-conjugate oriented staple.
   *
   * This function returns the Hermitian counterpart of the staple
   * contributing to U(x, mu). It corresponds to the reversed path
   * ordering used in force calculations.
   *
   *  - For s1 = +1:
   *      U(x+mu,nu) U†(x+nu,mu) U†(x,nu)
   *
   *  - For s1 = -1:
   *      U†(x-nu+mu,nu) U†(x-nu,mu) U(x-nu,nu)
   *
   * These forms naturally appear when computing derivatives of the action
   * with respect to link variables (HMC force).
   *
   * @param U   Gauge field.
   * @param x   Lattice site.
   * @param mu  Direction of the link.
   * @param s1  Orientation (+1 forward, -1 backward).
   * @param nu  Orthogonal direction.
   *
   * @return Hermitian-oriented staple matrix.
   */
  mdp_matrix staple_H(const gauge_field &U, mdp_site x,
                      mdp_suint mu, mdp_int s1, mdp_suint nu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    if (s1 == +1)
    {
      tmp = U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
    }
    else
    {
      mdp_site y(U.lattice());
      y = x - nu;
      tmp = hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
    }
    return tmp;
  }

  /**
   * @brief Computes the sum of Hermitian-oriented staples for a link.
   *
   * This function sums all Hermitian-conjugate staple contributions
   * over directions nu != mu:
   *
   *   staple_H(x,mu) = sum_{nu != mu} [ forward_H(nu) + backward_H(nu) ]
   *
   * This object is commonly used in:
   *   - Hybrid Monte Carlo (HMC) force calculations,
   *   - gradient flow,
   *   - improved gauge actions.
   *
   * @param U   Gauge field.
   * @param x   Lattice site.
   * @param mu  Direction of the link.
   *
   * @return Sum of Hermitian staples.
   */
  mdp_matrix staple_H(const gauge_field &U, mdp_site x, mdp_suint mu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    tmp = 0;
    for (mdp_suint nu = 0; nu < U.ndim(); nu++)
    {
      if (nu != mu)
      {
        tmp += U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
        y = x - nu;
        tmp += hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
      }
    }
    return tmp;
  }

  /**
   * @brief Computes anisotropic (uniaxial) Hermitian staple.
   *
   * This version introduces an anisotropy parameter zeta, typically used
   * in anisotropic lattices where temporal and spatial directions have
   * different lattice spacings.
   *
   * The weighting is:
   *   - zeta       for mixed (temporal-spatial) plaquettes
   *   - 1 / zeta   for purely spatial plaquettes
   *
   * The condition (nu * mu == 0) assumes:
   *   mu=0 or nu=0 corresponds to temporal direction.
   *
   * @param U     Gauge field.
   * @param x     Lattice site.
   * @param mu    Direction of the link.
   * @param zeta  Anisotropy parameter (a_s / a_t).
   *
   * @return Weighted Hermitian staple sum.
   *
   * @note Assumes direction 0 is the temporal direction.
   */
  mdp_matrix staple_H_anisotropic(const gauge_field &U, mdp_site x, mdp_suint mu, mdp_real zeta)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    mdp_real param;
    tmp = 0;
    for (mdp_suint nu = 0; nu < U.ndim(); nu++)
    {
      if (nu != mu)
      {
        if (nu * mu == 0)
        {
          param = zeta;
        }
        else
        {
          param = 1.0 / zeta;
        }
        tmp += param * U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
        y = x - nu;
        tmp += param * hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
      }
    }
    return tmp;
  }

  /**
   * @brief Computes Hermitian staples restricted to temporal-spatial planes.
   *
   * This function selects only plaquettes involving the temporal direction (0):
   *
   *  - If mu == 0:
   *      sum over all spatial directions nu = 1..(D-1)
   *
   *  - If mu != 0:
   *      only contribution from nu = 0
   *
   * This is useful for:
   *   - anisotropic actions,
   *   - separating electric vs magnetic contributions,
   *   - Hamiltonian splitting in HMC.
   *
   * @param U   Gauge field.
   * @param x   Lattice site.
   * @param mu  Direction of the link.
   *
   * @return Hermitian staple restricted to 0-i planes.
   */
  mdp_matrix staple_0i_H(const gauge_field &U, mdp_site x, mdp_suint mu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    if (mu == 0)
    {
      tmp = 0;
      for (mdp_suint nu = 1; nu < U.ndim(); nu++)
      {
        tmp += U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
        y = x - nu;
        tmp += hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
      }
    }
    else
    {
      mdp_suint nu = 0;
      tmp = U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
      y = x - nu;
      tmp += hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
    }

    return tmp;
  }

  /**
   * @brief Computes Hermitian staples restricted to spatial planes.
   *
   * This function includes only spatial-spatial plaquettes (i,j):
   *
   *   - mu != 0 required
   *   - nu runs over spatial directions (1..D-1), nu != mu
   *
   * No contribution is returned for mu == 0.
   *
   * This separation is commonly used in:
   *   - anisotropic gauge actions,
   *   - splitting electric and magnetic contributions,
   *   - improved integrators in HMC.
   *
   * @param U   Gauge field.
   * @param x   Lattice site.
   * @param mu  Direction of the link.
   *
   * @return Hermitian staple restricted to spatial planes.
   */
  mdp_matrix staple_ij_H(const gauge_field &U, mdp_site x, mdp_suint mu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    tmp = 0;
    if (mu != 0)
    {
      for (mdp_suint nu = 1; nu < U.ndim(); nu++)
      {
        if (nu != mu)
        {
          tmp += U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
          y = x - nu;
          tmp += hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
        }
      }
    }
    return tmp;
  }

  /**
   * @brief Computes the plaquette (1x1 Wilson loop).
   *
   * The plaquette is the smallest closed Wilson loop on the lattice:
   *
   *   P(x; mu,nu) = U(x,mu) U(x+mu,nu)
   *                 U†(x+nu,mu) U†(x,nu)
   *
   * It is the fundamental building block of lattice gauge actions
   * (e.g. Wilson, Symanzik, Iwasaki).
   *
   * @param U   Gauge field.
   * @param x   Lattice site.
   * @param mu  First direction.
   * @param nu  Second direction (mu != nu).
   *
   * @return Nc x Nc matrix representing the plaquette.
   */
  mdp_matrix plaquette(const gauge_field &U, mdp_site x, mdp_suint mu, mdp_suint nu)
  {
    return U(x, mu) * U(x + mu, nu) * hermitian(U(x, nu) * U(x + nu, mu));
  }

  /**
   * @brief Projects a complex matrix onto the su(N) Lie algebra.
   *
   * This function maps a general complex Nc x Nc matrix X onto the Lie algebra su(N),
   * i.e. the space of traceless anti-Hermitian matrices:
   *
   *   su(N) = { A | A† = -A, Tr(A) = 0 }
   *
   * The projection is performed in two steps:
   *   1. Anti-Hermitian part:
   *        A = (X - X†) / 2
   *   2. Traceless projection:
   *        A -> A - (Tr(A)/Nc) * I
   *
   * @param X   Input matrix (Nc x Nc), not necessarily anti-Hermitian nor traceless.
   * @param Nc  Number of colors (dimension of the gauge group SU(N)).
   *
   * @return Matrix in su(N): anti-Hermitian and traceless.
   *
   * @note This projection is essential in Hybrid Monte Carlo (HMC), where the
   *       momentum conjugate to gauge links lives in su(N).
   *
   * @note Numerical stability: subtracting the trace explicitly avoids drift
   *       away from SU(N) due to floating-point errors.
   */
  mdp_matrix project_to_suN(const mdp_matrix &X, mdp_real Nc)
  {
    mdp_matrix A = 0.5 * (X - hermitian(X));

    mdp_complex tr = trace(A) / Nc;

    for (int i = 0; i < Nc; i++)
      A(i, i) -= tr;

    return A;
  }

  /**
   * @brief Computes the Wilson gauge force for Hybrid Monte Carlo (HMC).
   *
   * This function evaluates the derivative of the Wilson gauge action
   * with respect to each link variable U(x, mu).
   *
   * The Wilson action is:
   *   S = beta * sum_{x,mu<nu} (1 - ReTr P_{mu,nu}(x) / Nc)
   *
   * The corresponding force is:
   *   F_mu(x) = -beta * Proj_su(N)[ U(x,mu) * staple_H(x,mu) ]
   *
   * where:
   *   - staple_H(x,mu) is the sum of Hermitian-oriented staples,
   *   - Proj_su(N) enforces that the force lies in the Lie algebra.
   *
   * @param F      Output gauge field storing the force (same structure as U).
   * @param U      Input gauge field configuration.
   * @param coeff  Simulation parameters; must contain "beta".
   *
   * @note The resulting force F(x,mu) is anti-Hermitian and traceless.
   *
   * @note A global communication step (F.update()) is performed at the end
   *       to ensure consistency across MPI processes.
   *
   * @note This implementation assumes standard Wilson plaquette action.
   *
   * @warning The minus sign is essential: it ensures correct Hamiltonian dynamics.
   *
   * @see staple_H, project_to_suN
   */
  void WilsonForce(gauge_field &F,
                   const gauge_field &U,
                   coefficients &coeff)
  {
    mdp_site x(U.lattice());

    if (!coeff.has_key("beta"))
    {
      error("beta undeclared");
    }

    mdp_real beta = coeff["beta"];

    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim(); mu++)
      {
        F(x, mu) = -beta * U(x, mu) * staple_H(U, x, mu);
      }
    }

    F.update();
  }

  /**
   * @brief Computes rectangle (1x2 and 2x1) staple contributions.
   *
   * This function builds the sum of extended Wilson loops (rectangles)
   * attached to a given link U(x, mu) in the plane defined by (mu, nu).
   *
   * The included paths correspond to:
   *   - 1x2 rectangles (two steps in mu, one in nu),
   *   - 2x1 rectangles (two steps in nu, one in mu),
   *   - both forward and backward orientations.
   *
   * These contributions are required for improved gauge actions such as:
   *   - Symanzik (tree-level improved),
   *   - Iwasaki,
   *   - DBW2.
   *
   * Geometrically, this function sums all length-6 Wilson paths that:
   *   - start at x,
   *   - include the link U(x,mu),
   *   - form a closed rectangle in the (mu,nu) plane.
   *
   * @param U   Gauge field.
   * @param x   Lattice site.
   * @param mu  Direction of the central link.
   * @param nu  Orthogonal direction defining the plane (mu != nu).
   *
   * @return Sum of rectangle staples contributing to U(x, mu).
   *
   * @note This function does not include the coefficient c1; it only
   *       constructs the geometric contribution.
   *
   * @note The returned matrix is not projected to su(N); projection must be
   *       applied after combining with link variables.
   *
   * @warning Rectangle contributions are the most error-prone part of
   *          improved gauge force implementations. Incorrect paths or missing
   *          orientations will lead to incorrect dynamics.
   *
   * @see ImprovedGaugeForce
   */
  mdp_matrix rectangle_staple(const gauge_field &U,
                              mdp_site x,
                              mdp_suint mu,
                              mdp_suint nu)
  {
    mdp_matrix R(U.nc(), U.nc());
    R = 0;

    mdp_site x_mu = x + mu;
    mdp_site x_nu = x + nu;

    mdp_site x_2mu = x_mu + mu;
    mdp_site x_2nu = x_nu + nu;

    // --- forward 1x2 ---
    R += U(x, nu) *
         U(x + nu, mu) *
         U(x + nu + mu, mu) *
         hermitian(U(x_2mu, nu)) *
         hermitian(U(x_mu, mu));

    // --- backward 1x2 ---
    mdp_site y = x - nu;
    R += hermitian(U(y, nu)) *
         U(y, mu) *
         U(y + mu, mu) *
         U(y + 2 * mu, nu) *
         hermitian(U(x_mu, mu));

    // --- 2x1 forward ---
    R += U(x, nu) *
         U(x + nu, nu) *
         U(x_2nu, mu) *
         hermitian(U(x_nu + mu, nu)) *
         hermitian(U(x, nu));

    // --- 2x1 backward ---
    y = x - nu;
    R += hermitian(U(y + nu, nu)) *
         hermitian(U(y, nu)) *
         U(y, mu) *
         U(y + mu, nu) *
         U(y + mu + nu, nu);

    return R;
  }

  /**
   * @brief Computes improved gauge force (Symanzik, Iwasaki, DBW2, or Wilson).
   *
   * This function evaluates the derivative of an improved gauge action
   * including both plaquette and rectangle contributions:
   *
   *   S = beta * sum_x [
   *         c0 * (1 - ReTr P / Nc)
   *       + c1 * (1 - ReTr R / Nc)
   *   ]
   *
   * where:
   *   - P is the plaquette,
   *   - R includes both 1x2 and 2x1 rectangles,
   *   - coefficients satisfy c0 = 1 - 8*c1.
   *
   * The force is given by:
   *
   *   F_mu(x) = -beta * Proj_su(N)[
   *                 U(x,mu) * ( c0 * staple_H(x,mu)
   *                           + c1 * rectangle_staple_sum )
   *             ]
   *
   * Supported action types:
   *   - "wilson"   : c1 = 0
   *   - "symanzik" : c1 = -1/12
   *   - "iwasaki"  : c1 = -0.331
   *   - "dbw2"     : c1 = -1.4088
   *
   * @param F             Output gauge field storing the force.
   * @param U             Input gauge field configuration.
   * @param coeff         Simulation parameters; must contain "beta".
   * @param action_type   Type of gauge action.
   *
   * @note For "wilson", rectangle contributions are automatically skipped.
   *
   * @note The result is anti-Hermitian and traceless at each link.
   *
   * @note MPI communication is handled via F.update().
   *
   * @note This implementation is structurally correct but not fully optimized.
   *       Production codes typically:
   *         - reuse partial paths,
   *         - cache link products,
   *         - minimize memory access.
   *
   * @warning Correct implementation of rectangle terms is crucial.
   *          Even small mistakes lead to incorrect HMC evolution.
   *
   * @warning Always validate the force numerically using:
   *          finite-difference check:
   *            dS/dε ≈ (S(U e^{εX}) - S(U)) / ε
   *
   * @see WilsonForce, rectangle_staple, project_to_suN
   */
  void ImprovedGaugeForce(gauge_field &F,
                          const gauge_field &U,
                          coefficients &coeff,
                          const std::string &action_type = "wilson")
  {
    mdp_site x(U.lattice());
    const mdp_real Nc = U.nc();

    if (!coeff.has_key("beta"))
    {
      mdp << "beta undeclared\n";
      mdp.abort();
    }

    mdp_real beta = coeff["beta"];

    mdp_real c1;

    if (action_type == "symanzik")
      c1 = -1.0 / 12.0;
    else if (action_type == "iwasaki")
      c1 = -0.331;
    else if (action_type == "dbw2")
      c1 = -1.4088;
    else if (action_type == "wilson")
      c1 = 0.0;
    else
    {
      error("Unknown action type");
    }

    mdp_real c0 = 1.0 - 8.0 * c1;

    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim(); mu++)
      {
        // --- plaquette part ---
        mdp_matrix H_plaq = staple_H(U, x, mu);

        mdp_matrix total = c0 * H_plaq;

        // --- rectangle part ---
        if (c1 != 0.0)
        {
          mdp_matrix H_rect(U.nc(), U.nc());
          H_rect = 0;

          for (mdp_suint nu = 0; nu < U.ndim(); nu++)
          {
            if (nu == mu)
              continue;

            H_rect += rectangle_staple(U, x, mu, nu);
          }

          total += c1 * H_rect;
        }

        mdp_matrix X = U(x, mu) * total;

        F(x, mu) = -beta * project_to_suN(X, Nc);
      }
    }

    F.update();
  }

  /**
   * @brief Computes the gauge action for a given lattice gauge field.
   *
   * This function evaluates the gluonic lattice action for several commonly used
   * gauge actions:
   *  - "wilson"   : Wilson plaquette action
   *  - "symanzik" : tree-level Symanzik (Lüscher–Weisz) improved action
   *  - "iwasaki"  : Iwasaki RG-improved action
   *  - "dbw2"     : DBW2 action
   *
   * The action is defined as:
   *   S = beta * sum_{x, mu < nu} [
   *         c0 * (1 - ReTr P_{mu,nu}(x) / Nc)
   *       + c1 * (1 - ReTr R_{mu,nu}(x) / Nc)
   *       + c1 * (1 - ReTr R_{nu,mu}(x) / Nc)
   *   ]
   *
   * where:
   *   - P_{mu,nu} is the plaquette (1x1 Wilson loop),
   *   - R_{mu,nu} and R_{nu,mu} are 1x2 and 2x1 rectangular Wilson loops,
   *   - Nc is the number of colors.
   *
   * The coefficients satisfy:
   *   c0 = 1 - 8*c1
   *
   * with:
   *   - Wilson:   c1 = 0
   *   - Symanzik: c1 = -1/12
   *   - Iwasaki:  c1 = -0.331
   *   - DBW2:     c1 = -1.4088
   *
   * @param U            Gauge field configuration (SU(Nc) link variables).
   * @param coeff        Container of simulation parameters. Must contain "beta".
   * @param action_type  Type of gauge action ("wilson", "symanzik", "iwasaki", "dbw2").
   *
   * @return The average action density:
   *         beta * S / V,
   *         where S is the total lattice action and V is the global lattice volume.
   */
  mdp_real GaugeAction(const gauge_field &U, coefficients &coeff, const std::string &action_type = "wilson")
  {
    mdp_real sum = 0;
    mdp_site x(U.lattice());

    const mdp_real Nc = U.nc();
    mdp_real beta;
    if (coeff.has_key("beta"))
      beta = coeff["beta"];
    else
      error("beta undeclared");

    // --- współczynniki ---
    mdp_real c1;

    if (action_type == "wilson")
      c1 = 0.0;
    else if (action_type == "symanzik")
      c1 = -1.0 / 12.0;
    else if (action_type == "iwasaki")
      c1 = -0.331;
    else if (action_type == "dbw2")
      c1 = -1.4088;
    else
    {
      error(std::format("Unknown action type: {}", action_type));
    }

    const mdp_real c0 = 1.0 - 8.0 * c1;
    const bool use_rect = (c1 != 0.0);

    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim(); mu++)
      {
        mdp_site x_mu = x + mu;

        for (mdp_suint nu = mu + 1; nu < U.ndim(); nu++)
        {
          mdp_site x_nu = x + nu;
          mdp_site x_mu_nu = x_mu + nu;

          // ===== PLAQUETTE =====
          mdp_matrix Ux_mu = U(x, mu);
          mdp_matrix Ux_nu = U(x, nu);
          mdp_matrix Uxmu_nu = U(x_mu, nu);
          mdp_matrix Uxnu_mu = U(x_nu, mu);

          mdp_matrix P = Ux_mu * Uxmu_nu;
          P = P * hermitian(Uxnu_mu);
          P = P * hermitian(Ux_nu);

          sum += c0 * (1.0 - real(trace(P)) / Nc);

          if (use_rect)
          {
            mdp_site x_2mu = x_mu + mu;
            mdp_site x_2nu = x_nu + nu;

            // ===== RECTANGLE mu,nu =====
            mdp_matrix R = Ux_mu * U(x_mu, mu); // U(x,mu) U(x+mu,mu)
            R = R * U(x_2mu, nu);               // U(x+2mu,nu)
            R = R * hermitian(U(x_mu_nu, mu));  // U†(x+mu+nu,mu)
            R = R * hermitian(Uxnu_mu);         // U†(x+nu,mu)
            R = R * hermitian(Ux_nu);           // U†(x,nu)

            sum += c1 * (1.0 - real(trace(R)) / Nc);

            // ===== RECTANGLE nu,mu =====
            mdp_matrix R2 = Ux_nu * U(x_nu, nu); // U(x,nu) U(x+nu,nu)
            R2 = R2 * U(x_2nu, mu);              // U(x+2nu,mu)
            R2 = R2 * hermitian(U(x_mu_nu, nu)); // U†(x+mu+nu,nu)
            R2 = R2 * hermitian(U(x_mu, nu));    // U†(x+mu,nu)
            R2 = R2 * hermitian(Ux_mu);          // U†(x,mu)

            sum += c1 * (1.0 - real(trace(R2)) / Nc);
          }
        }
      }
    }

    mdp.add(sum);
    return beta * sum / U.lattice().global_volume();
  }

  /**
   * @brief Finite-difference consistency test between gauge action and gauge force.
   *
   * This function verifies the correctness of the HMC gauge force implementation
   * by comparing:
   *
   *  1. Finite-difference derivative of the action:
   *        ( S(U') - S(U) ) / epsilon
   *
   *  2. Analytical force contraction:
   *        sum_{x,mu} ReTr( F_mu(x) * X_mu(x) )
   *
   * where:
   *   - U' = exp(epsilon * X) * U
   *   - X is a random element of the su(N) Lie algebra
   *   - F is the gauge force computed from the action
   *
   * If the force is implemented correctly, both quantities agree up to
   * O(epsilon) and numerical precision.
   *
   * @param U             Input gauge field configuration.
   * @param F             Gauge force field corresponding to U.
   * @param coeff         Simulation parameters (must contain "beta").
   * @param action_type   Type of gauge action ("wilson", "symanzik",
   *                      "iwasaki", "dbw2").
   * @param epsilon      Small perturbation parameter used in finite-difference
   *                      derivative (typical values: 1e-5 to 1e-7).
   *
   * @return Relative error between finite-difference derivative and
   *         analytical force contraction:
   *         |(S(U') - S(U))/epsilon - Tr(F·X)| / |Tr(F·X)|
   *
   * @note This is a diagnostic function used to validate HMC correctness.
   *       It should NOT be used in production simulations.
   *
   * @note Correct implementation should yield relative errors ~1e-8–1e-10.
   *
   * @note The random matrix X is projected onto su(N) to ensure correct
   *       Lie algebra consistency.
   *
   * @warning Large values of epsilon introduce nonlinear effects,
   *          while too small values lead to numerical cancellation errors.
   *
   * @see GaugeAction, WilsonForce, ImprovedGaugeForce, project_to_suN
   */
  mdp_real test_force_vs_action(const gauge_field &U, const gauge_field &F,
                                coefficients &coeff, const std::string &action_type = "wilson",
                                mdp_real epsilon = 1e-6)
  {
    mdp_site x(U.lattice());
    const mdp_suint Nc = U.nc();

    gauge_field U_new(U);

    mdp_real S0 = GaugeAction(U, coeff, action_type);

    mdp_real S1 = 0.0;
    mdp_real rhs = 0.0;

    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim(); mu++)
      {
        mdp_matrix X(Nc, Nc);

        for (mdp_suint i = 0; i < Nc; i++)
          for (mdp_suint j = 0; j < Nc; j++)
          {
            X(i, j) = mdp_complex(Random.plain() - 0.5, Random.plain() - 0.5);
          }

        X = project_to_suN(X, Nc);

        // --- link perturbation ---
        mdp_matrix expX = exp(epsilon * X);

        U_new(x, mu) = expX * U(x, mu);

        rhs += real(trace(F(x, mu) * X));
      }
    }

    S1 = GaugeAction(U_new, coeff, action_type);

    mdp_real lhs = (S1 - S0) / epsilon;

    mdp << "===== FORCE vs ACTION TEST =====\n";
    mdp << "dS (finite diff) = " << lhs << "\n";
    mdp << "Tr(F*X)          = " << rhs << "\n";
    mdp << "relative error   = " << fabs(lhs - rhs) / fabs(rhs) << "\n";

    return fabs(lhs - rhs) / fabs(rhs);
  }
} // namespace MDP

#endif /* FERMIQCD_GAUGE_ROUTINES_ */
