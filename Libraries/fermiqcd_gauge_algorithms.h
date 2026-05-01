/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_algorithms.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various stuff for gauge field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_GAUGE_ALGORITHMS_
#define FERMIQCD_GAUGE_ALGORITHMS_

#include <vector>
#include <initializer_list>
#include <cstddef>
#include "fermiqcd_gauge_field.h"
#include "fermiqcd_gauge_routines.h"

namespace MDP
{
  /** @brief make a cold gauge configuration
   */
  void set_cold(gauge_field &U)
  {
    begin_function("set_cold");
    mdp << "Creating a cold gauge configuration\n";
    mdp_site x(U.lattice());

    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim(); mu++)
        U(x, mu) = mdp_identity(U.nc());
    }
    U.update();
    end_function("set_cold");
  }

  /** @brief Make a hot gauge configuration
   */
  void set_hot(gauge_field &U)
  {
    begin_function("set_hot");
    mdp << "Creating a hot gauge configuration\n";
    mdp_site x(U.lattice());

    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim(); mu++)
        U(x, mu) = U.lattice().random(x).SU(U.nc());
    }
    U.update();
    end_function("set_hot");
  }

  /// Check that gauge field is unitary within precision
  void check_unitarity(const gauge_field &U, mdp_real precision = mdp_precision)
  {
    begin_function("check_unitarity");
    mdp_site x(U.lattice());
    mdp_int how_many = 0;

    forallsitesandcopies(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim(); mu++)
        if (max(inv(U(x, mu)) - hermitian(U(x, mu))) > precision)
          how_many++;
    }
    mdp.add(how_many);
    mdp << "Non unitary links found=" << how_many << "\n";
    end_function("check_unitarity");
  }

  /// Compute average plaquette on plane mu-nu
  mdp_real average_plaquette(const gauge_field &U, mdp_suint mu, mdp_suint nu)
  {
    mdp_real sum = 0;
    mdp_site x(U.lattice());

    forallsites(x)
    {
      sum += real(trace(plaquette(U, x, mu, nu)));
    }
    mdp.add(sum);
    return sum / (U.lattice().global_volume() * U.nc());
  }

  /// Compute average plaquette (all planes)
  mdp_real average_plaquette(const gauge_field &U)
  {
    mdp_real sum = 0;
    mdp_site x(U.lattice());

    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim() - 1; mu++)
        for (mdp_suint nu = mu + 1; nu < U.ndim(); nu++)
          sum += real(trace(plaquette(U, x, mu, nu)));
    }
    mdp.add(sum);
    return 2.0 * sum / (U.ndim() * (U.ndim() - 1) * U.lattice().global_volume() * U.nc());
  }

  /** @brief Compute average Time plaquette
   */
  mdp_real TimePlaquette(const gauge_field &U)
  {
    mdp_real sum = 0;
    mdp_site x(U.lattice());
    mdp_suint mu = 0;

    forallsites(x)
    {
      for (mdp_suint nu = mu + 1; nu < U.ndim(); nu++)
        sum += real(trace(plaquette(U, x, mu, nu)));
    }
    mdp.add(sum);
    return sum / (U.lattice().global_volume() * (U.ndim() - 1));
  }

  /** @brief Compute average Space plaquette
   */
  mdp_real SpacePlaquette(const gauge_field &U)
  {
    mdp_real sum = 0;
    mdp_site x(U.lattice());

    forallsites(x)
    {
      for (mdp_suint mu = 1; mu < U.ndim() - 1; mu++)
        for (mdp_suint nu = mu + 1; nu < U.ndim(); nu++)
          sum += real(trace(plaquette(U, x, mu, nu)));
    }
    mdp.add(sum);
    return 2 * sum / (U.lattice().global_volume() * (U.ndim() - 2) * (U.ndim() - 1));
  }

  /**
   * @brief Computes the local Polyakov loop matrix L(x) at spatial sites.
   *
   * The Polyakov loop is defined as an ordered product of temporal gauge links:
   * L(x) = \prod_{t=0}^{Nt-1} U_0(t, x)
   *
   * Only sites with x(0) == 0 (first time slice) store the result.
   *
   * @param U Gauge field
   *
   * @return Field of Polyakov loop matrices
   */
  mdp_matrix_field PolyakovField(const gauge_field &U)
  {
    mdp_site x(U.lattice());
    mdp_matrix_field L(U.lattice(), U.nc(), U.nc());

    const mdp_uint Nt = U.lattice().size(0);

    forallsites(x)
    {
      if (x(0) == 0)
      {
        mdp_site y = x;
        L(x) = 1;

        for (mdp_uint t = 0; t < Nt; t++)
        {
          L(x) *= U(y, 0);
          y = y + 0;
        }
      }
    }

    L.update();
    return L;
  }

  /**
   * @brief Computes the spatially averaged Polyakov loop.
   *
   * The observable is:
   *   <L> = (1 / (Ns^d * Nc)) * sum_x Tr L(x)
   *
   * where the sum runs over spatial sites only (t = 0 slice).
   *
   * @param U Gauge field
   *
   * @return Averaged Polyakov loop
   */
  mdp_complex PolyakovLoop(const gauge_field &U)
  {
    mdp_site x(U.lattice());
    mdp_complex sum = 0;

    mdp_matrix_field L = PolyakovField(U);
    mdp_uint spatial_volume = U.lattice().global_volume() / U.lattice().size(0);

    forallsites(x)
    {
      if (x(0) == 0)
        sum += trace(L(x));
    }

    mdp.add(sum);

    return sum / (spatial_volume * U.nc());
  }

  /** @brief Relaxation algorithm
   */
  void relaxation(gauge_field &U, mdp_uint n_iter = 1)
  {
    if (U.nc() == 1)
      error("fermiqcd_gauge_algorithms/heatbath(): U(1)? (use metropolis)");

    mdp_matrix M;
    mdp_complex a[4];
    mdp_site x(U.lattice());
    mdp_real dk;
    mdp_real e[4];

    for (mdp_uint iter = 0; iter < n_iter; iter++)
      for (mdp_parity parity : {EVEN, ODD})
        for (mdp_suint mu = 0; mu < U.ndim(); mu++)
        {
          forallsitesofparity(x, parity)
          {
            for (mdp_suint i = 0; i < U.nc() - 1; i++)
              for (mdp_suint j = i + 1; j < U.nc(); j++)
              {
                M = U(x, mu) * staple_H(U, x, mu);

                e[0] = real(M(i, i) + M(j, j));
                e[1] = imag(M(i, j) + M(j, i));
                e[2] = real(M(i, j) - M(j, i));
                e[3] = imag(M(i, i) - M(j, j));

                dk = 0;
                for (mdp_suint t = 0; t < 4; t++)
                {
                  dk += std::pow(e[t], 2);
                }
                dk = std::sqrt(dk);
                for (mdp_suint t = 0; t < 4; t++)
                {
                  e[t] /= dk;
                }

                a[0] = mdp_complex(e[0], -e[3]);
                a[1] = -mdp_complex(e[2], e[1]);
                a[2] = mdp_complex(e[2], -e[1]);
                a[3] = mdp_complex(e[0], e[3]);

                for (mdp_suint k = 0; k < U.nc(); k++)
                {
                  mdp_complex tmpUik = (a[0] * a[0] + a[2] * a[1]) * U(x, mu, i, k) + (a[1] * (a[0] + a[3])) * U(x, mu, j, k);
                  U(x, mu, j, k) = (a[2] * (a[0] + a[3])) * U(x, mu, i, k) + (a[1] * a[2] + a[3] * a[3]) * U(x, mu, j, k);
                  U(x, mu, i, k) = tmpUik;
                }
              }
          }
          U.update(parity, mu, U.nc() * U.nc());
        }
  }

  /// Unitarize SU(N) matrix
  /*
   *          [ e0 e1 e2 ]
   * SU(3) =  [ |  |  |  ]
   *          [ .  .  .  ]
   *
   *
   */
  mdp_real unitarize(gauge_field &U)
  {

    mdp_real stats;
    mdp_site x(U.lattice());
    mdp_real precision = 0;
    mdp_matrix e(U.nc(), U.nc());
    mdp_matrix m(U.nc() - 1, U.nc() - 1);
    mdp_real quadnorm;
    mdp_complex scalar;
    e = U(x, 0);

    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim(); mu++)
      {
        precision += abs(1.0 - det(U(x, mu)));

        // Gram-Schmidt's orthonormalization method for the 'U.nc' vectors (e1, e2, ..., en)
        for (mdp_suint i = 0; i < U.nc(); i++)
        {
          for (mdp_suint k = 0; k < i; k++)
          {
            scalar = 0;
            for (mdp_suint j = 0; j < U.nc(); j++)
            {
              scalar += conj(e(k, j)) * e(i, j);
            } //<ek|ei>
            for (mdp_suint j = 0; j < U.nc(); j++)
            {
              e(i, j) -= e(k, j) * scalar;
            } // ei = ei - ek<ek|ei>
          }
          quadnorm = 0;
          for (mdp_suint j = 0; j < U.nc(); j++)
            quadnorm += abs2(e(i, j)); // quadnorm = ||ei||^2

          for (mdp_suint j = 0; j < U.nc(); j++)
            e(i, j) /= std::sqrt(quadnorm); // ei = ei/||ei||
        }

        // The last vector is calculated as the cross product: ek = ei x ej
        mdp_suint i = U.nc() - 1;
        for (mdp_suint j = 0; j < U.nc(); j++)
        {
          for (mdp_suint s = 0; s < U.nc() - 1; s++)
          {
            for (mdp_suint k = 0; k < U.nc() - 1; k++)
            {
              if (k < j)
                m(s, k) = e(s, k);
              else
                m(s, k) = e(s, k + 1);
            }
          }
          e(i, j) = det(conj(m));
          if ((j + i) % 2 == 1)
            e(i, j) *= -1;
        }
        U(x, mu) = e;
      }
    }
    mdp.add(precision);
    stats = precision / (U.ndim() * U.lattice().global_volume());

    return stats;
  }

  /**
   * @brief Computes Polyakov loop correlation function C_mu(r)
   *        for a given spatial direction mu.
   *
   * Definition:
   *   C_mu(r) = < Tr L(x) Tr L^\dagger(x + r * e_mu) >
   *
   * Averaged over all spatial sites (t = 0 slice).
   *
   * @param U Gauge field
   * @param mu Spatial direction (1 <= mu < ndim)
   *
   * @return Correlation in mu direction as function of distance r
   */
  std::vector<mdp_real> PolyCor(const gauge_field &U, mdp_suint mu)
  {
    mdp_site x(U.lattice());
    mdp_matrix_field L = PolyakovField(U);

    const mdp_suint ndim = U.ndim();

    if (mu == 0 || mu >= ndim)
      error("PolyCor: mu must be a spatial direction (1 <= mu < ndim)");

    const mdp_uint Ns = U.lattice().size(mu);

    std::vector<mdp_complex> sum(Ns, 0);
    mdp_uint count = 0;

    forallsites(x)
    {
      if (x(0) == 0)
      {
        count++;

        mdp_vector coords;
        for (mdp_suint nu = 0; nu < ndim; nu++)
          coords[nu] = x(nu);

        mdp_site y = x;
        for (mdp_uint r = 0; r < Ns; r++)
        {
          coords[mu] = (x(mu) + r) % Ns;
          y.set(coords);

          sum[r] += trace(L(x)) * conj(trace(L(y)));
        }
      }
    }

    for (mdp_uint r = 0; r < Ns; r++)
      mdp.add(sum[r]);

    std::vector<mdp_real> correlation(Ns);
    for (mdp_uint r = 0; r < Ns; r++)
      correlation[r] = sum[r].real() / count;

    return correlation;
  }

  /**
   * @brief Computes Polyakov loop correlation functions for all spatial directions.
   *
   * For each spatial direction mu = 1,...,ndim-1:
   *   C_mu(r) = < Tr L(x) Tr L^\dagger(x + r * e_mu) >
   *
   * @param U Gauge field
   * @return Matrix of size (ndim-1) x Ns,
   *         where each row corresponds to a spatial direction.
   */
  mdp_matrix PolyCor(const gauge_field &U)
  {
    const mdp_suint ndim = U.ndim();
    const mdp_uint Ns = U.lattice().size(1);

    mdp_matrix correlation(ndim - 1, Ns);

    for (mdp_suint mu = 1; mu < ndim; mu++)
    {
      std::vector<mdp_real> corr_mu = PolyCor(U, mu);

      for (mdp_uint r = 0; r < Ns; r++)
      {
        correlation(mu - 1, r) = corr_mu[r];
      }
    }

    return correlation;
  }

  /// Given a field U compute the chromo-electro-magnetic field U.em
  void compute_em_field(gauge_field &U)
  {
    mdp_suint nc = U.nc();
    mdp_site x(U.lattice());
    // It is fine to use Nmdp_matrix even if there is twist .. how lucky!
    U.em.deallocate_field();
    U.em.allocate_em_field(U.lattice(), U.nc());
    mdp_matrix A(nc, nc);
    mdp_matrix b1(nc, nc), b2(nc, nc);
    /*
       Fast version of code for the clover term.
       A are the four clover leafs
    */
    forallsites(x)
    {
      for (mdp_suint mu = 0; mu < U.ndim() - 1; mu++)
        for (mdp_suint nu = mu + 1; nu < U.ndim(); nu++)
        {

          A =
              U(x, mu) * U(x + mu, nu) * hermitian(U(x, nu) * U(x + nu, mu)) +
              hermitian(U(x - nu, nu)) * U(x - nu, mu) *
                  U((x - nu) + mu, nu) * hermitian(U(x, mu)) +
              hermitian(U((x - mu) - nu, nu) * U(x - mu, mu)) * U((x - mu) - nu, mu) * U(x - nu, nu) +
              U(x, nu) * hermitian(U(x - mu, nu) * U((x + nu) - mu, mu)) * U(x - mu, mu);

          U.em(x, mu, nu) = (0.125) * (A - hermitian(A));
        }
    }
    U.em.update();
  }

  // /////////////////////////////////////////////////////
  // function to compute longlinks of V and attach them to
  // the handle1 of U
  // /////////////////////////////////////////////////////

  /// For use with asqtad staggered action
  /// Given field V makes a field U.long_links where (if length==2)
  /// @verbatim
  /// U.long_links(x,mu)=V(x,mu)*V(x+mu,mu);
  /// @endverbatim
  /// or (if length==3)
  /// @verbatim
  /// U.long_links(x,mu)=V(x,mu)*V(x+mu,mu)*V((x+mu)+mu,mu);
  /// @endverbatim
  void compute_long_links(gauge_field &U, gauge_field &V, mdp_sint length = 2)
  {
    if ((&(U.lattice()) != &(V.lattice())) || (U.nc() != V.nc()) || (U.ndim() != V.ndim()))
      error("fermiqcd_gauge_auxiliary_functions/compute_long_links: U and V are not compatible lattices");
    if (V.lattice().boundary_thickness() < length)
      error("fermiqcd_gauge_auxiliary_functions/compute_long_links: boundary thickness is not big enough");

    U.long_links.deallocate_field();
    U.long_links.allocate_mdp_nmatrix_field(V.lattice(), U.ndim(), U.nc(), U.nc());
    mdp_site x(U.lattice());

    if (length == 2)
      forallsites(x)
      {
        for (mdp_suint mu = 0; mu < V.ndim(); mu++)
          U.long_links(x, mu) = V(x, mu) * V(x + mu, mu);
      }
    if (length == 3)
      forallsites(x)
      {
        for (mdp_suint mu = 0; mu < V.ndim(); mu++)
          U.long_links(x, mu) = V(x, mu) * V(x + mu, mu) * V((x + mu) + mu, mu);
      }
    U.long_links.update();
  }

  // ////////////////////////////////////////////////////////////
  // ////////////////////////////////////////////////////////////
  // set phases for antiperiodic boundary conditions
  // ////////////////////////////////////////////////////////////
  // ////////////////////////////////////////////////////////////

  /// To set antiperiodic boundary conditions on in direction mu
  /// @verbatim
  ///    gauge_field U(lattice,nc);
  ///    // do heatbath on U
  ///    set_antiperiodic_phases(U,mu,true);
  ///    // use quarks (will have antiperiodic boundary conditions)
  ///    set_antiperiodic_phases(U,mu,false);
  /// @endverbatim
  void set_antiperiodic_phases(gauge_field &U, mdp_suint mu = 0, bool check = true)
  {
    begin_function("set_antiperiodic_phases");
    mdp_site x(U.lattice());

    if (check)
      mdp << "Setting antiperiodic boundary conditions on mu=" << mu << "\n";
    else
      mdp << "Removing antiperiodic boundary conditions on mu=" << mu << "\n";
    forallsites(x)
    {
      if (x(mu) == U.lattice().size(mu) - 1)
      {
        for (mdp_suint i = 0; i < U.nc(); i++)
          for (mdp_suint j = 0; j < U.nc(); j++)
            U(x, mu, i, j) *= -1;
      }
    }
    end_function("set_antiperiodic_phases");
  }

  /// takes a matrix M, performs a Cabibbo-Marinari cooling
  /// and returns the projected matrix
  mdp_matrix project_SU(mdp_matrix M, mdp_uint nstep = 1)
  {
    mdp_suint nc = M.rows();
    mdp_real e0, e1, e2, e3, dk, d;
    mdp_complex dc, u0, u1, u2, u3;
    mdp_matrix B(nc, nc);
    mdp_matrix C(nc, nc);
    mdp_matrix S(nc, nc);

    C = M;
    // /////////////////
    // preconditioning
    // /////////////////
    for (mdp_suint i = 0; i < nc; i++)
    {
      for (mdp_suint j = 0; j < i; j++)
      {
        dc = 0;
        for (mdp_suint k = 0; k < nc; k++)
          dc += conj(C(k, j)) * C(k, i);
        for (mdp_suint k = 0; k < nc; k++)
          C(k, i) -= dc * C(k, j);
      }
      d = 0;
      for (mdp_suint k = 0; k < nc; k++)
        d += std::pow(abs(C(k, i)), 2.0);
      d = std::sqrt(d);
      for (mdp_suint k = 0; k < nc; k++)
        C(k, i) /= d;
    }
    // ////////////////////////////
    // Cabibbo Marinari Projection
    // ////////////////////////////
    for (mdp_suint i = 0; i < nc; i++)
      for (mdp_suint j = 0; j < nc; j++)
        for (mdp_suint k = 0; k < nc; k++)
          B(i, j) += conj(M(k, i)) * C(k, j);
    for (mdp_uint step = 0; step < nstep; step++)
    {
      for (mdp_suint i = 0; i < nc - 1; i++)
        for (mdp_suint j = i + 1; j < nc; j++)
        {
          e0 = real(B(i, i)) + real(B(j, j));
          e1 = imag(B(i, j)) + imag(B(j, i));
          e2 = real(B(i, j)) - real(B(j, i));
          e3 = imag(B(i, i)) - imag(B(j, j));
          dk = std::sqrt(e0 * e0 + e1 * e1 + e2 * e2 + e3 * e3);
          u0 = mdp_complex(e0, -e3) / dk;
          u1 = mdp_complex(e2, -e1) / dk;
          u2 = mdp_complex(-e2, -e1) / dk;
          u3 = mdp_complex(e0, e3) / dk;
          // S=C;
          for (mdp_suint k = 0; k < nc; k++)
          {
            S(k, i) = C(k, i) * u0 + C(k, j) * u1;
            S(k, j) = C(k, i) * u2 + C(k, j) * u3;
          }
          if ((i == nc - 2) && (j == nc - 1))
            for (mdp_suint k = 0; k < nc; k++)
              for (mdp_suint l = 0; l < nc - 2; l++)
                S(k, l) = C(k, l);
          if ((i != nc - 2) || (j != nc - 1) || (step != nstep - 1))
            for (mdp_suint k = 0; k < nc; k++)
            {
              C(k, i) = B(k, i) * u0 + B(k, j) * u1;
              C(k, j) = B(k, i) * u2 + B(k, j) * u3;
              B(k, i) = C(k, i);
              B(k, j) = C(k, j);
              C(k, i) = S(k, i);
              C(k, j) = S(k, j);
            }
        }
    }
    return S;
  }

  class Path
  {
  public:
    // Each step is represented as a pair (orientation, mu).
    // Orientation +1 or -1 represents forward or backward step in mu direction.
    struct Step
    {
      int orientation; // either +1 or -1
      int direction;

      constexpr Step(int o = +1, int d = 0)
          : orientation(o >= 0 ? +1 : -1),
            direction(d)
      {
      }

      // [] accessors
      constexpr int &operator[](std::size_t i)
      {
        return (i == 0) ? orientation : direction;
      }

      constexpr const int &operator[](std::size_t i) const
      {
        return (i == 0) ? orientation : direction;
      }
    };

  private:
    std::vector<Step> m_steps;

  public:
    Path() = default;

    explicit Path(std::size_t length)
        : m_steps(length)
    {
    }

    constexpr Path(std::initializer_list<Step> init_steps)
        : m_steps(init_steps)
    {
    }

    constexpr std::size_t length() const
    {
      return m_steps.size();
    }

    constexpr Step &operator[](std::size_t i)
    {
      return m_steps[i];
    }

    constexpr const Step &operator[](std::size_t i) const
    {
      return m_steps[i];
    }

    constexpr auto begin() { return m_steps.begin(); }
    constexpr auto end() { return m_steps.end(); }

    constexpr auto begin() const { return m_steps.begin(); }
    constexpr auto end() const { return m_steps.end(); }

    constexpr void set_length(std::size_t n)
    {
      m_steps.resize(n);
    }

    // flip orientation
    static constexpr int flip(int o)
    {
      return -o;
    }

    constexpr void invert_path(int mu)
    {
      for (auto &step : m_steps)
      {
        if (step.direction == mu)
          step.orientation = flip(step.orientation);
      }
    }

    constexpr void rotate_path(int angle, int mu, int nu)
    {
      angle = (angle + 360) % 360;

      for (auto &step : m_steps)
      {
        if (step.direction != mu && step.direction != nu)
          continue;

        switch (angle)
        {
        case 90:
          if (step.direction == mu)
            step.direction = nu;
          else
          {
            step.direction = mu;
            step.orientation = flip(step.orientation);
          }
          break;

        case 180:
          step.orientation = flip(step.orientation);
          break;

        case 270:
          if (step.direction == mu)
          {
            step.direction = nu;
            step.orientation = flip(step.orientation);
          }
          else
            step.direction = mu;
          break;
        }
      }
    }
  };

  /// Takes a field U and path d of length and compute the average of
  /// the path on the entire lattice. Assumes computation can be done
  /// locally for each mdp_site
  ///
  /// Example:
  /// @verbatim
  ///   mdp_suint mu=0, nu=1;
  ///   gauge_field U(lattice,nc);
  ///   Path d = {{+1,mu},{+1,nu},{-1,mu},{-1,nu}};
  ///   mdp << "plaquette=" << average_path(U,d) << "\n";
  /// @endverbatim
  mdp_complex average_path(const gauge_field &U, const Path &d)
  {
    mdp_matrix_field psi1(U.lattice(), U.nc(), U.nc());
    mdp_matrix_field psi2(U.lattice(), U.nc(), U.nc());
    mdp_site x(U.lattice());
    mdp_complex sum = 0;

    for (size_t i = 0; i < d.length(); i++)
    {
      const auto &[orientation, direction] = d[i];

      if (i == 0)
        forallsites(x)
            psi1(x) = U(x, orientation, direction);
      else
        forallsites(x)
            psi1(x) = psi1(x) * U(x, orientation, direction);

      if (i < d.length() - 1)
      {
        psi1.update();
        // signs are correct this way, thanks J.Flynn
        if (orientation == +1)
          forallsites(x) psi2(x) = psi1(x - direction);
        else
          forallsites(x) psi2(x) = psi1(x + direction);
        psi1 = psi2;
      }
    }

    forallsites(x)
        sum += trace(psi1(x));

    return sum / (1.0 * U.lattice().global_volume() * U.nc());
  }

  /// Builds the ordered product of gauge links along the path d starting at x.
  /// The result is computed locally for each site.
  ///
  /// Example:
  /// @verbatim
  ///   mdp_suint mu=0, nu=1;
  ///   gauge_field U(lattice,nc);
  ///   Path d = {{+1,mu},{+1,nu},{-1,mu},{-1,nu}};
  ///   forallsites(x) {
  ///     mdp_matrix plaquette = build_path(U, x, d);
  ///     std::cout << "plaquette(x)=" << trace(plaquette) << std::endl;
  ///   }
  /// @endverbatim
  mdp_matrix build_path(const gauge_field &U, mdp_site x, const Path &d)
  {
    mdp_matrix tmp = mdp_identity(U.nc());

    for (const auto &[orientation, direction] : d)
    {
      tmp = tmp * U(x, orientation, direction);
      x = (orientation < 0) ? x - direction : x + direction;
    }

    return tmp;
  }
} // namespace MDP

#endif /* FERMIQCD_GAUGE_ALGORITHMS_ */
