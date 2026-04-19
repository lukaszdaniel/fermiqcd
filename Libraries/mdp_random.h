/////////////////////////////////////////////////////////////////
/// @file mdp_random.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Class mdp_random (the random number generator of MDP)
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_RANDOM_
#define MDP_RANDOM_

#include <cmath>
#include "mdp_global_vars.h"
#include "mdp_matrix.h"
#include "mdp_communicator.h"
#ifdef MDP_USE_SFMT
#include "mdp_prng_sfmt.h"
#else
#include "mdp_prng.h"
#endif

namespace MDP
{
  /// @brief Class to generate random numbers
  ///
  /// You should not instantiate this class because:
  /// - there is a global object mdp_global_random
  /// - each field "lattice" has a parallel generator "lattice.random(x)"
  /// Example:
  /// @verbatim
  ///    // print a uniform number in (0,1)
  ///    std::cout << mdp_global_random.plain() << std::endl;
  ///    // print a gaussian number
  ///    std::cout << mdp_global_random.gaussian() << std::endl;
  ///    // print a random SU(10) matrix
  ///    std::cout << mdp_global_random.SU(10) << std::endl;
  /// @endverbatim
  class mdp_random
  {
  private:
#ifdef MDP_USE_SFMT
    mdp_prng_sfmt m_prng;
#else
    mdp_prng m_prng;
#endif
    bool m_has_gauss = false;
    mdp_real m_gauss_cache = 0;

  public:
    mdp_random(mdp_int k = 0) : m_prng()
    {
      if (k == 0)
        m_prng.initialize(ME);
    }

    void initialize(mdp_int seed)
    {
      m_prng.initialize(seed);
    }

    /// return a uniform random number in (0,1)
    inline mdp_real plain() noexcept
    {
      return m_prng.plain();
    }

    inline mdp_real uniform() noexcept
    {
      return m_prng.plain();
    }

    /// returns a gaussian random number
    mdp_real gaussian(mdp_real sigma = 1.0) noexcept
    {
      if (m_has_gauss)
      {
        m_has_gauss = false;
        return m_gauss_cache * sigma;
      }

      mdp_real u1 = m_prng.plain();
      mdp_real u2 = m_prng.plain();

      mdp_real r = std::sqrt(-2.0f * std::log(u1));
      mdp_real theta = 2.0f * Pi * u2;

      m_gauss_cache = r * std::cos(theta);
      m_has_gauss = true;

      return sigma * r * std::sin(theta);
    }

    /// draws a random mdp_real in (0,1) from a distribution using accept-reject
    mdp_real distribution(mdp_real (*fp)(mdp_real, void *), void *a = nullptr)
    {
      mdp_real x, y;
      do
      {
        x = m_prng.plain();
        y = m_prng.plain();
      } while (fp(x, a) >= y);
      return x;
    }

    /// returns a random SU(n) matrix using Cabibbo-Marinari
    mdp_matrix SU(mdp_suint n) noexcept
    {
      if (n == 1)
      {
        mdp_matrix tmp(1, 1);
        mdp_real alpha = 2.0 * Pi * m_prng.plain();
        tmp(0, 0) = mdp_complex(std::cos(alpha), std::sin(alpha));
        return tmp;
      }

      mdp_matrix M = mdp_identity(n);

      for (mdp_suint i = 0; i < n - 1; ++i)
      {
        for (mdp_suint j = i + 1; j < n; ++j)
        {
          mdp_real alpha = Pi * m_prng.plain();
          mdp_real phi = 2.0f * Pi * m_prng.plain();

          mdp_real z = 2.0f * m_prng.plain() - 1.0f;
          mdp_real s = std::sqrt(1.0f - z * z);

          mdp_real sa, ca;
          mdp_real sp, cp;

#if 0 // defined(__GNUC__)
          __builtin_sincosf(alpha, &sa, &ca);
          __builtin_sincosf(phi, &sp, &cp);
#else
          sa = std::sin(alpha);
          ca = std::cos(alpha);
          sp = std::sin(phi);
          cp = std::cos(phi);
#endif

          mdp_real a0 = ca;
          mdp_real a1 = sa * s * cp;
          mdp_real a2 = sa * s * sp;
          mdp_real a3 = sa * z;

          mdp_complex *row_i = &M(i, 0);
          mdp_complex *row_j = &M(j, 0);

          mdp_complex A(a0, a3);
          mdp_complex B(a2, a1);
          mdp_complex C(-a2, a1);
          mdp_complex D(a0, -a3);

          for (mdp_suint k = 0; k < n; ++k)
          {
            mdp_complex xi = row_i[k];
            mdp_complex xj = row_j[k];

            row_i[k] = A * xi + B * xj;
            row_j[k] = C * xi + D * xj;
          }
        }
      }

      return M;
    }

    /// skip n numbers from the sequence
    void skip(mdp_uint n)
    {
      while (n--) [[likely]]
        m_prng.plain();
    }
  } mdp_global_random; /// the global random number generator
} // namespace MDP

#endif /* MDP_RANDOM_ */
