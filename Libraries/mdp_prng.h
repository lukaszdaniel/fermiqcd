/////////////////////////////////////////////////////////////////
/// @file mdp_prng.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Class mdp_prng (the random number generator of MDP)
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_PRNG_
#define MDP_PRNG_

#include <cmath>
#include "mdp_global_vars.h"
#include "mdp_matrix.h"
#include "mdp_communicator.h"

namespace MDP
{
  /// @brief Marsaglia's random number generator (same as UKQCD)
  ///
  /// You should not instantiate this class because:
  /// - there is a global object mdp_random
  /// - each field "lattice" has a parallel generator "lattice.random(x)"
  /// Example:
  /// @verbatim
  ///    // print a uniform number in (0,1)
  ///    cout << mdp_random.plain() << endl;
  ///    // print a gaussian number
  ///    cout << mdp_random.gaussian() << endl;
  ///    // print a random SU(10) matrix
  ///    cout << mdp_random.SU(10) << endl;
  /// @endverbatim
  class mdp_prng
  {
  private:
    static constexpr float cd = 7654321.0f / 16777216.0f;
    static constexpr float cm = 16777213.0f / 16777216.0f;
    static constexpr unsigned long N = 97;
    float m_c;
    int m_ui;
    int m_uj;
    std::array<float, N> m_u;
    bool m_has_gauss = false;
    float m_gauss_cache = 0;

  public:
    mdp_prng(mdp_int k = 0) : m_c(362436.0f / 16777216.0f), m_ui(N - 1), m_uj(32)
    {
      // c = 362436.0f / 16777216.0f;

      // ui = 97; /*  There is a bug in the original Fortran version */
      // uj = 33; /*  of UNI -- i and j should be SAVEd in UNI()     */

      if (k == 0)
        initialize(ME);
    }

    /// return a uniform random number in (0,1)
    inline float plain() noexcept
    {
      float luni = m_u[m_ui] - m_u[m_uj];
      if (luni < 0.0f)
        luni += 1.0f;

      m_u[m_ui] = luni;

      m_ui = (m_ui == 0) ? N - 1 : m_ui - 1;
      m_uj = (m_uj == 0) ? N - 1 : m_uj - 1;

      m_c -= cd;
      if (m_c < 0.0f)
        m_c += cm;

      luni -= m_c;
      if (luni < 0.0f)
        luni += 1.0f;

      return luni;
    }

    inline double uniform() noexcept
    {
      return static_cast<double>(plain());
    }

    ///////////////////////////////////////////////////////////////////////
    //      initialize: this takes a single integer in the range
    //    0 <= ijkl <= 900 000 000
    //  and produces four smaller integers needed for start. It is
    //  based on the ideas contained in the RMARIN subroutine in
    //    F. James, "A Review of Pseudorandom Number Generators",
    //      Comp. Phys. Commun. Oct 1990, p.340
    //  To reduce the modifications to the existing code, seed_uni now
    //  takes the role of a preprocessor for rstart.
    //
    //  This is useful for the parallel version of the code as James
    //  states that any integer ijkl will produce a statistically
    //  independent sequence of random numbers.
    //
    //     Very funny.
    //     If that statement was worth anything he would have provided
    //     a proof to go with it. spb 12/12/90
    ///////////////////////////////////////////////////////////////////////

    void initialize(mdp_int ijkl)
    {
      if ((ijkl < 0) || (ijkl > 900000000))
        error("Wrong initialization for random number generator");

      int ij = ijkl / 30082;
      int kl = ijkl % 30082;

      int i = ((ij / 177) % 177) + 2;
      int j = (ij % 177) + 2;
      int k = ((kl / 169) % 178) + 1;
      int l = kl % 169;
      if ((i <= 0) || (i > 178))
        error("Wrong initialization for random number generator");
      if ((j <= 0) || (j > 178))
        error("Wrong initialization for random number generator");
      if ((k <= 0) || (k > 178))
        error("Wrong initialization for random number generator");
      if ((l < 0) || (l > 168))
        error("Wrong initialization for random number generator");
      if (i == 1 && j == 1 && k == 1)
        error("Wrong initialization for random number generator");

      for (unsigned int ii = 0; ii < m_u.size(); ii++)
      {
        float s = 0.0f;
        float t = 0.5f;

        for (unsigned int jj = 0; jj < 24; jj++)
        {
          int m = (((i * j) % 179) * k) % 179;

          i = j;
          j = k;
          k = m;

          l = (53 * l + 1) % 169;

          if ((l * m) % 64 >= 32)
            s += t;

          t *= 0.5f;
        }

        m_u[ii] = s;
      }
    }

    /// returns a gaussian random number
    float gaussian(float sigma = 1.0f) noexcept
    {
      if (m_has_gauss)
      {
        m_has_gauss = false;
        return m_gauss_cache * sigma;
      }

      float u1 = plain();
      float u2 = plain();

      float r = std::sqrt(-2.0f * std::log(u1));
      float theta = 2.0f * Pi * u2;

      m_gauss_cache = r * std::cos(theta);
      m_has_gauss = true;

      return sigma * r * std::sin(theta);
    }

    /// draws a random float in (0,1) from a distribution using accept-reject
    double distribution(float (*fp)(float, void *), void *a = nullptr)
    {
      float x, y;
      do
      {
        x = plain();
        y = plain();
      } while ((*fp)(x, a) >= y);
      return x;
    }

    /// returns a random SU(n) matrix using Cabibbo-Marinari
    mdp_matrix SU(int n) noexcept
    {
      if (n == 1)
      {
        mdp_matrix tmp(1, 1);
        float alpha = 2.0 * Pi * plain();
        tmp(0, 0) = mdp_complex(std::cos(alpha), std::sin(alpha));
        return tmp;
      }

      mdp_matrix M = mdp_identity(n);

      for (int i = 0; i < n - 1; ++i)
      {
        for (int j = i + 1; j < n; ++j)
        {
          float alpha = Pi * plain();
          float phi = 2.0f * Pi * plain();

          float z = 2.0f * plain() - 1.0f;
          float s = std::sqrt(1.0f - z * z);

          float sa, ca;
          float sp, cp;

#if defined(__GNUC__)
          __builtin_sincosf(alpha, &sa, &ca);
          __builtin_sincosf(phi, &sp, &cp);
#else
          sa = std::sin(alpha);
          ca = std::cos(alpha);
          sp = std::sin(phi);
          cp = std::cos(phi);
#endif

          float a0 = ca;
          float a1 = sa * s * cp;
          float a2 = sa * s * sp;
          float a3 = sa * z;

          mdp_complex *row_i = &M(i, 0);
          mdp_complex *row_j = &M(j, 0);

          mdp_complex A(a0, a3);
          mdp_complex B(a2, a1);
          mdp_complex C(-a2, a1);
          mdp_complex D(a0, -a3);

          for (int k = 0; k < n; ++k)
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
    void skip(size_t n)
    {
      while (n--) [[likely]]
        plain();
    }
  } mdp_random; /// the global random number generator
} // namespace MDP

#endif /* MDP_PRNG_ */
