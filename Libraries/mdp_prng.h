/////////////////////////////////////////////////////////////////
/// @file mdp_prng.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Marsaglia's random number generator
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

  public:
    mdp_prng(mdp_int k = 0) : m_c(362436.0f / 16777216.0f), m_ui(N - 1), m_uj(32)
    {
      if (k == 0)
        initialize(ME);
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

    /**
     * @brief Returns a pseudorandom floating-point value in the range [0, 1].
     *
     * @return Random float in the range [0, 1].
     */
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
  };
} // namespace MDP

#endif /* MDP_PRNG_ */
