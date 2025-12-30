/////////////////////////////////////////////////////////////////
/// @file mdp_vector.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_vector
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_VECTOR_
#define MDP_VECTOR_

#include <array>
#include "mdp_global_vars.h"

namespace MDP
{
  /// @brief discrete vectors to navigate on a lattice
  ///
  class mdp_vector
  {
  private:
    std::array<int, 10> m_x;

  public:
    mdp_vector()
    {
      m_x.fill(0);
    }

    mdp_vector(int x0, int x1 = 0, int x2 = 0, int x3 = 0, int x4 = 0,
               int x5 = 0, int x6 = 0, int x7 = 0, int x8 = 0, int x9 = 0)
    {
      m_x[0] = x0;
      m_x[1] = x1;
      m_x[2] = x2;
      m_x[3] = x3;
      m_x[4] = x4;
      m_x[5] = x5;
      m_x[6] = x6;
      m_x[7] = x7;
      m_x[8] = x8;
      m_x[9] = x9;
    }

    // Overloading [] operator to access elements in array style
    int operator[](mdp_suint index) const
    {
      return m_x[index];
    }
  };

  mdp_vector binary2versor(mdp_int a)
  {
    mdp_vector v((a)&0x1,
                 (a >> 1) & 0x1,
                 (a >> 2) & 0x1,
                 (a >> 3) & 0x1,
                 (a >> 4) & 0x1,
                 (a >> 5) & 0x1,
                 (a >> 6) & 0x1,
                 (a >> 7) & 0x1,
                 (a >> 8) & 0x1,
                 (a >> 9) & 0x1);
    return v;
  }

  int versor2binary(int x0, int x1 = 0, int x2 = 0, int x3 = 0, int x4 = 0,
                           int x5 = 0, int x6 = 0, int x7 = 0, int x8 = 0, int x9 = 0)
  {
#ifdef CHECK_ALL
    if ((fabs(0.5 - x0) > 1) || (fabs(0.5 - x1) > 1) ||
        (fabs(0.5 - x2) > 1) || (fabs(0.5 - x3) > 1) ||
        (fabs(0.5 - x4) > 1) || (fabs(0.5 - x5) > 1) ||
        (fabs(0.5 - x6) > 1) || (fabs(0.5 - x7) > 1) ||
        (fabs(0.5 - x8) > 1) || (fabs(0.5 - x9) > 1))
      error("versor2binary");
#endif
    return x0 + 2 * x1 + 4 * x2 + 8 * x3 + 16 * x4 + 32 * x5 + 64 * x6 + 128 * x7 + 256 * x8 + 512 * x9;
  }

  mdp_int vector2binary(mdp_vector v)
  {
#ifdef CHECK_ALL
    if ((fabs(0.5 - v[0]) > 1) || (fabs(0.5 - v[1]) > 1) ||
        (fabs(0.5 - v[2]) > 1) || (fabs(0.5 - v[3]) > 1) ||
        (fabs(0.5 - v[4]) > 1) || (fabs(0.5 - v[5]) > 1) ||
        (fabs(0.5 - v[6]) > 1) || (fabs(0.5 - v[7]) > 1) ||
        (fabs(0.5 - v[8]) > 1) || (fabs(0.5 - v[9]) > 1) ||
        (fabs(0.5 - v[2]) > 1))
      error("vector2binary");
#endif
    return v[0] + 2 * v[1] + 4 * v[2] + 8 * v[3] + 16 * v[4] + 32 * v[5] + 64 * v[6] + 128 * v[7] + 256 * v[8] + 512 * v[9];
  }
} // namespace MDP

#endif /* MDP_VECTOR_ */
