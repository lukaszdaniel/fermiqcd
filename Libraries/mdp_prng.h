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
    float u[98];
    float c;
    float cd;
    float cm;
    int ui;
    int uj;

  public:
    /// return a uniform random number in (0,1)
    float plain()
    {
      float luni; /* local variable for Float */
      luni = u[ui] - u[uj];
      if (luni < 0.0)
        luni += 1.0;
      u[ui] = luni;
      if (--ui == 0)
        ui = 97;
      if (--uj == 0)
        uj = 97;
      if ((c -= cd) < 0.0)
        c += cm;
      if ((luni -= c) < 0.0)
        luni += 1.0;
      return ((float)luni);
    }

    ///////////////////////////////////////////////////////////////////////
    //      initializei: this takes a single integer in the range
    //		0 <= ijkl <= 900 000 000
    //	and produces the four smaller integers needed for start. It is
    //	based on the ideas contained in the RMARIN subroutine in
    //		F. James, "A Review of Pseudorandom Number Generators",
    //			Comp. Phys. Commun. Oct 1990, p.340
    //	To reduce the modifications to the existing code, seed_uni now
    //	takes the role of a preprocessor for rstart.
    //
    //	This is useful for the parallel version of the code as James
    //	states that any integer ijkl will produce a statistically
    //	independent sequence of random numbers.
    //
    //     Very funny.
    //     If that statement was worth anything he would have provided
    //     a proof to go with it. spb 12/12/90
    ///////////////////////////////////////////////////////////////////////

    void initialize(mdp_int ijkl)
    {
      int i, j, k, l, ij, kl;
      if ((ijkl < 0) || (ijkl > 900000000))
        error("Wrong initialization for random number generator");
      ij = ijkl / 30082;
      kl = ijkl - (30082 * ij);
      i = ((ij / 177) % 177) + 2;
      j = (ij % 177) + 2;
      k = ((kl / 169) % 178) + 1;
      l = kl % 169;
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
      int ii, jj, m;
      float s, t;
      for (ii = 1; ii <= 97; ii++)
      {
        s = 0.0;
        t = 0.5;
        for (jj = 1; jj <= 24; jj++)
        {
          m = ((i * j % 179) * k) % 179;
          i = j;
          j = k;
          k = m;
          l = (53 * l + 1) % 169;
          if (l * m % 64 >= 32)
            s += t;
          t *= 0.5;
        }
        u[ii] = s;
      }
      c = 362436.0 / 16777216.0;
      cd = 7654321.0 / 16777216.0;
      cm = 16777213.0 / 16777216.0;
      ui = 97; /*  There is a bug in the original Fortran version */
      uj = 33; /*  of UNI -- i and j should be SAVEd in UNI()     */
    }

    mdp_prng(mdp_int k = 0)
    {
      if (k == 0)
        initialize(ME);
    }

    /// returns a gaussian random number
    float gaussian(float sigma = 1)
    {
      static int i = 0;
      static float r = 0;

      if (i == 0)
      {
        float x = (float)std::sqrt(-2.0 * std::log(plain()));
        float y = (float)2.0 * Pi * plain();
        i = 1;
        r = sigma * x * ((float)std::cos(y));
        return sigma * x * ((float)std::sin(y));
      }
      else
      {
        i = 0;
        return r;
      }
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
    mdp_matrix SU(int n)
    {
      mdp_matrix tmp, small;
      int i, j;
      float alpha, sin_alpha;
      float a0, a1, a2, a3;
      float phi, cos_theta, sin_theta;

      if (n == 1)
      {
        tmp.dimension(1, 1);
        alpha = 2.0 * Pi * plain();
        tmp(0, 0) = mdp_complex(std::cos(alpha), std::sin(alpha));
        return tmp;
      }

      tmp = mdp_identity(n);

      for (i = 0; i < (n - 1); i++)
        for (j = i + 1; j < n; j++)
        {
          alpha = Pi * plain();
          phi = 2.0 * Pi * plain();
          cos_theta = 2.0 * plain() - 1.0;
          sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
          sin_alpha = std::sin(alpha);
          a0 = std::cos(alpha);
          a1 = sin_alpha * sin_theta * std::cos(phi);
          a2 = sin_alpha * sin_theta * std::sin(phi);
          a3 = sin_alpha * cos_theta;
          small = mdp_identity(n);
          small(i, i) = mdp_complex(a0, a3);
          small(i, j) = mdp_complex(a2, a1);
          small(j, i) = mdp_complex(-a2, a1);
          small(j, j) = mdp_complex(a0, -a3);
          tmp = small * tmp;
        }
      return tmp;
    }

    /// skip n numbers from the sequence
    void skip(int n)
    {
      for (int i = 0; i < n; i++)
        plain();
    }
  } mdp_random; /// the global random number generator
} // namespace MDP

#endif /* MDP_PRNG_ */
