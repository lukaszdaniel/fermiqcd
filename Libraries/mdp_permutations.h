/////////////////////////////////////////////////////////////////
/// @file mdp_permutations.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Functions to compute permutations
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_PERMUTATIONS_
#define MDP_PERMUTATIONS_

#include <memory>
#include <algorithm>
#include "mdp_global_vars.h"

// this is my favourite piece of code...

namespace MDP
{
  // this is just n!
  constexpr mdp_int mdp_permutations(int n)
  {
    if (n < 0)
      return 0;

    mdp_int a = 1;
    for (; n > 1; --n)
      a *= n;

    return a;
  }

  /// Returns i-th element of the k-th permutations of n numbers
  /// For example if n=4
  /// [0123] k=0
  /// [0132] k=1
  /// ...
  /// [3210] k=23
  /// Returns -1 on error when (i >= n || k >= n_permutations(n))
  int mdp_permutation(int n, int k, int i)
  {
    if (i >= n || k >= mdp_permutations(n))
      return -1;

    std::vector<mdp_int> fact(n + 1, 1);
    for (int j = 1; j <= n; ++j)
      fact[j] = fact[j - 1] * j;

    std::vector<bool> used(n, false);

    for (int pos = 0; pos < n; ++pos)
    {
      mdp_int f = fact[n - 1 - pos];
      int idx = k / f;
      k %= f;

      // find idx-th unused element
      int count = -1;
      for (int val = 0; val < n; ++val)
      {
        if (!used[val])
          ++count;

        if (count == idx)
        {
          used[val] = true;

          if (pos == i)
            return val;

          break;
        }
      }
    }

    return -1;
  }
} // namespace MDP

#endif /* MDP_PERMUTATIONS_ */
