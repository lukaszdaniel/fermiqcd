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

// this is my favourite piece of code...

namespace MDP
{
  // this is just n!
  mdp_int mdp_permutations(int n)
  {
    mdp_int a = 1;
    for (; n; n--)
      a *= n;
    return a;
  }

  // this sorts the first k elements of map[] assuming
  // the fisrt k-1 are already sorted
  void mdp_permutation_sort(int map[], int k)
  {
    for (int i = k - 1; i >= 0; i--)
      if (map[i] > map[i + 1])
      {
        std::swap(map[i], map[i + 1]);
      }
  }

  /// Returns j-th element of the k-th permutations of n numbers
  /// For example if n=4
  /// [0123] k=0
  /// [0132] k=1
  /// ...
  /// [3210] k=23
  /// Returns -1 on error when (i>n || k>n_permutations(n))
  int mdp_permutation(int n, int k, int i)
  {
    std::unique_ptr<int[]> map = std::make_unique<int[]>(i + 1);

    if (i > n || k > mdp_permutations(n))
      return -1;

    for (int j = 0; j <= i; j++)
    {
      map[j] = (k % mdp_permutations(n - j)) / mdp_permutations(n - 1 - j);
      mdp_permutation_sort(map.get(), j - 1);
      for (int l = 0; l < j; l++)
        if (map[l] <= map[j])
          map[j]++;
      if (i == j)
      {
        return map[j];
      }
    }
    return -1;
  }
} // namespace MDP

#endif /* MDP_PERMUTATIONS_ */
