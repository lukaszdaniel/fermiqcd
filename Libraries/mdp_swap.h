/////////////////////////////////////////////////////////////////
/// @file mdp_swap.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains swap function
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_SWAP_
#define MDP_SWAP_

#include <algorithm>
#include "mdp_global_vars.h"

namespace MDP
{
  template <class T>
  void swap(T *a, T *b, mdp_uint n)
  {
    for (mdp_uint i = 0; i < n; i++)
    {
      std::swap(a[i], b[i]);
    }
  }
} // namespace MDP

#endif /* MDP_SWAP_ */
