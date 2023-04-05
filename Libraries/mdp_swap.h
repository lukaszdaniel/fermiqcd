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

namespace MDP
{
  template <class T>
  void swap(T &a, T &b)
  {
    T c;
    c = a;
    a = b;
    b = c;
  }

  template <class T>
  void swap(T *a, T *b, int n)
  {
    int i;
    T c;
    for (i = 0; i < n; i++)
    {
      c = a[i];
      a[i] = b[i];
      b[i] = c;
    }
  }
} // namespace MDP

#endif /* MDP_SWAP_ */
