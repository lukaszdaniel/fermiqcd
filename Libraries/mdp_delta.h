/////////////////////////////////////////////////////////////////
/// @file mdp_delta.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration delta function
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_DELTA_
#define MDP_DELTA_

/// True if i==j, false otherwise
template <class T>
const bool delta(const T &i, const T &j)
{
  return (i == j) ? true : false;
}

#endif /* MDP_DELTA_ */
