/// @file mdp_epsilon.h
/// @version 2017-11-23
/// @author Lukasz Daniel <lukasz.daniel@gmail.com>
///
/// Contains declaration of the epsilon function
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MY_EPSILON_
#define MY_EPSILON_

namespace MDP
{
  int epsilon(const int i, const int j, const int k)
  {

    int sign = -1;

    if (i == k || i == j || j == k)
      sign = 0;
    else if ((i < j && j < k) || (k < i && i < j) || (j < k && k < i))
      sign = 1;

    return sign;
  }
} // namespace MDP

#endif /* MY_EPSILON_ */
