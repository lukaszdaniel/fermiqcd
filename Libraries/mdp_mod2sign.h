/////////////////////////////////////////////////////////////////
/// @file mdp_mod2sign.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains function mdp_mod2sign
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_MOD2SIGN_
#define MDP_MOD2SIGN_

/// Returns +1 is x%2==0 -1 otherwise
int mdp_mod2sign(int x)
{
  if ((x % 2) == 0)
    return +1;
  else
    return -1;
}

#endif /* MDP_MOD2SIGN_ */
