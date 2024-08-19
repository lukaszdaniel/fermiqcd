/////////////////////////////////////////////////////////////////
/// @file fermiqcd_check_differences.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Constants and parameters used by FermiQCD
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_CHECK_DIFFERENCES_
#define FERMIQCD_CHECK_DIFFERENCES_

#include "mdp_field.h"
namespace MDP
{
  // compares two fields and returns the max distance between
  // equivalent components.
  float check_differences(mdp_field<mdp_complex> &chi,
                          mdp_field<mdp_complex> &psi)
  {
    begin_function("check_differences");

    mdp_int i_min = psi.physical_local_start();
    mdp_int i_max = psi.physical_local_stop();
    float max = 0, tmp;
    if (&chi.lattice() != &psi.lattice())
      error("check_differences()\nFields defined on different lattices");
    for (mdp_int i = i_min; i < i_max; i++)
    {
      tmp = abs(chi[i] - psi[i]);
      if (tmp > max)
        max = tmp;
    }
    mdp << "Fields agree/disagree within precision=" << max << "\n";

    end_function("check_differences");
    return max;
  }
} // namespace MDP

#endif /* FERMIQCD_CHECK_DIFFERENCES_ */
