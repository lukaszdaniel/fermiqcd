/////////////////////////////////////////////////////////////////
/// @file fermiqcd_set_random.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Function to initialize fields
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_SET_RANDOM_
#define FERMIQCD_SET_RANDOM_

#include "mdp_field.h"

namespace MDP
{
  /**
   * @brief Initialize a complex field with Gaussian random numbers.
   *
   * Sets each complex component of the field @p chi to independent Gaussian
   * random values with mean 0 and standard deviation 1. This is commonly used
   * for stochastic propagators.
   *
   * The initialization can be restricted to a subset of lattice sites based on
   * their parity.
   *
   * @param chi Field to be initialized.
   * @param parity Selection of lattice sites to initialize:
   *        EVEN, ODD, or EVENODD (default: EVENODD).
   */
  void set_random(mdp_complex_field &chi,
                  mdp_parity parity = EVENODD)
  {
    mdp_uint i_max = chi.size_per_site();
    mdp_site x(chi.lattice());

    forallsitesofparity(x, parity)
    {
      for (mdp_uint i = 0; i < i_max; i++)
      {
        *chi.address(x, i) = chi.lattice().random(x).gaussian();
      }
    }
    chi.update();
  }

  /**
   * @brief Initialize a time-slice ("wall") of a complex field with Gaussian random numbers.
   *
   * Sets each complex component of the field @p chi to independent Gaussian
   * random values (mean 0, standard deviation 1), but only on the time-slice
   * specified by @p t.
   *
   * The operation can be restricted to lattice sites of a given parity.
   *
   * @warning This function does NOT reset or modify values on other time-slices.
   *          Existing data outside the selected wall remains unchanged.
   *
   * @param chi Field to be partially initialized.
   * @param t Time coordinate (typically x(0)) identifying the wall to initialize.
   * @param parity Selection of lattice sites to initialize:
   *        EVEN, ODD, or EVENODD (default: EVENODD).
   */
  void set_wall_random(mdp_complex_field &chi,
                       mdp_uint t = 0,
                       mdp_parity parity = EVENODD)
  {
    mdp_uint i_max = chi.size_per_site();
    mdp_site x(chi.lattice());

    forallsitesofparity(x, parity)
    {
      if (x(0) == t)
        for (mdp_uint i = 0; i < i_max; i++)
        {
          *chi.address(x, i) = chi.lattice().random(x).gaussian();
        }
    }
    chi.update();
  }

  /**
   * @brief Set a complex field to zero.
   *
   * Assigns zero to all components of the field @p chi on the selected
   * subset of lattice sites.
   *
   * The operation can be restricted to lattice sites of a given parity.
   *
   * @param chi Field to be cleared.
   * @param parity Selection of lattice sites to reset:
   *        EVEN, ODD, or EVENODD (default: EVENODD).
   */
  void set_zero(mdp_complex_field &chi,
                mdp_parity parity = EVENODD)
  {
    mdp_uint i_max = chi.size_per_site();
    mdp_site x(chi.lattice());
    forallsitesofparity(x, parity)
    {
      for (mdp_uint i = 0; i < i_max; i++)
      {
        *chi.address(x, i) = mdp_complex(0, 0);
      }
    }
    chi.update();
  }
} // namespace MDP

#endif /* FERMIQCD_SET_RANDOM_ */
