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

namespace MDP
{
  /// Set the complex field components of chi to be gaussian random numbers
  /// with mean=0 and sigma=1 (useful for stochastic propagators).
  /// can choose parity=EVEN, ODD or EVENODD
  void set_random(mdp_field<mdp_complex> &chi,
                  int parity = EVENODD)
  {
    int i;
    int i_max = chi.size_per_site();
    mdp_site x(chi.lattice());
    forallsitesofparity(x, parity)
    {
      for (i = 0; i < i_max; i++)
      {
        *chi.address(x, i) = chi.lattice().random(x).gaussian();
      }
    }
    chi.update();
  }

  /// Set the complex field components of chi to be gaussian random numbers
  /// on the wall identified by t
  /// with mean=0 and sigma=1 (useful for stochastic propagators).
  /// can choose parity=EVEN, ODD or EVENODD
  /// attention! does not set to zero other timeslices!!!
  void set_wall_random(mdp_field<mdp_complex> &chi,
                       int t = 0,
                       int parity = EVENODD)
  {
    int i;
    int i_max = chi.size_per_site();
    mdp_site x(chi.lattice());
    forallsitesofparity(x, parity)
    {
      if (x(0) == t)
        for (i = 0; i < i_max; i++)
        {
          *chi.address(x, i) = chi.lattice().random(x).gaussian();
        }
    }
    chi.update();
  }

  /// Set the complex field components of chi tozero.
  /// can choose parity=EVEN, ODD or EVENODD
  void set_zero(mdp_field<mdp_complex> &chi,
                int parity = EVENODD)
  {
    int i;
    int i_max = chi.size_per_site();
    mdp_site x(chi.lattice());
    forallsitesofparity(x, parity)
    {
      for (i = 0; i < i_max; i++)
      {
        *chi.address(x, i) = mdp_complex(0, 0);
      }
    }
    chi.update();
  }
} // namespace MDP

#endif /* FERMIQCD_SET_RANDOM_ */
