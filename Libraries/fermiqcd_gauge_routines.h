/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_routines.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various gauge multiplication routines
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_GAUGE_ROUTINES_
#define FERMIQCD_GAUGE_ROUTINES_

#include "mdp_matrix.h"
#include "fermiqcd_gauge_field.h"

namespace MDP
{
  // ///////////////////////////
  // If use class mdp_matrix
  // ///////////////////////////
  mdp_matrix staple(const gauge_field &U, mdp_site x,
                           int mu, int s1, int nu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    if (s1 == +1)
    {
      tmp = U(x, nu) * U(x + nu, mu) * hermitian(U(x + mu, nu));
    }
    else
    {
      mdp_site y(U.lattice());
      y = x - nu;
      tmp = hermitian(U(y, nu)) * U(y, mu) * U(y + mu, nu);
    }
    return tmp;
  }

  mdp_matrix staple(const gauge_field &U, mdp_site x, int mu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());

    tmp = 0;
    for (int nu = 0; nu < U.ndim(); nu++)
      if (nu != mu)
      {
        tmp += U(x, nu) * U(x + nu, mu) * hermitian(U(x + mu, nu));
        y = x - nu;
        tmp += hermitian(U(y, nu)) * U(y, mu) * U(y + mu, nu);
      }
    return tmp;
  }

  mdp_matrix staple_H(const gauge_field &U, mdp_site x,
                             int mu, int s1, int nu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    if (s1 == +1)
    {
      tmp = U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
    }
    else
    {
      mdp_site y(U.lattice());
      y = x - nu;
      tmp = hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
    }
    return tmp;
  }

  mdp_matrix staple_H(const gauge_field &U, mdp_site x, int mu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    int nu;
    tmp = 0;
    for (nu = 0; nu < U.ndim(); nu++)
      if (nu != mu)
      {
        tmp += U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
        y = x - nu;
        tmp += hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
      }
    return tmp;
  }

  mdp_matrix staple_H_unisotropic(const gauge_field &U, mdp_site x, int mu, mdp_real zeta)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    int nu;
    mdp_real param;
    tmp = 0;
    for (nu = 0; nu < U.ndim(); nu++)
      if (nu != mu)
      {
        if (nu * mu == 0)
        {
          param = zeta;
        }
        else
        {
          param = 1.0 / zeta;
        }
        tmp += param * U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
        y = x - nu;
        tmp += param * hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
      }
    return tmp;
  }

  mdp_matrix staple_0i_H(const gauge_field &U, mdp_site x, int mu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    int nu;
    if (mu == 0)
    {
      tmp = 0;
      for (nu = 1; nu < U.ndim(); nu++)
      {
        tmp += U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
        y = x - nu;
        tmp += hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
      }
    }
    else
    {
      nu = 0;
      tmp = U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
      y = x - nu;
      tmp += hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
    }

    return tmp;
  }

  mdp_matrix staple_ij_H(const gauge_field &U, mdp_site x, int mu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    int nu;
    tmp = 0;
    if (mu != 0)
      for (nu = 1; nu < U.ndim(); nu++)
        if (nu != mu)
        {
          tmp += U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
          y = x - nu;
          tmp += hermitian(U(y + mu, nu)) * hermitian(U(y, mu)) * U(y, nu);
        }
    return tmp;
  }

  // ///////////////////////////
  // plaquette
  // //////////////////////////
  mdp_matrix plaquette(const gauge_field &U, mdp_site x, int mu, int nu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    tmp = U(x, mu) * U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
    return tmp;
  }

  mdp_real SymanzikAction(const gauge_field &U)
  {
    mdp_real tmp = 0;
    mdp_site x(U.lattice());
    mdp_real c1 = -1.0 / 12; // Symanzik
    //  mdp_real c1=-0.331; //Iwasaki
    //  mdp_real c1=-1.4088; //DBW2

    // U.update();

    forallsites(x)
    {
      for (mdp_int mu = 0; mu < U.ndim() - 1; mu++)
      {
        for (mdp_int nu = mu + 1; nu < U.ndim(); nu++)
        {
          tmp += (1 - 8 * c1) * (1 - real(trace(plaquette(U, x, mu, nu))) / U.nc()) +
                 c1 * (2 - 1.0 / U.nc() * real(trace(U(x, mu) * (U(x + mu, mu) * U((x + mu) + mu, nu) * hermitian(U((x + mu) + nu, mu)) * hermitian(U(x + nu, mu)) * +U(x + mu, nu) * U((x + mu) + nu, nu) * hermitian(U((x + nu) + nu, mu)) * hermitian(U(x + nu, nu))) * hermitian(U(x, nu)))));
        }
      }
    }

    mdp.add(tmp);
    return tmp / (U.lattice().global_volume());
  }
} // namespace MDP

#endif /* FERMIQCD_GAUGE_ROUTINES_ */
