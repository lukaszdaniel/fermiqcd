/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_routines.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various gauge multiplication routines (some call SSE/SSE2 macors)
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
  mdp_matrix staple(gauge_field &U, mdp_site x,
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

  mdp_matrix staple(gauge_field &U, mdp_site x, int mu)
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

  mdp_matrix staple_H(gauge_field &U, mdp_site x,
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

  mdp_matrix staple_H(gauge_field &U, mdp_site x, int mu)
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

  mdp_matrix staple_H_unisotropic(gauge_field &U, mdp_site x, int mu, mdp_real zeta)
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

  mdp_matrix staple_0i_H(gauge_field &U, mdp_site x, int mu)
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

  mdp_matrix staple_ij_H(gauge_field &U, mdp_site x, int mu)
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
  mdp_matrix plaquette(gauge_field &U, mdp_site x, int mu, int nu)
  {
    mdp_matrix tmp(U.nc(), U.nc());
    tmp = U(x, mu) * U(x + mu, nu) * hermitian(U(x + nu, mu)) * hermitian(U(x, nu));
    return tmp;
  }

  mdp_real SymanzikAction(gauge_field &U)
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

#if 0
  // obsolete stuff

  void fast_mul_AB_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                        int &ni)
  {
    int ink, inj;
    for (int i = 0; i < ni; i++)
    {
      ink = i * ni;
      inj = i * ni;
      for (int k = 0; k < ni; k++)
      {
        c[ink + k] = a[inj] * b[k];
        for (int j = 1; j < ni; j++)
          c[ink + k] += a[inj + j] * b[j * ni + k];
      }
    }
  }

  void fast_mul_ABH_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                         int &ni)
  {
    int ink, inj;
    for (int i = 0; i < ni; i++)
    {
      ink = i * ni;
      inj = i * ni;
      for (int k = 0; k < ni; k++)
      {
        c[ink + k] = a[inj] * conj(b[k * ni]);
        for (int j = 1; j < ni; j++)
          c[ink + k] += a[inj + j] * conj(b[k * ni + j]);
      }
    }
  }

  void fast_mul_AHB_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                         int &ni)
  {
    int ink;
    for (int i = 0; i < ni; i++)
    {
      ink = i * ni;
      for (int k = 0; k < ni; k++)
      {
        c[ink + k] = conj(a[i]) * b[k];
        for (int j = 1; j < ni; j++)
          c[ink + k] += conj(a[j * ni + i]) * b[j * ni + k];
      }
    }
  }

  void fast_mul_AHBH_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                          int &ni)
  {
    int i, j, k, ink, inj;
    for (i = 0; i < ni; i++)
    {
      ink = i * ni;
      inj = i * ni;
      for (k = 0; k < ni; k++)
      {
        c[ink + k] = conj(a[i] * b[k * ni]);
        for (j = 1; j < ni; j++)
          c[ink + k] += conj(a[j * ni + i] * b[k * ni + j]);
      }
    }
  }

  void fast_mul_AB_addto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                           int &ni)
  {
    int i, j, k, ink, inj;
    for (i = 0; i < ni; i++)
    {
      ink = i * ni;
      inj = i * ni;
      for (k = 0; k < ni; k++)
      {
        c[ink + k] += a[inj] * b[k];
        for (j = 1; j < ni; j++)
          c[ink + k] += a[inj + j] * b[j * ni + k];
      }
    }
  }

  void fast_mul_ABH_addto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                            int &ni)
  {
    int i, j, k, ink, inj;
    for (i = 0; i < ni; i++)
    {
      ink = i * ni;
      inj = i * ni;
      for (k = 0; k < ni; k++)
      {
        c[ink + k] += a[inj] * conj(b[k * ni]);
        for (j = 1; j < ni; j++)
          c[ink + k] += a[inj + j] * conj(b[k * ni + j]);
      }
    }
  }

  void fast_mul_AHB_addto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                            int &ni)
  {
    int i, j, k, ink;
    for (i = 0; i < ni; i++)
    {
      ink = i * ni;
      for (k = 0; k < ni; k++)
      {
        c[ink + k] += conj(a[i]) * b[k];
        for (j = 1; j < ni; j++)
          c[ink + k] += conj(a[j * ni + i]) * b[j * ni + k];
      }
    }
  }

  void fast_add_AB_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                        int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] = a[i] + b[i];
  }

  void fast_add_AB_addto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                           int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] += a[i] + b[i];
  }

  void fast_add_aAbB_to_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
                          mdp_complex *c, int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] = c1 * a1[i] + c2 * a2[i];
  }

  void fast_add_aAbB_addto_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
                             mdp_complex *c, int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] += c1 * a1[i] + c2 * a2[i];
  }

  void fast_add_aAbB_to_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
                          mdp_complex c3, mdp_complex *a3, mdp_complex *c, int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] = c1 * a1[i] + c2 * a2[i] + c3 * a3[i];
  }

  void fast_add_aAbB_addto_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
                             mdp_complex c3, mdp_complex *a3, mdp_complex *c, int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] += c1 * a1[i] + c2 * a2[i] + c3 * a3[i];
  }

  void fast_scalar_mul_AB_to_C(mdp_complex a, mdp_complex *b, mdp_complex *c,
                               int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] = a * b[i];
  }

  void fast_scalar_mul_AB_addto_C(mdp_complex a, mdp_complex *b, mdp_complex *c,
                                  int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] += a * b[i];
  }

  void fast_mul_AB_subto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                           int &ni)
  {
    int i, j, k, ink, inj;
    for (i = 0; i < ni; i++)
    {
      ink = i * ni;
      inj = i * ni;
      for (k = 0; k < ni; k++)
      {
        c[ink + k] += a[inj] * b[k];
        for (j = 1; j < ni; j++)
          c[ink + k] -= a[inj + j] * b[j * ni + k];
      }
    }
  }

  void fast_mul_ABH_subto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                            int &ni)
  {
    int i, j, k, ink, inj;
    for (i = 0; i < ni; i++)
    {
      ink = i * ni;
      inj = i * ni;
      for (k = 0; k < ni; k++)
      {
        c[ink + k] += a[inj] * conj(b[k * ni]);
        for (j = 1; j < ni; j++)
          c[ink + k] -= a[inj + j] * conj(b[k * ni + j]);
      }
    }
  }

  void fast_mul_AHB_subto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                            int &ni)
  {
    int i, j, k, ink;
    for (i = 0; i < ni; i++)
    {
      ink = i * ni;
      for (k = 0; k < ni; k++)
      {
        c[ink + k] += conj(a[i]) * b[k];
        for (j = 1; j < ni; j++)
          c[ink + k] -= conj(a[j * ni + i]) * b[j * ni + k];
      }
    }
  }

  void fast_add_AB_subto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c,
                           int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] -= a[i] + b[i];
  }

  void fast_add_aAbB_subto_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
                             mdp_complex *c, int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] -= c1 * a1[i] + c2 * a2[i];
  }

  void fast_add_aAbB_subto_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
                             mdp_complex c3, mdp_complex *a3, mdp_complex *c, int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] -= c1 * a1[i] + c2 * a2[i] + c3 * a3[i];
  }

  void fast_scalar_mul_AB_subto_C(mdp_complex a, mdp_complex *b, mdp_complex *c,
                                  int &ni)
  {
    for (int i = 0; i < ni * ni; i++)
      c[i] -= a * b[i];
  }

  // ///////////////////////////
  // else use fast_mul functions
  // gains a factor 1.45 in time
  // ///////////////////////////

  mdp_matrix staple(gauge_field &U, mdp_site x,
                           int mu, int s1, int nu)
  {
    int nc = U.nc();
    mdp_matrix b1(nc, nc);
    mdp_site y(U.lattice());
    mdp_matrix tmp(nc, nc);
    if (s1 == +1)
    {
      fast_mul_AB_to_C(&U(x, nu, 0, 0), &U(x + nu, mu, 0, 0), &b1(0, 0), nc);
      fast_mul_ABH_to_C(&b1(0, 0), &U(x + mu, nu, 0, 0), &tmp(0, 0), nc);
    }
    else
    {
      y = x - nu;
      fast_mul_AB_to_C(&U(y, mu, 0, 0), &U(y + mu, nu, 0, 0), &b1(0, 0), nc);
      fast_mul_AHB_to_C(&U(y, nu, 0, 0), &b1(0, 0), &tmp(0, 0), nc);
    }
    return tmp;
  }

  mdp_matrix staple(gauge_field &U, mdp_site x, int mu)
  {
    int nc = U.nc();
    mdp_matrix b1(nc, nc);
    mdp_matrix tmp(nc, nc);
    mdp_site y(U.lattice());
    int nu;
    tmp = 0;
    for (nu = 0; nu < U.ndim(); nu++)
      if (nu != mu)
      {
        fast_mul_AB_to_C(&U(x, nu, 0, 0), &U(x + nu, mu, 0, 0), &b1(0, 0), nc);
        fast_mul_ABH_addto_C(&b1(0, 0), &U(x + mu, nu, 0, 0), &tmp(0, 0), nc);
        y = x - nu;
        fast_mul_AB_to_C(&U(y, mu, 0, 0), &U(y + mu, nu, 0, 0), &b1(0, 0), nc);
        fast_mul_AHB_addto_C(&U(y, nu, 0, 0), &b1(0, 0), &tmp(0, 0), nc);
      }
    return tmp;
  }

  mdp_matrix staple_H(gauge_field &U, mdp_site x,
                             int mu, int s1, int nu)
  {
    int nc = U.nc();
    mdp_matrix b1(nc, nc);
    mdp_matrix tmp(nc, nc);
    mdp_site y(U.lattice());
    if (s1 == +1)
    {
      fast_mul_ABH_to_C(&U(x + mu, nu, 0, 0), &U(x + nu, mu, 0, 0), &b1(0, 0), nc);
      fast_mul_ABH_to_C(&b1(0, 0), &U(x, nu, 0, 0), &tmp(0, 0), nc);
    }
    else
    {
      y = x - nu;
      fast_mul_AHB_to_C(&U(y, mu, 0, 0), &U(y, nu, 0, 0), &b1(0, 0), nc);
      fast_mul_AHB_to_C(&U(y + mu, nu, 0, 0), &b1(0, 0), &tmp(0, 0), nc);
    }
    return tmp;
  }

  mdp_matrix staple_H(gauge_field &U, mdp_site x, int mu)
  {
    int nc = U.nc();
    mdp_matrix b1(nc, nc);
    mdp_matrix tmp(nc, nc);
    mdp_site y(U.lattice());
    int nu;
    tmp = 0;
    for (nu = 0; nu < U.ndim(); nu++)
      if (nu != mu)
      {
        fast_mul_ABH_to_C(&U(x + mu, nu, 0, 0), &U(x + nu, mu, 0, 0), &b1(0, 0), nc);
        fast_mul_ABH_addto_C(&b1(0, 0), &U(x, nu, 0, 0), &tmp(0, 0), nc);
        y = x - nu;
        fast_mul_AHB_to_C(&U(y, mu, 0, 0), &U(y, nu, 0, 0), &b1(0, 0), nc);
        fast_mul_AHB_addto_C(&U(y + mu, nu, 0, 0), &b1(0, 0), &tmp(0, 0), nc);
      }
    return tmp;
  }

  mdp_matrix staple_0i_H(gauge_field &U, mdp_site x, int mu)
  {
    int nc = U.nc();
    mdp_matrix b1(nc, nc);
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    int nu;
    if (mu == 0)
    {
      tmp = 0;
      for (nu = 1; nu < U.ndim(); nu++)
      {
        fast_mul_ABH_to_C(&U(x + mu, nu, 0, 0), &U(x + nu, mu, 0, 0), &b1(0, 0), nc);
        fast_mul_ABH_addto_C(&b1(0, 0), &U(x, nu, 0, 0), &tmp(0, 0), nc);
        y = x - nu;
        fast_mul_AHB_to_C(&U(y, mu, 0, 0), &U(y, nu, 0, 0), &b1(0, 0), nc);
        fast_mul_AHB_addto_C(&U(y + mu, nu, 0, 0), &b1(0, 0), &tmp(0, 0), nc);
      }
    }
    else
    {
      nu = 0;
      fast_mul_ABH_to_C(&U(x + mu, nu, 0, 0), &U(x + nu, mu, 0, 0), &b1(0, 0), nc);
      fast_mul_ABH_to_C(&b1(0, 0), &U(x, nu, 0, 0), &tmp(0, 0), nc);
      y = x - nu;
      fast_mul_AHB_to_C(&U(y, mu, 0, 0), &U(y, nu, 0, 0), &b1(0, 0), nc);
      fast_mul_AHB_addto_C(&U(y + mu, nu, 0, 0), &b1(0, 0), &tmp(0, 0), nc);
    }
    return tmp;
  }

  mdp_matrix staple_ij_H(gauge_field &U, mdp_site x, int mu)
  {
    int nc = U.nc();
    mdp_matrix b1(nc, nc);
    mdp_matrix tmp(U.nc(), U.nc());
    mdp_site y(U.lattice());
    int nu;
    tmp = 0;
    if (mu != 0)
      for (nu = 1; nu < U.ndim(); nu++)
        if (nu != mu)
        {
          fast_mul_ABH_to_C(&U(x + mu, nu, 0, 0), &U(x + nu, mu, 0, 0), &b1(0, 0), nc);
          fast_mul_ABH_addto_C(&b1(0, 0), &U(x, nu, 0, 0), &tmp(0, 0), nc);
          y = x - nu;
          fast_mul_AHB_to_C(&U(y, mu, 0, 0), &U(y, nu, 0, 0), &b1(0, 0), nc);
          fast_mul_AHB_addto_C(&U(y + mu, nu, 0, 0), &b1(0, 0), &tmp(0, 0), nc);
        }
    return tmp;
  }

  // ///////////////////////////
  // plaquette
  // //////////////////////////
  mdp_matrix plaquette(gauge_field &U, mdp_site x, int mu, int nu)
  {
    int nc = U.nc();
    mdp_matrix b1(nc, nc);
    mdp_matrix b2(nc, nc);
    mdp_matrix tmp(U.nc(), U.nc());
    fast_mul_AB_to_C(&U(x, mu, 0, 0), &U(x + mu, nu, 0, 0), &b1(0, 0), nc);
    fast_mul_ABH_to_C(&b1(0, 0), &U(x + nu, mu, 0, 0), &b2(0, 0), nc);
    fast_mul_ABH_to_C(&b2(0, 0), &U(x, nu, 0, 0), &tmp(0, 0), nc);
    return tmp;
  }
#endif
} // namespace MDP

#endif /* FERMIQCD_GAUGE_ROUTINES_ */
