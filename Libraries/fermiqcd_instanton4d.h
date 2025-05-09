/////////////////////////////////////////////////////////////////
/// @file fermiqcd_instanton4d.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// class for building a single instanton gauge configuration
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_INSTANTON4D_
#define FERMIQCD_INSTANTON4D_

#include <vector>
#include "mdp_global_vars.h"
#include "mdp_epsilon.h"
#include "mdp_matrix.h"
#include "mdp_site.h"
#include "fermiqcd_gamma_matrices.h"

namespace MDP
{
  class Instanton4D
  {
  public:
    std::vector<mdp_real> p; // location of the instanton
    int nc;
    int sub_i, sub_j;
    mdp_real charge; // this is 1/g
    mdp_real lambda;
    mdp_matrix eta[4][4];
    Instanton4D(int nc_, int sub_i_, int sub_j_, mdp_real charge_, mdp_real lambda_, std::vector<mdp_real> &p_)
    {
      nc = nc_;
      sub_i = sub_i_;
      sub_j = sub_j_;
      lambda = lambda_;
      charge = charge_;
      p = p_;
      mdp_matrix T(nc, nc);
      mdp_matrix sigma_rot1[4];
      mdp_matrix sigma_rot2[4];
      mdp_matrix sigma_rot3[4];
      mdp_real alpha = 0.0, beta = 0.0, gamma = 0.0; // ???

      sigma_rot1[3] = sigma[3];
      sigma_rot1[1] = std::cos(alpha) * sigma[1] + std::sin(alpha) * sigma[2];
      sigma_rot1[2] = -std::sin(alpha) * sigma[1] + std::cos(alpha) * sigma[2];

      sigma_rot2[1] = sigma_rot1[1];
      sigma_rot2[2] = std::cos(beta) * sigma_rot1[2] + std::sin(beta) * sigma_rot1[3];
      sigma_rot2[3] = -std::sin(beta) * sigma_rot1[2] + std::cos(beta) * sigma_rot1[3];

      sigma_rot3[3] = sigma_rot2[3];
      sigma_rot3[1] = std::cos(gamma) * sigma_rot2[1] + std::sin(gamma) * sigma_rot2[2];
      sigma_rot3[2] = -std::sin(gamma) * sigma_rot2[1] + std::cos(gamma) * sigma_rot2[2];

      for (mdp_int mu = 0; mu < 4; mu++)
        for (mdp_int nu = 0; nu < 4; nu++)
        {
          eta[mu][nu].dimension(nc, nc);
          eta[mu][nu] = 0;
          for (int a = 1; a < 4; a++)
          {
            T = 0;
            T(sub_i, sub_i) = sigma_rot3[a](0, 0);
            T(sub_i, sub_j) = sigma_rot3[a](0, 1);
            T(sub_j, sub_i) = sigma_rot3[a](1, 0);
            T(sub_j, sub_j) = sigma_rot3[a](1, 1);
            if (a > 0 && mu > 0 && nu > 0)
            {
              eta[mu][nu] += epsilon(a, mu, nu) * T;
            }
            else if (mu == 0 && a > 0 && a == nu)
            {
              eta[mu][nu] -= T; // T[a]*(-1)*delta(a,nu)
            }
            else if (nu == 0 && a > 0 && a == mu)
            {
              eta[mu][nu] += T; // T[a]*delta(a,mu)
            }
          }
        }
      /*
      for (mdp_int mu = 0; mu < 4; mu++)
        for (mdp_int nu = 0; nu < 4; nu++)
        {
          eta[mu][nu].dimension(nc, nc);
          for (mdp_int i = 0; i < nc; i++)
            for (mdp_int j = 0; j < nc; j++)
              if (i == j && (i != sub_i && i != sub_j))
              {
                eta[mu][nu](i, j) = 1;
              }
              else if ((i != sub_i && i != sub_j) || (j != sub_i && j != sub_j))
              {
                eta[mu][nu](i, j) = 0;
              }
              else if (mu == 0 && nu == 0)
              {
                eta[mu][nu](i, j) = 0;
              }
              else
              {
                int i0 = (i == sub_i ? 0 : 1);
                int j0 = (j == sub_i ? 0 : 1);
                for (int a = 1; a < 4; a++)
                {
                  if (mu == 0)
                  {
                    // eta[mu][nu](i,j) += ((a==nu)?-1:0);
                    eta[mu][nu](i, j) += sigma[a](i0, j0) * ((a == nu) ? -1 : 0);
                  }
                  else if (nu == 0)
                  {
                    // eta[mu][nu](i,j) += ((a==mu)?+1:0);
                    eta[mu][nu](i, j) += sigma[a](i0, j0) * ((a == mu) ? +1 : 0);
                  }
                  else
                  {
                    // eta[mu][nu](i,j) += epsilon(a,mu,nu);
                    eta[mu][nu](i, j) += sigma[a](i0, j0) * epsilon(a, mu, nu);
                  }
                }
              }
        }
      */
    }
    mdp_matrix operator()(mdp_site &x, int mu)
    {
      int v[4];
      // mdp_lattice &lattice=x.lattice();
      for (int nu = 0; nu < 4; nu++)
        v[nu] = x(nu) - p[nu];
      float d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
      mdp_matrix A(nc, nc);
      A = 0;
      for (int nu = 0; nu < 4; nu++)
        A += eta[mu][nu] * v[nu];
      return (2.0 * charge / (d2 + lambda * lambda)) * A;
    }
  };
} // namespace MDP

#endif /* FERMIQCD_INSTANTON4D_ */
