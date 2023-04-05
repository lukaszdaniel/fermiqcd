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

namespace MDP
{
  class Instanton4D
  {
  private:
    static int epsilon123(int i, int j, int k)
    {
      if (i == j || j == k || i == k)
        return 0;
      if (i == 1 && j == 2 && k == 3)
        return 1;
      if (i == 1 && j == 3 && k == 2)
        return -1;
      if (i == 2 && j == 1 && k == 3)
        return -1;
      if (i == 2 && j == 3 && k == 1)
        return 1;
      if (i == 3 && j == 1 && k == 2)
        return 1;
      if (i == 3 && j == 2 && k == 1)
        return -1;
      return 0;
    }

  public:
    std::vector<mdp_real> p; // location of the instanton
    int nc;
    int sub_i, sub_j;
    mdp_real charge; // this is 1/g
    mdp_real lambda;
    mdp_matrix eta[4][4];
    Instanton4D(int nc, int sub_i, int sub_j, mdp_real charge, mdp_real lambda, std::vector<mdp_real> &p)
    {
      this->nc = nc;
      this->sub_i = sub_i;
      this->sub_j = sub_j;
      this->lambda = lambda;
      this->charge = charge;
      this->p = p;
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

      for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
        {
          this->eta[mu][nu].dimension(nc, nc);
          this->eta[mu][nu] = 0;
          for (int a = 1; a < 4; a++)
          {
            T = 0;
            T(sub_i, sub_i) = sigma_rot3[a](0, 0);
            T(sub_i, sub_j) = sigma_rot3[a](0, 1);
            T(sub_j, sub_i) = sigma_rot3[a](1, 0);
            T(sub_j, sub_j) = sigma_rot3[a](1, 1);
            if (a > 0 && mu > 0 && nu > 0)
            {
              this->eta[mu][nu] += epsilon123(a, mu, nu) * T;
            }
            else if (mu == 0 && a > 0 && a == nu)
            {
              this->eta[mu][nu] -= T; // T[a]*(-1)*delta(a,nu)
            }
            else if (nu == 0 && a > 0 && a == mu)
            {
              this->eta[mu][nu] += T; // T[a]*delta(a,mu)
            }
          }
        }
      /*
      for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
        {
          this->eta[mu][nu].dimension(nc, nc);
          for (int i = 0; i < nc; i++)
            for (int j = 0; j < nc; j++)
              if (i == j && (i != sub_i && i != sub_j))
              {
                this->eta[mu][nu](i, j) = 1;
              }
              else if ((i != sub_i && i != sub_j) || (j != sub_i && j != sub_j))
              {
                this->eta[mu][nu](i, j) = 0;
              }
              else if (mu == 0 && nu == 0)
              {
                this->eta[mu][nu](i, j) = 0;
              }
              else
              {
                int i0 = (i == sub_i ? 0 : 1);
                int j0 = (j == sub_i ? 0 : 1);
                for (int a = 1; a < 4; a++)
                {
                  if (mu == 0)
                  {
                    // this->eta[mu][nu](i,j) += ((a==nu)?-1:0);
                    this->eta[mu][nu](i, j) += sigma[a](i0, j0) * ((a == nu) ? -1 : 0);
                  }
                  else if (nu == 0)
                  {
                    // this->eta[mu][nu](i,j) += ((a==mu)?+1:0);
                    this->eta[mu][nu](i, j) += sigma[a](i0, j0) * ((a == mu) ? +1 : 0);
                  }
                  else
                  {
                    // this->eta[mu][nu](i,j) += epsilon123(a,mu,nu);
                    this->eta[mu][nu](i, j) += sigma[a](i0, j0) * epsilon123(a, mu, nu);
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
        v[nu] = x(nu) - this->p[nu];
      float d2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
      mdp_matrix A(this->nc, this->nc);
      A = 0;
      for (int nu = 0; nu < 4; nu++)
        A += this->eta[mu][nu] * v[nu];
      return (2.0 * this->charge / (d2 + this->lambda * this->lambda)) * A;
    }
  };
} // namespace MDP

#endif /* FERMIQCD_INSTANTON4D_ */
