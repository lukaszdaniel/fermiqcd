/////////////////////////////////////////////////////////////////
/// @file fermiqcd_su_generators.h
/// @version 2009-03-11
/// @author Simon Catterall and Massimo Di Pierro
///
/// Contains an SU generator class
///
/// Distributed under GPL2 license
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_SU_GENERATORS_
#define FERMIQCD_SU_GENERATORS_

#include <vector>
#include <cmath>
#include "mdp_matrix.h"

namespace MDP
{
  // Example:
  // #include "fermiqcd.h"
  // int main() {
  //  SU_Generators g(3);
  //  for(mdp_uint a=0; a<g.ngenerators; a++)
  //    std::cout << "g=" << g.lambda(a) << std::endl;
  //  return 0;
  // }
  class SU_Generators
  {
  public:
    const mdp_uint n;
    const mdp_uint ngenerators;

  private:
    std::vector<mdp_matrix> m_lambda;

    mdp_matrix build_matrix(mdp_uint a)
    {

      mdp_uint pos1 = (n * (n - 1) / 2);
      mdp_uint pos2 = (n * (n - 1));
      mdp_matrix mat(n, n);
      mdp_complex mult = 0;
      std::vector<mdp_complex> vec(ngenerators);
      for (mdp_uint i = 0; i < ngenerators; i++)
        vec[i] = (i == a) ? 1 : 0;

      for (mdp_uint i = 0, a = 0; i < n; i++)
      {
        mat(i, i) = 0;
        for (mdp_uint j = i + 1; j < n; j++)
        {
          mat(i, j) = 0.5 * (vec[a] - I * vec[pos1 + a]);
          mat(j, i) = 0.5 * (vec[a] + I * vec[pos1 + a]);
          a++;
        }
      }
      for (mdp_uint i = 0; i < n - 1; i++)
      {
        mult = vec[pos2 + i] * (1.0 / std::sqrt(2. + 2. / (1. + i)) / (1. + i));
        for (mdp_uint j = 0; j < i + 1; j++)
        {
          mat(j, j) += mult;
        }
        mat(i + 1, i + 1) -= (1.0 + i) * mult;
      }
      return I * std::sqrt(2) * mat;
    }

  public:
    SU_Generators(mdp_uint N) : n(N), ngenerators(n * n - 1), m_lambda(ngenerators)
    {
      for (mdp_uint a = 0; a < ngenerators; a++)
        m_lambda[a] = build_matrix(a);
    }

    const mdp_matrix &lambda(mdp_uint m) const
    {
      return m_lambda[m];
    }
  };
} // namespace MDP

#endif /* FERMIQCD_SU_GENERATORS_ */
