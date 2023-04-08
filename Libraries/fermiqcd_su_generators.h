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

namespace MDP
{
  // SU(n) generator
  // Example:
  // #include "fermiqcd.h"
  // int main() {
  //  SU_Generators g(3);
  //  for(int a=0; a<g.ngenerators; a++)
  //    cout << "g=" << g.lambda[a] << endl;
  //  return 0;
  // }
  class SU_Generators
  {
  public:
    std::vector<mdp_matrix> lambda;
    int n;
    int ngenerators;

    mdp_matrix build_matrix(int a)
    {

      int pos1 = (n * (n - 1) / 2);
      int pos2 = (n * (n - 1));
      mdp_matrix mat(n, n);
      mdp_complex mult = 0;
      std::vector<mdp_complex> vec(ngenerators);
      for (int i = 0; i < ngenerators; i++)
        vec[i] = (i == a) ? 1 : 0;

      for (int i = 0, a = 0; i < n; i++)
      {
        mat(i, i) = 0;
        for (int j = i + 1; j < n; j++)
        {
          mat(i, j) = 0.5 * (vec[a] - I * vec[pos1 + a]);
          mat(j, i) = 0.5 * (vec[a] + I * vec[pos1 + a]);
          a++;
        }
      }
      for (int i = 0; i < n - 1; i++)
      {
        mult = vec[pos2 + i] * (1.0 / std::sqrt(2. + 2. / (1. + i)) / (1. + i));
        for (int j = 0; j < i + 1; j++)
        {
          mat(j, j) += mult;
        }
        mat(i + 1, i + 1) -= (1 + i) * mult;
      }
      return I * std::sqrt(2) * mat;
    }

    SU_Generators(int N)
    {
      n = N;
      ngenerators = (n * n - 1);
      lambda.resize(ngenerators);

      for (int a = 0; a < ngenerators; a++)
        lambda[a] = build_matrix(a);
    }
  };
} // namespace MDP

#endif /* FERMIQCD_SU_GENERATORS_ */
