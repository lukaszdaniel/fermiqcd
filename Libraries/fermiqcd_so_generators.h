/////////////////////////////////////////////////////////////////
/// @file fermiqcd_so_generators.h
/// @version 2009-03-11
/// @author Simon Catterall and Massimo Di Pierro
///
/// Contains an SO generator class
///
/// Distributed under GPL2 license
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_SO_GENERATORS_
#define FERMIQCD_SO_GENERATORS_

#include <vector>
#include <cmath>
#include "mdp_matrix.h"

namespace MDP
{
  // Example:
  // #include "fermiqcd.h"
  // int main() {
  //  SO_Generators g(4);
  //  for(int a=0; a<g.ngenerators; a++)
  //    cout << "g=" << g.lambda(a) << endl;
  //  return 0;
  // }
  class SO_Generators
  {
  public:
    const int n;
    const int ngenerators;

  private:
    std::vector<mdp_matrix> m_lambda;

  public:
    SO_Generators(int N) : n(N), ngenerators(n * (n - 1) / 2), m_lambda(ngenerators)
    {
      mdp_matrix temp(n, n);
      mdp_complex z = 1.0 / std::sqrt(2.0);
      int k = 0;
      for (int j = 0; j < n - 1; j++)
        for (int i = j + 1; i < n; i++)
        {
          temp = 0;
          temp(i, j) = z;
          temp(j, i) = -z;
          m_lambda[k++] = temp;
        }
    }

    const mdp_matrix &lambda(int m) const
    {
      return m_lambda[m];
    }
  };
} // namespace MDP

#endif /* FERMIQCD_SO_GENERATORS_ */
