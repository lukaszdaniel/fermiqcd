#include "mdp.h"

using namespace MDP;

int main()
{
  mdp_prng_sfmt r;
  r.initialize(100);
  constexpr mdp_suint n = 30;
  mdp_suint c[n];
  mdp_uint sum1 = 0, sum2 = 0;
  mdp_uint i, counter = 0;
  mdp_real x;

  for (mdp_suint j = 0; j < n; j++)
    c[j] = 0;

  while (1)
  {
    x = r.plain();
    i = mdp_suint(30.0 * x);
    c[i]++;
    counter++;
    sum1 += 1;
    sum2 += 2 * c[i] - 1;
    mdp_real mean = (1.0 * sum1) / n;
    mdp_real sd = std::sqrt(std::abs(sum2 / n - mean * mean));
    std::cout << x << "\t" << counter << "\t" << mean << "\t" << sum2 / counter << "\t" << sd / mean << std::endl;
    if (counter == 20)
      break;
  }

  return 0;
}
