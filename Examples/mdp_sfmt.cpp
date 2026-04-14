#include "mdp.h"

using namespace MDP;

int main()
{
  mdp_prng_sfmt r;
  r.initialize(100);
  constexpr int n = 30;
  int c[n];
  long long sum1 = 0, sum2 = 0;
  float mean;
  float sd;
  int i, counter = 0;
  float x;
  for (int j = 0; j < n; j++)
    c[j] = 0;
  while (1)
  {
    x = r.plain();
    i = int(30.0 * x);
    c[i]++;
    counter++;
    sum1 += 1;
    sum2 += 2 * c[i] - 1;
    mean = (1.0 * sum1) / n;
    sd = std::sqrt(std::abs(sum2 / n - mean * mean));
    std::cout << x << "\t" << counter << "\t" << mean << "\t" << sum2 / counter << "\t" << sd / mean << std::endl;
    if (counter == 20)
      break;
  }

  return 0;
}
