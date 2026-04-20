// Program: example05.cpp
#include "mdp.h"

using namespace MDP;

int main()
{
  mdp_random random;
  mdp_uint bin[10];
  mdp_real x;
  for (mdp_suint n = 0; n < 10; n++)
    bin[n] = 0;
  for (mdp_suint i = 0; i < 1000; i++)
  {
    x = random.gaussian();
    for (mdp_suint n = 0; n < 10; n++)
      if ((x >= 0.5 * n) && (x < 0.5 * (n + 1)))
        bin[n]++;
  }
  for (mdp_suint n = 0; n < 10; n++)
    std::cout << "bin[" << n << "] = " << bin[n] << "\n";
  return 0;
}
