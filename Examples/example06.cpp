// Program: example06.cpp
#include "mdp.h"

using namespace MDP;

mdp_real Q(mdp_real x, [[maybe_unused]] void *a)
{
  return sin(Pi * x);
}

int main()
{
  mdp_random random;
  constexpr mdp_suint N = 100;
  constexpr mdp_real sigma = 0.3, a_bar = 1;
  mdp_real average = 0;
  for (mdp_suint i = 0; i < N; i++)
  {
    mdp_real a = (sigma * random.gaussian() + a_bar);
    mdp_real b = random.distribution(Q);
    average += a + b;
    std::cout << "average=" << average / (i + 1) << "\n";
  }
  return 0;
}
