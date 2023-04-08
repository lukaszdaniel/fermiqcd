// Program: example06.cpp
#include "mdp.h"

using namespace MDP;

float Q(float x, void *a)
{
   return sin(Pi * x);
}

int main()
{
   mdp_prng random;
   constexpr mdp_suint N = 100;
   float a, b, average = 0, sigma = 0.3, a_bar = 1;
   for (mdp_suint i = 0; i < N; i++)
   {
      a = (sigma * random.gaussian() + a_bar);
      b = random.distribution(Q);
      average += a + b;
      std::cout << "average=" << average / (i + 1) << "\n";
   }
   return 0;
}
