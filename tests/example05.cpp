// Program: example05.cpp
#include "mdp.h"

using namespace MDP;

int main()
{
   mdp_prng random;
   int i, n, bin[10];
   float x;
   for (n = 0; n < 10; n++)
      bin[n] = 0;
   for (i = 0; i < 1000; i++)
   {
      x = random.gaussian();
      for (n = 0; n < 10; n++)
         if ((x >= 0.5 * n) && (x < 0.5 * (n + 1)))
            bin[n]++;
   }
   for (n = 0; n < 10; n++)
      std::cout << "bin[" << n << "] = " << bin[n] << "\n";
   return 0;
}
