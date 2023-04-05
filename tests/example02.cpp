// Program: example02.cpp
#include "mdp.h"

using namespace MDP;

int main()
{
   mdp_matrix a(2, 2);
   a(0, 0) = 2;
   a(0, 1) = 3;
   a(1, 0) = 4;
   a(1, 1) = 5 * I;
   std::cout << inv(exp(a)) << "\n";
   return 0;
}
