// Program: example07.cpp
#include "mdp.h"

using namespace MDP;

constexpr int n = 100;

float f1(float *x, void *a) { return x[0] / x[1]; }

float f2(float *x, void *a) { return x[0] * x[1]; }

int main()
{
   mdp_matrix A;
   mdp_jackboot jb(n, 2);

   for (int i = 0; i < n; i++)
   {
      A = Random.SU(6) + Random.gaussian();
      jb(i, 0) = real(det(inv(A)));
      jb(i, 1) = real(det(A));
   }

   jb.f = f1;
   std::cout << "Result x[0]/x[1] = " << jb.mean() << "\n";
   std::cout << "Jackknife error  = " << jb.j_err() << "\n";
   std::cout << "Bootstrap error  = " << jb.b_err(100) << "\n";

   jb.f = f2;
   std::cout << "Result x[0]*x[1] = " << jb.mean() << "\n";
   std::cout << "Jackknife error  = " << jb.j_err() << "\n";
   std::cout << "Bootstrap error  = " << jb.b_err(100) << "\n";

   return 0;
}
