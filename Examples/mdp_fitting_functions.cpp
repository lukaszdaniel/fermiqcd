#include "mdp.h"

using namespace MDP;

// Example of a program that uses the Leverberger-Marquardt

// fitting function
mdp_real f(mdp_real x, mdp_real *a, [[maybe_unused]] mdp_int ma, [[maybe_unused]] void *junk)
{
  return a[0] * (exp(-a[1] * x) + exp(-a[1] * (10.0 - x))) * sin(a[2] * x);
}

int main()
{
  mdp_real x[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  mdp_measure y[10];
  constexpr mdp_int ma = 3;
  mdp_real a[3] = {1000.0, 0.1, 0.1};

  mdp_matrix covar(ma, ma);
  for (mdp_suint i = 0; i < 10; i++)
  {
    x[i] = i;
    y[i].set(f(x[i], a, ma, 0), 0.1);
  }

  a[0] = 1;
  a[1] = 1.0, a[2] = 0.2;

  BayesianLevenbergMarquardt(x, y, 0, 10, a, ma, covar, f);

  for (mdp_suint i = 0; i < 10; i++)
    printf("%f %f (%f) %f\n",
           x[i], y[i].getmean(), y[i].getmerr(), f(x[i], a, ma, 0));

  printf("%f %f %f\n", a[0], a[1], a[2]);
  std::cout << covar;
}