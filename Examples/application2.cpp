// Program: application2.cpp
#include "mdp.h"

using namespace MDP;

float resistance(int x0, int x1, int mu)
{
  return (Pi / 100 * x0);
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int mybox[] = {100, 20};
  mdp_lattice cylinder(2, mybox, default_partitioning<0>,
                       open_cylinder);
  mdp_real_scalar_field u(cylinder);
  mdp_real_vector_field r(cylinder, 2);
  mdp_site x(cylinder);
  mdp_site A(cylinder);
  mdp_site B(cylinder);
  float precision, old_u;
  float c = 0, J = 1, deltaJ, deltaR, Rtot;
  A.set(15, 7);
  B.set(62, 3);

  forallsitesandcopies(x)
  {
    u(x) = 0;
    r(x, 0) = resistance(x(0), x(1), 0);
    r(x, 1) = resistance(x(0), x(1), 1);
  }

  do
  {
    precision = 0;
    forallsites(x)
    {
      old_u = u(x);
      if (x == A)
      {
        c = +J;
        std::cout << ME << "\n";
      }

      if (x == B)
        c = -J;
      deltaJ = c;
      deltaR = 0;

      /* the next two lines take care of the finite cylinder */
      if (x + 0 != x)
      {
        deltaJ += u(x + 0) * r(x, 0);
        deltaR += r(x, 0);
      }

      if (x - 0 != x)
      {
        deltaJ += u(x - 0) * r(x - 0, 0);
        deltaR += r(x - 0, 0);
      }

      deltaJ += u(x + 1) * r(x, 1);
      deltaR += r(x, 1);
      deltaJ += u(x - 1) * r(x - 1, 1);
      deltaR += r(x - 1, 1);
      u(x) = deltaJ / deltaR;
      precision += std::pow(u(x) - old_u, 2);
    }
    u.update();
    mdp.add(precision);
  } while (sqrt(precision) > 0.00001);

  Rtot = 0;
  if (A.is_in())
    Rtot += u(A) / J;
  if (B.is_in())
    Rtot -= u(B) / J;

  mdp.add(Rtot);
  mdp << "R_A-R_B=" << Rtot << "\n";

  mdp.close_wormholes();
  return 0;
}
