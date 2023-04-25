// Program: example13.cpp
#include "mdp.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int mybox[] = {10, 10};
  mdp_lattice mylattice(2, mybox);
  mdp_site x(mylattice);
  int mu = 0;
  if (ME == 0)
  {
    x.set(0, 0);
    do
    {
      std::cout << "x=(" << x(0) << "," << x(1) << ")\n";
      if (x.is_in_boundary())
        error("I found the boundary");
      x = x + mu;
    } while (true);
  }
  mdp.close_wormholes();
}
