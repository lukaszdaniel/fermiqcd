// Program: example12.cpp
#include "mdp.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  constexpr Box mybox = {10, 10};
  mdp_lattice mylattice(mybox);
  mdp_site x(mylattice);
  forallsites(x)
  {
    std::cout << "Site=(" << x(0) << "," << x(1) << ") is stored by " << ME << "\n";
  }
  mdp.close_wormholes();
}
