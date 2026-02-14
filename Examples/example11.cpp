// Program: example11.cpp
#include "mdp.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  constexpr Box mybox = {8, 8};
  mdp_lattice mylattice(mybox);
  mdp.close_wormholes();
}
