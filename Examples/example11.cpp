// Program: example11.cpp
#include "mdp.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mpi.open_wormholes(argc, argv);
  int mybox[] = {8, 8};
  mdp_lattice mylattice(2, mybox);
  mpi.close_wormholes();
}
