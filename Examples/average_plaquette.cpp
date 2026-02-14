#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  constexpr Box L = {32, 24, 24, 24};
  mdp_lattice lattice(L, default_partitioning0, torus_topology, 1, 0, false);
  gauge_field U(lattice, 3);
  U.load(argv[1]);
  mdp_site x(lattice), y(lattice);

  x.set(0, 0, 0, 0);
  y.set(31, 23, 23, 23);
  for (int mu = 0; mu < 4; mu++)
    U(y, mu) = U(x, mu);
  std::cout << average_plaquette(U) << std::endl;

  mdp.close_wormholes();
  return 0;
}
