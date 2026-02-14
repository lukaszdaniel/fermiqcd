#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int nc = 3;
  int nt, nx;
  sscanf(argv[1], "%ix%i", &nt, &nx);
  const Box box = {nt, nx, nx, nx};
  mdp_lattice lattice(box, default_partitioning0, torus_topology, 0, 2, false);
  gauge_field U(lattice, nc);
  U.load(argv[2]);
  mdp_site x(lattice);
  x.set(0, 0, 0, 0);
  std::cout << U(x, 0) << std::endl;
  forallsites(x)
  {
    for (int mu = 0; mu < 4; mu++)
      if (max(U(x, mu) - 1) > 0)
      {
        std::cout << x << mu << "\n"
             << U(x, mu) << std::endl;
        break;
      }
  }
  mdp << "plaquette:" << average_plaquette(U) << "\n";
  mdp.close_wormholes();
  return 0;
}
