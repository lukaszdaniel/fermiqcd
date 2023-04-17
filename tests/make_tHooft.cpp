#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  define_base_matrices("FERMILAB");
  mdp.open_wormholes(argc, argv);
  int nc = 2;
  int L[] = {10, 40, 10, 10};
  mdp_lattice lattice(4, L, default_partitioning0, torus_topology, 0, 2, false);
  gauge_field U(lattice, nc);
  char filename[128];
  mdp_site x(lattice);
  set_cold(U);
  std::vector<mdp_real> p(4);

  for (int mu = 0; mu < 4; mu++)
    p[mu] = L[mu] / 2;

  int i = 0, j = 0;

  Instanton4D A(nc, i, j, 1.0, 3, p);

  forallsites(x)
  {
    for (int mu = 0; mu < U.ndim(); mu++)
      U(x, mu) = exp(-I * A(x, mu));
  }

  snprintf(filename, 128, "%s.vtk", argv[0]);
  mdp << "top=" << topological_charge_vtk(U, filename) << "\n";

  mdp.close_wormholes();
  return 0;
}
