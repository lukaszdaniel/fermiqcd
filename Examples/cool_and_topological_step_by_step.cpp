// #define BLOCKSITE 4
// #define TWISTED_BOUNDARY
#include "fermiqcd.h"

using namespace MDP;

void test_gauge(mdp_uint nt, mdp_uint nx, const char *filename)
{
  const Box box = {nt, nx, nx, nx};
  mdp_suint nc = 3;
  mdp_lattice lattice(box,
                      default_partitioning0,
                      torus_topology,
                      0, 1, false);
  gauge_field U(lattice, nc);
  std::string filename2;
  U.load(filename);
  // U.switch_endianess_4bytes();
  for (mdp_suint k = 0; k <= 20; k += 5)
  {
    filename2 = std::format("{}.topological_charge_{:03d}.vtk", filename, k);
    mdp_real tc = topological_charge_vtk(U, filename2, 0);
    mdp << "topological_charge=" << tc << "\n";
    ApeSmearing::smear(U, 0.7, 5, 10);
  }
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  mdp_uint nt, nx;
  sscanf(argv[1], "%ux%u", &nt, &nx);
  test_gauge(nt, nx, argv[2]);
  mdp.close_wormholes();
  return 0;
}
