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
  mdp << "success in allocating vector\n";
  std::string filename2;
  U.load(filename);
  // U.switch_endianess_4bytes();
  ApeSmearing::smear(U, 0.7, 20, 10);
  filename2 = std::format("{}.cooled20", filename);
  U.save(filename2);
  filename2 = std::format("{}.topological_charge_20.vtk", filename);
  float tc = topological_charge_vtk(U, filename2, 0);
  mdp << "topological_charge=" << tc << "\n";
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  mdp_field_file_header header = get_info(argv[1]);
  assert(header.ndim == 4);
  assert(header.box[2] == header.box[1]);
  assert(header.box[3] == header.box[1]);
  mdp_uint nt = header.box[0];
  mdp_uint nx = header.box[1];
  test_gauge(nt, nx, argv[1]);
  mdp.close_wormholes();
  return 0;
}
