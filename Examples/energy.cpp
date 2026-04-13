#include "fermiqcd.h"

using namespace MDP;

void compute_energy(mdp_uint nt, mdp_uint nx, const std::string &filename)
{
  const Box L = {nt, nx, nx, nx};
  mdp_suint nc = 3;
  std::string output;
  mdp_lattice lattice(L,
                      default_partitioning<1>,
                      torus_topology,
                      0, 1, false);
  gauge_field U(lattice, nc);
  mdp_real_scalar_field q(lattice);
  mdp_site x(lattice);

  U.load(filename);
  compute_em_field(U);
  forallsites(x)
  {
    if (x(0) == 0)
    {
      q(x) = 0;
      for (mdp_suint mu = 0; mu < 4; mu++)
        for (mdp_suint nu = mu + 1; nu < 4; nu++)
          q(x) += std::pow(abs(trace(U.em(x, mu, nu))), 2);
    }
  }

  /// FIX FILE NAMES depends on T, beta, and K
  output = std::format("{}.e2b2.vtk", filename);
  q.save_vtk(output, 0);
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
  compute_energy(nt, nx, argv[1]);
  mdp.close_wormholes();
  return 0;
}
