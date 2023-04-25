#include "fermiqcd.h"

using namespace MDP;

void compute_plaquette(mdp_int nt, mdp_int nx, const std::string &filename)
{
  mdp_int L[] = {nt, nx, nx, nx};
  mdp_int nc = 3;
  char output[512];
  mdp_lattice lattice(4, L,
                      default_partitioning<1>,
                      torus_topology,
                      0, 1, false);
  gauge_field U(lattice, nc);
  mdp_real_scalar_field q(lattice);
  mdp_site x(lattice);

  U.load(filename);
  forallsites(x)
  {
    if (x(0) == 0)
    {
      q(x) = 0;
      for (int mu = 0; mu < 4; mu++)
        for (int nu = mu + 1; nu < 4; nu++)
          q(x) += real(trace(plaquette(U, x, mu, nu)));
    }
  }

  /// FIX FILE NAMES depends on T, beta, and K
  snprintf(output, 512, "%s.plaquette.vtk", filename.c_str());
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
  mdp_int nt = header.box[0];
  mdp_int nx = header.box[1];
  compute_plaquette(nt, nx, argv[1]);
  mdp.close_wormholes();
  return 0;
}
