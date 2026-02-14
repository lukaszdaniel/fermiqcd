#include "fermiqcd.h"

using namespace MDP;

void compute_plaquette(mdp_int nt, mdp_int nx, const std::string &filename)
{
  const Box L = {nt, nx, nx, nx};
  mdp_int nc = 3;
  mdp_lattice lattice(L,
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
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  const std::string saved_field = "gauge.cold";
  mdp_field_file_header header = get_info(saved_field);
  assert(header.ndim == 4);
  assert(header.box[2] == header.box[1]);
  assert(header.box[3] == header.box[1]);
  mdp_int nt = header.box[0];
  mdp_int nx = header.box[1];
  compute_plaquette(nt, nx, saved_field);
  mdp.close_wormholes();
  return 0;
}
