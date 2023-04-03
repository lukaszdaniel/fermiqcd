#include "fermiqcd.h"

void compute_energy(int nt, int nx, std::string filename)
{
  int L[] = {nt, nx, nx, nx};
  int nc = 3;
  char output[512];
  mdp_lattice lattice(4, L,
                      default_partitioning<1>,
                      torus_topology,
                      0, 1, false);
  gauge_field U(lattice, nc);
  mdp_field<mdp_real> q(lattice);
  mdp_site x(lattice);

  U.load(filename);
  compute_em_field(U);
  forallsites(x) if (x(0) == 0)
  {
    q(x) = 0;
    for (int mu = 0; mu < 4; mu++)
      for (int nu = mu + 1; nu < 4; nu++)
        q(x) += pow(abs(trace(U.em(x, mu, nu))), 2);
  }

  /// FIX FILE NAMES depends on T, beta, and K
  snprintf(output, 512, "%s.e2b2.vtk", filename.c_str());
  q.save_vtk(output, 0);
}

int main(int argc, char **argv)
{
  mpi.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  mdp_field_file_header header = get_info(argv[1]);
  assert(header.ndim == 4);
  assert(header.box[2] == header.box[1]);
  assert(header.box[3] == header.box[1]);
  int nt = header.box[0];
  int nx = header.box[1];
  compute_energy(nt, nx, argv[1]);
  mpi.close_wormholes();
  return 0;
}
