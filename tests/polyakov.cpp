#include "fermiqcd.h"

using namespace MDP;

void polyakov(int nt, int nx, std::string filename)
{
  int L[] = {nt, nx, nx, nx};
  char output[512];
  mdp_lattice lattice(4, L,
                      default_partitioning<1>,
                      torus_topology,
                      0, 1, false);
  mdp_lattice lattice3d(3, L + 1,
                        default_partitioning<1>,
                        torus_topology,
                        0, 1, false);
  gauge_field U(lattice, 3);
  mdp_matrix_field V(lattice3d, 3, 3);
  mdp_real_vector_field q(lattice3d, 2);

  mdp_site x(lattice);
  mdp_site y(lattice3d);

  forallsites(y)
      V(y) = 1;

  U.load(filename);

  for (int t = 0; t < L[0]; t++)
  {
    forallsites(y)
    {
      x.set(t, y(0), y(1), y(2));
      V(y) = V(y) * U(x, 0);
    }
  }

  forallsites(y)
  {
    mdp_complex z = trace(V(y));
    q(y, 0) = real(z);
    q(y, 1) = imag(z);
  }
  /// FIX FILE NAMES depends on T, beta, and K
  snprintf(output, 512, "%s.polyakov.real.vtk", filename.c_str());
  q.save_vtk(output, -1, 0, 0, false);
  snprintf(output, 512, "%s.polyakov.real.vtk", filename.c_str());
  q.save_vtk(output, -1, 1, 0, false);
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  mdp_field_file_header header = get_info("gauge.cold");
  assert(header.ndim == 4);
  assert(header.box[2] == header.box[1]);
  assert(header.box[3] == header.box[1]);
  int nt = header.box[0];
  int nx = header.box[1];
  polyakov(nt, nx, "gauge.cold");
  mdp.close_wormholes();
  return 0;
}
