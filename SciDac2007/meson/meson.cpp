// #define SSE2
#include "../../Libraries/fermiqcd.h"
#include "../dump.h"

using namespace MDP;

int main(int argc, char **argv)
{
  define_base_matrices("FERMILAB");
  mdp.open_wormholes(argc, argv);
  int L[] = {20, 10, 10, 10};
  mdp_lattice lattice(4, L);
  mdp_site x(lattice);
  gauge_field U(lattice, 1);
  fermi_field psi(lattice, 1);
  fermi_field phi(lattice, 1);
  coefficients quark;
  quark["kappa"] = 0.12;
  mdp_lattice space(3, L + 1);
  mdp_real_vector_field s(space, L[0]);
  mdp_site x3(space);
  char filename[128];
  set_cold(U);

  for (int t = 0; t < L[0]; t++)
  {
    forallsites(x3)
    {
      s(x3, t) = 0;
    }
  }

  for (int a = 0; a < 4; a++)
  {
    x.set(10, 5, 5, 5);
    psi = 0;
    psi(x, a, 0) = 1;
    mul_invQ(phi, psi, U, quark, 1e-7);

    for (int t = 0; t < L[0]; t++)
    {
      forallsites(x3)
      {
        x.set(t, x3(0), x3(1), x3(2));
        for (int b = 0; b < 4; b++)
          s(x3, t) += real(phi(x, b, 0) * conj(phi(x, b, 0)));
      }
      snprintf(filename, 128, "meson.%.3i.vtk", t);
      dump(s, t, filename);
    }
  }
  mdp.close_wormholes();
  return 0;
}
