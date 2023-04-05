// #define SSE2
#define OSX
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
  staggered_field psi(lattice, 1);
  staggered_field phi(lattice, 1);
  coefficients quark;
  quark["mass"] = 0.3;

  set_cold(U);
  x.set(0, 5, 5, 5);
  psi(x, 0) = 1;

  mul_invQ(phi, psi, U, quark, 1e-7);

  mdp_lattice space(3, L + 1);
  mdp_field<float> s(space, 3);
  mdp_site x3(space);
  char filename[128];
  for (int t = 0; t < L[0] - 1; t++)
  {
    forallsites(x3)
    {
      x.set((L[0] / 2 + 1 + t) % L[0], x3(0), x3(1), x3(2));
      s(x3, 0) = real(phi(x, 0));
      s(x3, 1) = imag(phi(x, 0));
      s(x3, 2) = abs(phi(x, 0));
    }
    snprintf(filename, 128, "staggered.real.%.3i.vtk", t);
    dump(s, 0, filename);
    snprintf(filename, 128, "staggered.imag.%.3i.vtk", t);
    dump(s, 1, filename);
    snprintf(filename, 128, "staggered.abs.%.3i.vtk", t);
    dump(s, 2, filename);

    forallsites(x)
    {
      if (imag(phi(x, 0)) != 0)
        std::cout << x << std::endl;
    }
  }

  mdp.close_wormholes();
  return 0;
}
