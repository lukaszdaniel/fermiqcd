// #define SSE2
#include "../../Libraries/fermiqcd.h"
#include "../dump.h"

using namespace MDP;

int main(int argc, char **argv)
{
  define_base_matrices("FERMILAB");
  mdp.open_wormholes(argc, argv);
  constexpr Box L = {20, 10, 10, 10};
  constexpr Box L_space = {10, 10, 10};
  mdp_lattice lattice(L);
  mdp_site x(lattice);
  gauge_field U(lattice, 1);
  fermi_field psi(lattice, 1);
  fermi_field phi(lattice, 1);
  coefficients quark;
  quark["kappa"] = 0.115;

  set_cold(U);
  psi = 0;
  x.set(0, 5, 5, 5);
  psi(x, 0, 0) = 1;

  mul_invQ(phi, psi, U, quark, 1e-8);

  mdp_lattice space(L_space);
  mdp_real_vector_field s(space, 9);
  mdp_site x3(space);
  char filename[128];
  for (int t = 0; t < L[0] - 1; t++)
  {
    forallsites(x3)
    {
      x.set(t % L[0], x3(0), x3(1), x3(2));
      s(x3, 0) = real(phi(x, 0, 0));
      s(x3, 1) = imag(phi(x, 0, 0));
      s(x3, 2) = real(phi(x, 1, 0));
      s(x3, 3) = imag(phi(x, 1, 0));
      s(x3, 4) = real(phi(x, 2, 0));
      s(x3, 5) = imag(phi(x, 2, 0));
      s(x3, 6) = real(phi(x, 3, 0));
      s(x3, 7) = imag(phi(x, 3, 0));
      s(x3, 8) = std::sqrt(real(phi(x, 0, 0) * conj(phi(x, 0, 0)) + phi(x, 1, 0) * conj(phi(x, 1, 0)) + phi(x, 2, 0) * conj(phi(x, 2, 0)) + phi(x, 3, 0) * conj(phi(x, 3, 0))));
    }
    snprintf(filename, 128, "wilson.s0.real.%.3i.vtk", t);
    dump(s, 0, filename);
    snprintf(filename, 128, "wilson.s0.imag%.3i.vtk", t);
    dump(s, 1, filename);
    snprintf(filename, 128, "wilson.s1.real.%.3i.vtk", t);
    dump(s, 2, filename);
    snprintf(filename, 128, "wilson.s1.imag.%.3i.vtk", t);
    dump(s, 3, filename);
    snprintf(filename, 128, "wilson.s2.real.%.3i.vtk", t);
    dump(s, 4, filename);
    snprintf(filename, 128, "wilson.s2.imag.%.3i.vtk", t);
    dump(s, 5, filename);
    snprintf(filename, 128, "wilson.s3.real.%.3i.vtk", t);
    dump(s, 6, filename);
    snprintf(filename, 128, "wilson.s3.imag.%.3i.vtk", t);
    dump(s, 7, filename);
    snprintf(filename, 128, "wilson.squared.%.3i.vtk", t);
    dump(s, 8, filename);
  };

  mdp.close_wormholes();
  return 0;
}
