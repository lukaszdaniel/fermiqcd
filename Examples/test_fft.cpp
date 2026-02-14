#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int nc = 3;
  constexpr Box L = {8, 4, 4, 4};
  mdp_lattice lattice(L);
  fermi_field psi(lattice, nc);
  fermi_field phi(lattice, nc);
  fermi_field chi(lattice, nc);
  mdp_site x(lattice);
  set_random(psi);

  fermi_field_fft(phi, psi, 1, true);
  fermi_field_fft(chi, phi, -1, true);

  forallsites(x)
  {
    if (max(chi(x) - psi(x)) > 0.01)
    {
      mdp << psi(x) << "\n";
      mdp << chi(x) << "\n";
    }
  }

  mdp.close_wormholes();
  return 0;
}
