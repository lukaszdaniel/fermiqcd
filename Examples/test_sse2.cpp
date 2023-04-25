// #define BLOCKSITE 4
// #define TWISTED_BOUNDARY
#include "fermiqcd.h"

using namespace MDP;

void test_fermi()
{
  int box[] = {4, 4, 4, 4};
  int nc = 3;

  mdp << "\n\nTEST FERMI FIELDS\n\n";

  mdp_lattice lattice(4, box);
  gauge_field U(lattice, nc);
  fermi_field psi(lattice, nc);
  fermi_field chi1(lattice, nc);
  fermi_field chi2(lattice, nc);
  coefficients coeff;

  coeff["kappa_s"] = 0.1;
  coeff["kappa_t"] = 0.1;
  coeff["c_{sw}"] = 1.00;

  set_hot(U);
  compute_em_field(U);
  set_random(psi);

  default_fermi_action = FermiCloverActionSlow::mul_Q;
  mul_Q(chi1, psi, U, coeff);

  default_fermi_action = FermiCloverActionFast::mul_Q;
  mul_Q(chi2, psi, U, coeff);

  mdp << "\n\nCheching that CloverActionFast and CloverActionSlow agree\n\n";
  check_differences(chi1, chi2);

#ifdef SSE2
  default_fermi_inverter = MinimumResidueInverter<fermi_field, gauge_field>;
  default_fermi_action = FermiCloverActionSSE2::mul_Q;
  mul_Q(chi2, psi, U, coeff);

  mdp << "\n\nCheching that CloverActionSlow and CloverActionSSE2 agree\n\n";
  check_differences(chi1, chi2);
#endif
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  test_fermi();
  mdp.close_wormholes();
  return 0;
}
