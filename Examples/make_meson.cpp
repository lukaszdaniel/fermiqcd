// #define PARALLEL
#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  // Setting parameters here:
  coefficients gauge;
  gauge["beta"] = 6.0;

  coefficients light_quark;
  light_quark["kappa_s"] = 0.1;
  light_quark["kappa_t"] = 0.1;
  light_quark["r_s"] = 1;
  light_quark["r_t"] = 1;
  light_quark["c_{sw}"] = 0;
  light_quark["c_E"] = 1;
  light_quark["c_B"] = 1;

  coefficients heavy_quark;
  heavy_quark["kappa_s"] = 0.1;
  heavy_quark["kappa_t"] = 0.1;
  heavy_quark["r_s"] = 1;
  heavy_quark["r_t"] = 1;
  heavy_quark["c_{sw}"] = 0;
  heavy_quark["c_E"] = 1;

  double absolute_precision = 1e-23;
  double relative_precision = 1e-12;
#ifdef SSE2
  default_fermi_action = FermiCloverActionSSE2::mul_Q;
#else
  default_fermi_action = FermiCloverActionFast::mul_Q;
#endif
  default_fermi_inverter = MinimumResidueInverter<fermi_field, gauge_field>;

  int nt = 12;
  int nc = 3;
  // double alpha=0.0;

  int mybox[] = {nt, 4, 4, 4};

  // creating the fields
  mdp_lattice mylattice(4, mybox);
  gauge_field U(mylattice, nc);
  fermi_field psi(mylattice, nc);
  fermi_field phi(mylattice, nc);
  fermi_field chi(mylattice, nc);
  mdp_site x(mylattice);
  mdp_array<mdp_complex, 1> C2(nt);

  // creating one gauge configuration

  set_cold(U);
  // WilsonGaugeAction::heatbath(U,gauge,20); // 2 iterations

#ifdef USE_S

  fermi_propagator S(mylattice, nc);

  fermi_propagator::generate(S, U, light_quark);

  forallsites(x)
  {
    for (a = 0; a < 4; a++)
      for (b = 0; b < 4; b++)
        for (i = 0; i < 3; i++)
          for (j = 0; j < 3; j++)
          {
            C2(x(0)) += S(x, a, b, i, j) * conj(S(x, a, b, i, j));
          }
  }

#else

  for (int t = 0; t < nt; t++)
    C2(t) = 0;

  for (int a = 0; a < 4; a++)
    for (int i = 0; i < nc; i++)
    {
      psi = 0;
      x.set(0, 0, 0, 0);
      psi(x, a, i) = 1;
      psi.update();

#if 0
      // ******* smearing here *******************
      for (n = 0; n < 10; n++)
      {
        forallsites(x)
        {
          chi(x) = psi(x) + alpha * (psi(x + 0) + psi(x - 0) +
                                     psi(x + 1) + psi(x - 1) +
                                     psi(x + 2) + psi(x - 2) +
                                     psi(x + 3) + psi(x - 3));
        }
        forallsites(x)
        {
          psi(x) = chi(x);
        }
      }
      // **********end of smearing ***************** 
#endif

          // make light prop phi
          mul_invQ(phi, psi, U, light_quark, absolute_precision, relative_precision);
      // make heavy prop chi
      mul_invQ(chi, psi, U, heavy_quark, absolute_precision, relative_precision);

      forallsites(x)
      {
        for (int b = 0; b < 4; b++)
          for (int j = 0; j < nc; j++)
            C2(x(0)) += phi(x, b, j) * conj(chi(x, b, j));
      }
    }

#endif

  mdp.add(C2.address(), nt);
  for (int t = 0; t < nt; t++)
    printf("%i\t%e\t%e\n", t, real(C2(t)), imag(C2(t)));

  mdp.close_wormholes();
};
