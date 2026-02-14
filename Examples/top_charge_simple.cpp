#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv); // mpirun

  constexpr Box L = {4, 10, 10, 10};
  mdp_lattice lattice(L);
  gauge_field U(lattice, 3);
  coefficients gauge;
  gauge["beta"] = 5.0;
  // mdp_real_scalar_field Q(lattice);

  set_cold(U);

  for (int k = 0; k < 10; k++)
  {
    std::cout << k << std::endl;
    WilsonGaugeAction::heatbath(U, gauge, 1);
  }
  topological_charge_vtk(U, "top_charge_simple.vtk");

  mdp.close_wormholes();
  return 0;
}
