#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  constexpr mdp_suint nc = 4;
  constexpr mdp_suint Nconfig = 5;
  constexpr Box mybox = {2, 2, 2, 2, 2};
  mdp_lattice mylattice(mybox);
  gauge_field U(mylattice, nc);
  coefficients coeff;
  coeff["beta"] = 2.0;

  set_hot(U);

  for (mdp_suint config = 0; config < Nconfig; config++)
  {
    WilsonGaugeAction::heatbath(U, coeff, 1);
    mdp << "config=" << config
        << "plaquette="
        << average_plaquette(U) << "\n";
  }
  mdp.close_wormholes();
  return 0;
}
