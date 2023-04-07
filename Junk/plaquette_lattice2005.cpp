#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv); // START
  int n = 4;
  constexpr mdp_suint nconfig = 100;
  int L[] = {8, 8, 8, 8, 8};
  mdp_lattice lattice(5, L); // declare lattice
  gauge_field U(lattice, n); // declare fields
  coefficients gauge;
  gauge["beta"] = 6.0; // set parameters
  set_cold(U);
  for (mdp_suint k = 0; k < nconfig; k++)
  {
    WilsonGaugeAction::heatbath(U, gauge, 10);        // do heathbath
    U.save(std::string("gauge") + std::to_string(k)); // save gauge config
    mdp << average_plaquette(U) << "\n";              // print plaquette
  }
  mdp.close_wormholes(); // STOP
  return 0;
}
