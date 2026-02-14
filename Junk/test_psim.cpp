#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  constexpr Box box = {4, 4, 4, 4};
  int nc = 3;
  mdp_lattice lattice(box);
  gauge_field U(lattice, nc);
  coefficients coeff;
  coeff["beta"] = 6.0;
  set_hot(U);
  mdp << "average plaquette = " << average_plaquette(U) << "\n";
  mdp.close_wormholes();

  return 0;
}
