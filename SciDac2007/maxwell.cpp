#include "mdp.h"
#include "mdp_all.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  constexpr Box L = {20, 20, 20};
  mdp_lattice cube(L);
  mdp_site x(cube);
  mdp_complex_vector_field E(cube, 3);
  mdp_complex_vector_field B(cube, 3);
  mdp_complex_scalar_field q(cube);
  mdp_complex_vector_field j(cube, 3);

  forallsites(x)
  {
    E(x) = 0;
    B(x) = 0;
  }

  mdp.close_wormholes();
  return 0;
}
