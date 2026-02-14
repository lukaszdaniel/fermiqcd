#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int nc = 3;
  constexpr Box box = {8, 4, 4, 4};
  mdp_lattice lattice(box);

  gauge_field U(lattice, nc);
  set_hot(U);
  U.save("gauge.hot");
  mdp.close_wormholes();
  return 0;
}
