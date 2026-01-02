#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int nc = 3;
  int box[] = {8, 4, 4, 4};
  mdp_lattice lattice(4, box);
  mdp_site x(lattice);

  gauge_field U(lattice, nc);
  set_hot(U);
  U.save("gauge.hot");

  gauge_field V(lattice, nc);
  V.load("gauge.hot");

  x.set(3, 2, 1, 0);
  // std::cout << "U(x) = " << U(x, 1) << std::endl;
  // std::cout << "V(y) = " << V(x, 1) << std::endl;
  std::cout << "Matrices are equal = " << std::boolalpha << (U(x, 1) == V(x, 1)) << std::endl;

  mdp.close_wormholes();
  return 0;
}
