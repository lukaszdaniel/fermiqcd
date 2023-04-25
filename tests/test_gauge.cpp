#include "fermiqcd.h"

using namespace MDP;

void test_gauge()
{
  int box[] = {4, 4, 4, 4};
  int nc = 3;
  mdp_lattice lattice(4, box);
  gauge_field U(lattice, nc);
  set_hot(U);
  std::cout << average_plaquette(U, 0, 1) << std::endl;
  int path[][2] = {{+1, 0}, {+1, 1}, {-1, 0}, {-1, 1}};
  std::cout << real(average_path(U, 4, path)) << std::endl;
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  test_gauge();
  mdp.close_wormholes();
  return 0;
}
