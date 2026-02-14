#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  constexpr Box L = {10, 10, 8, 6};
  mdp_lattice lattice(L);
  mdp_site x(lattice);
  mdp_real_vector_field s(lattice, 2);
  forallsites(x)
  {
    s(x, 0) = std::sqrt(x(0) * x(0) + x(1) * x(1) + x(2) * x(2) + x(3) * x(3));
    s(x, 1) = std::sqrt((10 - x(0)) * (10 - x(0)) + x(1) * x(1) + x(2) * x(2) + x(3) * x(3));
  }
  s.save_vtk("test_save_vtk", 0);
  s.save_vtk("test_save_vtk");

  mdp.close_wormholes();
  return 0;
}
