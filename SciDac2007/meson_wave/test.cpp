#include "../../Libraries/fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  std::cout << Gamma[0] << std::endl;
  std::cout << Gamma5 << std::endl;
  std::cout << (1 + Gamma[0]) / 2 * Gamma5 << std::endl;
  return 0;
}
