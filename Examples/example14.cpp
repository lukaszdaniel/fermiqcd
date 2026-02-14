// Program: example14.cpp
#include "mdp.h"

using namespace MDP;

struct mystruct
{
  float value; /* or any other structure */
};

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  constexpr Box mybox = {10, 10};
  mdp_lattice mylattice(mybox);
  mdp_field<mystruct> myfield(mylattice);
  mdp_site x(mylattice);
  forallsites(x)
  {
    myfield(x).value = 0;
  }
  myfield.update();
  mdp.close_wormholes();
}
