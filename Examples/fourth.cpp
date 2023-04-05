#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv); // mpirun

  int L[] = {6, 6, 6, 6};
  mdp_lattice lattice(4, L);
  gauge_field U(lattice, 3);
  mdp_site x(lattice);
  coefficients gauge;
  gauge["beta"] = 6.0;

  forallsites(x)
  {
    for (int mu = 0; mu < 4; mu++)
      U(x, mu) = lattice.random(x).SU(3);
  }

  U.update();
  for (mdp_uint k = 0; k < 1000; k++)
  {
    WilsonGaugeAction::heatbath(U, gauge, 1);
    U.save("somename");
  }

  mdp << "hello world\n";

  mdp.close_wormholes();
  return 0;
}
