#include "mdp.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  constexpr Box L = {100};
  mdp_lattice line(L);
  mdp_int_scalar_field spin(line);
  mdp_site x(line);
  mdp_int dE = 0, H = L.volume(), dH = 0;
  mdp_real kappa = 2.0;
  forallsites(x)
  {
    spin(x) = +1;
  }

  for (mdp_suint i = 0; i < 100; i++)
  {
    dH = 0;
    for (mdp_parity parity : {EVEN, ODD})
    {
      forallsitesofparity(x, parity)
      {
        dE = 2 * spin(x) * (spin(x - 0) + spin(x + 0));
        if (exp(-kappa * dE) > mdp_global_random.plain())
        {
          spin(x) *= -1;
          dH = dH + 2 * spin(x);
        }
      }
      spin.update();
    }
    mdp.add(dH);
    H = H + dH;
    mdp << "magnetization=" << H << "\n";
  }
  mdp.close_wormholes();
  return 0;
}
