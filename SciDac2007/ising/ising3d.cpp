#include "../../Libraries/mdp.h"
// #include "../mdp_all.h"
#include "../dump.h"

#define X 0
#define Y 1
#define Z 2

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  constexpr Box L = {20, 20, 20};
  mdp_lattice cube(L);
  mdp_real_scalar_field spin(cube);
  mdp_site x(cube);
  mdp_int dE = 0, H = L.volume(), dH = 0;
  mdp_real kappa = 0.40;
  if (argc > 1)
    kappa = atof(argv[1]); // try 0.5 or 0.25

  forallsites(x)
  {
    spin(x) = (x(0) > L[0] / 4 && x(0) <= 3 * L[0] / 4) ? (+1) : (-1);
  }

  for (mdp_suint i = 0; i < 100; i++)
  {
    dH = 0;
    for (mdp_parity parity : {EVEN, ODD})
    {
      forallsitesofparity(x, parity)
      {
        dE = 2 * spin(x) * (spin(x - X) + spin(x + X) + spin(x - Y) + spin(x + Y) + spin(x - Z) + spin(x + Z));
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
    // spin.save_vtk(std::format("ising3d{:03d}.vtk", i), -1, -1, 0, true);
    dump(spin, 0, std::format("ising3d_{:03d}.vtk", i), true);
  }

  // spin.save_vtk("ising3d.vtk", -1, -1, 0, true);
  dump(spin);
  mdp.close_wormholes();
  return 0;
}
