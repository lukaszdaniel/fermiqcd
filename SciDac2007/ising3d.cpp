#include "mdp.h"

using namespace MDP;

#define X 0
#define Y 1
#define Z 2

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int L[] = {20, 20, 20};
  mdp_lattice cube(3, L);
  mdp_real_scalar_field spin(cube);
  mdp_site point(cube);
  int dE = 0, H = cube.size(), dH = 0;
  float kappa = 0.40;
  if (argc > 1)
    kappa = atof(argv[1]); // try 0.5 or 0.25

  forallsites(point)
  {
    spin(point) = (point(0) > L[0] / 4 && point(0) <= 3 * L[0] / 4) ? (+1) : (-1);
  }

  for (int i = 0; i < 100; i++)
  {
    dH = 0;
    for (int parity = 0; parity < 2; parity++)
    {
      forallsitesofparity(point, parity)
      {
        dE = 2 * spin(point) * (spin(point - X) + spin(point + X) + spin(point - Y) + spin(point + Y) + spin(point - Z) + spin(point + Z));
        if (exp(-kappa * dE) > mdp_random.plain())
        {
          spin(point) *= -1;
          dH = dH + 2 * spin(point);
        }
      }
      spin.update();
    }
    mdp.add(dH);
    H = H + dH;
    mdp << "magnetization=" << H << "\n";
  }
  spin.save_vtk("ising3d.vtk", -1, -1, 0, true);
  mdp.close_wormholes();
  return 0;
}
