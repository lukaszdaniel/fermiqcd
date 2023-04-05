#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv); // mpirun

  int L[] = {10, 10, 10, 10};
  mdp_lattice lattice(4, L);
  gauge_field U(lattice, 3);
  mdp_site x(lattice);
  coefficients gauge;
  gauge["beta"] = 6.0;

  int k;
  // int mu=0,nu=1;
  mdp_complex s = 0;

  forallsites(x)
  {
    for (int mu = 0; mu < 4; mu++)
      U(x, mu) = lattice.random(x).SU(3);
  }

  int d[4][2];

  U.update();
  for (k = 0; k < 100; k++)
  {
    WilsonGaugeAction::heatbath(U, gauge, 1);

    forallsites(x)
    {
      for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
          if (mu != nu)
          {
            // d={{1,mu},{1,nu},{-1,mu},{-1,nu}};
            d[0][0] = 1;
            d[0][1] = mu;
            d[1][0] = 1;
            d[1][1] = nu;
            d[2][0] = -1;
            d[2][1] = mu;
            d[3][0] = -1;
            d[3][1] = nu;
            // s+=trace(U(x,mu)*U(x+mu,nu)*hermitian(U(x+nu,mu))*hermitian(U(x,nu)));
            s += trace(build_path(U, x, 4, d)); // d={{1,mu},{1,nu},{-1,mu},{-1,nu}}
          }
    }
  }

  mdp << s / 100 / (6 * 6 * 6 * 6) / 12 << "\n";

  mdp.close_wormholes();
  return 0;
}
