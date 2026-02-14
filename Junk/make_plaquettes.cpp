#include "fermiqcd.h"

using namespace MDP;

mdp_real average_plaquette1(gauge_field &U, int mu, int nu)
{
  mdp_lattice &lattice = U.lattice();
  mdp_site x(lattice);
  mdp_real sum = 0.0;
  forallsites(x)
      sum += real(trace(U(x, mu) * U(x + mu, nu) *
                        hermitian(U(x, nu) * U(x + nu, mu))));
  mdp.add(sum);
  return sum / (lattice.size() * U.nc());
}

mdp_real average_plaquette2(gauge_field &U, int mu, int nu)
{
  mdp_lattice &lattice = U.lattice();
  mdp_site x(lattice);
  mdp_complex sum = 0.0;
  Path path = {{+1, mu}, {+1, nu}, {-1, mu}, {-1, nu}};
  sum = average_path(U, path);
  return real(sum);
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  int nc = 3;
  constexpr Box box = {8, 4, 4, 4};
  mdp_lattice lattice(box);
  gauge_field U(lattice, nc);
  int mu = 0, nu = 1;
  set_hot(U);
  mdp << "average_plaquette  : " << average_plaquette(U, mu, nu) << "\n";
  mdp << "average_plaquette1 : " << average_plaquette1(U, mu, nu) << "\n";
  mdp << "average_plaquette2 : " << average_plaquette2(U, mu, nu) << "\n";
  mdp.close_wormholes();
  return 0;
}
