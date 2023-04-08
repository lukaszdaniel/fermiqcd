// Program: application3.cpp
#include "mdp.h"

using namespace MDP;

// /////////////////////////////////////////////////////
class scalar_field : public mdp_real_scalar_field
{
public:
  const int ndim;

  scalar_field(mdp_lattice &a) : mdp_real_scalar_field(a), ndim(a.n_dimensions())
  {
  }
};

// /////////////////////////////////////////////////////
class ising_field : public mdp_int_scalar_field
{
public:
  const int ndim;
  mdp_real beta;
  mdp_real kappa;
  mdp_real magnetic_field;

  ising_field(mdp_lattice &a) : mdp_int_scalar_field(a), ndim(a.n_dimensions())
  {
  }
};

void set_cold(ising_field &S)
{
  mdp_site x(S.lattice());

  forallsites(x)
  {
    S(x) = 1;
  }
}

void montecarlo_multihit(ising_field &S, scalar_field &H, int n_iter = 1, int n_hits = 1)
{
  mdp_int parity;
  mdp_int new_spin;
  mdp_real delta_action;
  mdp_site x(S.lattice());

  for (int iter = 0; iter < n_iter; iter++)
  {
    for (parity = 0; parity <= 1; parity++)
    {
      forallsitesofparity(x, parity)
      {
        for (int hit = 0; hit < n_hits; hit++)
        {
          delta_action = S.kappa * (H(x) + S.magnetic_field);
          for (int mu = 0; mu < S.ndim; mu++)
            delta_action -= S(x + mu) + S(x - mu);
          new_spin = (S.lattice().random(x).plain() > 0.5) ? 1 : -1;
          delta_action *= S.beta * (new_spin - S(x));
          if (delta_action < S.lattice().random(x).plain())
            S(x) = new_spin;
        }
      }
    }
    S.update(parity);
  }
}

mdp_real average_spin(ising_field &S)
{
  mdp_real res = 0;
  mdp_site x(S.lattice());
  forallsites(x)
  {
    res += S(x);
  }
  mdp.add(res);
  return res / (S.lattice().global_volume());
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  constexpr int Nconfig = 10;
  int mybox[] = {128, 128};
  mdp_lattice mylattice(2, mybox);
  ising_field S(mylattice);  /* ths spin field +1 or -1     */
  scalar_field H(mylattice); /* the external magnetic field */
  mdp_site x(mylattice);
  mdp_jackboot jb(Nconfig, 1);
  jb.plain(0);

  S.beta = 0.5;         /* inverse square temperature             */
  S.kappa = 0.1;        /* coupling with the total extarnal field */
  S.magnetic_field = 1; /* an extra external magnetic field       */

  for (S.beta = 0.01; S.beta < 0.15; S.beta += 0.001)
  {
    forallsites(x)
    {
      S(x) = 1;
      H(x) = 0;
    }
    S.update();
    H.update();

    for (int conf = 0; conf < Nconfig; conf++)
    {
      montecarlo_multihit(S, H);
      jb(conf, 0) = average_spin(S);
    }

    mdp << "beta = " << S.beta << ", "
        << "average plaquette = " << jb.mean()
        << "(" << jb.j_err() << ")\n";
  }

  mdp.close_wormholes();

  return 0;
}
