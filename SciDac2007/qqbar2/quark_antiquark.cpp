#include "../../Libraries/fermiqcd.h"
#include "../mdp_all.h"
#include "../dump.h"

using namespace MDP;

class PunchedWilsonGaugeAction : public WilsonGaugeAction
{
public:
  static gauge_stats heatbath(gauge_field &U,
                              coefficients &coeff,
                              std::vector<mdp_site> &sites,
                              mdp_uint n_iter = 1)
  {
    begin_function("WilsonGaugeAction__heatbath");
    if (U.nc() == 1)
      error("fermiqcd_gauge_algorithms/heatbath(): U(1)? (use metropolis)");
    gauge_stats stats;
    mdp_real beta, zeta;

    if (coeff.has_key("beta"))
      beta = coeff["beta"];
    else
      error("beta undeclared");
    if (coeff.has_key("zeta"))
      zeta = coeff["zeta"];
    else
      zeta = 1;

    mdp_suint nc = U.nc();
    mdp_suint ndim = U.ndim();
    mdp_matrix M(nc, nc);
    mdp_site x(U.lattice());
    mdp_real time = mdp.time();
    mdp_complex a[4];

    mdp << coeff;

    for (mdp_uint iter = 0; iter < n_iter; iter++)
    {
      for (mdp_parity parity : {EVEN, ODD})
      {
        for (mdp_suint mu = 0; mu < ndim; mu++)
        {
          forallsitesofparity(x, parity)
          {
            if (mu == 0)
            {
              for (size_t q = 0; q < sites.size(); q++)
              {
                if (x(1) == sites[q](1) &&
                    x(2) == sites[q](2) &&
                    x(3) == sites[q](3))
                  continue;
              }
            }
            for (mdp_suint i = 0; i < nc - 1; i++)
            {
              for (mdp_suint j = i + 1; j < nc; j++)
              {
                if (zeta == 1)
                  M = U(x, mu) * staple_H(U, x, mu);
                else if (mu == 0)
                  M = zeta * U(x, 0) * staple_H(U, x, 0);
                else
                  M = (1.0 / zeta) * U(x, mu) * staple_H(U, x, mu);
                a[0] = M(i, i);
                a[1] = M(i, j);
                a[2] = M(j, i);
                a[3] = M(j, j);
                heatbath_SU2(U.lattice().random(x), beta / nc, a);
                for (mdp_suint k = 0; k < nc; k++)
                {
                  mdp_complex tmpUik = a[0] * U(x, mu, i, k) + a[1] * U(x, mu, j, k);
                  U(x, mu, j, k) = a[2] * U(x, mu, i, k) + a[3] * U(x, mu, j, k);
                  U(x, mu, i, k) = tmpUik;
                }
              }
            }
          }
          // The next command does all the communication!
          U.update(parity, mu, nc * nc);
        }
      }
    }
    mdp << "\t<stats>\n\t\t<time>" << mdp.time() - time << "</time>\n\t</stats>\n";
    end_function("WilsonGaugeAction__heatbath");
    return stats;
  }
};

void punched_ape_smearing(gauge_field &U,
                          std::vector<mdp_site> &sites,
                          mdp_real alpha = 0.7,
                          mdp_uint iterations = 20,
                          mdp_uint cooling_steps = 10)
{
  gauge_field V(U.lattice(), U.nc());
  mdp_site x(U.lattice());
  for (mdp_uint iter = 0; iter < iterations; iter++)
  {
    std::cout << "smearing step " << iter << "/" << iterations << std::endl;
    V = U;
    for (mdp_suint mu = 0; mu < U.ndim(); mu++)
    {
      forallsites(x)
      {
        if (mu == 0)
          for (size_t q = 0; q < sites.size(); q++)
          {
            if (x(1) == sites[q](1) &&
                x(2) == sites[q](2) &&
                x(3) == sites[q](3))
              continue;
          }
        U(x, mu) = (1.0 - alpha) * V(x, mu);
        for (mdp_suint nu = 0; nu < U.ndim(); nu++)
          if (nu != mu)
            U(x, mu) += (1.0 - alpha) / 6 *
                        (V(x, nu) * V(x + nu, mu) * hermitian(V(x + mu, nu)) +
                         hermitian(V(x - nu, nu)) * V(x - nu, mu) * V((x - nu) + mu, nu));
        U(x, mu) = project_SU(U(x, mu), cooling_steps);
      }
    }
    U.update();
  }
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);
  define_base_matrices("FERMILAB");
  mdp_uint N;
  constexpr mdp_suint nc = 2;
  constexpr Box L = {10, 16, 10, 10};
  constexpr Box L_space = {16, 10, 10};
  int idx[6][3] = {{0, 0, 1}, {1, 0, 2}, {2, 0, 3}, {3, 1, 2}, {4, 1, 3}, {5, 2, 3}};
  mdp_lattice lattice(L);
  mdp_lattice cube(L_space);

  gauge_field U(lattice, nc);
  gauge_field V(lattice, nc);
  gauge_field W(lattice, nc);
  mdp_complex_field PA(cube, 6);
  mdp_complex_field PB(cube, 6);
  mdp_complex P, PC = 0;
  mdp_real_scalar_field Q3(cube);
  mdp_site x(lattice), y(lattice);
  mdp_site x3(cube), x3b(cube);
  mdp_complex d;
  std::string filename;
  coefficients gauge;
  gauge["beta"] = 2.2;
  std::vector<mdp_site> sites;

  Path path = {{+1, 0}, {+1, 0}, {+1, 0}, {+1, 0}, {+1, 1}, {+1, 1}, {+1, 1}, {+1, 1}, {+1, 1}, {+1, 1}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 1}, {-1, 1}, {-1, 1}, {-1, 1}, {-1, 1}, {-1, 1}};

  forallsites(x3)
  {
    for (mdp_suint i = 0; i < 6; i++)
      PA(x3, i) = PB(x3, i) = 0;
  }
  PC = 0;

#if 0
  set_hot(V);
  for (int k = 0; k < 4000; k++)
  {
    std::cout << k << std::endl;
    WilsonGaugeAction::heatbath(V, gauge, 10);
  }
#endif

  for (mdp_suint conf = 0; conf < 1000; conf++)
  {
#if 0
    WilsonGaugeAction::heatbath(V, gauge, 50);
    filename = std::format("gauge_from_hot_10x16x10x10.{:03d}.fixed.mdp", conf);
    V.save(filename);

    U = V;
    ApeSmearing::smear(U, 0.7, 20, 10);
#endif
    filename = std::format("gauge_from_hot_10x16x10x10.{:03d}.fixed.mdp", conf);
    U.load(filename);

    // GaugeFixing::fix(U, GaugeFixing::Landau, 20);

#if 0
    for (mdp_uint shift0 = 0; shift0 < L[0]; shift0++)
    {
#endif
      for (mdp_uint shift1 = 0; shift1 < L[1]; shift1++)
      {
        for (mdp_uint shift2 = 0; shift2 < L[2]; shift2++)
        {
          for (mdp_uint shift3 = 0; shift3 < L[3]; shift3++)
          {

            x.set(3, 5, L[2] / 2, L[3] / 2);

            d = trace(build_path(U, x, path));
            PC += d;
            forallsites(x3)
            {
              x.set(L[0] / 2, x3(0), x3(1), x3(2));
              for (mdp_suint h = 0; h < 6; h++)
              {
                P = trace(plaquette(U, x, idx[h][1], idx[h][2]));
                PA(x3, idx[h][0]) += d * P;
                PB(x3, idx[h][0]) += P;
              }
            }

            forallsites(x)
            {
              for (mdp_suint mu = 0; mu < 4; mu++)
                W(x, mu) = U(x - 3, mu);
            }
            W.update();
            U = W;
          }

          forallsites(x)
          {
            for (mdp_suint mu = 0; mu < 4; mu++)
              W(x, mu) = U(x - 2, mu);
          }
          W.update();
          U = W;
        }

        forallsites(x)
        {
          for (mdp_suint mu = 0; mu < 4; mu++)
            W(x, mu) = U(x - 1, mu);
        }
        W.update();
        U = W;
      }

#if 0
      forallsites(x)
      {
        for (mdp_suint mu = 0; mu < 4; mu++)
          W(x, mu) = U(x - 0, mu);
      }
      W.update();
      U = W;
    }
#endif

    N = (conf + 1) * L[0];
    mdp_real m = 0;
    mdp_real Q;
    forallsites(x3) Q3(x3) = 0;
    forallsites(x3)
    {
      Q = 0;
      for (mdp_suint q = 0; q < 3; q++)
        Q += (real(PA(x3, q) / PC - PB(x3, q) / N));
      for (mdp_suint q = 3; q < 6; q++)
        Q += (real(PA(x3, q) / PC - PB(x3, q) / N));
      Q3(x3) += Q;
      // x3b.set(L[1] - x3(0), x3(1), x3(2));
      // Q3(x3b) += Q;
      // x3b.set(x3(0), L[2] - x3(1), x3(2));
      // Q3(x3b) += Q;
      // x3b.set(x3(0), x3(1), L[3] - x3(2));
      // Q3(x3b) += Q;
      // x3b.set(L[1] - x3(0), L[2] - x3(1), x3(2));
      // Q3(x3b) += Q;
      // x3b.set(L[1] - x3(0), x3(1), L[3] - x3(2));
      // Q3(x3b) += Q;
      // x3b.set(x3(0), L[2] - x3(1), L[3] - x3(2));
      // Q3(x3b) += Q;
      // x3b.set(L[1] - x3(0), L[2] - x3(1), L[3] - x3(2));
      // Q3(x3b) += Q;
    }

    // x3.set(0, 0, 0);
    // Q3(x3) = 0;

    forallsites(x3)
    {
      if (Q3(x3) < m)
        m = Q3(x3);
    }
    std::cout << "m=" << m << std::endl;
    forallsites(x3)
    {
      Q3(x3) -= m;
    }
    std::cout << "saving vtk file\n";
    filename = std::format("energy_density_E2B2_XYZ.{:03d}.vtk", conf);
    dump(Q3, 0, filename);
  }

  mdp.close_wormholes();
  return 0;
}
