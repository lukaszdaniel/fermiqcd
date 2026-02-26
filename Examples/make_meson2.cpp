// ////////////////////////////////////////////////////////
// Program C2.C written by Massimo Di Pierro @ July 2000
// ////////////////////////////////////////////////////////
// WORKING EXAMPLE of a Lattice QCD program to compute the
// Euclidean propagator of an Heavy-Light Meson, C2(t)
// To extract m_B and f_B fit output with
//
// C2(t) = 1/2 f_B^2 m_B exp(-m_B t) + ...
//
// and extrapolate to
// mq -> 0    (GeV)
// mh -> mb   (GeV)
// a  -> 0    (GeV^(-1))
//
// NOTE: This program works in parallel with MPI
//
// ////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////
// optional macros
// ////////////////////////////////////////////////////////
// #define USE_DOUBLE_PRECISION

// ////////////////////////////////////////////////////////
// open the libraries: (download from www.fermiqcd.net)
// ////////////////////////////////////////////////////////
#include "fermiqcd.h"

using namespace MDP;

// ////////////////////////////////////////////////////////
// useful macros
// ////////////////////////////////////////////////////////
#define GeV 1

// ////////////////////////////////////////////////////////
// main program
// ////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv); // open communications

  // //////////////////////////////////////////////////////
  // declare parameters of the simulation
  // //////////////////////////////////////////////////////
  int Nt = 12, Nx = 4, Ny = 4, Nz = 4; // lattice size
  int Nc = 3;                          // set colors, SU(Nc)
  // float        beta=5.7;      // set lattice spacing
  float mq = 0.2 * GeV; // set light quark pole-mass
  float mh = 0.7 * GeV; // set heavy quark pole-mass

  // //////////////////////////////////////////////////////
  // additional parameters (they only depend on beta!)
  // //////////////////////////////////////////////////////
  float a = 0.91 / GeV; // lattice spacing
  // float        cSW    =1.57;
  float kappa_c = 0.14315;

  // //////////////////////////////////////////////////////
  // setting the parameters for the light and heavy quarks
  // //////////////////////////////////////////////////////
  coefficients gauge_coeff;
  coefficients light_coeff;
  coefficients heavy_coeff;
  gauge_coeff["beta"] = 5.7;
  light_coeff["kappa"] = 1.0 / ((exp(mq * a) - 1.0) * 2.0 + 1.0 / kappa_c);
  heavy_coeff["kappa"] = 1.0 / ((exp(mh * a) - 1.0) * 2.0 + 1.0 / kappa_c);

  // //////////////////////////////////////////////////////
  // define gamma matrices in Euclidean space
  // //////////////////////////////////////////////////////
  define_base_matrices("FERMILAB");

  // //////////////////////////////////////////////////////
  // choose inversion method and action for computation
  // //////////////////////////////////////////////////////
  float absolute_precision = 1e-12;
  float relative_precision = 1e-8;
  default_fermi_inverter = MinimumResidueInverter<fermi_field, gauge_field>;
#ifdef SSE2
  default_fermi_action = FermiCloverActionSSE2::mul_Q;
#else
  default_fermi_action = FermiCloverActionFast::mul_Q;
#endif

  // //////////////////////////////////////////////////////
  // define the grid size on which the lattice is defined
  // //////////////////////////////////////////////////////
  const Box grid_size = {Nt, Nx, Ny, Nz};

  // //////////////////////////////////////////////////////
  // associate the lattice to the grid
  // //////////////////////////////////////////////////////
  mdp_lattice lattice(grid_size);

  // //////////////////////////////////////////////////////
  // define a gauge field
  // U(x,mu) = exp(i g A(x,mu) )
  // to the sites of the lattice
  // //////////////////////////////////////////////////////
  gauge_field U(lattice, Nc);

  // //////////////////////////////////////////////////////
  // define a light and a heavy propagator
  // S = <0| q(x) \bar q(0) |0>
  // to the sites of the lattice
  // //////////////////////////////////////////////////////
  fermi_propagator Sq(lattice, Nc);
  fermi_propagator Sh(lattice, Nc);

  // //////////////////////////////////////////////////////
  // define a variable site
  // to move on the lattice
  // //////////////////////////////////////////////////////
  mdp_site x(lattice);

  // //////////////////////////////////////////////////////
  // define the number of gauge configurations to be used
  // //////////////////////////////////////////////////////
  constexpr mdp_uint Nconfig = 100;

  // //////////////////////////////////////////////////////
  // define a container object for the propagator
  // this container computes bootstrap error automatically.
  // //////////////////////////////////////////////////////
  mdp_jackboot F(Nconfig, Nt);

  // //////////////////////////////////////////////////////
  // creating initial gauge configuration
  // //////////////////////////////////////////////////////
  set_hot(U);

  // //////////////////////////////////////////////////////
  // create and skip 100 gauge configuration
  // for each print average plaquette to check termalization
  // //////////////////////////////////////////////////////
  for (int i = 0; i < 100; i++)
  {
    WilsonGaugeAction::heatbath(U, gauge_coeff);
    mdp << "average_plaquette=" << average_plaquette(U) << "\n";
  }

  // //////////////////////////////////////////////////////
  // loop over the gauge configurations
  // //////////////////////////////////////////////////////
  for (mdp_uint config = 0; config < Nconfig; config++)
  {

    // skip 10 gauge configurations and
    // for each print average plaquette
    for (int i = 0; i < 10; i++)
    {
      WilsonGaugeAction::heatbath(U, gauge_coeff);
      mdp << "average_plaquette=" << average_plaquette(U) << "\n";
    }

    // compute electromagnetic-field
    compute_em_field(U);

    // compute light propagator Sq
    fermi_propagator::generate(Sq, U, light_coeff, absolute_precision, relative_precision);

    // compute heavy propagator Sh
    fermi_propagator::generate(Sh, U, heavy_coeff, absolute_precision, relative_precision);

    // ////////////////////////////////////////////////////
    // compute the pion propagator by
    // wick contracting
    //
    // C_2(t_x) =
    //   \sum_{x} \bar h(x) Gamma5 q(x) \bar q(0) Gamma5 h(0)
    //
    // as function of t = t_x
    // ////////////////////////////////////////////////////
    for (int t = 0; t < Nt; t++)
    {
      F(0, t) = 0;
    }
    forallsites(x)
    {
      int t = x(0);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          F(config, t) +=
              real(trace(Sq(x, i, j) * hermitian(Sh(x, i, j))));
    }

    // ////////////////////////////////////////////////////
    // for each t print out C_2(t) with the Bootstrap error
    // ////////////////////////////////////////////////////
    mdp << "\nRESULT FOR C2(t) (@ gauge = " << config << ")\n";
    mdp << "==================================\n";
    mdp << "t\tC2\t\t(error)\n";
    mdp << "==================================\n";
    for (int t = 0; t < Nt; t++)
    {
      mdp.barrier(); // syncronize output if in parallel!
      if (on_which_process(lattice, t) == ME)
      {
        F.plain(t);
        std::cout << t << "\t" << F.mean() << "\t" << F.b_err() << "\n";
      }
    }
    mdp << "==================================\n\n";
  }

  mdp.close_wormholes(); // close communications
  return 0;
}
