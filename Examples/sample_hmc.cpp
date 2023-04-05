#include <cmath>
#include <cstdlib>
#include "fermiqcd.h"
// #include "fermiqcd_cg_inverter.hpp"

using namespace MDP;

// #include <sstream>
// const char* output(std::string file, int num){
//   std::ostringstream ostr;
//   std::string os;
//   ostr << file << num << ends;
//   os=ostr.str();
//   return os.c_str();
// }

class parameters
{
public:
  int L, T, ndim, nspin, nc, read_in;
  int therm, sweeps, gap, seed, nf;
  mdp_real trajectory_length, beta, timestep, mass;
  mdp_real cg_absolute_precision, cg_relative_precision, cg_max_steps;
  std::string representation;
  parameters()
  {
    L = 6;
    T = 6;
    ndim = 4;
    nspin = 4;
    nc = 2;
    therm = 50;      // unused in tests
    sweeps = 100000; // unused in tests
    seed = 1;        // unused in tests
    nf = 2;
    gap = 1; // unused in tests
    timestep = 0.01;
    trajectory_length = 1;
    beta = 1.5;
    mass = 0.0;
    cg_absolute_precision = 1e-10;
    cg_relative_precision = 0;
    cg_max_steps = 1000;
    representation = "FUNDAMENTAL";
    // representation="SYMMETRIC";
    read_in = 0;
  }
  void read(std::string filename)
  {
    L = (int)val(prompt(filename, "L", "8"));
    T = (int)val(prompt(filename, "T", "16"));
    ndim = (int)val(prompt(filename, "NDIM", "4"));
    nspin = (int)val(prompt(filename, "NSPIN", "4"));
    nc = (int)val(prompt(filename, "NC", "2"));
    therm = (int)val(prompt(filename, "THERMALIZATION_STEPS", "100"));
    sweeps = (int)val(prompt(filename, "SWEEPS", "400"));
    seed = (int)val(prompt(filename, "SEED", "1"));
    nf = (int)val(prompt(filename, "NUMBER_FLAVORS", "2"));
    gap = (int)val(prompt(filename, "GAP", "10"));
    timestep = val(prompt(filename, "TIMESTEP", "0.01"));
    trajectory_length = val(prompt(filename, "TRAJECTORY_LENGTH", "100"));
    beta = val(prompt(filename, "BETA", "1.50"));
    mass = val(prompt(filename, "MASS", "0.00"));
    cg_absolute_precision = val(prompt(filename, "CG_ABSOLUTE_PRECISION", "1e-14"));
    cg_relative_precision = val(prompt(filename, "CG_RELATIVE_PRECISION", "0"));
    cg_max_steps = val(prompt(filename, "CG_MAX_STEPS", "1000"));
    // representation=prompt(filename,"REPRESENTATION","SYMMETRIC");
    representation = prompt(filename, "REPRESENTATION", "FUNDAMENTAL");
    read_in = (int)val(prompt(filename, "READIN", "0"));
  }
};

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  // double starttime;

  parameters param;
  // param.read("input");

  define_base_matrices("FERMILAB");
  int box[] = {param.T, param.L, param.L, param.L};
  // int sweep, number=0;
  mdp_lattice lattice(param.ndim, box);
  // unsigned int r=param.seed;

  gauge_field U(lattice, param.nc);
  int dimrep = 1;
  if (param.representation == "SYMMETRIC")
    dimrep = U.nc * (U.nc + 1) / 2;
  else if (param.representation == "FUNDAMENTAL")
    dimrep = U.nc;
  fermi_field psi(lattice, dimrep, param.nspin);
  // mdp_real pl;
  coefficients coeff;
  coeff["time"] = param.T;
  coeff["length"] = param.L;
  coeff["ndim"] = param.ndim;
  coeff["dimspin"] = param.nspin;
  coeff["mass"] = param.mass;
  coeff["kappa"] = 0.5 / (param.mass + param.ndim);
  coeff["beta"] = param.beta;
  coeff["dynamical_quarks"] = param.nf;
  coeff["measure"] = 0;
  if (param.representation == "SYMMETRIC")
    coeff["representation"] = 1;
  else if (param.representation == "FUNDAMENTAL")
    coeff["representation"] = 0;
  else
    throw std::string("representation not supported");
  coeff["timestep"] = param.timestep;
  coeff["trajectory_length"] = param.trajectory_length;
  coeff["cg_absolute_precision"] = param.cg_absolute_precision;
  coeff["cg_relative_precision"] = param.cg_relative_precision;
  coeff["cg_max_steps"] = param.cg_max_steps;

  default_fermi_action = FermiCloverActionFast::mul_Q;
#ifdef SSE2
  default_fermi_action = FermiCloverActionSSE2::mul_Q;
#endif

  set_hot(U);
  std::cout << "average_plaquette = " << average_plaquette(U) << std::endl;
  set_random(psi);

  HMC<gauge_field, fermi_field> hmc(U, psi, coeff);
  for (int i = 0; i < 100; i++)
  {
    hmc.step();
    mdp << "acceptance=" << (float)hmc.acceptance_rate() << "\n";
  }

  mdp.close_wormholes();
  return 0;
}
