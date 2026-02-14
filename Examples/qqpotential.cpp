#include "fermiqcd.h" // include FermiQCD libraries

using namespace MDP;

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv); // START
  constexpr Box L = {8, 8, 8, 8}; // lattice volume
  int n = 3;                      // SU(n) gauge group
  constexpr int N = 100;          // number of gauge configurations
  mdp_lattice lattice(L);         // make a 4D lattice
  gauge_field U(lattice, n);      // make a gauge field U
  coefficients gauge;             // set physical parameters
  gauge["beta"] = 6.0;            // beta=6/g^2 sets lattice spacing
  mdp_array<float, 1> P(4);       // make a 1D array of 4 float elements
  set_hot(U);                     // make a random gauge configuration
  // consider only square paths in the 1,2 plane
  // each path is an array of {verse,direction}
  Path path1x1 = {{1, 1}, {1, 2}, {-1, 1}, {-1, 2}};
  Path path2x2 = {{1, 1}, {1, 1}, {1, 2}, {1, 2}, {-1, 1}, {-1, 1}, {-1, 2}, {-1, 2}};
  Path path3x3 = {{1, 1}, {1, 1}, {1, 1}, {1, 2}, {1, 2}, {1, 2}, {-1, 1}, {-1, 1}, {-1, 1}, {-1, 2}, {-1, 2}, {-1, 2}};
  P(1) = P(2) = P(3) = 0; // zero the output
  for (int k = 0; k < N; k++)
  {                                              // loop over the MCMC
    WilsonGaugeAction::heatbath(U, gauge, 10);   // do 10 MCMC steps
    P(1) += -real(average_path(U, path1x1));  // compute 1x1 path
    P(2) += -real(average_path(U, path2x2));  // compute 2x2 path
    P(3) += -real(average_path(U, path3x3)); // compute 3x3 path
  }
  for (int r = 1; r < 4; r++) // print output
    mdp << "V(" << r << ")=" << -log(P(r) / N) / r << "\n";
  mdp.close_wormholes(); // STOP
  return 0;
}
