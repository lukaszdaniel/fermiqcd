/////////////////////////////////////////////////////////////////
/// @file mdp_field_test.cpp
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains a sample test (main) function
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////

#include <cassert>
#include "mdp.h"

/// For debugging purposes only
bool mdp_field_test(int argc, char **argv)
{
  mpi.open_wormholes(argc, argv);

  int box[] = {4, 4, 4, 4};
  mdp_lattice lattice(4, box);
  mdp_matrix_field M(lattice, 3, 3);
  mdp_site x(lattice);
  mdp_site y(lattice);
  y.set(0, 1, 0, 0);

  std::cout << "Initial M(0,1,0,0):\n";
  std::cout << M(y) << "\n";
  std::cout << "Final M(0,1,0,0):\n";
  M(y) = M.lattice().random(y).SU(3);
  std::cout << M(y) << "\n";

  double counter = 0;

  forallsites(x)
      M(x) = M.lattice().random(x).SU(3);

  forallsites(x)
      counter += real(trace(M(x) * inv(M(x)))) / 3;

  mpi.add(counter);
  counter /= lattice.size();

  assert((counter - 1) < mdp_precision);
  printf("lattice and field ...test passed.\n");

  mpi.close_wormholes();
  return (mdp_int)counter;
}

int main(int argc, char **argv)
{
  mdp_field_test(argc, argv);
  return 0;
}
