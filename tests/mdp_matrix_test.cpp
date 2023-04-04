/////////////////////////////////////////////////////////////////
/// @file mdp_matrix_test.cpp
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// For debugging only
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////

#include <cassert>
#include "mdp.h"

/// For debugging only
bool mdp_matrix_test()
{
  mdp_matrix a(3, 3), b, c;

  b = a;
  a = a + b; // Random.SU(3);

  assert(max(a * inv(a) - 1) < mdp_precision);
  printf("inversion                   ...test passed\n");

  assert(max(2 * a - mdp_complex(0, -1) * mdp_complex(0, 2) * a) < mdp_precision);
  printf("operator*                   ...test passed\n");

  assert(max(a + 3 * a - 4 * a) < mdp_precision);
  assert(max(a + mdp_complex(3, 2) * a - mdp_complex(4, 2) * a) < mdp_precision);
  printf("operator+ and operator-     ...test passed\n");

  assert(max(exp(mdp_complex(0, 1) * a) - cos(a) - mdp_complex(0, 1) * sin(a)) < mdp_precision);
  printf("exp, sin and cos            ...test passed\n");

  assert(max(a - log(exp(a))) < mdp_precision);
  printf("exp, log                    ...test passed\n");

  b = a;
  b *= mdp_complex(0, 2);
  b /= mdp_complex(0, 2);
  b += mdp_complex(2, 4);
  b -= mdp_complex(2, 4);
  assert(max(a - b) < mdp_precision);
  printf("*=, /=, +=, -= (mdp_complex)    ...test passed\n");

  b = a;
  b *= 2;
  b /= 2;
  b += 3;
  b -= 3;
  assert(max(a - b) < mdp_precision);
  printf("*=, /=, +=, -= (Real)       ...test passed\n");

  b = a;
  b *= a;
  b /= a;
  b += a;
  b -= a;
  assert(max(a - b) < 10.0 * mdp_precision);
  printf("*=, /=, +=, -= (mdp_matrix)     ...test passed\n");

  return 1;
}

void test_size()
{
  mdp_matrix a(3, 3);
  mdp_matrix b(2, 2);
  mdp_matrix c;

  b(0, 0) = 2;
  b(0, 1) = 1;
  b(1, 0) = 1;
  b(1, 1) = 4;

  std::cout << "sizeof(a) = " << sizeof(a) << "\n";
  std::cout << "sizeof(b) = " << sizeof(b) << "\n";
  std::cout << "sizeof(c) = " << sizeof(c) << "\n";
}

int main(int argc, char **argv)
{
  mdp_matrix_test();
  test_size();
  return 0;
}
