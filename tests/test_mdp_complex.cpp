#include <iostream>
#include <cassert>
#include <cmath>
// #include <iomanip>

#include "mdp_complex.h"

using namespace MDP;

static constexpr double EPS = 1e-10;

bool approx(double a, double b, double eps = EPS)
{
  return std::fabs(a - b) < eps;
}

bool approx_complex(const mdp_complex &a, const mdp_complex &b)
{
  return approx(a.real(), b.real()) && approx(a.imag(), b.imag());
}

void test_constructors()
{
  mdp_complex a;
  assert(a.real() == 0.0);
  assert(a.imag() == 0.0);

  mdp_complex b(3.0, 4.0);
  assert(b.real() == 3.0);
  assert(b.imag() == 4.0);

  mdp_complex c(b);
  assert(c == b);

  std::cout << "constructors OK\n";
}

void test_unary()
{
  mdp_complex a(2, -3);

  mdp_complex b = +a;
  assert(b == a);

  mdp_complex c = -a;
  assert(c.real() == -2);
  assert(c.imag() == 3);

  std::cout << "unary operators OK\n";
}

void test_add_sub()
{
  mdp_complex a(1, 2);
  mdp_complex b(3, 4);

  mdp_complex c = a + b;
  assert(c == mdp_complex(4, 6));

  mdp_complex d = b - a;
  assert(d == mdp_complex(2, 2));

  std::cout << "addition/subtraction OK\n";
}

void test_mul()
{
  mdp_complex a(1, 2);
  mdp_complex b(3, 4);

  mdp_complex c = a * b;
  assert(c == mdp_complex(-5, 10));

  std::cout << "multiplication OK\n";
}

void test_div()
{
  mdp_complex a(1, 2);
  mdp_complex b(3, 4);

  mdp_complex c = a / b;

  mdp_complex expected(0.44, 0.08);
  assert(approx_complex(c, expected));

  std::cout << "division OK\n";
}

void test_assignment()
{
  mdp_complex a(5, 6);
  mdp_complex b;

  b = a;
  assert(b == a);

  std::cout << "assignment OK\n";
}

void test_compound()
{
  mdp_complex a(1, 2);
  mdp_complex b(3, 4);

  mdp_complex c = a;
  c += b;
  assert(c == mdp_complex(4, 6));

  c = a;
  c -= b;
  assert(c == mdp_complex(-2, -2));

  c = a;
  c *= b;
  assert(c == mdp_complex(-5, 10));

  c = a;
  c /= b;
  assert(approx_complex(c, mdp_complex(0.44, 0.08)));

  std::cout << "compound operators OK\n";
}

void test_scalar_ops()
{
  mdp_complex a(2, 3);

  mdp_complex b = a + 5.0;
  assert(b == mdp_complex(7, 3));

  mdp_complex c = 5.0 + a;
  assert(c == mdp_complex(7, 3));

  mdp_complex d = a * 2.0;
  assert(d == mdp_complex(4, 6));

  mdp_complex e = 2.0 * a;
  assert(e == mdp_complex(4, 6));

  std::cout << "scalar operators OK\n";
}

void test_helpers()
{
  mdp_complex a(3, 4);

  assert(approx(abs(a), 5.0));
  assert(approx(abs2(a), 25.0));

  mdp_complex c = conj(a);
  assert(c == mdp_complex(3, -4));

  std::cout << "helper functions OK\n";
}

void test_functions()
{
  mdp_complex a(1, 1);

  mdp_complex e = exp(a);
  mdp_complex s = sin(a);
  mdp_complex c = cos(a);

  // std::cout << "exp(a) = " << std::setprecision(12) << e << "\n";
  // std::cout << "sin(a) = " << std::setprecision(12) << s << "\n";
  // std::cout << "cos(a) = " << std::setprecision(12) << c << "\n";
  assert(approx_complex(e, mdp_complex(1.46869393992, 2.28735528718)));
  assert(approx_complex(s, mdp_complex(1.29845758142, 0.634963914785)));
  assert(approx_complex(c, mdp_complex(0.833730025131, -0.988897705763)));

  std::cout << "math functions OK\n";
}

void test_i_helpers()
{
  mdp_complex a(2, 3);

  mdp_complex b = I * (a);
  assert(b == mdp_complex(-3, 2));

  mdp_complex c = -I * a;
  assert(c == mdp_complex(3, -2));

  std::cout << "i helpers OK\n";
}

int main()
{
  std::cout << "Testing mdp_complex...\n";

  test_constructors();
  test_unary();
  test_add_sub();
  test_mul();
  test_div();
  test_assignment();
  test_compound();
  test_scalar_ops();
  test_helpers();
  test_functions();
  test_i_helpers();

  std::cout << "\nALL TESTS PASSED\n";

  return 0;
}
