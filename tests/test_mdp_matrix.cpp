#include <iostream>
// #include <random>
#include <chrono>
#include <cassert>
#include <iomanip>
#include "mdp.h"

using namespace MDP;

// std::mt19937 rng(0);
// std::uniform_real_distribution<double> dist(-2.0, 2.0);

// mdp_complex rand_complex()
// {
//   return mdp_complex(dist(rng), dist(rng));
// }

void fill_random(mdp_matrix &A, bool su_matrix = false)
{
  if (su_matrix)
  {
    const int n = A.rows();
    A = Random.SU(n);
  }
  else
  {
    for (auto &x : A)
      x = mdp_complex(Random.gaussian(1.274), Random.gaussian(1.472)); // rand_complex();
  }
}

void print_test(const std::string &name)
{
  std::cout << "---- " << name << " ----\n";
}

void test_size(int n, bool su_matrix = false)
{
  std::cout << "\n===========================\n";
  std::cout << "Testing size " << n << "x" << n << "\n";
  std::cout << "===========================\n";

  mdp_matrix A(n, n);
  mdp_matrix B(n, n);

  fill_random(A, su_matrix);
  std::cout << "Created matrix A:\n"
            << A << "\n";
  fill_random(B, su_matrix);
  std::cout << "Created matrix B:\n"
            << B << "\n";

  print_test("Constructors");
  mdp_matrix C(A);
  std::cout << "C (copy of A):\n"
            << C << "\n";
  mdp_matrix D = B;
  std::cout << "D (copy of B):\n"
            << D << "\n";

  print_test("operator()");
  auto x = A(0, 0);
  std::cout << "A(0, 0): " << x << "\n";
  A(0, 0) = x;
  std::cout << "A(0, 0) after assignment: " << A(0, 0) << "\n";

  print_test("operator[]");
  auto y = A[0];
  std::cout << "A[0]: " << y << "\n";
  A[0] = y;
  std::cout << "A[0] after assignment: " << A[0] << "\n";

  print_test("Iterators");
  mdp_complex sum = 0;
  for (auto &v : A)
    sum += v;
  std::cout << "Sum of elements in A: " << sum << "\n";

  print_test("Unary operators");
  mdp_matrix U1 = +A;
  std::cout << "U1 (unary +A):\n"
            << U1 << "\n";
  mdp_matrix U2 = -A;
  std::cout << "U2 (unary -A):\n"
            << U2 << "\n";

  print_test("Matrix + -");
  mdp_matrix E = A + B;
  std::cout << "E (A + B):\n"
            << E << "\n";
  mdp_matrix F = A - B;
  std::cout << "F (A - B):\n"
            << F << "\n";

  print_test("Matrix *");
  mdp_matrix G = A * B;
  std::cout << "G (A * B):\n"
            << G << "\n";

  print_test("Matrix /");
  mdp_matrix H = A / B;
  std::cout << "H (A / B):\n"
            << H << "\n";

  print_test("Compound operators");
  mdp_matrix T = A;
  std::cout << "T (initially A):\n"
            << T << "\n";
  T += B;
  std::cout << "T after T += B:\n"
            << T << "\n";
  T -= B;
  std::cout << "T after T -= B:\n"
            << T << "\n";
  T *= B;
  std::cout << "T after T *= B:\n"
            << T << "\n";
  T /= B;
  std::cout << "T after T /= B:\n"
            << T << "\n";

  print_test("Scalar operations");

  mdp_complex s = mdp_complex(2.0);
  std::cout << "Scalar s: " << s << "\n";

  mdp_matrix S1 = A + s;
  std::cout << "S1 (A + s):\n"
            << S1 << "\n";
  mdp_matrix S2 = A - s;
  std::cout << "S2 (A - s):\n"
            << S2 << "\n";
  mdp_matrix S3 = A * s;
  std::cout << "S3 (A * s):\n"
            << S3 << "\n";
  mdp_matrix S4 = A / s;
  std::cout << "S4 (A / s):\n"
            << S4 << "\n";

  A += s;
  std::cout << "A after A += s:\n"
            << A << "\n";
  A -= s;
  std::cout << "A after A -= s:\n"
            << A << "\n";
  A *= s;
  std::cout << "A after A *= s:\n"
            << A << "\n";
  A /= s;
  std::cout << "A after A /= s:\n"
            << A << "\n";

  print_test("Global scalar ops");

  mdp_matrix GS1 = s + A;
  std::cout << "GS1 (s + A):\n"
            << GS1 << "\n";
  mdp_matrix GS2 = s - A;
  std::cout << "GS2 (s - A):\n"
            << GS2 << "\n";
  mdp_matrix GS3 = s * A;
  std::cout << "GS3 (s * A):\n"
            << GS3 << "\n";
  mdp_matrix GS4 = s / A;
  std::cout << "GS4 (s / A):\n"
            << GS4 << "\n";

  print_test("Identity / Zero");

  mdp_matrix I = mdp_identity(n);
  std::cout << "Identity matrix I:\n"
            << I << "\n";
  mdp_matrix Z = mdp_zero(n);
  std::cout << "Zero matrix Z:\n"
            << Z << "\n";

  print_test("Trace");
  auto tr = trace(A);
  std::cout << "Trace: " << tr << "\n";

  print_test("Determinant");
  auto d = det(A);
  std::cout << "Determinant: " << d << "\n";

  print_test("Inverse");
  mdp_matrix invA = inv(A);
  std::cout << "inv(A):\n"
            << invA << "\n";

  print_test("Transpose");
  mdp_matrix At = transpose(A);
  std::cout << "Transpose of A:\n"
            << At << "\n";

  print_test("Hermitian");
  mdp_matrix Ah = hermitian(A);
  std::cout << "Hermitian of A:\n"
            << Ah << "\n";

  print_test("Conjugate");
  mdp_matrix Ac = conj(A);
  std::cout << "Conjugate of A:\n"
            << Ac << "\n";

  print_test("Submatrix");
  if (n > 2)
  {
    auto sub = submatrix(A, 1, 1);
    std::cout << "Submatrix of A (removing row 1 and column 1):\n"
              << sub << "\n";
  }

  print_test("max()");
  auto m = max(A);
  std::cout << "max: " << m << "\n";

  print_test("Power");
  mdp_matrix P = pow(A, 2);
  std::cout << "A^2:\n"
            << P << "\n";

  print_test("exp");
  mdp_matrix ex = exp(A);
  std::cout << "exp(A):\n"
            << ex << "\n";

  print_test("log");
  mdp_matrix lg = log(I + A * 0.01);
  std::cout << "log(I + 0.01 * A):\n"
            << lg << "\n";

  print_test("sin");
  mdp_matrix sn = sin(A);
  std::cout << "sin(A):\n"
            << sn << "\n";

  print_test("cos");
  mdp_matrix cs = cos(A);
  std::cout << "cos(A):\n"
            << cs << "\n";

  print_test("Comparison");
  bool eq = (A == A);
  bool neq = (A != B);
  std::cout << "A == A: " << eq << "\n";
  std::cout << "A != B: " << neq << "\n";

  print_test("operator<<");
  std::cout << A << "\n";
}

bool mdp_matrix_test()
{
  mdp_matrix a(3, 3), b, c;

  a = Random.SU(3);
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

  // assert(max(a - log(exp(a))) < mdp_precision);
  // printf("exp, log                    ...test passed\n");

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

void test_object_size()
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

int main()
{
  std::cout << "Random matrix tests\n";
  test_size(2);
  test_size(3);
  test_size(4);
  test_size(8);

  std::cout << "\nRandom SU(n) matrix tests\n";
  test_size(2, true);
  test_size(3, true);
  test_size(4, true);
  test_size(8, true);

  std::cout << "\nOther tests\n";
  mdp_matrix_test();
  test_object_size();

  std::cout << "\nAll tests executed.\n";
}
