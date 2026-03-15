#include <iostream>
#include <cassert>
#include "mdp.h"

using namespace MDP;

bool approx(double a, double b, double eps = mdp_precision)
{
  return std::fabs(a - b) < eps;
}

bool approx_complex(const mdp_complex &a, const mdp_complex &b)
{
  return approx(a.real(), b.real()) && approx(a.imag(), b.imag());
}

bool compare_equal(const mdp_matrix &a, const mdp_matrix &b)
{
  if (a.rows() != b.rows() || a.cols() != b.cols())
    return false;

  for (mdp_uint i = 0; i < a.rows(); ++i)
    for (mdp_uint j = 0; j < b.cols(); ++j)
    {
      if (!approx_complex(a(i, j), b(i, j)))
        return false;
    }

  return true;
}

mdp_matrix exp_slow(const mdp_matrix &a)
{
#ifdef CHECK_ALL
  if (a.rows() != a.cols())
    error("exp(...)\nmdp_matrix is not squared");
#endif

  const mdp_uint n = a.rows();

  mdp_matrix result = mdp_identity(n);
  mdp_matrix term = a;

  result += term;

  for (mdp_uint i = 2; max(term) > mdp_precision; ++i)
  {
    term *= a;
    term *= 1.0 / i;
    result += term;
  }

  return result;
}

mdp_matrix log_slow(const mdp_matrix &a)
{
#ifdef CHECK_ALL
  if (a.rows() != a.cols())
    error("log(...)\nmdp_matrix is not squared");
#endif

  const mdp_uint n = a.rows();

  mdp_matrix b = a - mdp_identity(n);
  mdp_matrix term = b;
  mdp_matrix result = term;

  for (mdp_uint k = 2; max(term) > mdp_precision; ++k)
  {
    term *= b;

    mdp_real coef = ((k % 2) ? 1.0 : -1.0) / k;
    result += coef * term;
  }

  return result;
}

std::pair<mdp_matrix, mdp_matrix> sincos_slow(const mdp_matrix &a)
{
#ifdef CHECK_ALL
  if (a.rows() != a.cols())
    error("sincos(...)\nmdp_matrix is not squared");
#endif

  mdp_matrix a2 = a * a;

  mdp_matrix sin_sum = a;
  mdp_matrix cos_sum = mdp_identity(a.rows());

  mdp_matrix sin_t = a;
  mdp_matrix cos_t = mdp_identity(a.rows());

  mdp_uint i_sin = 1;
  mdp_uint i_cos = 0;

  do
  {
    // next sin term
    ++i_sin;
    sin_t = (-1.0 / i_sin) * sin_t * a2;
    ++i_sin;
    sin_t *= 1.0 / i_sin;
    sin_sum += sin_t;

    // next cos term
    ++i_cos;
    cos_t = (-1.0 / i_cos) * cos_t * a2;
    ++i_cos;
    cos_t *= 1.0 / i_cos;
    cos_sum += cos_t;

  } while (max(sin_t) > mdp_precision || max(cos_t) > mdp_precision);

  return {sin_sum, cos_sum};
}

mdp_matrix sin_slow(const mdp_matrix &a)
{
  return sincos_slow(a).first;
}

mdp_matrix cos_slow(const mdp_matrix &a)
{
  return sincos_slow(a).second;
}

void print_test(const std::string &name)
{
  std::cout << "---- " << name << " ----\n";
}

void test_size(const mdp_matrix a)
{
  const int n = a.rows();
  std::cout << "\n===========================\n";
  std::cout << "Testing size " << n << "x" << n << "\n";
  std::cout << "===========================\n";

  print_test("Constructors");
  mdp_matrix A(a);
  std::cout << "A (copy of a):\n"
            << A << "\n";

  mdp_matrix B(a);
  print_test("operator()");
  auto x = B(0, 0);
  std::cout << "B(0, 0): " << x << "\n";
  B(0, 0) = mdp_complex(0.629574, 0.162009);
  std::cout << "B(0, 0) after assignment of mdp_complex(0.629574, 0.162009): " << B(0, 0) << "\n";

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
  mdp_matrix ex_slow = exp_slow(A);
  if (!compare_equal(ex, ex_slow))
  {
    std::cout << "exp(A) differs\n\n";
  }
  std::cout << "exp(A):\n"
            << ex << "\n";

  print_test("log");
  mdp_matrix lg = log(I + A * 0.01);
  mdp_matrix lg_slow = log_slow(I + A * 0.01);
  if (!compare_equal(lg, lg_slow))
  {
    std::cout << "log(I + 0.01 * A) differs\n\n";
  }
  std::cout << "log(I + 0.01 * A):\n"
            << lg << "\n";

  print_test("sin");
  mdp_matrix sn = sin(A);
  mdp_matrix sn_slow = sin_slow(A);
  if (!compare_equal(sn, sn_slow))
  {
    std::cout << "sin(A) differs\n\n";
  }
  std::cout << "sin(A):\n"
            << sn << "\n";

  print_test("cos");
  mdp_matrix cs = cos(A);
  mdp_matrix cs_slow = cos_slow(A);
  if (!compare_equal(cs, cs_slow))
  {
    std::cout << "cos(A) differs\n\n";
  }
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
  mdp_matrix A(2, 2);
  A(0, 0) = mdp_complex(1.37706, 0.371378);
  A(0, 1) = mdp_complex(1.38901, 1.43178);
  A(1, 0) = mdp_complex(-0.462473, 0.494255);
  A(1, 1) = mdp_complex(-1.77315, -0.809862);

  mdp_matrix B(3, 3);
  B(0, 0) = mdp_complex(1.82862, -0.527034);
  B(0, 1) = mdp_complex(1.48035, -1.4386);
  B(0, 2) = mdp_complex(1.20364, -0.105568);
  B(1, 0) = mdp_complex(0.715518, 0.0819099);
  B(1, 1) = mdp_complex(0.328079, 0.882531);
  B(1, 2) = mdp_complex(1.03446, 0.149493);
  B(2, 0) = mdp_complex(-0.105598, -1.57637);
  B(2, 1) = mdp_complex(0.947673, -1.25467);
  B(2, 2) = mdp_complex(-1.45913, -1.1338);

  mdp_matrix C(4, 4);
  C(0, 0) = mdp_complex(0.611161, 1.8358);
  C(0, 1) = mdp_complex(1.9812, 0.540236);
  C(0, 2) = mdp_complex(-0.342526, 0.327401);
  C(0, 3) = mdp_complex(0.49404, -0.10121);
  C(1, 0) = mdp_complex(0.699009, -0.64797);
  C(1, 1) = mdp_complex(1.11338, -0.731193);
  C(1, 2) = mdp_complex(0.650107, 1.79828);
  C(1, 3) = mdp_complex(0.491384, -1.94571);
  C(2, 0) = mdp_complex(1.88778, 0.694639);
  C(2, 1) = mdp_complex(0.0384975, 1.51277);
  C(2, 2) = mdp_complex(-0.195363, -1.77714);
  C(2, 3) = mdp_complex(-0.233156, -1.92005);
  C(3, 0) = mdp_complex(-0.562222, 1.91835);
  C(3, 1) = mdp_complex(0.754645, -0.0764259);
  C(3, 2) = mdp_complex(1.67294, 1.5219);
  C(3, 3) = mdp_complex(0.260755, -1.13271);

  mdp_matrix D(8, 8);
  D(0, 0) = mdp_complex(-0.747127, -1.00081);
  D(0, 1) = mdp_complex(0.35386, 1.86166);
  D(0, 2) = mdp_complex(0.132825, 0.638674);
  D(0, 3) = mdp_complex(-0.420523, -1.07787);
  D(0, 4) = mdp_complex(-0.10053, 0.475234);
  D(0, 5) = mdp_complex(0.864298, -0.119471);
  D(0, 6) = mdp_complex(-0.466151, -0.848036);
  D(0, 7) = mdp_complex(1.51381, 0.996679);
  D(1, 0) = mdp_complex(-1.6305, -1.58855);
  D(1, 1) = mdp_complex(0.207265, -0.583813);
  D(1, 2) = mdp_complex(1.87585, -1.8655);
  D(1, 3) = mdp_complex(-1.11495, -0.716011);
  D(1, 4) = mdp_complex(-1.61096, -1.43494);
  D(1, 5) = mdp_complex(-0.958636, 1.93617);
  D(1, 6) = mdp_complex(-0.208295, 0.14809);
  D(1, 7) = mdp_complex(-0.590753, -1.60172);
  D(2, 0) = mdp_complex(1.36456, -0.123003);
  D(2, 1) = mdp_complex(-1.84976, 1.61859);
  D(2, 2) = mdp_complex(-1.33261, 0.0332618);
  D(2, 3) = mdp_complex(1.45973, 1.1162);
  D(2, 4) = mdp_complex(-1.44011, -0.354413);
  D(2, 5) = mdp_complex(1.9303, -1.86711);
  D(2, 6) = mdp_complex(-0.319699, -0.506837);
  D(2, 7) = mdp_complex(-0.538016, -1.79765);
  D(3, 0) = mdp_complex(-1.07703, -1.93349);
  D(3, 1) = mdp_complex(1.77649, 1.05965);
  D(3, 2) = mdp_complex(-0.642385, 0.999997);
  D(3, 3) = mdp_complex(-0.64406, -0.0418042);
  D(3, 4) = mdp_complex(-1.31605, -1.28204);
  D(3, 5) = mdp_complex(1.49829, -0.146196);
  D(3, 6) = mdp_complex(0.433011, 1.77648);
  D(3, 7) = mdp_complex(1.13458, 0.386622);
  D(4, 0) = mdp_complex(-1.79852, 0.000105175);
  D(4, 1) = mdp_complex(1.96959, 0.796392);
  D(4, 2) = mdp_complex(0.716362, -0.93095);
  D(4, 3) = mdp_complex(1.00338, 1.45713);
  D(4, 4) = mdp_complex(0.21697, 1.85796);
  D(4, 5) = mdp_complex(-1.11023, -1.15044);
  D(4, 6) = mdp_complex(0.278294, -1.125);
  D(4, 7) = mdp_complex(1.88095, -0.191564);
  D(5, 0) = mdp_complex(-1.65882, 0.722179);
  D(5, 1) = mdp_complex(-0.0486492, -1.77433);
  D(5, 2) = mdp_complex(1.90562, 1.52402);
  D(5, 3) = mdp_complex(0.169995, 0.470632);
  D(5, 4) = mdp_complex(0.975338, 1.41845);
  D(5, 5) = mdp_complex(0.708326, -0.0856147);
  D(5, 6) = mdp_complex(0.858788, 0.42818);
  D(5, 7) = mdp_complex(-0.175942, -0.122011);
  D(6, 0) = mdp_complex(-1.45112, 1.62567);
  D(6, 1) = mdp_complex(1.52634, -1.08312);
  D(6, 2) = mdp_complex(0.583138, 1.6177);
  D(6, 3) = mdp_complex(0.0788448, -0.701268);
  D(6, 4) = mdp_complex(-0.752559, -1.99978);
  D(6, 5) = mdp_complex(1.54135, -0.298194);
  D(6, 6) = mdp_complex(-0.175481, 0.719518);
  D(6, 7) = mdp_complex(1.15496, -0.0663655);
  D(7, 0) = mdp_complex(1.52119, -1.08223);
  D(7, 1) = mdp_complex(1.8298, -0.74523);
  D(7, 2) = mdp_complex(0.846335, -0.112994);
  D(7, 3) = mdp_complex(0.921769, -1.38522);
  D(7, 4) = mdp_complex(-1.14048, 0.585058);
  D(7, 5) = mdp_complex(1.23032, -1.25417);
  D(7, 6) = mdp_complex(0.699389, 0.988318);
  D(7, 7) = mdp_complex(-1.30036, -0.892425);

  std::cout << "Random matrix tests\n";
  test_size(A);
  test_size(B);
  test_size(C);
  test_size(D);

  A(0, 0) = mdp_complex(0.467606, -0.568824);
  A(0, 1) = mdp_complex(0.326484, 0.592614);
  A(1, 0) = mdp_complex(-0.326484, 0.592614);
  A(1, 1) = mdp_complex(0.467606, 0.568824);

  B(0, 0) = mdp_complex(-0.222349, -0.344606);
  B(0, 1) = mdp_complex(-0.165495, -0.153194);
  B(0, 2) = mdp_complex(-0.875952, 0.116874);
  B(1, 0) = mdp_complex(-0.236984, 0.722386);
  B(1, 1) = mdp_complex(-0.462023, 0.176074);
  B(1, 2) = mdp_complex(-0.215822, -0.361872);
  B(2, 0) = mdp_complex(-0.40385, -0.301181);
  B(2, 1) = mdp_complex(0.199874, 0.815308);
  B(2, 2) = mdp_complex(0.013522, -0.203313);

  C(0, 0) = mdp_complex(0.482691, -0.182589);
  C(0, 1) = mdp_complex(0.235091, -0.275466);
  C(0, 2) = mdp_complex(-0.150781, -0.406174);
  C(0, 3) = mdp_complex(-0.287186, -0.576484);
  C(1, 0) = mdp_complex(-0.198814, -0.16848);
  C(1, 1) = mdp_complex(0.235064, 0.768484);
  C(1, 2) = mdp_complex(-0.491037, -0.103469);
  C(1, 3) = mdp_complex(0.00481474, -0.18552);
  C(2, 0) = mdp_complex(0.691404, 0.0542527);
  C(2, 1) = mdp_complex(-0.11887, 0.3891);
  C(2, 2) = mdp_complex(0.114708, 0.508223);
  C(2, 3) = mdp_complex(-0.275961, 0.076722);
  C(3, 0) = mdp_complex(-0.428651, 0.0321493);
  C(3, 1) = mdp_complex(-0.157966, 0.180405);
  C(3, 2) = mdp_complex(0.528411, 0.0989837);
  C(3, 3) = mdp_complex(-0.479224, -0.488931);

  D(0, 0) = mdp_complex(0.0467698, -0.054417);
  D(0, 1) = mdp_complex(-0.00652483, -0.0111906);
  D(0, 2) = mdp_complex(0.0130637, -0.0136811);
  D(0, 3) = mdp_complex(-0.065583, 0.00840794);
  D(0, 4) = mdp_complex(-0.237801, -0.00795271);
  D(0, 5) = mdp_complex(-0.0270772, -0.00289287);
  D(0, 6) = mdp_complex(-0.153167, 0.786219);
  D(0, 7) = mdp_complex(0.50962, -0.176879);
  D(1, 0) = mdp_complex(0.314162, 0.157261);
  D(1, 1) = mdp_complex(-0.0563457, -0.107896);
  D(1, 2) = mdp_complex(0.0942334, 0.177524);
  D(1, 3) = mdp_complex(0.280064, -0.446779);
  D(1, 4) = mdp_complex(0.19632, 0.154332);
  D(1, 5) = mdp_complex(0.127958, 0.149584);
  D(1, 6) = mdp_complex(0.253371, -0.211605);
  D(1, 7) = mdp_complex(0.568894, 0.0979331);
  D(2, 0) = mdp_complex(-0.369263, 0.167828);
  D(2, 1) = mdp_complex(-0.025795, 0.0921083);
  D(2, 2) = mdp_complex(0.249461, 0.280965);
  D(2, 3) = mdp_complex(0.314787, -0.0810046);
  D(2, 4) = mdp_complex(0.58867, -0.299998);
  D(2, 5) = mdp_complex(0.0438447, 0.0518746);
  D(2, 6) = mdp_complex(-0.115724, 0.263086);
  D(2, 7) = mdp_complex(-0.138026, -0.19158);
  D(3, 0) = mdp_complex(-0.518813, 0.109263);
  D(3, 1) = mdp_complex(0.218869, -0.142409);
  D(3, 2) = mdp_complex(0.433253, -0.159455);
  D(3, 3) = mdp_complex(-0.00909684, -0.331774);
  D(3, 4) = mdp_complex(-0.286739, 0.434553);
  D(3, 5) = mdp_complex(0.00604235, 0.1133);
  D(3, 6) = mdp_complex(-0.169902, -0.0883752);
  D(3, 7) = mdp_complex(-0.0206458, -0.0799267);
  D(4, 0) = mdp_complex(0.169139, 0.385861);
  D(4, 1) = mdp_complex(0.366076, 0.0706944);
  D(4, 2) = mdp_complex(-0.264175, -0.143269);
  D(4, 3) = mdp_complex(0.41059, 0.125903);
  D(4, 4) = mdp_complex(-0.152028, -0.0567643);
  D(4, 5) = mdp_complex(-0.0941076, -0.0473815);
  D(4, 6) = mdp_complex(-0.247764, -0.198825);
  D(4, 7) = mdp_complex(0.0615746, -0.516331);
  D(5, 0) = mdp_complex(-0.0133888, -0.155558);
  D(5, 1) = mdp_complex(0.296539, -0.735655);
  D(5, 2) = mdp_complex(0.202315, -0.0648864);
  D(5, 3) = mdp_complex(0.0692384, 0.255873);
  D(5, 4) = mdp_complex(0.142289, -0.158811);
  D(5, 5) = mdp_complex(-0.0365092, -0.392316);
  D(5, 6) = mdp_complex(0.0982427, -0.0434956);
  D(5, 7) = mdp_complex(0.134994, 0.0247533);
  D(6, 0) = mdp_complex(-0.192755, 0.378598);
  D(6, 1) = mdp_complex(-0.215334, -0.189657);
  D(6, 2) = mdp_complex(-0.189217, 0.566869);
  D(6, 3) = mdp_complex(-0.215256, 0.0310179);
  D(6, 4) = mdp_complex(-0.28975, -0.0251216);
  D(6, 5) = mdp_complex(0.386708, -0.285543);
  D(6, 6) = mdp_complex(-0.0142949, -0.0832885);
  D(6, 7) = mdp_complex(0.016477, -0.0982466);
  D(7, 0) = mdp_complex(0.208386, 0.00238859);
  D(7, 1) = mdp_complex(-0.145974, -0.189477);
  D(7, 2) = mdp_complex(0.334975, -0.011556);
  D(7, 3) = mdp_complex(0.0747648, 0.440663);
  D(7, 4) = mdp_complex(-0.120851, -0.0494826);
  D(7, 5) = mdp_complex(0.429584, 0.60088);
  D(7, 6) = mdp_complex(-0.0945577, -0.0399801);
  D(7, 7) = mdp_complex(-0.0637924, -0.0999159);

  std::cout << "\nRandom SU(n) matrix tests\n";
  test_size(A);
  test_size(B);
  test_size(C);
  test_size(D);

  std::cout << "\nOther tests\n";
  mdp_matrix_test();
  test_object_size();

  std::cout << "\nAll tests executed.\n";
}
