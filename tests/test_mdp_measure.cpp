#include "mdp.h"

using namespace MDP;

void assert_near(double a, double b, double tol, const char *msg)
{
  if (std::fabs(a - b) > tol)
  {
    std::cerr << "FAIL: " << msg << " got=" << a << " expected=" << b << "\n";
    std::exit(1);
  }
}

// --- konstruktor i set/reset
void test_basic()
{
  mdp_measure m(5, 2, 3);
  assert_near(m.getmean(), 5, 1e-6, "mean");
  assert_near(m.getmerr(), 2, 1e-6, "error");
  assert_near(m.getnum(), 3, 1e-6, "num");

  m.reset();
  assert_near(m.getmean(), 0, 1e-6, "reset mean");
  assert_near(m.getmerr(), 0, 1e-6, "reset error");
  assert_near(m.getnum(), 0, 1e-6, "reset num");
}

// --- operator <<
void test_accumulation()
{
  mdp_measure m;
  m << 1.0f;
  m << 3.0f;

  assert_near(m.getmean(), 2.0f, 1e-5, "mean accumulation");
  assert_near(m.getmerr(), 1.0f, 1e-5, "error accumulation");
  assert_near(m.getnum(), 2, 1e-5, "num accumulation");
}

// --- + -
void test_add_sub()
{
  mdp_measure a(2, 0.5);
  mdp_measure b(3, 0.2);

  auto c = a + b;
  assert_near(c.getmean(), 5, 1e-6, "add mean");
  assert_near(c.getmerr(), 0.7, 1e-6, "add error");
  assert_near(c.getnum(), 1, 1e-6, "add num");

  auto d = a - b;
  assert_near(d.getmean(), -1, 1e-6, "sub mean");
  assert_near(d.getnum(), 1, 1e-6, "sub num");
}

// --- * /
void test_mul_div()
{
  mdp_measure a(2, 0.5);
  mdp_measure b(4, 0.2);

  auto c = a * b;
  assert_near(c.getmean(), 8, 1e-6, "mul mean");
  assert_near(c.getmerr(), 2.4, 1e-6, "mul error");
  assert_near(c.getnum(), 1, 1e-6, "mul num");

  auto d = a / b;
  assert_near(d.getmean(), 0.5, 1e-6, "div mean");
  assert_near(d.getmerr(), 0.15, 1e-6, "div error");
  assert_near(d.getnum(), 1, 1e-6, "div num");
}

// --- scalar ops
void test_scalar()
{
  mdp_measure a(2, 0.5);

  auto b = a * (-3);
  assert_near(b.getmean(), -6, 1e-6, "scalar mul mean");
  assert_near(b.getmerr(), 1.5, 1e-6, "scalar mul error");
  assert_near(b.getnum(), 1, 1e-6, "scalar mul num");

  auto c = a + 1;
  assert_near(c.getmean(), 3, 1e-6, "scalar add");
  assert_near(c.getmerr(), 0.5, 1e-6, "scalar add error");
  assert_near(c.getnum(), 1, 1e-6, "scalar add num");
}

// --- funkcje matematyczne
void test_functions()
{
  mdp_measure a(0, 0.1);

  auto s = a.sin();
  assert_near(s.getmean(), 0, 1e-6, "sin mean");
  assert_near(s.getmerr(), 0.1, 1e-6, "sin error");

  mdp_measure b(1, 0.1);
  auto e = b.exp();
  assert_near(e.getmean(), std::exp(1.0f), 1e-6, "exp mean");
  assert_near(e.getmerr(), 0.271828, 1e-6, "exp mean error");
}

// --- log
void test_log()
{
  mdp_measure a(2, 0.2);
  auto l = a.log();

  assert_near(l.getmean(), std::log(2.0f), 1e-6, "log mean");
  assert_near(l.getmerr(), 0.2f / 2.0f, 1e-6, "log error");
}

// --- operator >>
void test_sampling()
{
  mdp_measure a(5, 2);
  float x;
  a >> x;

  assert_near(x, 6.40774, 1e-6, "sampling (gaussian)");
}

// --- edge: division by small number
void test_edge_division()
{
  mdp_measure a(1, 0.1);
  mdp_measure b(1e-6, 0.1);

  auto c = a / b;
  assert(std::isfinite(c.getmean()));
  assert_near(c.getmean(), 1e6, 1e-6, "div mean");
  assert_near(c.getmerr(), 100000104448, 1e-6, "div error");
}

int main()
{
  std::cout << "Running mdp_measure tests...\n";

  test_basic();
  std::cout << "Testing basic\n";

  test_accumulation();
  std::cout << "Testing accumulation\n";

  test_add_sub();
  std::cout << "Testing add/sub\n";

  test_mul_div();
  std::cout << "Testing mul/div\n";

  test_scalar();
  std::cout << "Testing scalar\n";

  test_functions();
  std::cout << "Testing functions\n";

  test_log();
  std::cout << "Testing log\n";

  test_sampling();
  std::cout << "Testing sampling\n";

  test_edge_division();
  std::cout << "Testing edge cases\n";

  std::cout << "\nAll tests passed!\n";

  std::cout << "\nOther test output:\n";
  mdp_measure m;
  // store 10 measurements
  for (int i = 0; i < 10; i++)
    m << 3.0 + mdp_random.gaussian(2.0);

  std::cout << m.getmean() << " +/- " << m.getmerr() << std::endl;
  std::cout << m.getnum() << std::endl;

  m = sin(exp(m) + m);

  std::cout << m.getmean() << " +/- " << m.getmerr() << std::endl;
  std::cout << m.getnum() << std::endl;

  return 0;
}
