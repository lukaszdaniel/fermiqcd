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

float quad(float *x, mdp_int, void *)
{
  return (*x - 2.0f) * (*x - 2.0f);
}

float linear_model(float x, float *a, mdp_int, void *)
{
  return a[0] * x + a[1];
}

void test_linear_fit()
{
  float x[] = {0, 1, 2, 3};
  mdp_measure y[4];

  for (int i = 0; i < 4; i++)
  {
    y[i].setmean(2 * x[i] + 1);
    y[i].setmerror(0.1);
  }

  mdp_measure a[2];
  linear_fit(x, y, 0, 4, a);

  assert_near(a[0].getmean(), 2.0, 1e-5, "linear slope");
  assert_near(a[1].getmean(), 1.0, 1e-5, "linear intercept");
}

void test_golden()
{
  float xmin;
  float f = golden_rule(quad, xmin, 0, 2, 5);

  assert_near(xmin, 2.0, 1e-2, "golden xmin");
  assert_near(f, 0.0, 1e-3, "golden fmin");
}

void test_blm()
{
  float x[] = {0, 1, 2};
  mdp_measure y[3];

  for (int i = 0; i < 3; i++)
  {
    y[i].set(2 * x[i] + 1, 0.1, i);
  }

  float a[] = {0, 0};
  mdp_matrix cov(2, 2);

  BayesianLevenbergMarquardt(x, y, 0, 3, a, 2, cov, linear_model);

  assert_near(a[0], 2.0, 1e-1, "BLM slope");
  assert_near(a[1], 1.0, 1e-1, "BLM intercept");
}

int main()
{
  std::cout << "Running tests...\n";

  test_linear_fit();
  std::cout << "- linear_fit passed\n";

  test_golden();
  std::cout << "- golden_rule passed\n";

  test_blm();
  std::cout << "- BLM passed\n";

  std::cout << "All tests passed!\n";
  return 0;
}
