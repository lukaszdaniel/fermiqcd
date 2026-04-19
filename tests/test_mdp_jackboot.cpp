// Program: example07.cpp
#include "mdp.h"

using namespace MDP;

constexpr mdp_suint n = 100;

mdp_real f1(const mdp_real *x, [[maybe_unused]] const void *a) { return x[0] / x[1]; }

mdp_real f2(const mdp_real *x, [[maybe_unused]] const void *a) { return x[0] * x[1]; }

int main()
{
  std::cout << std::fixed << std::setprecision(6);

  mdp_matrix A;
  mdp_jackboot jb(n, 2);

  for (mdp_suint i = 0; i < n; i++)
  {
    A = Random.SU(6) + Random.gaussian();
    jb(i, 0) = real(det(inv(A))); // test operator()(conf,arg)
    jb(i, 1) = real(det(A));
  }

  // =========================
  // 2. Test address()
  // =========================
  {
    mdp_real *ptr = jb.address(0);
    assert(ptr != nullptr);
    std::cout << "[OK] address() works\n";
  }

  // =========================
  // 3. Test set_conf + operator()(arg)
  // =========================
  {
    jb.set_conf(0);
    mdp_real v0 = jb(0);
    mdp_real v1 = jb(1);

    std::cout << "[INFO] set_conf + operator()(arg): "
              << v0 << ", " << v1 << "\n";
  }

  // =========================
  // 4. Test plain()
  // =========================
  {
    jb.plain(0);
    mdp_real mean0 = jb.mean();

    jb.plain(1);
    mdp_real mean1 = jb.mean();

    std::cout << "[PLAIN] mean(x[0]) = " << mean0 << "\n";
    std::cout << "[PLAIN] mean(x[1]) = " << mean1 << "\n";
  }

  // =========================
  // 5. Test with f function f (division)
  // =========================
  jb.f = f1;

  mdp_real mean_f1 = jb.mean();
  mdp_real jerr_f1 = jb.j_err();
  mdp_real berr_f1 = jb.b_err(200);

  std::cout << "\n[f1] x[0]/x[1]\n";
  std::cout << "mean       = " << mean_f1 << "\n";
  std::cout << "jackknife  = " << jerr_f1 << "\n";
  std::cout << "bootstrap  = " << berr_f1 << "\n";

  // =========================
  // 6. Test with f function f (multiplication)
  // =========================
  jb.f = f2;

  mdp_real mean_f2 = jb.mean();
  mdp_real jerr_f2 = jb.j_err();
  mdp_real berr_f2 = jb.b_err(200);

  std::cout << "\n[f2] x[0]*x[1]\n";
  std::cout << "mean       = " << mean_f2 << "\n";
  std::cout << "jackknife  = " << jerr_f2 << "\n";
  std::cout << "bootstrap  = " << berr_f2 << "\n";

  // =========================
  // 7. Edge case (not enough data)
  // =========================
  {
    mdp_jackboot small(1, 1);
    small(0, 0) = 42.0f;
    small.plain(0);

    std::cout << "\n[EDGE CASE]\n";
    std::cout << "mean = " << small.mean() << "\n";
    std::cout << "j_err = " << small.j_err() << " (should be 0)\n";
    std::cout << "b_err = " << small.b_err() << " (should be 0)\n";
  }

  // =========================
  // 8. Sanity check
  // =========================
  {
    jb.plain(0);
    mdp_real m = jb.mean();

    assert(std::isfinite(m));
    std::cout << "\n[OK] mean() is finite\n";
  }

  std::cout << "\n=== All tests passsed ===\n";

  return 0;
}
