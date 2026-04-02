#include "mdp.h"

using namespace MDP;

// ===== Kolmogorov-Smirnov =====
double kolmogorov_smirnov_stat(std::vector<double> &data)
{
  std::sort(data.begin(), data.end());
  int n = data.size();

  double D_plus = 0.0;
  double D_minus = 0.0;

  for (int i = 0; i < n; ++i)
  {
    double Fi = (i + 1) / (double)n;
    double Dp = Fi - data[i];
    double Dm = data[i] - (i / (double)n);

    D_plus = std::max(D_plus, Dp);
    D_minus = std::max(D_minus, Dm);
  }

  return std::max(D_plus, D_minus);
}

// Approximated p-value
double ks_pvalue(double D, int n)
{
  double lambda = (std::sqrt(n) + 0.12 + 0.11 / std::sqrt(n)) * D;

  double sum = 0.0;
  for (int k = 1; k <= 100; ++k)
  {
    double term = std::exp(-2.0 * k * k * lambda * lambda);
    sum += (k % 2 == 1 ? 1 : -1) * term;
  }
  return 2.0 * sum;
}

// ===== Autocorrelation =====
double autocorrelation(const std::vector<double> &data, int lag)
{
  int n = data.size();

  double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;

  double num = 0.0;
  double denom = 0.0;

  for (int i = 0; i < n; ++i)
  {
    double diff = data[i] - mean;
    denom += diff * diff;
  }

  for (int i = 0; i < n - lag; ++i)
  {
    num += (data[i] - mean) * (data[i + lag] - mean);
  }

  return num / denom;
}

// ===== FUNKCJA TESTUJĄCA =====
template <typename RNG>
void run_tests(RNG &rng, const std::string &name, int N, int max_lag)
{
  std::vector<double> data(N);

  for (int i = 0; i < N; ++i)
  {
    data[i] = rng.plain();
  }

  std::cout << "\n=============================\n";
  std::cout << "Generator: " << name << "\n";

  std::vector<double> data_copy = data;

  double D = kolmogorov_smirnov_stat(data_copy);
  double p = ks_pvalue(D, N);

  std::cout << "\n[Kolmogorov-Smirnov]\n";
  std::cout << "D = " << D << "\n";
  std::cout << "p-value = " << p << "\n";

  if (p < 0.05)
    std::cout << "Result: REJECT (not uniform)\n";
  else
    std::cout << "Result: OK\n";

  std::cout << "\n[Autocorrelation]\n";

  double threshold = 2.0 / std::sqrt(N);

  for (int lag = 1; lag <= max_lag; ++lag)
  {
    double r = autocorrelation(data, lag);

    std::cout << "lag " << lag
              << ": r = " << r
              << " (+/-" << threshold << ")";

    if (std::abs(r) > threshold)
      std::cout << "  <-- suspicious";

    std::cout << "\n";
  }
}

int main()
{
  constexpr int N = 800000;
  constexpr int max_lag = 10;

  mdp_prng prng_random;
  mdp_prng_sfmt sfmt_random;

  prng_random.initialize(12345);
  sfmt_random.initialize(12345);

  run_tests(prng_random, "mdp_prng", N, max_lag);
  run_tests(sfmt_random, "mdp_prng_sfmt", N, max_lag);

  std::cout << "\nFirst 10 random numbers from each generator:\n";
  std::cout << "mdp_prng\tmdp_prng_sfmt\n";
  for (int n = 0; n < 10; n++)
  {
    std::cout << std::fixed << std::setprecision(6) << prng_random.plain() << "\t" << sfmt_random.plain() << "\n";
  }
  return 0;
}
