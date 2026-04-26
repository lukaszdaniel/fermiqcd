#include "fermiqcd.h"

using namespace MDP;

void test_ft(const mdp_suint n, const std::string &&type)
{
  std::vector<mdp_complex> in(n);  // input field
  std::vector<mdp_complex> out(n); // fourier transform of the input field
  std::vector<mdp_complex> chk(n); // inverse fourier transform (should be equal to input)

  for (mdp_suint i = 0; i < n; i++)
  {
    in[i] = exp(I * 1.37 * sin(i) + 0.287);
  }

  if (type == "dft")
  {
    dft(out, in, n, +1);  // fourier transform
    dft(chk, out, n, -1); // inverse fourier transform
  }
  else if (type == "fft")
  {
    fft(out, in, n, +1);  // fourier transform
    fft(chk, out, n, -1); // inverse fourier transform
  }
  else if (type == "ft_auto")
  {
    ft_auto(out, in, n, +1);  // fourier transform
    ft_auto(chk, out, n, -1); // inverse fourier transform
  }
  else
    error("Unknown type of Fourier transform: " + type);

  mdp_real accumulated_difference = 0.0;
  for (mdp_suint i = 0; i < n; i++)
  {
    mdp_real tmp = abs(chk[i] - in[i]);
    std::cout << "(" << real(in[i]) << ", " << imag(in[i]) << ")\n"
              << "(" << real(chk[i]) << ", " << imag(chk[i]) << ")" << (tmp > mdp_precision ? " difference at index " + std::to_string(i) : "") << "\n"
              << "(" << real(out[i]) << ", " << imag(out[i]) << ")\n\n";

    accumulated_difference += tmp;
  }
  std::cout << "Average difference: " << accumulated_difference / n << "\n\n";
}

void fermi_field_ft()
{

  constexpr Box box = {1, 16, 20, 32};
  mdp_lattice lattice(box);
  fermi_field psi(lattice, 3);
  fermi_field phi(lattice, 3);
  fermi_field chi(lattice, 3);

  set_random(psi);
  mdp_site x(lattice);

  for (mdp_uint t = 0; t < box[0]; t++)
  {
    if (on_which_process(lattice, t) == ME)
    {
      fermi_field_fft(t, phi, psi, 1);
      fermi_field_fft(t, chi, phi, -1);
    }
  }

  check_differences(psi, chi);
}

int main(int argc, char **argv)
{
  mdp.open_wormholes(argc, argv);

  std::cout << "Testing dft for n = 12:\n";
  test_ft(12, "dft");
  std::cout << "Testing dft for n = 32:\n";
  test_ft(32, "dft");
  std::cout << "Testing fft for n = 32:\n";
  test_ft(32, "fft");
  std::cout << "Testing ft_auto for n = 14:\n";
  test_ft(14, "ft_auto");
  std::cout << "Testing ft_auto for n = 16:\n";
  test_ft(16, "ft_auto");
  std::cout << "Testing fermi_field_fft:\n";

  fermi_field_ft();

  mdp.close_wormholes();

  return 0;
}
