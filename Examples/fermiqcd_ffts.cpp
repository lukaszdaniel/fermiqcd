#include "fermiqcd.h"

using namespace MDP;

void test_ft()
{
  constexpr mdp_suint n = 16;

  std::vector<mdp_complex> in(n), out(n), chk(n);

  for (mdp_suint i = 0; i < n; i++)
    in[i] = exp(I * 2.0 * sin(i));

  dft(out, in, n, +1);

  dft(chk, out, n, -1);

  mdp_real accumulated_difference = 0.0;
  for (mdp_suint i = 0; i < n; i++)
  {
    mdp_real tmp = abs(chk[i] - in[i]);
    if (tmp > mdp_precision)
      std::cout << "(" << real(in[i]) << ", " << imag(in[i]) << ") "
                << "(" << real(chk[i]) << ", " << imag(chk[i]) << ") "
                << "(" << real(out[i]) << ", " << imag(out[i]) << ")" << std::endl;

    accumulated_difference += tmp;
  }
  std::cout << "Average difference: " << accumulated_difference / n << "\n\n";
}

void fermi_field_ft()
{

  constexpr Box box = {1, 20, 20, 20};
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

  test_ft();
  fermi_field_ft();

  mdp.close_wormholes();

  return 0;
}
