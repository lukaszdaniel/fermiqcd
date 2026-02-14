#include "fermiqcd.h"

using namespace MDP;

int test1()
{
    mdp_int n = 16;

    double tmp = 0;

    std::vector<mdp_complex> in(n), out(n), chk(n);
    for (mdp_int i = 0; i < n; i++)
        in[i] = exp(I * 2.0 * sin(i));

    dft(out.data(), in.data(), n, +1);

    dft(chk.data(), out.data(), n, -1);

    for (mdp_int i = 0; i < n; i++)
    {
        if (tmp += abs(chk[i] - in[i]) > 0.00001)
            printf("(%f, %f) (%f, %f) (%f, %f)\n",
                   real(in[i]), imag(in[i]),
                   real(chk[i]), imag(chk[i]),
                   real(out[i]), imag(out[i]));
    }
    printf("%f\n", tmp / n);
    return 0;
};

void test2()
{

    constexpr Box box = {1, 20, 20, 20};
    mdp_lattice lattice(box);
    fermi_field psi(lattice, 3);
    fermi_field phi(lattice, 3);
    fermi_field chi(lattice, 3);

    set_random(psi);
    mdp_site x(lattice);

    for (int t = 0; t < box[0]; t++)
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
    test1();
    test2();
    mdp.close_wormholes();
    return 0;
}
