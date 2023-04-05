// #define BLOCKSITE 4
// #define TWISTED_BOUNDARY

#include "fermiqcd.h"

using namespace MDP;

// NOTES:
//   make meson does not work if you enable twist. The rest does work!

int main(int argc, char **argv)
{
    mpi.open_wormholes(argc, argv);
    define_base_matrices("FERMILAB");

    //  define_twist_matrices();

    printf("COMPUTING C2(t) FOR STAGGERED MESON: %s\n", argv[1]);

    float seed = 0;
    int t, Nt = 24, Ns = 8;
    int Nc = 3;
    int box[] = {Nt, Ns, Ns, Ns};

    mdp_real pU, u0 = 0.8629;
    mdp_array<mdp_real, 1> c;
    char s[256];

    mdp_lattice space_time(4, box,
                           default_partitioning<0>,
                           torus_topology,
                           seed, 3);

    int conf, nconf = (int)val(argv[2]);

    gauge_field U(space_time, Nc);
    gauge_field V(space_time, Nc);
    // eta_field eta(space_time);
    mdp_site x(space_time);
    mdp_matrix prop;
    coefficients mass_a;
    coefficients mass_b;

    mass_a["mass"] = 0.03;
    mass_b["mass"] = 0.03;

    default_staggered_action = StaggeredAsqtadActionFast::mul_Q;
    default_staggered_inverter = StaggeredBiCGUML::inverter;

    u0 = 0.8629;

    mdp_jackboot c2(nconf, Nt);

    for (conf = 0; conf < nconf; conf++)
    {
        c2.conf = conf;
        snprintf(s, 256, "/home/mdp/data/gauge_improved_b7.4_u0.8629/gauge32x08_b7.4_u0.8629_n%.6i", conf + 1);

        U.load(s);
        snprintf(s, 256, "/home/mdp/data/gauge_improved_b7.4_u0.8629/gauge24x08_b7.4_u0.8629_n%.6i", conf + 1);
        U.save(s);

        pU = pow(u0, 4);
        c = lepage_coefficients(pU, "Full");
        set_antiperiodic_phases(U, 0, true);
        lepage_improved_links(V, U, c);

        if (strcmp(argv[1], "5x5") == 0)
            prop = make_meson(U, V, Gamma5, Gamma5, mass_a, mass_b,
                              wall_source, wall_source, 1e-5);
        if (strcmp(argv[1], "05x05") == 0)
            prop = make_meson(U, V, Gamma[0] * Gamma5, Gamma[0] * Gamma5, mass_a, mass_b,
                              wall_source, wall_source, 1e-5);
        if (strcmp(argv[1], "5x35") == 0)
            prop = make_meson(U, V, Gamma5, Gamma[3] * Gamma5, mass_a, mass_b,
                              wall_source, wall_source, 1e-5);
        if (strcmp(argv[1], "05x12") == 0)
            prop = make_meson(U, V, Gamma[0] * Gamma5, Gamma[1] * Gamma[2], mass_a, mass_b,
                              wall_source, wall_source, 1e-5);
        if (strcmp(argv[1], "5x03") == 0)
            prop = make_meson(U, V, Gamma5, Gamma[0] * Gamma[3], mass_a, mass_b, wall_source,
                              wall_source, 1e-5);
        if (strcmp(argv[1], "05x3") == 0)
            prop = make_meson(U, V, Gamma[0] * Gamma5, Gamma[3], mass_a, mass_b,
                              wall_source, wall_source, 1e-5);
        if (strcmp(argv[1], "5x0") == 0)
            prop = make_meson(U, V, Gamma5, Gamma[0], mass_a, mass_b,
                              wall_source, wall_source, 1e-5);
        if (strcmp(argv[1], "05xI") == 0)
            prop = make_meson(U, V, Gamma[0] * Gamma5, Gamma1, mass_a, mass_b,
                              wall_source, wall_source, 1e-5);

        for (t = 0; t < Nt; t++)
            c2(t) = real(prop(0, t));

        float x, dx;

        if (ME == 0)
        {
            for (t = 0; t < Nt; t = t + 2)
            {
                c2.plain(t);
                x = c2.mean();
                dx = c2.j_err();
                printf("%i,\t%e,\t%e\n", t, log(x), dx / x);
            }
            printf("\n");
        }
    }
    mpi.close_wormholes();
    return 0;
}