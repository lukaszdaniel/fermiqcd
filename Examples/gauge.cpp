#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
    mdp.open_wormholes(argc, argv);

    bool hot = false;
    bool save_output = false;
    bool load_field = false;
    char input[1024] = "";
    char output[1024] = "";
    int nc = 3;
    mdp_real beta = 6.0;
    int nsteps = 0;
    int *L = nullptr;
    mdp_field_file_header header;

    // //////////////////////////////
    // Parsing command line arguments
    // //////////////////////////////
    for (int i = 1; i < argc; i++)
    {
        if (strncmp(argv[i], "-hot", 4) == 0)
        {
            hot = true;
        }
        else if (strncmp(argv[i], "-input", 6) == 0)
        {
            sscanf(argv[i], "-input=%s", input);
            load_field = true;
        }
        else if (strncmp(argv[i], "-output", 6) == 0)
        {
            sscanf(argv[i], "-output=%s", output);
            save_output = true;
        }
        else if (strncmp(argv[i], "-L", 2) == 0)
        {
            int box[4];
            sscanf(argv[i], "-L=%ix%ix%ix%i", box, box + 1, box + 2, box + 3);
            L = box;
        }
        else if (strncmp(argv[i], "-steps", 6) == 0)
        {
            sscanf(argv[i], "-steps=%i", &nsteps);
        }
        else if (strncmp(argv[i], "-nc", 3) == 0)
        {
            sscanf(argv[i], "-nc=%i", &nc);
        }
        else if (strncmp(argv[i], "-beta", 5) == 0)
        {
#ifndef USE_DOUBLE_PRECISION
            sscanf(argv[i], "-beta=%f", &beta);
#else
            sscanf(argv[i], "-beta=%lf", &beta);
#endif
        }
        else
            error("Wrong command line option");
    }

    if (load_field)
    {
        if (file_exists(input))
            header = get_info(input);

        if (header.ndim != 4)
            error("Sorry, mesons only in 4D");

        nc = (int)sqrt((double)header.bytes_per_site / (4 * sizeof(mdp_complex)));
        L = header.box;
    }

    if (!L)
        error("Lattice is not set");

    mdp_lattice lattice(4, L);

    gauge_field U(lattice, nc);
    coefficients coeff;
    coeff["beta"] = beta;

    if (load_field)
    {
        U.load(std::string(input));
    }
    else
    {
        if (hot)
        {
            set_hot(U);
        }
        else
        {
            set_cold(U);
        }
    }

    if (nsteps > 0)
    {
        ImprovedGaugeAction::heatbath(U, coeff, nsteps, "MILC");
    }

    if (save_output)
    {
        U.save(std::string(output));
    }

    mdp.close_wormholes();
    return 0;
}
