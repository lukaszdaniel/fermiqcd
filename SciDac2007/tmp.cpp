/*
    python3 iqcd.py code +cold +save_partitioning

    +cold warning: assuming default argument TxXxYxZ=16x4x4x4
    +cold warning: assuming default argument nc=3
    +save_partitioning warning: assuming default argument filename=partitioning
*/
#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
    mdp.open_wormholes(argc, argv);
    std::string filename;
    coefficients coeff;
    constexpr Box L = {16, 4, 4, 4};
    mdp_lattice spacetime(L);
    int nc = 3;
    gauge_field U(spacetime, nc);
    set_cold(U);

    save_partitioning_vtk(spacetime, "partitioning");

    mdp.close_wormholes();
    return 0;
}
