/*
    python qcd.py make +cold +save_partitioning_vtk

    +cold warning: assuming default argument TxXxYxZ=16x4x4x4
    +cold warning: assuming default argument nc=3
    +save_partitioning_vtk warning: assuming default argument filename=partitining
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
