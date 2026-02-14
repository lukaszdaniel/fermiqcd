#include "fermiqcd.h"

using namespace MDP;

int main(int argc, char **argv)
{
    mdp.open_wormholes(argc, argv);
    define_base_matrices("FERMILAB");
    constexpr Box L = {8, 32, 32, 32};
    mdp_lattice lattice(L);
    gauge_field U(lattice, 3);
    InstantonGenerator4D generator;
    for (int k = 0; k < 100; k++)
    {
        std::cout << k << std::endl;
        std::vector<SingleInstanton4D> instantons;
        instantons.push_back(SingleInstanton4D(7.5, 15.5, 15.5, 10 + 0.1 * k, 1, +1)); // instanton
        instantons.push_back(SingleInstanton4D(7.5, 15.5, 16.5, 20 - 0.1 * k, 1, -1)); // anti-instanton
        generator.generate(U, instantons);
        check_unitarity(U);
        mdp_site x(lattice);
        x.set(0, 0, 0, 0);
        std::cout << U(x, 0) << std::endl;
        topological_charge_vtk(U, std::string("engineered.topological.") + std::to_string(k) + ".vtk", 0);
    }
    mdp.close_wormholes();
    return 0;
}
