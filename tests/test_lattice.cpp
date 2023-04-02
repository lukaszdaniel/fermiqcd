#include "mdp.h"

/** @brief Test lattice volumes
 *
 * Sum of local volumes should add up to global volume.
 *
 * When running the program with more than 1 process,
 * the sum of enclosing volumes should be greater than the global volume
 * since each enclosing volume stores, apart from local sites,
 * copies of the neighbouring local lattices.
 *
 * @note For parallel processing run for example:
 * ./test_lattice.exe -PSIM_NPROCS=4
 */
int main(int argc, char **argv)
{
    mdp.open_wormholes(argc, argv);
    int mybox[] = {8, 8, 8};
    mdp_lattice mylattice(3, mybox);
    mdp_site x(mylattice);
    int sum_lattice_sites = 0;
    mdp_int global_vol = mylattice.global_volume();
    mdp_int enclosing_vol = mylattice.nvol;
    mdp_int local_vol = mylattice.local_volume();

    forallsites(x)
    {
        ++sum_lattice_sites;
    }
    mdp.add(sum_lattice_sites);
    mdp.add(enclosing_vol);
    mdp.add(local_vol);

    std::cout << "Number of lattice sites: " << sum_lattice_sites << "\n";
    std::cout << "Lattice global volume: " << global_vol << "\n";
    std::cout << "Lattice enclosing volume: " << enclosing_vol << "\n";
    std::cout << "Lattice local volume: " << local_vol << "\n";

    mdp.close_wormholes();
}
