#ifndef mdp_save_partitioning_vtk_
#define mdp_save_partitioning_vtk_

using namespace std;

void save_partitioning_vtk(mdp_lattice &lattice, string filename) {
	mdp_field<int> where(lattice);
	mdp_site x(lattice);
	forallsites(x)
	{
		where(x) = ME;
	}
	where.save_vtk(filename);
}

#endif /* mdp_save_partitioning_vtk_ */
