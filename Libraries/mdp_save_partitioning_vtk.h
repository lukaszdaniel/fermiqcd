#ifndef MDP_SAVE_PARTITIONING_VTK_
#define MDP_SAVE_PARTITIONING_VTK_

void save_partitioning_vtk(mdp_lattice &lattice, std::string filename)
{
  mdp_field<int> where(lattice);
  mdp_site x(lattice);
  forallsites(x)
  {
    where(x) = ME;
  }
  where.save_vtk(filename);
}

#endif /* MDP_SAVE_PARTITIONING_VTK_ */
