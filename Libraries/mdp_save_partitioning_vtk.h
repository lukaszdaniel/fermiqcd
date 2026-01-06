/////////////////////////////////////////////////////////////////
/// @file mdp_save_partitioning_vtk.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains a method to save partitioning to a VTK file
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_SAVE_PARTITIONING_VTK_
#define MDP_SAVE_PARTITIONING_VTK_


#include "mdp_lattice.h"
#include "mdp_field.h"

namespace MDP
{
  void save_partitioning_vtk(mdp_lattice &lattice, const std::string &filename)
  {
    mdp_field<int> where(lattice);
    mdp_site x(lattice);
    forallsites(x)
    {
      where(x) = ME;
    }
    where.save_vtk(filename);
  }
} // namespace MDP

#endif /* MDP_SAVE_PARTITIONING_VTK_ */
