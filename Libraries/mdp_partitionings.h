/////////////////////////////////////////////////////////////////
/// @file mdp_partitionings.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Example functions to do parallel partitioning of a lattice
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_PARTITIONINGS_
#define MDP_PARTITIONINGS_

#include <cmath>

namespace MDP
{
  /** @brief Partition function
   *
   * @param x Point x[] on the lattice
   * @param ndim Dimension of the box
   * @param nx Box containing the lattice
   *
   * @return Suggested process ID for point x[]
   *
   * @note Shorthand for generic default_partitioning<0>
   * function
   */
  int default_partitioning0(int *x,
                            int ndim,
                            int *nx)
  {
    default_partitioning<0>(x, ndim, nx);
  }

  /** @brief Partition function
   *
   * Function to calculate suggested process ID
   * based on the location of point x[], lattice size
   * and the number of available processes.
   *
   * Roughly speaking lattice of size S in dim
   * direction will be split into Nproc chunks.
   * Each chunk will belong to one process ID.
   *
   * @param x Point x[] on the lattice
   * @param ndim Dimension of the box
   * @param nx Box containing the lattice
   *
   * @return Suggested process ID for point x[]
   */
  template <int dim>
  int default_partitioning(int *x,
                           int ndim,
                           int *nx)
  {
    int partition_factor = std::ceil((1.0 * nx[dim]) / Nproc);
    return x[dim] / partition_factor;
  }
} // namespace MDP

#endif /* MDP_PARTITIONINGS_ */
