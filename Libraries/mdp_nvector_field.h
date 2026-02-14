/////////////////////////////////////////////////////////////////
/// @file mdp_nvector_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_nvector_field (deprecated)
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_NVECTOR_FIELD_
#define MDP_NVECTOR_FIELD_

#include "mdp_nmatrix_field.h"

namespace MDP
{
  /// @brief a field of vectors of vectors (DEPRECATED)
  ///
  /// Example:
  /// @verbatim
  ///    constexpr Box box = {10,10,10};
  ///    mdp_lattice lattice(box);
  ///    mdp_nvector_field h(lattice,10,2);
  ///    mdp_site x(lattice);
  ///    forallsites(x)
  ///      for(int i=0; i<10; i++)
  ///        h(x,i)=lattice.random(x).SU(3);
  /// @endverbatim
  class mdp_nvector_field : public mdp_nmatrix_field
  {
  public:
    mdp_nvector_field() : mdp_nmatrix_field()
    {
    }

    /** @brief declares a n-component vector field of i-component vectors at each site
     */
    mdp_nvector_field(mdp_lattice &a, mdp_uint n, mdp_uint i) : mdp_nmatrix_field(a, n, i, 1)
    {
    }

    mdp_nvector_field(const mdp_nvector_field &field) : mdp_nmatrix_field(field)
    {
    }

    /** @brief dynamically allocates a n-component vector field of i-component vectors at each site
     */
    void allocate_mdp_nvector_field(mdp_lattice &a, mdp_uint n, mdp_uint i)
    {
      allocate_mdp_nmatrix_field(a, n, i, 1);
    }

    /** @brief returns the n-th vector stored at site x
     */
    mdp_matrix operator()(mdp_site x, mdp_uint n)
    {
#ifdef CHECK_BOUNDARY
      if (n >= m_matrices)
      {
        error("field component can be indexed up to " + (m_matrices - 1));
      }
#endif
      return mdp_matrix(address(x, n * m_rows * m_columns), m_rows, m_columns);
    }

    /** @brief returns the i-th component of the n-th vector stored at site x
     */
    mdp_complex &operator()(mdp_site x, mdp_uint n, mdp_uint i)
    {
#ifdef CHECK_BOUNDARY
      if (n >= m_matrices)
      {
        error("field component can be indexed up to " + (m_matrices - 1));
      }
      if (i >= m_rows)
      {
        error("field rows can be indexed up to " + (m_rows - 1));
      }
#endif
      return address(x, n * m_rows * m_columns)[i * m_columns + 0];
    }

    /** @brief returns the i-th const component of the n-th vector stored at site x
     */
    const mdp_complex &operator()(mdp_site x, mdp_uint n, mdp_uint i) const
    {
#ifdef CHECK_BOUNDARY
      if (n >= m_matrices)
      {
        error("field component can be indexed up to " + (m_matrices - 1));
      }
      if (i >= m_rows)
      {
        error("field rows can be indexed up to " + (m_rows - 1));
      }
#endif
      return address(x, n * m_rows * m_columns)[i * m_columns + 0];
    }
  };
} // namespace MDP

#endif /* MDP_NVECTOR_FIELD_ */
