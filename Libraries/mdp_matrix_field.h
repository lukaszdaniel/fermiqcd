/////////////////////////////////////////////////////////////////
/// @file mdp_matrix_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_matrix_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_MATRIX_FIELD_
#define MDP_MATRIX_FIELD_

#include "mdp_nmatrix_field.h"

namespace MDP
{
  /// @brief a field of complex matrices
  ///
  /// Example:
  /// @verbatim
  ///    int box[]={10,10,10};
  ///    mdp_lattice lattice(3,box);
  ///    mdp_matrix_field h(lattice,5,5);
  ///    mdp_site x(lattice);
  ///    forallsites(x)
  ///       h(x)=lattice.random(x).SU(5);
  /// @endverbatim
  class mdp_matrix_field : public mdp_nmatrix_field
  {
  public:
    mdp_matrix_field() : mdp_nmatrix_field()
    {
    }

    /** @brief declares a field of ixj matrices at each site
     */
    mdp_matrix_field(mdp_lattice &a, mdp_uint i, mdp_uint j) : mdp_nmatrix_field(a, 1, i, j)
    {
    }

    mdp_matrix_field(const mdp_matrix_field &field) : mdp_nmatrix_field(field)
    {
    }

    /** @brief dynamically allocates a field of ixj matrices at each site
     */
    void allocate_mdp_matrix_field(mdp_lattice &a, mdp_uint i, mdp_uint j)
    {
      allocate_mdp_nmatrix_field(a, 1, i, j);
    }

    /** @brief returns the matrix stored at site x
     */
    mdp_matrix operator()(mdp_site x)
    {
      return mdp_matrix(address(x), m_rows, m_columns);
    }

    /** @brief returns the (i,j) component of the matrix stored at site x
     */
    mdp_complex &operator()(mdp_site x, mdp_uint i, mdp_uint j)
    {
#ifdef CHECK_BOUNDARY
      if (i >= m_rows)
      {
        error("field rows can be indexed up to " + (m_rows - 1));
      }
      if (j >= m_rows)
      {
        error("field columns can be indexed up to " + (m_columns - 1));
      }
#endif
      return address(x)[i * m_columns + j];
    }

    /** @brief returns the (i,j) const component of the matrix stored at site x
     */
    const mdp_complex &operator()(mdp_site x, mdp_uint i, mdp_uint j) const
    {
#ifdef CHECK_BOUNDARY
      if (i >= m_rows)
      {
        error("field rows can be indexed up to " + (m_rows - 1));
      }
      if (j >= m_rows)
      {
        error("field columns can be indexed up to " + (m_columns - 1));
      }
#endif
      return address(x)[i * m_columns + j];
    }
  };
} // namespace MDP

#endif /* MDP_MATRIX_FIELD_ */
