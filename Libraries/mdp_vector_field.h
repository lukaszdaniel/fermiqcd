/////////////////////////////////////////////////////////////////
/// @file mdp_vector_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_vector_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_VECTOR_FIELD_
#define MDP_VECTOR_FIELD_

#include "mdp_nvector_field.h"

namespace MDP
{
  /// @brief a field of vectors of complex numbers
  ///
  /// Example:
  /// @verbatim
  ///    constexpr Box box = {10,10,10};
  ///    mdp_lattice lattice(box);
  ///    mdp_vector_field h(lattice,10);
  ///    mdp_site x(lattice);
  ///    forallsites(x)
  ///      h(x)=0.0+0.0*I;
  /// @endverbatim
  class mdp_vector_field : public mdp_nvector_field
  {

  public:
    mdp_vector_field() : mdp_nvector_field()
    {
    }

    /** @brief declares a field of i-component vectors at each site
     */
    mdp_vector_field(mdp_lattice &a, int i) : mdp_nvector_field(a, 1, i)
    {
    }

    mdp_vector_field(const mdp_vector_field &field) : mdp_nvector_field(field)
    {
    }

    /** @brief dynamically allocates a field of i-component vectors at each site
     */
    void allocate_mdp_vector_field(mdp_lattice &a, int i)
    {
      allocate_mdp_nvector_field(a, 1, i);
    }

    /** @brief returns the matrix stored at site x
     */
    mdp_matrix operator()(mdp_site x)
    {
      return mdp_matrix(address(x), m_rows, m_columns);
    }

    /** @brief returns the i-th component of the vector stored at site x
     */
    mdp_complex &operator()(mdp_site x, mdp_uint i)
    {
#ifdef CHECK_BOUNDARY
      if (i >= m_rows)
      {
        error("field rows can be indexed up to " + (m_rows - 1));
      }
#endif
      return address(x)[i * m_columns + 0];
    }

    /** @brief returns the i-th const component of the vector stored at site x
     */
    const mdp_complex &operator()(mdp_site x, mdp_uint i) const
    {
#ifdef CHECK_BOUNDARY
      if (i >= m_rows)
      {
        error("field rows can be indexed up to " + (m_rows - 1));
      }
#endif
      return address(x)[i * m_columns + 0];
    }
  };
} // namespace MDP

#endif /* MDP_VECTOR_FIELD_ */
