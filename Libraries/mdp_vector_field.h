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

namespace MDP
{
  /// @brief a field of vectors of complex numbers
  ///
  /// Example:
  /// @verbatim
  ///    int box[]={10,10,10};
  ///    mdp_lattice lattice(3,box);
  ///    mdp_vector_field h(lattice,10);
  ///    mdp_site x(lattice);
  ///    forallsites(x)
  ///      h(x)=0.0+0.0*I;
  /// @endverbatim
  class mdp_vector_field : public mdp_field<mdp_complex>
  {
  private:
    mdp_int m_rows;
    mdp_int m_columns;
    mdp_int m_imax;

  public:
    mdp_vector_field() : mdp_field<mdp_complex>(), m_rows(0), m_columns(0), m_imax(0)
    {
    }

    mdp_vector_field(mdp_lattice &a, int i) : mdp_field<mdp_complex>(a, i), m_rows(i), m_columns(1), m_imax(i)
    {
    }

    mdp_vector_field(mdp_vector_field &field) : mdp_field<mdp_complex>(field), m_rows(field.m_rows), m_columns(field.m_columns), m_imax(field.m_imax)
    {
    }

    void allocate_mdp_vector_field(mdp_lattice &a, int i)
    {
      deallocate_field();
      m_rows = i;
      m_columns = 1;
      m_imax = i;
      allocate_field(a, m_imax);
    }

    mdp_matrix operator()(mdp_site x)
    {
      return mdp_matrix(address(x), m_rows, m_columns);
    }

    mdp_complex &operator()(mdp_site x, int i)
    {
      return address(x)[i];
    }

    const mdp_complex &operator()(mdp_site x, int i) const
    {
      return address(x)[i];
    }
  };
} // namespace MDP

#endif /* MDP_VECTOR_FIELD_ */
