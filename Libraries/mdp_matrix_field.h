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

namespace MDP
{
  /// @brief a field of matrices
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
  class mdp_matrix_field : public mdp_field<mdp_complex>
  {
  private:
    mdp_uint m_rows;
    mdp_uint m_columns;
    mdp_uint m_imax;

  public:
    mdp_matrix_field() : mdp_field<mdp_complex>(), m_rows(0), m_columns(0), m_imax(0)
    {
    }

    mdp_matrix_field(mdp_matrix_field &field) : mdp_field<mdp_complex>(field), m_rows(field.m_rows), m_columns(field.m_columns), m_imax(field.m_imax)
    {
    }

    /** @brief declares a field of ixj matrices at each site
     */
    mdp_matrix_field(mdp_lattice &a, mdp_uint i, mdp_uint j) : mdp_field<mdp_complex>(a, i * j), m_rows(i), m_columns(j), m_imax(i * j)
    {
    }

    /** @brief dynamically allocates a field of ixj matrices at each site
     */
    void allocate_mdp_matrix_field(mdp_lattice &a, mdp_uint i, mdp_uint j)
    {
      deallocate_field();
      m_rows = i;
      m_columns = j;
      m_imax = i * j;
      allocate_field(a, m_imax);
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
      return address(x)[i * m_columns + j];
    }

    /** @brief returns the (i,j) const component of the matrix stored at site x
     */
    const mdp_complex &operator()(mdp_site x, mdp_uint i, mdp_uint j) const
    {
      return address(x)[i * m_columns + j];
    }
  };
} // namespace MDP

#endif /* MDP_MATRIX_FIELD_ */
