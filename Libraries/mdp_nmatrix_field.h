/////////////////////////////////////////////////////////////////
/// @file mdp_nmatrix_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_nmatrix_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_NMATRIX_FIELD_
#define MDP_NMATRIX_FIELD_

namespace MDP
{
  /// @brief a field of vectors of complex matrices
  ///
  /// Example:
  /// @verbatim
  ///    int box[]={10,10,10};
  ///    mdp_lattice lattice(3,box);
  ///    mdp_nmatrix_field h(lattice,10,3,3);
  ///    mdp_site x(lattice);
  ///    forallsites(x)
  ///      for(int i=0; i<10; i++)
  ///        h(x,i)=lattice.random(x).SU(3);
  /// @endverbatim
  class mdp_nmatrix_field : public mdp_field<mdp_complex>
  {
  protected:
    mdp_uint m_matrices;
    mdp_uint m_rows;
    mdp_uint m_columns;

  public:
    mdp_nmatrix_field() : mdp_field<mdp_complex>(), m_matrices(0), m_rows(0), m_columns(0)
    {
    }

    /** @brief declares a n-component vector field of ixj matrices at each site
     */
    mdp_nmatrix_field(mdp_lattice &a, mdp_uint n, mdp_uint i, mdp_uint j) : mdp_field<mdp_complex>(a, n * i * j), m_matrices(n), m_rows(i), m_columns(j)
    {
    }

    mdp_nmatrix_field(const mdp_nmatrix_field &field) : mdp_field<mdp_complex>(field), m_matrices(field.m_matrices), m_rows(field.m_rows), m_columns(field.m_columns)
    {
    }

    /** @brief dynamically allocates a n-component vector field of ixj matrices at each site
     */
    void allocate_mdp_nmatrix_field(mdp_lattice &a, mdp_uint n, mdp_uint i, mdp_uint j)
    {
      m_matrices = n;
      m_rows = i;
      m_columns = j;
      allocate_field(a, m_matrices * m_rows * m_columns);
    }

    /** @brief returns the n-th matrix stored at site x
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

    /** @brief returns the (i,j) component of the n-th matrix stored at site x
     */
    mdp_complex &operator()(mdp_site x, mdp_uint n, mdp_uint i, mdp_uint j)
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
      if (j >= m_rows)
      {
        error("field columns can be indexed up to " + (m_columns - 1));
      }
#endif
      return address(x, n * m_rows * m_columns)[i * m_columns + j];
    }

    /** @brief returns the (i,j) const component of the n-th matrix stored at site x
     */
    const mdp_complex &operator()(mdp_site x, mdp_uint n, mdp_uint i, mdp_uint j) const
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
      if (j >= m_rows)
      {
        error("field columns can be indexed up to " + (m_columns - 1));
      }
#endif
      return address(x, n * m_rows * m_columns)[i * m_columns + j];
    }
  };
} // namespace MDP

#endif /* MDP_NMATRIX_FIELD_ */
