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
  /// @brief field of vectors of matrices
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
  private:
    mdp_uint m_matrices;
    mdp_uint m_rows;
    mdp_uint m_columns;
    mdp_uint m_imax;
    mdp_uint m_imax2;

  public:
    mdp_nmatrix_field() : mdp_field<mdp_complex>(), m_matrices(0), m_rows(0), m_columns(0), m_imax(0), m_imax2(0)
    {
    }

    /** @brief declares a n-component vector field of ixj matrices at each site
     */
    mdp_nmatrix_field(mdp_lattice &a, mdp_uint n, mdp_uint i, mdp_uint j) : mdp_field<mdp_complex>(a, i * j * n), m_matrices(n), m_rows(i), m_columns(j), m_imax(i * j * n), m_imax2(i * j)
    {
    }

    mdp_nmatrix_field(const mdp_nmatrix_field &field) : mdp_field<mdp_complex>(field), m_matrices(field.m_matrices), m_rows(field.m_rows), m_columns(field.m_columns), m_imax(field.m_imax), m_imax2(field.m_imax2)
    {
    }

    /** @brief dynamically allocates a n-component vector field of ixj matrices at each site
     */
    void allocate_mdp_nmatrix_field(mdp_lattice &a, mdp_uint n, mdp_uint i, mdp_uint j)
    {
      deallocate_field();
      m_rows = i;
      m_columns = j;
      m_matrices = n;
      m_imax = i * j * n;
      m_imax2 = i * j;
      allocate_field(a, m_imax);
    }

    /** @brief returns the n-th matrix stored at site x
     */
    mdp_matrix operator()(mdp_site x, mdp_uint n)
    {
      return mdp_matrix(address(x, n * m_imax2), m_rows, m_columns);
    }

    /** @brief returns the (i,j) component of the n-th matrix stored at site x
     */
    mdp_complex &operator()(mdp_site x, mdp_uint n, mdp_uint i, mdp_uint j)
    {
      return address(x, n * m_imax2)[i * m_columns + j];
    }

    /** @brief returns the (i,j) const component of the n-th matrix stored at site x
     */
    const mdp_complex &operator()(mdp_site x, mdp_uint n, mdp_uint i, mdp_uint j) const
    {
      return address(x, n * m_imax2)[i * m_columns + j];
    }
  };
} // namespace MDP

#endif /* MDP_NMATRIX_FIELD_ */
