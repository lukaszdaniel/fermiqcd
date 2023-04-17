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

namespace MDP
{
  /** @brief field of vectors of vectors (DEPRECATED)
   */
  class mdp_nvector_field : public mdp_field<mdp_complex>
  {
  private:
    mdp_uint m_rows;
    mdp_uint m_imax;
    mdp_uint m_imax2;

  public:
    mdp_nvector_field() : m_rows(0), m_imax(0), m_imax2(0)
    {
      reset_field();
    }

    mdp_nvector_field(mdp_nvector_field &field) : m_rows(field.m_rows), m_imax(field.m_imax), m_imax2(field.m_imax2)
    {
      allocate_field(field.lattice(), field.m_imax);
    }

    /** @brief declares a n-component vector field of i-component vectors at each site
     */
    mdp_nvector_field(mdp_lattice &a, mdp_uint n, mdp_uint i)
    {
      m_rows = i;
      m_imax = i * n;
      m_imax2 = i;
      allocate_field(a, m_imax);
    }

    /** @brief dynamically allocates a n-component vector field of i-component vectors at each site
     */
    void allocate_mdp_nvector_field(mdp_lattice &a, mdp_uint n, mdp_uint i)
    {
      deallocate_field();
      m_rows = i;
      m_imax = i * n;
      m_imax2 = i;
      allocate_field(a, m_imax);
    }

    /** @brief returns the n-th vector stored at site x
     */
    mdp_matrix operator()(mdp_site x, mdp_uint n)
    {
      return mdp_matrix(address(x, n * m_imax2), m_rows, 1);
    }

    /** @brief returns the i-th component of the n-th vector stored at site x
     */
    mdp_complex &operator()(mdp_site x, mdp_uint n, mdp_uint i)
    {
      return address(x, n * m_imax2)[i];
    }

    /** @brief returns the i-th const component of the n-th vector stored at site x
     */
    const mdp_complex &operator()(mdp_site x, mdp_uint n, mdp_uint i) const
    {
      return address(x, n * m_imax2)[i];
    }
  };
} // namespace MDP

#endif /* MDP_NVECTOR_FIELD_ */
