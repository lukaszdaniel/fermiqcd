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
  /// @brief field of vectors of vectors (DEPRECATED)
  class mdp_nvector_field : public mdp_field<mdp_complex>
  {
  public:
    mdp_uint rows;
    mdp_uint columns;
    mdp_uint imax;
    mdp_uint imax2;

    mdp_nvector_field()
    {
      rows = columns = imax = imax2 = 0;
      reset_field();
    }

    mdp_nvector_field(mdp_nvector_field &field)
    {
      rows = field.rows;
      columns = field.columns;
      imax = field.imax;
      imax2 = field.imax2;
      allocate_field(field.lattice(), field.imax);
    }

    /** @brief declares a n-component vector field of i-component vectors at each site
     */
    mdp_nvector_field(mdp_lattice &a, int n, int i)
    {
      rows = i;
      columns = 1;
      imax = i * n;
      imax2 = i;
      allocate_field(a, imax);
    }

    /** @brief dynamically allocates a n-component vector field of i-component vectors at each site
     */
    void allocate_mdp_nvector_field(mdp_lattice &a, int n, int i)
    {
      deallocate_field();
      rows = i;
      columns = 1;
      imax = i * n;
      imax2 = i;
      allocate_field(a, imax);
    }

    /** @brief returns the n-th vector stored at site x
     */
    mdp_matrix operator()(mdp_site x, int n)
    {
      return mdp_matrix(address(x, n * imax2), rows, columns);
    }

    /** @brief returns the i-th component of the n-th vector stored at site x
     */
    mdp_complex &operator()(mdp_site x, int n, int i)
    {
      return address(x, n * imax2)[i];
    }

    /** @brief returns the i-th const component of the n-th vector stored at site x
     */
    const mdp_complex &operator()(mdp_site x, int n, int i) const
    {
      return address(x, n * imax2)[i];
    }
  };
} // namespace MDP

#endif /* MDP_NVECTOR_FIELD_ */
