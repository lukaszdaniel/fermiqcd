/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Stuff for SSE/SSE2 compile with -DSSE2
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_STAGGERED_FIELD_
#define FERMIQCD_STAGGERED_FIELD_

#include "mdp_mod2sign.h"
#include "mdp_complex_field.h"

namespace MDP
{
  /// @brief staggered fermionic field
  ///
  /// Example:
  /// @verbatim
  /// staggered_field psi(lattice,nc);
  /// mdp_site x(lattice);
  /// forallsites(x)
  ///   for(int i=0; i<nc; i++)
  ///     psi(x,i)=0.0+0.0*I;
  /// @endverbatim
  class staggered_field : public mdp_complex_field
  {
  private:
    mdp_int m_nspin;
    mdp_int m_nc;

  public:
    staggered_field() : mdp_complex_field(), m_nspin(0), m_nc(0)
    {
    }

    staggered_field(mdp_lattice &a, int nc_, int nspin_ = 4) : mdp_complex_field(a, (nc_)), m_nspin(nspin_), m_nc(nc_)
    {
      // attention here that nspin_ is ignored in field allocation!
    }

    staggered_field(const staggered_field &chi) : mdp_complex_field(chi), m_nspin(chi.m_nspin), m_nc(chi.m_nc)
    {
    }

    mdp_int nspin() const
    {
      return m_nspin;
    }

    mdp_int nc() const
    {
      return m_nc;
    }

    staggered_field &operator=(const staggered_field &chi)
    {
      if (this == &chi)
        return *this;

      // base assignment
      mdp_complex_field::operator=(chi);

      // derived fields assignment
      m_nspin = chi.m_nspin;
      m_nc = chi.m_nc;

      return *this;
    }

    /** @brief returns the vector stored at site x
     */
    mdp_matrix operator()(mdp_site x)
    {
      return mdp_matrix(address(x), m_nc, 1);
    }

    /** @brief returns the i-th component of the vector stored at site x
     */
    mdp_complex &operator()(mdp_site x, int i)
    {
      return *(address(x, i));
    }

    /** @brief returns the i-th const component of the vector stored at site x
     */
    const mdp_complex &operator()(mdp_site x, int i) const
    {
      return *(address(x, i));
    }

    void operator=(mdp_complex a)
    {
      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] = a;
    }

    mdp_real component(mdp_site x, int mu)
    {
      return x(mu) % 2;
    }

    mdp_real eta(mdp_site x, int mu)
    {
#ifdef USE_GOLTERMAN
      mdp_int i_max = (mu + ndim() - 1) % ndim();
      int tmp = 0;
      for (mdp_int i = 1; i <= i_max; i++)
        tmp += x(i);
#else
      int tmp = 0;
      for (mdp_int i = 0; i < mu; i++)
        tmp += x(i);
#endif
      return mdp_mod2sign(tmp);
    }

    mdp_real eps(mdp_site x)
    {
      int tmp = x(0);
      for (mdp_int i = 1; i < ndim(); i++)
        tmp += x(i);
      return mdp_mod2sign(tmp);
    }

    mdp_real type(mdp_site x)
    {
      mdp_real tmp = x(0) % 2;
      for (mdp_int i = 1; i < ndim(); i++)
        tmp += (x(i) % 2) * std::pow(2.0, i);
      return tmp;
    }
  };
} // namespace MDP

#endif /* FERMIQCD_STAGGERED_FIELD_ */
