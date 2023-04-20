/////////////////////////////////////////////////////////////////
/// @file fermiqcd_sdwf_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// WORK IN PROGRESS
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_SDWF_FIELD_
#define FERMIQCD_SDWF_FIELD_

namespace MDP
{
  /** @brief field for domain wall staggered fermions
   */
  class sdwf_field : public mdp_complex_field
  {
  private:
    mdp_int m_nspin;
    mdp_int m_nc;
    mdp_int m_L5;

  public:
    sdwf_field() : mdp_complex_field(), m_nspin(0), m_nc(0), m_L5(0)
    {
    }

    sdwf_field(mdp_lattice &a, int L5_, int nc_, int nspin_ = 4) : mdp_complex_field(a, (L5_ * nc_)), m_nspin(nspin_), m_nc(nc_), m_L5(L5_)
    {
      // attention here that nspin_ is ignored in field allocation!
    }

    sdwf_field(const sdwf_field &chi) : mdp_complex_field(chi), m_nspin(chi.m_nspin), m_nc(chi.m_nc), m_L5(chi.m_L5)
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

    mdp_int L5() const
    {
      return m_L5;
    }

    /** @brief returns the \e x5 component of the vector stored at site x
     */
    mdp_matrix operator()(mdp_site x, int x5)
    {
      return mdp_matrix(address(x, x5 * m_nc), m_nc, 1);
    }

    /** @brief returns the \e x5 component of the vector of colour \e i stored at site x
     */
    mdp_complex &operator()(mdp_site x, int x5, int i)
    {
      return *(address(x, x5 * m_nc + i));
    }

    /** @brief returns the \e x5 const component of the vector of colour \e i stored at site x
     */
    const mdp_complex &operator()(mdp_site x, int x5, int i) const
    {
      return *(address(x, x5 * m_nc + i));
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
      int tmp = 0;
      int i_max = (mu + ndim() - 1) % ndim();

      for (int i = 1; i <= i_max; i++)
        tmp += x(i);
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
      mdp_real tmp;
      tmp = x(0) % 2;
      for (mdp_int i = 1; i < ndim(); i++)
        tmp += (x(i) % 2) * std::pow(2.0, i);
      return tmp;
    }

    mdp_site chiral_shift(mdp_site x)
    {
      for (mdp_int i = 0; i < ndim(); i++)
        if (x(i) % 2 == 1)
          x = x - i;
        else
          x = x + i;
      return x;
    }

    mdp_real chiral_phase(mdp_site x)
    { // (Gamma5 (x) 1)
      int tmp = ndim() / 2;
      for (int i = 1; i < ndim(); i += 2)
        tmp += x(i);
      return mdp_mod2sign(tmp);
    }

    mdp_real chiral_phase2(mdp_site x)
    { // (Gamma5 (x) Gamma5)
      int tmp = 0;
      for (int i = 0; i < ndim(); i++)
        tmp += x(i);
      return mdp_mod2sign(tmp);
    }
  };
} // namespace MDP

#endif /* FERMIQCD_SDWF_FIELD_ */
