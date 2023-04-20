/////////////////////////////////////////////////////////////////
/// @file fermiqcd_dwfermi_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains the class dwfermi_field fro domain wall fermions
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_DWFERMI_FIELD_
#define FERMIQCD_DWFERMI_FIELD_

namespace MDP
{
  /// @brief domain wall fermionic field
  ///
  /// Example:
  /// @verbatim
  /// int L5=10; // size in 5th dimension
  /// fermi_field psi(lattice,L5,nc);
  /// mdp_site x(lattice);
  /// forallsites(x)
  ///    for(int k=0; k<L5; k++)
  ///      for(int spin=0; spin<4; spin++)
  ///        for(int i=0; i<nc; i++)
  ///          psi(x,k,spin,i)=0.0+0.0*I;
  /// @endverbatim
  class dwfermi_field : public mdp_complex_field
  {
  private:
    mdp_int m_nspin;
    mdp_int m_nc;
    mdp_int m_L5;

  public:
    dwfermi_field() : mdp_complex_field(), m_nspin(0), m_nc(0), m_L5(0)
    {
    }

    dwfermi_field(mdp_lattice &a, int L5_, int nc_, int nspin_ = 4) : mdp_complex_field(a, (L5_ * nc_ * nspin_)), m_nspin(nspin_), m_nc(nc_), m_L5(L5_)
    {
    }

    dwfermi_field(const dwfermi_field &psi) : mdp_complex_field(psi), m_nspin(psi.m_nspin), m_nc(psi.m_nc), m_L5(psi.m_L5)
    {
    }

    void allocate_dwfermi_field(mdp_lattice &a, int L5_, int nc_, int nspin_ = 4)
    {
      deallocate_field();
      m_L5 = L5_;
      m_nc = nc_;
      m_nspin = nspin_;
      allocate_field(a, m_L5 * m_nspin * m_nc);
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

    /** @brief returns the matrix stored at site x and site \e L5_ in 5-th dimension
     */
    mdp_matrix operator()(mdp_site x, int L5_)
    {
      return mdp_matrix(address(x, L5_ * m_nc * m_nspin), m_nspin, m_nc);
    }

    /** @brief returns the vector of colour \e a stored at site x and site \e L5_ in 5-th dimension
     */
    mdp_matrix operator()(mdp_site x, int L5_, int a)
    {
      return mdp_matrix(address(x, (L5_ * m_nspin + a) * m_nc), m_nc, 1);
    }

    /** @brief returns the (a,i) component of the matrix stored at site x and site \e L5_ in 5-th dimension
     */
    mdp_complex &operator()(mdp_site x, int L5_, int a, int i)
    {
      return *(address(x, (L5_ * m_nspin + a) * m_nc + i));
    }

    /** @brief returns the (a,i) const component of the matrix stored at site x and site \e L5_ in 5-th dimension
     */
    const mdp_complex &operator()(mdp_site x, int L5_, int a, int i) const
    {
      return *(address(x, (L5_ * m_nspin + a) * m_nc + i));
    }

    void operator=(mdp_complex a)
    {
      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] = a;
    }
  };
} // namespace MDP

#endif /* FERMIQCD_DWFERMI_FIELD_ */
