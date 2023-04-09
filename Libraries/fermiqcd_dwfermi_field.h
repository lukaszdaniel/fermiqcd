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
  public:
    int nspin;
    int nc;
    int L5;

    dwfermi_field()
    {
      reset_field();
      L5 = 0;
    }

    dwfermi_field(mdp_lattice &a, int L5_, int nc_, int nspin_ = 4)
    {
      deallocate_field();
      L5 = L5_;
      nc = nc_;
      nspin = nspin_;
      allocate_field(a, L5 * nspin * nc);
    }

    dwfermi_field(const dwfermi_field &psi) : mdp_complex_field(psi)
    {
      L5 = psi.L5;
      nc = psi.nc;
      nspin = psi.nspin;
    }

    void allocate_dwfermi_field(mdp_lattice &a,
                                int L5_, int nc_, int nspin_ = 4)
    {
      deallocate_field();
      L5 = L5_;
      nc = nc_;
      nspin = nspin_;
      allocate_field(a, L5 * nspin * nc);
    }

    /** @brief returns the matrix stored at site x and site \e L5_ in 5-th dimension
     */
    mdp_matrix operator()(mdp_site x, int L5_)
    {
      return mdp_matrix(address(x, L5_ * nc * nspin), nspin, nc);
    }

    /** @brief returns the vector of colour \e a stored at site x and site \e L5_ in 5-th dimension
     */
    mdp_matrix operator()(mdp_site x, int L5_, int a)
    {
      return mdp_matrix(address(x, (L5_ * nspin + a) * nc), nc, 1);
    }

    /** @brief returns the (a,i) component of the matrix stored at site x and site \e L5_ in 5-th dimension
     */
    mdp_complex &operator()(mdp_site x, int L5_, int a, int i)
    {
      return *(address(x, (L5_ * nspin + a) * nc + i));
    }

    /** @brief returns the (a,i) const component of the matrix stored at site x and site \e L5_ in 5-th dimension
     */
    const mdp_complex &operator()(mdp_site x, int L5_, int a, int i) const
    {
      return *(address(x, (L5_ * nspin + a) * nc + i));
    }

    void operator=(mdp_complex a)
    {
      for (mdp_int i = 0; i < size; i++)
        m[i] = a;
    }
  };
} // namespace MDP

#endif /* FERMIQCD_DWFERMI_FIELD_ */
