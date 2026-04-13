/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains the class fermi_field
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_FERMI_FIELD_
#define FERMIQCD_FERMI_FIELD_

#include <iostream>
#include "mdp_field.h"

namespace MDP
{
  /// @brief wilson fermionic field
  ///
  /// Example:
  /// @verbatim
  /// fermi_field psi(lattice,nc);
  /// mdp_site x(lattice);
  /// forallsites(x)
  ///    for(mdp_suint spin=0; spin<4; spin++)
  ///      for(mdp_suint i=0; i<nc; i++)
  ///        psi(x,spin,i)=0.0+0.0*I;
  /// @endverbatim
  class fermi_field : public mdp_complex_field
  {
  private:
    mdp_suint m_nspin;
    mdp_suint m_nc;

  public:
    fermi_field() : mdp_complex_field(), m_nspin(0), m_nc(0)
    {
    }

    fermi_field(const mdp_lattice &a, mdp_suint nc_, mdp_suint nspin_ = 4) : mdp_complex_field(a, (nc_ * nspin_)), m_nspin(nspin_), m_nc(nc_)
    {
    }

    fermi_field(const fermi_field &chi) : mdp_complex_field(chi), m_nspin(chi.m_nspin), m_nc(chi.m_nc)
    {
    }

    void allocate_fermi_field(const mdp_lattice &a, mdp_suint nc_, mdp_suint nspin_ = 4)
    {
      deallocate_field();
      m_nspin = nspin_;
      m_nc = nc_;
      allocate_field(a, m_nspin * m_nc);
    }

    mdp_suint nspin() const
    {
      return m_nspin;
    }

    mdp_suint nc() const
    {
      return m_nc;
    }

    fermi_field &operator=(const fermi_field &chi)
    {
      if (this == &chi)
        return *this;

      m_nspin = chi.m_nspin;
      m_nc = chi.m_nc;

      return *this;
    }

    /** @brief returns the matrix stored at site x
     */
    mdp_matrix operator()(mdp_site x) const
    {
      return mdp_matrix(address(x), m_nspin, m_nc);
    }

    /** @brief returns the vector of spin \e a stored at site x
     */
    mdp_matrix operator()(mdp_site x, mdp_suint a) const
    {
      return mdp_matrix(address(x, a * m_nc), m_nc, 1);
    }

    /** @brief returns the (a,i) component of the matrix stored at site x
     */
    mdp_complex &operator()(mdp_site x, mdp_suint a, int i)
    {
      return *(address(x, a * m_nc + i));
    }

    /** @brief returns the (a,i) const component of the matrix stored at site x
     */
    const mdp_complex &operator()(mdp_site x, mdp_suint a, int i) const
    {
      return *(address(x, a * m_nc + i));
    }

    void operator=(mdp_complex a)
    {
      for (mdp_uint i = 0; i < m_size; i++)
        m_data[i] = a;
    }
  };

  // //////////////////////////////////////////////////////////////////////
  // function to print a quark field by components. It prompts for a site //
  // similar to the CANOPY one. I used this to test the converter        //
  /////////////////////////////////////////////////////////////////////////

  void print_fermi_field(fermi_field &psi)
  {
    begin_function("print_fermi_field");
    mdp_uint x0, x1, x2, x3;
    mdp_site x(psi.lattice());
    bool do_exit = false;
    do
    {
      mdp << "\nCheck point!\n";
      mdp << "Here you called the function to print the propagator\n";
      mdp << "Enter the coordinates (x0,x1,x2,x3 or 'quit' to end): ";
      if (isMainProcess())
      {
        std::string stringa;
        std::cin >> stringa;
        if (stringa == "quit")
          do_exit = true;
        else
        {
          std::stringstream ss(stringa);
          char comma;
          ss >> x0 >> comma >> x1 >> comma >> x2 >> comma >> x3;
        }
      }
      mdp.broadcast(do_exit, 0);
      if (do_exit)
      {
        mdp << "\n";
        break;
      }
      mdp.broadcast(x0, 0);
      mdp.broadcast(x1, 0);
      mdp.broadcast(x2, 0);
      mdp.broadcast(x3, 0);
      if (on_which_process(psi.lattice(), x0, x1, x2, x3) == ME)
      {
        x.set(x0, x1, x2, x3);
        mdp << psi(x);
      }
    } while (true);
    end_function("print_fermi_field");
  }
} // namespace MDP

#endif /* FERMIQCD_FERMI_FIELD_ */
