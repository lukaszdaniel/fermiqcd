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

namespace MDP
{
  /// @brief wilson fermionic field
  ///
  /// Example:
  /// @verbatim
  /// fermi_field psi(lattice,nc);
  /// mdp_site x(lattice);
  /// forallsites(x)
  ///    for(int spin=0; spin<4; spin++)
  ///      for(int i=0; i<nc; i++)
  ///        psi(x,spin,i)=0.0+0.0*I;
  /// @endverbatim
  class fermi_field : public mdp_complex_field
  {
  public:
    int nspin;
    int nc;

    fermi_field()
    {
      reset_field();
    }

    fermi_field(mdp_lattice &a, int nc_, int nspin_ = 4)
    {
      deallocate_field();
      nc = nc_;
      nspin = nspin_;
      allocate_field(a, nspin * nc);
    }

    void allocate_fermi_field(mdp_lattice &a, int nc_, int nspin_ = 4)
    {
      deallocate_field();
      nc = nc_;
      nspin = nspin_;
      allocate_field(a, nspin * nc);
    }

    fermi_field(const fermi_field &chi) : mdp_complex_field(chi)
    {
      nc = chi.nc;
      nspin = chi.nspin;
    }

    void operator=(const fermi_field &chi)
    {
      nc = chi.nc;
      nspin = chi.nspin;
      mdp_complex_field::operator=(chi);
    }

    /** @brief returns the matrix stored at site x
     */
    inline mdp_matrix operator()(mdp_site x)
    {
      return mdp_matrix(address(x), nspin, nc);
    }

    /** @brief returns the vector of spin \e a stored at site x
     */
    inline mdp_matrix operator()(mdp_site x, int a)
    {
      return mdp_matrix(address(x, a * nc), nc, 1);
    }

    /** @brief returns the (a,i) component of the matrix stored at site x
     */
    inline mdp_complex &operator()(mdp_site x, int a, int i)
    {
      return *(address(x, a * nc + i));
    }

    /** @brief returns the (a,i) const component of the matrix stored at site x
     */
    inline const mdp_complex &operator()(mdp_site x, int a, int i) const
    {
      return *(address(x, a * nc + i));
    }

    void operator=(mdp_complex a)
    {
      for (mdp_int i = 0; i < size; i++)
        m[i] = a;
    }
  };

  // //////////////////////////////////////////////////////////////////////
  // function to print a quark field by components. It prompts for a site //
  // similar to the CANOPY one. I used this to test the converter        //
  /////////////////////////////////////////////////////////////////////////

  void print_fermi_field(fermi_field &psi)
  {
    begin_function("print_fermi_field");
    int x0, x1, x2, x3;
    mdp_site x(psi.lattice());
    int do_exit = false;
    do
    {
      mdp << "\nCheck point!\n";
      mdp << "Here you called the function to print the propagator\n";
      mdp << "Enter the coordinates (x0,x1,x2,x3 or 'quit' to end): ";
      if (ME == 0)
      {
        std::string stringa;
        std::cin >> stringa;
        if (stringa == "quit")
          do_exit = true;
        else
          sscanf(stringa.c_str(), "%i,%i,%i,%i", &x0, &x1, &x2, &x3);
      }
      mpi.broadcast(do_exit, 0);
      if (do_exit == true)
      {
        mdp << "\n";
        break;
      }
      mpi.broadcast(x0, 0);
      mpi.broadcast(x1, 0);
      mpi.broadcast(x2, 0);
      mpi.broadcast(x3, 0);
      if (on_which_process(psi.lattice(), x0, x1, x2, x3) == ME)
      {
        x.set(x0, x1, x2, x3);
        mdp << psi(x);
      }
    } while (1);
    end_function("print_fermi_field");
  }
} // namespace MDP

#endif /* FERMIQCD_FERMI_FIELD_ */
