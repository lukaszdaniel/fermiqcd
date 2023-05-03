/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_propagator.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class fermi_propagator
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_FERMI_PROPAGATOR_
#define FERMIQCD_FERMI_PROPAGATOR_

namespace MDP
{
  /// @brief a Wilson/Clover quark propagator (all 12 components)
  ///
  /// Example of how to make a pion:
  /// @verbatim
  /// gauge_field U(lattice,nc);
  /// U.load("myfield");
  /// fermi_propagator S(lattice,nc);
  /// coefficients quark;
  /// quark["kappa"]=1.12;
  /// fermi_propagator::generate(S,U,quark);
  /// vector<float> sum(U.lattice.size(TIME));
  /// forallsites(x)
  ///   for(int alpha=0; alpha<4; alpha++)
  ///     for(int beta=0; beta<4; beta++)
  ///        sum(x(0))+=real(trace(S(x,alpha,beta)*
  ///                   hermitian(S(x,beta,alpha))));
  /// @endverbatim
  /// Note that S(x,alpha,beta,i,j) is
  /// \f$ \left<0|\bar q^i_\alpha(x), q^j_\beta(0)|\right> \f$
  class fermi_propagator : public mdp_complex_field
  {
  private:
    mdp_int m_nspin;
    mdp_int m_nc;

  public:
    fermi_propagator() : mdp_complex_field(), m_nspin(0), m_nc(0)
    {
    }

    fermi_propagator(mdp_lattice &a, int nc_, int nspin_ = 4) : mdp_complex_field(a, (nc_ * nc_ * nspin_ * nspin_)), m_nspin(nspin_), m_nc(nc_)
    {
    }

    fermi_propagator(const fermi_propagator &chi) : mdp_complex_field(chi), m_nspin(chi.m_nspin), m_nc(chi.m_nc)
    {
    }

    void allocate_fermi_propagator(mdp_lattice &a, int nc_, int nspin_ = 4)
    {
      deallocate_field();
      m_nspin = nspin_;
      m_nc = nc_;
      allocate_field(a, m_nspin * m_nspin * m_nc * m_nc);
    }

    mdp_int nspin() const
    {
      return m_nspin;
    }

    mdp_int nc() const
    {
      return m_nc;
    }

    /** @brief returns the matrix \e alpha, \e beta stored at site x
     */
    mdp_matrix operator()(mdp_site x, int a, int b)
    {
      mdp_matrix tmp(address(x, (a * m_nspin + b) * m_nc * m_nc), m_nc, m_nc);
      return tmp;
    }

    /** @brief returns the (i,j) component of the matrix \e alpha, \e beta stored at site x
     */
    mdp_complex &operator()(mdp_site x, int a, int b, int i, int j)
    {
      return *(address(x, ((a * m_nspin + b) * m_nc + i) * m_nc + j));
    }

    /// makes the quark propagator
    ///
    /// @param S the output propagator
    /// @param U the input gauge configuration
    /// @param coeff the parameters to be passed to the action
    /// @param absolute_precision the target absolute precision for inversion
    /// @param relative_precision the target relative precision for invcersion
    /// @param max_steps the max number of steps in inversion
    /// @param smf pointer to smearing function (smear sources)
    /// @param smear_coeff parameters for smearing
    static void generate(fermi_propagator &S, gauge_field &U,
                         coefficients &coeff,
                         mdp_real absolute_precision = fermi_inversion_precision,
                         mdp_real relative_precision = 0,
                         int max_steps = 2000,
                         void (*smf)(fermi_field &,
                                     gauge_field &,
                                     coefficients &) = nullptr,
                         coefficients smear_coeff = coefficients())
    {
      fermi_field psi(S.lattice(), S.nc(), S.nspin());
      fermi_field chi(S.lattice(), S.nc(), S.nspin());
      mdp_site x(S.lattice());

      double time = mdp.time();
      inversion_stats stats;
      begin_function("generate");
      mdp << "BEGIN Generating ordinary propagator\n";

      for (mdp_int b = 0; b < psi.nspin(); b++)
        for (mdp_int j = 0; j < psi.nc(); j++)
        {

          mdp << "Source: spin=" << b << ", color=" << j << "\n";

          forallsitesandcopies(x)
          {
            for (mdp_int a = 0; a < psi.nspin(); a++)
              for (mdp_int i = 0; i < psi.nc(); i++)
              {
                if ((x.is_equal(0)) && (a == b) && (i == j))
                  psi(x, a, i) = 1;
                else
                  psi(x, a, i) = 0;
              }
          }
          /*
            If a smearing function is passed (smf)
            the source is smeared before the inversion
            the sink must be smeared using smear_sink.
          */
          if (smf != nullptr)
            (*smf)(psi, U, smear_coeff);
          stats = mul_invQ(chi, psi, U, coeff,
                           absolute_precision, relative_precision, max_steps);

          forallsites(x)
          {
            for (mdp_int a = 0; a < psi.nspin(); a++)
              for (mdp_int i = 0; i < psi.nc(); i++)
              {
                S(x, a, b, i, j) = chi(x, a, i);
              }
          }
          std::cout << "Statistics: residue=" << stats.residue
                    << ", steps=" << stats.steps
                    << ", time=" << stats.time << "\n";
        }
      mdp << "END Generating ordinary propagator. ";
      mdp << "time=" << mdp.time() - time << "\n";
      end_function("generate");
    }
  };

  // //////////////////////////////////////////////////////////////////////
  // function to print a propagator by components. It prompts for a site //
  // similar to the CANOPY one. I used this to test the converter        //
  /////////////////////////////////////////////////////////////////////////

  void print_propagator(fermi_propagator &S)
  {
    begin_function("print_propagator");
    int x0, x1, x2, x3;
    mdp_site x(S.lattice());
    mdp_complex tmp;
    int do_exit = false;
    int nc = S.nc();
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
          sscanf(stringa.c_str(), "%i,%i,%i,%i", &x0, &x1, &x2, &x3);
      }
      mdp.broadcast(do_exit, 0);
      if (do_exit == true)
      {
        mdp << "\n";
        break;
      };
      mdp.broadcast(x0, 0);
      mdp.broadcast(x1, 0);
      mdp.broadcast(x2, 0);
      mdp.broadcast(x3, 0);
      if (on_which_process(S.lattice(), x0, x1, x2, x3) == ME)
      {
        x.set(x0, x1, x2, x3);
        for (mdp_int color_source = 0; color_source < nc; color_source++)
          for (mdp_int spin_source = 0; spin_source < 4; spin_source++)
          {
            mdp << "Source: spin=" << spin_source
                << ", color=" << color_source << "\n";
            for (mdp_int spin_sink = 0; spin_sink < 4; spin_sink++)
            {
              mdp << "[ ";
              for (mdp_int color_sink = 0; color_sink < nc; color_sink++)
              {
                tmp = S(x, spin_sink, spin_source, color_sink, color_source);
                mdp << tmp;
                if (color_sink < nc - 1)
                  mdp << ",\t";
              };
              mdp << " ]\n";
            };
          };
        fflush(stdout);
      };
    } while (1);
    begin_function("print_propagator");
  }

  void smear_propagator(fermi_propagator &S, gauge_field &U, int smear_steps = 10, float alpha = 1.0)
  {
    mdp_matrix_field V(U.lattice(), U.nc(), U.nc());
    mdp_site x(U.lattice());
    for (mdp_int n = 0; n < smear_steps; n++)
    {
      for (mdp_int a = 0; a < 4; a++)
        for (mdp_int b = 0; b < 4; b++)
        {
          forallsites(x)
          {
            V(x) = alpha * S(x, a, b);
            for (mdp_int mu = 0; mu < 4; mu++)
              V(x) += U(x, mu) * S(x + mu, a, b) + hermitian(U(x - mu, mu)) * S(x - mu, a, b);
          }
          V.update();
          forallsites(x)
          {
            for (mdp_int i = 0; i < U.nc(); i++)
              for (mdp_int j = 0; j < U.nc(); j++)
                S(x, a, b, i, j) = V(x, i, j) / (8.0 + alpha);
          }
        }
    }
  }
} // namespace MDP

#endif /* FERMIQCD_FERMI_PROPAGATOR_ */
