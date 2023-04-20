/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_propagator.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various stuff for staggered fermions
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_STAGGERED_PROPAGATOR_
#define FERMIQCD_STAGGERED_PROPAGATOR_

namespace MDP
{
  /// @brief staggared quark propagator
  ///
  /// On a (2n) dimensional lattice this makes 3*(2^n) sources at the
  /// vertices of the hypercube at the origin of the lattice
  /// and inverts the Staggered/Asqtad action on them.
  ///
  /// Example:
  /// @verbatim
  /// mdp_gauge U(lattice,nc);
  /// staggered_propagator S(lattice,nc);
  /// mdp_site x(lattice);
  /// mdp_site y(lattice);
  /// coefficients coeff;
  /// coeff["mass"]=1.0;
  /// generate(S,U,coeff);
  /// for(int i=0; i<(int) pow(2,lattice.ndim); i++) {
  ///    x=binary2versor(a);
  ///    cout << "source at:" << x << "\nprop:\n";
  ///    forallsites(y) cout << S(x,a) << endl;
  /// }
  /// @endverbatim
  class staggered_propagator : public mdp_complex_field
  {
  private:
    mdp_int m_nc;

  public:
    staggered_propagator() : mdp_complex_field()
    {
    }

    staggered_propagator(mdp_lattice &a, int nc_) : mdp_complex_field(a, a.ndim() * a.ndim() * nc_ * nc_), m_nc(nc_)
    {
    }

    staggered_propagator(const staggered_propagator &S) : mdp_complex_field(S), m_nc(S.m_nc)
    {
    }

    mdp_int nc() const
    {
      return m_nc;
    }

    /** @brief returns the matrix of colour \e a stored at site x
     */
    mdp_matrix operator()(mdp_site x, int a)
    {
      mdp_matrix tmp(address(x, a * m_nc * m_nc), m_nc, m_nc);
      return tmp;
    }

    /** @brief returns the (i,j) component of the matrix of colour \e a stored at site x
     */
    mdp_complex &operator()(mdp_site x, int a, int i, int j)
    {
      return *(address(x) + a * m_nc * m_nc + i * m_nc + j);
    }

    friend void generate(staggered_propagator &S, gauge_field &U,
                         coefficients &coeff,
                         mdp_real absolute_precision = fermi_inversion_precision,
                         mdp_real relative_precision = 0,
                         int max_steps = 2000,
                         void (*smf)(staggered_field &, gauge_field &) = nullptr,
                         int comp = 0)
    {
      staggered_field psi(S.lattice(), S.nc());
      staggered_field chi(S.lattice(), S.nc());
      mdp_site x(S.lattice());
      mdp_int ndim = S.lattice().n_dimensions();
      mdp_int nc = S.nc();

      double time = mpi.time();

      if (ME == 0 && shutup == false)
      {
        printf("BEGIN Generating ordinary propagator\n");
        fflush(stdout);
      }

      for (mdp_int a = 0; a < (1 << ndim); a++)
        for (mdp_int j = 0; j < nc; j++)
        {
          forallsitesandcopies(x)
          {
            for (mdp_int i = 0; i < nc; i++)
              psi(x, i) = 0;
          }

          x = binary2versor(a);
          if (ME == 0 && shutup == false)
          {
            printf("(source at (");

            for (mdp_int mu = 0; mu < ndim; mu++)
              printf("%i ", x(mu));

            printf("), Color: %i\n", j);
            fflush(stdout);
          }
          if (x.is_here())
            psi(x, j) = 1;

          /*
            If a smearing function is passed (smf)
            the source is smeared before the inversion
            the sink must be smeared using smear_sink.
          */
          if (smf != nullptr)
            (*smf)(psi, U);

          mul_invQ(chi, psi, U, coeff, absolute_precision, relative_precision, max_steps);

          forallsites(x)
          {
            for (mdp_int i = 0; i < nc; i++)
              S(x, a, i, j) = chi(x, i);
          }
        }

      if (ME == 0 && shutup == false)
      {
        printf("END Generating ordinary propagator. Time: %f (sec)\n",
               mpi.time() - time);
        fflush(stdout);
      }
    }
  };
} // namespace MDP

#endif /* FERMIQCD_STAGGERED_PROPAGATOR_ */
