/////////////////////////////////////////////////////////////////
/// @file fermiqcd_bicgstab_inverter_vtk.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains the stabilized biconjugate inverter
/// From hep-lat/9404013 (by A. Frommer et al.)
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_BICGSTAB_INVERTER_VTK_
#define FERMIQCD_BICGSTAB_INVERTER_VTK_

#include "mdp_field.h"
#include "fermiqcd_default_parameters.h"
#include "fermiqcd_minres_inverter.h"
#include "fermiqcd_minres_inverter_vtk.h"

namespace MDP
{
  /// @brief the stabilized biconjugate inverter
  ///
  /// It inverts mul_Q(psi_out,psi_in,U,coeff) iteratively
  /// @param psi_out the output field passed by reference
  /// @param psi_in the input field passed by reference
  /// @param U the gauge field to be passed to mul_Q
  /// @param coeff the gauge parameters to be passed to mul_Q
  /// @param absolute_precision the target absolute precision
  /// @param relative_precision the target relative precision
  /// @param max_steps the maximum number of steps
  ///
  /// Example:
  /// @verbatim
  /// gauge_field U(lattice,nc);
  /// fermi_field psi(lattice,nc);
  /// fermi_field chi(lattice,nc);
  /// coefficinets coeff;
  /// coeff["kappa"]=1.12;
  /// U.load("myfield");
  /// psi.load("myfield_psi");
  /// default_fermi_inverter=BiCGStab::inverter<fermi_field,gauge_field>;
  /// default_fermi_action=FermiCloverActionSlow::mul_Q;
  /// mul_invQ(chi,psi,U,coeff);
  /// chi.save("myfield_chi");
  /// @endverbatim
  /// Note that mul_invQ(chi,psi,U,coeff) reads
  /// \f$ \chi=(/\!\!\!D[U]+m)^{-1}\psi \f$
  class BiCGStabVtk
  {
  public:
    template <class fieldT, class fieldG>
    static inversion_stats inverter(fieldT &psi_out,
                                    fieldT &psi_in,
                                    fieldG &U,
                                    coefficients &coeff,
                                    mdp_real absolute_precision = mdp_precision,
                                    mdp_real relative_precision = 0,
                                    int max_steps = 2000)
    {
      mdp.begin_function("BiConugateGradientStabilizedInverter");
      const std::string filename_prefix = inversion_vtk_prefix;
      const int tc = 0;
      std::string filename1, filename2;
      int step = 0;
      fieldT p(psi_in);
      fieldT q(psi_in);
      fieldT r(psi_in);
      fieldT s(psi_in);
      fieldT t(psi_in);

      mdp_real_scalar_field sv(psi_in.lattice());
      mdp_site x(psi_in.lattice());
      double residue, rresidue = -1, old_rresidue;
      mdp_complex alpha, beta, rho, rho_old, omega;
      double time = mdp.time();
      inversion_stats stats;

      mdp << "\tstep\tresidue\t\ttime (sec)\n";

      // Initial conditions
      if (BiCGStabRestart == false)
        psi_out = psi_in; // first guess and a stupid guess!

      psi_out.update();
      mul_Q(r, psi_out, U, coeff);
      r *= -1;
      r += psi_in;
      q = r;

      p = 0.0;
      s = 0.0;

      rho_old = alpha = omega = 1;

      mdp << "\t<target>\n"
          << "\t\t<max_steps>" << max_steps << "</max_steps>\n"
          << "\t\t<absolute_precision>" << absolute_precision << "</absolute_precision>\n"
          << "\t\t<relative_precision>" << relative_precision << "</relative_precision>\n"
          << "\t</target>\n";

      do
      {

        rho = q * r;
        beta = (rho / rho_old) * (alpha / omega);
        rho_old = rho;
        p *= beta;
        p += r;
        mdp_add_scaled_field(p, -beta * omega, s);
        p.update();
        mul_Q(s, p, U, coeff);
        alpha = rho / (q * s);
        mdp_add_scaled_field(r, -alpha, s);
        r.update();
        mul_Q(t, r, U, coeff);
        omega = t * r;
        omega /= norm_square(t);
        mdp_add_scaled_field(psi_out, omega, r);
        mdp_add_scaled_field(psi_out, alpha, p);

        // computation of residue
        residue = norm_square(r);
        residue = std::sqrt(residue / r.global_size());

        // computation of rresidue
        old_rresidue = rresidue;
        rresidue = relative_residue(r, psi_out);

        // make VTK files
        for (mdp_int a = 0; a < 4; a++)
        {
          forallsites(x)
          {
            sv(x) = 0.0;
            for (mdp_int k = 0; k < psi_in.nc(); k++)
              sv(x) += std::sqrt(real(psi_out(x, a, k) * conj(psi_out(x, a, k))));
          }
          filename1 = filename_prefix + ".field" + std::to_string(a) + "." + std::to_string(step) + ".vtk";
          sv.save_vtk(filename1, tc);
        }
        forallsites(x)
        {
          for (mdp_int a = 0; a < 4; a++)
          {
            sv(x) = 0.0;
            for (mdp_int k = 0; k < psi_in.nc(); k++)
              sv(x) += std::log(real(r(x, a, k) * conj(r(x, a, k))) + PRECISION);
          }
        }
        filename2 = filename_prefix + ".residue." + std::to_string(step) + ".vtk";
        sv.save_vtk(filename2, tc);

        mdp_add_scaled_field(r, -omega, t);

        mdp << "\t\t<step>" << step << "</step>\n"
            << "\t\t<residue>" << residue << "</residue>\n"
            << "\t\t<relative_residue>" << rresidue << "</relative_residue>\n"
            << "\t\t<time>" << mdp.time() - time << "</time>\n\n";

        if ((step > 10) && (rresidue == old_rresidue))
          error("not converging");
        step++;

      } while (residue > absolute_precision &&
               rresidue > relative_precision &&
               step < max_steps);

      psi_out.update();

      stats.target_absolute_precision = absolute_precision;
      stats.target_relative_precision = relative_precision;
      stats.max_steps = max_steps;
      stats.absolute_precision = residue;
      stats.relative_precision = rresidue;
      stats.residue = residue;
      stats.steps = step;
      stats.mul_Q_steps = 2 * step + 1;
      stats.time = mdp.time() - time;

      mdp << "\t<stats>\n"
          << "\t\t<max_steps>" << step << "</max_steps>\n"
          << "\t\t<absolute_precision>" << residue << "</absolute_precision>\n"
          << "\t\t<relative_precision>" << rresidue << "</relative_precision>\n"
          << "\t\t<time>" << stats.time << "</time>\n"
          << "\t</stats>\n";

      mdp.end_function("BiConugateGradientStabilizedInverter");
      return stats;
    }
  };

  template <class fieldT, class fieldG>
  inversion_stats BiConjugateGradientStabilizedInverterVtk(fieldT &psi_out,
                                                           fieldT &psi_in,
                                                           fieldG &U,
                                                           coefficients &coeff,
                                                           mdp_real absolute_precision = mdp_precision,
                                                           mdp_real relative_precision = 0,
                                                           int max_steps = 2000)
  {
    return BiCGStabVtk::inverter(psi_out, psi_in, U, coeff,
                                 absolute_precision,
                                 relative_precision,
                                 max_steps);
  }
} // namespace MDP

#endif /* FERMIQCD_BICGSTAB_INVERTER_VTK_ */
