/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_rotation.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains methods for field rotation
///
/// Distributed under GPL2 license
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_FERMI_ROTATION_
#define FERMIQCD_FERMI_ROTATION_

#include "fermiqcd_fermi_field.h"
#include "fermiqcd_gauge_field.h"
#include "fermiqcd_gamma_matrices.h"
#include "fermiqcd_coefficients.h"

namespace MDP
{
  /**
   * @brief Applies a rotation (finite-difference derivative-like transformation)
   *        to a fermionic field using the gauge field and given coefficients.
   *
   * This function modifies the input fermion field `psi` in place by adding
   * a contribution proportional to a covariant nearest-neighbor difference.
   * The transformation involves gamma matrices and gauge links, resembling
   * a discretized derivative operator acting on the field.
   *
   * @param psi Fermionic field to be rotated. The result is accumulated directly
   *   into this field.
   *
   * @param U Gauge field providing link variables used to transport fermion
   *   fields between neighboring lattice sites.
   *
   * @param coeff Container of coefficients. The parameter `"rotate_field:d1"`
   *   controls the strength of the rotation term.
   */
  void rotate_field(fermi_field &psi, const gauge_field &U, coefficients &coeff)
  {
    fermi_field psi_rot(psi);
    mdp_real d1 = coeff["rotate_field:d1"];
    mdp_site x(psi.lattice());

    forallsites(x)
    {
      for (mdp_suint mu = 1; mu < U.ndim(); mu++)
      {
        for (mdp_suint a = 0; a < psi.nspin(); a++)
        {
          for (mdp_suint b = 0; b < psi.nspin(); b++)
          {
            psi(x, a) += d1 * 0.5 * Gamma[mu](a, b) * (U(x, mu) * psi_rot(x + mu, b) - (U(x, -1, mu)) * psi_rot(x - mu, b));
          }
        }
      }
    }
    psi.update();
  }
} // namespace MDP

#endif /* FERMIQCD_FERMI_ROTATION_ */
