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

namespace MDP
{
  void rotate_field(fermi_field &psi, gauge_field &U, coefficients &coeff)
  {
    fermi_field psi_rot(psi.lattice(), psi.nspin, psi.nc);
    mdp_real d1 = coeff["rotate_field:d1"];
    psi_rot = psi;
    mdp_site x(psi.lattice());
    forallsites(x)
    {
      for (int mu = 1; mu < U.ndim(); mu++)
        for (int a = 0; a < psi.nspin; a++)
          for (int b = 0; b < psi.nspin; b++)
          {
            psi(x, a) += d1 * 0.5 * Gamma[mu](a, b) * (U(x, mu) * psi_rot(x + mu, b) - (U(x, -1, mu)) * psi_rot(x - mu, b));
          }
    }
    psi.update();
  }
} // namespace MDP

#endif /* FERMIQCD_FERMI_ROTATION_ */
