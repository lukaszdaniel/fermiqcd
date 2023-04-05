/////////////////////////////////////////////////////////////////
/// @file fermiqcd_default_parameters.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Constants and parameters used by FermiQCD
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_DEFAULT_PARAMETERS_
#define FERMIQCD_DEFAULT_PARAMETERS_

namespace MDP
{
#define TIME 0
#define SPACE_X 1
#define SPACE_Y 2
#define SPACE_Z 3

#ifndef DAGGER
#define DAGGER -1
#endif

    mdp_real fermi_inversion_precision = 1e-6;
    mdp_real staggered_inversion_precision = 1e-6;
    mdp_real dwfermi_inversion_precision = 1e-6;
    mdp_real sdwf_inversion_precision = 1e-6;

    /// Set this to true to run BuCGStab with restart
    bool BiCGStabRestart = false;

#ifdef BLOCKSITE
    mdp_matrix OmegaTwist[BLOCKSITE];
#endif
} // namespace MDP

#endif /* FERMIQCD_DEFAULT_PARAMETERS_ */