/////////////////////////////////////////////////////////////////
/// @file fermiqcd_selector.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains a method to set field action and associated field inverter
///
/// Distributed under GPL2 license
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_SELECTOR_
#define FERMIQCD_SELECTOR_

#include <string>

namespace MDP
{
  void select_action_and_inverter(std::string quark_action, std::string inverter)
  {
    if (quark_action == "clover_fast")
    {
      default_fermi_action = FermiCloverActionFast::mul_Q;
    }
    else if (quark_action == "clover_slow")
    {
      default_fermi_action = FermiCloverActionSlow::mul_Q;
    }
#ifdef SSE2
    else if (quark_action == "clover_sse2")
    {
      default_fermi_action = FermiCloverActionSSE2::mul_Q;
    }
#endif
    else if (quark_action == "staggered_fast")
    {
      default_staggered_action = StaggeredAsqtadActionFast::mul_Q;
    }
    else if (quark_action == "staggered_slow")
    {
      default_staggered_action = StaggeredAsqtadActionSlow::mul_Q;
    }
#ifdef SSE2
    else if (quark_action == "staggered_sse2")
    {
      default_staggered_action = StaggeredAsqtadActionSSE2::mul_Q;
    }
#endif
    else
      mdp.error_message("quark action not supported");

    if (inverter == "minres")
    {
      default_fermi_inverter = MinRes::inverter<fermi_field, gauge_field>;
    }
    else if (inverter == "bicgstab")
    {
      default_fermi_inverter = BiCGStab::inverter<fermi_field, gauge_field>;
    }
    else if (inverter == "minres-vtk")
    {
      default_fermi_inverter = MinResVtk::inverter<fermi_field, gauge_field>;
    }
    else if (inverter == "bicgstab-vtk")
    {
      default_fermi_inverter = BiCGStabVtk::inverter<fermi_field, gauge_field>;
    }
    else if (inverter == "bicguml")
    {
      default_staggered_inverter = StaggeredBiCGUML::inverter;
    }
    else if (inverter == "bicgstab-staggered")
    {
      default_staggered_inverter = BiCGStab::inverter<staggered_field, gauge_field>;
    }
    else if (inverter == "minres-staggered")
    {
      default_staggered_inverter = MinRes::inverter<staggered_field, gauge_field>;
    }
    else
      mdp.error_message("quark inverter not supported");
  }
} // namespace MDP

#endif /* FERMIQCD_SELECTOR_ */
