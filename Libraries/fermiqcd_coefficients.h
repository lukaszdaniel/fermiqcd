/////////////////////////////////////////////////////////////////
/// @file fermiqcd_coefficients.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Container for action parameters
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
///
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////
#ifndef FERMIQCD_COEFFICIENTS_
#define FERMIQCD_COEFFICIENTS_

#include <map>
#include <string>
#include "mdp_global_vars.h"
#include "mdp_communicator.h"

namespace MDP
{
  /// @brief container for action parameters
  ///
  /// All FermiQCD actions are classe and share the same prototype.
  /// Parameters are passed to the action via coefficients objects which
  /// are nothing more than hash tables.
  ///
  /// Example:
  /// @verbatim
  ///    gauge_field U(lattice,nc);
  ///    coefficients gauge;
  ///    gauge["beta"]=6.0;
  ///    WilsonGaugeAction::heatbath(U,gauge);
  /// @endverbatim
  /// Please check the spelling of the variables you store into the
  /// coefficients object (each action has its own coefficients).
  ///
  /// Why? This allows the creating of new actions while reusing inverters
  /// and simplify passing parameters to the action.
  class coefficients : public std::map<std::string, mdp_real>
  {
  public:
    bool has_key(const std::string &s) const
    {
      return (find(s) != end());
    }
  };

  void dagger(coefficients &coeff)
  {
    if (!coeff.has_key("sign"))
      coeff["sign"] = -1;
    else
      coeff["sign"] = -coeff["sign"];
  }

  std::ostream &operator<<(std::ostream &os, const coefficients &coeff)
  {
    begin_function("print_coefficients");

    for (const auto &[coefficient, value] : coeff)
    {
      std::cout << "<coefficient name=\"" << coefficient << "\">" << value
                << "</coefficient>"
                << "\n";
    }
    end_function("print_coefficients");
    return os;
  }
} // namespace MDP

#endif /* FERMIQCD_COEFFICIENTS_ */
