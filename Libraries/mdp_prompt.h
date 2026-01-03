/////////////////////////////////////////////////////////////////
/// @file mdp_prompt.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Functions to parse user input of parameters
/// in a way safe to parallel programs
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_PROMPT_
#define MDP_PROMPT_

#include <string>
#include <string_view>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "mdp_macros.h"
#include "mdp_communicator.h"

namespace MDP
{
  inline constexpr std::string_view STD_INPUT = "";
  inline constexpr std::string_view STD_INPUT_FILE = "<stdin>";

  /// Converts string to float
  double val(const std::string &s)
  {
    return std::stod(s);
  }

  static std::string trim(const std::string &s)
  {
    const auto begin = s.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
      return "";
    const auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(begin, end - begin + 1);
  }

  /// Try prompt("<stdin>","VALUE","4.0")
  /// It will prompt the user for variable VALUE and take 4.0 as default
  std::string prompt(const std::string &filename,
                     const std::string &variable,
                     const std::string &def_val = "0.0",
                     int p = 0)
  {
    std::string response = def_val;

#ifdef PARALLEL
    mdp.barrier();
#endif

    if (isSubProcess(p))
    {
      if (filename == STD_INPUT)
      {
        std::cout << "Input value of " << variable << " (default is '" << def_val << "'): ";
        std::cin >> response;
      }
      else
      {
        std::istream *in = nullptr;
        std::ifstream file;

        if (filename == STD_INPUT_FILE)
        {
          std::cout << "Input variable (default is '" << variable << " " << def_val << "')...\n";
          in = &std::cin;
        }
        else
        {
          file.open(filename);
          if (!file)
            error("Unable to open file");
          in = &file;
        }

        std::string line;
        while (std::getline(*in, line))
        {
          // ignore comments
          line = trim(line);
          if (line.empty() || line[0] == '#')
            continue;

          if (auto pos = line.find('#'); pos != std::string::npos)
            line = trim(line.substr(0, pos));

          std::istringstream iss(line);
          std::string name, value, eq;

          if (!(iss >> name))
            continue;

          if (name != variable)
            continue;

          // handle: name value  OR  name = value
          if (iss >> eq && eq == "=")
            iss >> value;
          else
            value = eq;

          if (!value.empty())
          {
            response = value;
            break;
          }
        }
      }

      std::cout << "... Adopting " << variable << " equal to \"" << response << "\"\n\n";
    }

#ifdef PARALLEL
    mdp.broadcast(response.data(), response.size() + 1, p);
#endif

    return response;
  }
} // namespace MDP

#endif /* MDP_PROMPT_ */
