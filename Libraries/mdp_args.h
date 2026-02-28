/////////////////////////////////////////////////////////////////
/// @file mdp_args.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Other junk that did not fit anywhere else
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////
#ifndef MDP_ARGS_
#define MDP_ARGS_

#include <string>
#include <vector>
#include "mdp_global_vars.h"
#include "mdp_utils.h"

namespace MDP
{
  class mdp_args
  {
  private:
    std::vector<std::string> m_args;

  public:
    mdp_args(int argc, char **argv)
    {
      for (int i = 1; i < argc; i++)
        m_args.push_back(argv[i]);
    }

    size_t length() const
    {
      return m_args.size();
    }

    bool have(const std::string &name) const
    {
      for (size_t i = 0; i < m_args.size(); i++)
      {
        if (m_args[i] == name || startswith(m_args[i], name + ":"))
          return true;
      }
      return false;
    }

    float get(const std::string &name, const std::string &key, float value = 0.0) const
    {
      for (size_t i = 0; i < m_args.size(); i++)
      {
        if (startswith(m_args[i], name + ":"))
          return parse_float(m_args[i], key, value);
      }
      return value;
    }

    float get(const std::string &name, const std::string &key, double value = 0.0) const
    {
      for (size_t i = 0; i < m_args.size(); i++)
      {
        if (startswith(m_args[i], name + ":"))
        {
          value = parse_float(m_args[i], key, value);
          break;
        }
      }
      mdp << "INPUT " << name << ":" << key << "=" << value << "\n";
      return value;
    }

    int get(const std::string &name, const std::string &key, int value = 0) const
    {
      for (size_t i = 0; i < m_args.size(); i++)
      {
        if (startswith(m_args[i], name + ":"))
        {
          value = parse_int(m_args[i], key, value);
          break;
        }
      }
      mdp << "INPUT " << name << ":" << key << "=" << value << "\n";
      return value;
    }

    mdp_uint get(const std::string &name, const std::string &key, mdp_uint value = 0) const
    {
      for (size_t i = 0; i < m_args.size(); i++)
      {
        if (startswith(m_args[i], name + ":"))
        {
          value = parse_uint(m_args[i], key, value);
          break;
        }
      }
      mdp << "INPUT " << name << ":" << key << "=" << value << "\n";
      return value;
    }

    std::string get(const std::string &name, const std::string &key, std::string value = "") const
    {
      int i = value.find('|');
      if (i >= 0)
        value = value.substr(0, i);
      for (size_t i = 0; i < m_args.size(); i++)
      {
        if (startswith(m_args[i], name + ":"))
        {
          value = parse_string(m_args[i], key, value);
          break;
        }
      }
      mdp << "INPUT " << name << ":" << key << "=" << value << "\n";
      return value;
    }
  };
} // namespace MDP

#endif /* MDP_ARGS_ */
