/////////////////////////////////////////////////////////////////
/// @file mdp_log.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_log
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_LOG_
#define MDP_LOG_

#include <vector>
#include <fstream>
#include <ostream>
#include <iostream>
#include <string>

namespace MDP
{
  /// @brief base class of class mdp_communicator (DO NOT INSTANTIATE)
  /// @see class mdp_communicator
  class mdp_log
  {
  private:
    int m_level;
    int m_max_level;
    std::vector<std::string> m_level_tag;
    std::ostream *m_os;

  protected:
    bool m_print;

  public:
    mdp_log() : m_level(0), m_max_level(100000), m_print(true)
    {
      connect(std::cout);
    }

    void abort()
    {
      exit(-1);
    }

    void set_level(int i)
    {
      m_max_level = i;
    }

    void connect(std::ostream &os1)
    {
      m_os = &os1;
    }

    void connect(std::ofstream &os2)
    {
      m_os = &os2; // is this correct? I think so!
    }

    bool printing() const
    {
      return m_print;
    }

    void disablePrinting()
    {
      m_print = false;
    }

    void enablePrinting()
    {
      m_print = true;
    }

    void restorePrinting(bool value)
    {
      m_print = value;
    }

    void error_message(std::string s, std::string file = "unkown", int line = 0)
    {
      if (m_print)
      {
        begin_function("error");
        *m_os << "In file \"" << file;
        *m_os << "\", before line " << line;
        *m_os << ", this error occurred: " << s << "\n";
        for (; m_level; m_level--)
          if (m_level < m_max_level)
            *m_os << "</" << m_level_tag[m_level - 1] << ">\n";
      }
      throw s;
    }

    void begin_function(std::string s)
    {
      m_level_tag.resize(++m_level);
      m_level_tag[m_level - 1] = s;
      if (m_print && m_level < m_max_level)
        *m_os << "<" << s << ">\n";
    }

    void end_function(std::string s)
    {
      if (m_level_tag[m_level - 1] == s)
      {
        if (m_print && m_level < m_max_level)
          *m_os << "</" << m_level_tag[m_level - 1] << ">\n";
        m_level--;
      }
      else
        error_message("missing end_function()", "unkown", 0);
    }

    template <class T>
    mdp_log &operator<<(const T x)
    {
      if (m_print && m_level < m_max_level)
        *m_os << x;
      return (*this);
    }
  };
} // namespace MDP

#endif /* MDP_LOG_ */
