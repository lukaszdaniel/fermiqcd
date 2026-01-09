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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include <iomanip>

namespace MDP
{
  /// @brief base class of class mdp_communicator (DO NOT INSTANTIATE)
  /// @see class mdp_communicator
  class mdp_log
  {
  private:
    bool m_printEnabled;
    bool m_footerEnabled;
    bool m_cpuValid;

    std::vector<std::string> m_level_tag;

    std::ostream *m_os;
    std::ofstream m_file;

    std::chrono::system_clock::time_point m_timeStart;
    std::clock_t m_cpuStart;

    mdp_log(const mdp_log &) = delete;
    mdp_log &operator=(const mdp_log &) = delete;
    mdp_log(mdp_log &&) = default;
    mdp_log &operator=(mdp_log &&) = default;

    static std::string formatTime(const std::chrono::system_clock::time_point &tp)
    {
      std::time_t t = std::chrono::system_clock::to_time_t(tp);
      std::tm tm{};

      if (auto *ptm = std::localtime(&t))
        tm = *ptm;
      else
        return "time not available";

      std::ostringstream oss;
      oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
      return oss.str();
    }

    void openTarget(const std::string &target)
    {
      if (target.empty() || target == "std::cout")
      {
        m_os = &std::cout;
      }
      else if (target == "std::stderr")
      {
        m_os = &std::cerr;
      }
      else
      {
        m_file.open(target);
        if (!m_file)
          throw std::runtime_error("cannot open log file: " + target);
        m_os = &m_file;
      }
    }

    void closeTarget()
    {
      if (m_file.is_open())
        m_file.close();
    }

    // --- nagłówek / stopka --------------------------------------------------

    void writeHeader(int argc = 0, char *argv[] = nullptr, const std::string &id = "")
    {
      m_footerEnabled = true;

      m_timeStart = std::chrono::system_clock::now();
      m_cpuStart = std::clock();
      m_cpuValid = (m_cpuStart != static_cast<std::clock_t>(-1));

      (*m_os) << "# program    " << (argv ? argv[0] : "unknown") << "   " << id << "\n";
      (*m_os) << "# time start " << formatTime(m_timeStart) << "\n";

      if (argv)
      {
        for (int i = 0; i < argc; ++i)
          (*m_os) << "# argv " << std::setw(2) << i << " " << argv[i] << "\n";
      }

      (*m_os) << "# header end\n";
      m_os->flush();
    }

    void writeFooter()
    {
      if (m_cpuValid)
      {
        double cpuTime = static_cast<double>(std::clock() - m_cpuStart) / CLOCKS_PER_SEC;

        (*m_os) << "# CPU time  " << std::setw(12) << cpuTime << " sec\n";
      }

      auto timeEnd = std::chrono::system_clock::now();
      double wallTime = std::chrono::duration<double>(timeEnd - m_timeStart).count();

      (*m_os) << "# user time " << std::setw(12) << wallTime << " sec\n";
      (*m_os) << "# time stop " << formatTime(timeEnd) << "\n";

      m_os->flush();
    }

  public:
    class function_scope
    {
    public:
      function_scope(mdp_log &log, std::string name)
          : m_log(log), m_name(std::move(name))
      {
        m_log.begin_function(m_name);
      }

      ~function_scope() noexcept
      {
        m_log.end_function(m_name);
      }

      function_scope(const function_scope &) = delete;
      function_scope &operator=(const function_scope &) = delete;

    private:
      mdp_log &m_log;
      std::string m_name;
    };

    mdp_log() : m_printEnabled(true), m_footerEnabled(false), m_cpuValid(false), m_os(nullptr)
    {
      openTarget("std::cout");
    }

    mdp_log(const std::string &target, int argc, char *argv[], const std::string &id = "") : m_printEnabled(true), m_footerEnabled(false), m_cpuValid(false), m_os(nullptr)
    {
      openTarget(target);
      writeHeader(argc, argv, id);
    }

    mdp_log(const std::string &target) : mdp_log(target, 0, nullptr, "")
    {
    }

    ~mdp_log()
    {
      if (m_footerEnabled)
      {
        writeFooter();
        closeTarget();
      }
    }

    void enablePrinting() { m_printEnabled = true; }
    void disablePrinting() { m_printEnabled = false; }
    bool printingEnabled() const { return m_printEnabled; }

    void setPrinting(bool value)
    {
      m_printEnabled = value;
    }

    [[noreturn]]
    void error_message(const std::string &message,
                       const std::string &file = "unknown",
                       int line = 0)
    {
      if (m_printEnabled)
      {
        begin_function("error");

        (*m_os)
            << "In file \"" << file
            << "\", before line " << line
            << ", this error occurred: "
            << message << "\n";

        while (m_level_tag.size() > 0)
        {
          (*m_os) << "</" << m_level_tag.back() << ">\n";
          m_level_tag.pop_back();
        }

        m_os->flush();
      }

      throw std::runtime_error(message);
    }

    void begin_function(const std::string &name)
    {
      m_level_tag.push_back(name);

      if (m_printEnabled)
        (*m_os) << "<" << name << ">\n";
      m_os->flush();
    }

    void end_function(const std::string &name)
    {
      if (m_level_tag.size() == 0)
        return;

      if (m_level_tag.back() != name)
        error_message("missing end_function()", "unknown", 0);

      if (m_printEnabled)
        (*m_os) << "</" << name << ">\n";

      m_level_tag.pop_back();
      m_os->flush();
    }

    void flush()
    {
      if (m_printEnabled)
        m_os->flush();
    }

    template <typename T>
    mdp_log &operator<<(const T &value)
    {
      if (m_printEnabled)
        (*m_os) << value;
      return *this;
    }
  };
} // namespace MDP

#endif /* MDP_LOG_ */
