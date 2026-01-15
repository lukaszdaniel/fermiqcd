/////////////////////////////////////////////////////////////////
/// @file mdp_utils.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Other junk that did not fit anywhere else
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_UTILS_
#define MDP_UTILS_

#include <functional>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#ifndef _WIN32
#include <glob.h> // for glob_t
#endif
#include "mdp_field.h"

namespace fs = std::filesystem;

namespace MDP
{
  bool file_exists(const std::string &filename)
  {
    return fs::exists(filename);
  }

  std::vector<std::string> glob(const std::string &pattern)
  {
    std::vector<std::string> v;

#ifndef _WIN32
    glob_t pglob;
    pglob.gl_offs = 2;

    if (glob(pattern.c_str(), 0, nullptr, &pglob) != 0)
    {
      v.push_back("?");
    }
    else
    {
      for (size_t i = 0; i < pglob.gl_pathc; ++i)
        v.emplace_back(pglob.gl_pathv[i]);
    }

    globfree(&pglob);

#else
    fs::path p(pattern);
    fs::path dir = p.parent_path();
    std::string mask = p.filename().string();

    if (dir.empty())
      dir = ".";

    if (!fs::exists(dir) || !fs::is_directory(dir))
    {
      v.push_back("?");
      return v;
    }

    std::function<bool(const char *, const char *)> match = [&](const char *str, const char *pat) -> bool
    {
      while (*pat)
      {
        if (*pat == '*')
        {
          ++pat;
          if (!*pat)
            return true;
          while (*str)
            if (match(str++, pat))
              return true;
          return false;
        }
        else if (*pat == '?')
        {
          if (!*str)
            return false;
          ++str;
          ++pat;
        }
        else
        {
          if (*str != *pat)
            return false;
          ++str;
          ++pat;
        }
      }
      return *str == '\0';
    };

    bool found = false;
    for (const auto &entry : fs::directory_iterator(dir))
    {
      const std::string name = entry.path().filename().string();
      if (match(name.c_str(), mask.c_str()))
      {
        v.push_back(entry.path().string());
        found = true;
      }
    }

    if (!found)
      v.push_back("?");

#endif

    return v;
  }

  std::string latest_file(const std::string &pattern)
  {
    std::vector<std::string> v = glob(pattern);
    return v[v.size() - 1];
  }

  std::string next_to_latest_file(std::string pattern)
  {
    int i = pattern.find("*");
    if (i < 0)
      return pattern;
    std::vector<std::string> v = glob(pattern);
    std::string latest = v[v.size() - 1];
    if (latest == "?")
      return pattern.replace(i, 1, std::to_string(0));

    latest = latest.substr(i, 5);
    int k = 1 + atoi(latest.c_str());

    return pattern.replace(i, 1, std::to_string(k).c_str());
  }

  mdp_field_file_header get_info(const std::string &filename, int proc = 0)
  {
    mdp_field_file_header myheader;
    if (isSubProcess(proc))
    {
      std::ifstream in(filename, std::ios::binary);
      if (!in)
        error("Unable to open file");

      in.read(reinterpret_cast<char *>(&myheader), sizeof(mdp_field_file_header));

      if (!in)
        error("Error while reading file");

      mdp_field_file_header::switch_header_endianess(myheader);
    }
    mdp.broadcast(myheader, proc);
    return myheader;
  }

  int mail(const std::string &email, const std::string &message)
  {
    static int ret;
    if (isMainProcess())
    {
      std::string s = "echo '" + message + "' | mail -s 'MDP MESSAGE' " + email;
      ret = system(s.c_str());
    }
    mdp.broadcast(ret, 0);
    return ret;
  }

  int mail_file(const std::string &email, const std::string &filename)
  {
    static int ret;
    if (isMainProcess())
    {
      std::string s = "more " + filename + " | mail -s 'MDP MESSAGE' " + email;
      ret = system(s.c_str());
    }
    mdp.broadcast(ret, 0);
    return ret;
  }

  bool startswith(const std::string &a, const std::string &b)
  {
    return (a.substr(0, b.length()) == b);
  }

  bool endswith(const std::string &a, const std::string &b)
  {
    int i = a.length();
    int j = b.length();
    if (i < j)
      return false;
    return (a.substr(i - j, j) == b);
  }

  int parse_int(const std::string &a, const std::string &b, int value = 0)
  {
    int i = a.find(std::string(":") + b);
    if (i < 0)
      return value;
    else
    {
      i += b.length() + 2;
      int j = a.find(":", i);
      if (j < 0)
        j = a.length();
      sscanf(a.substr(i, j - i).c_str(), "%i", &value);
      return value;
    }
  }

  mdp_uint parse_uint(const std::string &a, const std::string &b, mdp_uint value = 0)
  {
    int i = a.find(std::string(":") + b);
    if (i < 0)
      return value;
    else
    {
      i += b.length() + 2;
      int j = a.find(":", i);
      if (j < 0)
        j = a.length();
      sscanf(a.substr(i, j - i).c_str(), "%u", &value);
      return value;
    }
  }

  float parse_float(const std::string &a, const std::string &b, float value = 0.0)
  {
    int i = a.find(std::string(":") + b);
    if (i < 0)
      return value;
    else
    {
      i += b.length() + 2;
      int j = a.find(":", i);
      if (j < 0)
        j = a.length();
      sscanf(a.substr(i, j - i).c_str(), "%f", &value);
      return value;
    }
  }

  std::string parse_string(const std::string &a, const std::string &b, const std::string &value = "")
  {
    int i = a.find(std::string(":") + b);
    if (i < 0)
      return value;
    else
    {
      i += b.length() + 2;
      int j = a.find(":", i);
      if (j < 0)
        j = a.length();
      return a.substr(i, j - i);
    }
  }
} // namespace MDP

#endif /* MDP_UTILS_ */
